#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine/engine.h"

#include <GL/gl.h>
#include <GL/glut.h>

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

static int s_width = 800;
static int s_height = 600;

static bool s_paused = false;
static float s_averageDt = 0.f;
static float s_slow = 1.0f;
static float s_speed = 1.0f; // speed of bone

static double averageTime()
{
    static bool init = false;
    static struct timeval lastTime;

    if (!init) { gettimeofday(&lastTime,NULL); init = true; }
    
    struct timeval now;
 
    long long sec;
    long long microsec;
    
    gettimeofday(&now,NULL);
    
    sec = now.tv_sec - lastTime.tv_sec;
    microsec = now.tv_usec - lastTime.tv_usec;
    
    lastTime = now;
   
    // Filter the average time to smooth it out

    double newDt = (double)sec + ( (double)microsec / (1000.0*1000.0) );

    return s_averageDt += 0.5*(newDt - s_averageDt);
}

static void wait(long long delay) // in micro secs
{
    struct timeval last;
    struct timeval now;
    
    gettimeofday(&last,NULL);
    gettimeofday(&now,NULL);

    while((now.tv_usec - last.tv_usec) < delay)
    {
	gettimeofday(&now,NULL);
	if (now.tv_sec > last.tv_sec)
	{
	    now.tv_usec += 1000000;
	}
    }
}

//=================================================================================================================

struct link
{
    u32 v1;
    u32 v2;
    float restLen;
};

struct VertInfo
{
    float maxDist;  // Max dist to travel to bone
    float lastMaxDist;
    vec3 defaultPos;
    float dist0;
    float dotA;
};


struct Bone
{
    float radius;
    float radius2;
    vec3 p;
    vec3 lastP;
};

static const u32 N = 20;

struct Sleeve
{
    float radius;
    vec3 verts[N];
    vec3 lastVerts[N];
    VertInfo vertInfos[N];
    link links[N];
    u32 numVerts;
    u32 numLinks;
};

struct Arm
{
    Bone bone;
    Sleeve sleeve;
};


static const u32 s_numArms = 50;
static Arm s_arm[s_numArms];


static void setupLink(link* ln, u32 v1, u32 v2, vec3* verts)
{
    ln->v1 = v1;
    ln->v2 = v2;
    ln->restLen = vecdist(&verts[v1], &verts[v2]);
}

static void start()
{
    for (u32 a = 0; a < s_numArms; a++)
    {
	Sleeve& s = s_arm[a].sleeve;
	Bone& b = s_arm[a].bone;

	s.radius = 0.3f;
	s.numVerts = N;
	s.numLinks = N;

	b.radius = 0.05f;
	b.radius2 = SQR(b.radius);

	float dAlpha = PIx2 / (float)N;
	float a = 0.0f;
	for (u32 i=0; i<s.numVerts; i++)
	{
	    vecset(&s.verts[i], sinf(a), cosf(a), 0.f);
	    vecscale(&s.verts[i], s.radius);
	    veccpy(&s.lastVerts[i], &s.verts[i]);
	    a += dAlpha;
	}

	for (u32 i=0; i<(N-1); i++)
	{
	    setupLink(&s.links[i], i, i+1, s.verts);
	}
	setupLink(&s.links[N-1], N-1, 0, s.verts);

	for (u32 i=0; i<N; i++)
	{
	    VertInfo& vi = s.vertInfos[i];

	    vi.maxDist = (s.radius - b.radius) * 0.9f;
	    vi.lastMaxDist = vi.maxDist;
	    vecsub(&vi.defaultPos, &s.verts[i], &b.p);
	    vi.dist0 = vecsize(&vi.defaultPos);
	    vi.dotA = 1.0f;
	}
    }
}

static void keyboard(unsigned char key, int x, int y)
{
    if (key == 'p') s_paused = !s_paused;
    if (key == '<') s_slow *= 0.8f;
    if (key == '>') s_slow *= 1.2f;

    if (s_slow < 0.1f) s_slow = 0.1f;
    if (s_slow > 2.f) s_slow = 2.f;
   
    if (key == ',') s_speed *= 0.8f;
    if (key == '.') s_speed *= 1.2f;
   
    if (s_speed < 0.01f) s_speed= 0.01f;
    if (s_speed > 4.f) s_speed = 4.f;
}

static void constrainLink(link* ln, vec3* verts, float weight)
{
    vec3* v1 = &verts[ ln->v1 ];
    vec3* v2 = &verts[ ln->v2 ];

    vec3 diff, diff0, offset;

    vecsub(&diff, v2, v1);
    
    float size=vecsize(&diff);	// Try to remove this sqrtf
    vecscale(&diff0, &diff, ln->restLen/size);
    vecsub(&offset, &diff, &diff0);
     
    vecaddscale(v1, v1, &offset, weight);
    vecsubscale(v2, v2, &offset, weight);
}


static void positionConstrain(vec3* vert, vec3* pos, float maxDist)
{
    vec3 diff, diff0;

    vecsub(&diff, vert, pos);
    
    float size=vecsize(&diff);	// Try to remove this sqrtf
    if (size > maxDist)
    {
	vecscale(&diff0, &diff, maxDist/size);
	vecadd(vert, pos, &diff0);
    }
}

//static float distSqPointToLine(vec3* linepos, vec3* dir, vec3* point)
//{
//    vec3 toPoint, pointOnLine;
//    vecsub(&toPoint, point, linepos);
//    float dot=vecdot(&toPoint, dir)/vecsizesq(dir);
//    vecadd(&pointOnLine, linepos, dir);
//    return vecdistsq(point, &pointOnLine);
//}

static vec3 s_grav = {0.f, -2.f, 0.f};

#define EPS_SQ (0.000001f)

// New Method (no sqrts!!)
static void doPositionConstraints2(Arm& arm)
{
    Sleeve& s = arm.sleeve;
    Bone& b = arm.bone;
    
    // Position Constraints
    for (u32 i=0; i<s.numVerts; i++)
    {
	VertInfo& vi = s.vertInfos[i];
	float maxDist = vi.maxDist;
	vec3* v = &s.verts[i];
	vec3* lastv = &s.lastVerts[i];
	vec3& relPos0 = vi.defaultPos;
	
	vec3 relPos;
	vecsub(&relPos, v, &b.p);
	
	vec3 diff;
	vecsub(&diff, &relPos, &relPos0);

	//float dot = vecdot(&diff, &relPos0) / vi.dist0 / vecsize(&diff);
	float dot = vecdot(&diff, &relPos0) / vecsizesq(&diff) * 1.2f;  // This approx eliminates the sqrtf
	
	if (dot < -1.f) dot = -1.f;

	vi.dotA = 2.f;
	float p = 1.f + vi.dotA + vi.dotA*dot;
	
	//vi.dotA = 4.f;
	//float p = 1.f + vi.dotA * (dot+1)*(dot+1);
	
	//float p = 1.f + vi.dotA + vi.dotA * dot * fabsf(dot);
	maxDist*=p;
//	maxDist += 0.9f*(maxDist*p - vi.lastMaxDist);
	
	vec3 pos0;
	vecadd(&pos0, &b.p, &relPos0);
	positionConstrain(v, &pos0, maxDist);

//	vi.lastMaxDist = maxDist;
    }
}

static void animateScene ()
{
    if (s_paused) return;

    float dt = averageTime();
    printf("fps = %3.1f\n", 1/dt);
    dt = 1.f/30.f;
    
    // Move the bones
    {
	static float time = PI_D2;
	time+=dt*s_speed*0.5f;

	float radius = 0.5f;

	for (u32 a=0; a<s_numArms; a++)
	{
	    Sleeve& s = s_arm[a].sleeve;
	    Bone& b = s_arm[a].bone;

	    float dx = 0.2f*cos(PIx2 * a / s_numArms);
	    float dy = 0.2f*sin(PIx2 * a / s_numArms);

	    b.lastP = b.p;

	    b.p.x = 2.0f*radius*cosf(time);
	    b.p.y = 2.0f*radius*sinf(2*time);

	    b.p.x *= cosf(5.f*time);
	    b.p.y *= sinf(5.f*time);

	    b.p.x += dx;
	    b.p.y += dy;
	}
    }

    static float remaining = 0.f;
    float subdt = 1.f/240.f;

    remaining += dt;
	
    // Store the last positions
    for (u32 a=0; a<s_numArms; a++)
    {
	Sleeve& s = s_arm[a].sleeve;
	for (u32 i=0; i<s.numVerts; i++)
	{
	    s.lastVerts[i] = s.verts[i];
	}
    }

    while (remaining > 0.f)
    {
	remaining -= subdt;

	vec3 grav;
	vecscale(&grav, &s_grav, subdt);
	vecscale(&grav, &s_grav, subdt);

	for (u32 a=0; a<s_numArms; a++)
	{
	    Sleeve& s = s_arm[a].sleeve;
	    Bone& b = s_arm[a].bone;

	    // Integrate
	    for (u32 i=0; i<s.numVerts; i++)
	    {
		vec3 tmp = {0.f, 0.f, 0.f};
		vec3* v = &s.verts[i];
		vec3* old = &s.lastVerts[i];

		//vecsub(&tmp, v, old);
		//vecscale(&tmp, 0.95f);    // Time varying stablisation

		vecadd(&tmp, &grav);
		vecadd(v, &tmp);
	    }
	    doPositionConstraints2(s_arm[a]);
	    for (u32 repeat=0; repeat<5; repeat++)
	    {
		// Link Constraints
		for (u32 i=0; i<s.numLinks; i++)
		{
		    constrainLink(&s.links[i], s.verts, 0.5f);
		}
	    }
	}
    }

    wait(16000);
    //wait((long long)(15000.0/s_slow));
    
    glutPostRedisplay();
}

static void drawCircle(float x, float y, float r, bool filled = true)
{
    u32 N = 20;
    float dAlpha = PIx2 / (float)N;
    float a = 0.0f;
    filled ? glBegin(GL_POLYGON) : glBegin(GL_LINE_STRIP);
	for (u32 i=0; i<(N+1); i++)
	{
	    glVertex3f(x+r*cosf(a), y+r*sinf(a), 0.f);
	    a += dAlpha;
	}
    glEnd();
}

static const float drawEps = 0.003f;

static void drawPoint(float x, float y, float size = drawEps)
{
    glBegin(GL_QUADS);
	glVertex3f(x-size, y-size, 0.f);
	glVertex3f(x-size, y+size, 0.f);
	glVertex3f(x+size, y+size, 0.f);
	glVertex3f(x+size, y-size, 0.f);
    glEnd();
}

static void drawCloth(Sleeve& s)
{
    vec3* v = s.verts;
    
    for (u32 i=0; i < s.numVerts; i++)
    {
	drawPoint(v[i].x, v[i].y);
    }
    
    glBegin(GL_LINE_STRIP);
    for (u32 i=0; i < (s.numLinks + 1); i++)
    {
	vec3* v1 = &v[ s.links[i % s.numLinks].v1 ];
	glVertex3f(v1->x, v1->y, 0.f); 
    }
    glEnd();
}

static void display()
{
    // clear the window
    glClearColor (0.5,0.5,0.5,0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // go to GL_MODELVIEW matrix mode and set the camera
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    // leave openGL in a known state - flat shaded white, no textures
    glDisable (GL_TEXTURE_2D);
    glShadeModel (GL_FLAT);
    glDisable (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);
    glColor4f (1.f,1.f,1.f,1.f);

    static float lastViewX = 0.f;
    static float lastViewY = 0.f;

    //lastViewX += 0.05f*(b.p.x - lastViewX);
    //lastViewY += 1.f*(b.p.y - lastViewY);
    
    glTranslatef(-lastViewX, -lastViewY, -1.5f);

    // Draw the orgin
    glColor3f(0.f, 0.f, 1.f);
    drawCircle( 0.0f, -0.1f, 0.01f);
    drawCircle( 0.0f, +0.1f, 0.01f);
    drawCircle(+0.1f,  0.0f, 0.01f);
    drawCircle(-0.1f,  0.0f, 0.01f);
    drawCircle( 0.0f,  0.0f, 0.02f, false);

    // Draw the bone(s)
    for (u32 a=0; a<s_numArms; a++)
    {
	if (a > 0) break;
	Sleeve& s = s_arm[a].sleeve;
	Bone& b = s_arm[a].bone;



	glColor3f(0.f, 0.f, 0.f);
	drawCircle(b.p.x, b.p.y, b.radius);

	// Draw the cloth
	glColor3f(1.f, 1.f, 1.f);
	drawCloth(s);
	// Debug
	glColor3f(0.f, 1.0f, 1.f);
    }
    

    glutSwapBuffers();
}

static void reshape(GLint width, GLint height)
{
    s_width = width;
    s_height = height;
    
    // setup viewport
    glViewport (0,0,s_width,s_height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    const float vnear = 0.1f;
    const float vfar = 100.0f;
    const float k = 0.8f;     // view scale, 1 = +/- 45 degrees
    if (s_width >= s_height) {
	float k2 = float(s_height)/float(s_width);
	glFrustum (-vnear*k,vnear*k,-vnear*k*k2,vnear*k*k2,vnear,vfar);
    }
    else {
	float k2 = float(s_width)/float(s_height);
	glFrustum (-vnear*k*k2,vnear*k*k2,-vnear*k,vnear*k,vnear,vfar);
    }

}

static void initGraphics()
{
}

static void mouseButton(int button, int state, int x, int y)
{
}

static void mouseMotion(int x, int y)
{
}


int main (int argc, char **argv)
{
    start();

    // GLUT Window Initialization:
    glutInit (&argc, argv);
    glutInitWindowSize (s_width, s_height);
    glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow ("CS248 GLUT example");

    // Initialize OpenGL graphics state
    initGraphics();

    // Register callbacks:
    glutDisplayFunc (display);
    glutReshapeFunc (reshape);
    glutKeyboardFunc (keyboard);
    glutMouseFunc (mouseButton);
    glutMotionFunc (mouseMotion);
    glutIdleFunc (animateScene);

    //BuildPopupMenu ();
    //glutAttachMenu (GLUT_RIGHT_BUTTON);

    // Turn the flow of control over to GLUT
    glutMainLoop ();

    printf("average dt was = %f, and FPS = %f\n", s_averageDt, 1.f/s_averageDt);

    return 0;
}


