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
};


struct Bone
{
    float radius;
    float radius2;
    mtx pose;
    vec3 p;
};

static const u32 N = 30;
static float s_sleeveRadius = 0.1f;
static u32 s_numVerts = N;
static u32 s_numLinks = N;
static vec3 s_verts[N];
static vec3 s_vertsDiff[N];
static vec3 s_lastVerts[N];
static vec3 s_debugVerts[N];
static VertInfo s_vertInfos[N];
static link s_links[N];
static Bone s_bone;
static vec3 s_debugPos = {0.f, 0.f, 0.f};

static void setupLink(link* ln, u32 v1, u32 v2)
{
    ln->v1 = v1;
    ln->v2 = v2;
    ln->restLen = vecdist(&s_verts[v1], &s_verts[v2]);
}

static void start()
{
    s_bone.radius = 0.05f;
    s_bone.radius2 = SQR(s_bone.radius);
    matrixIdent(&s_bone.pose);
    veczero(&s_bone.p);

    float dAlpha = PIx2 / (float)s_numVerts;
    float a = 0.0f;
    for (u32 i=0; i<s_numVerts; i++)
    {
	vecset(&s_verts[i], sinf(a), cosf(a), 0.f);
	vecscale(&s_verts[i], s_sleeveRadius);
	veczero(&s_vertsDiff[i]);
	veccpy(&s_lastVerts[i], &s_verts[i]);
	veccpy(&s_debugVerts[i], &s_verts[i]);
	a += dAlpha;
    }

    for (u32 i=0; i<(N-1); i++)
    {
	setupLink(&s_links[i], i, i+1);
    }
    setupLink(&s_links[N-1], N-1, 0);

    for (u32 i=0; i<N; i++)
    {
	VertInfo& vi = s_vertInfos[i];
	
	vi.maxDist = (s_sleeveRadius - s_bone.radius) * 0.8f;
	vi.lastMaxDist = vi.maxDist;
	veccpy(&vi.defaultPos, &s_verts[i]);
    }

}

static void keyboard(unsigned char key, int x, int y)
{
    if (key == 'p') s_paused = !s_paused;
    if (key == '<') s_slow *= 0.8f;
    if (key == '>') s_slow *= 1.2f;

    if (s_slow < 0.1f) s_slow = 0.1f;
    if (s_slow > 2.f) s_slow = 2.f;
}

static void constrainLink(link* ln, float weight)
{
    vec3* v1 = &s_verts[ ln->v1 ];
    vec3* v2 = &s_verts[ ln->v2 ];

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

static float distSqPointToLine(vec3* linepos, vec3* dir, vec3* point)
{
    vec3 toPoint, pointOnLine;
    vecsub(&toPoint, point, linepos);
    float dot=vecdot(&toPoint, dir)/vecsizesq(dir);
    vecadd(&pointOnLine, linepos, dir);
    return vecdistsq(point, &pointOnLine);
}

static vec3 s_grav = {0.f, -2.f, 0.f};

static void doPositionConstraints()
{
    // Position Constraints
    vec3 avpos = {0.f, 0.f, 0.f};
    for (u32 i=0; i<s_numVerts; i++)
    {
        VertInfo& vi = s_vertInfos[i];
        vec3* v = &s_verts[i];
        vecadd(&avpos, v);
    }
    vecscale(&avpos, 1.0f/(float)s_numVerts);

    vec3 up;
    vecsub(&up, &s_bone.p, &avpos);

    if (vecsizesq(&up) > 0.00001f)
    {
        vecnormalise(&up);
    }

    s_debugPos = avpos;
    
    // Position Constraints
    for (u32 i=0; i<s_numVerts; i++)
    {
	VertInfo& vi = s_vertInfos[i];
	float maxDist = vi.maxDist;
	vec3* v = &s_verts[i];
	vec3 diff;
	veccpy(&diff, &vi.defaultPos);
	vecnormalise(&diff);

	float f = vecdot(&diff, &up);
	float a = 2.3f;
	float p = 1+a - a*f;//*f*f;
	if (f > 0.5f)
	{
	    maxDist *= p;
	    maxDist = 0.7f*maxDist + 0.3f*vi.lastMaxDist;
	}
	else
	{
	    maxDist = maxDist*10.f;
	}

	vec3 dpos;
	vecadd(&dpos, &s_bone.p, &vi.defaultPos);
	positionConstrain(v, &dpos, maxDist);

	vi.lastMaxDist = maxDist;
    }
    
}

static void doPositionConstraints2()
{
    // Position Constraints
    vec3 avpos = {0.f, 0.f, 0.f};
    for (u32 i=0; i<s_numVerts; i++)
    {
        VertInfo& vi = s_vertInfos[i];
        vec3* v = &s_verts[i];
        vecadd(&avpos, v);
    }
    vecscale(&avpos, 1.0f/(float)s_numVerts);

    vec3 up;
    vecsub(&up, &s_bone.p, &avpos);

    if (vecsizesq(&up) > 0.00001f)
    {
        vecnormalise(&up);
    }

    s_debugPos = avpos;
    
    // Position Constraints
    for (u32 i=0; i<s_numVerts; i++)
    {
	VertInfo& vi = s_vertInfos[i];
	float maxDist = vi.maxDist;
	vec3* v = &s_verts[i];
	vec3 diff;
	veccpy(&diff, &vi.defaultPos);
	vecnormalise(&diff);

	float f = vecdot(&diff, &up);
	float a = 2.3f;
	float p = 1+a - a*f;//*f*f;
	if (f > 0.5f)
	{
	    maxDist *= p;
	    maxDist = 0.7f*maxDist + 0.3f*vi.lastMaxDist;
	}
	else
	{
	    maxDist = maxDist*10.f;
	}

	vec3 dpos;
	vecadd(&dpos, &s_bone.p, &vi.defaultPos);
	positionConstrain(v, &dpos, maxDist);

	vi.lastMaxDist = maxDist;
    }
    
}

static void animateScene ()
{
    if (s_paused) return;

    static float remaining = 0.f;
    float dt = averageTime() * s_slow;
    printf("fps = %3.1f\n", 1/dt);

    float subdt = 1.f/240.f;

    remaining += dt;

    // Move the bone
    {
	static float time = PI_D2;
	float speed = 0.5f;
	time+=dt*speed;

	float radius = 0.3f;

	s_bone.p.x = 3.0f*radius*cosf(time);
	s_bone.p.y = 1.0f*radius*sinf(2*time);
	
	s_bone.p.x *= cosf(5.f*time);
	s_bone.p.y *= sinf(5.f*time);
    }
    while (remaining > 0.f)
    {
	remaining -= subdt;
    
	// Make a tmp copy of verts (needed for last position)
	vec3 tmp[N];
	memcpy(tmp, s_verts, sizeof(vec3)*N);

	vec3 grav;
	vecscale(&grav, &s_grav, subdt);
	vecscale(&grav, &s_grav, subdt);

	// Vertlet
	for (u32 i=0; i<s_numVerts; i++)
	{
	    vec3 tmp = {0.f, 0.f, 0.f};
	    vec3* v = &s_verts[i];
	    vec3* old = &s_lastVerts[i];

	    //vecsub(&tmp, v, old);
	    //vecscale(&tmp, 0.95f);    // Time varying stablisation

	    vecadd(&tmp, &grav);
	    vecadd(v, &tmp);

	    veccpy(&s_debugVerts[i], v);
	    vec3* vel = &s_vertsDiff[i];
	}

	for (u32 repeat=0; repeat<3; repeat++)
	{
	    // Link Constraints
	    for (u32 i=0; i<s_numLinks; i++)
	    {
		constrainLink(&s_links[i], 0.5f);
	    }

	    doPositionConstraints2();

	}

	// Make a copy of the last position
	memcpy(s_lastVerts, tmp, sizeof(vec3)*N);

    }


    wait(15000);
    //wait(16000/s_slow);
    
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

static void drawCloth(vec3* v)
{
    for (u32 i=0; i < s_numVerts; i++)
    {
	drawPoint(v[i].x, v[i].y);
    }
    
    glBegin(GL_LINE_STRIP);
    for (u32 i=0; i < (s_numLinks + 1); i++)
    {
	vec3* v1 = &v[ s_links[i % s_numLinks].v1 ];
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

    //lastViewX += 0.05f*(s_bone.p.x - lastViewX);
    //lastViewY += 1.f*(s_bone.p.y - lastViewY);
    
    glTranslatef(-lastViewX, -lastViewY, -1.5f);

    // Draw the orgin
    glColor3f(0.f, 0.f, 1.f);
    drawCircle( 0.0f, -0.1f, 0.01f);
    drawCircle( 0.0f, +0.1f, 0.01f);
    drawCircle(+0.1f,  0.0f, 0.01f);
    drawCircle(-0.1f,  0.0f, 0.01f);
    drawCircle( 0.0f,  0.0f, 0.02f, false);

    // Draw the bone
    glColor3f(0.f, 0.f, 0.f);
    drawCircle(s_bone.p.x, s_bone.p.y, s_bone.radius);

    // Draw the cloth
    glColor3f(1.f, 1.f, 1.f);
    drawCloth(s_verts);
    
    // Debug
    glColor3f(0.f, 1.0f, 1.f);
    //drawCloth(s_debugVerts);

    drawPoint(s_debugPos.x, s_debugPos.y, 0.02f);
    
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


