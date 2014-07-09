#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine/engine.h"

#include <GL/gl.h>
#include <GL/glut.h>

#include "obj.cpp"

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

static int s_width = 800;
static int s_height = 600;

static bool s_paused = false;
static float s_averageDt = 0.f;
static float s_slow = 1.0f;
static float s_speed = 1.0f; // speed of bone
static bool s_mouseDown = false;
static int s_lastMouseX=0;
static int s_lastMouseY=0;
static float s_elev=0.f;
static float s_yaw=0.f;
static float s_camRadius=10.f;

static float min(float x, float y)
{
	return x<y ? x : y;
}

static float max(float x, float y)
{
	return x>y ? x : y;
}

static float sqr(float x)
{
	return x*x;
}

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

struct Wheel
{
	vec3 pos;
	float radius;
	float width;
	float distortion;
};

static Wheel s_wheel;



static void start()
{
	Wheel* w = &s_wheel;
	w->radius = 1.f;
	w->width = 0.25f;
	w->distortion = 0.2f;
	vecset(&w->pos, 0.f, w->radius, 0.f);
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

	if (key =='a' ) s_camRadius += 0.5f;
	if (key =='z' ) s_camRadius -= 0.5f;

	s_camRadius = min(10.f, s_camRadius);
	s_camRadius = max(1.f, s_camRadius);
}


//static float distSqPointToLine(vec3* linepos, vec3* dir, vec3* point)
//{
//    vec3 toPoint, pointOnLine;
//    vecsub(&toPoint, point, linepos);
//    float dot=vecdot(&toPoint, dir)/vecsizesq(dir);
//    vecadd(&pointOnLine, linepos, dir);
//    return vecdistsq(point, &pointOnLine);
//}

#define EPS_SQ (0.000001f)

static float s_time = 0.f;

static void animateScene ()
{
    if (s_paused) return;

    float dt = averageTime() * 2.f;
    //printf("fps = %3.1f\n", 1/dt);
    //dt = 1.f/30.f;
   
	static float time = 0.f;

	s_time += dt*s_speed;

	Wheel* w=&s_wheel;
	w->distortion = sqr(sinf(s_time*1.f));
	if (w->distortion>0.f)
		w->distortion = sqrtf(w->distortion);

	w->distortion *= 0.1f;
	if (w->distortion>0.3f)
		w->distortion = 0.3f;

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

static void drawPoint(vec3 v, float size = drawEps)
{
	float x = v.x, y = v.y, z = v.z;
    glBegin(GL_QUADS);

	glVertex3f(x-size, y-size, z+size);
	glVertex3f(x-size, y+size, z+size);
	glVertex3f(x+size, y+size, z+size);
	glVertex3f(x+size, y-size, z+size);
	
	glVertex3f(x-size, y-size, z-size);
	glVertex3f(x-size, y+size, z-size);
	glVertex3f(x+size, y+size, z-size);
	glVertex3f(x+size, y-size, z-size);
	
	glVertex3f(x-size, y+size, z-size);
	glVertex3f(x-size, y+size, z+size);
	glVertex3f(x+size, y+size, z+size);
	glVertex3f(x+size, y+size, z-size);
	
	glVertex3f(x-size, y-size, z-size);
	glVertex3f(x-size, y-size, z+size);
	glVertex3f(x+size, y-size, z+size);
	glVertex3f(x+size, y-size, z-size);
	
	glVertex3f(x+size, y-size, z-size);
	glVertex3f(x+size, y-size, z+size);
	glVertex3f(x+size, y+size, z+size);
	glVertex3f(x+size, y+size, z-size);
	
	glVertex3f(x-size, y-size, z-size);
	glVertex3f(x-size, y-size, z+size);
	glVertex3f(x-size, y+size, z+size);
	glVertex3f(x-size, y+size, z-size);

    glEnd();
}

void dlDrawLine(const vec3* from, const vec3* to)
{
	glBegin(GL_LINES);
	glVertex3f(from->x, from->y, from->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}

void dlDrawLine(const vec3* from, const vec3* to, const vec4* colour)
{
	glBegin(GL_LINES);
	glColor3f(colour->x, colour->y, colour->z);
	glVertex3f(from->x, from->y, from->z);
	glColor3f(colour->x, colour->y, colour->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}

static void drawWheelCentre(float x, float y, float radius)
{
    // Draw the orgin
    glColor3f(0.f, 0.f, 1.f);
    drawCircle( x, y-radius, 0.01f);
    drawCircle( x, y+radius, 0.01f);
    drawCircle( x+radius, y, 0.01f);
    drawCircle( x-radius, y, 0.01f);
    drawCircle( x, y, 0.02f, false);
}

struct DistortionParams
{
	vec3 distortionNormal;	// only need x,y of this
	float distortion;
	float innerRadius;
	float outerRadius;
	float invRadialDifference;
	float lateralDistortion;
	float lateralDistortion2;
};

static void DistortTyre(vec3* verts, int numVerts, const DistortionParams* params)
{
	const vec3 distortionNormal = params->distortionNormal;
	const float distortion = params->distortion;
	const float innerRadius = params->innerRadius;
	const float outerRadius = params->outerRadius;
	const float invRadialDifference = params->invRadialDifference;
	const float lateralDistortion = params->lateralDistortion;
	const float lateralDistortion2 = params->lateralDistortion2;

	for (int i=0; i<numVerts; i++)
	{
		vec3* p0 = &verts[i];

		float y = (p0->x * distortionNormal.x) + (p0->y * distortionNormal.y);
		if (y>=0.f)
			continue;

		float x = -(p0->x * distortionNormal.y) + (p0->y * distortionNormal.x);

		float outX = x;
		float outY = y;

		// Radial offset of vertex
		float vr = sqrtf(sqr(x) + sqr(y));
		if (vr < innerRadius)
			continue;

		float ivr = 1.f/vr;

		// Distortion scale, based on radial position between innerRadius and outerRadius
		//float overall = min(1.f, 1.f - (outerRadius-vr)*invRadialDifference);
		float overall = 1.f - (outerRadius-vr)*invRadialDifference;

		// Normalise the vertex position
		float xx = fabsf(x) * ivr;
		float yy = fabsf(y) * ivr;
		float excess = yy - (1.f-distortion);	// excess +ve for vertices underneath the ground
		float add;
		if (excess>0.f)
		{
			//float s = 2.f-max(yy,0.f*xx)-distortion; // approx of: float s = (1-distortion)/yy;
			float s = 1.f - overall * max(excess,0.f);
			outX *= s + overall * max(0.f,excess)*distortion*xx;
			outY *= s;
		}
		else
		{
			float s = 1.f + overall * sqr(distortion)*(excess - (distortion - 1.f));
			outX *= s;
			outY *= 1.f - overall * sqr(distortion*excess);
		}

		// divide :(
		float bell = 1.f/(1.f+lateralDistortion2*sqr(sqr(xx))) - 1.f/(1.f+lateralDistortion2);

		{
			float yyy1 = 0.5f*(yy + fabsf(y/outerRadius));
			float yyy2 = min(1.f, fabsf((vr-innerRadius)*invRadialDifference));// * (fabsf(y)+fabsf(x))/outerRadius;
			float yyy3 = min(1.f, fabsf(y)/(1.f - distortion));

			float yyy = (yyy1 + yyy2 + yyy3)*0.333f;

			float y2 = sqr(yyy);
			float y4 = sqr(y2);
			float bulge = y2 - y4;
			p0->z *= 1.f + lateralDistortion*distortion*bell*bulge;
		}

		p0->x = outY * distortionNormal.x - outX * distortionNormal.y;
		p0->y = outX * distortionNormal.x + outY * distortionNormal.y;
	}
}


static void drawWheel(Wheel* w, float z)
{
	DistortionParams params;

	//vecset(&params.distortionNormal, 1.f, 1.f, 0.f);
	vecset(&params.distortionNormal, 0.f, 1.f, 0.f);
	vecnormalise(&params.distortionNormal, &params.distortionNormal);
	params.distortion = w->distortion;
	params.innerRadius = w->radius * 0.6f;
	params.outerRadius = w->radius;
	params.invRadialDifference = 1.f/(params.outerRadius-params.innerRadius);
	params.lateralDistortion = 10.f;
	params.lateralDistortion2 = sqr(7.f);

	// Generate wheel geometry
    const u32 N = 64;
	vec3 verts[N*5];
    float dAlpha = PIx2 / (float)N;
    float angle = s_time*0.2f;
	for (u32 i=0; i<N; i++)
	{
	    vecset(&verts[i], params.outerRadius*cosf(angle), params.outerRadius*sinf(angle), z);
		vecscale(&verts[i + N], &verts[i], 3.2f/5.f);
		vecscale(&verts[i + 2*N], &verts[i], 3.7f/5.f);
		vecscale(&verts[i + 3*N], &verts[i], 4.f/5.f);
		vecscale(&verts[i + 4*N], &verts[i], 4.5f/5.f);
		verts[i + N].z = z;
		verts[i + 2*N].z = z;
		verts[i + 3*N].z = z;
		verts[i + 4*N].z = z;
	    angle += dAlpha;
	}
			
	// Vertex shader distortion
	if (w->distortion>0.f)
	{
		DistortTyre(verts, N*5, &params);
	}
    
	// DRAW
	const float centreX = w->pos.x - params.distortionNormal.x * params.distortion * w->radius;
	const float centreY = w->pos.y - params.distortionNormal.y * params.distortion * w->radius;

	drawWheelCentre(centreX,centreY,0.1f);
		
	for (u32 i=0; i<N*5; i++)
	{
		verts[i].x+=centreX; verts[i].y+=centreY;
	}

	for (u32 i=0; i<N; i++)
	{
		vec3 p0 = verts[i];
		vec3 p1 = verts[(i+1)%N];
		dlDrawLine(&p0, &p1);
		drawPoint(p0,0.01f);
	}
	
	for (u32 i=N; i<5*N; i++)
	{
		vec3 p0 = verts[i];
		drawPoint(p0,0.01f);
	}
}

static void display()
{
    // clear the window
    glClearColor (0.5,0.5,0.5,0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // go to GL_MODELVIEW matrix mode and set the camera
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	float posX = s_wheel.radius*s_camRadius*cosf(s_yaw)*cosf(s_elev);
	float posY = s_wheel.radius*s_camRadius*sinf(s_elev);
	float posZ = s_wheel.radius*s_camRadius*sinf(s_yaw)*cosf(s_elev);
	vec3 up = {0.f, 1.f, 0.f};

    gluLookAt(posX, posY, posZ, 0, 0, 0, up.x, up.y, up.z);


    // leave openGL in a known state - flat shaded white, no textures
    glDisable (GL_TEXTURE_2D);
    glShadeModel (GL_FLAT);
    glDisable (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);
    glColor4f (1.f,1.f,1.f,1.f);

    static float lastViewX = 0.f;
    static float lastViewY = 1.f;

	// Draw ground line
	vec3 p0 = {-5.f, 0.f, 0.f};
	vec3 p1 = {+5.f, 0.f, 0.f};
	vec4 col = {1.f,1.f,1.f,1.f};
	dlDrawLine(&p0, &p1, &col);

	drawWheel(&s_wheel, -0.22f);
	drawWheel(&s_wheel, +0.22f);


	// glColor3f(0.f, 0.f, 0.f);
	// drawCircle(b.p.x, b.p.y, b.radius);

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
    const float k = 0.3f;     // view scale, 1 = +/- 45 degrees
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
    // 0 == down, 1 == up
    int mouseDown = !state;

    if (!s_mouseDown && mouseDown)
    {
		s_lastMouseX = x;
		s_lastMouseY = y;
    }
    s_mouseDown = mouseDown;
}

static void mouseMotion(int x, int y)
{
    if (!s_mouseDown) return;

    float dx = (float)(x - s_lastMouseX)/(float)s_width;
    float dy = (float)(y - s_lastMouseY)/(float)s_height;

    s_lastMouseX = x;
    s_lastMouseY = y;
    
    s_yaw += dx*1.0f;
    s_elev += dy*1.0f;

    if (s_elev > 0.9*PI_D2) s_elev = 0.9*PI_D2;
    if (s_elev < -0.9*PI_D2) s_elev = -0.9*PI_D2;

    glutPostRedisplay();
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

    //printf("average dt was = %f, and FPS = %f\n", s_averageDt, 1.f/s_averageDt);

    return 0;
}


