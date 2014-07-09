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


static float function(float x, float y, float a, float b)
{
    return y + x*x*a + b;
}

static void derivative(float x, y

//=================================================================================================================


static void start()
{
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



static void animateScene ()
{
    if (s_paused) return;

    float dt = averageTime();
    printf("fps = %3.1f\n", 1/dt);
    dt = 1.f/30.f;

    wait(16000);
    
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

//static void drawCloth(Sleeve& s)
//{
//    vec3* v = s.verts;
//    
//    for (u32 i=0; i < s.numVerts; i++)
//    {
//	drawPoint(v[i].x, v[i].y);
//    }
//    
//    glBegin(GL_LINE_STRIP);
//    for (u32 i=0; i < (s.numLinks + 1); i++)
//    {
//	vec3* v1 = &v[ s.links[i % s.numLinks].v1 ];
//	glVertex3f(v1->x, v1->y, 0.f); 
//    }
//    glEnd();
//}

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

#if 0
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
#endif

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


