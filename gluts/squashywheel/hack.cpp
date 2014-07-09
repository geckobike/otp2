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

    float dt = averageTime();
    //printf("fps = %3.1f\n", 1/dt);
    dt = 1.f/30.f;
   
	static float time = 0.f;

	s_time += dt*s_speed;

	Wheel* w=&s_wheel;
	w->distortion = sqr(sinf(s_time*1.f));
	if (w->distortion>0.f)
		w->distortion = sqrtf(w->distortion);

	w->distortion *= 0.4f;

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
	glVertex3f(x-size, y-size, z);
	glVertex3f(x-size, y+size, z);
	glVertex3f(x+size, y+size, z);
	glVertex3f(x+size, y-size, z);
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

static void drawWheel(Wheel* w, float z)
{
	float distortion = w->distortion;

	const float x = w->pos.x;
	const float y = w->pos.y - distortion;
	const float r = w->radius;

    // Draw the orgin
	drawWheelCentre(x,y,0.1f);

    const u32 N = 64;
	vec3 verts[N*5];
    float dAlpha = PIx2 / (float)N;
    float angle = s_time;
	for (u32 i=0; i<N; i++)
	{
	    vecset(&verts[i], r*cosf(angle), r*sinf(angle), z);
		vecscale(&verts[i + N], &verts[i], 3.f/5.f);
		vecscale(&verts[i + 2*N], &verts[i], 3.5f/5.f);
		vecscale(&verts[i + 3*N], &verts[i], 4.f/5.f);
		vecscale(&verts[i + 4*N], &verts[i], 4.5f/5.f);
		verts[i + N].z = z;
		verts[i + 2*N].z = z;
		verts[i + 3*N].z = z;
		verts[i + 4*N].z = z;
	    angle += dAlpha;
	}
			
	// Radial distortion
	if (distortion>0.f)
	{
		for (u32 i=0; i<N*5; i++)
		{
			vec3* p0 = &verts[i];
			if (p0->y>0.f)
				continue;



#if 0		// Simple squash - based on x coord :(
			float xx = fabsf(p0->x)/r;
			float yy = fabsf(p0->y)/r;
			p0->y += (1.f-xx*xx*xx)*distortion*r;
#endif

#if 0 		// Simple squash based on depth - reasonable
			float xx = fabsf(p0->x)/r;
			float yy = fabsf(p0->y)/r;
			if (yy > (1.f-distortion))
			{
				float s = 2.f-yy-distortion;
				//float s = (1-distortion)/yy;
				p0->x *= s;
				p0->y *= s;
			}
#endif

#if 1 		// Simple squash based on depth - reasonable
			float vr = sqrtf(sqr(p0->x) + sqr(p0->y));

			float overall = max(0.f, 1.f - (r-vr)/(r-0.6f));


			float xx = fabsf(p0->x) / vr;
			float yy = fabsf(p0->y) / vr;
			float excess = yy - (1.f-distortion);
			if (excess>0.f)
			{
				//float s = 2.f-max(yy,0.f*xx)-distortion; // approx of: float s = (1-distortion)/yy;
				float s = 1.f - overall * max(excess,0.f);
				p0->x *= s + overall * max(0.f,excess)*distortion*xx;
				p0->y *= s;
			}
			else
			{
				float s = 1.f + overall * sqr(distortion)*(excess - (distortion - 1.f));
				p0->x *= s;
				p0->y *= 1.f - overall * sqr(distortion*excess);
			}

			{
				float bell = 1.f/(1.f+10.f*sqr(sqr(xx))) - 1.f/11.f;

				float yyy = min(1.f, yy/(1.f - distortion));
				float y2 = sqr(yyy);
				float y4 = sqr(y2);
				float y8 = sqr(y4);
				float bulge = y2 - y8;

				float s = 1.f + 15.f*overall * sqr(distortion)*(excess - (distortion - 1.f))*bulge;

				//p0->z *= 1.f + 5.f*distortion * overall * bulge * bell;
				p0->z *= s;
			}
#endif
	
		}
	}
		
	for (u32 i=0; i<N*5; i++)
	{
		verts[i].x+=x; verts[i].y+=y;
	}

	for (u32 i=0; i<N; i++)
	{
		vec3 p0 = verts[i];
		vec3 p1 = verts[(i+1)%N];
		dlDrawLine(&p0, &p1);
		drawPoint(p0,0.02f);
	}
	
	for (u32 i=N; i<5*N; i++)
	{
		vec3 p0 = verts[i];
		drawPoint(p0,0.02f);
	}

	//const TriMesh m = GetWheelObj1();
	//for (int i=0; i<min(10,m.numTris); i++)
	//{
	//	const int i0 = m.indices[i*3+0] - 1; // obj file, so index from 1
	//	const int i1 = m.indices[i*3+1] - 1;
	//	const int i2 = m.indices[i*3+2] - 1;
	//	vec3 v0 = m.verts[i0];
	//	vec3 v1 = m.verts[i1];
	//	vec3 v2 = m.verts[i2];
	//	dlDrawLine(&v0,&v1);
	//	dlDrawLine(&v1,&v2);
	//	dlDrawLine(&v2,&v0);
	//	//printf("%f,%f,%f\n", v0.x,v0.y,v0.z);
	//}

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

	float posX = 10.f*cosf(s_yaw)*cosf(s_elev);
	float posY = 10.f*sinf(s_elev);
	float posZ = 10.f*sinf(s_yaw)*cosf(s_elev);
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

    //lastViewX += 0.05f*(b.p.x - lastViewX);
    //lastViewY += 1.f*(b.p.y - lastViewY);
    
  //  glTranslatef(-lastViewX, -lastViewY, -10.f);
//	glRotatef(s_yaw, 0,1,0);
//	glRotatef(s_elev, 1,0,0);

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


