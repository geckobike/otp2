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

#define rnd rand
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define assert(x) { if (!(x)) { printf("ABORT: from file %s, line %d: expr failed: %s\n", __FILE__, __LINE__, #x); fflush(stdout); exit(0); } }

#define veclist(v) (v)[0], (v)[1], (v)[2], (v)[3]

typedef unsigned short u16;
typedef unsigned char u8;

static int s_width = 800;
static int s_height = 600;

static bool s_paused = false;
static float s_averageDt = 0.f;

#if 1
static float s_elev = 0.858333;
static float s_yaw = -0.017500;
static float s_cameraRadius = 2.899999;
static vec3 s_cameraPos = {0.899089, 0.700000, 0.040486};
static vec3 s_cameraDir;
static vec3 s_cameraForward = {-0.017499, 0.000000, -0.999847};
static vec3 s_cameraLeft = {-0.999847, 0.000000, 0.017499};
#else
static float s_elev=0.f;
static float s_yaw=0.f;
static float s_cameraRadius=1.f;
static vec3 s_cameraPos={0.f, 1.f, 0.f};
static vec3 s_cameraDir;
static vec3 s_cameraForward = {0.f, 0.f, 1.f};	// normalised direction in XZ plane
static vec3 s_cameraLeft = {-1.f, 0.f, 0.f};	// normalised direction in XZ plane
#endif
static int s_mouseDown = false;
static int s_lastMouseX=0;
static int s_lastMouseY=0;




#define ARRAY_SIZEOF(a) ((int)( sizeof((a)) / sizeof((a)[0]) ))

/*
=======================================================================================
    Utils
=======================================================================================
*/

template <class T>
static void sortQuickSort(T* array, int l, int r)	// use: sortQuickSort(array, 0, N-1), need to implement "operator < ()"
{
	int i,j, mi;
	T k, m;

	i = l; j = r; mi = (l + r)/2;
	m = array[mi];
	while (i <= j) {
		while((i<=j) && array[i] < m) i++;
		while((j>l) && m < array[j]) j--;
		if (i <= j) {
			k = array[i]; array[i] = array[j]; array[j] = k;
			i++; j--;
		}
	}
	if (l < j) sortQuickSort<T>(array, l, j);
	if (i < r) sortQuickSort<T>(array, i, r);
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

void dlDrawLine(const vec3* from, const vec3* to)
{
	glBegin(GL_LINES);
	glVertex3f(from->x, from->y, from->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}


static long long getTimeLong()
{
    struct timeval now;
    gettimeofday(&now,NULL);
    return  now.tv_sec*1000000 + now.tv_usec;
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
	glVertex3f(x-size, 0.f, y-size);
	glVertex3f(x-size, 0.f, y+size);
	glVertex3f(x+size, 0.f, y+size);
	glVertex3f(x+size, 0.f, y-size);
    glEnd();
}


void updateCamera(const vec3* zz, const vec3* yy)
{
    vec3 up={0.f, 1.f, 0.f};
    vec3 forward={0.f, 0.f, -1.f};

    mtx rot;
    matrixRotY(&rot, s_yaw);
    vec3mtx33mulvec3(&s_cameraForward, &rot, &forward);

    s_cameraForward.y = 0.f;
    s_cameraForward.z = -cosf(s_yaw);
    s_cameraForward.x = +sinf(s_yaw);

    s_cameraDir = s_cameraForward;

    veccross(&s_cameraLeft, &up, &s_cameraForward);
    matrixRot(&rot, &s_cameraLeft, s_elev);
    vec3mtx33mulvec3(&s_cameraDir, &rot, &s_cameraForward);

    vec3 target = s_cameraPos;
    vec3 pos;
    vecsubscale(&pos, &target, &s_cameraDir, s_cameraRadius);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(pos.x, pos.y, pos.z, target.x, target.y, target.z, up.x, up.y, up.z);
}


/*
========================================================================================================
	Poly
========================================================================================================
*/
union vec2
{
    struct
    {
	float x,y;
    };
    struct
    {
	float v[2];
    };
};

typedef struct Poly
{
    //float scale;
    //float invScale;
    int n;
    vec2* verts;
    vec2* diffs;
    float* k;
};

static vec2 s_verts[] =
{
    {0.f, 0.f},
    {-0.5f, 1.f},
    {+1.5f, 1.2f},
    {+2.0f, -0.5f},
};

static Poly s_poly;

static float cross2(const vec2* a, const vec2* b)
{
    return (a->x*b->y) - (a->y*b->x);
}

bool intersectPointInConcavePoly2D_nobranch(const vec2 *point, vec2 *const vertexes, const int numvertexes)
{
    int i;
    const vec2 *v1 = &vertexes[numvertexes - 1];
    boolean inside = FALSE;

    for (i = 0; i < numvertexes; i++)
    {
	vec2 *v2 = &vertexes[i];
	float d = (point->y - v1->y)*(v2->x - v1->x);
	float tmp = (point->x - v1->x)*(v2->y - v1->y);
	inside ^= point->y < v2->y & point->y > v1->y & d > tmp;
	inside ^= point->y > v2->y & point->y < v1->y & d < tmp;
	v1 = v2;
    }

    return inside;
}

boolean intersectPointInConcavePoly2D(const vec2 *point, vec2 *const vertexes, const int numvertexes)
{
    int i;
    const vec2 *v1 = &vertexes[numvertexes - 1];
    boolean inside = FALSE;

    for (i = 0; i < numvertexes; i++)
    {
	vec2 *v2 = &vertexes[i];

	if (point->y < v2->y)
	{
	    if (v1->y <= point->y)
	    {
		if ((point->y - v1->y)*(v2->x - v1->x) > (point->x - v1->x)*(v2->y - v1->y))
		{
		    inside = !inside;
		}
	    }
	}
	else
	if (point->y < v1->y)
	{
	    if ((point->y - v1->y)*(v2->x - v1->x) < (point->x - v1->x)*(v2->y - v1->y))
	    {
		inside = !inside;
	    }
	}

	v1 = v2;
    }

    return inside;
}

static bool inside_nob(float x, float y, Poly* poly)
{
    vec2 p = {x,y};
    return intersectPointInConcavePoly2D(&p, poly->verts, poly->n);
}


static bool inside(float x, float y, Poly* poly)
{
    for (int i=0; i<poly->n; i++)
    {
	const vec2* v0 = poly->verts + (i+0)%poly->n;
	const vec2* v1 = poly->verts + (i+1)%poly->n;
	const vec2* d = poly->diffs + i;

	float dx = d->x;
	float dy = d->y;

	if ((((x*dy - y*dx) < poly->k[i]))) return false;
    }
    return true;
}

static void drawPoly()
{
    glDisable(GL_DEPTH_TEST);
    glColor3f(0.f, 0.6f, 0.f);
    
    int n = ARRAY_SIZEOF(s_verts);

    for (int i=0; i<n; i++)
    {
	vec2* v0 = s_poly.verts + i;
	vec2* v1 = s_poly.verts + (i+1)%n;
	vec3 p0 = {v0->x, 0.f, v0->y};
	vec3 p1 = {v1->x, 0.f, v1->y};
	dlDrawLine(&p0, &p1);
    }

    float dy = 0.03f;
    float dx = 0.03f;
    float y=-0.6f;
    for ( ; y< 1.5f; y+=dy)
    {
	float x=-0.6f;
	for ( ; x< 1.5f; x+=dx)
	{
	    if (inside(x,y, &s_poly))
	    {
		glColor3f(1.f, 1.f, 0.f);
	    }
	    else
	    {
		glColor3f(0.f, 0.f, 1.f);
	    }
	    drawPoint(x,y, 0.01f);
	}
    }

    long long t1 = getTimeLong();
    
    for (u32 rep=0; rep<1; rep++)
    {
	float dy = 0.03f;
	float dx = 0.03f;
	float y=-0.6f;
	for ( ; y< 1.5f; y+=dy)
	{
	    float x=-1.f;
	    for ( ; x< 1.5f; x+=dx)
	    {
		bool d = inside(x,y, &s_poly);
	    }
	}
    }

    long long t2 = getTimeLong();
    float delta = (float)(t2-t1);

    printf("time = %f\n", delta);

    {
	long long t1 = getTimeLong();

	for (u32 rep=0; rep<100; rep++)
	{
	    float dy = 0.03f;
	    float dx = 0.03f;
	    float y=-0.6f;
	    for ( ; y< 1.5f; y+=dy)
	    {
		float x=-1.f;
		for ( ; x< 1.5f; x+=dx)
		{
		    bool d = inside_nob(x,y, &s_poly);
		}
	    }
	}

	long long t2 = getTimeLong();
	float delta = (float)(t2-t1);

	printf("nob time = %f\n", delta);
    }
}

static void initPoly()
{
    int n = ARRAY_SIZEOF(s_verts);
    s_poly.n=n;
    s_poly.verts = s_verts;
    s_poly.diffs = (vec2*)malloc(n*sizeof(s_poly.diffs[0]));
    s_poly.k = (float*)malloc(n*sizeof(s_poly.k[0]));

    for (int i=0; i<n; i++)
    {
	s_poly.diffs[i].x = s_verts[(i+1)%n].x - s_verts[i].x;
	s_poly.diffs[i].y = s_verts[(i+1)%n].y - s_verts[i].y;
	s_poly.k[i] = cross2(s_verts+i, s_poly.diffs+i);
    }
}


/*
========================================================================================================
	Demo
========================================================================================================
*/

void render()
{
    updateCamera(&s_cameraPos, &s_cameraDir);
    
    u32 N = 10;
    float width=7.f;
    float dx = width/(float)(N-1);
    float dy = width/(float)(N-1);
    float y=-width*0.5f;
    float x=-width*0.5f;
    
    glColor3f(0.2f, 0.2f, 0.2f);
    for (u32 i=0; i<N; i++)
    {
	vec3 from = {x, 0.f, -y};
	vec3 to = {x, 0.f, +y};
	dlDrawLine(&from, &to);
	x += dx;
    }
    
    y=-width*0.5f;
    x=-width*0.5f;
    for (u32 i=0; i<N; i++)
    {
	vec3 from = {+x, 0.f, y};
	vec3 to = {-x, 0.f, y};
	dlDrawLine(&from, &to);
	y += dy;
    }


    drawPoly();
}


static void start()
{
    initPoly();

    GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    /*	light_position is NOT default value	*/
    GLfloat light_position0[] = { 1.0, 10.0, 1.0, 0.0 };
    GLfloat light_position1[] = { -1.0, -10.0, -1.0, 0.0 };
  
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glClearColor(0.3,0.3,0.3,0);

    //  glEnable(GL_CULL_FACE);
    //  glCullFace(GL_BACK);
}

static void keyboard(unsigned char key, int x, int y)
{
    if (key==27) exit(0);

    vec3 up={0.f, 1.f, 0.f};

    if (key=='q') vecaddscale(&s_cameraPos, &s_cameraPos, &up, -0.1f);
    if (key=='e') vecaddscale(&s_cameraPos, &s_cameraPos, &up, +0.1f);

    if (key=='w')
    {
	s_cameraRadius -= 0.1f;
    }
    if (key=='s')
    {
	s_cameraRadius += 0.1f;
    }

    if (key=='a')
    {
	vecaddscale(&s_cameraPos, &s_cameraPos, &s_cameraLeft, +0.1f);
    }
    if (key=='d')
    {
	vecaddscale(&s_cameraPos, &s_cameraPos, &s_cameraLeft, -0.1f);
    }

    glutPostRedisplay();

    printf("static float s_elev = %f;\n", s_elev);
    printf("static float s_yaw = %f;\n", s_yaw);
    printf("static float s_cameraRadius = %f;\n", s_cameraRadius);
    printf("static vec3 s_cameraPos = {%f, %f, %f};\n", vec3list(s_cameraPos));
    printf("static vec3 s_cameraDir;\n");
    printf("static vec3 s_cameraForward = {%f, %f, %f};\n", vec3list(s_cameraForward));
    printf("static vec3 s_cameraLeft = {%f, %f, %f};\n", vec3list(s_cameraLeft));
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

    render();

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
    const float vnear = 0.01f;
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
    glutPostRedisplay();
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
    glutCreateWindow ("triangles to quad mesh");

    // Initialize OpenGL graphics state
    initGraphics();

    // Register callbacks:
    glutDisplayFunc (display);
    glutReshapeFunc (reshape);
    glutKeyboardFunc (keyboard);
    glutMouseFunc (mouseButton);
    glutMotionFunc (mouseMotion);
    //glutIdleFunc (animateScene);

    // Turn the flow of control over to GLUT
    glutMainLoop ();

    //printf("average dt was = %f, and FPS = %f\n", s_averageDt, 1.f/s_averageDt);

    return 0;
}


