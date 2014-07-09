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

static float s_elev=0.f;
static float s_yaw=0.f;
static float s_cameraRadius=1.f;
static vec3 s_cameraPos={1.f, 1.f, 2.f};
static vec3 s_cameraDir;
static vec3 s_cameraForward = {0.f, 0.f, 1.f};	// normalised direction in XZ plane
static vec3 s_cameraLeft = {-1.f, 0.f, 0.f};	// normalised direction in XZ plane
static int s_mouseDown = false;
static int s_lastMouseX=0;
static int s_lastMouseY=0;

#define EXTERN_INLINE

EXTERN_INLINE vec3 operator + (const vec3& a, const vec3& b)
{
    vec3 tmp;
    tmp.x = a.x + b.x;
    tmp.y = a.y + b.y;
    tmp.z = a.z + b.z;
    return tmp;
}

EXTERN_INLINE vec3 operator - (const vec3& a, const vec3& b)
{
    vec3 tmp;
    tmp.x = a.x - b.x;
    tmp.y = a.y - b.y;
    tmp.z = a.z - b.z;
    return tmp;
}

EXTERN_INLINE vec3 operator * (float s, const vec3& a)
{
    vec3 tmp;
    tmp.x = s * a.x;
    tmp.y = s * a.y;
    tmp.z = s * a.z;
    return tmp;
}

EXTERN_INLINE vec3 operator * (const vec3& a, float s)
{
    vec3 tmp;
    tmp.x = s * a.x;
    tmp.y = s * a.y;
    tmp.z = s * a.z;
    return tmp;
}

EXTERN_INLINE vec3 operator / (const vec3& a, float s)
{
    float f = (1.f/s);
    vec3 tmp;
    tmp.x = a.x * f;
    tmp.y = a.y * f;
    tmp.z = a.z * f;
    return tmp;
}


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
	glVertex3f(x-size, y-size, 0.f);
	glVertex3f(x-size, y+size, 0.f);
	glVertex3f(x+size, y+size, 0.f);
	glVertex3f(x+size, y-size, 0.f);
    glEnd();
}

static void drawPoint3d(float x, float y, float z, float size = drawEps)
{
    glBegin(GL_QUADS);
	glVertex3f(x-size, y-size, z+size);
	glVertex3f(x-size, y+size, z+size);
	glVertex3f(x+size, y+size, z+size);
	glVertex3f(x+size, y-size, z+size);

	glVertex3f(x-size, y-size, z-size);
	glVertex3f(x-size, y+size, z-size);
	glVertex3f(x+size, y+size, z-size);
	glVertex3f(x+size, y-size, z-size);

	glVertex3f(x-size, y-size, z-size);
	glVertex3f(x-size, y-size, z+size);
	glVertex3f(x-size, y+size, z+size);
	glVertex3f(x-size, y+size, z-size);

	glVertex3f(x+size, y-size, z-size);
	glVertex3f(x+size, y-size, z+size);
	glVertex3f(x+size, y+size, z+size);
	glVertex3f(x+size, y+size, z-size);

	glVertex3f(x-size, y+size, z-size);
	glVertex3f(x-size, y+size, z+size);
	glVertex3f(x+size, y+size, z+size);
	glVertex3f(x+size, y+size, z-size);

	glVertex3f(x-size, y-size, z-size);
	glVertex3f(x-size, y-size, z+size);
	glVertex3f(x+size, y-size, z+size);
	glVertex3f(x+size, y-size, z-size);
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


bool s_active = true;

static void animateScene ()
{
    glutPostRedisplay();
}


////////////////////////////////////////


static const int nrows = 10;
static const int ncols = 10;
static float s_heights[nrows][ncols] = {0.f};

////////////////////////////////////////

/*
    
    0---------> X axis
    |
    |
    |
    |  Z axis
   \./ 

    a---D1->---b
    |         /|
    |       /  |
    D2    /    |
   \|/  /      |
    | /        |
    c----------d	a,b,c,d are the heights at those grid positions

*/

static bool linetestA(vec3* out, const vec3* start, const vec3* dir, float a, float b, float c)
{
    //vec3 normal;
    //vec3 d1 = {1, b-a, 0.f};
    //vec3 d2 = {0.f, c-a, 1.f};
    //veccross(&normal, &d2, &d1);

    vec3 normal = {a-b, 1.f, a-c};
    float t = vecdot(dir, &normal);
    if (t<0.f)
    {
	vec3 diff = {start->x, start->y - a, start->z};
	float t2 = -vecdot(&normal, &diff);
	if (t2>=t)   // equiv: t <= 1.f
	{
	    t2 = t2 / t;
	    *out = *start + *dir * t2;

	    if ((out->z<(1.f-out->x)) && (out->x > 0.f) && (out->x <= 1.f) && (out->z > 0.f) && (out->z <= 1.f))
	    {
		return true;
	    }
	    return false;
	}
    }
    return false;
}



void render()
{
    animateScene();
    
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

    static float time = 0.f;
    time += 1.f/60.f;

    s_heights[0][0] = 0.5f;
    s_heights[0][1] = 0.0;//5f*sinf(0.01f*time);

    {
	vec3 from;
	vec3 to;

	glColor3f(1.0f, 1.0f, 1.0f);
	
	float a = s_heights[0][0];
	float b = s_heights[0][1];
	float c = s_heights[1][0];
	float d = s_heights[1][1];

	vecset(&from, 0.f, a, 0.f);
	vecset(&to, 1.f, b, 0.f);
	dlDrawLine(&from, &to);

	vecset(&from, 0.f, a, 0.f);
	vecset(&to, 0.f, c, 1.f);
	dlDrawLine(&from, &to);

	vecset(&from, 1.f, b, 0.f);
	vecset(&to, 0.f, c, 1.f);
	dlDrawLine(&from, &to);

	vecset(&from, 1.f, d, 1.f);
	vecset(&to, 0.f, c, 1.f);
	dlDrawLine(&from, &to);

	vecset(&from, 1.f, d, 1.f);
	vecset(&to, 1.f, b, 0.f);
	dlDrawLine(&from, &to);

	vec3 out;
	vec3 start = {2.f, 0.124999f, 0.5f};
	vec3 dir = {-2.901f, 0.f, 0.f};
	to = start + 1.8f*dir;

	bool hit = linetestA(&out, &start, &dir, a, b, c);
	
	if (hit)
	{
	    glColor3f(1.0f, 0.0f, 0.0f);
	    dlDrawLine(&start, &to);
	    drawPoint3d(out.x, out.y, out.z, 0.01f);
	}
	else
	{
	    glColor3f(0.0f, 1.0f, 0.0f);
	    dlDrawLine(&start, &to);
	}
    }

}


static void start()
{
    //createLinks(&gMesh, &gEdges);

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

    if (key==' ')
    {
	s_active = !s_active;
    }

    glutPostRedisplay();
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
    glutIdleFunc (animateScene);

    // Turn the flow of control over to GLUT
    glutMainLoop ();

    //printf("average dt was = %f, and FPS = %f\n", s_averageDt, 1.f/s_averageDt);

    return 0;
}


