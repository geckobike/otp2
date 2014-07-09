#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ENABLE_TIMER
#include "engine/engine.h"

#include <GL/gl.h>
#include <GL/glut.h>

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

typedef unsigned short u16;
typedef unsigned char u8;

static int s_width = 800;
static int s_height = 600;

static bool s_paused = false;
static float s_averageDt = 0.f;

static float s_elev=0.f;
static float s_yaw=0.f;
static float s_cameraRadius=1.f;
static vec3 s_cameraPos={0.f, 1.f, 2.5f};
static vec3 s_cameraDir;
static vec3 s_cameraForward = {0.f, 0.f, -1.f};	// normalised direction in XZ plane
static vec3 s_cameraLeft = {-1.f, 0.f, 0.f};	// normalised direction in XZ plane
static Timer g_time;
static int s_mouseDown = false;
static int s_lastMouseX=0;
static int s_lastMouseY=0;
static float g_angle = 0.f;

// Approximation of exp(-x)
inline float approxExp(float x)
{
	return 1.f/(1.f+x);
}
	
// Approximation of 1.f - exp(-x)
inline float approxOneExp(float x)
{
	return x/(1.f+x);
}

float fsgnf(float x)
{
	return x>=0.f ? +1.f : -1.f;
}

#define ARRAY_SIZEOF(a) ((int)( sizeof((a)) / sizeof((a)[0]) ))


/*
=======================================================================================
    Utils
=======================================================================================
*/


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

void dlDrawCircle(float x, float y, float r, bool filled = true)
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

void dlMultMatrix(const mtx *mat)
{
    glMultMatrixf((float *)mat);
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


// static void dlDrawCross(float size, vec4* colour)
// {
// 	vec3 p0 = {+size,0.f,0.f};
// 	vec3 p1 = {-size,0.f,0.f};
// 	vec3 p2 = {0.f,+size,0.f};
// 	vec3 p3 = {0.f,-size,0.f};
// 	vec3 p4 = {0.f,0.f,+size};
// 	vec3 p5 = {0.f,0.f,-size};
// 	dlDrawLine(&p0, &p1, colour);
// 	dlDrawLine(&p2, &p3, colour);
// 	dlDrawLine(&p4, &p5, colour);
// }
// 
// static void BodyDraw(Body* body, vec4* colour)
// {
// 	mtx m;
// 	quaternionToMatrix(&m, &body->rot);
// 	veccpy((vec3*)&m.v[3], &body->pos);
// 	glPushMatrix();
// 	dlMultMatrix(&m);
// 	dlDrawCross(0.02f, colour);
// 	glPopMatrix();
// }
// 
// 
// static void CapsuleBodyDraw(CapsuleBody* cap, vec4* colour)
// {
// 	mtx m;
// 	quaternionToMatrix(&m, &cap->body.rot);
// 	veccpy((vec3*)&m.v[3], &cap->body.pos);
// 	glPushMatrix();
// 	dlMultMatrix(&m);
// 	{
// 	vec3 pos0 = {+cap->radius, -0.5f*cap->length, 0.f};
// 	vec3 pos1 = {+cap->radius, +0.5f*cap->length, 0.f};
// 	vec3 pos2 = {-cap->radius, -0.5f*cap->length, 0.f};
// 	vec3 pos3 = {-cap->radius, +0.5f*cap->length, 0.f};
// 	vec3 pos4 = {0.f, -0.5f*cap->length - cap->radius, 0.f};
// 	vec3 pos5 = {0.f, +0.5f*cap->length + cap->radius, 0.f};
// 	dlDrawLine(&pos0, &pos1, colour);
// 	dlDrawLine(&pos1, &pos5, colour);
// 	dlDrawLine(&pos5, &pos3, colour);
// 	dlDrawLine(&pos3, &pos2, colour);
// 	dlDrawLine(&pos2, &pos4, colour);
// 	dlDrawLine(&pos4, &pos0, colour);
// 	}
// 	{
// 	vec3 pos0 = {0.f, -0.5f*cap->length, +cap->radius};
// 	vec3 pos1 = {0.f, +0.5f*cap->length, +cap->radius};
// 	vec3 pos2 = {0.f, -0.5f*cap->length, -cap->radius};
// 	vec3 pos3 = {0.f, +0.5f*cap->length, -cap->radius};
// 	vec3 pos4 = {0.f, -0.5f*cap->length - cap->radius, 0.f};
// 	vec3 pos5 = {0.f, +0.5f*cap->length + cap->radius, 0.f};
// 	dlDrawLine(&pos0, &pos1, colour);
// 	dlDrawLine(&pos1, &pos5, colour);
// 	dlDrawLine(&pos5, &pos3, colour);
// 	dlDrawLine(&pos3, &pos2, colour);
// 	dlDrawLine(&pos2, &pos4, colour);
// 	dlDrawLine(&pos4, &pos0, colour);
// 	}
// 	dlDrawCross(0.02f, colour);
// 	glPopMatrix();
// }

struct TrafficLight
{
	enum
	{
		k_red=0, k_yellow, k_green, k_yellow2, k_numStates
	};

	Vec3 m_pos;
	float m_timer;
	int m_state;
	float m_stateTimes[k_numStates];

	void Init(Vec3 pos, float timeRed, float timeGreen)
	{
		const float lightChangeTime = 2.0f;
		m_pos = pos;
		m_state = k_red;
		m_stateTimes[k_red] = timeRed;
		m_stateTimes[k_yellow] = lightChangeTime;
		m_stateTimes[k_green] = timeGreen;
		m_stateTimes[k_yellow2] = lightChangeTime;
		m_timer = 0.f;
	}

	void Update(float dt)
	{
		m_timer += dt;
		if (m_timer >= m_stateTimes[m_state]) // Change State
		{
			m_timer -= m_stateTimes[m_state];
			m_state = (m_state + 1) % k_numStates;
		}
	}

	int GetState() { return m_state; }
};

struct Car
{
	Vec3 pos;
	Vec3 vel;
};


struct Simulation
{
public:
	void Init();
	void Draw();
	void Update(float dt);
public:
	TrafficLight m_light;
};


void Simulation::Init()
{
	Vec3 p0 = {0.f, 0.f, 0.f};
	m_light.Init(p0, 60.0, 60.0);
}

void Simulation::Draw()
{
	vec4 colour = {1.f, 1.f, 1.f, 1.f};
	vec4 red = {1.f, 0.f, 0.f, 1.f};
	vec4 green = {0.f, 1.f, 0.f, 1.f};
	vec4 blue = {0.f, 0.f, 1.f, 1.f};
	vec3 a0 = {0.f, 0.f, 0.f};
	vec3 ax = {1.f, 0.f, 0.f};
	vec3 ay = {0.f, 1.f, 0.f};
	vec3 az = {0.f, 0.f, 1.f};
	dlDrawLine(&a0, &ax, &red);
	dlDrawLine(&a0, &ay, &green);
	dlDrawLine(&a0, &az, &blue);

	vec4 cols[4] = 
	{
		{1.f, 0.f, 0.f, 1.f}, // red
		{1.f, 1.f, 1.f, 1.f},
		{0.f, 1.f, 0.f, 1.f}, // green
		{0.f, 1.f, 1.f, 1.f},
	};

}


void Simulation::Update(float dt)
{
	m_light.Update(dt);
}


Simulation g_simulation;

static void start()
{
	g_simulation.Init();

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

    // if (key=='w')
    // {
	// s_cameraRadius -= 0.1f;
    // }
    // if (key=='s')
    // {
	// s_cameraRadius += 0.1f;
    // }

	vec3 forward = s_cameraForward;
	forward.y = 0.f;
	vecnormalise(&forward, &forward);
    if (key=='w')
		vecaddscale(&s_cameraPos, &s_cameraPos, &forward, +0.1f);
    if (key=='s')
		vecaddscale(&s_cameraPos, &s_cameraPos, &forward, -0.1f);

    if (key=='a')
    {
	vecaddscale(&s_cameraPos, &s_cameraPos, &s_cameraLeft, +0.1f);
    }
    if (key=='d')
    {
	vecaddscale(&s_cameraPos, &s_cameraPos, &s_cameraLeft, -0.1f);
    }

    glutPostRedisplay();
}

static void animateScene()
{
	float dt = timerGetTickSinceLastUUpdate(&g_time);
	if (dt < 10e-3) return;
	timerUpdate(&g_time);
	
	dt = MIN(50e-3f, dt);
	
	static float remaining = 0.f;
	remaining += dt;
	while (remaining>0.f)
	{
		float subDt = 0.01f;
		g_simulation.Update(subDt);
		remaining -= subDt;
	}

	// keep ticking
    glutPostRedisplay();
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

	g_simulation.Draw();
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
	
	timerUpdate(&g_time);

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


