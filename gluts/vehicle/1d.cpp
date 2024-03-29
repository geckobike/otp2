#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ENABLE_TIMER
#include "engine/engine.h"

#include <GL/gl.h>
#include <GL/glut.h>

#include "springs.h"

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#define clampf(x,min,max) ((x) < (min)) ? (min) : ( (x) > (max) ? (max) : x )

static int s_width = 480 * 2;
static int s_height = 320 * 2;

static bool s_paused = false;
static float s_averageDt = 0.f;

#define XYZ(v) (v).x, (v).y, (v).z
#define XYZp(v) (v)->x, (v)->y, (v)->z

static void drawCircle(float x, float y, float r, bool filled = true);

//======================================================================================================================================================
//======================================================================================================================================================

static Timer g_time;

static void dlDrawLine(const vec3* from, const vec3* to)
{
	glBegin(GL_LINES);
	glVertex3f(from->x, from->y, from->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}

static void drawCircle(float x, float y, float r, bool filled)
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



//======================================================================================================================================================
//======================================================================================================================================================

struct Wheel;
struct Chassis;


struct Wheel
{
    Chassis* chassis;
	vec3 offset;

    float mass;

	float spring;	// This will be calculated from the rest len
	float damp;

	float lowTravel;
	float highTravel;
	float restLen;	// This is the natural unconstrained length, if travel is not restricted (will be clamped to >=travel)

	vec3 vel;
	vec3 pos;
};

struct Chassis
{
    Wheel* wheel;

    float mass;
    vec3 pos;
	vec3 vel;
};


const int numVehicles = 1;
static Chassis s_chassis[numVehicles];
static Wheel s_wheel[numVehicles];
static float gravity = 10.f;


static void vehicleDraw()
{
	for (int i=0; i<numVehicles; i++)
	{
		Chassis* c = &s_chassis[i];

		float cx = c->pos.x;
		float cy = c->pos.y;
		float cz = c->pos.z;

		glColor4f(1,1,1,1);
		drawCircle(cx, cz, 0.1f);
		glBegin(GL_LINES);
		glVertex3f(cx - 0.5f, cz, 0.f);
		glVertex3f(cx + 0.5f, cz, 0.f);
		glEnd();

		Wheel* w = c->wheel;
		drawCircle(w->pos.x, w->pos.z, 0.05f);

		// Draw the travel limits
		glColor4f(1,0,0,1);
		glBegin(GL_LINES);
		glVertex3f(cx - 0.1f, cz + w->offset.z + w->highTravel, 0.f);
		glVertex3f(cx + 0.1f, cz + w->offset.z + w->highTravel, 0.f);
		glVertex3f(cx - 0.1f, cz + w->offset.z - w->lowTravel, 0.f);
		glVertex3f(cx + 0.1f, cz + w->offset.z - w->lowTravel, 0.f);
		glEnd();
	}
}

static void wheelTick(Wheel* wheel, float dt)
{
}

static void chassisTick(Chassis* chassis, float dt)
{
}

static void vehicleSubTick(Chassis* c, float dt)
{
	const int numIterations = 1;
	const float invIterations = 1.f/(float)(numIterations);
	
	Wheel* w = c->wheel;

	float diff = w->pos.z - (c->pos.z + w->offset.z);
	float compression = diff + w->restLen;

	// Add gravity
	c->vel.z -= gravity * dt;
	w->vel.z -= gravity * dt;

	// Calc the force/velocity that the spring will apply
	float dv = w->vel.z - c->vel.z;
	float requiredVelocity = physicsSpringGetVelocity(compression, dv, dt, w->spring, w->damp);	// This is not technically correct, but will do!

	// Do we hit the ground? Work out the time spent in contact with the ground, and the time spent without contact with the ground
	float collisionDistance = w->pos.z;
	if (collisionDistance < 0.f) collisionDistance = 0.f;
	float wheelSpeed = w->vel.z + requiredVelocity - dv;
	float timeNoCollision;
	if (wheelSpeed >= 0.f)					// There's a lot of ifs here!
	{
		if (collisionDistance<=0.01f)
		{
			timeNoCollision = 0.f;
		}
		else
		{
			timeNoCollision = dt;
		}
	}
	else
	{
		timeNoCollision = -(collisionDistance / wheelSpeed);
	}
	if (timeNoCollision > dt) timeNoCollision = dt;
	const float timeInCollision = dt - timeNoCollision;

	// Finally work out how much to push the wheel and chassis by
	const float correctionSpeed = (requiredVelocity - dv);
	const float velocityCorrectionChassis = -correctionSpeed * invIterations * (timeInCollision/dt);
	const float velocityCorrectionWheel = correctionSpeed * invIterations * (timeNoCollision/dt);

	// Calculate the min/max compression lengths
	const float minCompression = w->restLen - w->lowTravel;
	const float maxCompression = w->restLen + w->highTravel;

	// Iterate
	float collisionError = 0.f;
	for (int repeat = 0; repeat < numIterations; repeat++)
	{
		// Apply the spring
		if (1)
		{
			c->vel.z += velocityCorrectionChassis;
			w->vel.z += velocityCorrectionWheel;
		}

		// Clamp the spring
		float compressionError = 0.f;
		if (1)
		{
			float compression2 = compression - (c->vel.z * dt) + (w->vel.z * dt);
			if (compression2>maxCompression)
			{
				compressionError = (compression2 - maxCompression) / dt;
			}
			if (compression2<minCompression)
			{
				compressionError = (compression2 - minCompression) / dt;
			}

			// Add to chassis, subtract from wheel
			// dv = impulse / m;
			float denom = 1.f/c->mass + 1.f/w->mass;
			float correctionImpulse = compressionError / denom;

			c->vel.z += correctionImpulse / c->mass;
			w->vel.z -= correctionImpulse / w->mass;
		}

		// Ground collision (again)
		if (1)
		{
			collisionError = 0.f - w->pos.z - w->vel.z*dt;
			if (collisionError < 0.f) collisionError = 0.f;

			collisionError = collisionError / dt;

			w->vel.z += collisionError;
			if (compressionError>0.f)
			{
				c->vel.z += collisionError;
			}
		}
	}

	// Finally, integrate
	{
		c->pos.z += c->vel.z * dt;
		w->pos.z += w->vel.z * dt;
	}
}

static void vehicleSubTick(float dt)
{
	for (int i=0; i<numVehicles; i++)
	{
		vehicleSubTick(&s_chassis[i], dt);
	}
}

const float subDt = 1.0/30.f;//0.01f;

static void vehicleTick(float dt)
{
	// Fixed time step integration
	static float remaining = 0.f;
	remaining += dt;
	while (remaining>0.f)
	{
		remaining -= subDt;
		vehicleSubTick(subDt);
	}
}

static void vehicleReset(Chassis* chassis)
{
	Chassis* c = chassis;

	c->mass = 1000.f;
	vecset(&c->pos, 0.f, 0.f, 1.f);
	veczero(&c->vel);

	Wheel* w = c->wheel;
	w->chassis = c;
	w->mass = 1.0f;
	w->spring = 1.f;
	w->damp = 10.f;
	w->lowTravel = 0.1f;	// 20 cm up and down from offset
	w->highTravel = 0.1f;	// 20 cm up and down from offset
	w->restLen = 0.3f;
	veczero(&w->vel);
	vecset(&w->offset, 0.f, 0.f, -1.f);
	vecadd(&w->pos, &c->pos, &w->offset);

	// Auto Calc suspension
	if (w->lowTravel < 0.05f) w->lowTravel = 0.05f;
	if (w->highTravel < 0.01f) w->highTravel = 0.01f;
	if (w->restLen < w->lowTravel) w->restLen = w->lowTravel;

	// To pre calculate spring we need the subdt that we will work with! Pooh!
	w->spring = gravity * (1.f + subDt*w->damp) / (w->restLen - gravity * subDt * subDt);
}

static void vehicleReset()
{
	s_chassis[0].wheel = &s_wheel[0];

	vehicleReset(&s_chassis[0]);

	s_chassis[0].pos.x = -1.f;
	s_wheel[0].pos.x = -1.f;
}

static void vehicleInit()
{
	vehicleReset();
}



//======================================================================================================================================================
//======================================================================================================================================================

static void keyboard(unsigned char key, int x, int y)
{
	if (key=='r') vehicleReset();
	if (key==27) exit(0);
	if (key==' ')
	{
		for (int i=0; i<numVehicles; i++)
		{
			Chassis* c = &s_chassis[i];
			Wheel* w= c->wheel;
			const float jumpSpeed = 5.f;
			if (c->vel.z < jumpSpeed)
			{
				c->vel.z += jumpSpeed;
				w->vel.z += jumpSpeed;
			}
		}
	}
}

static void animateScene ()
{
	float dt = timerGetTickSinceLastUUpdate(&g_time);
	timerUpdate(&g_time);

	vehicleTick(dt);

    glutPostRedisplay();
}

static void display()
{
    // clear the window
    glClearColor (0.3,0.3,0.3,0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // go to GL_MODELVIEW matrix mode and set the camera
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    // leave openGL in a known state - flat shaded white, no textures
    glDisable (GL_TEXTURE_2D);
    glShadeModel (GL_FLAT);
    glDisable (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);
    
    glTranslatef(0.f, -2.f, -5.5f);

	// Draw the Ground Plane (2D)
    glColor4f (1.f,1.f,1.f,1.f);
	glBegin(GL_LINES);
		glVertex3f(-20.f,0.f,0.f);
		glVertex3f(+20.f,0.f,0.f);
	glEnd();

	vehicleDraw();
    
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
	timerUpdate(&g_time);
    vehicleInit();

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

    return 0;
}


