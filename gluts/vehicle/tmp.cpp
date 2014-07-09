#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ENABLE_TIMER
#include "engine/engine.h"
#include <assert.h>

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


static inline void addImpulseAtOffset(vec3* vel, vec3* angVel, float invMass, float invInertia, const vec3* offset, const vec3* impulse)
{
    vec3 tmp;
	vecaddscale(vel, vel, impulse, invMass);
    veccross(&tmp, offset, impulse);
	vecscale(&tmp, &tmp, invInertia);
	vecadd(angVel, angVel, &tmp);
}

static inline float computeDenominator(float invMass, float invInertia, const vec3* offset, const vec3* norm)
{
    // If you apply an impulse of 1.0f in the direction of 'norm'
    // at position specified by 'offset' then the point will change
    // velocity by the amount calculated here
	vec3 cross;
	veccross(&cross, offset, norm);
	vecscale(&cross, &cross, invInertia);
	veccross(&cross, &cross, offset);
	return vecdot(norm, &cross) + invMass;
}
static inline float getSpringForce(float x, float vel, float dt, float mass, float spring, float damp)
{
	float a = spring * dt;
	return - /*mass* */ (spring * x + (a + damp) * vel);// / (mass  + dt * (a + damp) );
}

static inline float getRestingSpring(float x, float gravity, float dt, float mass, float damp)
{
	// Calc spring constant
	return (mass * gravity + damp * gravity * dt) / (gravity*dt*dt + x);
}

//======================================================================================================================================================
//======================================================================================================================================================

struct Wheel;
struct Chassis;

struct Suspension
{
    Chassis* chassis;
	Wheel* wheel;

	vec3 offset;			// Local to the chassis (this is the top of the suspension spring, and should be inside the chassis' collision geom)

	float spring;			// This will be calculated from the rest len
	float damp;

	float startingCompression;
	float minCompression;
	float maxCompression;
	float linetestLength;

	float vel;
	float compression;			// Current compression value
	float newCompression;

	// Calculated each frame
	vec3 axis;

	// Linetest
	float hitDistance;
	vec3 hitPos;
	vec3 hitNorm;
	float hitPlaneD;

	// Calculated each frame
	vec3 worldOffset;
	vec3 worldDefaultPos;
};


struct Wheel
{
	Suspension* suspension;

    float mass;					// Doesn't really matter!

	// Determined each frame
	vec3 pos;
};

struct Chassis
{
	Suspension* suspension[2];

    float mass;
	float inertia;

	mtx pose;
	vec3 vel;
	float angVel;
};

const int numWheels = 2;
static Chassis s_chassis;
static Wheel s_wheel[numWheels];
static Suspension s_suspension[numWheels];
static float gravity = 10.f;

static void suspensionDraw(Suspension* s)
{
	Chassis* c = s->chassis;
	const vec3* x = &c->pose.v[0].v3;
	const vec3* y = &c->pose.v[1].v3;
	const vec3* z = &c->pose.v[2].v3;
	const vec3* chassisPos = &c->pose.v[3].v3;

	float cx = chassisPos->x;
	float cy = chassisPos->y;
	float cz = chassisPos->z;

	Wheel* w = s->wheel;
	float wx = w->pos.x;
	float wy = w->pos.y;
	float wz = w->pos.z;

	vec3 point;
	vec3mtx43mulvec3(&point, &c->pose, &s->offset);
	
	glColor4f(0,0,1,1);
	drawCircle(point.x, point.z, 0.02f);

	vec3 minPoint, maxPoint;
	vecaddscale(&minPoint, &point, z, -s->startingCompression + s->minCompression);
	vecaddscale(&maxPoint, &point, z, -s->startingCompression + s->maxCompression);

	glColor4f(1,0,0,1);
	drawCircle(minPoint.x, minPoint.z, 0.02f);
	glColor4f(1,0,0,1);
	drawCircle(maxPoint.x, maxPoint.z, 0.02f);

	glColor4f(1,1,1,1);
	drawCircle(wx, wz, 0.05f);
}

static void vehicleDraw()
{
	Chassis* c = &s_chassis;

	vec3* x = &c->pose.v[0].v3;
	vec3* y = &c->pose.v[1].v3;
	vec3* z = &c->pose.v[2].v3;
	vec3* chassisPos = &c->pose.v[3].v3;

	float cx = chassisPos->x;
	float cy = chassisPos->y;
	float cz = chassisPos->z;

	vec3 a, b;
	vecaddscale(&a, chassisPos, x, -0.5f);
	vecaddscale(&b, chassisPos, x, +0.5f);


	glColor4f(1,1,1,1);
	drawCircle(cx, cz, 0.1f);
	glBegin(GL_LINES);
	glVertex3f(a.x, a.z, 0.f);
	glVertex3f(b.x, b.z, 0.f);
	glEnd();

	suspensionDraw(c->suspension[0]);
	suspensionDraw(c->suspension[1]);
}

static vec3 getPointVel(Chassis* c, Suspension* s)
{
	vec3 v;
	vec3 angVel = {0.f, c->angVel, 0.f};
	vec3 offset;
	vec3mtx33mulvec3(&offset, &c->pose, &s->offset);
	veccross(&v, &angVel, &offset);
	vecadd(&v, &v, &c->vel);
	return v;
}



static float linetest(const vec3* pos, const vec3* dir, vec3* hitPos, vec3* normal)
{
	vecset(normal, 0.f, 0.f, 1.f);
	float dot = -dir->z;
	if (dot>0.f)
	{
		float distance = pos->z / dot;
		vecaddscale(hitPos, pos, dir, distance);
		return distance;
	}
	else
	{
		vecset(hitPos, 0.f, 0.f, -1000.f);
		return 1000.f;
	}
}

static float g_vel = 0.f;
static float g_angVel = 0.f;

static bool g_step = true;
static bool g_pauseEachFrame = false;


static void vehicleSubTick(Chassis* c, float dt)
{
	//g_step = true;
	if (!g_step) return;
	if (g_pauseEachFrame)
	{
		g_step = false;
	}

	const int numIterations = 8;
	float invIterations = 0.1f;
	
	s_chassis.vel.x = 0.f;	// HACK

	// This bit is done by the physics engine
	if(1)
	{
		//s_chassis.vel.z -= gravity*dt;
		{
			// Run Solver 
		}
		vec3* chassisPos = &s_chassis.pose.v[3].v3;
		vec3* x = &s_chassis.pose.v[0].v3;
		vec3* y = &s_chassis.pose.v[1].v3;
		vec3* z = &s_chassis.pose.v[2].v3;
		vecaddscale(chassisPos, chassisPos, &s_chassis.vel, dt);

		// Rotate
		float angVel = s_chassis.angVel;
		float c = cosf(dt*angVel);
		float s = sinf(dt*angVel);
		vec3 tmpX = *x;
		vec3 tmpY = *y;
		vec3 tmpZ = *z;
		vecscale(z, &tmpZ, c);
		vecaddscale(z, z, &tmpX, s);	// this doesn't seem the correct way around!
		vecscale(x, &tmpX, c);
		vecaddscale(x, x, &tmpZ, -s);

		//s_chassis.angVel = 0.f;
	}

	if (g_vel!=0.f)
	{
		s_chassis.vel.z = g_vel;
		//for (int i=0; i<numWheels; i++)
		//{
		//	Suspension* s = c->suspension[i];
		//	s->vel -= g_vel * s->axis.z;
		//}
	}
	if (g_angVel!=0.f)
	{
		s_chassis.angVel = g_angVel;
	}
		
	// Add our own gravity to help with collision calculations
	// Remove gravity after the process
	s_chassis.vel.z -= gravity*dt;

	vec3 extraVel={0.f, 0.f, 0.f};
	vec3 extraAngVel={0.f, 0.f, 0.f};

	for (int i=0; i<numWheels; i++)
	{
		Suspension* s = c->suspension[i];
		Wheel* w = s->wheel;
		
		// Update the suspension axis
		s->axis = c->pose.v[2].v3;

		s->compression = s->newCompression;
		
		// Calculate the world position and offset of the suspension point
		vec3mtx33mulvec3(&s->worldOffset, &c->pose, &s->offset);
		vec3mtx43mulvec3(&s->worldDefaultPos, &c->pose, &s->offset);

		// Linetest
		vec3 pos;
		vec3 dir;
		vecneg(&dir, &s->axis);
		vecaddscale(&pos, &s->worldDefaultPos, &s->axis, s->linetestLength - s->startingCompression);
		s->hitDistance = s->linetestLength - linetest(&pos, &dir, &s->hitPos, &s->hitNorm);		// Note hit distance measured from the lower end tip of the suspension, like compression
		s->hitPlaneD = -vecdot(&s->hitPos, &s->hitNorm);		// Plane equation, A*x + B*x + C*y + D = 0;

		vec3 pointVel = getPointVel(c, s);
		float chassisInlineSpeed = vecdot(&s->axis, &pointVel);

		// Convert suspension wheel speed to world space
		s->vel = s->vel + chassisInlineSpeed;

		// Add gravity
		//s->vel -= (gravity) * s->axis.z * dt;		// NB this is the world velocity projected onto the axis

		// Add spring force if contact!
		if ((s->compression + s->vel * dt) < s->hitDistance)
		{
			float damp = s->damp * c->mass;
			float spring = getRestingSpring(s->startingCompression, gravity, dt, c->mass, damp);
			float x = s->compression;
			float force = -getSpringForce(x, -vecdot(&s->axis, &pointVel), dt, c->mass, spring, damp) * (1.f/(float)(numWheels));
			vec3 impulse;
			vecscale(&impulse, &s->axis, force*dt);
			addImpulseAtOffset(&extraVel, &extraAngVel, 1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &impulse);
		}
	}

	vecadd(&c->vel, &c->vel, &extraVel);
	c->angVel += extraAngVel.y;

	// TODO: will need to optimise all the chassisInlineSpeed calculations!
	
	for (int repeat=0; repeat<numIterations; repeat++)
	{
		for (int i=0; i<numWheels; i++)
		{
			float collisionError = 0.f;
			Suspension* s = c->suspension[i];
			Wheel* w = s->wheel;

			// Side ground collision
			if (1)
			{
				float dot = vecdot(&s->axis, &s->hitNorm);
				if (dot<0.99f && dot>0.01f)
				{
					vec3 pointVel = getPointVel(c, s);
					float chassisInlineSpeed = vecdot(&s->hitNorm, &pointVel);
					vec3 pos;
					vecaddscale(&pos, &s->worldDefaultPos, &s->axis, s->maxCompression - s->startingCompression);
					vecaddscale(&pos, &pos, &s->hitNorm, dt*chassisInlineSpeed);
					float penetrationError = -vecdot(&pos, &s->hitNorm) - s->hitPlaneD + (1.f - dot*dot)*0.1f;		// Positive penetration means below the hit plane
					if (penetrationError>0.f)
					{
						vec3 angVel = {0.f, c->angVel, 0.f};
						float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &s->hitNorm);
						vec3 correctionImpulse;
						vecscale(&correctionImpulse, &s->hitNorm, invIterations*(penetrationError - 0.05*chassisInlineSpeed) / (denom * dt));
						addImpulseAtOffset(&c->vel, &angVel, 1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &correctionImpulse);
						c->angVel = angVel.y;
						printf("penetrationError = %f, impulse = %f %f %f, angVel = %f\n", penetrationError, XYZ(correctionImpulse), c->angVel);
						g_pauseEachFrame = false;
					}
				}
			}

			// Clamp the low limit of spring
			if (1)
			{
				// This affects only the wheel
				vec3 pointVel = getPointVel(c, s);
				float chassisInlineSpeed = vecdot(&s->axis, &pointVel);
				float compression2 = s->compression + (s->vel - chassisInlineSpeed) * dt;
				if (compression2<s->minCompression)
				{
					// This is where you would add any downward force to stop the vehicle from tipping! (but only if there is ground collision)
					float compressionError = (compression2 - s->minCompression) / dt;
					s->vel = s->vel - compressionError * invIterations;
				}
				if (compression2>s->maxCompression)
				{
					float compressionError = (compression2 - s->maxCompression) / dt;
					s->vel = s->vel - compressionError * invIterations;
				}
			}

			// Ground collision
			if (1)
			{
				collisionError = s->hitDistance - s->compression - s->vel*dt;
				if (collisionError < 0.f) collisionError = 0.f;
				collisionError = collisionError / dt;
				s->vel += collisionError * invIterations;
			}

			// Clamp the upper limit of spring
			if (1)
			{
				vec3 pointVel = getPointVel(c, s);
				float chassisInlineSpeed = vecdot(&s->axis, &pointVel);
				float compression2 = s->compression + (s->vel - chassisInlineSpeed) * dt;
				if (compression2>s->maxCompression)
				{
					float compressionError = (compression2 - s->maxCompression) / dt;
					if (0)
					{
						// This affects the chassis when the wheel collided
						float denom = 1.f/c->mass;
						float correctionImpulse = invIterations * compressionError / denom;
						vecaddscale(&c->vel, &c->vel, &s->axis, correctionImpulse / c->mass);
					}
					else
					{
						vec3 angVel = {0.f, c->angVel, 0.f};
						float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &s->axis);
						float correctionImpulse = invIterations * compressionError / denom;
						vec3 impulse;
						vecscale(&impulse, &s->axis, correctionImpulse);
						addImpulseAtOffset(&c->vel, &angVel, 1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &impulse);
						c->angVel = angVel.y;
					}
				}
			}
		}
		invIterations = 0.5f + 0.5f * invIterations;
	}

	g_angVel = 0.f;
	g_vel = 0.f;

	for (int i=0; i<numWheels; i++)
	{
		Suspension* s = c->suspension[i];
		Wheel* w = s->wheel;
		
		vec3 pointVel = getPointVel(c, s);
		float chassisInlineSpeed = vecdot(&s->axis, &pointVel);
		
		// Convert suspension velocity back to chassis space
		s->vel = s->vel - chassisInlineSpeed;
		s->newCompression = s->compression + s->vel * dt;

		s->vel = -5.0f;

		// Update the wheel pos
		vecaddscale(&w->pos, &s->worldDefaultPos, &s->axis, s->compression - s->startingCompression);
		//printf("compression = %f\n", s->compression);
	}


	// Subtract gravity, since the physics engine will do it as well
	//s_chassis.vel.z += gravity*dt;
}

static void vehicleSubTick(float dt)
{
	vehicleSubTick(&s_chassis, dt);
}

const float subDt = 0.02;//1.0/30.f;//0.01f;

static void vehicleTick(float dt)
{
	// Make this more like a physics engine with a "passive"
	// callback function

	// Fixed time step integration
	static float remaining = 0.f;
	remaining += dt;
	while (remaining>0.f)
	{
		remaining -= subDt;
		vehicleSubTick(subDt);
	}
}

static void wheelReset(Wheel* w, const vec3* pos)
{
	w->mass = 1.f;
	w->pos = *pos;
}

static void suspensionReset(Suspension* s, const vec3* offset)
{
	s->offset.x = s->wheel->pos.x;
	s->offset.y = s->wheel->pos.y;
	s->offset.z = offset->z;
	
	vec3* chassisPos = &s->chassis->pose.v[3].v3;

	/*
	   -----------
	         ||
	         ||
			---- ------------------------------------- 
	         ||                                      ^
	         ||                                      |
	         ||                                      |  maxCompression
	         ||                                      |
	        <==>  Wheel Starting pos --------------------------
	         ||                                      |  ^
	         ||                                      |  |
	         ||                                      |  |
			---- -------------------------------^-   |  | compression 
	         ||                  minCompression |    |  |
	         ||                                 |    |  |
	   ----------- --------------------------------------------------------------- unconstrained position of wheel

	   compresssion is measured upwards from the point where the spring would be relaxed.
	*/

	const float travelMax = 0.2f;
	const float travelMin = 0.2f;

	s->startingCompression = 0.2f;
	if (s->startingCompression < travelMin) s->startingCompression = travelMin;
	s->minCompression = s->startingCompression - travelMin;
	s->maxCompression = s->startingCompression + travelMax;

	s->spring = 0.f;
	s->damp = 3.1f;
	s->vel = 0.f;
	s->compression = s->startingCompression;
	s->newCompression = s->startingCompression;

	s->linetestLength = 1.f;
	s->hitDistance = 10000.f;

	if (s->linetestLength < s->maxCompression*1.2f) s->linetestLength = s->maxCompression*1.2f;

	vecset(&s->axis, 0.f, 0.f, 1.f);
	
	//=====================
	// Auto calc stiffness
	//=====================
 	// To pre calculate spring we need the subdt that we will work with! Pooh!
	// This calculation is made in isolation of other suspension
 	s->spring = gravity * (1.f + subDt*s->damp) / (s->compression - gravity * subDt * subDt);

	printf("s->spring = %f\n", s->spring);
	printf("s->damp= %f\n", s->damp);
	printf("s->startingCompression = %f\n", s->startingCompression);
	printf("s->compression = %f\n", s->compression);
}

static void vehicleReset(Chassis* chassis)
{
	Chassis* c = chassis;
	
	const float wheelOffset = 0.5f;

	//================
	// Reset chassis
	//================
	c->mass = 10.f;
	c->inertia = c->mass;
	matrixIdent(&c->pose);
	vecset(&c->pose.v[3].v3, 0.f, 0.f, wheelOffset);
	veczero(&c->vel);
	c->angVel = 0.f;
	
	//=================
	// Position Wheels
	//=================
	vec3 wpos0 = {+2.f, 0.f, 0.f};
	vec3 wpos1 = {-2.f, 0.f, 0.f};
	wheelReset(c->suspension[0]->wheel, &wpos0);
	wheelReset(c->suspension[1]->wheel, &wpos1);
	
	//==================
	// Reset suspension
	//==================
	vec3 offset0 = {+2.f, 0.f, -wheelOffset};
	vec3 offset1 = {-2.f, 0.f, -wheelOffset};
	
	suspensionReset(c->suspension[0], &offset0);
	suspensionReset(c->suspension[1], &offset1);

	// Set to new position
	matrixRotY(&c->pose, 1.5f);
	c->pose.v[3].v3.z = 4.0f;

	//=======================
	// Auto balance stiffness
	//=======================

	for (int i=0; i<numWheels; i++)
	{
		Suspension* s = c->suspension[i];
		s->spring *= 1.f/(float)numWheels;
	}

	vehicleTick(subDt);
}

static void vehicleReset()
{
	// Set up the pointers
	s_chassis.suspension[0] = &s_suspension[0];
	s_chassis.suspension[1] = &s_suspension[1];
	s_suspension[0].chassis = &s_chassis;
	s_suspension[1].chassis = &s_chassis;

	s_suspension[0].wheel = &s_wheel[0];
	s_suspension[1].wheel = &s_wheel[1];
	s_wheel[0].suspension = &s_suspension[0];
	s_wheel[1].suspension = &s_suspension[1];

	vehicleReset(&s_chassis);
}

static void vehicleInit()
{
	vehicleReset();
}



//======================================================================================================================================================
//======================================================================================================================================================

static void keyboard(unsigned char key, int x, int y)
{
	Chassis* c = &s_chassis;
	if (key=='r')
	{
		vehicleReset();
		g_pauseEachFrame = false;
	}
	if (key==27) exit(0);
	if (key==' ')
	{
		g_vel = 5.f;
	}
	if (key=='z')
	{
		g_angVel = 1.0f;
	}
	if (key=='s')
	{
		g_step = true;
		g_pauseEachFrame = true;
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


