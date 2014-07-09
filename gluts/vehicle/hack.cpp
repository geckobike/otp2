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
	return - mass * (spring * x + (a + damp) * vel) / (mass  + dt * (a + damp) );
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

	float compression;			// Current compression value

	// Calculated each frame
	vec3 axis;

	// Linetest
	vec3 hitPos;
	vec3 unfilteredHitNorm;
	vec3 hitNorm;
	float hitPlaneD;
	bool collision;

	// Calculated each frame
	vec3 worldOffset;
	vec3 worldDefaultPos;
};


struct Wheel
{
	Suspension* suspension;

    float mass;					// Doesn't really matter!
	float invMass;

	// Determined each frame
	vec3 pos;
	vec3 vel;
};

struct Chassis
{
	Suspension* suspension[2];

    float mass;
	float invMass;
	float inertia;

	mtx pose;
	vec3 vel;
	vec3 angVel;
};

const int numWheels = 2;
static Chassis s_chassis;
static Wheel s_wheel[numWheels];
static Suspension s_suspension[numWheels];
static float gravity = 10.f;
const float subDt = 0.02;//1.0/30.f;//0.01f;

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

	vec3 wheelVel;
	vecaddscale(&wheelVel, &w->pos, &w->vel, subDt);
	glColor4f(0,0,1,1);
	glBegin(GL_LINES);
	glVertex3f(wx, wz, 0.f);
	glVertex3f(wheelVel.x, wheelVel.z, 0.f);
	glEnd();


	// Linetest hit norm
	vec3 norm;
	vecaddscale(&norm, &w->pos, &s->unfilteredHitNorm, 1.f);
	glColor4f(1,0,0,1);
	glBegin(GL_LINES);
	glVertex3f(wx, wz, 0.f);
	glVertex3f(norm.x, norm.z, 0.f);
	glEnd();

	vecaddscale(&norm, &w->pos, &s->hitNorm, 1.f);
	glColor4f(1,1,0,1);
	glBegin(GL_LINES);
	glVertex3f(wx, wz, 0.f);
	glVertex3f(norm.x, norm.z, 0.f);
	glEnd();
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
	
	vec3 vel;
	vecaddscale(&vel, chassisPos, &c->vel, subDt);
	glColor4f(0,1,0,1);
	glBegin(GL_LINES);
	glVertex3f(cx, cz, 0.f);
	glVertex3f(vel.x, vel.z, 0.f);
	glEnd();


	for (int i=0; i<numWheels; i++)
	{
		suspensionDraw(c->suspension[i]);
	}
}

static vec3 getPointVel(Chassis* c, const vec3* offset)
{
	// offset is in world space

	vec3 v;
	veccross(&v, &c->angVel, offset);
	vecadd(&v, &v, &c->vel);
	return v;
}

static vec3 getPointVel(Chassis* c, Suspension* s)
{
	vec3 offset;
	vec3mtx33mulvec3(&offset, &c->pose, &s->offset);
	return getPointVel(c, &offset);
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

static void doCorrection(Chassis* c, const vec3* worldOffset, const vec3* axis, float requiredVelocityChange)
{
	float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, worldOffset, axis);
	float correction = requiredVelocityChange / denom;
	vec3 impulse;
	vecscale(&impulse, axis, correction);
	addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, worldOffset, &impulse);
}

static void doCorrection2(Chassis* c, const vec3* worldOffset, const vec3* axis, float requiredVelocityChange)
{
	const float t = 0.5f;
	const float u = 1.f - t;
	// Add half to the linear
    vec3 tmp;
	vecaddscale(&c->vel, &c->vel, axis, requiredVelocityChange*t);

	// Add the other half to the angular
	vec3 cross;
	veccross(&cross, worldOffset, axis);
	veccross(&cross, &cross, worldOffset);
	float angularResponse = vecdot(axis, &cross);
	float angVelChange = requiredVelocityChange / angularResponse;

	vec3 impulse;
	vecscale(&impulse, axis, angVelChange*u);
	veccross(&tmp, worldOffset, &impulse);
	vecadd(&c->angVel, &c->angVel, &tmp);
}

static void doCorrection3(Chassis* c, const vec3* worldOffset, const vec3* axis, float requiredVelocityChange, float linearRatio)
{
	const float u = 1.f - linearRatio;
	// Add half to the linear
    vec3 tmp;
	vecaddscale(&c->vel, &c->vel, axis, requiredVelocityChange*linearRatio);

	// Add the other half to the angular
	vec3 cross;
	veccross(&cross, worldOffset, axis);
	veccross(&cross, &cross, worldOffset);
	float angularResponse = vecdot(axis, &cross);
	float angVelChange = requiredVelocityChange / angularResponse;

	vec3 impulse;
	vecscale(&impulse, axis, angVelChange*u);
	veccross(&tmp, worldOffset, &impulse);
	vecadd(&c->angVel, &c->angVel, &tmp);
}

static int g_step = 1;		// 1 = pause each frame, 2 = on

static const int numIterations = 5;
	

static void vehicleSubTick(Chassis* c, float dt)
{
	if (g_step==0) return;
	if (g_step&1) g_step = 0;

	vec3* chassisPos = &s_chassis.pose.v[3].v3;
	vec3* x = &s_chassis.pose.v[0].v3;
	vec3* y = &s_chassis.pose.v[1].v3;
	vec3* z = &s_chassis.pose.v[2].v3;

	// This bit is done by the physics engine
	if(1)
	{
		vecaddscale(chassisPos, chassisPos, &s_chassis.vel, dt);

		// Rotate
		float angVel = s_chassis.angVel.y;
		float c = cosf(dt*angVel);
		float s = sinf(dt*angVel);
		vec3 tmpX = *x;
		vec3 tmpY = *y;
		vec3 tmpZ = *z;
		vecscale(z, &tmpZ, c);
		vecaddscale(z, z, &tmpX, s);	// this doesn't seem the correct way around!
		vecscale(x, &tmpX, c);
		vecaddscale(x, x, &tmpZ, -s);
	}

	// Damp
	vecscale(&c->vel, &c->vel, expf(-dt*1.f));
	vecscale(&c->angVel, &c->angVel, expf(-dt*1.f));

	if (fabsf(c->angVel.x)<0.01f) c->angVel.x = 0.f;
	if (fabsf(c->angVel.y)<0.01f) c->angVel.y = 0.f;
	if (fabsf(c->angVel.z)<0.01f) c->angVel.z = 0.f;

	bool addGravity = true;

	vec3 springImpulse[numWheels] = {0};
	float springVelocity[numWheels];
	float springMaxImpulse[numWheels]; 
	float impulseApplied[numWheels]={0};

	// Prepare
	for (int i=0; i<numWheels; i++)
	{
		Suspension* s = c->suspension[i];
		Wheel* w = s->wheel;

		s->collision = false;
		
		// Update the suspension axis
		s->axis = c->pose.v[2].v3;

		// Calculate the world position and offset of the suspension point
		vec3mtx33mulvec3(&s->worldOffset, &c->pose, &s->offset);
		vec3mtx43mulvec3(&s->worldDefaultPos, &c->pose, &s->offset);

		vecaddscale(&w->pos, &s->worldDefaultPos, &s->axis, s->compression - s->startingCompression);

		// Linetest
		vec3 pos;
		vec3 dir;
		vecneg(&dir, &s->axis);
		vecaddscale(&pos, &s->worldDefaultPos, &s->axis, s->linetestLength - s->startingCompression);
		s->linetestLength - linetest(&pos, &dir, &s->hitPos, &s->hitNorm);		// Note hit distance measured from the lower end tip of the suspension, like compression
		s->unfilteredHitNorm = s->hitNorm;
	
		if (0)
		{
			// Create bastard normal!
			vecaddscale(&s->hitNorm, &s->hitNorm, x, 3.f);
			vecnormalise(&s->hitNorm);
			s->hitNorm = *x;
			s->unfilteredHitNorm = s->hitNorm;

			// Filter
			vecadd(&s->hitNorm, &s->hitNorm, &s->axis);
			vecnormalise(&s->hitNorm);
		}

		s->hitPlaneD = -vecdot(&s->hitPos, &s->hitNorm);		// Plane equation, A*x + B*x + C*y + D = 0;

		
		vec3 pointVel = getPointVel(c, &s->worldOffset);

		// Convert suspension wheel speed to world space
		vecadd(&w->vel, &w->vel, &pointVel);
		
		if(1)
		{
			float damp = s->damp * w->mass;
			float spring = getRestingSpring(s->startingCompression, gravity, dt, c->mass, damp);
			float x = s->compression;
			float v = vecdot(&s->axis, &w->vel) - vecdot(&s->axis, &pointVel);
			float force = -getSpringForce(x, v, dt, c->mass, spring, damp) * (1.f/(float)(numWheels));

			vecscale(&springImpulse[i], &s->axis, force*dt/(float)numIterations);

			force = -force;
			springMaxImpulse[i] = -force*dt;
			springVelocity[i] = v - springMaxImpulse[i]*(1.f/w->mass + 1.f/c->mass);

			springMaxImpulse[i] = fabsf(springMaxImpulse[i]);
			impulseApplied[i] = 0.f;
		}

		// Add gravity after the spring (this is important!)
		if (addGravity) w->vel.z -= gravity * dt;
	}

	// Add our own gravity
	if (addGravity) s_chassis.vel.z -= gravity*dt;

	if (0)
	{
		// Test downward hack force 
		static float time = 0.f;
		time += dt;
		//s_chassis.vel.z += sinf(5.f*time) * 1.1f*gravity*dt;
		s_chassis.vel.z += 1.f*gravity*dt;
	}

	//=============
	//   VERBOSE
	//=============
	#define verbose false
	#define dump if (verbose) printf

	printf("==========================================\n");
	printf("START ITERATION\n");
	printf("==========================================\n");

	float solverERP = numIterations>1 ? 0.1f : 1.f;
	float changeSolverERP = numIterations>1 ? (1.f - solverERP)/(float)(numIterations-1) : 0.f;

	for (int repeat=0; repeat<numIterations; repeat++)
	{
		printf(" == Start Iter == \n");

			for (int i=0; i<numWheels; i++)
			{
				if (i==1) continue;
				Suspension* s = c->suspension[i];
				Wheel* w = s->wheel;

				if (1)	// Add Spring
				{
					addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &springImpulse[i]);
					vecaddscale(&w->vel, &w->vel, &springImpulse[i], -1.f/w->mass);
				}
				else
				{
					vec3 pointVel = getPointVel(c, &s->worldOffset);
					float v = vecdot(&s->axis, &w->vel) - vecdot(&s->axis, &pointVel);
					float error = v - springVelocity[i];
					//float denom = c->invMass + w->invMass;
					float denom = computeDenominator(c->invMass, 1.f/c->inertia, &s->worldOffset, &s->axis) + w->invMass;
					float impulse = solverERP * error / denom;

					printf("[%d] error = %f denom = %f, impulse = %f applied = %f, springMaxImpulse = %f\n", i, error, denom, impulse, impulseApplied[i], springMaxImpulse[i]);


					float total = impulseApplied[i] + impulse;
					if (total>+springMaxImpulse[i]) impulse -= total - springMaxImpulse[i];
					if (total<-springMaxImpulse[i]) impulse -= total + springMaxImpulse[i];
					impulseApplied[i] += impulse;

					printf("[%d] applying impulse = %f\n", i, impulse);

					//vecscale(&error, &error, -solverERP/denom);
					//addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &offset, &error);
					//vecaddscale(&w->vel, &w->vel, &error, -solverERP/w->mass);

					vecaddscale(&c->vel, &c->vel, &s->axis, +impulse/c->mass);
					vecaddscale(&w->vel, &w->vel, &s->axis, -impulse/w->mass);
				}
			}
		

		dump("after springs: c->angVel = %f %f %f\n", XYZ(c->angVel));

		for (int i=0; i<numWheels; i++)
		{
			Suspension* s = c->suspension[i];
			Wheel* w = s->wheel;
		
			const bool groundCollision = true;
			const bool axisError = true;
			const bool wheelLimits = true;
			const bool friction = true;

			float penetration = -0.001f;

			if (groundCollision)	// Ground Collision
			{
				penetration = -vecdot(&w->pos, &s->hitNorm) - s->hitPlaneD - vecdot(&s->hitNorm, &w->vel) * dt;
				penetration = penetration / dt;
				if (penetration > 0.f)
				{
					vecaddscale(&w->vel, &w->vel, &s->hitNorm, penetration*solverERP);

					// Friction
					if (friction)
					{
						vec3 lateral;
						vec3 lateralVel;
						vecaddscale(&lateralVel, &w->vel, &s->hitNorm, -vecdot(&s->hitNorm, &w->vel));

						//float mu = 1.0f;

						float sizeSq = vecsizesq(&lateralVel);
						if (sizeSq>0.001f)
						{
							vecscale(&lateral, &lateralVel, 1.f/sqrtf(sizeSq));
							float v = vecdot(&lateralVel, &lateral);
							//if (fabsf(v)<0.01f)
							//{
							//	mu = 1.f;
							//}

							//float error = -mu*v;
							//vecscale(&lateralVel, &lateralVel, 1.f/v);

							vecsubscale(&w->vel, &w->vel, &lateralVel, solverERP);

							//====================================
							// THIS DOESN'T LOOK RIGHT?
							//====================================
							vecsubscale(&c->vel, &c->vel, &lateralVel, solverERP);

							// float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &s->offset, &lateral);
							// vec3 impulse;
							// vecscale(&impulse, &lateral, -v/denom);
							// addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &s->offset, &impulse);

							//{
							//	vec3 pointVel = getPointVel(c, s);
							//	vec3 lateralVel;
							//	vecaddscale(&lateralVel, &pointVel, &s->hitNorm, -vecdot(&s->hitNorm, &pointVel));
							//	printf("adding vel [%d] , lateralVel (after) = %f %f %f\n", i, XYZ(lateralVel));
							//}
						}

					}	
					dump("gound collision errorV = %f, vel of wheel after = %f\n", penetration, vecdot(&w->vel, &s->hitNorm));
				}
			}

			if (axisError)	// Axis Error
			{
				vec3 pointvel = getPointVel(c, s);
				vec3 error;
				vecsub(&error, &pointvel, &w->vel);
				vecaddscale(&error, &error, &s->axis, -vecdot(&error, &s->axis));

				vec3 norm;
				if (vecsizesq(&error)>0.001f)
				{
					dump("axis error %f\n", vecsize(&error));
					vecnormalise(&norm, &error);

					vec3 offset;
					vecsub(&offset, &w->pos, chassisPos);
					float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &offset, &norm) + 1.f/w->mass;
					vecscale(&error, &error, -solverERP/denom);
					addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &offset, &error);
					vecaddscale(&w->vel, &w->vel, &error, -solverERP/w->mass);
				}
				dump("axis error vel of wheel after = %f, inline = %f\n", vecdot(&w->vel, &s->hitNorm), vecdot(&w->vel, &s->axis));
			}

			if (wheelLimits)
			{
				vec3 pointvel = getPointVel(c, &s->worldOffset);
				vec3 vel;
				vecsub(&vel, &w->vel, &pointvel);
				float speed = vecdot(&vel, &s->axis);
				float compression2 = s->compression + speed * dt;

				// Lower limit error
				if (compression2 < s->minCompression)
				{
					if (penetration>0.f)
					{
						#if 0
						//===========================
						// LOOKS LIKE THIS IS BROKEN!
						//===========================
						// Probably wont need this
						float error = (s->minCompression - compression2)/dt;
						float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &s->axis);
						vec3 impulse;
						vecscale(&impulse, &s->axis, -solverERP*error/denom);
						addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &impulse);
						#endif
					}
					{
						dump("lower limit hit\n");
						float error = (s->minCompression - compression2)/dt;
						// float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &s->axis) + 1.f/w->mass;
						// vec3 impulse;
						// vecscale(&impulse, &s->axis, -solverERP*error/denom);
						// addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &s->worldOffset, &impulse);

						if (1)
						{
							float denom = 1.f/c->mass + 1.f/w->mass;
							vec3 impulse;
							vecscale(&impulse, &s->axis, -error/denom);
							vecaddscale(&c->vel, &c->vel, &impulse, +solverERP/c->mass);
							vecaddscale(&w->vel, &w->vel, &impulse, -solverERP/w->mass);
							dump("c->angVel = %f %f %f\n", XYZ(c->angVel));
						}

					}
				}
				if (compression2 > s->maxCompression)
				{
					float error = (s->maxCompression - compression2)/dt;
					vecaddscale(&w->vel, &w->vel, &s->axis, solverERP*error);

					dump("upper limit hit: axis error vel of wheel after = %f\n", vecdot(&w->vel, &s->hitNorm));
				}
			}

			if (groundCollision)	// Upper limit (for chassis)
			{
				vec3 pos, offset;
				vecaddscale(&pos, &s->worldDefaultPos, &s->axis, s->maxCompression-s->startingCompression);
				vecsub(&offset, &pos, chassisPos);
				vec3 pointvel = getPointVel(c, &offset);
				float penetration = -vecdot(&pos, &s->hitNorm) - s->hitPlaneD - vecdot(&s->hitNorm, &pointvel) * dt;
				penetration = penetration / dt;
				if (penetration > 0.f)
				{
					float denom = 1.f/c->mass;
					
					dump("upper limit hit chassis\n");

					if (0)
					{
						vec3 impulse;
						vecscale(&impulse, &s->hitNorm, penetration);
						vecadd(&c->vel, &c->vel, &impulse);
						//vecadd(&w->vel, &w->vel, &impulse);
					}
					if (1)
					{
						// Should probably fix axis error here directly!
						float denom = computeDenominator(1.f/c->mass, 1.f/c->inertia, &s->offset, &s->hitNorm);
						vec3 impulse;
						vecscale(&impulse, &s->hitNorm, solverERP*penetration/denom);
						addImpulseAtOffset(&c->vel, &c->angVel, 1.f/c->mass, 1.f/c->inertia, &offset, &impulse);
					}

					if (0)
					{
						// Not Working!
						doCorrection2(c, &offset, &s->hitNorm, penetration);
					}
					dump("upper chassis limit hit, vel of wheel after = %f\n", vecdot(&w->vel, &s->hitNorm));
				}
			}

		}
		solverERP += changeSolverERP;
	}

	for (int i=0; i<numWheels; i++)
	{
		Suspension* s = c->suspension[i];
		Wheel* w = s->wheel;
		vec3 pointVel = getPointVel(c, s);
		
		// Convert suspension wheel speed back to car space
		vecsub(&w->vel, &w->vel, &pointVel);

		s->compression = s->compression + vecdot(&s->axis, &w->vel) * dt;

		// Damp wheel speed
		//vecscale(&w->vel, &w->vel, expf(-dt));
	}
}

static void vehicleSubTick(float dt)
{
	vehicleSubTick(&s_chassis, dt);
}


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

static void suspensionReset(Suspension* s, const vec3* offset)
{
	s->offset = *offset;
	
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
	s->damp = 10.0f;
	s->compression = s->startingCompression;

	s->linetestLength = 10.f;

	veczero(&s->wheel->vel);

	if (s->linetestLength < s->maxCompression*1.2f) s->linetestLength = s->maxCompression*1.2f;

	vecset(&s->axis, 0.f, 0.f, 1.f);

	// Limits!
	float smallestMass = s->chassis->mass * 0.05f;
	if (s->wheel->mass < smallestMass) s->wheel->mass = smallestMass;
	if (s->damp>10.f) s->damp = 10.f;
	s->wheel->invMass = 1.f/s->wheel->mass;
}

static void vehicleReset(Chassis* chassis)
{
	Chassis* c = chassis;
	
	const float wheelOffset = 0.3f;

	//================
	// Reset chassis
	//================
	c->mass = 10.f;
	c->inertia = c->mass * 0.2;
	c->invMass = 1.f/c->mass;
	matrixIdent(&c->pose);
	vecset(&c->pose.v[3].v3, 0.f, 0.f, wheelOffset);
	veczero(&c->vel);
	veczero(&c->angVel);
	
	//==================
	// Reset suspension
	//==================
	vec3 offset0 = {+2.f, 0.f, -wheelOffset};
	vec3 offset1 = {-2.f, 0.f, -wheelOffset};
	
	suspensionReset(c->suspension[0], &offset0);
	suspensionReset(c->suspension[1], &offset1);

	// Set to new position
	if (0)
	{
		matrixRotY(&c->pose, 1.0f);
	}

	if(0)
	{
		vecset(&c->pose.v[0].v3, -1.f, 0.f, 0.f);
		vecset(&c->pose.v[2].v3, 0.f, 0.f, -1.f);
	}

	vec3 offset;
	vec3mtx43mulvec3(&offset, &c->pose, &c->suspension[0]->offset);
	c->pose.v[3].v3.z -= offset.z;
	c->pose.v[3].v3.z += 1.f;
	
	// Set the correct wheel position
	for (int i=0; i<numWheels; i++)
	{
		Suspension* s = c->suspension[i];
		Wheel* w = s->wheel;
		s->axis = c->pose.v[2].v3;
		vec3mtx33mulvec3(&s->worldOffset, &c->pose, &s->offset);
		vec3mtx43mulvec3(&s->worldDefaultPos, &c->pose, &s->offset);
		vecaddscale(&w->pos, &s->worldDefaultPos, &s->axis, s->compression - s->startingCompression);
	}

}

static void vehicleReset()
{
	// Set up the pointers
	for (int i=0; i<numWheels; i++)
	{
		s_chassis.suspension[i] = &s_suspension[i];
		s_suspension[i].chassis = &s_chassis;
		s_suspension[i].wheel = &s_wheel[i];
		s_wheel[i].suspension = &s_suspension[i];
	}

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
		g_step = 2;
	}
	if (key==27) exit(0);
	if (key=='s')
	{
		g_step = g_step^1;
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
    
    glTranslatef(0.f, -1.f, -5.5f);

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


