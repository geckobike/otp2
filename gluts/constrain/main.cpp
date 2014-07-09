#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ENABLE_TIMER
#include "engine/engine.h"
#include "matrixmath.h"

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


static inline float estimateAngle(float dot)
{
	return (dot<1.f) ? acosf(dot) : 0.f;
	return (dot<1.f) ? sqrtf((2.f - 2.f*dot)) : 0.f;
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

/*
=================================================================
   Constraints
=================================================================
*/

struct Matrix2x2
{
	float a[2][2];
};

struct Body
{
	quaternion rot;
	vec3 pos;
	vec3 v;
	vec3 w;
	float mass;
	float inertia;
	float invMass;
	float invInertia;
	float damping;
	float angDamping;
};

void BodyInit(Body* body, float mass=0.f, float inertia=0.f);
void BodySetMass(Body* body, float mass=0.f, float inertia=0.f);
void BodyUpdate(Body* body, float dt);
void BodyApplyImpulse(Body* b, const vec3* impulse, const vec3* offset);

struct ConstraintAngSpring
{
	vec3 axis;
	float targetSpeed;
	float denominator;
	float denominator2;
};

struct Constraint
{
	vec3 offset[2];  // relative offset of the point constraint, in the space of each body
	Body* body[2];

	float spring;
	float damp;

	vec3 woffset[2]; // world offset
	vec3 werror; // world positional error
	Matrix2x2 response2d;
	mtx response3d;
};

// A helper class for drawing purposes
struct CapsuleBody
{
	Body body;
	float radius;
	float length;
};

static void dlDrawCross(float size, vec4* colour)
{
	vec3 p0 = {+size,0.f,0.f};
	vec3 p1 = {-size,0.f,0.f};
	vec3 p2 = {0.f,+size,0.f};
	vec3 p3 = {0.f,-size,0.f};
	vec3 p4 = {0.f,0.f,+size};
	vec3 p5 = {0.f,0.f,-size};
	dlDrawLine(&p0, &p1, colour);
	dlDrawLine(&p2, &p3, colour);
	dlDrawLine(&p4, &p5, colour);
}

static void BodyDraw(Body* body, vec4* colour)
{
	mtx m;
	quaternionToMatrix(&m, &body->rot);
	veccpy((vec3*)&m.v[3], &body->pos);
	glPushMatrix();
	dlMultMatrix(&m);
	dlDrawCross(0.02f, colour);
	glPopMatrix();
}


static void CapsuleBodyDraw(CapsuleBody* cap, vec4* colour)
{
	mtx m;
	quaternionToMatrix(&m, &cap->body.rot);
	veccpy((vec3*)&m.v[3], &cap->body.pos);
	glPushMatrix();
	dlMultMatrix(&m);
	{
	vec3 pos0 = {+cap->radius, -0.5f*cap->length, 0.f};
	vec3 pos1 = {+cap->radius, +0.5f*cap->length, 0.f};
	vec3 pos2 = {-cap->radius, -0.5f*cap->length, 0.f};
	vec3 pos3 = {-cap->radius, +0.5f*cap->length, 0.f};
	vec3 pos4 = {0.f, -0.5f*cap->length - cap->radius, 0.f};
	vec3 pos5 = {0.f, +0.5f*cap->length + cap->radius, 0.f};
	dlDrawLine(&pos0, &pos1, colour);
	dlDrawLine(&pos1, &pos5, colour);
	dlDrawLine(&pos5, &pos3, colour);
	dlDrawLine(&pos3, &pos2, colour);
	dlDrawLine(&pos2, &pos4, colour);
	dlDrawLine(&pos4, &pos0, colour);
	}
	{
	vec3 pos0 = {0.f, -0.5f*cap->length, +cap->radius};
	vec3 pos1 = {0.f, +0.5f*cap->length, +cap->radius};
	vec3 pos2 = {0.f, -0.5f*cap->length, -cap->radius};
	vec3 pos3 = {0.f, +0.5f*cap->length, -cap->radius};
	vec3 pos4 = {0.f, -0.5f*cap->length - cap->radius, 0.f};
	vec3 pos5 = {0.f, +0.5f*cap->length + cap->radius, 0.f};
	dlDrawLine(&pos0, &pos1, colour);
	dlDrawLine(&pos1, &pos5, colour);
	dlDrawLine(&pos5, &pos3, colour);
	dlDrawLine(&pos3, &pos2, colour);
	dlDrawLine(&pos2, &pos4, colour);
	dlDrawLine(&pos4, &pos0, colour);
	}
	dlDrawCross(0.02f, colour);
	glPopMatrix();
}

void BodyInit(Body* body, float mass, float inertia)
{
	veczero(&body->pos);
	quaternionIdent(&body->rot);
	veczero(&body->v);
	veczero(&body->w);
	BodySetMass(body, mass, inertia);
	body->damping = 1.1f;
	body->angDamping = 1.1f;
}

void BodySetMass(Body* body, float mass, float inertia)
{
	body->mass = mass;
	body->inertia = inertia;
	body->invMass = mass>0.f ? 1.f/mass : 0.f;
	body->invInertia = inertia>0.f ? 1.f/inertia : 0.f;
}

void UtilGetNewRotation(quaternion* qNew, const quaternion* q, const vec3* w, float dt)
{
   	if (1)
	{
		// Seem to get much better results with this _very_ approximat integration method
		quaternion dq;
		vecscale(&dq.v, w, 0.5f*dt);
		dq.w = 1.f - 0.5f*vecsizesq(&dq.v);
		// float size = 0.5f*dt*vecnormaliseSize(&q.v, angvel);
		// vecscale(&q.v, &q.v, sinf(size));
		// q.w = cosf(size);
		quaternionMul(qNew, &dq, q);
		quaternionNormalise(qNew);
	}
}


void BodyUpdate(Body* body, float dt)
{
	vecaddscale(&body->pos, &body->pos, &body->v, dt); 

	UtilGetNewRotation(&body->rot, &body->rot, &body->w, dt);
#if 0
	if (1)
	{
		// Seem to get much better results with this _very_ approximat integration method
		quaternion q;
		vecscale(&q.v, &body->w, 0.5f*dt);
		q.w = 1.f - 0.5f*vecsizesq(&q.v);
		// float size = 0.5f*dt*vecnormaliseSize(&q.v, angvel);
		// vecscale(&q.v, &q.v, sinf(size));
		// q.w = cosf(size);
		quaternionMul(&body->rot, &q, &body->rot);
		quaternionNormalise(&body->rot);
	}
	else
	{
		/*
		 * From thesis: Brian Vincent Mirtich
		 *
		 *        | qs |         |  -x  -y  -z  | 
		 * dq/dt =| qx | =  1/2  |  +s  -z  +y  |  w
		 *        | qy |         |  +z  +s  -x  |  
		 *        | qz |         |  -y  +x  +s  |
		*/
		vec3 w = body->w;
		quaternion q = body->rot;
		quaternion dq;
		dq.w = (-q.x*w.x - q.y*w.y - q.z*w.z) * 0.5f * dt;
		dq.x = (+q.w*w.x - q.z*w.y + q.y*w.z) * 0.5f * dt;
		dq.y = (+q.z*w.x + q.w*w.y - q.x*w.z) * 0.5f * dt;
		dq.z = (-q.y*w.x + q.x*w.y + q.w*w.z) * 0.5f * dt;
		quaternionAdd(&body->rot, &body->rot, &dq);
		quaternionNormalise(&body->rot);

		// And is the same as:
		// q.w -= (w*q.v)*dt*0.5f;
		// q.v += ((w^q.v)+w*q.w)*(dt*0.5f);
	}
#endif
	float damping = MAX(0.f, 1.f - dt*body->damping);
	float angDamping = MAX(0.f, 1.f - dt*body->angDamping);
	vecscale(&body->v, &body->v, damping);
	vecscale(&body->w, &body->w, angDamping);
}

void BodyApplyImpulse(Body* b, const vec3* impulse, const vec3* offset)
{
	vec3 torque;
	vecaddscale(&b->v, &b->v, impulse, b->invMass);
	veccross(&torque, offset, impulse);
	vecaddscale(&b->w, &b->w, &torque, b->invInertia);	// Treat inertia as a scalar (for now)
}

void BodyApplyAngImpulse(Body* b, const vec3* angImpulse)
{
	vecaddscale(&b->w, &b->w, angImpulse, b->invInertia);	// Treat inertia as a scalar (for now)
}

void BodyAddToResponseMatrix(Matrix2x2* response, Body* b, float rx, float ry)
{
	// rx, ry, is the world offset to the point where the impulse would be applied
	// in the literature this is JMJ^T
	response->a[0][0] += b->invMass + b->invInertia * SQR(ry);
	response->a[0][1] += -b->invInertia * rx*ry;
	response->a[1][0] += -b->invInertia * rx*ry;
	response->a[1][1] += b->invMass + b->invInertia * SQR(rx);
}

void BodyAddToResponseMatrix(mtx* response, Body* b, const vec3* woffset)
{
	// Linear response to an unity impulse
	response->f[0][0] += b->invMass; 
	response->f[1][1] += b->invMass; 
	response->f[2][2] += b->invMass; 

	// Torque response to an unity impulse at woffset
	vec3 xAxis = {1.f,0.f,0.f};
	vec3 yAxis = {0.f,1.f,0.f};
	vec3 zAxis = {0.f,0.f,1.f};
	vec3 torque, v0, v1, v2, w;
	veccross(&torque, woffset, &xAxis);
	vecscale(&w, &torque, b->invInertia);	// Treat inertia as a scalar (for now)
	veccross(&v0, &w, woffset);
	veccross(&torque, woffset, &yAxis);
	vecscale(&w, &torque, b->invInertia);	// Treat inertia as a scalar (for now)
	veccross(&v1, &w, woffset);
	veccross(&torque, woffset, &zAxis);
	vecscale(&w, &torque, b->invInertia);	// Treat inertia as a scalar (for now)
	veccross(&v2, &w, woffset);

	response->f[0][0] += v0.x;
	response->f[0][1] += v0.y;
	response->f[0][2] += v0.z;
	response->f[1][0] += v1.x;
	response->f[1][1] += v1.y;
	response->f[1][2] += v1.z;
	response->f[2][0] += v2.x;
	response->f[2][1] += v2.y;
	response->f[2][2] += v2.z;
}

void BodyGetVelocityAtOffset(vec3* pointvel, Body* body, const vec3* woffset)
{
	veccross(pointvel, &body->w, woffset);
	vecadd(pointvel, pointvel, &body->v);
}


void Matrix2x2Invert(Matrix2x2* out, Matrix2x2* in)
{
	const float a11 = in->a[0][0];
	const float a12 = in->a[0][1];
	const float a21 = in->a[1][0];
	const float a22 = in->a[1][1];

	float det = a11*a22 - a12*a21;
	if (fabsf(det)<1e-6f)
	{
		out->a[0][0] = 0.f;
		out->a[0][1] = 0.f;
		out->a[1][0] = 0.f;
		out->a[1][1] = 0.f;
	}
	else
	{
		out->a[0][0] = +a22/det;
		out->a[0][1] = -a12/det;
		out->a[1][0] = -a21/det;
		out->a[1][1] = +a11/det;
	}
}

void Matrix2x2Mul(float out[2], const Matrix2x2* m, float in[2])
{
	float a = m->a[0][0] * in[0] + m->a[0][1] * in[1];
	float b = m->a[1][0] * in[0] + m->a[1][1] * in[1];
	out[0] = a;
	out[1] = b;
}


void GetErrorVelocity(vec3* errorVel, Constraint* c)
{
	// In the literature this is dC/dt or J*v
	float erp = 0.8f; // position error reduction parameter
	vec3 v0, v1;
	BodyGetVelocityAtOffset(&v0, c->body[0], &c->woffset[0]);
	BodyGetVelocityAtOffset(&v1, c->body[1], &c->woffset[1]);
	vecsub(errorVel, &v1, &v0);
	vecaddscale(errorVel, errorVel, &c->werror, erp);
}

void ConstraintDraw(Constraint* c)
{
	mtx bodyMat;
	mtx m;
    matrixIdent(&m);

	for (int i=0; i<2; i++)
	{
		mtx bodyMat;
		quaternionToMatrix(&bodyMat, &c->body[i]->rot);
		veccpy((vec3*)&bodyMat.v[3], &c->body[i]->pos);
		glPushMatrix();
		dlMultMatrix(&bodyMat);
		veccpy((vec3*)&m.v[3], &c->offset[i]);
		glPushMatrix();
		dlMultMatrix(&m);
		dlDrawCircle(0.f, 0.f, 0.01f, false);
		glPopMatrix();
		glPopMatrix();
	}
}

float ConstraintCalcRotationalDenominator(Body* b0, Body* b1, const vec3* axis)
{
	vec3 response0, response1;
	vecscale(&response0, axis, b0->invInertia);
	vecscale(&response1, axis, b1->invInertia);
	return vecdot(axis, &response0) + vecdot(axis, &response1);
}

void ConstraintPrepare(Constraint* c, float dt)
{
	vec3 p0, p1;

	// Calculate the world offsets
	quaternionRotateVector(&c->woffset[0], &c->body[0]->rot, &c->offset[0]);
	quaternionRotateVector(&c->woffset[1], &c->body[1]->rot, &c->offset[1]);
	
	// Calculate the world positional error
	vecadd(&p0, &c->body[0]->pos, &c->woffset[0]);
	vecadd(&p1, &c->body[1]->pos, &c->woffset[1]);
	vecsub(&c->werror, &p1, &p0);
	vecscale(&c->werror, &c->werror, 1.f/dt);
	
	// Calculate the response matrix
	memset(&c->response2d, 0, sizeof(c->response2d));
	BodyAddToResponseMatrix(&c->response2d, c->body[0], c->woffset[0].x, c->woffset[0].y);
	BodyAddToResponseMatrix(&c->response2d, c->body[1], c->woffset[1].x, c->woffset[1].y);
	Matrix2x2Invert(&c->response2d, &c->response2d);
	
	memset(&c->response3d, 0, sizeof(c->response3d));
	BodyAddToResponseMatrix(&c->response3d, c->body[0], &c->woffset[0]);
	BodyAddToResponseMatrix(&c->response3d, c->body[1], &c->woffset[1]);
	matrixInvert(&c->response3d, &c->response3d);
}

// What is JMJ^T ? Its just the k matrix
// i.e. if you poke the constraint with an impulse of (1,0,0,...) you get this error out
// the exercise is to eliminate the error


void ConstraintSolve(vec3* solution, Constraint* c)
{
	vec3 errorVel;
	GetErrorVelocity(&errorVel, c);
	//Matrix2x2Mul((float*)solution, &c->response2d, (float*)&errorVel);
	vec3mtx33mulvec3(solution, &c->response3d, &errorVel);
}

void ConstraintAngSpringSolve(vec3* solution, Constraint* c, ConstraintAngSpring* angSpring)
{
	vec3 relangvel;
	vecsub(&relangvel, &c->body[1]->w, &c->body[0]->w);
	float relSpeed = vecdot(&relangvel, &angSpring->axis);
	float impulseSize = (relSpeed-angSpring->targetSpeed)/angSpring->denominator;
	vecscale(solution, &angSpring->axis, impulseSize);
	
	// This is wrong!
	//const vec3* yaxis = &c->body[0]->rot.v[1].v3;
	//relSpeed = vecdot(&relangvel, yaxis);
	//impulseSize = relSpeed/angSpring->denominator2;
	//vecaddscale(solution, solution, yaxis, 0.1f*impulseSize);
}

void ConstraintAngLimitSolve(vec3* solution, Constraint* c, float maxAngle, float dt)
{
	quaternion q0, q1;
	UtilGetNewRotation(&q0, &c->body[0]->rot, &c->body[0]->w, dt);
	UtilGetNewRotation(&q1, &c->body[1]->rot, &c->body[1]->w, dt);
    
	vec3 axesA; quaternionGetYBasis(&axesA, &q0);
    vec3 axesB; quaternionGetYBasis(&axesB, &q1);

	float angle = estimateAngle(vecdot(&axesA, &axesB));
	if (angle>(maxAngle+DEG2RAD(0.1f)))
	{
		vec3 axis;
		veccross(&axis, &axesA, &axesB);    // NB: positive rotation about axis rotates A towards B
		vecnormalise(&axis);
		float denominator = ConstraintCalcRotationalDenominator(c->body[0], c->body[1], &axis);
		vecscale(solution, &axis, MIN(DEG2RAD(90.f/dt), angle - maxAngle) / (dt * denominator));
	}
	else
	{
		veczero(solution);
	}
}

void ConstraintAngLimitSolve2(vec3* solution, Constraint* c, float minAngle, float dt)
{
	quaternion q0, q1;
	UtilGetNewRotation(&q0, &c->body[0]->rot, &c->body[0]->w, dt);
	UtilGetNewRotation(&q1, &c->body[1]->rot, &c->body[1]->w, dt);
    
	vec3 axesA; quaternionGetYBasis(&axesA, &q0);
    vec3 axesB; quaternionGetYBasis(&axesB, &q1);

	float angle = estimateAngle(vecdot(&axesA, &axesB));
	if (angle<(minAngle-DEG2RAD(0.1f)))
	{
		vec3 axis;
		veccross(&axis, &axesA, &axesB);    // NB: positive rotation about axis rotates A towards B
		vecnormalise(&axis);
		float denominator = ConstraintCalcRotationalDenominator(c->body[0], c->body[1], &axis);
		vecscale(solution, &axis, MIN(DEG2RAD(90.f/dt), angle - minAngle) / (dt * denominator));
	}
	else
	{
		veczero(solution);
	}
}

struct Simulation
{
public:
	void Init();
	void Draw();
	void Update(float dt);

	void GetConstraintSpring(ConstraintAngSpring* angSpring, Constraint* c, float dt);

public:
	enum { k_numBodies=4 };
	enum { k_numConstraints=k_numBodies };
	Body m_static;
	CapsuleBody m_cap[k_numBodies];
	Constraint m_constraint[k_numConstraints];

	float m_topOfChainY;
	float m_staticBodyPos0;
	bool m_update;
};


void Simulation::Init()
{
	m_update = true;
	m_topOfChainY = 1.5f;
	m_staticBodyPos0 = 1.f;
	const float radius = 0.05f;
	const float length = 0.3f;
	const float mass1 = 1.f;
	const float mass2 = 5.f;
	const float inertia = mass1*length*length*0.25f * 0.1f;
	const float spring = 1.0f;
	const float spring2 = 100.f;
	const float damp = 0.f;
	const float damp2 = 0.0f;
	const float z = 0.f;
	const vec3 topOfChain = {0.f, m_topOfChainY, z};
	const vec3 capsuleOffset = {0.f, radius + 0.5f*length, 0.f};
	
	// Static/Kinematic world body
	BodyInit(&m_static);
	m_static.pos.y = m_staticBodyPos0;
	
	// Body 0
	BodyInit(&m_cap[0].body, mass1, inertia);
	vecsub(&m_cap[0].body.pos, &topOfChain, &capsuleOffset);
	vecset(&m_cap[0].body.v, 0.f, 0.f, 0.f);
	m_cap[0].body.w.z = 0.f;
	m_cap[0].radius = radius;
	m_cap[0].length = length;

	/////////////////////////////////////////////
	// Constraint 1 - attach body to static world
	/////////////////////////////////////////////
	m_constraint[0].body[0] = &m_static;
	m_constraint[0].offset[0] = topOfChain;
	vecsub(&m_constraint[0].offset[0], &m_constraint[0].offset[0], &m_static.pos);
	m_constraint[0].body[1] = &m_cap[0].body;
	m_constraint[0].offset[1] = capsuleOffset;
	m_constraint[0].spring = spring;
	m_constraint[0].damp = damp;
	
	/////////////////////////////////////////////
	// Make the chain
	/////////////////////////////////////////////
	for (int i=1; i<k_numBodies; i++)
	{
		float mass = mass1 + (mass2 - mass1)*(float)(i+1)/(float)k_numBodies;
		BodyInit(&m_cap[i].body, mass, inertia);
		vecsubscale(&m_cap[i].body.pos, &m_cap[i-1].body.pos, &capsuleOffset, 2.f);
		vecset(&m_cap[i].body.v, 0.f, 0.f, 0.f);
		m_cap[i].body.w.z = 0.f;
		m_cap[i].radius = radius;
		m_cap[i].length = length;
		
		m_constraint[i].body[0] = &m_cap[i-1].body;
		m_constraint[i].offset[0] = vecneg(capsuleOffset);
		m_constraint[i].body[1] = &m_cap[i].body;
		m_constraint[i].offset[1] = capsuleOffset;
		m_constraint[i].spring = spring + (spring2 - spring)*(float)(i+1)/(float)k_numBodies;
		m_constraint[i].damp = damp + (damp2 - damp)*(float)(i+1)/(float)k_numBodies;
	}
}

void Simulation::Draw()
{
	vec4 colour = {1.f, 1.f, 1.f, 1.f};
	BodyDraw(&m_static, &colour);
	for (int i=0; i<k_numBodies; i++)
		CapsuleBodyDraw(&m_cap[i], &colour);

	for (int c=0; c<k_numConstraints; c++)
		ConstraintDraw(&m_constraint[c]);

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
}

inline float physicsSpringGetVelocity(float x, float velocity, float dt, float k, float b)
{
    float a = dt*k;
    return (velocity - a*x) / (1.f + dt*(a + b));
}

void Simulation::GetConstraintSpring(ConstraintAngSpring* angSpring, Constraint* c, float dt)
{
	// This is a hack and assumes that the y-axis of the bodies is to be aligned
    const float smallvalue = 1e-6f;
    vec3 axesA; quaternionGetYBasis(&axesA, &c->body[0]->rot);
    vec3 axesB; quaternionGetYBasis(&axesB, &c->body[1]->rot);
	veccross(&angSpring->axis, &axesA, &axesB);	// NB: a positive rotation about this axis, rotates axisA towards axisB
	float axisLenSq = vecsizesq(&angSpring->axis);
	if (axisLenSq>smallvalue)
	{
		float dot = vecdot(&axesA, &axesB);
		float angle = estimateAngle(dot);
		vec3 relangvel;
		vecsub(&relangvel, &c->body[1]->w, &c->body[0]->w);
		vecscale(&angSpring->axis, &angSpring->axis, 1.f/sqrtf(axisLenSq));
		float relSpeed = vecdot(&relangvel, &angSpring->axis);
		angSpring->targetSpeed = -approxOneExp(dt*c->spring)*angle/dt;//physicsSpringGetVelocity(angle, relSpeed, dt, c->spring, c->damp);
		//angSpring->targetSpeed = physicsSpringGetVelocity(angle, relSpeed, dt, c->spring, c->damp);
		angSpring->denominator = ConstraintCalcRotationalDenominator(c->body[0], c->body[1], &angSpring->axis);
	
		// off axis response
		//vecscale(&response0, &c->body[0]->rot.v[1].v3, c->body[0]->invInertia);
		//vecscale(&response1, &c->body[0]->rot.v[1].v3, c->body[1]->invInertia);
		angSpring->denominator2 = 1.f;//vecdot(&c->body[0]->rot.v[1].v3, &response0) + vecdot(&c->body[0]->rot.v[1].v3, &response1);
	}
	else
	{
		veczero(&angSpring->axis);
		angSpring->denominator = 1.f;
		angSpring->denominator2 = 1.f;
	}
}

float sin_clipped(float x, float t=0.8f)
{
	float s = sinf(x);
	if (fabsf(s)<t)
	{
		return (float)fsgnf(s)*fabsf(s)*(1.f/t);
	}
	return fsgnf(s);
}

//float sin_clipped(float x, float a=0.1f)
//{
//	float s = sinf(x);
//	float as = fabsf(s);
//	return fsgnf(s) * as / ( a + as) * (1.f + a);
//}

void Simulation::Update(float dt)
{
	vec3 pos0 = m_static.pos;
	if (m_update)
	{
		static float ftime = 0.f;
		m_static.pos.y = m_staticBodyPos0 + 0.08f*sin_clipped(3.5f*ftime);
		m_static.pos.x = 0.3f*sin_clipped(2.f*ftime, 1.f);
		m_static.pos.z = 0.17f*sin_clipped(1.f*ftime, 0.7f);
		m_static.w.z = 1.0f * sin_clipped(2.f*ftime);
		m_static.w.x = 2.0f * sin_clipped(3.f*ftime);
		m_static.w.y = 0.0f;
		ftime += 2.f*dt;
	}
	else
	{
		veczero(&m_static.v);
		veczero(&m_static.w);
	}

	// Add gravity
	for (int i=0; i<k_numBodies; i++)
	{
		if (m_cap[i].body.mass>0)
		{
			Body* b = &m_cap[i].body;
			b->v.y -= 10.f * dt;
			vec3 yaxis;
			quaternionGetYBasis(&yaxis, &b->rot);
			vecsubscale(&b->w, &b->w, &yaxis, vecdot(&b->w, &yaxis));
		}
	}
	
	// Randomly poke the body
	static float time = 0.f;
	time +=dt;
	//if (time > 3.f)
	//{
	//	time = 0.f;
	//	float v = (frnd()>0.5f) ? 2.f : - 1.f;
	//	for (int i=0; i<k_numBodies; i++)
	//		m_cap[i].body.v.x += v;
	//}

	// Springs
	ConstraintAngSpring springs[k_numConstraints];
	for (int c=0; c<k_numConstraints; c++)
	{
		GetConstraintSpring(&springs[c], &m_constraint[c], dt);
	}

	for (int c=0; c<k_numConstraints; c++)
		ConstraintPrepare(&m_constraint[c], dt);

	const int numIterations = 16;
	float subDt = dt / (float)numIterations;
	for (int repeat=0; repeat<numIterations; repeat++)
	{
		for (int c=0; c<k_numConstraints; c++)
		{
			Constraint* con = &m_constraint[c];
			if(0){
				vec3 angSpringImpulse;
				ConstraintAngSpringSolve(&angSpringImpulse, con, &springs[c]);
				vecscale(&angSpringImpulse, &angSpringImpulse, c==0 ? 0.f : 0.8f);
				BodyApplyAngImpulse(con->body[0], &angSpringImpulse);
				angSpringImpulse = vecneg(angSpringImpulse);
				BodyApplyAngImpulse(con->body[1], &angSpringImpulse);
			}
			if(0){
				vec3 angLimitImpulse;
				ConstraintAngLimitSolve(&angLimitImpulse, con, DEG2RAD((float)c*30.f), dt);
				vecscale(&angLimitImpulse, &angLimitImpulse, c==0 ? 0.8f : 0.8f);
				BodyApplyAngImpulse(con->body[0], &angLimitImpulse);
				angLimitImpulse = vecneg(angLimitImpulse);
				BodyApplyAngImpulse(con->body[1], &angLimitImpulse);
			}
			{
				vec3 impulse;
				float erp = repeat == (numIterations-1) ? 1.f : 0.65f;
				ConstraintSolve(&impulse, &m_constraint[c]);
				vecscale(&impulse, &impulse, erp);
				BodyApplyImpulse(m_constraint[c].body[0], &impulse, &m_constraint[c].woffset[0]);
				impulse = vecneg(impulse);
				BodyApplyImpulse(m_constraint[c].body[1], &impulse, &m_constraint[c].woffset[1]);
			}
		}
	}

	BodyUpdate(&m_static, dt);
	for (int i=0; i<k_numBodies; i++)
		BodyUpdate(&m_cap[i].body, dt);
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

	if (key=='#') g_simulation.m_update = !g_simulation.m_update;

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
#if 0
	mtx rot;
	matrixRotZ(&rot, DEG2RAD(45.f));

	float dataA[] = {
		1,2,3,
		2,1,6,
		3,6,1,
	};

	float dataB[] = {
		0,
		10,
		0,
	};
	
	MMatrix* matA = MMatrixCreate(3,3,0);
	MMatrix* matB = MMatrixCreate(3,1,0);
	MMatrix* matC = MMatrixCreate(3,1,0);
	
	MMatrixInit(matA, dataA);
	MMatrixInit(matB, dataB);

	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			float* v = (float*)&rot.v[i];
			matA->data[i*3+j] = v[j];
		}
	}
	
	printf("A:\n");
	MMatrixDump(matA);
	// printf("B:\n");
	// MMatrixDump(matB);

	// MMatrixMul(matC, matA, matB);
	// 
	// printf("C:\n");
	// MMatrixDump(matC);

	MMatrix* invMatA = MMatrixCreate(3,3,0);
	MMatrixGaussJordanInvert(invMatA, matA);
	printf("invA:\n");
	MMatrixDump(invMatA);

	MMatrixDelete(matA);
	MMatrixDelete(invMatA);
	MMatrixDelete(matB);
	MMatrixDelete(matC);

	exit(0);
#endif
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


