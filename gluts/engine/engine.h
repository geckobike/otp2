#ifndef __ENGINE_H__
#define __ENGINE_H__

typedef short s16;
typedef unsigned short u16;
typedef unsigned int u32;
typedef int s32;

#define PI	    (3.1415926535897931159979635f)
#define PIx2	    (6.2831853071795862319959270f)
#define PI_D2	    (1.57079632679489661923f)
#define SQRT1_2	    (0.707106781186547524401f)
#define DEG2RAD(a)  ((a)*PI/180.0f)
#define RAD2DEG(a)  ((a)*180.0f/PI)

#define FLOAT_MAX	(3.4028e+38f)
#define FLOAT_MIN	(-FLOAT_MAX)
#define FLOAT_EPSILON	(1.1920929e-7f)

#define SQR(x) ((x)*(x))

#define VECALIGN 0
#define MTXALIGN 0
#define __ALIGN(x,y) x

union vec3_u
{
    __ALIGN(float v[3], VECALIGN);

    struct
    {
	float x;
	float y;
	float z;
    };
};

typedef vec3_u vec3;

union vec4_u
{
    __ALIGN(float v[4], VECALIGN);
	
    struct
    {
	union
	{
	    vec3 v3;

	    struct
	    {
		float x;
		float y;
		float z;
	    };
	};

        float w;
    };

    struct
    {
	float r;
	float g;
	float b;
	float a;
    };

#if PLATFROM_CELL
    vec_float4 v128;
#endif
};

typedef vec4_u vec4;

union mtx_u
{
    __ALIGN(float f[4][4], MTXALIGN);
    __ALIGN(vec4 v[4], MTXALIGN);
};

typedef mtx_u mtx;

union quaternion_u
{
    __ALIGN(float q[4], QUATALIGN);

    struct
    {
	union
	{
	    vec3 v;

	    struct
	    {
		float x;
		float y;
		float z;
	    };
	};

	float w;
    };
};

typedef union quaternion_u quaternion;

typedef int boolean;

#define FALSE 0
#define TRUE (!FALSE)


extern void vecaddscale(vec3 *v, const vec3 *v1, const float f);
extern void vecaddscale(vec3 *v, const vec3 *v1, const vec3 *v2, const float f);
extern float vecsizesq(const vec3 *v);
extern float vecdistsq(const vec3 *v1, const vec3 *v2);
extern float vecsize(const vec3 *v);
extern float vecdist(const vec3 *v1, const vec3 *v2);
extern void veccpy(vec3* a, const vec3* b);
extern void vecscale(vec3* a, float s);
extern void vecscale(vec3* a, const vec3* b, float s);
extern void vecadd(vec3* v, const vec3* a);
extern void vecadd(vec3* v, const vec3* a, const vec3* b);
extern void vecsub(vec3* v, const vec3* a);
extern void vecsub(vec3* v, const vec3* a, const vec3* b);
extern void vecsubscale(vec3 *v, const vec3 *v1, const vec3 *v2, const float f);
extern float vecdot(const vec3 *v1, const vec3 *v2);
extern void vecnormalise(vec3* v);
extern void vecnormalise(vec3* v, const vec3* v1);
extern void vecset(vec3* v, float x, float y, float z);
extern void veczero(vec3* a);
extern float sgn(float f);
extern void veccpy4(vec4 *v, const vec4 *v1);

static inline void vecmidpoint(vec3 *v, const vec3 *v1, const vec3 *v2)
{
    v->x = (v1->x + v2->x) * 0.5f;
    v->y = (v1->y + v2->y) * 0.5f;
    v->z = (v1->z + v2->z) * 0.5f;
}

static void vecneg(vec3 *v, const vec3 *v1)
{
    v->x = -v1->x;
    v->y = -v1->y;
    v->z = -v1->z;
}

static void vecnegcpy(vec3 *v, const vec3 *v1)
{
    v->x = -v1->x;
    v->y = -v1->y;
    v->z = -v1->z;
}

static void veccross(vec3 *v, const vec3 *v1, const vec3 *v2)
{
    float x = (v1->y*v2->z) - (v1->z*v2->y);
    float y = (v1->z*v2->x) - (v1->x*v2->z);
    float z = (v1->x*v2->y) - (v1->y*v2->x);
    v->x = x;
    v->y = y;
    v->z = z;
}

static void vecset4(vec4 *v, const float x, const float y, const float z, float w)
{
    v->x = x;
    v->y = y;
    v->z = z;
    v->w = w;
}

inline vec3 operator + (const vec3& a, const vec3& b)
{
	vec3 out = {{(a.x + b.x), (a.y + b.y), (a.z + b.z)}};
	return out;
}

inline vec3 operator - (const vec3& a, const vec3& b)
{
	vec3 out = {{(a.x - b.x), (a.y - b.y), (a.z - b.z)}};
	return out;
}

inline vec3 vecneg(const vec3& a)
{
	vec3 out = {{(-a.x), (-a.y), (-a.z)}};
	return out;
}

#define fapprox(f1,f2, epsilon)	    (fabsf((f2)-(f1))<(epsilon))
#define fapprox0(f, epsilon)	    (fabsf(f)<(epsilon))
#define fapx	    fapprox
#define fapx0	    fapprox0

#define vecisunit(v)	    (fapx(vecsizesq(v),1.f,0.005f))

#define XYZ(v) (v).x, (v).y, (v).z
#define veclist(v) (v)[0], (v)[1], (v)[2], (v)[3]

#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))

extern void vec3mtx43mulvec3(vec3 *v, const mtx *m, const vec3 *v1);
extern void vec4mtx44mulvec3(vec4 *v, const mtx *m, const vec3 *v1);
extern void vec3mtx33mulvec3(vec3 *v, const mtx *m, const vec3 *v1);
extern void vec3tranmtx33mulvec3(vec3 *v, const mtx *m, const vec3 *v1);
extern void vec3tranmtx43mulvec3(vec3 *v, const mtx *m, const vec3 *v1);
extern void vec4mtx44mulvec4(vec4 *v, const mtx *m, const vec4 *v1);
extern void vec4mtx44mulvec4pre(vec4 *v, const mtx *m, const vec4 *v1);

#define vec3list(v) (v).x, (v).y, (v).z
#define vec3listp(v) (v)->x, (v)->y, (v)->z

extern void matrixIdent(mtx *m);
extern void matrixZero(mtx *m);
extern void matrixAbs(mtx *m1, const mtx *m0);
extern boolean matrixCompare(const mtx *m0, const mtx *m1);
extern void matrixScale(mtx *m, float scale);
extern void matrixScaleXYZ(mtx *m, float sx, float sy, float sz);
extern void matrixDoScaleXYZ(mtx *m1, const mtx *m0, float sx, float sy, float sz);
extern float matrixGetDeterminant(const mtx *m);
extern void matrixCopy(mtx *m, const mtx *m1);
extern void matrixCopy33(mtx *outMtx, const mtx *inMtx);
extern void matrixInvert(mtx *m, const mtx *m1);
extern void matrixInvert44(mtx *m, const mtx *m1);
extern void matrixTranspose(mtx *m, const mtx *m1);
extern void matrixTransposeWithTrans(mtx *outMtx, const mtx *inMtx);
extern void matrixLook(mtx *m, const vec3 *pos, const vec3 *dir, const vec3 *up, const vec3 *right);
extern void matrixLookAt(mtx *m, const vec3 *pos, const vec3 *centre, const vec3 *up);
extern void matrixLookNonTransposed(mtx *m, const vec3 *pos, const vec3 *dir, const vec3 *up, const vec3 *right);
extern void matrixLookAtNonTransposed(mtx *m, const vec3 *pos, const vec3 *lookatpoint, const vec3 *up);
extern void matrixPerspective(mtx *m, float fovy, float aspect, float zn, float zf);
extern void matrixFrustum(mtx *m, float left, float right, float bottom, float top, float zNear, float zFar);
extern void matrixOrtho(mtx *mat,  float left,  float right,  float bottom,  float top,  float zNear,  float zFar);
extern void matrixRot(mtx *m, const vec3 *u, float angle);
extern void matrixRotCS(mtx *m, vec3 *u, float cosAngle, float sinAngle);
extern mtx* matrixRotX(mtx *m, float angle);
extern mtx* matrixRotY(mtx *m, float angle);
extern mtx* matrixRotZ(mtx *m, float angle);
extern void matrixRotXYZ(mtx *m, float anglex, float angley, float anglez);
extern void matrixRotZXY(mtx *m, float anglex, float angley, float anglez);
extern void matrixRotZYX(mtx *m, float anglex, float angley, float anglez);
extern void matrixTransRotXYZ(mtx *m, float tx, float ty, float tz, float anglex, float angley, float anglez);
extern void matrixTransRotZXY(mtx *m, float tx, float ty, float tz, float anglex, float angley, float anglez);
extern void matrixTrans(mtx *m, float x, float y, float z);
extern void matrixTransv(mtx *m, const vec3* vec);
extern void matrixSetTrans(mtx *m, float x, float y, float z);
extern void matrixSetTransv(mtx *outMtx, const vec3 *inTrans);
extern const vec3* matrixGetTrans(const mtx *m);
extern void matrixSetRightv(mtx *outMtx, const vec3 *inRight);
extern const vec3* matrixGetRight(const mtx *inMtx);
extern void matrixSetUpv(mtx *outMtx, const vec3 *inUp);
extern const vec3* matrixGetUp(const mtx *inMtx);
extern void matrixSetDirv(mtx *outMtx, const vec3 *inDir);
extern const vec3* matrixGetDir(const mtx *inMtx);

#if PLATFORM_PC || PLATFORM_LINUX
#define matrixMultiply(a,b,c) matrixMultiplyAligned(a,b,c)
#else
extern void matrixMultiply(mtx *m, const mtx *m1, const mtx *m2);
#endif
extern void matrixMultiplyAligned(mtx *m, const mtx *m1, const mtx *m2);
extern void matrix33Multiply33(mtx *m, const mtx *m1, const mtx *m2);
extern void matrix33Multiply33Tran1(mtx *m, const mtx *m1, const mtx *m2);
extern void matrixFrom2Axes(mtx *m, vec3 *v1, vec3 *v2);
extern void matrixToRotationZXY(const mtx *m, float *r1, float *r2, float *r3);
extern void matrixVecRot(const mtx *m, vec3 *vec);
extern void matrixVecRotCpy(vec3 *out, const mtx *m, const vec3 *vec);
extern void matrixVecRotTranspose(const mtx *m, vec3 *v);
extern void matrixVecMulCpy(vec3 *outVec, const mtx *m, const vec3 *v);
extern void matrixVecMul(const mtx *m, vec3 *vec);
extern void matrixVecMulTranspose(const mtx *m, vec3 *v);
extern void matrixVec3Mul4(const mtx *m, const vec3 *vec, vec4 *newvec);
extern void matrixFrom3Axis(mtx *outMtx, const vec3 *inX, const vec3 *inY, const vec3 *inZ);

extern void matrixPrintfTranspose(const mtx *m);
extern void matrixPrintf(const mtx *m);
extern void matrixFromXAxis( mtx *m, const vec3 *xaxis );
extern void matrixFromYAxis( mtx *m, const vec3 *yaxis );
extern void matrixFromYAxisWithPos( mtx *m, const vec3 *yaxis, const vec3 *inPos );
extern void matrixFromZAxis( mtx *m, const vec3 *zaxis );
extern void matrixFromZAxisWithPos( mtx *m, const vec3 *zaxis, const vec3 *inPos );

extern void matrixroll(union mtx_u *newmtx, vec3 *axis, float axialrot);
extern void matrixReNormalise(mtx* m);


extern void quaternionRotateVector(vec3 *v, const quaternion *q1, const vec3 *v1);
extern void quaternionSet(quaternion *q, float x, float y, float z, float w);
extern void quaternionIdent(quaternion *q);
extern void quaternionCopy(quaternion *q, const quaternion *q1);
extern void quaternionNegate(quaternion *q, const quaternion *q1);
extern void quaternionAdd(quaternion *q, const quaternion *q1, const quaternion *q2);
extern void quaternionSub(quaternion *q, const quaternion *q1, const quaternion *q2);
extern void quaternionToMatrix(mtx* m, const quaternion* q);
extern void quaternionNormalise(quaternion *q);
extern void quaternionMul(quaternion *q, const quaternion *q1, const quaternion *q2);
extern void quaternionGetXBasis(vec3* out, const quaternion* q);
extern void quaternionGetYBasis(vec3* out, const quaternion* q);
extern void quaternionGetZBasis(vec3* out, const quaternion* q);

#define rnd rand
inline float frnd()
{
	return ((float)rand()*(1.f/RAND_MAX));
}


/// TIMER
#ifdef ENABLE_TIMER

#include <time.h>
#include <sys/time.h>

struct Timer
{
	timeval time;
};

static void timerUpdate(Timer* timer)
{
	gettimeofday(&timer->time,NULL);
}

static double timerGetTickSinceLastUUpdate(Timer* timer)
{
    struct timeval now;
 
    long long sec;
    long long microsec;
    
    gettimeofday(&now,NULL);
    
    sec = now.tv_sec - timer->time.tv_sec;
    microsec = now.tv_usec - timer->time.tv_usec;
   
    return (double)sec + ( (double)microsec / (1e6) );
}
#endif



#define assert(x) { if (!(x)) { printf("ABORT: from file %s, line %d: expr failed: %s\n", __FILE__, __LINE__, #x); fflush(stdout); exit(0); } }



#endif
