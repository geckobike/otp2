#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <time.h>
//#include <sys/time.h>

#include "engine.h"
#ifndef assert
#define assert(expr)		    (void)sizeof(expr)
#endif
#define criticalAssert(expr)	    (void)sizeof(expr)
#define assertMsg(expr,msg)	    (void)sizeof(expr,msg)
#define criticalAssertMsg(expr,msg) (void)sizeof(expr,msg)
#define assertVerify(condition, msg)	(condition)
#define assertInit()
#define assertDisable(b)		    
#define assertSetCallback(f,d)	    (void)sizeof(f)
#define assertGetCallback(f,d)	    (void)sizeof(f)
#define assertIf(condition)	    if (condition)
#define assertIfMsg(condition,m)    if (condition)
#define assertNot(condition)	    if (condition)
#define assertNotMsg(condition,m)   if (condition)
#define assertElse(condition)	    if (!(condition))
#define assertElseMsg(condition,m)  if (!(condition))
#define assertSetProcessName(n)	    (void)sizeof(n)
#define assertGetProcessName()	    ""
#define assertIsUnitVector(v)	    (void)sizeof(v)

void vecaddscale(vec3 *v, const vec3 *v1, const float f)
{
    v->x += v1->x;
    v->y += v1->y;
    v->z += v1->z;
}

void vecaddscale(vec3 *v, const vec3 *v1, const vec3 *v2, const float f)
{
    v->x = v1->x + f*v2->x;
    v->y = v1->y + f*v2->y;
    v->z = v1->z + f*v2->z;
}

float vecsizesq(const vec3 *v)
{
    return (SQR(v->x) + SQR(v->y) + SQR(v->z));
}

float vecdistsq(const vec3 *v1, const vec3 *v2)
{
    return SQR(v1->x - v2->x) + SQR(v1->y - v2->y) + SQR(v1->z - v2->z);
}

float vecsize(const vec3 *v)
{
    return sqrtf(vecsizesq(v));
}

float vecdist(const vec3 *v1, const vec3 *v2)
{
    return sqrtf(vecdistsq(v1, v2));
}

void veccpy(vec3* a, const vec3* b)
{
    *a = *b;
}

void vecscale(vec3* a, float s)
{
    a->x*=s;
    a->y*=s;
    a->z*=s;
}

void vecscale(vec3* a, const vec3* b, float s)
{
    a->x = b->x * s;
    a->y = b->y * s;
    a->z = b->z * s;
}

void vecadd(vec3* v, const vec3* a)
{
    v->x += a->x;
    v->y += a->y;
    v->z += a->z;
}

void vecadd(vec3* v, const vec3* a, const vec3* b)
{
    v->x = a->x + b->x;
    v->y = a->y + b->y;
    v->z = a->z + b->z;
}

void vecsub(vec3* v, const vec3* a)
{
    v->x -= a->x;
    v->y -= a->y;
    v->z -= a->z;
}

void vecsub(vec3* v, const vec3* a, const vec3* b)
{
    v->x = a->x - b->x;
    v->y = a->y - b->y;
    v->z = a->z - b->z;
}

void vecsubscale(vec3 *v, const vec3 *v1, const vec3 *v2, const float f)
{
    v->x = v1->x - f*v2->x;
    v->y = v1->y - f*v2->y;
    v->z = v1->z - f*v2->z;
}

float vecdot(const vec3 *v1, const vec3 *v2)
{
    return (v1->x*v2->x) + (v1->y*v2->y) + (v1->z*v2->z);
}

void vecnormalise(vec3* v)
{
    float size = 1.0f/vecsize(v);
    vecscale(v, size);
}

void vecnormalise(vec3* v, const vec3* v1)
{
    float size = 1.0f/vecsize(v1);
    vecscale(v, v1, size);
}

void vecset(vec3* v, float x, float y, float z)
{
    v->x = x;
    v->y = y;
    v->z = z;
}
void veczero(vec3* a)
{
    a->x=0.f;
    a->y=0.f;
    a->z=0.f;
}

float sgn(float f)
{
    if(f<0.f)	return -1.f;
    if(f>0.f)	return 1.f;
    return(0.f);
}

void veccpy4(vec4 *v, const vec4 *v1)
{
    v->x = v1->x;
    v->y = v1->y;
    v->z = v1->z;
    v->w = v1->w;
}

void vec3mtx43mulvec3(vec3 *v, const mtx *m, const vec3 *v1)
{
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;
    
    v->x = x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0] + m->f[3][0];
    v->y = x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1] + m->f[3][1];
    v->z = x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2] + m->f[3][2];
}

void vec4mtx44mulvec3(vec4 *v, const mtx *m, const vec3 *v1)
{
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;

    v->x = x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0] + m->f[3][0];
    v->y = x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1] + m->f[3][1];
    v->z = x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2] + m->f[3][2];
    v->w = x*m->f[0][3] + y*m->f[1][3] + z*m->f[2][3] + m->f[3][3];
}

void vec3mtx33mulvec3(vec3 *v, const mtx *m, const vec3 *v1)
{
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;

    v->x = x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0];
    v->y = x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1];
    v->z = x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2];
}

void vec3tranmtx33mulvec3(vec3 *v, const mtx *m, const vec3 *v1)
{
    // 'tran' meaning transpose
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;

    v->x = x*m->f[0][0] + y*m->f[0][1] + z*m->f[0][2];
    v->y = x*m->f[1][0] + y*m->f[1][1] + z*m->f[1][2];
    v->z = x*m->f[2][0] + y*m->f[2][1] + z*m->f[2][2];
}

void vec3tranmtx43mulvec3(vec3 *v, const mtx *m, const vec3 *v1)
{
    // 'tran' meaning transpose
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;
    
    v->x = x*m->f[0][0] + y*m->f[0][1] + z*m->f[0][2] + m->f[0][3];
    v->y = x*m->f[1][0] + y*m->f[1][1] + z*m->f[1][2] + m->f[1][3];
    v->z = x*m->f[2][0] + y*m->f[2][1] + z*m->f[2][2] + m->f[2][3];
}


void vec4mtx44mulvec4(vec4 *v, const mtx *m, const vec4 *v1)
{
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;
    const float w = v1->w;

    v->x = x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0] + w*m->f[3][0];
    v->y = x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1] + w*m->f[3][1];
    v->z = x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2] + w*m->f[3][2];
    v->w = x*m->f[0][3] + y*m->f[1][3] + z*m->f[2][3] + w*m->f[3][3];
}

void vec4mtx44mulvec4pre(vec4 *v, const mtx *m, const vec4 *v1)
{
    const float x = v1->x;
    const float y = v1->y;
    const float z = v1->z;
    const float w = v1->w;

    v->x = x*m->f[0][0] + y*m->f[0][1] + z*m->f[0][2] + w*m->f[0][3];
    v->y = x*m->f[1][0] + y*m->f[1][1] + z*m->f[1][2] + w*m->f[1][3];
    v->z = x*m->f[2][0] + y*m->f[2][1] + z*m->f[2][2] + w*m->f[2][3];
    v->w = x*m->f[3][0] + y*m->f[3][1] + z*m->f[3][2] + w*m->f[3][3];
}

//================================================================================================================

void matrixIdent(mtx *m)
{
    m->f[0][0] = 1.0f;
    m->f[0][1] = 0.0f;
    m->f[0][2] = 0.0f;
    m->f[0][3] = 0.0f;

    m->f[1][0] = 0.0f;
    m->f[1][1] = 1.0f;
    m->f[1][2] = 0.0f;
    m->f[1][3] = 0.0f;

    m->f[2][0] = 0.0f;
    m->f[2][1] = 0.0f;
    m->f[2][2] = 1.0f;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}

void matrixZero(mtx *m)
{
    m->f[0][0] = 0.0f;
    m->f[0][1] = 0.0f;
    m->f[0][2] = 0.0f;
    m->f[0][3] = 0.0f;

    m->f[1][0] = 0.0f;
    m->f[1][1] = 0.0f;
    m->f[1][2] = 0.0f;
    m->f[1][3] = 0.0f;

    m->f[2][0] = 0.0f;
    m->f[2][1] = 0.0f;
    m->f[2][2] = 0.0f;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 0.0f;
}

void matrixAbs(mtx *m1, const mtx *m0)
{
    m1->f[0][0] = fabsf(m0->f[0][0]);
    m1->f[0][1] = fabsf(m0->f[0][1]);
    m1->f[0][2] = fabsf(m0->f[0][2]);
    m1->f[1][0] = fabsf(m0->f[1][0]);
    m1->f[1][1] = fabsf(m0->f[1][1]);
    m1->f[1][2] = fabsf(m0->f[1][2]);
    m1->f[2][0] = fabsf(m0->f[2][0]);
    m1->f[2][1] = fabsf(m0->f[2][1]);
    m1->f[2][2] = fabsf(m0->f[2][2]);
    m1->f[3][0] = 0.0f;
    m1->f[3][1] = 0.0f;
    m1->f[3][2] = 0.0f;
}

void matrixScale(mtx *m, float scale)
{
    m->f[0][0] = scale;
    m->f[0][1] = 0.0f;
    m->f[0][2] = 0.0f;
    m->f[0][3] = 0.0f;

    m->f[1][0] = 0.0f;
    m->f[1][1] = scale;
    m->f[1][2] = 0.0f;
    m->f[1][3] = 0.0f;

    m->f[2][0] = 0.0f;
    m->f[2][1] = 0.0f;
    m->f[2][2] = scale;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}

void matrixScaleXYZ(mtx *m, float sx, float sy, float sz)
{
    m->f[0][0] = sx;
    m->f[0][1] = 0.0f;
    m->f[0][2] = 0.0f;
    m->f[0][3] = 0.0f;

    m->f[1][0] = 0.0f;
    m->f[1][1] = sy;
    m->f[1][2] = 0.0f;
    m->f[1][3] = 0.0f;

    m->f[2][0] = 0.0f;
    m->f[2][1] = 0.0f;
    m->f[2][2] = sz;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}

void matrixDoScaleXYZ(mtx *m1, const mtx *m0, float sx, float sy, float sz)
{
    vecscale(&m1->v[0].v3, &m0->v[0].v3, sx);
    vecscale(&m1->v[1].v3, &m0->v[1].v3, sy);
    vecscale(&m1->v[2].v3, &m0->v[2].v3, sz);
}

float matrixGetDeterminant(const mtx *m)
{
    return m->f[0][0] * (m->f[1][1]*m->f[2][2] - m->f[1][2]*m->f[2][1])
	- m->f[0][1] * (m->f[1][0]*m->f[2][2] - m->f[1][2]*m->f[2][0])
	+ m->f[0][2] * (m->f[1][0]*m->f[2][1] - m->f[1][1]*m->f[2][0]);
}

void matrixCopy(mtx *m, const mtx *m1)
{
    memcpy(m, m1, sizeof(mtx));
}

// copy 3x3 rotation submatrix
void matrixCopy33(mtx *outMtx, const mtx *inMtx)
{
    veccpy4(&outMtx->v[0],&inMtx->v[0]);
    veccpy4(&outMtx->v[1],&inMtx->v[1]);
    veccpy4(&outMtx->v[2],&inMtx->v[2]);
}

// m = output
// m1 = input
void matrixInvert(mtx *m, const mtx *m1)
{
    int i, j;
    float rdet;
    mtx tmp;
    mtx *m0 = (m == m1) ? & tmp : m;

    rdet = matrixGetDeterminant(m1);

//    if (fabsf(rdet) < FLOAT_EPSILON)
    if (fabsf(rdet) == 0.0f)
    {
	return;
    }

    rdet = 1.0f/rdet;

    // Get adjoint of rotational submatrix
    m0->f[0][0] = m1->f[1][1]*m1->f[2][2] - m1->f[1][2]*m1->f[2][1];
    m0->f[0][1] = m1->f[0][2]*m1->f[2][1] - m1->f[0][1]*m1->f[2][2];
    m0->f[0][2] = m1->f[0][1]*m1->f[1][2] - m1->f[0][2]*m1->f[1][1];
    m0->f[1][0] = m1->f[1][2]*m1->f[2][0] - m1->f[1][0]*m1->f[2][2];
    m0->f[1][1] = m1->f[0][0]*m1->f[2][2] - m1->f[0][2]*m1->f[2][0];
    m0->f[1][2] = m1->f[0][2]*m1->f[1][0] - m1->f[0][0]*m1->f[1][2];
    m0->f[2][0] = m1->f[1][0]*m1->f[2][1] - m1->f[1][1]*m1->f[2][0];
    m0->f[2][1] = m1->f[0][1]*m1->f[2][0] - m1->f[0][0]*m1->f[2][1];
    m0->f[2][2] = m1->f[0][0]*m1->f[1][1] - m1->f[0][1]*m1->f[1][0];

    // Scale by reciprocal of determinant
    for (i = 0; i < 3; i++)
    {
	for (j = 0; j < 3; j++)
	{
	    m0->f[i][j] *= rdet;
	}
    }

    // Calculate translation row vector
    m0->f[3][0] = -(m1->f[3][0]*m0->f[0][0] + m1->f[3][1]*m0->f[1][0] + m1->f[3][2]*m0->f[2][0]);
    m0->f[3][1] = -(m1->f[3][0]*m0->f[0][1] + m1->f[3][1]*m0->f[1][1] + m1->f[3][2]*m0->f[2][1]);
    m0->f[3][2] = -(m1->f[3][0]*m0->f[0][2] + m1->f[3][1]*m0->f[1][2] + m1->f[3][2]*m0->f[2][2]);

    m0->f[0][3] =
    m0->f[1][3] =
    m0->f[2][3] = 0.0f;
    m0->f[3][3] = 1.0f;

    if (m == m1)
    {
	matrixCopy(m, m0);
    }
}

static inline float det2x2(float a, float b, float c, float d)
{
    float ans;
    ans = a * d - b * c;
    return ans;
}

static float det3x3(float  a1, float a2, float a3, float b1, float b2, float b3, float c1, float c2, float c3)
{
    float ans;
    ans = a1 * det2x2( b2, b3, c2, c3 )
        - b1 * det2x2( a2, a3, c2, c3 )
        + c1 * det2x2( a2, a3, b2, b3 );
    return ans;
}

static float det4x4( const mtx *m )
{
    float ans;
    float a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

    // assign to individual variable names to aid selecting correct elements
	
    a1 = m->f[0][0]; b1 = m->f[0][1]; 
    c1 = m->f[0][2]; d1 = m->f[0][3];

    a2 = m->f[1][0]; b2 = m->f[1][1]; 
    c2 = m->f[1][2]; d2 = m->f[1][3];

    a3 = m->f[2][0]; b3 = m->f[2][1]; 
    c3 = m->f[2][2]; d3 = m->f[2][3];

    a4 = m->f[3][0]; b4 = m->f[3][1]; 
    c4 = m->f[3][2]; d4 = m->f[3][3];

    ans = a1 * det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4)
        - b1 * det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4)
        + c1 * det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4)
        - d1 * det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);
    return ans;
}

static void adjoint(const mtx *in, mtx *out)
{
    float a1, a2, a3, a4, b1, b2, b3, b4;
    float c1, c2, c3, c4, d1, d2, d3, d4;

    // assign to individual variable names to aid selecting correct values

    a1 = in->f[0][0]; b1 = in->f[0][1]; 
    c1 = in->f[0][2]; d1 = in->f[0][3];

    a2 = in->f[1][0]; b2 = in->f[1][1]; 
    c2 = in->f[1][2]; d2 = in->f[1][3];

    a3 = in->f[2][0]; b3 = in->f[2][1];
    c3 = in->f[2][2]; d3 = in->f[2][3];

    a4 = in->f[3][0]; b4 = in->f[3][1]; 
    c4 = in->f[3][2]; d4 = in->f[3][3];

    // row column labeling reversed since we transpose rows & columns

    out->f[0][0]  =   det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
    out->f[1][0]  = - det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
    out->f[2][0]  =   det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
    out->f[3][0]  = - det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);
        
    out->f[0][1]  = - det3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4);
    out->f[1][1]  =   det3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4);
    out->f[2][1]  = - det3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4);
    out->f[3][1]  =   det3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4);
        
    out->f[0][2]  =   det3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4);
    out->f[1][2]  = - det3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4);
    out->f[2][2]  =   det3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4);
    out->f[3][2]  = - det3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4);
        
    out->f[0][3]  = - det3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3);
    out->f[1][3]  =   det3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3);
    out->f[2][3]  = - det3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3);
    out->f[3][3]  =   det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

void matrixInvert44(mtx *m, const mtx *min)
{
    // Full inverse of 4x4 matrix 'min' designed for accuracy rather than speed
    int i, j;
    float det;
    mtx	temp;
    
    // calculate the adjoint matrix
    adjoint(min, &temp);

    //  calculate the 4x4 determinant
    //  if the determinant is zero, 
    //  then the inverse matrix is not unique.
    det = det4x4(min);

    if (fabsf(det) < FLOAT_EPSILON)
    {
        assertMsg(0, "Non-singular matrix, no inverse!");
	matrixCopy(m, min);
	return;
    }

    // scale the adjoint matrix to get the inverse
    for (i=0; i<4; i++)
        for(j=0; j<4; j++)
	    m->f[i][j] = temp.f[i][j] / det;
}

// transpose the rotation part of the mtx, then invert the translation part
// this will produce the same results as a matrixInvert() for an orthonormalised mtx
void matrixTransposeWithTrans(mtx *outMtx, const mtx *inMtx)
{
    assert(outMtx!=inMtx);	// add new func (preferable) or extend this func to deal with it

    // transpose rotation submtx
    outMtx->f[0][0]=inMtx->f[0][0];
    outMtx->f[0][1]=inMtx->f[1][0];
    outMtx->f[0][2]=inMtx->f[2][0];
    outMtx->f[0][3]=0.f;
    outMtx->f[1][0]=inMtx->f[0][1];
    outMtx->f[1][1]=inMtx->f[1][1];
    outMtx->f[1][2]=inMtx->f[2][1];
    outMtx->f[1][3]=0.0f;
    outMtx->f[2][0]=inMtx->f[0][2];
    outMtx->f[2][1]=inMtx->f[1][2];
    outMtx->f[2][2]=inMtx->f[2][2];
    outMtx->f[2][3]=0.0f;
    outMtx->f[3][3]=1.0f;

    // invert position
    vec3mtx33mulvec3((vec3*)&outMtx->v[3],outMtx,(const vec3*)&inMtx->v[3]);
    vecneg(&outMtx->v[3].v3, &outMtx->v[3].v3);
}

void matrixTranspose(mtx *m, const mtx *m1)
{
    if (m == m1)
    {
#define SWAPF(a,b) { float tmp = (a); (a) = (b); (b) = tmp; }
    SWAPF(m->f[1][0], m->f[0][1]);
    SWAPF(m->f[2][0], m->f[0][2]);
    SWAPF(m->f[3][0], m->f[0][3]);

    SWAPF(m->f[2][1], m->f[1][2]);
    SWAPF(m->f[3][1], m->f[1][3]);

    SWAPF(m->f[3][2], m->f[2][3]);
#undef SWAPF
    }
    else
    {
	m->f[0][0] = m1->f[0][0];
	m->f[1][0] = m1->f[0][1];
	m->f[2][0] = m1->f[0][2];
	m->f[3][0] = m1->f[0][3];
	m->f[0][1] = m1->f[1][0];
	m->f[1][1] = m1->f[1][1];
	m->f[2][1] = m1->f[1][2];
	m->f[3][1] = m1->f[1][3];
	m->f[0][2] = m1->f[2][0];
	m->f[1][2] = m1->f[2][1];
	m->f[2][2] = m1->f[2][2];
	m->f[3][2] = m1->f[2][3];
	m->f[0][3] = m1->f[3][0];
	m->f[1][3] = m1->f[3][1];
	m->f[2][3] = m1->f[3][2];
	m->f[3][3] = m1->f[3][3];
    }
}


// the mtx produced by this is transposed, ie the x unit axis is in [0][0] [1][0] [2][0]
void matrixLook(mtx *m, const vec3 *pos, const vec3 *dir, const vec3 *up, const vec3 *right)
{
    m->f[0][0] = -right->x;
    m->f[0][1] = up->x;
    m->f[0][2] = -dir->x;
    m->f[0][3] = 0.0f;

    m->f[1][0] = -right->y;
    m->f[1][1] = up->y;
    m->f[1][2] = -dir->y;
    m->f[1][3] = 0.0f;

    m->f[2][0] = -right->z;
    m->f[2][1] = up->z;
    m->f[2][2] = -dir->z;
    m->f[2][3] = 0.0f;

    m->f[3][0] = ((pos->x*right->x) + (pos->y*right->y) + (pos->z*right->z));
    m->f[3][1] = -((pos->x*up->x)    + (pos->y*up->y)    + (pos->z*up->z));
    m->f[3][2] = ((pos->x*dir->x)   + (pos->y*dir->y)   + (pos->z*dir->z));
    m->f[3][3] = 1.0f;
}

void	matrixLookAt(mtx *m, const vec3 *pos, const vec3 *centre, const vec3 *up)
{
    vec3    dir;
    vec3    right;
    vec3    up2;
    vecsub(&dir, centre, pos);
    vecnormalise(&dir, &dir);
    veccross(&right, up, &dir);
    vecnormalise(&right, &right);
    veccross(&up2, &dir, &right);
    vecnormalise(&up2, &up2);
    matrixLook(m, pos, &dir, &up2, &right);
}

void matrixLookNonTransposed(mtx *m, const vec3 *pos, const vec3 *dir, const vec3 *up, const vec3 *right)
{
    m->f[0][0] = right->x;
    m->f[0][1] = right->y;
    m->f[0][2] = right->z;
    m->f[0][3] = 0.0f;

    m->f[1][0] = up->x;
    m->f[1][1] = up->y;
    m->f[1][2] = up->z;
    m->f[1][3] = 0.0f;

    m->f[2][0] = dir->x;
    m->f[2][1] = dir->y;
    m->f[2][2] = dir->z;
    m->f[2][3] = 0.0f;

    m->f[3][0] = pos->x;
    m->f[3][1] = pos->y;
    m->f[3][2] = pos->z;
    m->f[3][3] = 1.0f;
}

// generates a mtx with the z axis pointing from pos towards lookatpoint
void	matrixLookAtNonTransposed(mtx *m, const vec3 *pos, const vec3 *lookatpoint, const vec3 *up)
{
    vec3    dir;
    vec3    right;
    vec3    up2;
    vecsub(&dir, lookatpoint, pos);
    vecnormalise(&dir, &dir);
    veccross(&right, up, &dir);
    vecnormalise(&right, &right);
    veccross(&up2, &dir, &right);
    vecnormalise(&up2, &up2);
    matrixLookNonTransposed(m, pos, &dir, &up2, &right);
}

void matrixFrustum(mtx *m, float left, float right, float bottom, float top, float zNear, float zFar)
{
    float tmp;

    tmp     =  1.0f / (right - left);
    m->f[0][0] =  (2*zNear) * tmp;
    m->f[1][0] =  0.0f;
    m->f[2][0] =  (right + left) * tmp;
    m->f[3][0] =  0.0f;
              
    tmp  =     1.0f / (top - bottom);
    m->f[0][1] =  0.0f;
    m->f[1][1] =  (2*zNear) * tmp;
    m->f[2][1] =  (top + bottom) * tmp;
    m->f[3][1] =  0.0f;
              
    m->f[0][2] =  0.0f;
    m->f[1][2] =  0.0f;
              
    tmp  =     1.0f / (zFar - zNear);
              
    // scale z to (-w, w) range
    m->f[2][2] = -(zFar + zNear) * tmp;
    m->f[3][2] = -(2*zFar*zNear) * tmp;
              
    m->f[0][3] =  0.0f;
    m->f[1][3] =  0.0f;
    m->f[2][3] = -1.0f;
    m->f[3][3] =  0.0f;
}


void matrixPerspective(mtx *m, float fovy, float aspect, float zn, float zf)
{
/*    float invtan;
    float dz;

    fovy *= 0.5f;
    invtan = cosf(fovy)/sinf(fovy);
    dz = zf - zn;

    m->f[0][0] = invtan/aspect;
    m->f[0][1] = 0.0f;
    m->f[0][2] = 0.0f;
    m->f[0][3] = 0.0f;

    m->f[1][0] = 0.0f;
    m->f[1][1] = invtan;
    m->f[1][2] = 0.0f;
    m->f[1][3] = 0.0f;

    m->f[2][0] = 0.0f;
    m->f[2][1] = 0.0f;
    m->f[2][2] = -1.0f;
    m->f[2][3] = -1.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = -2.0f*zn;
    m->f[3][3] = 0.0f;
*/

    float invtan;
    float dz;

    fovy *= 0.5f;
    invtan = cosf(fovy)/sinf(fovy);
    dz = zf - zn;

    m->f[0][0] = invtan/aspect;
    m->f[0][1] = 0.0f;
    m->f[0][2] = 0.0f;
    m->f[0][3] = 0.0f;

    m->f[1][0] = 0.0f;
    m->f[1][1] = invtan;
    m->f[1][2] = 0.0f;
    m->f[1][3] = 0.0f;

    m->f[2][0] = 0.0f;
    m->f[2][1] = 0.0f;
    m->f[2][2] = -(zf + zn)/dz;
    m->f[2][3] = -1.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = -(2.0f*zf*zn)/dz;
    m->f[3][3] = 0.0f;
}

void matrixOrtho(mtx *mat,  float left,  float right,  float bottom,  float top,  float zNear,  float zFar)
{
    // Column major as gl
    mat->f[0][0]    =	2.f/(right-left);
    mat->f[0][1]    =   0.f;
    mat->f[0][2]    =   0.f;
    mat->f[0][3]    =   0.f;

    mat->f[1][0]    =   0.f;
    mat->f[1][1]    =   2.f/(top-bottom);
    mat->f[1][2]    =   0.f;
    mat->f[1][3]    =   0.f;

    mat->f[2][0]    =   0.f;
    mat->f[2][1]    =   0.f;
    mat->f[2][2]    =   -2.f/(zFar-zNear);
    mat->f[2][3]    =   0.f;

    mat->f[3][0]    =   -(right+left)/(right-left);
    mat->f[3][1]    =   -(top+bottom)/(top-bottom);
    mat->f[3][2]    =   -(zFar+zNear)/(zFar-zNear);
    mat->f[3][3]    =	1.0f;
}

void matrixRot(mtx *m, const vec3 *u, float angle)
{
    float c = cosf(angle);
    float s = sinf(angle);
    float v = 1.0f - c;

    float vxx = v*u->x*u->x;
    float vxy = v*u->x*u->y;
    float vxz = v*u->x*u->z;
    float vyy = v*u->y*u->y;
    float vyz = v*u->y*u->z;
    float vzz = v*u->z*u->z;

    float sx = s*u->x;
    float sy = s*u->y;
    float sz = s*u->z;

    m->f[0][0] = vxx + c;
    m->f[0][1] = vxy + sz;
    m->f[0][2] = vxz - sy;
    m->f[0][3] = 0.0f;

    m->f[1][0] = vxy - sz;
    m->f[1][1] = vyy + c;
    m->f[1][2] = vyz + sx;
    m->f[1][3] = 0.0f;

    m->f[2][0] = vxz + sy;
    m->f[2][1] = vyz - sx;
    m->f[2][2] = vzz + c;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}

void matrixRotCS(mtx *m, vec3 *u, float cosAngle, float sinAngle)   // Use this version if you have already calculated cosAngle and sinAngle
{
    float v = 1.0f - cosAngle;

    float vxx = v*u->x*u->x;
    float vxy = v*u->x*u->y;
    float vxz = v*u->x*u->z;
    float vyy = v*u->y*u->y;
    float vyz = v*u->y*u->z;
    float vzz = v*u->z*u->z;

    float sx = sinAngle*u->x;
    float sy = sinAngle*u->y;
    float sz = sinAngle*u->z;

    m->f[0][0] = vxx + cosAngle;
    m->f[0][1] = vxy + sz;
    m->f[0][2] = vxz - sy;
    m->f[0][3] = 0.0f;

    m->f[1][0] = vxy - sz;
    m->f[1][1] = vyy + cosAngle;
    m->f[1][2] = vyz + sx;
    m->f[1][3] = 0.0f;

    m->f[2][0] = vxz + sy;
    m->f[2][1] = vyz - sx;
    m->f[2][2] = vzz + cosAngle;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}

mtx* matrixRotX(mtx *m, float angle)
{
    float s,c;

    matrixIdent(m);
    s=sinf(angle);
    c=cosf(angle);
    m->f[1][1]=c;
    m->f[1][2]=s;
    m->f[2][1]=-s;
    m->f[2][2]=c;

    return m;
}

mtx* matrixRotY(mtx *m, float angle)
{
    float s,c;

    matrixIdent(m);
    s=sinf(angle);
    c=cosf(angle);
    m->f[2][2]=c;
    m->f[2][0]=s;
    m->f[0][2]=-s;
    m->f[0][0]=c;

    return m;
}

mtx* matrixRotZ(mtx *m, float angle)
{
    float s,c;

    matrixIdent(m);
    s=sinf(angle);
    c=cosf(angle);
    m->f[0][0]=c;
    m->f[0][1]=s;
    m->f[1][0]=-s;
    m->f[1][1]=c;

    return m;
}

#if 0
void matrixRotXYZ(mtx *m, float anglex, float angley, float anglez)
{
    // this is code for factoring as Ry Rx Rz
    float sx=sinf(anglex);
    float cx=cosf(anglex);
    float sy=sinf(angley);
    float cy=cosf(angley);
    float sz=sinf(anglez);
    float cz=cosf(anglez);

    float sxsy=sx*sy;
    float sxcy=sx*cy;
    
    m->f[0][0]=(cz*cy)+(sz*sxsy);
    m->f[0][1]=(sz*cx);
    m->f[0][2]=(cz*-sy)+(sz*sxcy);
    m->f[0][3]=0.0f;
    m->f[1][0]=(-sz*cy)+(cz*sxsy);
    m->f[1][1]=(cz*cx);
    m->f[1][2]=(sz*sy)+(cz*sxcy);
    m->f[1][3]=0.0f;
    m->f[2][0]=cx*sy;
    m->f[2][1]=-sx;
    m->f[2][2]=cx*cy;
    m->f[2][3]=0.0f;
    m->f[3][0]=0.0f;
    m->f[3][1]=0.0f;
    m->f[3][2]=0.0f;
    m->f[3][3]=1.0f;
}
#endif

void matrixRotXYZ(mtx *m, float anglex, float angley, float anglez)
{
    // this is code for factoring as Rz Ry Rx
    float sx=sinf(anglex);
    float cx=cosf(anglex);
    float sy=sinf(angley);
    float cy=cosf(angley);
    float sz=sinf(anglez);
    float cz=cosf(anglez);

    m->f[0][0] = cy*cz;
    m->f[0][1] = cy*sz;
    m->f[0][2] = -sy;
    m->f[0][3] = 0.0f;

    m->f[1][0] = cz*sx*sy - cx*sz;
    m->f[1][1] = cx*cz + sx*sy*sz;
    m->f[1][2] = cy*sx;
    m->f[1][3] = 0.0f;

    m->f[2][0] = cx*cz*sy + sx*sz;
    m->f[2][1] = -cz*sx + cx*sy*sz;
    m->f[2][2] = cx*cy;
    m->f[2][3] = 0.0f;

    m->f[3][0]=0.0f;
    m->f[3][1]=0.0f;
    m->f[3][2]=0.0f;
    m->f[3][3]=1.0f;
}

void matrixRotZYX(mtx *m, float anglex, float angley, float anglez)
{
    // this is code for factoring as Rx Ry Rz
    float sx=sinf(anglex);
    float cx=cosf(anglex);
    float sy=sinf(angley);
    float cy=cosf(angley);
    float sz=sinf(anglez);
    float cz=cosf(anglez);

    m->f[0][0] = cy*cz;
    m->f[0][1] = cz*sx*sy + cx*sz;
    m->f[0][2] = -cx*cz*sy + sx*sz;
    m->f[0][3] = 0.0f;

    m->f[1][0] = -cy*sz;
    m->f[1][1] = cx*cz - sx*sy*sz;
    m->f[1][2] = cz*sx + cx*sy*sz;
    m->f[1][3] = 0.0f;

    m->f[2][0] = sy;
    m->f[2][1] = -cy*sx;
    m->f[2][2] = cx*cy;
    m->f[2][3] = 0.0f;

    m->f[3][0]=0.0f;
    m->f[3][1]=0.0f;
    m->f[3][2]=0.0f;
    m->f[3][3]=1.0f;
}

void matrixRotZXY(mtx *m, float anglex, float angley, float anglez)
{
    // this is code for factoring as Ry Rx Rz
    float sx=sinf(anglex);
    float cx=cosf(anglex);
    float sy=sinf(angley);
    float cy=cosf(angley);
    float sz=sinf(anglez);
    float cz=cosf(anglez);

    float sxsy=sx*sy;
    float sxcy=sx*cy;

	// These are needed because of a compiler bug, optimisation was screwing up the
	// the floating point stack. Putting these values into intermediates fixs the problem
	float sxsysz = sz * sxsy;
	float szsxcy = sz * sxcy;
	float czsxsy = cz * sxsy;
	float czsxcy = cz * sxcy;

    m->f[0][0]=(cz*cy)+sxsysz;
    m->f[0][1]=(sz*cx);
    m->f[0][2]=(cz*-sy)+(szsxcy);
    m->f[0][3]=0.0f;
    
    m->f[1][0]=(-sz*cy)+(czsxsy);
    m->f[1][1]=(cz*cx);
    m->f[1][2]=(sz*sy)+(czsxcy);
    m->f[1][3]=0.0f;
    
    m->f[2][0]=cx*sy;
    m->f[2][1]=-sx;
    m->f[2][2]=cx*cy;
    m->f[2][3]=0.0f;

    m->f[3][0]=0.0f;
    m->f[3][1]=0.0f;
    m->f[3][2]=0.0f;
    m->f[3][3]=1.0f;
}

void matrixTransRotXYZ(mtx *m, float tx, float ty, float tz, float anglex, float angley, float anglez)
{
    float sx=sinf(anglex);
    float cx=cosf(anglex);
    float sy=sinf(angley);
    float cy=cosf(angley);
    float sz=sinf(anglez);
    float cz=cosf(anglez);
    float sxsy=sx*sy;
    float sxcy=sx*cy;
    
    m->f[0][0] = cy*cz;
    m->f[0][1] = cy*sz;
    m->f[0][2] = -sy;
    m->f[0][3] = 0.0f;

    m->f[1][0] = cz*sx*sy - cx*sz;
    m->f[1][1] = cx*cz + sx*sy*sz;
    m->f[1][2] = cy*sx;
    m->f[1][3] = 0.0f;

    m->f[2][0] = cx*cz*sy + sx*sz;
    m->f[2][1] = -cz*sx + cx*sy*sz;
    m->f[2][2] = cx*cy;
    m->f[2][3] = 0.0f;

    m->f[3][0]=tx;
    m->f[3][1]=ty;
    m->f[3][2]=tz;
    m->f[3][3]=1.0f;
}

void matrixTransRotZXY(mtx *m, float tx, float ty, float tz, float anglex, float angley, float anglez)
{
    float sx=sinf(anglex);
    float cx=cosf(anglex);
    float sy=sinf(angley);
    float cy=cosf(angley);
    float sz=sinf(anglez);
    float cz=cosf(anglez);
    float sxsy=sx*sy;
    float sxcy=sx*cy;
    
    m->f[0][0]=(cz*cy)+(sz*sxsy);
    m->f[0][1]=(sz*cx);
    m->f[0][2]=(cz*-sy)+(sz*sxcy);
    m->f[0][3]=0.0f;
    m->f[1][0]=(-sz*cy)+(cz*sxsy);
    m->f[1][1]=(cz*cx);
    m->f[1][2]=(sz*sy)+(cz*sxcy);
    m->f[1][3]=0.0f;
    m->f[2][0]=cx*sy;
    m->f[2][1]=-sx;
    m->f[2][2]=cx*cy;
    m->f[2][3]=0.0f;
    m->f[3][0]=tx;
    m->f[3][1]=ty;
    m->f[3][2]=tz;
    m->f[3][3]=1.0f;
}

void matrixTrans(mtx *m, float x, float y, float z)
{
    matrixIdent(m);
    vecset(&m->v[3].v3, x, y, z);
}

void matrixTransv(mtx *m, const vec3* vec)
{
    matrixIdent(m);
    vecset(&m->v[3].v3, vec->x, vec->y, vec->z);
}

void matrixSetTrans(mtx *m, float x, float y, float z)
{
    vecset(&m->v[3].v3, x, y, z);
}

const vec3* matrixGetTrans(const mtx *m)
{
    return &m->v[3].v3;
}

void matrixSetTransv(mtx *outMtx, const vec3 *inTrans)
{
    veccpy(&outMtx->v[3].v3, inTrans);
}

void matrixSetRightv(mtx *outMtx, const vec3 *inRight)
{
    veccpy(&outMtx->v[0].v3, inRight);
}

const vec3* matrixGetRight(const mtx *inMtx)
{
    return &inMtx->v[0].v3;
}

void matrixSetUpv(mtx *outMtx, const vec3 *inUp)
{
    veccpy(&outMtx->v[1].v3, inUp);
}

const vec3* matrixGetUp(const mtx *inMtx)
{
    return &inMtx->v[1].v3;
}

void matrixSetDirv(mtx *outMtx, const vec3 *inDir)
{
    veccpy(&outMtx->v[2].v3, inDir);    
}

const vec3* matrixGetDir(const mtx *inMtx)
{
    return &inMtx->v[2].v3;
}

#if ENABLE_SSE
#if 0
void matrixMultiplyAligned(mtx *m, const mtx *m1, const mtx *m2)
{
    __m128 a, b, c, d;
    
    assert( (((u32)(m->v)) & 16) == 0 );
    assert( (((u32)(m1->v)) & 16) == 0 );
    
    a = _mm_add_ps(
		    _mm_add_ps(
		    _mm_add_ps(
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[0], m1->sm[0],_MM_SHUFFLE(0,0,0,0)), m2->sm[0]), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[0], m1->sm[0],_MM_SHUFFLE(1,1,1,1)), m2->sm[1])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[0], m1->sm[0],_MM_SHUFFLE(2,2,2,2)), m2->sm[2])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[0], m1->sm[0],_MM_SHUFFLE(3,3,3,3)), m2->sm[3]));
    b = _mm_add_ps(
		    _mm_add_ps(
		    _mm_add_ps(
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[1], m1->sm[1],_MM_SHUFFLE(0,0,0,0)), m2->sm[0]), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[1], m1->sm[1],_MM_SHUFFLE(1,1,1,1)), m2->sm[1])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[1], m1->sm[1],_MM_SHUFFLE(2,2,2,2)), m2->sm[2])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[1], m1->sm[1],_MM_SHUFFLE(3,3,3,3)), m2->sm[3]));
    c = _mm_add_ps(
		    _mm_add_ps(
		    _mm_add_ps(
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[2], m1->sm[2],_MM_SHUFFLE(0,0,0,0)), m2->sm[0]), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[2], m1->sm[2],_MM_SHUFFLE(1,1,1,1)), m2->sm[1])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[2], m1->sm[2],_MM_SHUFFLE(2,2,2,2)), m2->sm[2])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[2], m1->sm[2],_MM_SHUFFLE(3,3,3,3)), m2->sm[3]));
    d = _mm_add_ps(
		    _mm_add_ps(
		    _mm_add_ps(
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[3], m1->sm[3],_MM_SHUFFLE(0,0,0,0)), m2->sm[0]), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[3], m1->sm[3],_MM_SHUFFLE(1,1,1,1)), m2->sm[1])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[3], m1->sm[3],_MM_SHUFFLE(2,2,2,2)), m2->sm[2])), 
		    _mm_mul_ps(_mm_shuffle_ps(m1->sm[3], m1->sm[3],_MM_SHUFFLE(3,3,3,3)), m2->sm[3]));
    
    m->sm[0] = a;
    m->sm[1] = b;
    m->sm[2] = c;
    m->sm[3] = d;
}
#endif
#else
void matrixMultiplyAligned(mtx *m, const mtx *m1, const mtx *m2)
{
#if PLATFORM_CELL
#if 0
    VmathMatrix4 sm __attribute__((aligned(16)));
    VmathMatrix4 sm1 __attribute__((aligned(16)));
    VmathMatrix4 sm2 __attribute__((aligned(16)));
    
    sm1 = sceM4MakeFromCols(sceV4MakeFromElems(m1->f[0][0], m1->f[0][1], m1->f[0][2], m1->f[0][3]),
			    sceV4MakeFromElems(m1->f[1][0], m1->f[1][1], m1->f[1][2], m1->f[1][3]),
			    sceV4MakeFromElems(m1->f[2][0], m1->f[2][1], m1->f[2][2], m1->f[2][3]),
			    sceV4MakeFromElems(m1->f[3][0], m1->f[3][1], m1->f[3][2], m1->f[3][3]));

    sm2 = sceM4MakeFromCols(sceV4MakeFromElems(m2->f[0][0], m2->f[0][1], m2->f[0][2], m2->f[0][3]),
			    sceV4MakeFromElems(m2->f[1][0], m2->f[1][1], m2->f[1][2], m2->f[1][3]),
			    sceV4MakeFromElems(m2->f[2][0], m2->f[2][1], m2->f[2][2], m2->f[2][3]),
			    sceV4MakeFromElems(m2->f[3][0], m2->f[3][1], m2->f[3][2], m2->f[3][3]));

    sm = sceM4Mul(sm2, sm1);

    m->f[0][0] = sceM4GetElem(sm, 0, 0); m->f[0][1] = sceM4GetElem(sm, 0, 1); m->f[0][2] = sceM4GetElem(sm, 0, 2); m->f[0][3] = sceM4GetElem(sm, 0, 3);
    m->f[1][0] = sceM4GetElem(sm, 1, 0); m->f[1][1] = sceM4GetElem(sm, 1, 1); m->f[1][2] = sceM4GetElem(sm, 1, 2); m->f[1][3] = sceM4GetElem(sm, 1, 3);
    m->f[2][0] = sceM4GetElem(sm, 2, 0); m->f[2][1] = sceM4GetElem(sm, 2, 1); m->f[2][2] = sceM4GetElem(sm, 2, 2); m->f[2][3] = sceM4GetElem(sm, 2, 3);
    m->f[3][0] = sceM4GetElem(sm, 3, 0); m->f[3][1] = sceM4GetElem(sm, 3, 1); m->f[3][2] = sceM4GetElem(sm, 3, 2); m->f[3][3] = sceM4GetElem(sm, 3, 3);
#else

    VmathMatrix4 sm1 __attribute__((aligned(16)));
    VmathMatrix4 sm2 __attribute__((aligned(16)));
    
    sm1 = m1->csm;
    sm2 = m2->csm;

    m->csm = vmathM4Mul_V(sm2, sm1);
    //m->sm = sceM4Mul(m1->sm, m2->sm);
//    assert((((u64)(&m->sm)) & 0xf) == 0);
//    assert((((u64)(&m1->sm)) & 0xf) == 0);
//    assert((((u64)(&m2->sm)) & 0xf) == 0);
//    assert((((u64)(&sm1)) & 0xf) == 0);
//    assert((((u64)(&sm2)) & 0xf) == 0);
#endif
#else
    mtx tmp;
    mtx *cm = m;
    boolean copy = FALSE;
    
    if ((m == m1) || (m == m2))
    {
	copy = TRUE;
	cm = &tmp;
    }

    cm->f[0][0] = m1->f[0][0]*m2->f[0][0] + m1->f[0][1]*m2->f[1][0] + m1->f[0][2]*m2->f[2][0] + m1->f[0][3]*m2->f[3][0];
    cm->f[0][1] = m1->f[0][0]*m2->f[0][1] + m1->f[0][1]*m2->f[1][1] + m1->f[0][2]*m2->f[2][1] + m1->f[0][3]*m2->f[3][1];
    cm->f[0][2] = m1->f[0][0]*m2->f[0][2] + m1->f[0][1]*m2->f[1][2] + m1->f[0][2]*m2->f[2][2] + m1->f[0][3]*m2->f[3][2];
    cm->f[0][3] = m1->f[0][0]*m2->f[0][3] + m1->f[0][1]*m2->f[1][3] + m1->f[0][2]*m2->f[2][3] + m1->f[0][3]*m2->f[3][3];

    cm->f[1][0] = m1->f[1][0]*m2->f[0][0] + m1->f[1][1]*m2->f[1][0] + m1->f[1][2]*m2->f[2][0] + m1->f[1][3]*m2->f[3][0];
    cm->f[1][1] = m1->f[1][0]*m2->f[0][1] + m1->f[1][1]*m2->f[1][1] + m1->f[1][2]*m2->f[2][1] + m1->f[1][3]*m2->f[3][1];
    cm->f[1][2] = m1->f[1][0]*m2->f[0][2] + m1->f[1][1]*m2->f[1][2] + m1->f[1][2]*m2->f[2][2] + m1->f[1][3]*m2->f[3][2];
    cm->f[1][3] = m1->f[1][0]*m2->f[0][3] + m1->f[1][1]*m2->f[1][3] + m1->f[1][2]*m2->f[2][3] + m1->f[1][3]*m2->f[3][3];

    cm->f[2][0] = m1->f[2][0]*m2->f[0][0] + m1->f[2][1]*m2->f[1][0] + m1->f[2][2]*m2->f[2][0] + m1->f[2][3]*m2->f[3][0];
    cm->f[2][1] = m1->f[2][0]*m2->f[0][1] + m1->f[2][1]*m2->f[1][1] + m1->f[2][2]*m2->f[2][1] + m1->f[2][3]*m2->f[3][1];
    cm->f[2][2] = m1->f[2][0]*m2->f[0][2] + m1->f[2][1]*m2->f[1][2] + m1->f[2][2]*m2->f[2][2] + m1->f[2][3]*m2->f[3][2];
    cm->f[2][3] = m1->f[2][0]*m2->f[0][3] + m1->f[2][1]*m2->f[1][3] + m1->f[2][2]*m2->f[2][3] + m1->f[2][3]*m2->f[3][3];

    cm->f[3][0] = m1->f[3][0]*m2->f[0][0] + m1->f[3][1]*m2->f[1][0] + m1->f[3][2]*m2->f[2][0] + m1->f[3][3]*m2->f[3][0];
    cm->f[3][1] = m1->f[3][0]*m2->f[0][1] + m1->f[3][1]*m2->f[1][1] + m1->f[3][2]*m2->f[2][1] + m1->f[3][3]*m2->f[3][1];
    cm->f[3][2] = m1->f[3][0]*m2->f[0][2] + m1->f[3][1]*m2->f[1][2] + m1->f[3][2]*m2->f[2][2] + m1->f[3][3]*m2->f[3][2];
    cm->f[3][3] = m1->f[3][0]*m2->f[0][3] + m1->f[3][1]*m2->f[1][3] + m1->f[3][2]*m2->f[2][3] + m1->f[3][3]*m2->f[3][3];

    if (copy)
    {
	matrixCopy(m, cm);
    }
#endif
}
#endif


#if !(PLATFORM_PC || PLATFORM_LINUX)
void matrixMultiply(mtx *m, const mtx *m1, const mtx *m2)
{
#if PLATFORM_CELL
#if 0
    VmathMatrix4 sm __attribute__((aligned(16)));
    VmathMatrix4 sm1 __attribute__((aligned(16)));
    VmathMatrix4 sm2 __attribute__((aligned(16)));
    
    sm1 = sceM4MakeFromCols(sceV4MakeFromElems(m1->f[0][0], m1->f[0][1], m1->f[0][2], m1->f[0][3]),
			    sceV4MakeFromElems(m1->f[1][0], m1->f[1][1], m1->f[1][2], m1->f[1][3]),
			    sceV4MakeFromElems(m1->f[2][0], m1->f[2][1], m1->f[2][2], m1->f[2][3]),
			    sceV4MakeFromElems(m1->f[3][0], m1->f[3][1], m1->f[3][2], m1->f[3][3]));

    sm2 = sceM4MakeFromCols(sceV4MakeFromElems(m2->f[0][0], m2->f[0][1], m2->f[0][2], m2->f[0][3]),
			    sceV4MakeFromElems(m2->f[1][0], m2->f[1][1], m2->f[1][2], m2->f[1][3]),
			    sceV4MakeFromElems(m2->f[2][0], m2->f[2][1], m2->f[2][2], m2->f[2][3]),
			    sceV4MakeFromElems(m2->f[3][0], m2->f[3][1], m2->f[3][2], m2->f[3][3]));

    sm = sceM4Mul(sm2, sm1);

    m->f[0][0] = sceM4GetElem(sm, 0, 0); m->f[0][1] = sceM4GetElem(sm, 0, 1); m->f[0][2] = sceM4GetElem(sm, 0, 2); m->f[0][3] = sceM4GetElem(sm, 0, 3);
    m->f[1][0] = sceM4GetElem(sm, 1, 0); m->f[1][1] = sceM4GetElem(sm, 1, 1); m->f[1][2] = sceM4GetElem(sm, 1, 2); m->f[1][3] = sceM4GetElem(sm, 1, 3);
    m->f[2][0] = sceM4GetElem(sm, 2, 0); m->f[2][1] = sceM4GetElem(sm, 2, 1); m->f[2][2] = sceM4GetElem(sm, 2, 2); m->f[2][3] = sceM4GetElem(sm, 2, 3);
    m->f[3][0] = sceM4GetElem(sm, 3, 0); m->f[3][1] = sceM4GetElem(sm, 3, 1); m->f[3][2] = sceM4GetElem(sm, 3, 2); m->f[3][3] = sceM4GetElem(sm, 3, 3);
#else

    VmathMatrix4 sm1 __attribute__((aligned(16)));
    VmathMatrix4 sm2 __attribute__((aligned(16)));
    
    sm1 = m1->csm;
    sm2 = m2->csm;

    m->csm = vmathM4Mul_V(sm2, sm1);
    //m->sm = sceM4Mul(m1->sm, m2->sm);
//    assert((((u64)(&m->sm)) & 0xf) == 0);
//    assert((((u64)(&m1->sm)) & 0xf) == 0);
//    assert((((u64)(&m2->sm)) & 0xf) == 0);
//    assert((((u64)(&sm1)) & 0xf) == 0);
//    assert((((u64)(&sm2)) & 0xf) == 0);
#endif
#else

    mtx tmp;
    mtx *cm = m;
    boolean copy = FALSE;
    
    if ((m == m1) || (m == m2))
    {
	copy = TRUE;
	cm = &tmp;
    }

    cm->f[0][0] = m1->f[0][0]*m2->f[0][0] + m1->f[0][1]*m2->f[1][0] + m1->f[0][2]*m2->f[2][0] + m1->f[0][3]*m2->f[3][0];
    cm->f[0][1] = m1->f[0][0]*m2->f[0][1] + m1->f[0][1]*m2->f[1][1] + m1->f[0][2]*m2->f[2][1] + m1->f[0][3]*m2->f[3][1];
    cm->f[0][2] = m1->f[0][0]*m2->f[0][2] + m1->f[0][1]*m2->f[1][2] + m1->f[0][2]*m2->f[2][2] + m1->f[0][3]*m2->f[3][2];
    cm->f[0][3] = m1->f[0][0]*m2->f[0][3] + m1->f[0][1]*m2->f[1][3] + m1->f[0][2]*m2->f[2][3] + m1->f[0][3]*m2->f[3][3];

    cm->f[1][0] = m1->f[1][0]*m2->f[0][0] + m1->f[1][1]*m2->f[1][0] + m1->f[1][2]*m2->f[2][0] + m1->f[1][3]*m2->f[3][0];
    cm->f[1][1] = m1->f[1][0]*m2->f[0][1] + m1->f[1][1]*m2->f[1][1] + m1->f[1][2]*m2->f[2][1] + m1->f[1][3]*m2->f[3][1];
    cm->f[1][2] = m1->f[1][0]*m2->f[0][2] + m1->f[1][1]*m2->f[1][2] + m1->f[1][2]*m2->f[2][2] + m1->f[1][3]*m2->f[3][2];
    cm->f[1][3] = m1->f[1][0]*m2->f[0][3] + m1->f[1][1]*m2->f[1][3] + m1->f[1][2]*m2->f[2][3] + m1->f[1][3]*m2->f[3][3];

    cm->f[2][0] = m1->f[2][0]*m2->f[0][0] + m1->f[2][1]*m2->f[1][0] + m1->f[2][2]*m2->f[2][0] + m1->f[2][3]*m2->f[3][0];
    cm->f[2][1] = m1->f[2][0]*m2->f[0][1] + m1->f[2][1]*m2->f[1][1] + m1->f[2][2]*m2->f[2][1] + m1->f[2][3]*m2->f[3][1];
    cm->f[2][2] = m1->f[2][0]*m2->f[0][2] + m1->f[2][1]*m2->f[1][2] + m1->f[2][2]*m2->f[2][2] + m1->f[2][3]*m2->f[3][2];
    cm->f[2][3] = m1->f[2][0]*m2->f[0][3] + m1->f[2][1]*m2->f[1][3] + m1->f[2][2]*m2->f[2][3] + m1->f[2][3]*m2->f[3][3];

    cm->f[3][0] = m1->f[3][0]*m2->f[0][0] + m1->f[3][1]*m2->f[1][0] + m1->f[3][2]*m2->f[2][0] + m1->f[3][3]*m2->f[3][0];
    cm->f[3][1] = m1->f[3][0]*m2->f[0][1] + m1->f[3][1]*m2->f[1][1] + m1->f[3][2]*m2->f[2][1] + m1->f[3][3]*m2->f[3][1];
    cm->f[3][2] = m1->f[3][0]*m2->f[0][2] + m1->f[3][1]*m2->f[1][2] + m1->f[3][2]*m2->f[2][2] + m1->f[3][3]*m2->f[3][2];
    cm->f[3][3] = m1->f[3][0]*m2->f[0][3] + m1->f[3][1]*m2->f[1][3] + m1->f[3][2]*m2->f[2][3] + m1->f[3][3]*m2->f[3][3];

    if (copy)
    {
	matrixCopy(m, cm);
    }
#endif
}
#endif

void matrix33Multiply33(mtx *m, const mtx *m1, const mtx *m2)
{
    mtx tmp;
    mtx *cm = m;
    boolean copy = FALSE;
    
    if ((m == m1) || (m == m2))
    {
	copy = TRUE;
	cm = &tmp;
    }

    cm->f[0][0] = m1->f[0][0]*m2->f[0][0] + m1->f[0][1]*m2->f[1][0] + m1->f[0][2]*m2->f[2][0];
    cm->f[0][1] = m1->f[0][0]*m2->f[0][1] + m1->f[0][1]*m2->f[1][1] + m1->f[0][2]*m2->f[2][1];
    cm->f[0][2] = m1->f[0][0]*m2->f[0][2] + m1->f[0][1]*m2->f[1][2] + m1->f[0][2]*m2->f[2][2];
    cm->f[0][3] = 0.0f;

    cm->f[1][0] = m1->f[1][0]*m2->f[0][0] + m1->f[1][1]*m2->f[1][0] + m1->f[1][2]*m2->f[2][0];
    cm->f[1][1] = m1->f[1][0]*m2->f[0][1] + m1->f[1][1]*m2->f[1][1] + m1->f[1][2]*m2->f[2][1];
    cm->f[1][2] = m1->f[1][0]*m2->f[0][2] + m1->f[1][1]*m2->f[1][2] + m1->f[1][2]*m2->f[2][2];
    cm->f[1][3] = 0.0f;

    cm->f[2][0] = m1->f[2][0]*m2->f[0][0] + m1->f[2][1]*m2->f[1][0] + m1->f[2][2]*m2->f[2][0];
    cm->f[2][1] = m1->f[2][0]*m2->f[0][1] + m1->f[2][1]*m2->f[1][1] + m1->f[2][2]*m2->f[2][1];
    cm->f[2][2] = m1->f[2][0]*m2->f[0][2] + m1->f[2][1]*m2->f[1][2] + m1->f[2][2]*m2->f[2][2];
    cm->f[2][3] = 0.0f;

    cm->f[3][0] = 0.0f;
    cm->f[3][1] = 0.0f;
    cm->f[3][2] = 0.0f;
    cm->f[3][3] = 1.0f;

    if (copy)
    {
	matrixCopy(m, cm);
    }
}

void matrix33Multiply33Tran1(mtx *m, const mtx *m1, const mtx *m2)
{
    mtx tmp;
    mtx *cm = m;
    boolean copy = FALSE;
    
    if ((m == m1) || (m == m2))
    {
	copy = TRUE;
	cm = &tmp;
    }

    cm->f[0][0] = m1->f[0][0]*m2->f[0][0] + m1->f[1][0]*m2->f[1][0] + m1->f[2][0]*m2->f[2][0];
    cm->f[0][1] = m1->f[0][0]*m2->f[0][1] + m1->f[1][0]*m2->f[1][1] + m1->f[2][0]*m2->f[2][1];
    cm->f[0][2] = m1->f[0][0]*m2->f[0][2] + m1->f[1][0]*m2->f[1][2] + m1->f[2][0]*m2->f[2][2];
    cm->f[0][3] = 0.0f;

    cm->f[1][0] = m1->f[0][1]*m2->f[0][0] + m1->f[1][1]*m2->f[1][0] + m1->f[2][1]*m2->f[2][0];
    cm->f[1][1] = m1->f[0][1]*m2->f[0][1] + m1->f[1][1]*m2->f[1][1] + m1->f[2][1]*m2->f[2][1];
    cm->f[1][2] = m1->f[0][1]*m2->f[0][2] + m1->f[1][1]*m2->f[1][2] + m1->f[2][1]*m2->f[2][2];
    cm->f[1][3] = 0.0f;

    cm->f[2][0] = m1->f[0][2]*m2->f[0][0] + m1->f[1][2]*m2->f[1][0] + m1->f[2][2]*m2->f[2][0];
    cm->f[2][1] = m1->f[0][2]*m2->f[0][1] + m1->f[1][2]*m2->f[1][1] + m1->f[2][2]*m2->f[2][1];
    cm->f[2][2] = m1->f[0][2]*m2->f[0][2] + m1->f[1][2]*m2->f[1][2] + m1->f[2][2]*m2->f[2][2];
    cm->f[2][3] = 0.0f;

    cm->f[3][0] = 0.0f;
    cm->f[3][1] = 0.0f;
    cm->f[3][2] = 0.0f;
    cm->f[3][3] = 1.0f;

    if (copy)
    {
	matrixCopy(m, cm);
    }
}

void matrixFrom2Axes(mtx *m, vec3 *v1, vec3 *v2)
{
    float l, k;
    vec3 a, b;

    l = 1.0f/vecsize(v1);
    vecscale(&a, v1, l);

    k = vecdot(&a, v2);
    vecsubscale(&b, v2, &a, k);

    l = 1.0f/vecsize(&b);
    vecscale(&b, &b, l);

    m->f[0][0] = a.x;
    m->f[0][1] = a.y;
    m->f[0][2] = a.z;
    m->f[0][3] = 0.0f;

    m->f[1][0] = b.x;
    m->f[1][1] = b.y;
    m->f[1][2] = b.z;
    m->f[1][3] = 0.0f;

    m->f[2][0] = a.y*b.z - b.y*a.z;
    m->f[2][1] = a.z*b.x - b.z*a.x;
    m->f[2][2] = a.x*b.y - b.x*a.y;
    m->f[2][3] = 0.0f;

    m->f[3][0] = 0.0f;
    m->f[3][1] = 0.0f;
    m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}


void matrixPrintfTranspose(const mtx *m)
{
    printf("%f %f %f %f\n", m->f[0][0], m->f[1][0], m->f[2][0], m->f[3][0]);
    printf("%f %f %f %f\n", m->f[0][1], m->f[1][1], m->f[2][1], m->f[3][1]);
    printf("%f %f %f %f\n", m->f[0][2], m->f[1][2], m->f[2][2], m->f[3][2]);
    printf("%f %f %f %f\n", m->f[0][3], m->f[1][3], m->f[2][3], m->f[3][3]);
}

void matrixPrintf(const mtx *m)
{
    printf("%f %f %f %f\n", m->f[0][0], m->f[0][1], m->f[0][2], m->f[0][3]);
    printf("%f %f %f %f\n", m->f[1][0], m->f[1][1], m->f[1][2], m->f[1][3]);
    printf("%f %f %f %f\n", m->f[2][0], m->f[2][1], m->f[2][2], m->f[2][3]);
    printf("%f %f %f %f\n", m->f[3][0], m->f[3][1], m->f[3][2], m->f[3][3]);
}

void matrixToRotationZXY(const mtx *m, float *r1, float *r2, float *r3)
{
#if 0
    // This code was broken for the case of making a matrix out of 
    // rx = 90; ry = 0; rz = 0
    // and then extracting the rotations using this code would result in
    // rx = 360; ry = 0; rz = 0;
    
#define EPSILON	    (1e-5f) // should be fine accuracy

    float cx, rx, ry, rz;
    
    rx = atan2f( -m->f[2][1], sqrtf( SQR(m->f[2][0]) + SQR(m->f[2][2]) ) );
    cx = cosf(rx);
    
    if( fabsf(cx) > EPSILON )
    {
	cx = 1.0f/cx;
	rz = atan2f( m->f[0][1]*cx , m->f[1][1]*cx );
        ry = atan2f( m->f[2][0]*cx , m->f[2][2]*cx );
    }
    else
    {
	if( fabsf(PI*2.0f - rx) < EPSILON )
	{
	    rx = -PI*2.0f;
            ry = 0.0f;
            rz = atan2f( -m->f[0][2], m->f[0][0] );
	}
        else
	{
	    rx = PI*2.0f;
            ry = 0.0f;
            rz = atan2f( m->f[0][2], m->f[0][0] );
        }
    }

    if( r1 ) *r1=rx;
    if( r2 ) *r2=ry;
    if( r3 ) *r3=rz;
#undef EPSILON
#else
    float rx, ry, rz;
    
    // this is code for factoring as Ry Rx Rz
    // JH: But still doesn't work - generates inconsistent angles for continuously rotating matrices
    rx = asinf(-m->f[2][1]);
    if (rx < PI_D2)
    {
	if (rx > -PI_D2)
	{
	    ry = atan2f(m->f[2][0], m->f[2][2]);
	    rz = atan2f(m->f[0][1], m->f[1][1]);
	}
	else
	{
	    // not a unique solution (ry+rz constant)
	    ry = atan2f(-m->f[1][0], m->f[0][0]);
	    rz = 0.0f;
	}
    }
    else
    {
	// not a unique solution (ry-rz constant)
	ry = atan2f(m->f[1][0], m->f[0][0]);
	rz = 0.0f;
    }

    if( r1 ) *r1=rx;
    if( r2 ) *r2=ry;
    if( r3 ) *r3=rz;
#endif
}

// Rotate a vector with the rotation components of a matrix. (from SS).
void matrixVecRot(const mtx *m, vec3 *vec)
{
    const float x = vec->x;
    const float y = vec->y;
    const float z = vec->z;
    const float xOut = x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0];
    const float yOut = x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1];
    const float zOut = x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2];

    vec->x = xOut;
    vec->y = yOut;
    vec->z = zOut;
}

void matrixVecRotCpy(vec3 *out, const mtx *m, const vec3 *vec)
{
    const float x = vec->x;
    const float y = vec->y;
    const float z = vec->z;

    out->x = x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0];
    out->y = x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1];
    out->z = x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2];
}

/*
** matrixVecRotTranspose - Multiply a vector by the rotation components of the transpose of a matrix
*/
void matrixVecRotTranspose(const mtx *m, vec3 *v)
{
    const float x = v->v[0];
    const float y = v->v[1];
    const float z = v->v[2];

    v->v[0] = (x*m->f[0][0]) + (y*m->f[0][1]) + (z*m->f[0][2]);
    v->v[1] = (x*m->f[1][0]) + (y*m->f[1][1]) + (z*m->f[1][2]);
    v->v[2] = (x*m->f[2][0]) + (y*m->f[2][1]) + (z*m->f[2][2]);
}

void matrixVecMulCpy(vec3 *outVec, const mtx *m, const vec3 *v)
{
    const float x = v->x;
    const float y = v->y;
    const float z = v->z;

    vecset(outVec, x*m->f[0][0] + y*m->f[1][0] + z*m->f[2][0] + m->f[3][0],
		   x*m->f[0][1] + y*m->f[1][1] + z*m->f[2][1] + m->f[3][1],
		   x*m->f[0][2] + y*m->f[1][2] + z*m->f[2][2] + m->f[3][2]);
}

void matrixVecMul(const mtx *m, vec3 *vec)
{
    matrixVecMulCpy(vec, m, vec);
}

/*
** matrixVecMulTranspose - Multiply a vector by the transpose of a matrix
*/
void matrixVecMulTranspose(const mtx *m, vec3 *v)
{
    const float x = v->v[0]-m->f[3][0];
    const float y = v->v[1]-m->f[3][1];
    const float z = v->v[2]-m->f[3][2];

    v->v[0] = (x*m->f[0][0]) + (y*m->f[0][1]) + (z*m->f[0][2]);
    v->v[1] = (x*m->f[1][0]) + (y*m->f[1][1]) + (z*m->f[1][2]);
    v->v[2] = (x*m->f[2][0]) + (y*m->f[2][1]) + (z*m->f[2][2]);
}

/*
** matrixVec3Mul4 - Multiply a vector by a matrix
*/
void matrixVec3Mul4(const mtx *m, const vec3 *vec, vec4 *newvec)
{
    const float x=vec->x;
    const float y=vec->y;
    const float z=vec->z;
    
    newvec->x = (x*m->f[0][0]) + (y*m->f[1][0]) + (z*m->f[2][0]) + m->f[3][0];
    newvec->y = (x*m->f[0][1]) + (y*m->f[1][1]) + (z*m->f[2][1]) + m->f[3][1];
    newvec->z = (x*m->f[0][2]) + (y*m->f[1][2]) + (z*m->f[2][2]) + m->f[3][2];
    newvec->w = (x*m->f[0][3]) + (y*m->f[1][3]) + (z*m->f[2][3]) + m->f[3][3];
}

// zaxis passed in should be normalised
// creates an orthonormal basis with the given zaxis with the yaxis in the vertical plane
void matrixFromZAxis( mtx *m, const vec3 *zaxis )
{
    matrixFromZAxisWithPos(m,zaxis,NULL);
}
    
// as above but you can explictly specify the position that should be set for the mtx
// zaxis should be normalised
void matrixFromZAxisWithPos( mtx *m, const vec3 *zaxis, const vec3 *inPos )
{
    vec3    right;
    vec3    up;
    vec3    up2;
    
    vecset( &up2, 0.f, 1.f, 0.f );
    if ( fabsf(vecdot(&up2, zaxis)) > 0.99f )
    {
	vecset( &up2, -1.f, 0.f, 0.f );
    }
    veccross(&right, &up2, zaxis);
    vecnormalise( &right, &right );
    veccross(&up, zaxis, &right);
//    vecnormalise( &up, &up );		    // technically not necessary as crossing two normalised vecs will result in another normalised vec
    assert(vecisunit(&up));		    // asserted here - also will trip if zaxis wasn't normalised already
    
    vecset4( &m->v[0], right.x, right.y, right.z, 0.f );
    vecset4( &m->v[1], up.x, up.y, up.z, 0.f );
    vecset4( &m->v[2], zaxis->x, zaxis->y, zaxis->z, 0.f );
    if (inPos)
    {
	vecset4( &m->v[3], vec3listp(inPos), 1.f );
    }
    else
    {
	vecset4( &m->v[3], 0.f, 0.f, 0.f, 1.f );
    }

    assert( matrixGetDeterminant( m ) > 0.f );
}

// yaxis passed in should be normalised
// creates an orthonormal basis with the given yaxis 
void matrixFromYAxis( mtx *m, const vec3 *yaxis )
{
    matrixFromYAxisWithPos(m,yaxis,NULL);
}
    
// as above but you can explictly specify the position that should be set for the mtx
// yaxis should be normalised
void matrixFromYAxisWithPos( mtx *m, const vec3 *yaxis, const vec3 *inPos )
{
    vec3    right;
    vec3    right2;
    vec3    forward;
    
    vecset( &right2, 1.f, 0.f, 0.f );
    if ( fabsf(vecdot(&right2, yaxis)) > 0.99f )
    {
	vecset( &right2, 0.f, 1.f, 0.f );
    }
    veccross(&forward, &right2, yaxis);
    vecnormalise( &forward, &forward );	    // not neccessary either surely!!
    veccross(&right, yaxis, &forward);
//    vecnormalise( &right, &right );	    // technically not necessary as crossing two normalised vecs will result in another normalised vec
    assert(vecisunit(&right));		    // asserted here - also will trip if zaxis wasn't normalised already
    
    vecset4( &m->v[0], right.x, right.y, right.z, 0.f );
    vecset4( &m->v[1], yaxis->x, yaxis->y, yaxis->z, 0.f );
    vecset4( &m->v[2], forward.x, forward.y, forward.z, 0.f );
    if (inPos)
    {
	vecset4( &m->v[3], vec3listp(inPos), 1.f );
    }
    else
    {
	vecset4( &m->v[3], 0.f, 0.f, 0.f, 1.f );
    }

    assert( matrixGetDeterminant( m ) > 0.f );
}


void matrixFromXAxis( mtx *m, const vec3 *xaxis )
{
    vec3    right;
    vec3    up;
    vec3    up2;
    
    vecset( &up2, 0.f, 1.f, 0.f );
    if ( fabsf(vecdot(&up2, xaxis)) > 0.99f )
    {
	vecset( &up2, 1.f, 0.f, 0.f );
    }
    veccross(&right, &up2, xaxis);
    vecnormalise( &right, &right );
    veccross(&up, &right, xaxis);
    vecnormalise( &up, &up );
    
    vecset4( &m->v[0], xaxis->x, xaxis->y, xaxis->z, 0.f );
    vecset4( &m->v[1], up.x, up.y, up.z, 0.f );
    vecset4( &m->v[2], right.x, right.y, right.z, 0.f );
    vecset4( &m->v[3], 0.f, 0.f, 0.f, 1.f );
    
    assert( matrixGetDeterminant( m ) > 0.f );
}

// set a matrix from 3 unit vectors
void matrixFrom3Axis(mtx *outMtx, const vec3 *inX, const vec3 *inY, const vec3 *inZ)
{
    veccpy((vec3*)&outMtx->v[0],inX);	    outMtx->v[0].w=0.0f;
    veccpy((vec3*)&outMtx->v[1],inY);	    outMtx->v[1].w=0.0f;
    veccpy((vec3*)&outMtx->v[2],inZ);	    outMtx->v[2].w=0.0f;
    veczero((vec3*)&outMtx->v[3]);	    outMtx->v[3].w=1.0f;
}

void matrixroll(union mtx_u *newmtx, vec3 *axis, float axialrot)
{
    float s,c,t,x,y,z;

    // Rotate about axis.
    s=sinf(axialrot);
    c=cosf(axialrot);
    t=1.0f-c;

    x=axis->x;
    y=axis->y;
    z=axis->z;

    newmtx->f[0][0]=(t*(x*x))+c;
    newmtx->f[0][1]=(t*x*y)-(s*z);
    newmtx->f[0][2]=(t*x*z)+(s*y);
    newmtx->f[1][0]=(t*x*y)+(s*z);
    newmtx->f[1][1]=(t*(y*y))+c;
    newmtx->f[1][2]=(t*y*z)-(s*x);
    newmtx->f[2][0]=(t*x*z)-(s*y);
    newmtx->f[2][1]=(t*y*z)+(s*x);
    newmtx->f[2][2]=(t*(z*z))+c;
    newmtx->f[3][0]=newmtx->f[3][1]=newmtx->f[3][2]=newmtx->f[0][3]=newmtx->f[1][3]=newmtx->f[2][3]=0.0f;
    newmtx->f[3][3]=1.0f;
}

void matrixReNormalise(mtx* m)
{
    vec3* x = &m->v[0].v3;
    vec3* y = &m->v[1].v3;
    vec3* z = &m->v[2].v3;
    
    vecnormalise(x, x);
    veccross(z, x, y);
    vecnormalise(z, z);
    veccross(y, z, x);
}

void quaternionRotateVector(vec3 *v, const quaternion *q1, const vec3 *v1)
{
    float x, y, z, w;

    // Calculate quat * vector * conjugate of quat

    x = (v1->x*q1->w) - (v1->y*q1->z) + (v1->z*q1->y);
    y = (v1->y*q1->w) - (v1->z*q1->x) + (v1->x*q1->z);
    z = (v1->z*q1->w) - (v1->x*q1->y) + (v1->y*q1->x);
    w = (v1->x*q1->x) + (v1->y*q1->y) + (v1->z*q1->z);

    v->x = (q1->w*x) + (q1->x*w) + (q1->y*z) - (q1->z*y);
    v->y = (q1->w*y) + (q1->y*w) + (q1->z*x) - (q1->x*z);
    v->z = (q1->w*z) + (q1->z*w) + (q1->x*y) - (q1->y*x);
}


void quaternionSet(quaternion *q, float x, float y, float z, float w)
{
    q->x = x;
    q->y = y;
    q->z = z;
    q->w = w;
}

void quaternionIdent(quaternion *q)
{
    q->x = 0.0f;
    q->y = 0.0f;
    q->z = 0.0f;
    q->w = 1.0f;
}

void quaternionCopy(quaternion *q, const quaternion *q1)
{
    q->x = q1->x;
    q->y = q1->y;
    q->z = q1->z;
    q->w = q1->w;
}

void quaternionNegate(quaternion *q, const quaternion *q1)
{
    q->x = -q1->x;
    q->y = -q1->y;
    q->z = -q1->z;
    q->w = -q1->w;
}

void quaternionAdd(quaternion *q, const quaternion *q1, const quaternion *q2)
{
    q->x = q1->x + q2->x;
    q->y = q1->y + q2->y;
    q->z = q1->z + q2->z;
    q->w = q1->w + q2->w;
}

void quaternionSub(quaternion *q, const quaternion *q1, const quaternion *q2)
{
    q->x = q1->x - q2->x;
    q->y = q1->y - q2->y;
    q->z = q1->z - q2->z;
    q->w = q1->w - q2->w;
}

void quaternionToMatrix(mtx* m, const quaternion* q)
{
    float wx,wy,wz,xx,xy,xz,yy,yz,zz;

    float xs = q->x + q->x;
    float ys = q->y + q->y;
    float zs = q->z + q->z;

    wx=q->w*xs; wy=q->w*ys; wz=q->w*zs;
    xx=q->x*xs; xy=q->x*ys; xz=q->x*zs;
    yy=q->y*ys; yz=q->y*zs; zz=q->z*zs;

    m->f[0][0]=1.0f-(yy+zz);
    m->f[0][1]=xy+wz;
    m->f[0][2]=xz-wy;

    m->f[1][0]=xy-wz;
    m->f[1][1]=1.0f-(xx+zz);
    m->f[1][2]=yz+wx;

    m->f[2][0]=xz+wy;
    m->f[2][1]=yz-wx;
    m->f[2][2]=1.0f-(xx+yy);

    m->f[0][3] = m->f[1][3] = m->f[2][3] = 0.0f;
    m->f[3][0] = m->f[3][1] = m->f[3][2] = 0.0f;
    m->f[3][3] = 1.0f;
}

void quaternionGetXBasis(vec3* out, const quaternion* q)
{
    float wx,wy,wz,xx,xy,xz,yy,yz,zz;

    float xs = q->x + q->x;
    float ys = q->y + q->y;
    float zs = q->z + q->z;

    wx=q->w*xs; wy=q->w*ys; wz=q->w*zs;
    xx=q->x*xs; xy=q->x*ys; xz=q->x*zs;
    yy=q->y*ys; yz=q->y*zs; zz=q->z*zs;

    out->x=1.0f-(yy+zz);
    out->y=xy+wz;
    out->z=xz-wy;
}

void quaternionGetYBasis(vec3* out, const quaternion* q)
{
    float wx,wy,wz,xx,xy,xz,yy,yz,zz;

    float xs = q->x + q->x;
    float ys = q->y + q->y;
    float zs = q->z + q->z;

    wx=q->w*xs; wy=q->w*ys; wz=q->w*zs;
    xx=q->x*xs; xy=q->x*ys; xz=q->x*zs;
    yy=q->y*ys; yz=q->y*zs; zz=q->z*zs;

    out->x=xy-wz;
    out->y=1.0f-(xx+zz);
    out->z=yz+wx;

}

void quaternionGetZBasis(vec3* out, const quaternion* q)
{
    float wx,wy,wz,xx,xy,xz,yy,yz,zz;

    float xs = q->x + q->x;
    float ys = q->y + q->y;
    float zs = q->z + q->z;

    wx=q->w*xs; wy=q->w*ys; wz=q->w*zs;
    xx=q->x*xs; xy=q->x*ys; xz=q->x*zs;
    yy=q->y*ys; yz=q->y*zs; zz=q->z*zs;

    out->x=xz+wy;
    out->y=yz-wx;
    out->z=1.0f-(xx+yy);
}


void quaternionNormalise(quaternion *q)
{
    float mag;
    mag=sqrtf(1.0f/((q->x*q->x)+(q->y*q->y)+(q->z*q->z)+(q->w*q->w)));
    q->x=q->x*mag;
    q->y=q->y*mag;
    q->z=q->z*mag;
    q->w=q->w*mag;
}



// A WORD TO THE WISE : quaternionMul tkaes it's args the other way around from matrixMultiply() - eg quaternionMul(out,a,b) == matrixMultiply(out,b,a). Somebody didn't think that one through. MT
// q1*q2 = (s1s2 - v1.v2, s1v2+s2v1+v1 x v2)
void quaternionMul(quaternion *q, const quaternion *q1, const quaternion *q2)
{
    float x = (q1->w*q2->x) + (q1->x*q2->w) + (q1->y*q2->z) - (q1->z*q2->y);
    float y = (q1->w*q2->y) + (q1->y*q2->w) + (q1->z*q2->x) - (q1->x*q2->z);
    float z = (q1->w*q2->z) + (q1->z*q2->w) + (q1->x*q2->y) - (q1->y*q2->x);
    float w = (q1->w*q2->w) - (q1->x*q2->x) - (q1->y*q2->y) - (q1->z*q2->z);
    q->x = x;
    q->y = y;
    q->z = z;
    q->w = w;
}
