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

#define fsel(test, a, b)    (((test)>=0.f) ? (a) : (b))
#define CLOTH_COLLISION_EPSILON 1e-6f

//#define print 0 &&
#define print printf

bool smallerQuadraticSolution(float a, float b, float c, float tMin, float tMax, float &t)
{
    if (fabsf(a) < CLOTH_COLLISION_EPSILON)
    {
	if (b == 0.0f)
	{
	    if (c != 0.0f) return false;
	    else t = 0.0f;
	}
	else 
	{
	    t = -c/b;
	}
    }
    else
    {
	b *= 0.5f;
	float d = b*b - a*c;
	if (d < 0.0f) return false;
	d = sqrtf(d);
	if (a < 0.0f) { a = -a; b = -b; }	// make a positive
	a = 1.0f / a;
	t = (-b-d)*a;	// this is the smaller solution -- only if a is positive!!
	if (t > tMax) return false;
	if (t < tMin) 
	{
	    t = (-b+d)*a;
	}
    }
    return (tMin <= t && t <= tMax);
}

bool cpSolveQuadForSphere(float a, float b, float c, float *t) // assumes a is positive and only positive t is valid
{
    if (b>=0.f)
    {
	return false;
    }
    else
    {
	float r;
	float inner;
	b = b * -0.5f;
	inner=SQR(b)-a*c;
	if (inner<0.0f)
	{
	    return false;
	}
	// inner will never be larger than fabs(b), therefore no further checks needed
	inner=sqrtf(inner);
	*t=(b-inner)/a;
	return true;
    }
}

/* dir, start */
bool cpSolveSphereLineTest(vec3 *a, vec3 *b, float r, float *t, int *inside)
{
    float A,B,C;
    C=vecsizesq(b)-SQR(r);
    if (C<0.0f) // already inside (flipin' 'eck!)
    {
	*t=0.0f;
	*inside=TRUE;
	print("ret true (4) ");
	return true;
    }
    A=vecsizesq(a);
    B=2.0f*vecdot(a,b);
    if (A<CLOTH_COLLISION_EPSILON) // hardly moving so ignore the At2 part
    {
	*t=fsel(fabsf(C),-C/B,0.0f); // OK to be NaN as below comparison should fail
	print("ret blah ");
	return (*t>=0.0f && *t<=1.0f);
    }
    return cpSolveQuadForSphere(A,B,C,t) && *t<=1.0f;
}


int main (int argc, char **argv)
{
#if 0
    float radius = 1.f;
    vec3 dir = {3.0f, 0.f, 0.f};
    float y = radius + 1.f;
    float dy = 0.05001f;
    for ( ; y >= 0.f; y -= dy)
    {
	int inside;
	float t = -1.f;
	vec3 pos = {0.8f, y, 0.0f};
	bool hit = cpSolveSphereLineTest(&dir, &pos, radius, &t, &inside);
	printf("hit = %d, y = %f, t = %f\n", hit, pos.y, t);
    }
#endif
    u32 skip = 0;
    u32 i = 3;
    u32 discrete = 1;
    u32 d = discrete<<i;
    bool add, add2;

#define DOIT(_i, _discrete) \
    discrete = _discrete; \
    i = _i;\
    d = discrete<<i;\
    add = (skip&1<<i)==0 || discrete;\
    add2 = (~skip | d) & (1<<i) ; \
    printf("idx = %d, should skip = %d, new_collision_is_discrete = %d, can add = %d, %d\n", i, bool(skip & 1<<i), discrete, add, add2);\
    skip |= d;\

    DOIT(3, 1);
    DOIT(3, 0);
    DOIT(3, 1);
    DOIT(3, 0);
    

}


