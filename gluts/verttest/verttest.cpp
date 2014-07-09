#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine/engine.h"

void vecRotateAboutVec(vec3* v, const vec3* axis, float angle)
{
    float axisSizeSq = vecsizesq(axis);
    if (axisSizeSq < 1e-6f) return;
    
    float dot=vecdot(v, axis) / axisSizeSq;
    vec3 a, c;
    vecscale(&c, axis, dot);
    vecsub(&a, v, &c);
    vec3 b;
    veccross(&b, axis, &a);
    vecscale(&b, &b, 1.0f/sqrtf(axisSizeSq)); // crap need to normalise

    printf("a = %f %f %f (%f)\n", vec3list(a), vecsize(&a));
    printf("b = %f %f %f (%f)\n", vec3list(b), vecsize(&b));
    printf("c = %f %f %f (%f)\n", vec3list(c), vecsize(&c));

    // v = a*cos + b*sin + c
    
    vec3 tmp1;
    vecscale(&tmp1, &a, cosf(angle));
    
    vec3 tmp2;
    vecscale(&tmp2, &b, sinf(angle));

    vecadd(v, &c, &tmp1);
    vecadd(v, v, &tmp2);
    
    printf("v = %f %f %f (%f)\n", vec3listp(v), vecsize(v));
}
    
void printvec(const char* msg, const vec3* v)
{
    printf("%s %f %f %f (%f)\n", msg, vec3listp(v), vecsize(v));
}
#define PRIINTVEC(x) printvec(#x, &x)
#define PRIINTVECP(x) printvec(#x, x)

void vecRotateAboutVec2(vec3* v, const vec3* axis, float angle)
{
    
    float axisSizeSq = vecsizesq(axis);
    if (axisSizeSq < 1e-6f) return;
    
    float dot=vecdot(v, axis) / axisSizeSq;
    vec3 a, c;
    vecscale(&c, axis, dot);
    vecsub(&a, v, &c);

    vec3 b;
    veccross(&b, axis, v);
    printvec("b", &b);

    vec3 ab;
    veccross(&ab, &a, &b);
    PRIINTVEC(ab);

    float s = sqrtf(axisSizeSq);
    printf("(a.a)*|axis| = %f\n", vecsizesq(&a)*s);
    vecscale(&ab, &ab, 1.f/vecdot(v,axis) );
    PRIINTVEC(ab);
}

void vecRotateAboutVec3(vec3* p, const vec3* axis, float angle)
{
    float x = p->x;
    float y = p->y;
    float z = p->z;
    
    float u = axis->x;
    float v = axis->y;
    float w = axis->z;
    
    float ux=u*x;
    float uy=u*y;
    float uz=u*z;
    float vx=v*x;
    float vy=v*y;
    float vz=v*z;
    float wx=w*x;
    float wy=w*y;
    float wz=w*z;
    
    float sa=sinf(angle);
    float ca=cosf(angle);
    
    p->x=u*(ux+vy+wz)+(x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa;
    p->y=v*(ux+vy+wz)+(y*(u*u+w*w)-v*(ux+wz))*ca+(wx-uz)*sa;
    p->z=w*(ux+vy+wz)+(z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa;
}

struct linedef
{
    vec3 point;
    vec3 dir;
};


float distanceLineLineSq(const linedef *line1, const linedef *line2, float *t1, float *t2)
{
    vec3 u;
    float a, b, c, d, e, f, det;
    float s, t;
    float distsq;

    vecsub(&u, &line1->point, &line2->point);
    a = vecdot(&line1->dir, &line1->dir);
    b = -vecdot(&line1->dir, &line2->dir);
    c = vecdot(&line2->dir, &line2->dir);
    d = vecdot(&line1->dir, &u);
    e = -vecdot(&line2->dir, &u);
    f = vecdot(&u, &u);
    det = fabsf(a*c - b*b);

    if (a==0.f)
    {
	s = 0.f;
	t = 0.f;
	distsq = FLOAT_MAX;
    }
    else if (det <= 0.0f)
    {
        // lines are parallel, select any closest pair of points
        s = -d/a;
        t = 0.0f;
        distsq = d*s + f;
    }
    else
    {
        // lines are not parallel
        s = (b*e - c*d)/det;
        t = (b*d - a*e)/det;
        distsq = s*(a*s + b*t + 2.0f*d) + t*(b*s + c*t + 2.0f*e) + f;
    }

    if (t1)
    {
	*t1 = s;
    }

    if (t2)
    {
	*t2 = t;
    }

    return fabsf(distsq);
}



int main (int argc, char **argv)
{
    vec3 axis = {0.f, 0.f, 2.0f};
    vec3 v = {0.1f, 0.3f, 0.5f};

    //vecRotateAboutVec(&v, &axis, PI*0.5f);
    //vecRotateAboutVec2(&v, &axis, PI*0.5f);
    //PRIINTVEC(v);
    //vecRotateAboutVec3(&v, &axis, PI*0.5f);
    //PRIINTVEC(v);
    //vecscale(&v, &v, 1.f/vecsize(&axis));
    //PRIINTVEC(v);

    linedef line1;
    vecset(&line1.point, 1.f, 0.f, 0.f);
    vecset(&line1.dir, -1.f, 0.f, 0.f);
    
    linedef line2;
    vecset(&line2.point, -1.f, 0.f, 0.f);
    vecset(&line2.dir, 1.f, 0.0f, 0.f);

    float t1=99.f, t2=99.f;
    float d = distanceLineLineSq(&line1, &line2, &t1, &t2);

    printf("t1 = %f, t2 = %f, d = %f\n", t1, t2, d);


    return 0;
}


