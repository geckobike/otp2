#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine/engine.h"

typedef unsigned long long u64;
typedef unsigned short u16;

static double getTime()
{
    struct timeval now;
    long long sec;
    long long microsec;
    gettimeofday(&now,NULL);
    sec = now.tv_sec;
    microsec = now.tv_usec;
    return (double)sec + ( (double)microsec * 1e-6 );
}

u64 newrndseedval = 0x0000000021454a53;

u32 rnd(void)
{
    return (u32)(((newrndseedval = ((u64)214013)*newrndseedval + ((u64)2531011))>>16)&0x7fff);
}

float rndf(void)
{
    return ((float)rnd()*(1.0f/32767.0f));
}


float rndrangef(float min, float max)
{
    return (min + (rndf()*(max - min)));
}

#define randomrange rndrangef

// Binary search for a value that satisfies: array[idx] <= value < array[idx+1]
// Array must be at least 2 values long

template <class T>
int search(T* array, T value, int N)
{
    int	lo, hi, mid;

    if ((value < array[0]) || (value >= array[N-1])) return -1;
    
    lo=0;
    hi=N-1;
    
    while (lo!=hi)
    {
	mid=(hi+lo)/2;
	
	if (value < array[mid])
	{
	    hi=mid;
	}
	else
	{
	    if (value < array[mid+1]) return mid; // Found the right interval
	    lo = mid+1;
	}
    }
    return lo;
}

int main (int argc, char **argv)
{
    const u32 N = 250;
    u32 indices[N];
    float maxDiff = 10250.f;
    //float maxDiff = 30.f;
    
    indices[0] = 0;
    
    for (u32 i=0;i<(N-1);i++)
    {
	u16 diff = (u16)randomrange(30.f, maxDiff);
	indices[i+1] = indices[i] + diff;
	maxDiff = (maxDiff - 30.f) * 0.98f + 30.f;
	//maxDiff += 130.f;
    }
    printf("Range is 0 to %d\n--------------------------\n", indices[N-1]);

    u32 lookfor = -1;
   

    double time = getTime();
    int idx;
    for (u32 i=0; i<1; i++)
    {
	idx = search(indices, lookfor, N);
    }
    printf("time taken = %f\n", (getTime() - time)*1e3);


    indices[0] = 0;
    indices[1] = 100;
    lookfor = -1;
    idx = search(indices, lookfor, 2);
    

    if (idx >=0)
    {
	printf("%d is located within %d and %d (idx == %d)\n", lookfor, indices[idx], indices[idx+1], idx);
    }
    else
    {
	printf("Not Found\n");
    }
    
    return 0;
}


