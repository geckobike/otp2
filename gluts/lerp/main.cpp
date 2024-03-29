#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ENABLE_TIMER
#include "engine/engine.h"
#include <assert.h>



int main (int argc, char **argv)
{
	float dts[2] = {0.1f, 0.01f};

	FILE * file[2];

	file[0] = fopen("file0.txt", "w");
	file[1] = fopen("file1.txt", "w");

	for (int i=0; i<2; i++)
	{
		float t = 0.f;
		float dt = dts[i];
		float target = -2.f;
		float x = 0.f;
		float spring = 3.f;

		while (t<3.f)
		{
			t = t + dt;
			x += (target-x)*(1.f - exp(-dt*spring));
			fprintf(file[i], "%f, %f\n", t, x);
		}
	}

	fclose(file[0]);
	fclose(file[1]);

	return 0;
}


