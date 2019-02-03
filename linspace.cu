//nvcc linspace.cu -o linspace -lglut -lGL -lm; ./'linspace'
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void linspace(float* vec, float start, float stop, int num, int useEndpoint)
{
/*	Create 'num' evenly spaced samples, calculated over the interval ['start', 'stop'].
 *	The endpoint of the interval can be excluded.
 */

	int q = 0;	
	float div;
	float delta;	
	float step;
		
	if(num > 1)
	{
		if(useEndpoint == 1)
		{
			q = 1;
		}

		div = num - q;
		delta = stop - start;
		step = float(delta)/div;

		for(int i = 0;i<num;i++)
		{
			vec[i] = i*step + start;
		}		
	}
}

int main( void )
{
	int n = 5;
	int R = 1;
	float* arr = (float*)malloc(n*sizeof(float));
	linspace( arr,
			  -1.0, 
			  1.0,
			  n, 1);

	printf("\n\n");
	
    for(int i = 0; i < n ; i++)
    {
        printf("%f\n", arr[i]);
    }
	printf("\n\n");

	free(arr);
}