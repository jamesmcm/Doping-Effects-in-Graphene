#include <stdio.h>
#include <stdlib.h>

#define P 500
#define bins 35

int main()
{
	double input[P] = {0.0};		/* this is the input array */
	double minval = 30;			/* minimum value of the minimum bin */
	double maxval = 70;			/* maximum value of the maximum bin */
	int output[bins] = {0};
	double binsize = 0.0;
	int i, j, k, m;				/* not sure if necessary, but my compiler seemed to need this */

	binsize = (maxval - minval) / bins;

	for (i=0; i<P; ++i)			/* for testing purposes generates a random array */
		input[i] = rand() % 100;

	for (j=0; j<P; ++j)
	{
		for (k=0; k<bins; ++k)
		{
			if (input[j] >= (minval + k*binsize) && input[j] < (minval + (k+1)*binsize))
				++output[k];
		}
	}
	
	for (m=0; m<bins; ++m)			/* just prints out the number of elements in each bin, histogram generating function goes here */
	{
		printf("%d\n",output[m]);
	}

	return 0;
}

