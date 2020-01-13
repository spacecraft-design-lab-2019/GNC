/* main.c */

#include "foo.h"
#include <stdio.h>


int main() {
	int size;
	size = 3;
	double vec1[size], vec2[size], total_C[size];	



	for (int i = 0; i < size; ++i)
	{
		vec1[i] = 1.1;
		vec2[i] = 2.2;
	}

	vec_sum(vec1, vec2, total_C);

	for (int i = 0; i < size; ++i)
	{
		printf("%f\n", total_C[i] );
	}
	

	return 0;
}