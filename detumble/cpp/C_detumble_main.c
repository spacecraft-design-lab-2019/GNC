/* main.c */

#include "detumble_algorithms.h"
#include <stdio.h>


int main() {
	int size;
	size = 3;
	double B_dot[size], max_dipoles[size], commanded_dipole[size];	



	for (int i = 0; i < size; ++i)
	{
		B_dot[i] = -1.0;
		max_dipoles[i] = 3.0;
	}

    detumble_B_dot_C(B_dot, max_dipoles, commanded_dipole);

	for (int i = 0; i < size; ++i)
	{
		printf("%f\n", commanded_dipole[i] );
	}
	

	return 0;
}