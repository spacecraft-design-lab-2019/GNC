/* C_linkage_detumble.cpp 
Example archetype script for wrapping C++ functions in C to make them accessible to C (and therefore CircuitPython).
Paul DeTrempe
detrempe@stanford.edu

*/


#include "C_linkage_detumble.h"
#include "detumble_algorithms.h"
#include <iostream>
#include "../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;

extern "C" {

	void detumble_B_dot(double* omega, double* B, double* max_dipoles, double* commanded_dipole){
		// Takes in two pointers to C/C++ arrays and one pointer to where you would like the answer stored
		// Updates the contents of the answer array

  // ----------------------Try to add two vectors using Eigen (map to Eigen, add, map back)-------

		// Create Eigen vectors from C/C++ arrays
		Vector3d w_vec = Map<Vector3d>(omega);
		Vector3d B_vec = Map<Vector3d>(B);
		Vector3d max_dipole_vec = Map<Vector3d>(max_dipoles);
		Vector3d commanded_dipole_vec = Map<Vector3d>(commanded_dipole);

		// Perform operation on vector/matrix (using Eigen notation)
		commanded_dipole_vec = detumble_B_cross_bang_bang(w_vec, B_vec, 1.0, max_dipole_vec);


		// Copy contents back into array preallocated for holding answer used at higher level by C calling function
		for (int i = 0; i < 3; ++i)
		{
			commanded_dipole[i] = commanded_dipole_vec(i);
		}

	}



} /* extern "C" */