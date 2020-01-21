/* foo.cpp */

#include "foo.h"
#include <iostream>
#include "../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;

extern "C" {

	void vec_sum(double* arg1, double* arg2, double* total){
		// Takes in two pointers to C/C++ arrays and one pointer to where you would like the answer stored
		// Updates the contents of the answer array


  // ----------------------Try to add two 3x1 vectors using arrays in C----------------------
  // int size;

  // double sum[size];
  // for (int i = 0; i < 3; ++i)
  // {
  // 	total[i] = arg1[i] + arg2[i];
  // }

  // ----------------------Try to add two vectors using Eigen (map to Eigen, add, map back)-------

		// Create Eigen vectors from C/C++ arrays
		Vector3d vector1 = Map<Vector3d>(arg1);
		Vector3d vector2 = Map<Vector3d>(arg2);
		Vector3d added_vectors = Map<Vector3d>(total);

		// Perform operation on vector/matrix (using Eigen notation)
		added_vectors = vector1 + vector2;

		// Copy contents back into array preallocated for holding answer used at higher level by C calling function
		for (int i = 0; i < 3; ++i)
		{
			total[i] = added_vectors(i);
		}

	}



} /* extern "C" */