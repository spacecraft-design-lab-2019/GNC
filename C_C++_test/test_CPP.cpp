// C++ functions to be called from C
// Paul DeTrempe
// 1/10/2019

#include "../eigen-git-mirror/Eigen/Dense"
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// constants for default dipole moments
const double dipole_x = 8.8e-3;
const double dipole_y = 1.372e-2;
const double dipole_z = 8.2e-3;

Vector3d get_B_dot_bang_bang( Vector3d B_dot, Vector3d max_dipoles( double dipole_x, double dipole_y, double dipole_z) );

int main(){
	return 0;
}

// C wrapper function exposing Eigen inputs and outputs
extern "C"{


	double* test_B_dot(double *B_dot, double *max_dipoles)
	{
    // Change array pointers (because arrays are 2nd class citizens in C) to Eigen Vectors (get it? eigenvectors?)
		Map<Vector3d> B_dot_Eigen(B_dot);
		Map<Vector3d> max_dipoles_Eigen(B_dot);
		Vector3d vec_sum;

	// Vector3d commanded_dipoles_Eigen;

	// First try outputting the addition of two vectors
		vec_sum = B_dot_Eigen + max_dipoles_Eigen;

    // Call the C++ function
    // commanded_dipoles_Eigen = detumble_B_dot_bang_bang( B_dot_Eigen, max_dipoles_Eigen);

    // Change Eigen vector to array pointer for output
		double* vec_sum_array = vec_sum.data();

		return vec_sum_array;
	}
}

// C++ control law
Vector3d detumble_B_dot_bang_bang(Vector3d B_dot, Vector3d max_dipoles){
    /*
    !-----------------This function replaces detumble_B_dot (11/17/2019)-----------------------!

    Takes in the rate of change of magnetic field as a 3x1 vector, [nanoTesla/s](though it should be noted that units don't
    actually matter for this function, as it uses the sign of the dot product of the magnetic field rate of change with the
    principal axes of the spacecraft, not the magnitude), and an optional input for the maximum dipole as a 3x1 vector,
    [Amp-m^2].

    It returns a dipole acting to maximally oppose the rate of change of the magnetic field.
    */

    Vector3d M;             // dipole moment [Amp-m^2]

    if( signbit(B_dot(0)) ){        // signbit returns true is a number is negative
    	M(0) = max_dipoles(0);
    }
    else{
    	M(0) = -max_dipoles(0);
    }

    if( signbit(B_dot(1)) ){        // signbit returns true is a number is negative
    	M(1) = max_dipoles(1);
    }
    else{
    	M(1) = -max_dipoles(1);
    }

    if( signbit(B_dot(2)) ){        // signbit returns true is a number is negative
    	M(2) = max_dipoles(2);
    }
    else{
    	M(2) = -max_dipoles(2);
    }


    return M;

}

