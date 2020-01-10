// Paul DeTrempe
// Script to test calling C++ from C

#include <stdio.h>
#include "test_CPP.h"

// constants for default dipole moments
const double dipole_x = 8.8e-3;
const double dipole_y = 1.372e-2;
const double dipole_z = 8.2e-3;

int main()
{

    printf("%s", "Hello World");

    // Initialize array and call C++ function
    double B_dot[3] = { 1.0, 1.0, 1.0};
    double max_dipoles[3] = {dipole_x, dipole_y, dipole_z};
    double *dipole_command;

    // B_dot = { 1.0, 1.0, 1.0};
    // max_dipoles = {dipole_x, dipole_y, dipole_z};

    dipole_command = test_B_dot(B_dot, max_dipoles);


    for(int i = 0; i < 3; i++)
      printf("%f ", dipole_command[i]);

    return 0;
}