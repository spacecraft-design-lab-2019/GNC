/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ilqrCar.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 01-Mar-2020 20:37:52
 */

#ifndef ILQRCAR_H
#define ILQRCAR_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "ilqrCar_types.h"

#include <stdbool.h>

/* Function Declarations */
void ilqrCar(const double x0[2004], const double xg[4], const double u0[1000],
             const double u_lims[4], double x[2004], double u[1000], double K
             [4000], bool *result);

#endif

/*
 * File trailer for ilqrCar.h
 *
 * [EOF]
 */
