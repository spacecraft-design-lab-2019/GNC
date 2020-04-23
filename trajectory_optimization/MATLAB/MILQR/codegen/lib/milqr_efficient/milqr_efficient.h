/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr_efficient.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 16:18:57
 */

#ifndef MILQR_EFFICIENT_H
#define MILQR_EFFICIENT_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "milqr_efficient_types.h"

/* Function Declarations */
void milqr_efficient(const double x0[7], const double xg[7], const double u0[297],
                     const double u_lims[6], double dt, const double B_ECI[297],
                     double x[700], double u[297], double K[1782], bool *result);
void milqr_efficient_initialize(void);
void milqr_efficient_terminate(void);

#endif

/*
 * File trailer for milqr_efficient.h
 *
 * [EOF]
 */
