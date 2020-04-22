/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr_efficient.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 20-Apr-2020 22:26:32
 */

#ifndef MILQR_EFFICIENT_H
#define MILQR_EFFICIENT_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
//#include "rtwtypes.h"
#include "milqr_efficient_types.h"

/* Function Declarations */
void milqr_efficient(const float x0[700], const float xg[7], const float u0[297],
                     const float u_lims[6], float dt, const float B_ECI[300],
                     float x[700], float u[297], float K[1782], boolean_T
                     *result);
void milqr_efficient_initialize(void);
void milqr_efficient_terminate(void);

#endif

/*
 * File trailer for milqr_efficient.h
 *
 * [EOF]
 */
