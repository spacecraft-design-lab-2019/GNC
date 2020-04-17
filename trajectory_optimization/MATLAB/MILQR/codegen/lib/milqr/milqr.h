/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 16-Apr-2020 13:51:41
 */

#ifndef MILQR_H
#define MILQR_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "milqr_types.h"

/* Function Declarations */
void milqr(const float x0[7000], const float xg[7], const float u0[2997], const
           float u_lims[6], const float B_ECI[3000], float x[7000], float u[2997],
           float K[17982], bool *result);
void milqr_initialize(void);
void milqr_terminate(void);

#endif

/*
 * File trailer for milqr.h
 *
 * [EOF]
 */
