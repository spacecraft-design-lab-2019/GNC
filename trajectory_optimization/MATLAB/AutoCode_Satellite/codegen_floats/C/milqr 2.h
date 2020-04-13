/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 19-Mar-2020 19:23:47
 */

#ifndef MILQR_H
#define MILQR_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>

/* Function Declarations */
void milqr(const float x0[3500], const float xg[7], const float u0[1497], const
           float u_lims[6], float x[3500], float u[1497], float K[8982],
           bool *result);

#endif

/*
 * File trailer for milqr.h
 *
 * [EOF]
 */
