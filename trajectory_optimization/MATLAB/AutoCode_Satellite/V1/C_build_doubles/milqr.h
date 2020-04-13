/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 19-Mar-2020 19:51:09
 */

#ifndef MILQR_H
#define MILQR_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>

typedef bool boolean_T;

/* Function Declarations */
void milqr(const double x0[3500], const double xg[7], const double u0[1497],
           const double u_lims[6], double x[3500], double u[1497], double K[8982],
           boolean_T *result);

#endif

/*
 * File trailer for milqr.h
 *
 * [EOF]
 */
