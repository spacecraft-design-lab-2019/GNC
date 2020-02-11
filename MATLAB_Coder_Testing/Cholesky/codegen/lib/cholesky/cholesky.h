/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: cholesky.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 20:17:58
 */

#ifndef CHOLESKY_H
#define CHOLESKY_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "cholesky_types.h"

/* Function Declarations */
extern void cholesky(const double A_data[], const int A_size[2], double R_data[],
                     int R_size[2]);
extern void cholesky_initialize(void);
extern void cholesky_terminate(void);

#endif

/*
 * File trailer for cholesky.h
 *
 * [EOF]
 */
