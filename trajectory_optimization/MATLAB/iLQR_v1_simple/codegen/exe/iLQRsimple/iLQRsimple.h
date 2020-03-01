/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: iLQRsimple.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 20:32:25
 */

#ifndef ILQRSIMPLE_H
#define ILQRSIMPLE_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "iLQRsimple_types.h"

/* Function Declarations */
extern void iLQRsimple(const double x0[2], const double xg[2], const double
  utraj0[399], const double Q[4], double R, const double Qf[4], double dt,
  double tol, double xtraj[800], double utraj[399], double K[798]);
extern void iLQRsimple_initialize(void);
extern void iLQRsimple_terminate(void);

#endif

/*
 * File trailer for iLQRsimple.h
 *
 * [EOF]
 */
