/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: iLQRv1.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 07-Feb-2020 14:37:48
 */

#ifndef ILQRV1_H
#define ILQRV1_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "iLQRv1_types.h"

/* Function Declarations */
extern void iLQRv1(const real_T x0[2], const real_T xg[2], const real_T utraj0
                   [399], const real_T Q[4], real_T R, const real_T Qf[4],
                   real_T dt, real_T tol, real_T xtraj[800], real_T utraj[399],
                   real_T K[798]);
extern void iLQRv1_initialize(void);
extern void iLQRv1_terminate(void);

#endif

/*
 * File trailer for iLQRv1.h
 *
 * [EOF]
 */
