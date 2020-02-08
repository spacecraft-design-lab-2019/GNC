/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: iLQRv1.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 06-Feb-2020 15:05:27
 */

#ifndef ILQRV1_H
#define ILQRV1_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "iLQRv1_types.h"

/* Function Declarations */
extern void iLQRv1(const double x0[2], const double xg[2], const double
                   utraj0_data[], const int utraj0_size[2], const double Q[4],
                   double R, const double Qf[4], double dt, double tol, double
                   xtraj_data[], int xtraj_size[2], double utraj_data[], int
                   utraj_size[2], double K_data[], int K_size[3]);
extern void iLQRv1_initialize(void);
extern void iLQRv1_terminate(void);

#endif

/*
 * File trailer for iLQRv1.h
 *
 * [EOF]
 */
