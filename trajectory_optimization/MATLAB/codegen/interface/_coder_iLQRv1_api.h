/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_iLQRv1_api.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 07-Feb-2020 14:37:48
 */

#ifndef _CODER_ILQRV1_API_H
#define _CODER_ILQRV1_API_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void iLQRv1(real_T x0[2], real_T xg[2], real_T utraj0[399], real_T Q[4],
                   real_T R, real_T Qf[4], real_T dt, real_T tol, real_T xtraj
                   [800], real_T utraj[399], real_T K[798]);
extern void iLQRv1_api(const mxArray * const prhs[8], int32_T nlhs, const
  mxArray *plhs[3]);
extern void iLQRv1_atexit(void);
extern void iLQRv1_initialize(void);
extern void iLQRv1_terminate(void);
extern void iLQRv1_xil_shutdown(void);
extern void iLQRv1_xil_terminate(void);

#endif

/*
 * File trailer for _coder_iLQRv1_api.h
 *
 * [EOF]
 */
