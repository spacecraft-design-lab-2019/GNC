/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_milqr_api.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 17-Mar-2020 19:06:03
 */

#ifndef _CODER_MILQR_API_H
#define _CODER_MILQR_API_H

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
void milqr(real_T x0[3500], real_T xg[7], real_T u0[1497], real_T u_lims[6],
           real_T x[3500], real_T u[1497], real_T K[8982], boolean_T *result);
void milqr_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray *plhs
               [4]);
void milqr_atexit(void);
void milqr_initialize(void);
void milqr_terminate(void);
void milqr_xil_shutdown(void);
void milqr_xil_terminate(void);

#endif

/*
 * File trailer for _coder_milqr_api.h
 *
 * [EOF]
 */
