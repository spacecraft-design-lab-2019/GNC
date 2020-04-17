/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_milqr_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 16-Apr-2020 13:51:41
 */

#ifndef _CODER_MILQR_API_H
#define _CODER_MILQR_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_milqr_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
void milqr(real32_T x0[7000], real32_T xg[7], real32_T u0[2997], real32_T
           u_lims[6], real32_T B_ECI[3000], real32_T x[7000], real32_T u[2997],
           real32_T K[17982], boolean_T *result);
void milqr_api(const mxArray * const prhs[5], int32_T nlhs, const mxArray *plhs
               [4]);
void milqr_atexit(void);
void milqr_initialize(void);
void milqr_terminate(void);
void milqr_xil_terminate(void);

#endif

/*
 * File trailer for _coder_milqr_api.h
 *
 * [EOF]
 */
