/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_functionA_api.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 15:18:21
 */

#ifndef _CODER_FUNCTIONA_API_H
#define _CODER_FUNCTIONA_API_H

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
extern void functionA(real_T A[16], real_T B[240], real_T C[4], real_T Z[16],
                      real_T Y[4], real_T X[12], real_T W[16], real_T V[16]);
extern void functionA_api(const mxArray * const prhs[3], int32_T nlhs, const
  mxArray *plhs[5]);
extern void functionA_atexit(void);
extern void functionA_initialize(void);
extern void functionA_terminate(void);
extern void functionA_xil_shutdown(void);
extern void functionA_xil_terminate(void);

#endif

/*
 * File trailer for _coder_functionA_api.h
 *
 * [EOF]
 */
