/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_cholesky_api.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 20:17:58
 */

#ifndef _CODER_CHOLESKY_API_H
#define _CODER_CHOLESKY_API_H

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
extern void cholesky(real_T A_data[], int32_T A_size[2], real_T R_data[],
                     int32_T R_size[2]);
extern void cholesky_api(const mxArray * const prhs[1], int32_T nlhs, const
  mxArray *plhs[1]);
extern void cholesky_atexit(void);
extern void cholesky_initialize(void);
extern void cholesky_terminate(void);
extern void cholesky_xil_shutdown(void);
extern void cholesky_xil_terminate(void);

#endif

/*
 * File trailer for _coder_cholesky_api.h
 *
 * [EOF]
 */
