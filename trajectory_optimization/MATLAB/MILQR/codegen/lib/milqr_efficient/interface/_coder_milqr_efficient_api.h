/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_milqr_efficient_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 20-Apr-2020 22:26:32
 */

#ifndef _CODER_MILQR_EFFICIENT_API_H
#define _CODER_MILQR_EFFICIENT_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_milqr_efficient_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
void milqr_efficient(real32_T x0[700], real32_T xg[7], real32_T u0[297],
                     real32_T u_lims[6], real32_T dt, real32_T B_ECI[300],
                     real32_T x[700], real32_T u[297], real32_T K[1782],
                     boolean_T *result);
void milqr_efficient_api(const mxArray * const prhs[6], int32_T nlhs, const
  mxArray *plhs[4]);
void milqr_efficient_atexit(void);
void milqr_efficient_initialize(void);
void milqr_efficient_terminate(void);
void milqr_efficient_xil_terminate(void);

#endif

/*
 * File trailer for _coder_milqr_efficient_api.h
 *
 * [EOF]
 */
