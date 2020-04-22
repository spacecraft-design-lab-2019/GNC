/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_get_magnetic_field_series_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 12:32:15
 */

#ifndef _CODER_GET_MAGNETIC_FIELD_SERIES_API_H
#define _CODER_GET_MAGNETIC_FIELD_SERIES_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_get_magnetic_field_series_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
void get_magnetic_field_series(real32_T x0[6], real32_T t[100], real32_T
  B_eci_vec[297], real32_T X[600]);
void get_magnetic_field_series_api(const mxArray * const prhs[2], int32_T nlhs,
  const mxArray *plhs[2]);
void get_magnetic_field_series_atexit(void);
void get_magnetic_field_series_initialize(void);
void get_magnetic_field_series_terminate(void);
void get_magnetic_field_series_xil_terminate(void);

#endif

/*
 * File trailer for _coder_get_magnetic_field_series_api.h
 *
 * [EOF]
 */
