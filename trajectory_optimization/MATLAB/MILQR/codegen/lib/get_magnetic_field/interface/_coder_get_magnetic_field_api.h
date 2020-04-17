/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_get_magnetic_field_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 11-Apr-2020 23:47:25
 */

#ifndef _CODER_GET_MAGNETIC_FIELD_API_H
#define _CODER_GET_MAGNETIC_FIELD_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_get_magnetic_field_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
void get_magnetic_field(real32_T lat, real32_T lon, real32_T alt, real32_T year,
  real32_T order, real32_T B_NED[3]);
void get_magnetic_field_api(const mxArray * const prhs[5], int32_T nlhs, const
  mxArray *plhs[1]);
void get_magnetic_field_atexit(void);
void get_magnetic_field_initialize(void);
void get_magnetic_field_terminate(void);
void get_magnetic_field_xil_terminate(void);

#endif

/*
 * File trailer for _coder_get_magnetic_field_api.h
 *
 * [EOF]
 */
