/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_ilqrCar_api.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 01-Mar-2020 20:37:52
 */

#ifndef _CODER_ILQRCAR_API_H
#define _CODER_ILQRCAR_API_H

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
void ilqrCar(real_T x0[2004], real_T xg[4], real_T u0[1000], real_T u_lims[4],
             real_T x[2004], real_T u[1000], real_T K[4000], boolean_T *result);
void ilqrCar_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray
                 *plhs[4]);
void ilqrCar_atexit(void);
void ilqrCar_initialize(void);
void ilqrCar_terminate(void);
void ilqrCar_xil_shutdown(void);
void ilqrCar_xil_terminate(void);

#endif

/*
 * File trailer for _coder_ilqrCar_api.h
 *
 * [EOF]
 */
