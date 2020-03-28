/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_ilqrCar_mex.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 01-Mar-2020 20:37:52
 */

/* Include Files */
#include "_coder_ilqrCar_mex.h"
#include "_coder_ilqrCar_api.h"

/* Function Declarations */
MEXFUNCTION_LINKAGE void ilqrCar_mexFunction(int32_T nlhs, mxArray *plhs[4],
  int32_T nrhs, const mxArray *prhs[4]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[4]
 *                int32_T nrhs
 *                const mxArray *prhs[4]
 * Return Type  : void
 */
void ilqrCar_mexFunction(int32_T nlhs, mxArray *plhs[4], int32_T nrhs, const
  mxArray *prhs[4])
{
  const mxArray *outputs[4];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4, 7,
                        "ilqrCar");
  }

  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 7,
                        "ilqrCar");
  }

  /* Call the function. */
  ilqrCar_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(ilqrCar_atexit);

  /* Module initialization. */
  ilqrCar_initialize();

  /* Dispatch the entry-point. */
  ilqrCar_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  ilqrCar_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_ilqrCar_mex.c
 *
 * [EOF]
 */
