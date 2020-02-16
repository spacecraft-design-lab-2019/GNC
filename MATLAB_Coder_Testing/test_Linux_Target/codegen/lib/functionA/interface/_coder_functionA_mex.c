/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_functionA_mex.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 15:18:21
 */

/* Include Files */
#include "_coder_functionA_mex.h"
#include "_coder_functionA_api.h"

/* Function Declarations */
MEXFUNCTION_LINKAGE void functionA_mexFunction(int32_T nlhs, mxArray *plhs[5],
  int32_T nrhs, const mxArray *prhs[3]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[5]
 *                int32_T nrhs
 *                const mxArray *prhs[3]
 * Return Type  : void
 */
void functionA_mexFunction(int32_T nlhs, mxArray *plhs[5], int32_T nrhs, const
  mxArray *prhs[3])
{
  const mxArray *outputs[5];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4, 9,
                        "functionA");
  }

  if (nlhs > 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 9,
                        "functionA");
  }

  /* Call the function. */
  functionA_api(prhs, nlhs, outputs);

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
  mexAtExit(functionA_atexit);

  /* Module initialization. */
  functionA_initialize();

  /* Dispatch the entry-point. */
  functionA_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  functionA_terminate();
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
 * File trailer for _coder_functionA_mex.c
 *
 * [EOF]
 */
