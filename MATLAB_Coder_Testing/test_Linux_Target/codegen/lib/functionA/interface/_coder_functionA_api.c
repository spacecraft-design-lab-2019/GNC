/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_functionA_api.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 15:18:21
 */

/* Include Files */
#include "_coder_functionA_api.h"
#include "_coder_functionA_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131483U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "functionA",                         /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[16];
static const mxArray *b_emlrt_marshallOut(const real_T u[4]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *B, const
  char_T *identifier))[240];
static const mxArray *c_emlrt_marshallOut(const real_T u[12]);
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[240];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *C, const
  char_T *identifier))[4];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *A, const
  char_T *identifier))[16];
static const mxArray *emlrt_marshallOut(const real_T u[16]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4];
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[16];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[240];
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4];

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[16]
 */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[16]
{
  real_T (*y)[16];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const real_T u[4]
 * Return Type  : const mxArray *
 */
  static const mxArray *b_emlrt_marshallOut(const real_T u[4])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[1] = { 0 };

  static const int32_T iv1[1] = { 4 };

  y = NULL;
  m = emlrtCreateNumericArray(1, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 1);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *B
 *                const char_T *identifier
 * Return Type  : real_T (*)[240]
 */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *B, const
  char_T *identifier))[240]
{
  real_T (*y)[240];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(B), &thisId);
  emlrtDestroyArray(&B);
  return y;
}
/*
 * Arguments    : const real_T u[12]
 * Return Type  : const mxArray *
 */
  static const mxArray *c_emlrt_marshallOut(const real_T u[12])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 4, 3 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[240]
 */
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[240]
{
  real_T (*y)[240];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *C
 *                const char_T *identifier
 * Return Type  : real_T (*)[4]
 */
  static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *C,
  const char_T *identifier))[4]
{
  real_T (*y)[4];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(C), &thisId);
  emlrtDestroyArray(&C);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *A
 *                const char_T *identifier
 * Return Type  : real_T (*)[16]
 */
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *A, const
  char_T *identifier))[16]
{
  real_T (*y)[16];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(A), &thisId);
  emlrtDestroyArray(&A);
  return y;
}
/*
 * Arguments    : const real_T u[16]
 * Return Type  : const mxArray *
 */
  static const mxArray *emlrt_marshallOut(const real_T u[16])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 4, 4 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[4]
 */
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4]
{
  real_T (*y)[4];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[16]
 */
  static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[16]
{
  real_T (*ret)[16];
  static const int32_T dims[2] = { 4, 4 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[16])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[240]
 */
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[240]
{
  real_T (*ret)[240];
  static const int32_T dims[3] = { 4, 3, 20 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 3U, dims);
  ret = (real_T (*)[240])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[4]
 */
  static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4]
{
  real_T (*ret)[4];
  static const int32_T dims[1] = { 4 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray * const prhs[3]
 *                int32_T nlhs
 *                const mxArray *plhs[5]
 * Return Type  : void
 */
void functionA_api(const mxArray * const prhs[3], int32_T nlhs, const mxArray
                   *plhs[5])
{
  real_T (*Z)[16];
  real_T (*Y)[4];
  real_T (*X)[12];
  real_T (*W)[16];
  real_T (*V)[16];
  real_T (*A)[16];
  real_T (*B)[240];
  real_T (*C)[4];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  Z = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  Y = (real_T (*)[4])mxMalloc(sizeof(real_T [4]));
  X = (real_T (*)[12])mxMalloc(sizeof(real_T [12]));
  W = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  V = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));

  /* Marshall function inputs */
  A = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "A");
  B = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "B");
  C = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "C");

  /* Invoke the target function */
  functionA(*A, *B, *C, *Z, *Y, *X, *W, *V);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*Z);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(*Y);
  }

  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(*X);
  }

  if (nlhs > 3) {
    plhs[3] = emlrt_marshallOut(*W);
  }

  if (nlhs > 4) {
    plhs[4] = emlrt_marshallOut(*V);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void functionA_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  functionA_xil_terminate();
  functionA_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void functionA_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void functionA_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_functionA_api.c
 *
 * [EOF]
 */
