/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_get_magnetic_field_series_api.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 12:32:15
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_get_magnetic_field_series_api.h"
#include "_coder_get_magnetic_field_series_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131467U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "get_magnetic_field_series",         /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[6]);
static const mxArray *b_emlrt_marshallOut(const real32_T u[600]);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *t, const
  char_T *identifier, real32_T y[100]);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[100]);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[6]);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *x0, const
  char_T *identifier, real32_T y[6]);
static const mxArray *emlrt_marshallOut(const real32_T u[297]);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[100]);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[6]
 * Return Type  : void
 */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[6])
{
  e_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const real32_T u[600]
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(const real32_T u[600])
{
  const mxArray *y;
  const mxArray *m1;
  static const int32_T iv1[2] = { 6, 100 };

  real32_T *pData;
  int32_T i1;
  int32_T i;
  int32_T b_i;
  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv1, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m1);
  i1 = 0;
  for (i = 0; i < 100; i++) {
    for (b_i = 0; b_i < 6; b_i++) {
      pData[i1] = u[b_i + 6 * i];
      i1++;
    }
  }

  emlrtAssign(&y, m1);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *t
 *                const char_T *identifier
 *                real32_T y[100]
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *t, const
  char_T *identifier, real32_T y[100])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(t), &thisId, y);
  emlrtDestroyArray(&t);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[100]
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[100])
{
  f_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[6]
 * Return Type  : void
 */
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[6])
{
  static const int32_T dims[1] = { 6 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 1U, dims);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *x0
 *                const char_T *identifier
 *                real32_T y[6]
 * Return Type  : void
 */
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *x0, const
  char_T *identifier, real32_T y[6])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(x0), &thisId, y);
  emlrtDestroyArray(&x0);
}

/*
 * Arguments    : const real32_T u[297]
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real32_T u[297])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 3, 99 };

  real32_T *pData;
  int32_T i0;
  int32_T i;
  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m0);
  i0 = 0;
  for (i = 0; i < 99; i++) {
    pData[i0] = u[3 * i];
    i0++;
    pData[i0] = u[1 + 3 * i];
    i0++;
    pData[i0] = u[2 + 3 * i];
    i0++;
  }

  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[100]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[100])
{
  static const int32_T dims[2] = { 1, 100 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 2U, dims);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray * const prhs[2]
 *                int32_T nlhs
 *                const mxArray *plhs[2]
 * Return Type  : void
 */
void get_magnetic_field_series_api(const mxArray * const prhs[2], int32_T nlhs,
  const mxArray *plhs[2])
{
  real32_T x0[6];
  real32_T t[100];
  real32_T B_eci_vec[297];
  real32_T X[600];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "x0", x0);
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t", t);

  /* Invoke the target function */
  get_magnetic_field_series(x0, t, B_eci_vec, X);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(B_eci_vec);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(X);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void get_magnetic_field_series_atexit(void)
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
  get_magnetic_field_series_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void get_magnetic_field_series_initialize(void)
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
void get_magnetic_field_series_terminate(void)
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
 * File trailer for _coder_get_magnetic_field_series_api.c
 *
 * [EOF]
 */
