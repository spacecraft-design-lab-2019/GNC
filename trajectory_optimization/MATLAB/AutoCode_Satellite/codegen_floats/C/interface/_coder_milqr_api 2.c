/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_milqr_api.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 19-Mar-2020 19:23:47
 */

/* Include Files */
#include "_coder_milqr_api.h"
#include "_coder_milqr_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131483U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "milqr",                             /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[3500]);
static const mxArray *b_emlrt_marshallOut(const real32_T u[1497]);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *xg, const
  char_T *identifier, real32_T y[7]);
static const mxArray *c_emlrt_marshallOut(const real32_T u[8982]);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[7]);
static const mxArray *d_emlrt_marshallOut(const boolean_T u);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u0, const
  char_T *identifier, real32_T y[1497]);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *x0, const
  char_T *identifier, real32_T y[3500]);
static const mxArray *emlrt_marshallOut(const real32_T u[3500]);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[1497]);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u_lims,
  const char_T *identifier, real32_T y[6]);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[6]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[3500]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[7]);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[1497]);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[6]);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[3500]
 * Return Type  : void
 */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[3500])
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const real32_T u[1497]
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(const real32_T u[1497])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 3, 499 };

  real32_T *pData;
  int32_T i;
  int32_T b_i;
  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m);
  i = 0;
  for (b_i = 0; b_i < 499; b_i++) {
    pData[i] = u[3 * b_i];
    i++;
    pData[i] = u[3 * b_i + 1];
    i++;
    pData[i] = u[3 * b_i + 2];
    i++;
  }

  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *xg
 *                const char_T *identifier
 *                real32_T y[7]
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *xg, const
  char_T *identifier, real32_T y[7])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(xg), &thisId, y);
  emlrtDestroyArray(&xg);
}

/*
 * Arguments    : const real32_T u[8982]
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrt_marshallOut(const real32_T u[8982])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[3] = { 3, 6, 499 };

  real32_T *pData;
  int32_T i;
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  y = NULL;
  m = emlrtCreateNumericArray(3, iv, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m);
  i = 0;
  for (b_i = 0; b_i < 499; b_i++) {
    for (c_i = 0; c_i < 6; c_i++) {
      i1 = 3 * c_i + 18 * b_i;
      pData[i] = u[i1];
      i++;
      pData[i] = u[i1 + 1];
      i++;
      pData[i] = u[i1 + 2];
      i++;
    }
  }

  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[7]
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[7])
{
  j_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const boolean_T u
 * Return Type  : const mxArray *
 */
static const mxArray *d_emlrt_marshallOut(const boolean_T u)
{
  const mxArray *y;
  const mxArray *m;
  y = NULL;
  m = emlrtCreateLogicalScalar(u);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u0
 *                const char_T *identifier
 *                real32_T y[1497]
 * Return Type  : void
 */
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u0, const
  char_T *identifier, real32_T y[1497])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(u0), &thisId, y);
  emlrtDestroyArray(&u0);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *x0
 *                const char_T *identifier
 *                real32_T y[3500]
 * Return Type  : void
 */
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *x0, const
  char_T *identifier, real32_T y[3500])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(x0), &thisId, y);
  emlrtDestroyArray(&x0);
}

/*
 * Arguments    : const real32_T u[3500]
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real32_T u[3500])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 7, 500 };

  real32_T *pData;
  int32_T i;
  int32_T b_i;
  int32_T c_i;
  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m);
  i = 0;
  for (b_i = 0; b_i < 500; b_i++) {
    for (c_i = 0; c_i < 7; c_i++) {
      pData[i] = u[c_i + 7 * b_i];
      i++;
    }
  }

  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[1497]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[1497])
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u_lims
 *                const char_T *identifier
 *                real32_T y[6]
 * Return Type  : void
 */
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u_lims,
  const char_T *identifier, real32_T y[6])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  h_emlrt_marshallIn(sp, emlrtAlias(u_lims), &thisId, y);
  emlrtDestroyArray(&u_lims);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[6]
 * Return Type  : void
 */
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[6])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[3500]
 * Return Type  : void
 */
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[3500])
{
  static const int32_T dims[2] = { 7, 500 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 2U, dims);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[7]
 * Return Type  : void
 */
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[7])
{
  static const int32_T dims[1] = { 7 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 1U, dims);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[1497]
 * Return Type  : void
 */
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[1497])
{
  static const int32_T dims[2] = { 3, 499 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 2U, dims);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[6]
 * Return Type  : void
 */
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[6])
{
  static const int32_T dims[2] = { 3, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 2U, dims);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray * const prhs[4]
 *                int32_T nlhs
 *                const mxArray *plhs[4]
 * Return Type  : void
 */
void milqr_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray *plhs
               [4])
{
  static real32_T x0[3500];
  real32_T xg[7];
  real32_T u0[1497];
  real32_T u_lims[6];
  static real32_T x[3500];
  real32_T u[1497];
  static real32_T K[8982];
  boolean_T result;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "x0", x0);
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "xg", xg);
  e_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "u0", u0);
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "u_lims", u_lims);

  /* Invoke the target function */
  milqr(x0, xg, u0, u_lims, x, u, K, &result);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(x);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(u);
  }

  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(K);
  }

  if (nlhs > 3) {
    plhs[3] = d_emlrt_marshallOut(result);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void milqr_atexit(void)
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
  milqr_xil_terminate();
  milqr_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void milqr_initialize(void)
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
void milqr_terminate(void)
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
 * File trailer for _coder_milqr_api.c
 *
 * [EOF]
 */
