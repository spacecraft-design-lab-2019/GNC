/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: functionA.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 15:18:21
 */

/* Include Files */
#include "functionA.h"
#include <math.h>

/* Function Definitions */

/*
 * Testing MATLAB  coder for different matlab functionalities
 *  Will try to run the tests using a .mex file
 *
 *  * Creating a matrix inside the function
 *  * Backlash (linear solvers)
 *  * Matrix Slicing
 *  * 3D matrices
 *  * Matrix multiplication
 *  * Matrix transposes
 *  * Returning multiple values from a function
 * Arguments    : const double A[16]
 *                const double B[240]
 *                const double C[4]
 *                double Z[16]
 *                double Y[4]
 *                double X[12]
 *                double W[16]
 *                double V[16]
 * Return Type  : void
 */
void functionA(const double A[16], const double B[240], const double C[4],
               double Z[16], double Y[4], double X[12], double W[16], double V
               [16])
{
  int i;
  static const signed char iv[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 1 };

  double b_A[16];
  signed char ipiv[4];
  int j;
  int mmj_tmp;
  int b;
  int jj;
  int jp1j;
  int iy;
  int jA;
  int ix;
  double smax;
  int k;
  double s;
  double d;
  static const signed char b_b[16] = { 1, 0, 0, 0, 3, 1, 0, 0, 0, 0, 1, 0, 0, 5,
    0, 1 };

  /*  A = (4 x 4) matrix */
  /*  B = (4 x 3 x 20) matrix */
  /*  C = (4 x 1) vector */
  /*  Test 0: Creating a matrix in the function */
  for (i = 0; i < 16; i++) {
    Z[i] = iv[i];
    b_A[i] = A[i];
  }

  Z[4] = 3.0;
  Z[13] = 5.0;

  /*  Test 1: Solving linear systems */
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (j = 0; j < 3; j++) {
    mmj_tmp = 2 - j;
    b = j * 5;
    jj = j * 5;
    jp1j = b + 2;
    iy = 4 - j;
    jA = 0;
    ix = b;
    smax = fabs(b_A[jj]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = fabs(b_A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (b_A[jj + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = (signed char)(iy + 1);
        smax = b_A[j];
        b_A[j] = b_A[iy];
        b_A[iy] = smax;
        ix = j + 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
      }

      i = (jj - j) + 4;
      for (iy = jp1j; iy <= i; iy++) {
        b_A[iy - 1] /= b_A[jj];
      }
    }

    iy = b + 4;
    jA = jj;
    for (k = 0; k <= mmj_tmp; k++) {
      smax = b_A[iy];
      if (b_A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 6;
        jp1j = (jA - j) + 8;
        for (b = i; b <= jp1j; b++) {
          b_A[b - 1] += b_A[ix] * -smax;
          ix++;
        }
      }

      iy += 4;
      jA += 4;
    }
  }

  Y[0] = C[0];
  Y[1] = C[1];
  Y[2] = C[2];
  Y[3] = C[3];
  if (ipiv[0] != 1) {
    iy = ipiv[0] - 1;
    Y[0] = Y[iy];
    Y[iy] = C[0];
  }

  if (ipiv[1] != 2) {
    smax = Y[1];
    iy = ipiv[1] - 1;
    Y[1] = Y[iy];
    Y[iy] = smax;
  }

  if (ipiv[2] != 3) {
    smax = Y[2];
    iy = ipiv[2] - 1;
    Y[2] = Y[iy];
    Y[iy] = smax;
  }

  if (Y[0] != 0.0) {
    for (iy = 2; iy < 5; iy++) {
      Y[iy - 1] -= Y[0] * b_A[iy - 1];
    }
  }

  if (Y[1] != 0.0) {
    for (iy = 3; iy < 5; iy++) {
      Y[iy - 1] -= Y[1] * b_A[iy + 3];
    }
  }

  if (Y[2] != 0.0) {
    for (iy = 4; iy < 5; iy++) {
      Y[3] -= Y[2] * b_A[11];
    }
  }

  if (Y[3] != 0.0) {
    Y[3] /= b_A[15];
    for (iy = 0; iy < 3; iy++) {
      Y[iy] -= Y[3] * b_A[iy + 12];
    }
  }

  if (Y[2] != 0.0) {
    Y[2] /= b_A[10];
    for (iy = 0; iy < 2; iy++) {
      Y[iy] -= Y[2] * b_A[iy + 8];
    }
  }

  if (Y[1] != 0.0) {
    Y[1] /= b_A[5];
    for (iy = 0; iy < 1; iy++) {
      Y[0] -= Y[1] * b_A[4];
    }
  }

  if (Y[0] != 0.0) {
    Y[0] /= b_A[0];
  }

  /*  Test 2: Slicing & 3D matrices */
  for (i = 0; i < 3; i++) {
    iy = i << 2;
    X[iy] = B[iy + 48];
    X[iy + 1] = B[iy + 49];
    X[iy + 2] = B[iy + 50];
    X[iy + 3] = B[iy + 51];
  }

  /*  (4 x 3) */
  /*  Test 3: Multiplication */
  /*  Test 3: Transpose */
  for (i = 0; i < 4; i++) {
    smax = A[i + 4];
    s = A[i + 8];
    d = A[i + 12];
    for (jp1j = 0; jp1j < 4; jp1j++) {
      iy = jp1j << 2;
      W[i + iy] = ((A[i] * (double)b_b[iy] + smax * (double)b_b[iy + 1]) + s *
                   (double)b_b[iy + 2]) + d * (double)b_b[iy + 3];
      V[jp1j + (i << 2)] = C[jp1j] * C[i];
    }
  }

  /*  Transpose */
}

/*
 * File trailer for functionA.c
 *
 * [EOF]
 */
