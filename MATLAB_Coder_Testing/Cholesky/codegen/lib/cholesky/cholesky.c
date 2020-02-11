/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: cholesky.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 20:17:58
 */

/* Include Files */
#include "cholesky.h"
#include <math.h>
#include <string.h>

/* Function Definitions */

/*
 * Returns the Cholesky decomposition of a Hermitian
 *  positive-definite matrix A
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                double R_data[]
 *                int R_size[2]
 * Return Type  : void
 */
void cholesky(const double A_data[], const int A_size[2], double R_data[], int
              R_size[2])
{
  int jmax;
  int n;
  int info;
  int j;
  boolean_T exitg1;
  int idxAjj;
  double ssq;
  int ix;
  int iy;
  int i;
  int idxAjp1j;
  int b_i;
  int iac;
  double c;
  int i1;
  int ia;
  R_size[0] = A_size[0];
  R_size[1] = A_size[1];
  jmax = A_size[0] * A_size[1];
  if (0 <= jmax - 1) {
    memcpy(&R_data[0], &A_data[0], jmax * sizeof(double));
  }

  n = A_size[1];
  if (A_size[1] != 0) {
    info = 0;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j <= n - 1)) {
      idxAjj = j + j * n;
      ssq = 0.0;
      if (j >= 1) {
        ix = j;
        iy = j;
        for (jmax = 0; jmax < j; jmax++) {
          ssq += R_data[ix] * R_data[iy];
          ix += n;
          iy += n;
        }
      }

      ssq = R_data[idxAjj] - ssq;
      if (ssq > 0.0) {
        ssq = sqrt(ssq);
        R_data[idxAjj] = ssq;
        if (j + 1 < n) {
          jmax = (n - j) - 1;
          i = j + 2;
          idxAjp1j = idxAjj + 2;
          if ((jmax != 0) && (j != 0)) {
            ix = j;
            b_i = (j + n * (j - 1)) + 2;
            for (iac = i; n < 0 ? iac >= b_i : iac <= b_i; iac += n) {
              c = -R_data[ix];
              iy = idxAjj + 1;
              i1 = (iac + jmax) - 1;
              for (ia = iac; ia <= i1; ia++) {
                R_data[iy] += R_data[ia - 1] * c;
                iy++;
              }

              ix += n;
            }
          }

          ssq = 1.0 / ssq;
          b_i = idxAjj + jmax;
          for (jmax = idxAjp1j; jmax <= b_i + 1; jmax++) {
            R_data[jmax - 1] *= ssq;
          }
        }

        j++;
      } else {
        R_data[idxAjj] = ssq;
        info = j + 1;
        exitg1 = true;
      }
    }

    if (info == 0) {
      jmax = A_size[1];
    } else {
      jmax = info - 1;
    }

    for (j = 2; j <= jmax; j++) {
      for (i = 0; i <= j - 2; i++) {
        R_data[i + R_size[0] * (j - 1)] = 0.0;
      }
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void cholesky_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void cholesky_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for cholesky.c
 *
 * [EOF]
 */
