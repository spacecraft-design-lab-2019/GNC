/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 16-Apr-2020 13:51:41
 */

/* Include Files */
#include <math.h>
#include <string.h>
#include "milqr.h"

/* Variable Definitions */
static const signed char iv0[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

static const signed char iv1[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
  0, 0, 1 };

static const signed char iv2[49] = { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 1 };

static const float fv0[18] = { 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
  0.0F, 0.01F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F, 0.01F };

static const float fv1[9] = { 0.01F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F,
  0.01F };

static const float fv2[9] = { 0.05F, 0.0F, 0.0F, 0.0F, 0.05F, 0.0F, 0.0F, 0.0F,
  0.05F };

/* Function Declarations */
static void b_abs(const float x[2997], float y[2997]);
static void b_forwardRollout(const float x[7000], const float xg[7], const float
  u[2997], const float l[2997], const float K[17982], float alpha, const float
  u_lims[6], const float B_ECI[3000], float xnew[7000], float unew[2997], float
  fx[35964], float fu[17982], float cx[6000], float cu[2997], float cxx[36000],
  float cuu[8991], float *cost);
static void b_sign(float *x);
static void backwardPass(const float fx[35964], const float fu[17982], const
  float cx[6000], const float cu[2997], const float cxx[36000], const float cuu
  [8991], float lambda, const float u_lims[6], const float u[2997], float l[2997],
  float K[17982], float dV[2], bool *diverged);
static void boxQPsolve(const float Quu[9], const float Qu[3], const float
  lower_lim[3], const float upper_lim[3], const float u0[3], float u[3], float
  *result, float Luu_data[], int Luu_size[2], bool b_free[3]);
static void chol_free(const float A_data[], const int A_size[2], float L_data[],
                      int L_size[2], float *fail);
static void forwardRollout(const float x[7000], const float xg[7], const float
  u[2997], const float K[17982], const float u_lims[6], const float B_ECI[3000],
  float xnew[7000], float unew[2997], float fx[35964], float fu[17982], float
  cx[6000], float cu[2997], float cxx[36000], float cuu[8991], float *cost);
static float mean(const float x[999]);
static void power(float y[11]);
static void satellite_dynamics(const float x[7], const float u[3], const float
  B_ECI[3], float xdot[7], float dxdot[70]);
static void updateLambda(float *lambda, float direction);

/* Function Definitions */

/*
 * Arguments    : const float x[2997]
 *                float y[2997]
 * Return Type  : void
 */
static void b_abs(const float x[2997], float y[2997])
{
  int k;
  for (k = 0; k < 2997; k++) {
    y[k] = fabsf(x[k]);
  }
}

/*
 * Uses an rk method to roll out a trajectory
 *  Returns the new trajectory, cost and the derivatives along the trajectory
 * Arguments    : const float x[7000]
 *                const float xg[7]
 *                const float u[2997]
 *                const float l[2997]
 *                const float K[17982]
 *                float alpha
 *                const float u_lims[6]
 *                const float B_ECI[3000]
 *                float xnew[7000]
 *                float unew[2997]
 *                float fx[35964]
 *                float fu[17982]
 *                float cx[6000]
 *                float cu[2997]
 *                float cxx[36000]
 *                float cuu[8991]
 *                float *cost
 * Return Type  : void
 */
static void b_forwardRollout(const float x[7000], const float xg[7], const float
  u[2997], const float l[2997], const float K[17982], float alpha, const float
  u_lims[6], const float B_ECI[3000], float xnew[7000], float unew[2997], float
  fx[35964], float fu[17982], float cx[6000], float cu[2997], float cxx[36000],
  float cuu[8991], float *cost)
{
  int i11;
  int k;
  int dx_tmp;
  float dx[6];
  int b_dx_tmp;
  int c_dx_tmp;
  float q_idx_0;
  int q_idx_1_tmp;
  int q_idx_2_tmp;
  int q_idx_3_tmp;
  float b_xnew;
  float b_xg;
  float out;
  int c_sign;
  float z1[3];
  float q[16];
  int xnew_tmp;
  int b_k;
  float q_error[4];
  float fv11[9];
  int unew_tmp;
  float d_sign[12];
  float b[7];
  float dxdot1[70];
  float c_xnew[7];
  float dxdot2[70];
  float x1[7];
  float fx_tmp[49];
  float b_fx_tmp[42];
  float fv12[49];
  float d_xnew[42];
  float c_fx_tmp[42];
  int b_xnew_tmp;
  float fv13[21];

  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0], 0, 7000U * sizeof(float));
  memset(&unew[0], 0, 2997U * sizeof(float));
  *cost = 0.0F;
  for (i11 = 0; i11 < 7; i11++) {
    xnew[i11] = x[i11];
  }

  for (k = 0; k < 999; k++) {
    /*  Find the state error vector dx */
    dx_tmp = 7 * k + 4;
    dx[3] = xnew[dx_tmp] - x[dx_tmp];
    b_dx_tmp = 7 * k + 5;
    dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
    c_dx_tmp = 7 * k + 6;
    dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
    q_idx_0 = x[7 * k];
    q_idx_1_tmp = 1 + 7 * k;
    q_idx_2_tmp = 2 + 7 * k;
    q_idx_3_tmp = 3 + 7 * k;

    /*  conjugate */
    /*  Forms the "left matrix" for quaternion multiplication */
    /*  where q*q1 = L(q)q1 */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Input */
    /*  q - the quaternion to build the left matrix from */
    /*  Output */
    /*  L - the left matrix */
    /*  L = [s, -v'; */
    /*       v, sI+skew(v)] */
    /* --------------------------------------------------- */
    q[0] = q_idx_0;
    q[4] = x[q_idx_1_tmp];
    q[8] = x[q_idx_2_tmp];
    q[12] = x[q_idx_3_tmp];
    q[1] = -x[q_idx_1_tmp];
    q[5] = q_idx_0;
    q[9] = x[q_idx_3_tmp];
    q[13] = -x[q_idx_2_tmp];
    q[2] = -x[q_idx_2_tmp];
    q[6] = -x[q_idx_3_tmp];
    q[10] = q_idx_0;
    q[14] = x[q_idx_1_tmp];
    q[3] = -x[q_idx_3_tmp];
    q[7] = x[q_idx_2_tmp];
    q[11] = -x[q_idx_1_tmp];
    q[15] = q_idx_0;
    q_idx_0 = 0.0F;
    for (i11 = 0; i11 < 4; i11++) {
      q_error[i11] = 0.0F;
      b_xg = ((q[i11] * xnew[7 * k] + q[i11 + 4] * xnew[q_idx_1_tmp]) + q[i11 +
              8] * xnew[q_idx_2_tmp]) + q[i11 + 12] * xnew[q_idx_3_tmp];
      q_error[i11] = b_xg;
      q_idx_0 += b_xg * b_xg;
    }

    q_idx_0 = sqrtf(q_idx_0);
    q_error[0] /= q_idx_0;
    q_error[1] /= q_idx_0;
    q_error[2] /= q_idx_0;
    q_error[3] /= q_idx_0;

    /*  re-normalize */
    /*  inverse Cayley Map */
    dx[0] = q_error[1] / q_error[0];
    dx[1] = q_error[2] / q_error[0];
    dx[2] = q_error[3] / q_error[0];

    /*  Find the new control and ensure it is within the limits */
    for (b_k = 0; b_k < 3; b_k++) {
      b_xg = 0.0F;
      for (i11 = 0; i11 < 6; i11++) {
        b_xg += K[(b_k + 3 * i11) + 18 * k] * dx[i11];
      }

      unew_tmp = b_k + 3 * k;
      unew[unew_tmp] = fminf(u_lims[3 + b_k], fmaxf(u_lims[b_k], (u[unew_tmp] -
        alpha * l[unew_tmp]) - b_xg));
    }

    /*  Step the dynamics forward */
    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Step dynamics */
    /*  Explicit midpoint step from x_k to x_{k+1} */
    satellite_dynamics(*(float (*)[7])&xnew[7 * k], *(float (*)[3])&unew[3 * k],
                       *(float (*)[3])&B_ECI[3 * k], b, dxdot1);
    for (i11 = 0; i11 < 7; i11++) {
      c_xnew[i11] = xnew[i11 + 7 * k] + 0.015F * b[i11];
    }

    satellite_dynamics(c_xnew, *(float (*)[3])&unew[3 * k], *(float (*)[3])&
                       B_ECI[3 * k], b, dxdot2);
    for (i11 = 0; i11 < 7; i11++) {
      x1[i11] = xnew[i11 + 7 * k] + 0.03F * b[i11];
    }

    /*  Re-normalize the quaternion */
    q_idx_0 = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] *
                    x1[3]);

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    for (i11 = 0; i11 < 7; i11++) {
      for (b_k = 0; b_k < 7; b_k++) {
        unew_tmp = b_k + 7 * i11;
        fx_tmp[unew_tmp] = 0.00045F * dxdot2[unew_tmp];
      }
    }

    fv11[0] = 0.0F;
    fv11[3] = -x1[3];
    fv11[6] = x1[2];
    fv11[1] = x1[3];
    fv11[4] = 0.0F;
    fv11[7] = -x1[1];
    fv11[2] = -x1[2];
    fv11[5] = x1[1];
    fv11[8] = 0.0F;
    for (i11 = 0; i11 < 3; i11++) {
      b_fx_tmp[i11] = -x1[1 + i11];
      b_fx_tmp[i11 + 3] = 0.0F;
      unew_tmp = 6 * (i11 + 1);
      b_fx_tmp[unew_tmp] = x1[0] * (float)iv0[i11] + fv11[i11];
      b_fx_tmp[unew_tmp + 3] = 0.0F;
      b_fx_tmp[1 + unew_tmp] = x1[0] * (float)iv0[i11 + 3] + fv11[i11 + 3];
      b_fx_tmp[unew_tmp + 4] = 0.0F;
      b_fx_tmp[2 + unew_tmp] = x1[0] * (float)iv0[i11 + 6] + fv11[i11 + 6];
      b_fx_tmp[unew_tmp + 5] = 0.0F;
      for (b_k = 0; b_k < 6; b_k++) {
        b_fx_tmp[b_k + 6 * (i11 + 4)] = iv1[i11 + 3 * b_k];
      }
    }

    for (i11 = 0; i11 < 7; i11++) {
      for (b_k = 0; b_k < 7; b_k++) {
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += fx_tmp[i11 + 7 * xnew_tmp] * dxdot1[xnew_tmp + 7 * b_k];
        }

        xnew_tmp = i11 + 7 * b_k;
        fv12[xnew_tmp] = ((float)iv2[xnew_tmp] + 0.03F * dxdot2[xnew_tmp]) +
          b_xg;
      }
    }

    b_xnew = xnew[7 * k];
    fv11[0] = 0.0F;
    fv11[3] = -xnew[3 + 7 * k];
    fv11[6] = xnew[2 + 7 * k];
    fv11[1] = xnew[3 + 7 * k];
    fv11[4] = 0.0F;
    fv11[7] = -xnew[1 + 7 * k];
    fv11[2] = -xnew[2 + 7 * k];
    fv11[5] = xnew[1 + 7 * k];
    fv11[8] = 0.0F;
    for (i11 = 0; i11 < 6; i11++) {
      for (b_k = 0; b_k < 7; b_k++) {
        unew_tmp = i11 + 6 * b_k;
        c_fx_tmp[unew_tmp] = 0.0F;
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += b_fx_tmp[i11 + 6 * xnew_tmp] * fv12[xnew_tmp + 7 * b_k];
        }

        c_fx_tmp[unew_tmp] = b_xg;
      }
    }

    for (i11 = 0; i11 < 3; i11++) {
      d_xnew[7 * i11] = -xnew[(i11 + 7 * k) + 1];
      xnew_tmp = 7 * (i11 + 3);
      d_xnew[xnew_tmp] = 0.0F;
      d_xnew[7 * i11 + 1] = b_xnew * (float)iv0[3 * i11] + fv11[3 * i11];
      d_xnew[xnew_tmp + 1] = 0.0F;
      b_xnew_tmp = 1 + 3 * i11;
      d_xnew[7 * i11 + 2] = b_xnew * (float)iv0[b_xnew_tmp] + fv11[b_xnew_tmp];
      d_xnew[xnew_tmp + 2] = 0.0F;
      b_xnew_tmp = 2 + 3 * i11;
      d_xnew[7 * i11 + 3] = b_xnew * (float)iv0[b_xnew_tmp] + fv11[b_xnew_tmp];
      d_xnew[xnew_tmp + 3] = 0.0F;
    }

    for (i11 = 0; i11 < 6; i11++) {
      d_xnew[7 * i11 + 4] = iv1[3 * i11];
      d_xnew[7 * i11 + 5] = iv1[1 + 3 * i11];
      d_xnew[7 * i11 + 6] = iv1[2 + 3 * i11];
    }

    for (i11 = 0; i11 < 6; i11++) {
      for (b_k = 0; b_k < 6; b_k++) {
        unew_tmp = (i11 + 6 * b_k) + 36 * k;
        fx[unew_tmp] = 0.0F;
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += c_fx_tmp[i11 + 6 * xnew_tmp] * d_xnew[xnew_tmp + 7 * b_k];
        }

        fx[unew_tmp] = b_xg;
      }
    }

    for (i11 = 0; i11 < 7; i11++) {
      for (b_k = 0; b_k < 3; b_k++) {
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += fx_tmp[i11 + 7 * xnew_tmp] * dxdot1[xnew_tmp + 7 * (7 + b_k)];
        }

        fv13[i11 + 7 * b_k] = 0.03F * dxdot2[i11 + 7 * (7 + b_k)] + b_xg;
      }
    }

    for (i11 = 0; i11 < 6; i11++) {
      for (b_k = 0; b_k < 3; b_k++) {
        unew_tmp = (i11 + 6 * b_k) + 18 * k;
        fu[unew_tmp] = 0.0F;
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += b_fx_tmp[i11 + 6 * xnew_tmp] * fv13[xnew_tmp + 7 * b_k];
        }

        fu[unew_tmp] = b_xg;
      }
    }

    xnew_tmp = 7 * (k + 1);
    xnew[xnew_tmp] = x1[0] / q_idx_0;
    xnew[1 + xnew_tmp] = x1[1] / q_idx_0;
    xnew[2 + xnew_tmp] = x1[2] / q_idx_0;
    xnew[3 + xnew_tmp] = x1[3] / q_idx_0;
    xnew[xnew_tmp + 4] = x1[4];
    xnew[xnew_tmp + 5] = x1[5];
    xnew[xnew_tmp + 6] = x1[6];

    /*  Calculate the cost */
    /*  Calculates the cost contribution of a given state and control */
    /*  Utilizes a standard quadratic cost function for the angular velocity */
    /*  and a geodesic cost for the attitude quaternion */
    /*  Also returns the 2nd order expansion of the cost function */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Inputs */
    /* ===================================== */
    /*  x        - [quaternion; omega] (7x1) */
    /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
    /*  terminal - integer 0 or 1 */
    /* --------------------------------------------------- */
    /*  control hessian */
    /*  cumulative angular velocity cost hessian */
    /*  cumulative geodesic cost weighting */
    /*  Finds the geodesic quaternion-error cost */
    /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
    /*  Also records the sign (+ or -) which minimizes quat_cost */
    /*  this is used when calculating the Jacobain and Hessian */
    q_idx_0 = xnew[7 * k];
    b_xnew = xg[0] * q_idx_0;
    b_xg = ((b_xnew + xg[1] * xnew[q_idx_1_tmp]) + xg[2] * xnew[q_idx_2_tmp]) +
      xg[3] * xnew[q_idx_3_tmp];
    if (1.0F + b_xg < 1.0F - b_xg) {
      out = 1.0F + b_xg;
      c_sign = 1;
    } else {
      out = 1.0F - b_xg;
      c_sign = -1;
    }

    z1[0] = xnew[dx_tmp] - xg[4];
    z1[1] = xnew[b_dx_tmp] - xg[5];
    z1[2] = xnew[c_dx_tmp] - xg[6];

    /*  State cost Hessian */
    b_xg = ((b_xnew + xg[1] * xnew[1 + 7 * k]) + xg[2] * xnew[2 + 7 * k]) + xg[3]
      * xnew[3 + 7 * k];
    for (i11 = 0; i11 < 3; i11++) {
      b_k = 6 * i11 + 36 * k;
      cxx[b_k] = (float)(-c_sign * iv0[3 * i11]) * b_xg;
      unew_tmp = 6 * (i11 + 3) + 36 * k;
      cxx[unew_tmp] = 0.0F;
      cxx[b_k + 1] = (float)(-c_sign * iv0[1 + 3 * i11]) * b_xg;
      cxx[unew_tmp + 1] = 0.0F;
      cxx[b_k + 2] = (float)(-c_sign * iv0[2 + 3 * i11]) * b_xg;
      cxx[unew_tmp + 2] = 0.0F;
    }

    for (i11 = 0; i11 < 6; i11++) {
      b_k = 6 * i11 + 36 * k;
      cxx[b_k + 3] = fv0[3 * i11];
      cxx[b_k + 4] = fv0[1 + 3 * i11];
      cxx[b_k + 5] = fv0[2 + 3 * i11];
    }

    /*  State cost Jacobian */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    fv11[0] = 0.0F;
    fv11[3] = -xnew[3 + 7 * k];
    fv11[6] = xnew[2 + 7 * k];
    fv11[1] = xnew[3 + 7 * k];
    fv11[4] = 0.0F;
    fv11[7] = -xnew[1 + 7 * k];
    fv11[2] = -xnew[2 + 7 * k];
    fv11[5] = xnew[1 + 7 * k];
    fv11[8] = 0.0F;
    for (i11 = 0; i11 < 3; i11++) {
      d_sign[i11] = (float)c_sign * -xnew[(i11 + 7 * k) + 1];
      unew_tmp = 3 * (i11 + 1);
      d_sign[unew_tmp] = (float)c_sign * (q_idx_0 * (float)iv0[i11] + fv11[i11]);
      d_sign[1 + unew_tmp] = (float)c_sign * (q_idx_0 * (float)iv0[i11 + 3] +
        fv11[i11 + 3]);
      d_sign[2 + unew_tmp] = (float)c_sign * (q_idx_0 * (float)iv0[i11 + 6] +
        fv11[i11 + 6]);
    }

    /*  Control cost Hessian & Jacobian */
    b_xg = 0.0F;
    for (i11 = 0; i11 < 3; i11++) {
      unew_tmp = i11 + 6 * k;
      cx[unew_tmp] = ((d_sign[i11] * xg[0] + d_sign[i11 + 3] * xg[1]) +
                      d_sign[i11 + 6] * xg[2]) + d_sign[i11 + 9] * xg[3];
      cx[unew_tmp + 3] = (fv1[i11] * z1[0] + fv1[i11 + 3] * z1[1]) + fv1[i11 + 6]
        * z1[2];
      unew_tmp = i11 + 3 * k;
      cu[unew_tmp] = 0.0F;
      xnew_tmp = 3 * i11 + 9 * k;
      cuu[xnew_tmp] = fv2[3 * i11];
      b_k = 1 + 3 * i11;
      cuu[xnew_tmp + 1] = fv2[b_k];
      b_xnew_tmp = 2 + 3 * i11;
      cuu[xnew_tmp + 2] = fv2[b_xnew_tmp];
      cu[unew_tmp] = (fv2[i11] * unew[3 * k] + fv2[i11 + 3] * unew[1 + 3 * k]) +
        fv2[i11 + 6] * unew[2 + 3 * k];
      b_xg += ((0.5F * z1[0] * fv1[3 * i11] + 0.5F * z1[1] * fv1[b_k]) + 0.5F *
               z1[2] * fv1[b_xnew_tmp]) * z1[i11];
    }

    q_idx_0 = 0.0F;
    for (i11 = 0; i11 < 3; i11++) {
      q_idx_0 += ((0.5F * unew[3 * k] * fv2[3 * i11] + 0.5F * unew[1 + 3 * k] *
                   fv2[1 + 3 * i11]) + 0.5F * unew[2 + 3 * k] * fv2[2 + 3 * i11])
        * unew[i11 + 3 * k];
    }

    *cost += (out + b_xg) + q_idx_0;
  }

  /*  Final cost */
  /*  Calculates the cost contribution of a given state and control */
  /*  Utilizes a standard quadratic cost function for the angular velocity */
  /*  and a geodesic cost for the attitude quaternion */
  /*  Also returns the 2nd order expansion of the cost function */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Inputs */
  /* ===================================== */
  /*  x        - [quaternion; omega] (7x1) */
  /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
  /*  terminal - integer 0 or 1 */
  /* --------------------------------------------------- */
  /*  control hessian */
  /*  terminal angular velocity cost hessian */
  /*  terminal geodesic cost weight */
  /*  Finds the geodesic quaternion-error cost */
  /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
  /*  Also records the sign (+ or -) which minimizes quat_cost */
  /*  this is used when calculating the Jacobain and Hessian */
  b_xnew = xg[0] * xnew[6993];
  b_xg = ((b_xnew + xg[1] * xnew[6994]) + xg[2] * xnew[6995]) + xg[3] * xnew
    [6996];
  if (1.0F + b_xg < 1.0F - b_xg) {
    out = 1.0F + b_xg;
    c_sign = 1;
  } else {
    out = 1.0F - b_xg;
    c_sign = -1;
  }

  z1[0] = xnew[6997] - xg[4];
  z1[1] = xnew[6998] - xg[5];
  z1[2] = xnew[6999] - xg[6];

  /*  State cost Hessian */
  xnew_tmp = -c_sign * 10;
  b_xg = ((b_xnew + xg[1] * xnew[6994]) + xg[2] * xnew[6995]) + xg[3] * xnew
    [6996];
  for (i11 = 0; i11 < 3; i11++) {
    cxx[35964 + 6 * i11] = (float)(xnew_tmp * iv0[3 * i11]) * b_xg;
    b_k = 6 * (i11 + 3);
    cxx[35964 + b_k] = 0.0F;
    cxx[6 * i11 + 35965] = (float)(xnew_tmp * iv0[1 + 3 * i11]) * b_xg;
    cxx[b_k + 35965] = 0.0F;
    cxx[6 * i11 + 35966] = (float)(xnew_tmp * iv0[2 + 3 * i11]) * b_xg;
    cxx[b_k + 35966] = 0.0F;
  }

  for (i11 = 0; i11 < 6; i11++) {
    cxx[6 * i11 + 35967] = iv1[3 * i11];
    cxx[6 * i11 + 35968] = iv1[1 + 3 * i11];
    cxx[6 * i11 + 35969] = iv1[2 + 3 * i11];
  }

  /*  State cost Jacobian */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  xnew_tmp = c_sign * 10;
  fv11[0] = 0.0F;
  fv11[3] = -xnew[6996];
  fv11[6] = xnew[6995];
  fv11[1] = xnew[6996];
  fv11[4] = 0.0F;
  fv11[7] = -xnew[6994];
  fv11[2] = -xnew[6995];
  fv11[5] = xnew[6994];
  fv11[8] = 0.0F;
  for (i11 = 0; i11 < 3; i11++) {
    d_sign[i11] = (float)xnew_tmp * -xnew[i11 + 6994];
    unew_tmp = 3 * (i11 + 1);
    d_sign[unew_tmp] = (float)xnew_tmp * (xnew[6993] * (float)iv0[i11] +
      fv11[i11]);
    d_sign[1 + unew_tmp] = (float)xnew_tmp * (xnew[6993] * (float)iv0[i11 + 3] +
      fv11[i11 + 3]);
    d_sign[2 + unew_tmp] = (float)xnew_tmp * (xnew[6993] * (float)iv0[i11 + 6] +
      fv11[i11 + 6]);
  }

  /*  Control cost Hessian & Jacobian */
  b_xg = 0.0F;
  for (i11 = 0; i11 < 3; i11++) {
    cx[5994 + i11] = ((d_sign[i11] * xg[0] + d_sign[i11 + 3] * xg[1]) +
                      d_sign[i11 + 6] * xg[2]) + d_sign[i11 + 9] * xg[3];
    cx[i11 + 5997] = ((float)iv0[i11] * z1[0] + (float)iv0[i11 + 3] * z1[1]) +
      (float)iv0[i11 + 6] * z1[2];
    b_xg += ((0.5F * z1[0] * (float)iv0[3 * i11] + 0.5F * z1[1] * (float)iv0[1 +
              3 * i11]) + 0.5F * z1[2] * (float)iv0[2 + 3 * i11]) * z1[i11];
  }

  *cost += 10.0F * out + b_xg;
}

/*
 * Arguments    : float *x
 * Return Type  : void
 */
static void b_sign(float *x)
{
  if (*x < 0.0F) {
    *x = -1.0F;
  } else {
    if (*x > 0.0F) {
      *x = 1.0F;
    }
  }
}

/*
 * Perfoms the LQR backward pass to find the optimal controls
 *  Solves a quadratic program (QP) at each timestep for the optimal
 *  controls given the control limits
 * Arguments    : const float fx[35964]
 *                const float fu[17982]
 *                const float cx[6000]
 *                const float cu[2997]
 *                const float cxx[36000]
 *                const float cuu[8991]
 *                float lambda
 *                const float u_lims[6]
 *                const float u[2997]
 *                float l[2997]
 *                float K[17982]
 *                float dV[2]
 *                bool *diverged
 * Return Type  : void
 */
static void backwardPass(const float fx[35964], const float fu[17982], const
  float cx[6000], const float cu[2997], const float cxx[36000], const float cuu
  [8991], float lambda, const float u_lims[6], const float u[2997], float l[2997],
  float K[17982], float dV[2], bool *diverged)
{
  int i3;
  float Vx[6];
  int k;
  int i4;
  bool exitg1;
  int b_k;
  float Vxx[36];
  float a;
  float Qx_tmp[36];
  float Qu_tmp[18];
  float Qu[3];
  int i5;
  float b_u_lims[3];
  float Quu[9];
  float b_Quu[9];
  float c_u_lims[3];
  float Quu_tmp[18];
  float f1;
  float fv10[3];
  int i6;
  int u_lims_tmp_tmp;
  float Qux[18];
  int b_u_lims_tmp_tmp;
  float lk[3];
  float result;
  float Luu_data[9];
  int Luu_size[2];
  bool b_free[3];
  float Kk[18];
  bool out;
  bool exitg2;
  signed char tmp_data[3];
  float b_Kk[18];
  float b_cx[6];
  float b_Qx_tmp[36];
  float b_cxx[36];
  int value;
  int n;
  int Y_size_idx_0;
  int loop_ub;
  float Y_data[18];
  float coeffs_data[3];
  int a_size_idx_1;
  int i;
  int b_loop_ub;
  int c_loop_ub;
  int X_size_idx_0;
  float a_data[3];
  int j;
  int LT_size_idx_0;
  float b_data[3];
  float LT_data[9];
  int d_loop_ub;
  int e_loop_ub;
  int b_i;
  signed char b_tmp_data[3];

  /*  Initialize matrices (for C code, not needed in MATLAB) */
  memset(&l[0], 0, 2997U * sizeof(float));
  memset(&K[0], 0, 17982U * sizeof(float));

  /*  Change in cost */
  dV[0] = 0.0F;
  dV[1] = 0.0F;

  /*  Set cost-to-go Jacobian and Hessian equal to final costs */
  for (i3 = 0; i3 < 6; i3++) {
    Vx[i3] = cx[5994 + i3];
    for (i4 = 0; i4 < 6; i4++) {
      b_k = i4 + 6 * i3;
      Vxx[b_k] = cxx[35964 + b_k];
    }
  }

  *diverged = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 999)) {
    /*  Define cost gradients and hessians */
    for (i3 = 0; i3 < 6; i3++) {
      for (i4 = 0; i4 < 6; i4++) {
        Qx_tmp[i4 + 6 * i3] = fx[(i3 + 6 * i4) + 36 * (998 - k)];
      }

      b_k = i3 + 18 * (998 - k);
      Qu_tmp[3 * i3] = fu[b_k];
      Qu_tmp[1 + 3 * i3] = fu[b_k + 6];
      Qu_tmp[2 + 3 * i3] = fu[b_k + 12];
    }

    for (i3 = 0; i3 < 3; i3++) {
      a = 0.0F;
      for (i4 = 0; i4 < 6; i4++) {
        i5 = i3 + 3 * i4;
        a += Qu_tmp[i5] * Vx[i4];
        Quu_tmp[i5] = 0.0F;
        f1 = 0.0F;
        for (i6 = 0; i6 < 6; i6++) {
          f1 += Qu_tmp[i3 + 3 * i6] * Vxx[i6 + 6 * i4];
        }

        Quu_tmp[i5] = f1;
      }

      Qu[i3] = cu[i3 + 3 * (998 - k)] + a;
      for (i4 = 0; i4 < 3; i4++) {
        a = 0.0F;
        for (i5 = 0; i5 < 6; i5++) {
          a += Quu_tmp[i3 + 3 * i5] * fu[(i5 + 6 * i4) + 18 * (998 - k)];
        }

        b_k = i3 + 3 * i4;
        b_Quu[b_k] = cuu[b_k + 9 * (998 - k)] + a;
      }

      for (i4 = 0; i4 < 6; i4++) {
        b_k = i3 + 3 * i4;
        Qux[b_k] = 0.0F;
        a = 0.0F;
        for (i5 = 0; i5 < 6; i5++) {
          a += Quu_tmp[i3 + 3 * i5] * fx[(i5 + 6 * i4) + 36 * (998 - k)];
        }

        Qux[b_k] = a;
      }
    }

    /*  Regularization (for Cholesky positive definiteness) */
    /*  Solve the Quadratic program with control limits */
    for (i3 = 0; i3 < 9; i3++) {
      Quu[i3] = b_Quu[i3] + (float)iv0[i3] * lambda;
    }

    b_u_lims[0] = u_lims[0] - u[3 * (998 - k)];
    c_u_lims[0] = u_lims[3] - u[3 * (998 - k)];
    i3 = 3 * ((int)fminf(999.0F, (999.0F + -(float)k) + 1.0F) - 1);
    fv10[0] = -l[i3];
    u_lims_tmp_tmp = 1 + 3 * (998 - k);
    b_u_lims[1] = u_lims[1] - u[u_lims_tmp_tmp];
    c_u_lims[1] = u_lims[4] - u[u_lims_tmp_tmp];
    fv10[1] = -l[1 + i3];
    b_u_lims_tmp_tmp = 2 + 3 * (998 - k);
    b_u_lims[2] = u_lims[2] - u[b_u_lims_tmp_tmp];
    c_u_lims[2] = u_lims[5] - u[b_u_lims_tmp_tmp];
    fv10[2] = -l[2 + i3];
    boxQPsolve(Quu, Qu, b_u_lims, c_u_lims, fv10, lk, &result, Luu_data,
               Luu_size, b_free);
    if (result < 2.0F) {
      *diverged = true;

      /*  fprintf('\nDiverged with lambda = %f\n',lambda); */
      exitg1 = true;
    } else {
      /*  Solve for feedback gains in non-clamped rows of u */
      /*  (using cholesky factor of Quu) */
      memset(&Kk[0], 0, 18U * sizeof(float));
      out = false;
      b_k = 0;
      exitg2 = false;
      while ((!exitg2) && (b_k < 3)) {
        if (b_free[b_k]) {
          out = true;
          exitg2 = true;
        } else {
          b_k++;
        }
      }

      if (out) {
        b_k = 0;
        if (b_free[0]) {
          tmp_data[0] = 1;
          b_k = 1;
        }

        if (b_free[1]) {
          tmp_data[b_k] = 2;
          b_k++;
        }

        if (b_free[2]) {
          tmp_data[b_k] = 3;
        }

        /*  Solves the linear system AX = B for the unknown X using the Cholesky */
        /*  decomposition L of the matrix A. Where LL' = A */
        /*  X can be a vector or a matrix of size n x m */
        /*  Copyright (c) 2020 Robotic Exploration Lab */
        /*   */
        /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
        /*  of this software and associated documentation files (the "Software"), to deal */
        /*  in the Software without restriction, including without limitation the rights */
        /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
        /*  copies of the Software, and to permit persons to whom the Software is */
        /*  furnished to do so, subject to the following conditions: */
        /*   */
        /*  The above copyright notice and this permission notice shall be included in all */
        /*  copies or substantial portions of the Software. */
        /*   */
        /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
        /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
        /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
        /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
        /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
        /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
        /*  SOFTWARE. */
        /*  Solution is found in O(nm) time using back-substitution */
        /*  This implementation only works for lower triangular factorisations */
        value = Luu_size[0];
        n = Luu_size[0];

        /*  Check sizes match and L is lower-triangular */
        /*  Solve LY = B for Y */
        /* ======================= */
        Y_size_idx_0 = Luu_size[0];
        loop_ub = Luu_size[0] * 6;
        if (0 <= loop_ub - 1) {
          memset(&Y_data[0], 0, (unsigned int)(loop_ub * (int)sizeof(float)));
        }

        if (0 <= Luu_size[0] - 1) {
          memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                  (float)));
        }

        i3 = Luu_size[0];
        if (0 <= Luu_size[0] - 1) {
          a_size_idx_1 = Luu_size[0];
          b_loop_ub = Luu_size[0];
          c_loop_ub = Luu_size[0];
        }

        for (i = 0; i < i3; i++) {
          if (1 + i != 1) {
            for (i4 = 0; i4 < i; i4++) {
              coeffs_data[i4] = Luu_data[i + Luu_size[0] * i4];
            }
          }

          if (0 <= b_loop_ub - 1) {
            memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(b_loop_ub * (int)
                    sizeof(float)));
          }

          for (j = 0; j < 6; j++) {
            for (i4 = 0; i4 < c_loop_ub; i4++) {
              b_data[i4] = Y_data[i4 + Y_size_idx_0 * j];
            }

            if ((a_size_idx_1 == 1) || (Y_size_idx_0 == 1)) {
              a = 0.0F;
              for (i4 = 0; i4 < a_size_idx_1; i4++) {
                a += a_data[i4] * b_data[i4];
              }
            } else {
              a = 0.0F;
              for (i4 = 0; i4 < a_size_idx_1; i4++) {
                a += a_data[i4] * b_data[i4];
              }
            }

            Y_data[i + Y_size_idx_0 * j] = (Qux[(tmp_data[i] + 3 * j) - 1] - a) /
              Luu_data[i + Luu_size[0] * i];
          }
        }

        /*  Solve L'X = Y for X */
        /* ======================= */
        X_size_idx_0 = Luu_size[0];
        loop_ub = Luu_size[0] * 6;
        if (0 <= loop_ub - 1) {
          memset(&Qu_tmp[0], 0, (unsigned int)(loop_ub * (int)sizeof(float)));
        }

        LT_size_idx_0 = Luu_size[1];
        loop_ub = Luu_size[0];
        for (i3 = 0; i3 < loop_ub; i3++) {
          b_k = Luu_size[1];
          for (i4 = 0; i4 < b_k; i4++) {
            LT_data[i4 + LT_size_idx_0 * i3] = Luu_data[i3 + Luu_size[0] * i4];
          }
        }

        if (0 <= Luu_size[0] - 1) {
          memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                  (float)));
        }

        i3 = (int)((1.0F + (-1.0F - (float)Luu_size[0])) / -1.0F);
        if (0 <= i3 - 1) {
          a_size_idx_1 = Luu_size[0];
          d_loop_ub = Luu_size[0];
          e_loop_ub = Luu_size[0];
        }

        for (i = 0; i < i3; i++) {
          b_i = (value - i) - 1;
          if (b_i + 1 != value) {
            if (b_i + 2 > n) {
              i4 = -1;
              i5 = 0;
              i6 = -1;
            } else {
              i4 = b_i;
              i5 = value;
              i6 = b_i;
            }

            loop_ub = i5 - i4;
            for (i5 = 0; i5 <= loop_ub - 2; i5++) {
              coeffs_data[(i6 + i5) + 1] = LT_data[b_i + LT_size_idx_0 * ((i4 +
                i5) + 1)];
            }
          }

          if (0 <= d_loop_ub - 1) {
            memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(d_loop_ub * (int)
                    sizeof(float)));
          }

          for (j = 0; j < 6; j++) {
            for (i4 = 0; i4 < e_loop_ub; i4++) {
              b_data[i4] = Qu_tmp[i4 + X_size_idx_0 * j];
            }

            if ((a_size_idx_1 == 1) || (X_size_idx_0 == 1)) {
              a = 0.0F;
              for (i4 = 0; i4 < a_size_idx_1; i4++) {
                a += a_data[i4] * b_data[i4];
              }
            } else {
              a = 0.0F;
              for (i4 = 0; i4 < a_size_idx_1; i4++) {
                a += a_data[i4] * b_data[i4];
              }
            }

            Qu_tmp[b_i + X_size_idx_0 * j] = (Y_data[b_i + Y_size_idx_0 * j] - a)
              / LT_data[b_i + LT_size_idx_0 * b_i];
          }
        }

        b_k = 0;
        if (b_free[0]) {
          b_tmp_data[0] = 1;
          b_k = 1;
        }

        if (b_free[1]) {
          b_tmp_data[b_k] = 2;
          b_k++;
        }

        if (b_free[2]) {
          b_tmp_data[b_k] = 3;
        }

        loop_ub = Luu_size[0];
        for (i3 = 0; i3 < 6; i3++) {
          for (i4 = 0; i4 < loop_ub; i4++) {
            Kk[(b_tmp_data[i4] + 3 * i3) - 1] = -Qu_tmp[i4 + X_size_idx_0 * i3];
          }
        }
      }

      /*  Update Cost to Go Jacobian and Hessian */
      for (i3 = 0; i3 < 3; i3++) {
        for (i4 = 0; i4 < 6; i4++) {
          b_Kk[i4 + 6 * i3] = Kk[i3 + 3 * i4];
        }
      }

      memcpy(&Qu_tmp[0], &b_Kk[0], 18U * sizeof(float));
      for (i3 = 0; i3 < 6; i3++) {
        for (i4 = 0; i4 < 3; i4++) {
          b_k = i3 + 6 * i4;
          b_Kk[b_k] = 0.0F;
          b_Kk[b_k] = (Qu_tmp[i3] * b_Quu[3 * i4] + Qu_tmp[i3 + 6] * b_Quu[1 + 3
                       * i4]) + Qu_tmp[i3 + 12] * b_Quu[2 + 3 * i4];
        }
      }

      memcpy(&Quu_tmp[0], &b_Kk[0], 18U * sizeof(float));
      for (i3 = 0; i3 < 3; i3++) {
        for (i4 = 0; i4 < 6; i4++) {
          b_Kk[i4 + 6 * i3] = Qux[i3 + 3 * i4];
        }
      }

      for (i3 = 0; i3 < 6; i3++) {
        a = 0.0F;
        for (i4 = 0; i4 < 6; i4++) {
          a += Qx_tmp[i3 + 6 * i4] * Vx[i4];
        }

        b_cx[i3] = ((cx[i3 + 6 * (998 - k)] + a) + ((Quu_tmp[i3] * lk[0] +
          Quu_tmp[i3 + 6] * lk[1]) + Quu_tmp[i3 + 12] * lk[2])) + ((Qu_tmp[i3] *
          Qu[0] + Qu_tmp[i3 + 6] * Qu[1]) + Qu_tmp[i3 + 12] * Qu[2]);
      }

      for (i3 = 0; i3 < 6; i3++) {
        Vx[i3] = b_cx[i3] + ((b_Kk[i3] * lk[0] + b_Kk[i3 + 6] * lk[1]) + b_Kk[i3
                             + 12] * lk[2]);
        for (i4 = 0; i4 < 6; i4++) {
          b_k = i3 + 6 * i4;
          b_Qx_tmp[b_k] = 0.0F;
          a = 0.0F;
          for (i5 = 0; i5 < 6; i5++) {
            a += Qx_tmp[i3 + 6 * i5] * Vxx[i5 + 6 * i4];
          }

          b_Qx_tmp[b_k] = a;
        }

        for (i4 = 0; i4 < 6; i4++) {
          a = 0.0F;
          for (i5 = 0; i5 < 6; i5++) {
            a += b_Qx_tmp[i3 + 6 * i5] * fx[(i5 + 6 * i4) + 36 * (998 - k)];
          }

          b_k = i3 + 6 * i4;
          b_cxx[b_k] = cxx[b_k + 36 * (998 - k)] + a;
        }
      }

      for (i3 = 0; i3 < 6; i3++) {
        for (i4 = 0; i4 < 6; i4++) {
          i5 = 1 + 3 * i4;
          i6 = 2 + 3 * i4;
          b_k = i3 + 6 * i4;
          b_Qx_tmp[b_k] = (b_cxx[b_k] + ((Quu_tmp[i3] * Kk[3 * i4] + Quu_tmp[i3
            + 6] * Kk[i5]) + Quu_tmp[i3 + 12] * Kk[i6])) + ((Qu_tmp[i3] * Qux[3 *
            i4] + Qu_tmp[i3 + 6] * Qux[i5]) + Qu_tmp[i3 + 12] * Qux[i6]);
        }
      }

      for (i3 = 0; i3 < 6; i3++) {
        for (i4 = 0; i4 < 6; i4++) {
          b_k = i3 + 6 * i4;
          Qx_tmp[b_k] = 0.0F;
          Qx_tmp[b_k] = (b_Kk[i3] * Kk[3 * i4] + b_Kk[i3 + 6] * Kk[1 + 3 * i4])
            + b_Kk[i3 + 12] * Kk[2 + 3 * i4];
        }
      }

      for (i3 = 0; i3 < 36; i3++) {
        Vxx[i3] = b_Qx_tmp[i3] + Qx_tmp[i3];
      }

      for (i3 = 0; i3 < 6; i3++) {
        for (i4 = 0; i4 < 6; i4++) {
          b_k = i4 + 6 * i3;
          Qx_tmp[b_k] = 0.5F * (Vxx[b_k] + Vxx[i3 + 6 * i4]);
        }
      }

      memcpy(&Vxx[0], &Qx_tmp[0], 36U * sizeof(float));

      /*  Ensure Hessian is symmetric */
      /*  Record control cost change to check convergence */
      a = 0.0F;
      for (i3 = 0; i3 < 3; i3++) {
        a += ((0.5F * lk[0] * b_Quu[3 * i3] + 0.5F * lk[1] * b_Quu[1 + 3 * i3])
              + 0.5F * lk[2] * b_Quu[2 + 3 * i3]) * lk[i3];
      }

      dV[0] += (lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2];
      dV[1] += a;

      /*  Update Control Vectors */
      l[3 * (998 - k)] = -lk[0];
      l[u_lims_tmp_tmp] = -lk[1];
      l[b_u_lims_tmp_tmp] = -lk[2];
      for (i3 = 0; i3 < 6; i3++) {
        b_k = 3 * i3 + 18 * (998 - k);
        K[b_k] = -Kk[3 * i3];
        K[b_k + 1] = -Kk[1 + 3 * i3];
        K[b_k + 2] = -Kk[2 + 3 * i3];
      }

      k++;
    }
  }
}

/*
 * Finds the optimal control within limits to minimize a quadratic cost
 *  Minimizes 0.5*u'*Quu*u + u'*Qu  s.t. lower_lim <= u <= upper_lim
 *
 *  Inputs:
 *  ==========================================
 *  Quu       - control cost Hessian (positive definite)   (m, m)
 *  Qu        - control cost Jacobian                      (m)
 *  lower     - control lower limit                        (m)
 *  upper     - control upper limit                        (m)
 *  u0        - initial control input for warm-start       (m)
 *
 *  Outputs:
 *  =====================================
 *  u         - optimal feed-forward control                (m)
 *  result    - gives exit criterion (explained below)
 *  Luu       - cholesky factor                             (m, m)
 *  free      - set of free dimensions                      (m)
 * Arguments    : const float Quu[9]
 *                const float Qu[3]
 *                const float lower_lim[3]
 *                const float upper_lim[3]
 *                const float u0[3]
 *                float u[3]
 *                float *result
 *                float Luu_data[]
 *                int Luu_size[2]
 *                bool b_free[3]
 * Return Type  : void
 */
static void boxQPsolve(const float Quu[9], const float Qu[3], const float
  lower_lim[3], const float upper_lim[3], const float u0[3], float u[3], float
  *result, float Luu_data[], int Luu_size[2], bool b_free[3])
{
  bool clamped[3];
  int i7;
  float old_cost;
  float scale;
  float z1[3];
  float absxk;
  float t;
  float cost;
  int iter;
  int b_iter;
  bool exitg1;
  int i;
  bool out;
  float grad[3];
  int k;
  bool prev_clamped[3];
  bool exitg2;
  bool b_clamped;
  bool factorize;
  bool guard1 = false;
  int trueCount;
  signed char tmp_data[3];
  signed char b_tmp_data[3];
  int Quu_size[2];
  float Quu_data[9];
  float indef;
  int i8;
  float y;
  float grad_norm;
  float grad_clamped[3];
  signed char c_tmp_data[3];
  int value;
  int n;
  float Y_data[3];
  float coeffs_data[3];
  int a_size_idx_1;
  int loop_ub;
  int b_loop_ub;
  float a_data[3];
  int LT_size_idx_0;
  int c_loop_ub;
  float b_data[3];
  float LT_data[9];
  int d_loop_ub;
  int e_loop_ub;
  int b_i;
  float delta_u;
  float b_delta_u[3];
  float step;
  float u_c_idx_0;
  float u_c_idx_1;
  float u_c_idx_2;
  float cost_c;

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Results */
  /*  =========================== */
  /*   0: No descent direction found */
  /*   1: Hessian is not positive definite */
  /*   2: Maximum main iterations exceeded        */
  /*   3: Maximum line-search iterations exceeded  */
  /*  */
  /*   4: Cost reduction smaller than tolerance      */
  /*   5: Gradient smaller than tolerance     */
  /*   6: All controls are clamped  */
  /*  Initialize arrays */
  /*  Indicies of clamped controls */
  clamped[0] = false;
  b_free[0] = true;
  clamped[1] = false;
  b_free[1] = true;
  clamped[2] = false;
  b_free[2] = true;
  Luu_size[0] = 3;
  Luu_size[1] = 3;
  for (i7 = 0; i7 < 9; i7++) {
    Luu_data[i7] = 0.0F;
  }

  /*  Placeholder to return if Luu not assigned */
  /*  Initialize scalars */
  old_cost = 0.0F;
  *result = 0.0F;

  /*  Solver Options */
  /*  max iterations */
  /*  min norm of non-clamped gradient */
  /*  min relative improvement */
  /*  factor for decreasing stepsize */
  /*  min stepsize for linesearch */
  /*  Armijo tolerance (fraction of linear improvement required) */
  /*  Initial controls */
  /*  Returns array x with all values clamped between lower and upper */
  /*  Initial cost value */
  scale = fmaxf(lower_lim[0], fminf(upper_lim[0], u0[0]));
  z1[0] = scale;
  u[0] = scale;
  absxk = scale * Qu[0];
  scale = fmaxf(lower_lim[1], fminf(upper_lim[1], u0[1]));
  z1[1] = scale;
  u[1] = scale;
  absxk += scale * Qu[1];
  scale = fmaxf(lower_lim[2], fminf(upper_lim[2], u0[2]));
  z1[2] = scale;
  u[2] = scale;
  absxk += scale * Qu[2];
  t = 0.0F;
  for (i7 = 0; i7 < 3; i7++) {
    t += ((0.5F * z1[0] * Quu[3 * i7] + 0.5F * z1[1] * Quu[1 + 3 * i7]) + 0.5F *
          scale * Quu[2 + 3 * i7]) * z1[i7];
  }

  cost = absxk + t;

  /*  Start optimisation */
  iter = 1;
  b_iter = 1;
  exitg1 = false;
  while ((!exitg1) && (b_iter - 1 < 100)) {
    iter = b_iter;
    if (*result != 0.0F) {
      exitg1 = true;
    } else {
      /*  Check relative cost change for convergence */
      if ((b_iter > 1) && (old_cost - cost < 1.0E-8F * fabsf(old_cost))) {
        *result = 4.0F;
        exitg1 = true;
      } else {
        old_cost = cost;

        /*  Gradient of cost function */
        /*  Find clamped controls */
        for (i = 0; i < 3; i++) {
          scale = Qu[i] + ((Quu[i] * u[0] + Quu[i + 3] * u[1]) + Quu[i + 6] * u
                           [2]);
          grad[i] = scale;
          prev_clamped[i] = clamped[i];
          out = false;
          b_clamped = false;
          if ((u[i] == lower_lim[i]) && (scale > 0.0F)) {
            out = true;
            b_clamped = true;
          }

          if ((u[i] == upper_lim[i]) && (scale < 0.0F)) {
            out = true;
            b_clamped = true;
          }

          b_free[i] = !out;
          clamped[i] = b_clamped;
        }

        /*  Check if all controls clamped */
        out = true;
        k = 0;
        exitg2 = false;
        while ((!exitg2) && (k < 3)) {
          if (!clamped[k]) {
            out = false;
            exitg2 = true;
          } else {
            k++;
          }
        }

        if (out) {
          *result = 6.0F;
          exitg1 = true;
        } else {
          /*  Cholesky factorize if clamped controls have changed */
          if (b_iter == 1) {
            factorize = true;
          } else {
            factorize = false;
            k = 0;
            exitg2 = false;
            while ((!exitg2) && (k < 3)) {
              if (prev_clamped[k] != clamped[k]) {
                factorize = true;
                exitg2 = true;
              } else {
                k++;
              }
            }
          }

          /*  Cholesky (check for non-PD) */
          guard1 = false;
          if (factorize) {
            trueCount = 0;
            if (b_free[0]) {
              trueCount = 1;
            }

            if (b_free[1]) {
              trueCount++;
            }

            if (b_free[2]) {
              trueCount++;
            }

            k = 0;
            if (b_free[0]) {
              tmp_data[0] = 1;
              k = 1;
            }

            if (b_free[1]) {
              tmp_data[k] = 2;
              k++;
            }

            if (b_free[2]) {
              tmp_data[k] = 3;
            }

            Quu_size[0] = trueCount;
            Quu_size[1] = trueCount;
            for (i7 = 0; i7 < trueCount; i7++) {
              for (i8 = 0; i8 < trueCount; i8++) {
                Quu_data[i8 + trueCount * i7] = Quu[(tmp_data[i8] + 3 *
                  (tmp_data[i7] - 1)) - 1];
              }
            }

            chol_free(Quu_data, Quu_size, Luu_data, Luu_size, &indef);
            if (indef != 0.0F) {
              *result = 1.0F;
              exitg1 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            /*  check gradient-norm of free controls */
            trueCount = 0;
            if (b_free[0]) {
              trueCount = 1;
            }

            if (b_free[1]) {
              trueCount++;
            }

            if (b_free[2]) {
              trueCount++;
            }

            k = 0;
            if (b_free[0]) {
              b_tmp_data[0] = 1;
              k = 1;
            }

            if (b_free[1]) {
              b_tmp_data[k] = 2;
              k++;
            }

            if (b_free[2]) {
              b_tmp_data[k] = 3;
            }

            if (trueCount == 0) {
              grad_norm = 0.0F;
            } else {
              y = 0.0F;
              if (trueCount == 1) {
                grad_norm = fabsf(grad[b_tmp_data[0] - 1]);
              } else {
                scale = 1.29246971E-26F;
                for (k = 0; k < trueCount; k++) {
                  absxk = fabsf(grad[b_tmp_data[k] - 1]);
                  if (absxk > scale) {
                    t = scale / absxk;
                    y = 1.0F + y * t * t;
                    scale = absxk;
                  } else {
                    t = absxk / scale;
                    y += t * t;
                  }
                }

                grad_norm = scale * sqrtf(y);
              }
            }

            if (grad_norm < 1.0E-8F) {
              *result = 5.0F;
              exitg1 = true;
            } else {
              /*  get search direction */
              scale = u[0] * (float)clamped[0];
              absxk = u[1] * (float)clamped[1];
              t = u[2] * (float)clamped[2];
              k = 0;
              for (i = 0; i < 3; i++) {
                grad_clamped[i] = Qu[i] + ((Quu[i] * scale + Quu[i + 3] * absxk)
                  + Quu[i + 6] * t);
                if (b_free[i]) {
                  c_tmp_data[k] = (signed char)(i + 1);
                  k++;
                }
              }

              /*  Solves the linear system AX = B for the unknown X using the Cholesky */
              /*  decomposition L of the matrix A. Where LL' = A */
              /*  X can be a vector or a matrix of size n x m */
              /*  Copyright (c) 2020 Robotic Exploration Lab */
              /*   */
              /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
              /*  of this software and associated documentation files (the "Software"), to deal */
              /*  in the Software without restriction, including without limitation the rights */
              /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
              /*  copies of the Software, and to permit persons to whom the Software is */
              /*  furnished to do so, subject to the following conditions: */
              /*   */
              /*  The above copyright notice and this permission notice shall be included in all */
              /*  copies or substantial portions of the Software. */
              /*   */
              /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
              /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
              /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
              /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
              /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
              /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
              /*  SOFTWARE. */
              /*  Solution is found in O(nm) time using back-substitution */
              /*  This implementation only works for lower triangular factorisations */
              value = Luu_size[0];
              n = Luu_size[0];

              /*  Check sizes match and L is lower-triangular */
              /*  Solve LY = B for Y */
              /* ======================= */
              if (0 <= Luu_size[0] - 1) {
                memset(&Y_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                        (float)));
              }

              if (0 <= Luu_size[0] - 1) {
                memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)
                        sizeof(float)));
              }

              i7 = Luu_size[0];
              if (0 <= Luu_size[0] - 1) {
                a_size_idx_1 = Luu_size[0];
                loop_ub = Luu_size[0];
                b_loop_ub = Luu_size[0];
              }

              for (i = 0; i < i7; i++) {
                if (1 + i != 1) {
                  for (i8 = 0; i8 < i; i8++) {
                    coeffs_data[i8] = Luu_data[i + Luu_size[0] * i8];
                  }
                }

                if (0 <= loop_ub - 1) {
                  memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(loop_ub *
                          (int)sizeof(float)));
                }

                if (0 <= b_loop_ub - 1) {
                  memcpy(&b_data[0], &Y_data[0], (unsigned int)(b_loop_ub * (int)
                          sizeof(float)));
                }

                if ((a_size_idx_1 == 1) || (Luu_size[0] == 1)) {
                  scale = 0.0F;
                  for (i8 = 0; i8 < a_size_idx_1; i8++) {
                    scale += a_data[i8] * b_data[i8];
                  }
                } else {
                  scale = 0.0F;
                  for (i8 = 0; i8 < a_size_idx_1; i8++) {
                    scale += a_data[i8] * b_data[i8];
                  }
                }

                Y_data[i] = (grad_clamped[c_tmp_data[i] - 1] - scale) /
                  Luu_data[i + Luu_size[0] * i];
              }

              /*  Solve L'X = Y for X */
              /* ======================= */
              if (0 <= Luu_size[0] - 1) {
                memset(&z1[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof(float)));
              }

              LT_size_idx_0 = Luu_size[1];
              c_loop_ub = Luu_size[0];
              for (i7 = 0; i7 < c_loop_ub; i7++) {
                k = Luu_size[1];
                for (i8 = 0; i8 < k; i8++) {
                  LT_data[i8 + LT_size_idx_0 * i7] = Luu_data[i7 + Luu_size[0] *
                    i8];
                }
              }

              if (0 <= Luu_size[0] - 1) {
                memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)
                        sizeof(float)));
              }

              i7 = (int)((1.0F + (-1.0F - (float)Luu_size[0])) / -1.0F);
              if (0 <= i7 - 1) {
                a_size_idx_1 = Luu_size[0];
                d_loop_ub = Luu_size[0];
                e_loop_ub = Luu_size[0];
              }

              for (i = 0; i < i7; i++) {
                b_i = (value - i) - 1;
                if (b_i + 1 != value) {
                  if (b_i + 2 > n) {
                    i8 = -1;
                    k = 0;
                    trueCount = -1;
                  } else {
                    i8 = b_i;
                    k = value;
                    trueCount = b_i;
                  }

                  c_loop_ub = k - i8;
                  for (k = 0; k <= c_loop_ub - 2; k++) {
                    coeffs_data[(trueCount + k) + 1] = LT_data[b_i +
                      LT_size_idx_0 * ((i8 + k) + 1)];
                  }
                }

                if (0 <= d_loop_ub - 1) {
                  memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(d_loop_ub *
                          (int)sizeof(float)));
                }

                if (0 <= e_loop_ub - 1) {
                  memcpy(&b_data[0], &z1[0], (unsigned int)(e_loop_ub * (int)
                          sizeof(float)));
                }

                if ((a_size_idx_1 == 1) || (Luu_size[0] == 1)) {
                  scale = 0.0F;
                  for (i8 = 0; i8 < a_size_idx_1; i8++) {
                    scale += a_data[i8] * b_data[i8];
                  }
                } else {
                  scale = 0.0F;
                  for (i8 = 0; i8 < a_size_idx_1; i8++) {
                    scale += a_data[i8] * b_data[i8];
                  }
                }

                z1[b_i] = (Y_data[b_i] - scale) / LT_data[b_i + LT_size_idx_0 *
                  b_i];
              }

              c_loop_ub = Luu_size[0];
              for (i7 = 0; i7 < c_loop_ub; i7++) {
                b_data[i7] = -z1[i7] - u[c_tmp_data[i7] - 1];
              }

              k = 0;

              /*  cholesky solver */
              /*  check projected change in cost is a reduction */
              delta_u = 0.0F;
              if (b_free[0]) {
                delta_u = b_data[0];
                k = 1;
              }

              z1[0] = delta_u * grad[0];
              b_delta_u[0] = delta_u;
              delta_u = 0.0F;
              if (b_free[1]) {
                delta_u = b_data[k];
                k++;
              }

              z1[1] = delta_u * grad[1];
              b_delta_u[1] = delta_u;
              delta_u = 0.0F;
              if (b_free[2]) {
                delta_u = b_data[k];
              }

              y = (z1[0] + z1[1]) + delta_u * grad[2];
              if (y >= 0.0F) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0F;

                /*  Returns array x with all values clamped between lower and upper */
                scale = fmaxf(lower_lim[0], fminf(upper_lim[0], u[0] +
                  b_delta_u[0]));
                z1[0] = scale;
                u_c_idx_0 = scale;
                absxk = scale * Qu[0];
                scale = fmaxf(lower_lim[1], fminf(upper_lim[1], u[1] +
                  b_delta_u[1]));
                z1[1] = scale;
                u_c_idx_1 = scale;
                absxk += scale * Qu[1];
                scale = fmaxf(lower_lim[2], fminf(upper_lim[2], u[2] + delta_u));
                z1[2] = scale;
                u_c_idx_2 = scale;
                absxk += scale * Qu[2];
                t = 0.0F;
                for (i7 = 0; i7 < 3; i7++) {
                  t += ((0.5F * z1[0] * Quu[3 * i7] + 0.5F * z1[1] * Quu[1 + 3 *
                         i7]) + 0.5F * scale * Quu[2 + 3 * i7]) * z1[i7];
                }

                cost_c = absxk + t;
                exitg2 = false;
                while ((!exitg2) && ((cost_c - cost) / (step * y) < 0.1F)) {
                  step *= 0.6F;

                  /*  Returns array x with all values clamped between lower and upper */
                  scale = fmaxf(lower_lim[0], fminf(upper_lim[0], u[0] + step *
                    b_delta_u[0]));
                  z1[0] = scale;
                  u_c_idx_0 = scale;
                  absxk = scale * Qu[0];
                  scale = fmaxf(lower_lim[1], fminf(upper_lim[1], u[1] + step *
                    b_delta_u[1]));
                  z1[1] = scale;
                  u_c_idx_1 = scale;
                  absxk += scale * Qu[1];
                  scale = fmaxf(lower_lim[2], fminf(upper_lim[2], u[2] + step *
                    delta_u));
                  z1[2] = scale;
                  u_c_idx_2 = scale;
                  absxk += scale * Qu[2];
                  t = 0.0F;
                  for (i7 = 0; i7 < 3; i7++) {
                    t += ((0.5F * z1[0] * Quu[3 * i7] + 0.5F * z1[1] * Quu[1 + 3
                           * i7]) + 0.5F * scale * Quu[2 + 3 * i7]) * z1[i7];
                  }

                  cost_c = absxk + t;
                  if (step < 1.0E-20F) {
                    *result = 3.0F;
                    exitg2 = true;
                  }
                }

                /*  Update the control */
                u[0] = u_c_idx_0;
                u[1] = u_c_idx_1;
                u[2] = u_c_idx_2;
                cost = cost_c;
                b_iter++;
              }
            }
          }
        }
      }
    }
  }

  if (iter >= 100) {
    *result = 2.0F;
  }
}

/*
 * Wrapper for MATLAB chol for use with auto coder
 *  Inputs:
 * ===========
 *  A - positive semi-definite matrix
 * Arguments    : const float A_data[]
 *                const int A_size[2]
 *                float L_data[]
 *                int L_size[2]
 *                float *fail
 * Return Type  : void
 */
static void chol_free(const float A_data[], const int A_size[2], float L_data[],
                      int L_size[2], float *fail)
{
  int A_size_idx_0;
  int nmj;
  int k;
  float b_A_data[9];
  int n;
  int info;
  int b_n;
  int j;
  bool exitg1;
  int jj;
  float ajj;
  int ix;
  int iy;
  int i9;
  int i10;
  float c_A_data[9];
  int jp1j;
  int iac;
  float c;
  int ia;
  A_size_idx_0 = A_size[0];
  nmj = A_size[1];
  k = A_size[0] * A_size[1];
  if (0 <= k - 1) {
    memcpy(&b_A_data[0], &A_data[0], (unsigned int)(k * (int)sizeof(float)));
  }

  n = A_size[1];
  info = 0;
  if (A_size[1] != 0) {
    b_n = A_size[0];
    info = 0;
    if (A_size[0] != 0) {
      j = 0;
      exitg1 = false;
      while ((!exitg1) && (j <= b_n - 1)) {
        jj = j + j * n;
        ajj = 0.0F;
        if (j >= 1) {
          ix = j;
          iy = j;
          for (k = 0; k < j; k++) {
            ajj += b_A_data[ix] * b_A_data[iy];
            ix += n;
            iy += n;
          }
        }

        ajj = b_A_data[jj] - ajj;
        if (ajj > 0.0F) {
          ajj = sqrtf(ajj);
          b_A_data[jj] = ajj;
          if (j + 1 < b_n) {
            nmj = (b_n - j) - 1;
            k = j + 2;
            jp1j = jj + 2;
            if ((nmj == 0) || (j == 0)) {
            } else {
              ix = j;
              i9 = (j + n * (j - 1)) + 2;
              for (iac = k; n < 0 ? iac >= i9 : iac <= i9; iac += n) {
                c = -b_A_data[ix];
                iy = jj + 1;
                i10 = (iac + nmj) - 1;
                for (ia = iac; ia <= i10; ia++) {
                  b_A_data[iy] += b_A_data[ia - 1] * c;
                  iy++;
                }

                ix += n;
              }
            }

            ajj = 1.0F / ajj;
            i9 = jj + nmj;
            for (k = jp1j; k <= i9 + 1; k++) {
              b_A_data[k - 1] *= ajj;
            }
          }

          j++;
        } else {
          b_A_data[jj] = ajj;
          info = j + 1;
          exitg1 = true;
        }
      }
    }

    if (info == 0) {
      nmj = A_size[1];
    } else {
      nmj = info - 1;
    }

    for (j = 2; j <= nmj; j++) {
      for (k = 0; k <= j - 2; k++) {
        b_A_data[k + A_size_idx_0 * (j - 1)] = 0.0F;
      }
    }

    if (1 > nmj) {
      k = 0;
      nmj = 0;
    } else {
      k = nmj;
    }

    for (i9 = 0; i9 < nmj; i9++) {
      for (i10 = 0; i10 < k; i10++) {
        c_A_data[i10 + k * i9] = b_A_data[i10 + A_size_idx_0 * i9];
      }
    }

    A_size_idx_0 = k;
    k *= nmj;
    if (0 <= k - 1) {
      memcpy(&b_A_data[0], &c_A_data[0], (unsigned int)(k * (int)sizeof(float)));
    }
  }

  *fail = (float)info;
  L_size[0] = A_size_idx_0;
  L_size[1] = nmj;
  k = A_size_idx_0 * nmj;
  if (0 <= k - 1) {
    memcpy(&L_data[0], &b_A_data[0], (unsigned int)(k * (int)sizeof(float)));
  }
}

/*
 * Uses an rk method to roll out a trajectory
 *  Returns the new trajectory, cost and the derivatives along the trajectory
 * Arguments    : const float x[7000]
 *                const float xg[7]
 *                const float u[2997]
 *                const float K[17982]
 *                const float u_lims[6]
 *                const float B_ECI[3000]
 *                float xnew[7000]
 *                float unew[2997]
 *                float fx[35964]
 *                float fu[17982]
 *                float cx[6000]
 *                float cu[2997]
 *                float cxx[36000]
 *                float cuu[8991]
 *                float *cost
 * Return Type  : void
 */
static void forwardRollout(const float x[7000], const float xg[7], const float
  u[2997], const float K[17982], const float u_lims[6], const float B_ECI[3000],
  float xnew[7000], float unew[2997], float fx[35964], float fu[17982], float
  cx[6000], float cu[2997], float cxx[36000], float cuu[8991], float *cost)
{
  int i0;
  int k;
  int dx_tmp;
  float dx[6];
  int b_dx_tmp;
  int c_dx_tmp;
  float q_idx_0;
  int q_idx_1_tmp;
  int q_idx_2_tmp;
  int q_idx_3_tmp;
  float b_xnew;
  float b_xg;
  float out;
  int c_sign;
  float z1[3];
  float q[16];
  int xnew_tmp;
  int b_k;
  float q_error[4];
  float fv5[9];
  int unew_tmp;
  float d_sign[12];
  float b[7];
  float dxdot1[70];
  float c_xnew[7];
  float dxdot2[70];
  float x1[7];
  float fx_tmp[49];
  float b_fx_tmp[42];
  float fv6[49];
  float d_xnew[42];
  float c_fx_tmp[42];
  int b_xnew_tmp;
  float fv7[21];

  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0], 0, 7000U * sizeof(float));
  memset(&unew[0], 0, 2997U * sizeof(float));
  *cost = 0.0F;
  for (i0 = 0; i0 < 7; i0++) {
    xnew[i0] = x[i0];
  }

  for (k = 0; k < 999; k++) {
    /*  Find the state error vector dx */
    dx_tmp = 7 * k + 4;
    dx[3] = xnew[dx_tmp] - x[dx_tmp];
    b_dx_tmp = 7 * k + 5;
    dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
    c_dx_tmp = 7 * k + 6;
    dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
    q_idx_0 = x[7 * k];
    q_idx_1_tmp = 1 + 7 * k;
    q_idx_2_tmp = 2 + 7 * k;
    q_idx_3_tmp = 3 + 7 * k;

    /*  conjugate */
    /*  Forms the "left matrix" for quaternion multiplication */
    /*  where q*q1 = L(q)q1 */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Input */
    /*  q - the quaternion to build the left matrix from */
    /*  Output */
    /*  L - the left matrix */
    /*  L = [s, -v'; */
    /*       v, sI+skew(v)] */
    /* --------------------------------------------------- */
    q[0] = q_idx_0;
    q[4] = x[q_idx_1_tmp];
    q[8] = x[q_idx_2_tmp];
    q[12] = x[q_idx_3_tmp];
    q[1] = -x[q_idx_1_tmp];
    q[5] = q_idx_0;
    q[9] = x[q_idx_3_tmp];
    q[13] = -x[q_idx_2_tmp];
    q[2] = -x[q_idx_2_tmp];
    q[6] = -x[q_idx_3_tmp];
    q[10] = q_idx_0;
    q[14] = x[q_idx_1_tmp];
    q[3] = -x[q_idx_3_tmp];
    q[7] = x[q_idx_2_tmp];
    q[11] = -x[q_idx_1_tmp];
    q[15] = q_idx_0;
    q_idx_0 = 0.0F;
    for (i0 = 0; i0 < 4; i0++) {
      q_error[i0] = 0.0F;
      b_xg = ((q[i0] * xnew[7 * k] + q[i0 + 4] * xnew[q_idx_1_tmp]) + q[i0 + 8] *
              xnew[q_idx_2_tmp]) + q[i0 + 12] * xnew[q_idx_3_tmp];
      q_error[i0] = b_xg;
      q_idx_0 += b_xg * b_xg;
    }

    q_idx_0 = sqrtf(q_idx_0);
    q_error[0] /= q_idx_0;
    q_error[1] /= q_idx_0;
    q_error[2] /= q_idx_0;
    q_error[3] /= q_idx_0;

    /*  re-normalize */
    /*  inverse Cayley Map */
    dx[0] = q_error[1] / q_error[0];
    dx[1] = q_error[2] / q_error[0];
    dx[2] = q_error[3] / q_error[0];

    /*  Find the new control and ensure it is within the limits */
    for (b_k = 0; b_k < 3; b_k++) {
      b_xg = 0.0F;
      for (i0 = 0; i0 < 6; i0++) {
        b_xg += K[(b_k + 3 * i0) + 18 * k] * dx[i0];
      }

      unew_tmp = b_k + 3 * k;
      unew[unew_tmp] = fminf(u_lims[3 + b_k], fmaxf(u_lims[b_k], u[unew_tmp] -
        b_xg));
    }

    /*  Step the dynamics forward */
    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Step dynamics */
    /*  Explicit midpoint step from x_k to x_{k+1} */
    satellite_dynamics(*(float (*)[7])&xnew[7 * k], *(float (*)[3])&unew[3 * k],
                       *(float (*)[3])&B_ECI[3 * k], b, dxdot1);
    for (i0 = 0; i0 < 7; i0++) {
      c_xnew[i0] = xnew[i0 + 7 * k] + 0.015F * b[i0];
    }

    satellite_dynamics(c_xnew, *(float (*)[3])&unew[3 * k], *(float (*)[3])&
                       B_ECI[3 * k], b, dxdot2);
    for (i0 = 0; i0 < 7; i0++) {
      x1[i0] = xnew[i0 + 7 * k] + 0.03F * b[i0];
    }

    /*  Re-normalize the quaternion */
    q_idx_0 = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] *
                    x1[3]);

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    for (i0 = 0; i0 < 7; i0++) {
      for (b_k = 0; b_k < 7; b_k++) {
        unew_tmp = b_k + 7 * i0;
        fx_tmp[unew_tmp] = 0.00045F * dxdot2[unew_tmp];
      }
    }

    fv5[0] = 0.0F;
    fv5[3] = -x1[3];
    fv5[6] = x1[2];
    fv5[1] = x1[3];
    fv5[4] = 0.0F;
    fv5[7] = -x1[1];
    fv5[2] = -x1[2];
    fv5[5] = x1[1];
    fv5[8] = 0.0F;
    for (i0 = 0; i0 < 3; i0++) {
      b_fx_tmp[i0] = -x1[1 + i0];
      b_fx_tmp[i0 + 3] = 0.0F;
      unew_tmp = 6 * (i0 + 1);
      b_fx_tmp[unew_tmp] = x1[0] * (float)iv0[i0] + fv5[i0];
      b_fx_tmp[unew_tmp + 3] = 0.0F;
      b_fx_tmp[1 + unew_tmp] = x1[0] * (float)iv0[i0 + 3] + fv5[i0 + 3];
      b_fx_tmp[unew_tmp + 4] = 0.0F;
      b_fx_tmp[2 + unew_tmp] = x1[0] * (float)iv0[i0 + 6] + fv5[i0 + 6];
      b_fx_tmp[unew_tmp + 5] = 0.0F;
      for (b_k = 0; b_k < 6; b_k++) {
        b_fx_tmp[b_k + 6 * (i0 + 4)] = iv1[i0 + 3 * b_k];
      }
    }

    for (i0 = 0; i0 < 7; i0++) {
      for (b_k = 0; b_k < 7; b_k++) {
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += fx_tmp[i0 + 7 * xnew_tmp] * dxdot1[xnew_tmp + 7 * b_k];
        }

        xnew_tmp = i0 + 7 * b_k;
        fv6[xnew_tmp] = ((float)iv2[xnew_tmp] + 0.03F * dxdot2[xnew_tmp]) + b_xg;
      }
    }

    b_xnew = xnew[7 * k];
    fv5[0] = 0.0F;
    fv5[3] = -xnew[3 + 7 * k];
    fv5[6] = xnew[2 + 7 * k];
    fv5[1] = xnew[3 + 7 * k];
    fv5[4] = 0.0F;
    fv5[7] = -xnew[1 + 7 * k];
    fv5[2] = -xnew[2 + 7 * k];
    fv5[5] = xnew[1 + 7 * k];
    fv5[8] = 0.0F;
    for (i0 = 0; i0 < 6; i0++) {
      for (b_k = 0; b_k < 7; b_k++) {
        unew_tmp = i0 + 6 * b_k;
        c_fx_tmp[unew_tmp] = 0.0F;
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += b_fx_tmp[i0 + 6 * xnew_tmp] * fv6[xnew_tmp + 7 * b_k];
        }

        c_fx_tmp[unew_tmp] = b_xg;
      }
    }

    for (i0 = 0; i0 < 3; i0++) {
      d_xnew[7 * i0] = -xnew[(i0 + 7 * k) + 1];
      xnew_tmp = 7 * (i0 + 3);
      d_xnew[xnew_tmp] = 0.0F;
      d_xnew[7 * i0 + 1] = b_xnew * (float)iv0[3 * i0] + fv5[3 * i0];
      d_xnew[xnew_tmp + 1] = 0.0F;
      b_xnew_tmp = 1 + 3 * i0;
      d_xnew[7 * i0 + 2] = b_xnew * (float)iv0[b_xnew_tmp] + fv5[b_xnew_tmp];
      d_xnew[xnew_tmp + 2] = 0.0F;
      b_xnew_tmp = 2 + 3 * i0;
      d_xnew[7 * i0 + 3] = b_xnew * (float)iv0[b_xnew_tmp] + fv5[b_xnew_tmp];
      d_xnew[xnew_tmp + 3] = 0.0F;
    }

    for (i0 = 0; i0 < 6; i0++) {
      d_xnew[7 * i0 + 4] = iv1[3 * i0];
      d_xnew[7 * i0 + 5] = iv1[1 + 3 * i0];
      d_xnew[7 * i0 + 6] = iv1[2 + 3 * i0];
    }

    for (i0 = 0; i0 < 6; i0++) {
      for (b_k = 0; b_k < 6; b_k++) {
        unew_tmp = (i0 + 6 * b_k) + 36 * k;
        fx[unew_tmp] = 0.0F;
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += c_fx_tmp[i0 + 6 * xnew_tmp] * d_xnew[xnew_tmp + 7 * b_k];
        }

        fx[unew_tmp] = b_xg;
      }
    }

    for (i0 = 0; i0 < 7; i0++) {
      for (b_k = 0; b_k < 3; b_k++) {
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += fx_tmp[i0 + 7 * xnew_tmp] * dxdot1[xnew_tmp + 7 * (7 + b_k)];
        }

        fv7[i0 + 7 * b_k] = 0.03F * dxdot2[i0 + 7 * (7 + b_k)] + b_xg;
      }
    }

    for (i0 = 0; i0 < 6; i0++) {
      for (b_k = 0; b_k < 3; b_k++) {
        unew_tmp = (i0 + 6 * b_k) + 18 * k;
        fu[unew_tmp] = 0.0F;
        b_xg = 0.0F;
        for (xnew_tmp = 0; xnew_tmp < 7; xnew_tmp++) {
          b_xg += b_fx_tmp[i0 + 6 * xnew_tmp] * fv7[xnew_tmp + 7 * b_k];
        }

        fu[unew_tmp] = b_xg;
      }
    }

    xnew_tmp = 7 * (k + 1);
    xnew[xnew_tmp] = x1[0] / q_idx_0;
    xnew[1 + xnew_tmp] = x1[1] / q_idx_0;
    xnew[2 + xnew_tmp] = x1[2] / q_idx_0;
    xnew[3 + xnew_tmp] = x1[3] / q_idx_0;
    xnew[xnew_tmp + 4] = x1[4];
    xnew[xnew_tmp + 5] = x1[5];
    xnew[xnew_tmp + 6] = x1[6];

    /*  Calculate the cost */
    /*  Calculates the cost contribution of a given state and control */
    /*  Utilizes a standard quadratic cost function for the angular velocity */
    /*  and a geodesic cost for the attitude quaternion */
    /*  Also returns the 2nd order expansion of the cost function */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Inputs */
    /* ===================================== */
    /*  x        - [quaternion; omega] (7x1) */
    /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
    /*  terminal - integer 0 or 1 */
    /* --------------------------------------------------- */
    /*  control hessian */
    /*  cumulative angular velocity cost hessian */
    /*  cumulative geodesic cost weighting */
    /*  Finds the geodesic quaternion-error cost */
    /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
    /*  Also records the sign (+ or -) which minimizes quat_cost */
    /*  this is used when calculating the Jacobain and Hessian */
    q_idx_0 = xnew[7 * k];
    b_xnew = xg[0] * q_idx_0;
    b_xg = ((b_xnew + xg[1] * xnew[q_idx_1_tmp]) + xg[2] * xnew[q_idx_2_tmp]) +
      xg[3] * xnew[q_idx_3_tmp];
    if (1.0F + b_xg < 1.0F - b_xg) {
      out = 1.0F + b_xg;
      c_sign = 1;
    } else {
      out = 1.0F - b_xg;
      c_sign = -1;
    }

    z1[0] = xnew[dx_tmp] - xg[4];
    z1[1] = xnew[b_dx_tmp] - xg[5];
    z1[2] = xnew[c_dx_tmp] - xg[6];

    /*  State cost Hessian */
    b_xg = ((b_xnew + xg[1] * xnew[1 + 7 * k]) + xg[2] * xnew[2 + 7 * k]) + xg[3]
      * xnew[3 + 7 * k];
    for (i0 = 0; i0 < 3; i0++) {
      b_k = 6 * i0 + 36 * k;
      cxx[b_k] = (float)(-c_sign * iv0[3 * i0]) * b_xg;
      unew_tmp = 6 * (i0 + 3) + 36 * k;
      cxx[unew_tmp] = 0.0F;
      cxx[b_k + 1] = (float)(-c_sign * iv0[1 + 3 * i0]) * b_xg;
      cxx[unew_tmp + 1] = 0.0F;
      cxx[b_k + 2] = (float)(-c_sign * iv0[2 + 3 * i0]) * b_xg;
      cxx[unew_tmp + 2] = 0.0F;
    }

    for (i0 = 0; i0 < 6; i0++) {
      b_k = 6 * i0 + 36 * k;
      cxx[b_k + 3] = fv0[3 * i0];
      cxx[b_k + 4] = fv0[1 + 3 * i0];
      cxx[b_k + 5] = fv0[2 + 3 * i0];
    }

    /*  State cost Jacobian */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    fv5[0] = 0.0F;
    fv5[3] = -xnew[3 + 7 * k];
    fv5[6] = xnew[2 + 7 * k];
    fv5[1] = xnew[3 + 7 * k];
    fv5[4] = 0.0F;
    fv5[7] = -xnew[1 + 7 * k];
    fv5[2] = -xnew[2 + 7 * k];
    fv5[5] = xnew[1 + 7 * k];
    fv5[8] = 0.0F;
    for (i0 = 0; i0 < 3; i0++) {
      d_sign[i0] = (float)c_sign * -xnew[(i0 + 7 * k) + 1];
      unew_tmp = 3 * (i0 + 1);
      d_sign[unew_tmp] = (float)c_sign * (q_idx_0 * (float)iv0[i0] + fv5[i0]);
      d_sign[1 + unew_tmp] = (float)c_sign * (q_idx_0 * (float)iv0[i0 + 3] +
        fv5[i0 + 3]);
      d_sign[2 + unew_tmp] = (float)c_sign * (q_idx_0 * (float)iv0[i0 + 6] +
        fv5[i0 + 6]);
    }

    /*  Control cost Hessian & Jacobian */
    b_xg = 0.0F;
    for (i0 = 0; i0 < 3; i0++) {
      unew_tmp = i0 + 6 * k;
      cx[unew_tmp] = ((d_sign[i0] * xg[0] + d_sign[i0 + 3] * xg[1]) + d_sign[i0
                      + 6] * xg[2]) + d_sign[i0 + 9] * xg[3];
      cx[unew_tmp + 3] = (fv1[i0] * z1[0] + fv1[i0 + 3] * z1[1]) + fv1[i0 + 6] *
        z1[2];
      unew_tmp = i0 + 3 * k;
      cu[unew_tmp] = 0.0F;
      xnew_tmp = 3 * i0 + 9 * k;
      cuu[xnew_tmp] = fv2[3 * i0];
      b_k = 1 + 3 * i0;
      cuu[xnew_tmp + 1] = fv2[b_k];
      b_xnew_tmp = 2 + 3 * i0;
      cuu[xnew_tmp + 2] = fv2[b_xnew_tmp];
      cu[unew_tmp] = (fv2[i0] * unew[3 * k] + fv2[i0 + 3] * unew[1 + 3 * k]) +
        fv2[i0 + 6] * unew[2 + 3 * k];
      b_xg += ((0.5F * z1[0] * fv1[3 * i0] + 0.5F * z1[1] * fv1[b_k]) + 0.5F *
               z1[2] * fv1[b_xnew_tmp]) * z1[i0];
    }

    q_idx_0 = 0.0F;
    for (i0 = 0; i0 < 3; i0++) {
      q_idx_0 += ((0.5F * unew[3 * k] * fv2[3 * i0] + 0.5F * unew[1 + 3 * k] *
                   fv2[1 + 3 * i0]) + 0.5F * unew[2 + 3 * k] * fv2[2 + 3 * i0]) *
        unew[i0 + 3 * k];
    }

    *cost += (out + b_xg) + q_idx_0;
  }

  /*  Final cost */
  /*  Calculates the cost contribution of a given state and control */
  /*  Utilizes a standard quadratic cost function for the angular velocity */
  /*  and a geodesic cost for the attitude quaternion */
  /*  Also returns the 2nd order expansion of the cost function */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Inputs */
  /* ===================================== */
  /*  x        - [quaternion; omega] (7x1) */
  /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
  /*  terminal - integer 0 or 1 */
  /* --------------------------------------------------- */
  /*  control hessian */
  /*  terminal angular velocity cost hessian */
  /*  terminal geodesic cost weight */
  /*  Finds the geodesic quaternion-error cost */
  /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
  /*  Also records the sign (+ or -) which minimizes quat_cost */
  /*  this is used when calculating the Jacobain and Hessian */
  b_xnew = xg[0] * xnew[6993];
  b_xg = ((b_xnew + xg[1] * xnew[6994]) + xg[2] * xnew[6995]) + xg[3] * xnew
    [6996];
  if (1.0F + b_xg < 1.0F - b_xg) {
    out = 1.0F + b_xg;
    c_sign = 1;
  } else {
    out = 1.0F - b_xg;
    c_sign = -1;
  }

  z1[0] = xnew[6997] - xg[4];
  z1[1] = xnew[6998] - xg[5];
  z1[2] = xnew[6999] - xg[6];

  /*  State cost Hessian */
  xnew_tmp = -c_sign * 10;
  b_xg = ((b_xnew + xg[1] * xnew[6994]) + xg[2] * xnew[6995]) + xg[3] * xnew
    [6996];
  for (i0 = 0; i0 < 3; i0++) {
    cxx[35964 + 6 * i0] = (float)(xnew_tmp * iv0[3 * i0]) * b_xg;
    b_k = 6 * (i0 + 3);
    cxx[35964 + b_k] = 0.0F;
    cxx[6 * i0 + 35965] = (float)(xnew_tmp * iv0[1 + 3 * i0]) * b_xg;
    cxx[b_k + 35965] = 0.0F;
    cxx[6 * i0 + 35966] = (float)(xnew_tmp * iv0[2 + 3 * i0]) * b_xg;
    cxx[b_k + 35966] = 0.0F;
  }

  for (i0 = 0; i0 < 6; i0++) {
    cxx[6 * i0 + 35967] = iv1[3 * i0];
    cxx[6 * i0 + 35968] = iv1[1 + 3 * i0];
    cxx[6 * i0 + 35969] = iv1[2 + 3 * i0];
  }

  /*  State cost Jacobian */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  xnew_tmp = c_sign * 10;
  fv5[0] = 0.0F;
  fv5[3] = -xnew[6996];
  fv5[6] = xnew[6995];
  fv5[1] = xnew[6996];
  fv5[4] = 0.0F;
  fv5[7] = -xnew[6994];
  fv5[2] = -xnew[6995];
  fv5[5] = xnew[6994];
  fv5[8] = 0.0F;
  for (i0 = 0; i0 < 3; i0++) {
    d_sign[i0] = (float)xnew_tmp * -xnew[i0 + 6994];
    unew_tmp = 3 * (i0 + 1);
    d_sign[unew_tmp] = (float)xnew_tmp * (xnew[6993] * (float)iv0[i0] + fv5[i0]);
    d_sign[1 + unew_tmp] = (float)xnew_tmp * (xnew[6993] * (float)iv0[i0 + 3] +
      fv5[i0 + 3]);
    d_sign[2 + unew_tmp] = (float)xnew_tmp * (xnew[6993] * (float)iv0[i0 + 6] +
      fv5[i0 + 6]);
  }

  /*  Control cost Hessian & Jacobian */
  b_xg = 0.0F;
  for (i0 = 0; i0 < 3; i0++) {
    cx[5994 + i0] = ((d_sign[i0] * xg[0] + d_sign[i0 + 3] * xg[1]) + d_sign[i0 +
                     6] * xg[2]) + d_sign[i0 + 9] * xg[3];
    cx[i0 + 5997] = ((float)iv0[i0] * z1[0] + (float)iv0[i0 + 3] * z1[1]) +
      (float)iv0[i0 + 6] * z1[2];
    b_xg += ((0.5F * z1[0] * (float)iv0[3 * i0] + 0.5F * z1[1] * (float)iv0[1 +
              3 * i0]) + 0.5F * z1[2] * (float)iv0[2 + 3 * i0]) * z1[i0];
  }

  *cost += 10.0F * out + b_xg;
}

/*
 * Arguments    : const float x[999]
 * Return Type  : float
 */
static float mean(const float x[999])
{
  float y;
  int k;
  y = x[0];
  for (k = 0; k < 998; k++) {
    y += x[k + 1];
  }

  y /= 999.0F;
  return y;
}

/*
 * Arguments    : float y[11]
 * Return Type  : void
 */
static void power(float y[11])
{
  int k;
  for (k = 0; k < 11; k++) {
    y[k] = powf(10.0F, -0.3F * (float)k);
  }
}

/*
 * Calculates the continuous time state derivative and Jacobians
 * Arguments    : const float x[7]
 *                const float u[3]
 *                const float B_ECI[3]
 *                float xdot[7]
 *                float dxdot[70]
 * Return Type  : void
 */
static void satellite_dynamics(const float x[7], const float u[3], const float
  B_ECI[3], float xdot[7], float dxdot[70])
{
  float c;
  float q_out[4];
  float B_B_idx_0;
  float B_B_idx_1;
  float B_B_idx_2;
  float qdot_tmp[9];
  int i1;
  float wdot_tmp[9];
  float B_B[3];
  int i2;
  float fv8[12];
  float q1[3];
  int a_tmp;
  static const signed char a[9] = { 100, 0, 0, 0, 100, 0, 0, 0, 100 };

  int wdot_tmp_tmp;
  float b_wdot_tmp[9];
  float fv9[9];
  float b_a[9];
  static const short c_a[9] = { -200, 0, 0, 0, -200, 0, 0, 0, -200 };

  int dxdot_tmp;
  int b_dxdot_tmp;
  static const signed char d_a[9] = { -100, 0, 0, 0, -100, 0, 0, 0, -100 };

  /*  Inputs */
  /* ======================= */
  /*  x - the current state */
  /*  u - the current control (magentic moment vector) */
  /*  B - Earth magnetic field vector in ECI co-ordinates (3x1) */
  /*  TODO: ideally, we'd like to pass in a B_B (body frame magnetic field */
  /*  vector) */
  /*  kgm^2 */
  /*  Angular velocity */
  /*  Quaternion */
  /*  Magnetic field section  */
  /*  given B_ECI (ECI magnetic field at the time step) */
  /*  this function multiplie s a vector by a quaternion rotation */
  /*  this function multiplies quaternions */
  c = (B_ECI[0] * x[1] + B_ECI[1] * x[2]) + B_ECI[2] * x[3];
  q_out[1] = x[0] * B_ECI[0] + (B_ECI[1] * x[3] - B_ECI[2] * x[2]);
  q_out[2] = x[0] * B_ECI[1] + (B_ECI[2] * x[1] - B_ECI[0] * x[3]);
  q_out[3] = x[0] * B_ECI[2] + (B_ECI[0] * x[2] - B_ECI[1] * x[1]);

  /*  this function multiplies quaternions */
  B_B_idx_0 = (x[0] * q_out[1] + (0.0F - c) * -x[1]) + (-x[2] * q_out[3] - -x[3]
    * q_out[2]);
  B_B_idx_1 = (x[0] * q_out[2] + (0.0F - c) * -x[2]) + (-x[3] * q_out[1] - -x[1]
    * q_out[3]);
  B_B_idx_2 = (x[0] * q_out[3] + (0.0F - c) * -x[3]) + (-x[1] * q_out[2] - -x[2]
    * q_out[1]);

  /*  Non-linear dynamics */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  qdot_tmp[0] = 0.0F;
  qdot_tmp[3] = -x[3];
  qdot_tmp[6] = x[2];
  qdot_tmp[1] = x[3];
  qdot_tmp[4] = 0.0F;
  qdot_tmp[7] = -x[1];
  qdot_tmp[2] = -x[2];
  qdot_tmp[5] = x[1];
  qdot_tmp[8] = 0.0F;
  for (i1 = 0; i1 < 9; i1++) {
    qdot_tmp[i1] += x[0] * (float)iv0[i1];
  }

  /*  wdot = Jinv*(u - skew_mat(w)*J*w); */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  wdot_tmp[0] = 0.0F;
  wdot_tmp[3] = -x[6];
  wdot_tmp[6] = x[5];
  wdot_tmp[1] = x[6];
  wdot_tmp[4] = 0.0F;
  wdot_tmp[7] = -x[4];
  wdot_tmp[2] = -x[5];
  wdot_tmp[5] = x[4];
  wdot_tmp[8] = 0.0F;

  /*  with magnetorquers */
  B_B[0] = -(B_B_idx_1 * u[2] - B_B_idx_2 * u[1]);
  B_B[1] = -(B_B_idx_2 * u[0] - B_B_idx_0 * u[2]);
  B_B[2] = -(B_B_idx_0 * u[1] - B_B_idx_1 * u[0]);
  for (i1 = 0; i1 < 3; i1++) {
    i2 = i1 << 2;
    fv8[i2] = 0.5F * -x[1 + i1];
    q1[i1] = 0.0F;
    for (a_tmp = 0; a_tmp < 3; a_tmp++) {
      wdot_tmp_tmp = i1 + 3 * a_tmp;
      b_wdot_tmp[wdot_tmp_tmp] = 0.0F;
      c = (wdot_tmp[i1] * fv1[3 * a_tmp] + wdot_tmp[i1 + 3] * fv1[1 + 3 * a_tmp])
        + wdot_tmp[i1 + 6] * fv1[2 + 3 * a_tmp];
      b_wdot_tmp[wdot_tmp_tmp] = c;
      fv8[(a_tmp + i2) + 1] = 0.5F * qdot_tmp[a_tmp + 3 * i1];
      q1[i1] += c * x[4 + a_tmp];
    }

    B_B[i1] -= q1[i1];
  }

  for (i1 = 0; i1 < 4; i1++) {
    q_out[i1] = 0.0F;
    q_out[i1] = (fv8[i1] * x[4] + fv8[i1 + 4] * x[5]) + fv8[i1 + 8] * x[6];
  }

  for (i1 = 0; i1 < 3; i1++) {
    q1[i1] = 0.0F;
    q1[i1] = ((float)a[i1] * B_B[0] + (float)a[i1 + 3] * B_B[1]) + (float)a[i1 +
      6] * B_B[2];
  }

  xdot[0] = q_out[0];
  xdot[1] = q_out[1];
  xdot[2] = q_out[2];
  xdot[3] = q_out[3];

  /*  Jacobians  */
  for (i1 = 0; i1 < 3; i1++) {
    xdot[i1 + 4] = q1[i1];
    q1[i1] = (fv1[i1] * x[4] + fv1[i1 + 3] * x[5]) + fv1[i1 + 6] * x[6];
  }

  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  B = [zeros(4,3); */
  /*       Jinv]; */
  /*  with magnetorquers: */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  fv9[0] = 0.0F;
  fv9[3] = -q1[2];
  fv9[6] = q1[1];
  fv9[1] = q1[2];
  fv9[4] = 0.0F;
  fv9[7] = -q1[0];
  fv9[2] = -q1[1];
  fv9[5] = q1[0];
  fv9[8] = 0.0F;
  for (i1 = 0; i1 < 9; i1++) {
    b_wdot_tmp[i1] -= fv9[i1];
  }

  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      a_tmp = i1 + 3 * i2;
      b_a[a_tmp] = 0.0F;
      b_a[a_tmp] = ((float)c_a[i1] * b_wdot_tmp[3 * i2] + (float)c_a[i1 + 3] *
                    b_wdot_tmp[1 + 3 * i2]) + (float)c_a[i1 + 6] * b_wdot_tmp[2
        + 3 * i2];
    }
  }

  fv9[0] = 0.0F;
  fv9[3] = -B_B_idx_2;
  fv9[6] = B_B_idx_1;
  fv9[1] = B_B_idx_2;
  fv9[4] = 0.0F;
  fv9[7] = -B_B_idx_0;
  fv9[2] = -B_B_idx_1;
  fv9[5] = B_B_idx_0;
  fv9[8] = 0.0F;
  dxdot[0] = 0.0F;
  for (i1 = 0; i1 < 3; i1++) {
    c = x[4 + i1];
    dxdot_tmp = 7 * (i1 + 1);
    dxdot[dxdot_tmp] = 0.5F * -c;
    b_dxdot_tmp = 7 * (i1 + 4);
    dxdot[b_dxdot_tmp] = 0.5F * -x[1 + i1];
    dxdot[i1 + 1] = 0.5F * c;
    for (i2 = 0; i2 < 3; i2++) {
      wdot_tmp_tmp = i1 + 3 * i2;
      b_wdot_tmp[wdot_tmp_tmp] = 0.0F;
      b_wdot_tmp[wdot_tmp_tmp] = ((float)d_a[i1] * fv9[3 * i2] + (float)d_a[i1 +
        3] * fv9[1 + 3 * i2]) + (float)d_a[i1 + 6] * fv9[2 + 3 * i2];
      a_tmp = i2 + 3 * i1;
      dxdot[(i2 + dxdot_tmp) + 1] = 0.5F * -wdot_tmp[a_tmp];
      dxdot[(i2 + b_dxdot_tmp) + 1] = 0.5F * qdot_tmp[a_tmp];
    }
  }

  for (i1 = 0; i1 < 4; i1++) {
    dxdot[7 * i1 + 4] = 0.0F;
    dxdot[7 * i1 + 5] = 0.0F;
    dxdot[7 * i1 + 6] = 0.0F;
  }

  for (i1 = 0; i1 < 3; i1++) {
    dxdot_tmp = 7 * (i1 + 4);
    dxdot[dxdot_tmp + 4] = 0.5F * b_a[3 * i1];
    b_dxdot_tmp = 1 + 3 * i1;
    dxdot[dxdot_tmp + 5] = 0.5F * b_a[b_dxdot_tmp];
    a_tmp = 2 + 3 * i1;
    dxdot[dxdot_tmp + 6] = 0.5F * b_a[a_tmp];
    dxdot_tmp = 7 * (i1 + 7);
    dxdot[dxdot_tmp] = 0.0F;
    dxdot[1 + dxdot_tmp] = 0.0F;
    dxdot[2 + dxdot_tmp] = 0.0F;
    dxdot[3 + dxdot_tmp] = 0.0F;
    dxdot[dxdot_tmp + 4] = b_wdot_tmp[3 * i1];
    dxdot[dxdot_tmp + 5] = b_wdot_tmp[b_dxdot_tmp];
    dxdot[dxdot_tmp + 6] = b_wdot_tmp[a_tmp];
  }
}

/*
 * Increases or decreases the regularization parameter according
 *  to a non-linear scaling regime.
 * Arguments    : float *lambda
 *                float direction
 * Return Type  : void
 */
static void updateLambda(float *lambda, float direction)
{
  if (direction == 1.0F) {
    /*  increase lambda */
    *lambda = fmaxf(*lambda * 1.6F, 1.0E-6F);
  } else {
    /*  decrease lambda */
    *lambda = *lambda * 0.625F * (float)(*lambda > 1.0E-6F);

    /*  set = 0 if lambda too small */
  }
}

/*
 * Solves finite horizon optimal control problem using a
 *  multiplicative iterative linear quadratic regulator
 * Arguments    : const float x0[7000]
 *                const float xg[7]
 *                const float u0[2997]
 *                const float u_lims[6]
 *                const float B_ECI[3000]
 *                float x[7000]
 *                float u[2997]
 *                float K[17982]
 *                bool *result
 * Return Type  : void
 */
void milqr(const float x0[7000], const float xg[7], const float u0[2997], const
           float u_lims[6], const float B_ECI[3000], float x[7000], float u[2997],
           float K[17982], bool *result)
{
  float alphas[11];
  float lambda;
  static float l[2997];
  float dV[2];
  static float fv3[17982];
  static float fx[35964];
  static float fu[17982];
  static float cx[6000];
  static float cu[2997];
  static float cxx[36000];
  static float cuu[8991];
  float cost;
  float cost_n;
  static float x_n[7000];
  static float u_n[2997];
  static float fx_n[35964];
  static float fu_n[17982];
  static float cx_n[6000];
  static float cu_n[2997];
  static float cxx_n[36000];
  static float cuu_n[8991];
  int iter;
  int b_iter;
  bool exitg1;
  bool backPassCheck;
  int exitg2;
  bool diverged;
  static float z[2997];
  static float fv4[2997];
  int j;
  float maxval[999];
  float f0;
  bool lineSearchCheck;
  bool exitg3;
  float expected_change;
  float dcost;
  float c_ratio;

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Note that this solver will not work without control limits */
  /*  Inputs */
  /*  =========================================== */
  /*  x0 - The intial trajectory (n, N) */
  /*  */
  /*  xg - The goal state (n, 1) */
  /*  */
  /*  u0 - The initial control sequeunce (m, N-1) */
  /*  */
  /*  u_lims - The control limits (m, 2) (lower, upper) */
  /*  */
  /*  B_ECI - A time sequence of Earth magnetic field vectors in ECI (3, N) */
  /*  Outputs */
  /*  =========================================== */
  /*  x - Final nominal trajectory (n, N) */
  /*  */
  /*  u - Final open-loop controls (m, N-1) */
  /*  */
  /*  K - Feedback control gains (n-1, m, N-1)  */
  /*  */
  /*  result - Indicates convergence (boolean) */
  /*  Options (pass in as array) */
  /*  Timestep (Should be highest we can get away with) */
  /*  maximum iterations */
  /*  cost reduction exit tolerance */
  /*  feedforward control change exit criterion */
  /*  max regularizaion param allowed for exit */
  /*  minimum accepted cost reduction ratio */
  /*  maximum regularization parameter */
  /*  set lambda = 0 below this value */
  /*  amount to scale dlambda by */
  /*  Init optimisation params */
  power(alphas);

  /*  line search param vector */
  lambda = 1.0F;

  /*  error state size (3 param. error representation for attitude) */
  /*  Initial Forward rollout */
  memset(&l[0], 0, 2997U * sizeof(float));
  memset(&K[0], 0, 17982U * sizeof(float));
  dV[0] = 0.0F;
  dV[1] = 0.0F;
  memset(&fv3[0], 0, 17982U * sizeof(float));
  forwardRollout(x0, xg, u0, fv3, u_lims, B_ECI, x, u, fx, fu, cx, cu, cxx, cuu,
                 &cost);

  /*  Convergence check params */
  /*  Expected cost change */
  /*  Ratio of cost change to expected cost change */
  *result = false;

  /*  variable initialization */
  cost = 0.0F;
  cost_n = 0.0F;
  memcpy(&x_n[0], &x0[0], 7000U * sizeof(float));
  memcpy(&u_n[0], &u0[0], 2997U * sizeof(float));
  memset(&fx_n[0], 0, 35964U * sizeof(float));
  memset(&fu_n[0], 0, 17982U * sizeof(float));
  memset(&cx_n[0], 0, 6000U * sizeof(float));
  memset(&cu_n[0], 0, 2997U * sizeof(float));
  memset(&cxx_n[0], 0, 36000U * sizeof(float));
  memset(&cuu_n[0], 0, 8991U * sizeof(float));

  /*  fprintf("\n=====Running MILQR Optimisation====\n"); */
  iter = 1;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter < 300)) {
    iter = b_iter + 1;

    /*      fprintf("\n---New Iteration---"); */
    /*      fprintf("\n lambda: "); */
    /*      fprintf(string(lambda)); */
    /*      if exist('dcost') */
    /*          fprintf("\n cost change: "); */
    /*          fprintf(string(dcost)); */
    /*      end */
    /*      fprintf("\n"); */
    /*  Backward Pass */
    /* ======================================= */
    backPassCheck = false;
    do {
      exitg2 = 0;
      if (!backPassCheck) {
        backwardPass(fx, fu, cx, cu, cxx, cuu, lambda, u_lims, u, l, K, dV,
                     &diverged);
        if (diverged) {
          /*  fprintf("---Warning: Cholesky factorizaton failed---\n"); */
          /*  Increase regularization parameter (lambda) */
          updateLambda(&lambda, 1.0F);
          if (lambda > 1.0E+7F) {
            exitg2 = 1;
          }
        } else {
          backPassCheck = true;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);

    /*  Check relative change of feedforward control */
    /*  Terminate if sufficiently small (success) */
    b_abs(l, z);
    b_abs(u, fv4);
    for (j = 0; j < 2997; j++) {
      z[j] /= fv4[j] + 1.0F;
    }

    for (j = 0; j < 999; j++) {
      maxval[j] = z[3 * j];
      f0 = z[3 * j + 1];
      if (maxval[j] < f0) {
        maxval[j] = f0;
      }

      f0 = z[3 * j + 2];
      if (maxval[j] < f0) {
        maxval[j] = f0;
      }
    }

    /*  Avg over time of max change */
    if ((lambda < 1.0E-5F) && (mean(maxval) < 0.0001F)) {
      /*  fprintf("\n---Success: Control change decreased below tolerance---\n"); */
      *result = true;
      exitg1 = true;
    } else {
      /*  Forward Line-Search */
      /* =========================================== */
      lineSearchCheck = false;
      if (backPassCheck) {
        j = 0;
        exitg3 = false;
        while ((!exitg3) && (j < 11)) {
          b_forwardRollout(x, xg, u, l, K, alphas[j], u_lims, B_ECI, x_n, u_n,
                           fx_n, fu_n, cx_n, cu_n, cxx_n, cuu_n, &cost_n);
          expected_change = alphas[j] * dV[0] + powf(alphas[j], 2.0F) * dV[1];
          if (expected_change < 0.0F) {
            c_ratio = (cost_n - cost) / expected_change;
          } else {
            /*  Non-positive expected cost reduction */
            /*  actual cost change must be negative to accept the step */
            f0 = cost_n - cost;
            b_sign(&f0);
            c_ratio = -f0;
          }

          if (c_ratio > 0.0F) {
            lineSearchCheck = true;
            exitg3 = true;
          } else {
            j++;
          }
        }
      }

      /*  Parameter Updates */
      /* ============================================= */
      if (lineSearchCheck) {
        /*  Decrease Lambda */
        updateLambda(&lambda, -1.0F);
        dcost = cost - cost_n;

        /*  Update the trajectory and controls */
        memcpy(&x[0], &x_n[0], 7000U * sizeof(float));
        memcpy(&u[0], &u_n[0], 2997U * sizeof(float));
        memcpy(&fx[0], &fx_n[0], 35964U * sizeof(float));
        memcpy(&fu[0], &fu_n[0], 17982U * sizeof(float));
        memcpy(&cx[0], &cx_n[0], 6000U * sizeof(float));
        memcpy(&cu[0], &cu_n[0], 2997U * sizeof(float));
        memcpy(&cxx[0], &cxx_n[0], 36000U * sizeof(float));
        memcpy(&cuu[0], &cuu_n[0], 8991U * sizeof(float));
        cost = cost_n;

        /*  Change in cost small enough to terminate? */
        if (dcost < 1.0E-7F) {
          *result = true;

          /*  fprintf('\n---Success: cost change < tolerance---\n'); */
          exitg1 = true;
        } else {
          b_iter++;
        }
      } else {
        /*  No cost reduction (based on cost change ratio) */
        /*  Increase lambda */
        updateLambda(&lambda, 1.0F);
        if (lambda > 1.0E+7F) {
          /*  fprintf("\n---Diverged: new lambda > lambda_max---\n"); */
          exitg1 = true;
        } else {
          b_iter++;
        }
      }
    }
  }

  if (iter == 300) {
    /*  Ddin't converge completely */
    *result = false;

    /*      % fprintf("\n---Warning: Max iterations exceeded---\n"); */
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void milqr_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void milqr_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for milqr.c
 *
 * [EOF]
 */
