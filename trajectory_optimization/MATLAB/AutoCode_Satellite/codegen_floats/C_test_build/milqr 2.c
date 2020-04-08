/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 19-Mar-2020 19:23:47
 */

/* Include Files */
#include "milqr.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static const signed char iv[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

static const signed char iv1[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
    0, 0, 1 };

static const float fv[18] = { 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
    0.0F, 0.01F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F, 0.01F };

static const float fv1[9] = { 0.01F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F,
    0.01F };

static const float fv2[9] = { 0.05F, 0.0F, 0.0F, 0.0F, 0.05F, 0.0F, 0.0F, 0.0F,
    0.05F };

/* Function Declarations */
static void b_forwardRollout(const float x[3500], const float xg[7], const float
    u[1497], const float l[1497], const float K[8982], float alpha, const float
    u_lims[6], float xnew[3500], float unew[1497], float fx[17964], float fu
    [8982], float cx[3000], float cu[1497], float cxx[18000], float cuu[4491],
    float *cost);
static void backwardPass(const float fx[17964], const float fu[8982], const
    float cx[3000], const float cu[1497], const float cxx[18000], const float
    cuu[4491], float lambda, const float u_lims[6], const float u[1497], float
    l[1497], float K[8982], float dV[2], bool *diverge);
static void boxQPsolve(const float Quu[9], const float Qu[3], const float lower
                       [3], const float upper[3], const float u0[3], float u[3],
                       float *result, float Luu_data[], int Luu_size[2],
                       bool b_free[3]);
static void chol_free(const float A_data[], const int A_size[2], float L_data[],
                      int L_size[2], float *fail);
static void forwardRollout(const float x[3500], const float xg[7], const float
    u[1497], const float K[8982], const float u_lims[6], float xnew[3500], float
    unew[1497], float fx[17964], float fu[8982], float cx[3000], float cu[1497],
    float cxx[18000], float cuu[4491], float *cost);
static void satellite_step(const float x0[7], const float u0[3], float x[7],
    float fx[36], float fu[18]);

/* Function Definitions */

/*
 * Arguments    : const float x[3500]
 *                const float xg[7]
 *                const float u[1497]
 *                const float l[1497]
 *                const float K[8982]
 *                float alpha
 *                const float u_lims[6]
 *                float xnew[3500]
 *                float unew[1497]
 *                float fx[17964]
 *                float fu[8982]
 *                float cx[3000]
 *                float cu[1497]
 *                float cxx[18000]
 *                float cuu[4491]
 *                float *cost
 * Return Type  : void
 */
static void b_forwardRollout(const float x[3500], const float xg[7], const float
    u[1497], const float l[1497], const float K[8982], float alpha, const float
    u_lims[6], float xnew[3500], float unew[1497], float fx[17964], float fu
    [8982], float cx[3000], float cu[1497], float cxx[18000], float cuu[4491],
    float *cost)
{
    int i;
    float b_fv[9];
    int k;
    int dx_tmp;
    float dx[6];
    int b_dx_tmp;
    int c_dx_tmp;
    float q_idx_0;
    int q_idx_1_tmp;
    int q_idx_2_tmp;
    int q_idx_3_tmp;
    float xg_tmp;
    float out;
    int b_sign;
    float c[3];
    float q[16];
    int b_k;
    float c_sign[12];
    float f;
    float f1;
    float q_error[4];
    int i1;
    float b_xnew[7];
    int unew_tmp;
    float b_fv1[3];

    /*  Uses an rk method to roll out a trajectory */
    /*  Returns the new trajectory, cost and the derivatives along the trajectory */
    /*  If feed-forward and feed-back controls l, K are non-zero */
    /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
    /*  Sizes */
    /*  Initialize outputs (Don't overwrite x or u) */
    memset(&xnew[0], 0, 3500U * sizeof(float));
    memset(&unew[0], 0, 1497U * sizeof(float));
    *cost = 0.0F;
    for (i = 0; i < 7; i++) {
        xnew[i] = x[i];
    }

    b_fv[0] = 0.0F;
    b_fv[4] = 0.0F;
    b_fv[8] = 0.0F;
    for (k = 0; k < 499; k++) {
        /*  Find the state error vector dx */
        dx_tmp = 7 * k + 4;
        dx[3] = xnew[dx_tmp] - x[dx_tmp];
        b_dx_tmp = 7 * k + 5;
        dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
        c_dx_tmp = 7 * k + 6;
        dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

        /*  Calculate error between qnew and q_nom */
        /*  Defined as conj(q_nom)*qnew */
        q_idx_0 = x[7 * k];
        q_idx_1_tmp = 7 * k + 1;
        q_idx_2_tmp = 7 * k + 2;
        q_idx_3_tmp = 7 * k + 3;

        /*  conjugate */
        /*  Forms the "left matrix" for quaternion multiplication */
        /*  where q*q1 = L(q)q1 */
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
        f = xnew[7 * k];
        for (i = 0; i < 4; i++) {
            f1 = ((q[i] * f + q[i + 4] * xnew[q_idx_1_tmp]) + q[i + 8] *
                  xnew[q_idx_2_tmp]) + q[i + 12] * xnew[q_idx_3_tmp];
            q_error[i] = f1;
            q_idx_0 += f1 * f1;
        }

        q_idx_0 = sqrtf(q_idx_0);
        f = q_error[0] / q_idx_0;
        f1 = q_error[1] / q_idx_0;
        xg_tmp = q_error[2] / q_idx_0;
        q_idx_0 = q_error[3] / q_idx_0;

        /*  re-normalize */
        /*  inverse Cayley Map */
        dx[0] = f1 / f;
        dx[1] = xg_tmp / f;
        dx[2] = q_idx_0 / f;

        /*  Find the new control and ensure it is within the limits */
        for (b_k = 0; b_k < 3; b_k++) {
            f = 0.0F;
            for (i = 0; i < 6; i++) {
                f += K[(b_k + 3 * i) + 18 * k] * dx[i];
            }

            unew_tmp = b_k + 3 * k;
            unew[unew_tmp] = fminf(u_lims[b_k + 3], fmaxf(u_lims[b_k],
                                    (u[unew_tmp] - alpha * l[unew_tmp]) - f));
        }

        /*  Step the dynamics forward */
        for (i1 = 0; i1 < 7; i1++) {
            b_xnew[i1] = xnew[i1 + 7 * k];
        }

        satellite_step(b_xnew, *(float (*)[3])&unew[3 * k], *(float (*)[7])&
                       xnew[7 * (k + 1)], *(float (*)[36])&fx[36 * k], *(float (*)
                        [18])&fu[18 * k]);

        /*  Calculate the cost */
        /*  Calculates the cost contribution of a given state and control  */
        /*  Also calculates the 2nd order expansion of the cost function */
        /*  Inputs */
        /* ===================================== */
        /*  x        - [quaternion; omega] (7x1) */
        /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
        /*  terminal - int 0 or 1 */
        /* --------------------------------------------------- */
        /*  control hessian */
        /*  angular velocity hessian */
        /*  Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q)) */
        /*  Also record the sign for use in the backward pass */
        q_idx_0 = xnew[7 * k];
        xg_tmp = ((xg[0] * q_idx_0 + xg[1] * xnew[q_idx_1_tmp]) + xg[2] *
                  xnew[q_idx_2_tmp]) + xg[3] * xnew[q_idx_3_tmp];
        if (xg_tmp + 1.0F < 1.0F - xg_tmp) {
            out = xg_tmp + 1.0F;
            b_sign = 1;
        } else {
            out = 1.0F - xg_tmp;
            b_sign = -1;
        }

        c[0] = xnew[dx_tmp] - xg[4];
        c[1] = xnew[b_dx_tmp] - xg[5];
        c[2] = xnew[c_dx_tmp] - xg[6];

        /*  State cost Hessian */
        for (i = 0; i < 3; i++) {
            b_k = 6 * i + 36 * k;
            cxx[b_k] = (float)(-b_sign * iv[3 * i]) * xg_tmp;
            dx_tmp = 6 * (i + 3) + 36 * k;
            cxx[dx_tmp] = 0.0F;
            cxx[b_k + 1] = (float)(-b_sign * iv[3 * i + 1]) * xg_tmp;
            cxx[dx_tmp + 1] = 0.0F;
            cxx[b_k + 2] = (float)(-b_sign * iv[3 * i + 2]) * xg_tmp;
            cxx[dx_tmp + 2] = 0.0F;
        }

        for (i = 0; i < 6; i++) {
            b_k = 6 * i + 36 * k;
            cxx[b_k + 3] = fv[3 * i];
            cxx[b_k + 4] = fv[3 * i + 1];
            cxx[b_k + 5] = fv[3 * i + 2];
        }

        /*  State cost Jacobian */
        /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
        b_fv[3] = -xnew[q_idx_3_tmp];
        b_fv[6] = xnew[q_idx_2_tmp];
        b_fv[1] = xnew[q_idx_3_tmp];
        b_fv[7] = -xnew[q_idx_1_tmp];
        b_fv[2] = -xnew[q_idx_2_tmp];
        b_fv[5] = xnew[q_idx_1_tmp];
        for (i = 0; i < 3; i++) {
            c_sign[i] = (float)b_sign * -xnew[(i + 7 * k) + 1];
            b_k = 3 * (i + 1);
            c_sign[b_k] = (float)b_sign * (q_idx_0 * (float)iv[i] + b_fv[i]);
            c_sign[b_k + 1] = (float)b_sign * (q_idx_0 * (float)iv[i + 3] +
                b_fv[i + 3]);
            c_sign[b_k + 2] = (float)b_sign * (q_idx_0 * (float)iv[i + 6] +
                b_fv[i + 6]);
        }

        /*  Control cost Hessian & Jacobian */
        f = 0.0F;
        f1 = unew[3 * k];
        i = 3 * k + 1;
        b_dx_tmp = 3 * k + 2;
        for (c_dx_tmp = 0; c_dx_tmp < 3; c_dx_tmp++) {
            b_k = c_dx_tmp + 6 * k;
            cx[b_k] = ((c_sign[c_dx_tmp] * xg[0] + c_sign[c_dx_tmp + 3] * xg[1])
                       + c_sign[c_dx_tmp + 6] * xg[2]) + c_sign[c_dx_tmp + 9] *
                xg[3];
            cx[b_k + 3] = (fv1[c_dx_tmp] * c[0] + fv1[c_dx_tmp + 3] * c[1]) +
                fv1[c_dx_tmp + 6] * c[2];
            b_k = 3 * c_dx_tmp + 9 * k;
            q_idx_0 = fv2[3 * c_dx_tmp];
            cuu[b_k] = q_idx_0;
            dx_tmp = 3 * c_dx_tmp + 1;
            cuu[b_k + 1] = fv2[dx_tmp];
            unew_tmp = 3 * c_dx_tmp + 2;
            cuu[b_k + 2] = fv2[unew_tmp];
            cu[c_dx_tmp + 3 * k] = (fv2[c_dx_tmp] * f1 + fv2[c_dx_tmp + 3] *
                                    unew[i]) + fv2[c_dx_tmp + 6] * unew[b_dx_tmp];
            f += ((0.5F * c[0] * fv1[3 * c_dx_tmp] + 0.5F * c[1] * fv1[dx_tmp])
                  + 0.5F * c[2] * fv1[unew_tmp]) * c[c_dx_tmp];
            b_fv1[c_dx_tmp] = (0.5F * f1 * q_idx_0 + 0.5F * unew[i] * fv2[dx_tmp])
                + 0.5F * unew[b_dx_tmp] * fv2[unew_tmp];
        }

        *cost += (out + f) + ((b_fv1[0] * f1 + b_fv1[1] * unew[3 * k + 1]) +
                              b_fv1[2] * unew[3 * k + 2]);
    }

    /*  Final cost */
    /*  Calculates the cost contribution of a given state and control  */
    /*  Also calculates the 2nd order expansion of the cost function */
    /*  Inputs */
    /* ===================================== */
    /*  x        - [quaternion; omega] (7x1) */
    /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
    /*  terminal - int 0 or 1 */
    /* --------------------------------------------------- */
    /*  control hessian */
    /*  terminal angular velocity hessian  */
    /*  Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q)) */
    /*  Also record the sign for use in the backward pass */
    xg_tmp = ((xg[0] * xnew[3493] + xg[1] * xnew[3494]) + xg[2] * xnew[3495]) +
        xg[3] * xnew[3496];
    if (xg_tmp + 1.0F < 1.0F - xg_tmp) {
        out = xg_tmp + 1.0F;
        b_sign = 1;
    } else {
        out = 1.0F - xg_tmp;
        b_sign = -1;
    }

    c[0] = xnew[3497] - xg[4];
    c[1] = xnew[3498] - xg[5];
    c[2] = xnew[3499] - xg[6];

    /*  State cost Hessian */
    for (i = 0; i < 3; i++) {
        cxx[6 * i + 17964] = (float)(-b_sign * iv[3 * i]) * xg_tmp;
        b_k = 6 * (i + 3);
        cxx[b_k + 17964] = 0.0F;
        cxx[6 * i + 17965] = (float)(-b_sign * iv[3 * i + 1]) * xg_tmp;
        cxx[b_k + 17965] = 0.0F;
        cxx[6 * i + 17966] = (float)(-b_sign * iv[3 * i + 2]) * xg_tmp;
        cxx[b_k + 17966] = 0.0F;
    }

    for (i = 0; i < 6; i++) {
        cxx[6 * i + 17967] = iv1[3 * i];
        cxx[6 * i + 17968] = iv1[3 * i + 1];
        cxx[6 * i + 17969] = iv1[3 * i + 2];
    }

    /*  State cost Jacobian */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    b_fv[0] = 0.0F;
    b_fv[3] = -xnew[3496];
    b_fv[6] = xnew[3495];
    b_fv[1] = xnew[3496];
    b_fv[4] = 0.0F;
    b_fv[7] = -xnew[3494];
    b_fv[2] = -xnew[3495];
    b_fv[5] = xnew[3494];
    b_fv[8] = 0.0F;
    for (i = 0; i < 3; i++) {
        c_sign[i] = (float)b_sign * -xnew[i + 3494];
        b_k = 3 * (i + 1);
        c_sign[b_k] = (float)b_sign * (xnew[3493] * (float)iv[i] + b_fv[i]);
        c_sign[b_k + 1] = (float)b_sign * (xnew[3493] * (float)iv[i + 3] +
            b_fv[i + 3]);
        c_sign[b_k + 2] = (float)b_sign * (xnew[3493] * (float)iv[i + 6] +
            b_fv[i + 6]);
    }

    /*  Control cost Hessian & Jacobian */
    f = 0.0F;
    for (i = 0; i < 3; i++) {
        cx[i + 2994] = ((c_sign[i] * xg[0] + c_sign[i + 3] * xg[1]) + c_sign[i +
                        6] * xg[2]) + c_sign[i + 9] * xg[3];
        cx[i + 2997] = ((float)iv[i] * c[0] + (float)iv[i + 3] * c[1]) + (float)
            iv[i + 6] * c[2];
        f += ((0.5F * c[0] * (float)iv[3 * i] + 0.5F * c[1] * (float)iv[3 * i +
               1]) + 0.5F * c[2] * (float)iv[3 * i + 2]) * c[i];
    }

    *cost += out + f;
}

/*
 * Arguments    : const float fx[17964]
 *                const float fu[8982]
 *                const float cx[3000]
 *                const float cu[1497]
 *                const float cxx[18000]
 *                const float cuu[4491]
 *                float lambda
 *                const float u_lims[6]
 *                const float u[1497]
 *                float l[1497]
 *                float K[8982]
 *                float dV[2]
 *                bool *diverge
 * Return Type  : void
 */
static void backwardPass(const float fx[17964], const float fu[8982], const
    float cx[3000], const float cu[1497], const float cxx[18000], const float
    cuu[4491], float lambda, const float u_lims[6], const float u[1497], float
    l[1497], float K[8982], float dV[2], bool *diverge)
{
    int i;
    float Vx[6];
    int k;
    int i1;
    bool exitg1;
    int b_k;
    float Vxx[36];
    float f;
    float Qx_tmp[36];
    float Qu_tmp[18];
    float Qu[3];
    int i2;
    int u_lims_tmp;
    float Quu[9];
    float b_Quu[9];
    float b_u_lims[3];
    float coeffs;
    float c_u_lims[3];
    int i3;
    float b_fv[3];
    float Quu_tmp[18];
    int b_u_lims_tmp;
    float Qux[18];
    int c_u_lims_tmp;
    float lk[3];
    float result;
    float Luu_data[9];
    int Luu_size[2];
    bool b_free[3];
    float Kk[18];
    bool out;
    bool exitg2;
    signed char tmp_data[3];
    float dV_idx_0;
    float b_Kk[18];
    int value;
    int Y_size_idx_0;
    int loop_ub;
    float Y_data[18];
    int coeffs_size_idx_0;
    float b_cx[6];
    float b_Quu_tmp[6];
    float coeffs_data[3];
    float b_Qx_tmp[36];
    int b_i;
    int X_size_idx_0;
    int j;
    float f1;
    float b_cxx[36];
    int LT_size_idx_0;
    int c_i;
    float LT_data[9];
    signed char b_tmp_data[3];

    /*  Perfoms the LQR backward pass to find the optimal controls */
    /*  Solves a quadratic program (QP) at each timestep for the optimal */
    /*  controls given the control limits */
    /*  Initialize matrices (for C) */
    memset(&l[0], 0, 1497U * sizeof(float));
    memset(&K[0], 0, 8982U * sizeof(float));

    /*  Change in cost */
    dV[0] = 0.0F;
    dV[1] = 0.0F;

    /*  Set cost-to-go Jacobain and Hessian equal to final costs */
    for (i = 0; i < 6; i++) {
        Vx[i] = cx[i + 2994];
        for (i1 = 0; i1 < 6; i1++) {
            b_k = i1 + 6 * i;
            Vxx[b_k] = cxx[b_k + 17964];
        }
    }

    *diverge = false;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 499)) {
        /*  Define cost gradients */
        for (i = 0; i < 6; i++) {
            for (i1 = 0; i1 < 6; i1++) {
                Qx_tmp[i1 + 6 * i] = fx[(i + 6 * i1) + 36 * (498 - k)];
            }

            b_k = i + 18 * (498 - k);
            Qu_tmp[3 * i] = fu[b_k];
            Qu_tmp[3 * i + 1] = fu[b_k + 6];
            Qu_tmp[3 * i + 2] = fu[b_k + 12];
        }

        for (i = 0; i < 3; i++) {
            f = 0.0F;
            for (i1 = 0; i1 < 6; i1++) {
                i2 = i + 3 * i1;
                f += Qu_tmp[i2] * Vx[i1];
                coeffs = 0.0F;
                for (i3 = 0; i3 < 6; i3++) {
                    coeffs += Qu_tmp[i + 3 * i3] * Vxx[i3 + 6 * i1];
                }

                Quu_tmp[i2] = coeffs;
            }

            Qu[i] = cu[i + 3 * (498 - k)] + f;
            for (i1 = 0; i1 < 3; i1++) {
                f = 0.0F;
                for (i2 = 0; i2 < 6; i2++) {
                    f += Quu_tmp[i + 3 * i2] * fu[(i2 + 6 * i1) + 18 * (498 - k)];
                }

                b_k = i + 3 * i1;
                b_Quu[b_k] = cuu[b_k + 9 * (498 - k)] + f;
            }

            for (i1 = 0; i1 < 6; i1++) {
                f = 0.0F;
                for (i2 = 0; i2 < 6; i2++) {
                    f += Quu_tmp[i + 3 * i2] * fx[(i2 + 6 * i1) + 36 * (498 - k)];
                }

                Qux[i + 3 * i1] = f;
            }
        }

        /*  Regularization (for Cholesky positive definiteness) */
        /*  Solve the Quadratic program with control limits */
        for (i = 0; i < 9; i++) {
            Quu[i] = b_Quu[i] + (float)iv[i] * lambda;
        }

        u_lims_tmp = 3 * (498 - k);
        b_u_lims[0] = u_lims[0] - u[u_lims_tmp];
        c_u_lims[0] = u_lims[3] - u[u_lims_tmp];
        i = 3 * ((int)fminf(499.0F, (-(float)k + 499.0F) + 1.0F) - 1);
        b_fv[0] = -l[i];
        b_u_lims_tmp = u_lims_tmp + 1;
        b_u_lims[1] = u_lims[1] - u[b_u_lims_tmp];
        c_u_lims[1] = u_lims[4] - u[b_u_lims_tmp];
        b_fv[1] = -l[i + 1];
        c_u_lims_tmp = u_lims_tmp + 2;
        b_u_lims[2] = u_lims[2] - u[c_u_lims_tmp];
        c_u_lims[2] = u_lims[5] - u[c_u_lims_tmp];
        b_fv[2] = -l[i + 2];
        boxQPsolve(Quu, Qu, b_u_lims, c_u_lims, b_fv, lk, &result, Luu_data,
                   Luu_size, b_free);
        if (result < 1.0F) {
            *diverge = true;
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
                /*  Solution is found in O(nm) time using back-substitution */
                /*  This implementation only works for lower triangular factorisations */
                value = Luu_size[0] - 1;

                /*  Solve LY = B for Y */
                /* ======================= */
                Y_size_idx_0 = Luu_size[0];
                loop_ub = Luu_size[0] * 6;
                if (0 <= loop_ub - 1) {
                    memset(&Y_data[0], 0, loop_ub * sizeof(float));
                }

                coeffs_size_idx_0 = Luu_size[0];
                loop_ub = Luu_size[0];
                if (0 <= loop_ub - 1) {
                    memset(&coeffs_data[0], 0, loop_ub * sizeof(float));
                }

                i = Luu_size[0];
                for (b_i = 0; b_i < i; b_i++) {
                    if (b_i + 1 != 1) {
                        for (i1 = 0; i1 < b_i; i1++) {
                            coeffs_data[i1] = Luu_data[b_i + Luu_size[0] * i1];
                        }
                    }

                    for (j = 0; j < 6; j++) {
                        if (((signed char)coeffs_size_idx_0 == 1) ||
                                (Y_size_idx_0 == 1)) {
                            coeffs = 0.0F;
                            for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                coeffs += coeffs_data[i1] * Y_data[i1 +
                                    Y_size_idx_0 * j];
                            }
                        } else {
                            coeffs = 0.0F;
                            for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                coeffs += coeffs_data[i1] * Y_data[i1 +
                                    Y_size_idx_0 * j];
                            }
                        }

                        Y_data[b_i + Y_size_idx_0 * j] = (Qux[(tmp_data[b_i] + 3
                            * j) - 1] - coeffs) / Luu_data[b_i + Luu_size[0] *
                            b_i];
                    }
                }

                /*  Solve L'X = Y for X */
                /* ======================= */
                X_size_idx_0 = Luu_size[0];
                loop_ub = Luu_size[0] * 6;
                if (0 <= loop_ub - 1) {
                    memset(&Qu_tmp[0], 0, loop_ub * sizeof(float));
                }

                LT_size_idx_0 = Luu_size[1];
                loop_ub = Luu_size[0];
                coeffs_size_idx_0 = Luu_size[0];
                for (i = 0; i < loop_ub; i++) {
                    b_k = Luu_size[1];
                    for (i1 = 0; i1 < b_k; i1++) {
                        LT_data[i1 + LT_size_idx_0 * i] = Luu_data[i + Luu_size
                            [0] * i1];
                    }

                    coeffs_data[i] = 0.0F;
                }

                i = (int)(((-1.0F - (float)Luu_size[0]) + 1.0F) / -1.0F);
                for (b_i = 0; b_i < i; b_i++) {
                    c_i = value - b_i;
                    if (c_i + 1 != value + 1) {
                        if (c_i + 2 > Luu_size[0]) {
                            i1 = 0;
                            i2 = -1;
                            i3 = 0;
                        } else {
                            i1 = c_i + 1;
                            i2 = value;
                            i3 = c_i + 1;
                        }

                        loop_ub = i2 - i1;
                        for (i2 = 0; i2 <= loop_ub; i2++) {
                            coeffs_data[i3 + i2] = LT_data[c_i + LT_size_idx_0 *
                                (i1 + i2)];
                        }
                    }

                    for (j = 0; j < 6; j++) {
                        if (((signed char)coeffs_size_idx_0 == 1) ||
                                (X_size_idx_0 == 1)) {
                            coeffs = 0.0F;
                            for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                coeffs += coeffs_data[i1] * Qu_tmp[i1 +
                                    X_size_idx_0 * j];
                            }
                        } else {
                            coeffs = 0.0F;
                            for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                coeffs += coeffs_data[i1] * Qu_tmp[i1 +
                                    X_size_idx_0 * j];
                            }
                        }

                        Qu_tmp[c_i + X_size_idx_0 * j] = (Y_data[c_i +
                            Y_size_idx_0 * j] - coeffs) / LT_data[c_i +
                            LT_size_idx_0 * c_i];
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

                for (i = 0; i < 6; i++) {
                    for (i1 = 0; i1 < X_size_idx_0; i1++) {
                        Kk[(b_tmp_data[i1] + 3 * i) - 1] = -Qu_tmp[i1 +
                            X_size_idx_0 * i];
                    }
                }
            }

            /*  Update Cost to Go */
            f = 0.0F;
            for (i = 0; i < 3; i++) {
                f += ((0.5F * lk[0] * b_Quu[3 * i] + 0.5F * lk[1] * b_Quu[3 * i
                       + 1]) + 0.5F * lk[2] * b_Quu[3 * i + 2]) * lk[i];
            }

            dV_idx_0 = dV[0] + ((lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2]);
            coeffs = dV[1] + f;
            dV[0] = dV_idx_0;
            dV[1] = coeffs;
            for (i = 0; i < 3; i++) {
                for (i1 = 0; i1 < 6; i1++) {
                    b_Kk[i1 + 6 * i] = Kk[i + 3 * i1];
                }
            }

            memcpy(&Qu_tmp[0], &b_Kk[0], 18U * sizeof(float));
            for (i = 0; i < 6; i++) {
                f = Qu_tmp[i + 6];
                coeffs = Qu_tmp[i + 12];
                for (i1 = 0; i1 < 3; i1++) {
                    b_Kk[i + 6 * i1] = (Qu_tmp[i] * b_Quu[3 * i1] + f * b_Quu[3 *
                                        i1 + 1]) + coeffs * b_Quu[3 * i1 + 2];
                }
            }

            memcpy(&Quu_tmp[0], &b_Kk[0], 18U * sizeof(float));
            for (i = 0; i < 3; i++) {
                for (i1 = 0; i1 < 6; i1++) {
                    b_Kk[i1 + 6 * i] = Qux[i + 3 * i1];
                }
            }

            for (i = 0; i < 6; i++) {
                f = 0.0F;
                for (i1 = 0; i1 < 6; i1++) {
                    f += Qx_tmp[i + 6 * i1] * Vx[i1];
                }

                b_cx[i] = ((cx[i + 6 * (498 - k)] + f) + ((Quu_tmp[i] * lk[0] +
                             Quu_tmp[i + 6] * lk[1]) + Quu_tmp[i + 12] * lk[2]))
                    + ((Qu_tmp[i] * Qu[0] + Qu_tmp[i + 6] * Qu[1]) + Qu_tmp[i +
                       12] * Qu[2]);
                b_Quu_tmp[i] = (b_Kk[i] * lk[0] + b_Kk[i + 6] * lk[1]) + b_Kk[i
                    + 12] * lk[2];
            }

            for (i = 0; i < 6; i++) {
                Vx[i] = b_cx[i] + b_Quu_tmp[i];
                for (i1 = 0; i1 < 6; i1++) {
                    f = 0.0F;
                    for (i2 = 0; i2 < 6; i2++) {
                        f += Qx_tmp[i + 6 * i2] * Vxx[i2 + 6 * i1];
                    }

                    b_Qx_tmp[i + 6 * i1] = f;
                }

                for (i1 = 0; i1 < 6; i1++) {
                    f = 0.0F;
                    for (i2 = 0; i2 < 6; i2++) {
                        f += b_Qx_tmp[i + 6 * i2] * fx[(i2 + 6 * i1) + 36 * (498
                            - k)];
                    }

                    b_k = i + 6 * i1;
                    b_cxx[b_k] = cxx[b_k + 36 * (498 - k)] + f;
                    Qx_tmp[b_k] = (Quu_tmp[i] * Kk[3 * i1] + Quu_tmp[i + 6] *
                                   Kk[3 * i1 + 1]) + Quu_tmp[i + 12] * Kk[3 * i1
                        + 2];
                }

                f = Qu_tmp[i + 6];
                coeffs = Qu_tmp[i + 12];
                dV_idx_0 = b_Kk[i + 6];
                f1 = b_Kk[i + 12];
                for (i1 = 0; i1 < 6; i1++) {
                    i2 = 3 * i1 + 1;
                    i3 = 3 * i1 + 2;
                    b_k = i + 6 * i1;
                    b_Qx_tmp[b_k] = (b_cxx[b_k] + Qx_tmp[b_k]) + ((Qu_tmp[i] *
                        Qux[3 * i1] + f * Qux[i2]) + coeffs * Qux[i3]);
                    Qx_tmp[b_k] = (b_Kk[i] * Kk[3 * i1] + dV_idx_0 * Kk[i2]) +
                        f1 * Kk[i3];
                }
            }

            for (i = 0; i < 36; i++) {
                Vxx[i] = b_Qx_tmp[i] + Qx_tmp[i];
            }

            for (i = 0; i < 6; i++) {
                for (i1 = 0; i1 < 6; i1++) {
                    b_k = i1 + 6 * i;
                    Qx_tmp[b_k] = 0.5F * (Vxx[b_k] + Vxx[i + 6 * i1]);
                }
            }

            memcpy(&Vxx[0], &Qx_tmp[0], 36U * sizeof(float));

            /*  Make sure Hessian is symmetric */
            /*  Update Control Vectors */
            l[u_lims_tmp] = -lk[0];
            l[b_u_lims_tmp] = -lk[1];
            l[c_u_lims_tmp] = -lk[2];
            for (i = 0; i < 6; i++) {
                b_k = 3 * i + 18 * (498 - k);
                K[b_k] = -Kk[3 * i];
                K[b_k + 1] = -Kk[3 * i + 1];
                K[b_k + 2] = -Kk[3 * i + 2];
            }

            k++;
        }
    }
}

/*
 * Arguments    : const float Quu[9]
 *                const float Qu[3]
 *                const float lower[3]
 *                const float upper[3]
 *                const float u0[3]
 *                float u[3]
 *                float *result
 *                float Luu_data[]
 *                int Luu_size[2]
 *                bool b_free[3]
 * Return Type  : void
 */
static void boxQPsolve(const float Quu[9], const float Qu[3], const float lower
                       [3], const float upper[3], const float u0[3], float u[3],
                       float *result, float Luu_data[], int Luu_size[2],
                       bool b_free[3])
{
    bool clamped[3];
    int i;
    float oldvalue;
    float scale;
    float z1[3];
    float absxk;
    float t;
    float value;
    int iter;
    int b_iter;
    bool exitg1;
    int k;
    bool out;
    float grad[3];
    bool prev_clamped[3];
    bool exitg2;
    bool factorize;
    bool guard1 = false;
    int trueCount;
    signed char tmp_data[3];
    signed char b_tmp_data[3];
    int Quu_size[2];
    float Quu_data[9];
    float indef;
    int i1;
    float fcnOutput;
    float gnorm;
    float grad_clamped[3];
    float deltaX[3];
    signed char c_tmp_data[3];
    int loop_ub;
    int coeffs_size_idx_0;
    float Y_data[3];
    float coeffs_data[3];
    int X_size_idx_0;
    int LT_size_idx_0;
    int b_i;
    float d_tmp_data[3];
    int i2;
    int i3;
    float LT_data[9];
    float step;
    float uc_idx_0;
    float uc_idx_1;
    float uc_idx_2;
    float vc;

    /*  Finds the optimal control with limits to minimize a quadratic cost */
    /*  Minimize 0.5*u'*Quu*u + u'*Qu  s.t. lower <= u <= upper */
    /*  */
    /*   inputs: */
    /*      Quu       - positive definite matrix              (m * m) */
    /*      Qu        - bias vector                           (m) */
    /*      lower     - lower bounds                          (m) */
    /*      upper     - upper bounds                          (m) */
    /*      u0        - initial control input for warm-start  (m) */
    /*  */
    /*   outputs: */
    /*      u         - solution                   (m) */
    /*      result    - result type (roughly, higher is better, see below) */
    /*      Luu       - cholesky factor            (m * m) */
    /*      free      - set of free dimensions     (m) */
    /*  Initialize arrays */
    clamped[0] = false;
    b_free[0] = true;
    clamped[1] = false;
    b_free[1] = true;
    clamped[2] = false;
    b_free[2] = true;
    Luu_size[0] = 3;
    Luu_size[1] = 3;
    for (i = 0; i < 9; i++) {
        Luu_data[i] = 0.0F;
    }

    /*  Placeholder to return if Luu not assigned */
    /*  Initialize scalars */
    oldvalue = 0.0F;
    *result = 0.0F;

    /*  options */
    /*  maximum number of iterations */
    /*  minimum norm of non-fixed gradient */
    /*  minimum relative improvement */
    /*  factor for decreasing stepsize */
    /*  minimal stepsize for linesearch */
    /*  Armijo parameter (fraction of linear improvement required) */
    /*  initial state */
    /*  Returns array x with all values clamped between lower and upper */
    /*  initial objective value */
    scale = fmaxf(lower[0], fminf(upper[0], u0[0]));
    z1[0] = scale;
    u[0] = scale;
    absxk = scale * Qu[0];
    scale = fmaxf(lower[1], fminf(upper[1], u0[1]));
    z1[1] = scale;
    u[1] = scale;
    absxk += scale * Qu[1];
    scale = fmaxf(lower[2], fminf(upper[2], u0[2]));
    z1[2] = scale;
    u[2] = scale;
    absxk += scale * Qu[2];
    t = 0.0F;
    for (i = 0; i < 3; i++) {
        t += ((0.5F * z1[0] * Quu[3 * i] + 0.5F * z1[1] * Quu[3 * i + 1]) + 0.5F
              * scale * Quu[3 * i + 2]) * z1[i];
    }

    value = absxk + t;

    /*  main loop */
    iter = 1;
    b_iter = 1;
    exitg1 = false;
    while ((!exitg1) && (b_iter - 1 < 100)) {
        iter = b_iter;
        if (*result != 0.0F) {
            exitg1 = true;
        } else {
            /*  check relative improvement */
            if ((b_iter > 1) && (oldvalue - value < 1.0E-8F * fabsf(oldvalue)))
            {
                *result = 3.0F;
                exitg1 = true;
            } else {
                oldvalue = value;

                /*  get gradient */
                /*  find clamped dimensions */
                for (k = 0; k < 3; k++) {
                    scale = Qu[k] + ((Quu[k] * u[0] + Quu[k + 3] * u[1]) + Quu[k
                                     + 6] * u[2]);
                    grad[k] = scale;
                    prev_clamped[k] = clamped[k];
                    out = false;
                    clamped[k] = false;
                    if ((u[k] == lower[k]) && (scale > 0.0F)) {
                        out = true;
                        clamped[k] = true;
                    }

                    if ((u[k] == upper[k]) && (scale < 0.0F)) {
                        out = true;
                        clamped[k] = true;
                    }

                    b_free[k] = !out;
                }

                /*  check for all clamped */
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
                    *result = 5.0F;
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

                    /*  Cholesky (check for non PD) */
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
                        for (i = 0; i < trueCount; i++) {
                            for (i1 = 0; i1 < trueCount; i1++) {
                                Quu_data[i1 + trueCount * i] = Quu[(tmp_data[i1]
                                    + 3 * (tmp_data[i] - 1)) - 1];
                            }
                        }

                        chol_free(Quu_data, Quu_size, Luu_data, Luu_size, &indef);
                        if (indef != 0.0F) {
                            *result = -1.0F;
                            exitg1 = true;
                        } else {
                            guard1 = true;
                        }
                    } else {
                        guard1 = true;
                    }

                    if (guard1) {
                        /*  check gradient norm */
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
                            gnorm = 0.0F;
                        } else {
                            fcnOutput = 0.0F;
                            if (trueCount == 1) {
                                gnorm = fabsf(grad[b_tmp_data[0] - 1]);
                            } else {
                                scale = 1.29246971E-26F;
                                for (k = 0; k < trueCount; k++) {
                                    absxk = fabsf(grad[b_tmp_data[k] - 1]);
                                    if (absxk > scale) {
                                        t = scale / absxk;
                                        fcnOutput = fcnOutput * t * t + 1.0F;
                                        scale = absxk;
                                    } else {
                                        t = absxk / scale;
                                        fcnOutput += t * t;
                                    }
                                }

                                gnorm = scale * sqrtf(fcnOutput);
                            }
                        }

                        if (gnorm < 1.0E-8F) {
                            *result = 4.0F;
                            exitg1 = true;
                        } else {
                            /*  get search direction */
                            scale = u[0] * (float)clamped[0];
                            absxk = u[1] * (float)clamped[1];
                            t = u[2] * (float)clamped[2];
                            for (k = 0; k < 3; k++) {
                                grad_clamped[k] = Qu[k] + ((Quu[k] * scale +
                                    Quu[k + 3] * absxk) + Quu[k + 6] * t);
                                deltaX[k] = 0.0F;
                            }

                            k = 0;
                            if (b_free[0]) {
                                c_tmp_data[0] = 1;
                                k = 1;
                            }

                            if (b_free[1]) {
                                c_tmp_data[k] = 2;
                                k++;
                            }

                            if (b_free[2]) {
                                c_tmp_data[k] = 3;
                            }

                            /*  Solves the linear system AX = B for the unknown X using the Cholesky */
                            /*  decomposition L of the matrix A. Where LL' = A */
                            /*  X can be a vector or a matrix of size n x m */
                            /*  Solution is found in O(nm) time using back-substitution */
                            /*  This implementation only works for lower triangular factorisations */
                            trueCount = Luu_size[0] - 1;

                            /*  Solve LY = B for Y */
                            /* ======================= */
                            loop_ub = Luu_size[0];
                            coeffs_size_idx_0 = Luu_size[0];
                            if (0 <= loop_ub - 1) {
                                memset(&Y_data[0], 0, loop_ub * sizeof(float));
                                memset(&coeffs_data[0], 0, loop_ub * sizeof
                                       (float));
                            }

                            i = Luu_size[0];

                            /*  Solve L'X = Y for X */
                            /* ======================= */
                            X_size_idx_0 = Luu_size[0];
                            LT_size_idx_0 = Luu_size[1];
                            for (k = 0; k < i; k++) {
                                if (k + 1 != 1) {
                                    for (i1 = 0; i1 < k; i1++) {
                                        coeffs_data[i1] = Luu_data[k + Luu_size
                                            [0] * i1];
                                    }
                                }

                                if (((signed char)coeffs_size_idx_0 == 1) ||
                                        (Luu_size[0] == 1)) {
                                    scale = 0.0F;
                                    for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                        scale += coeffs_data[i1] * Y_data[i1];
                                    }
                                } else {
                                    scale = 0.0F;
                                    for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                        scale += coeffs_data[i1] * Y_data[i1];
                                    }
                                }

                                Y_data[k] = (grad_clamped[c_tmp_data[k] - 1] -
                                             scale) / Luu_data[k + Luu_size[0] *
                                    k];
                                z1[k] = 0.0F;
                                loop_ub = Luu_size[1];
                                for (i1 = 0; i1 < loop_ub; i1++) {
                                    LT_data[i1 + LT_size_idx_0 * k] = Luu_data[k
                                        + Luu_size[0] * i1];
                                }
                            }

                            coeffs_size_idx_0 = Luu_size[0];
                            loop_ub = Luu_size[0];
                            if (0 <= loop_ub - 1) {
                                memset(&coeffs_data[0], 0, loop_ub * sizeof
                                       (float));
                            }

                            i = (int)(((-1.0F - (float)Luu_size[0]) + 1.0F) /
                                      -1.0F);
                            for (k = 0; k < i; k++) {
                                b_i = trueCount - k;
                                if (b_i + 1 != trueCount + 1) {
                                    if (b_i + 2 > Luu_size[0]) {
                                        i1 = 0;
                                        i2 = -1;
                                        i3 = 0;
                                    } else {
                                        i1 = b_i + 1;
                                        i2 = trueCount;
                                        i3 = b_i + 1;
                                    }

                                    loop_ub = i2 - i1;
                                    for (i2 = 0; i2 <= loop_ub; i2++) {
                                        coeffs_data[i3 + i2] = LT_data[b_i +
                                            LT_size_idx_0 * (i1 + i2)];
                                    }
                                }

                                if (((signed char)coeffs_size_idx_0 == 1) ||
                                        (X_size_idx_0 == 1)) {
                                    scale = 0.0F;
                                    for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                        scale += coeffs_data[i1] * z1[i1];
                                    }
                                } else {
                                    scale = 0.0F;
                                    for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                                        scale += coeffs_data[i1] * z1[i1];
                                    }
                                }

                                z1[b_i] = (Y_data[b_i] - scale) / LT_data[b_i +
                                    LT_size_idx_0 * b_i];
                            }

                            for (i = 0; i < X_size_idx_0; i++) {
                                d_tmp_data[i] = -z1[i] - u[c_tmp_data[i] - 1];
                            }

                            k = 0;

                            /*  cholesky solver */
                            /*  check for descent direction */
                            if (b_free[0]) {
                                deltaX[0] = d_tmp_data[0];
                                k = 1;
                            }

                            if (b_free[1]) {
                                deltaX[1] = d_tmp_data[k];
                                k++;
                            }

                            if (b_free[2]) {
                                deltaX[2] = d_tmp_data[k];
                            }

                            fcnOutput = (deltaX[0] * grad[0] + deltaX[1] * grad
                                         [1]) + deltaX[2] * grad[2];
                            if (fcnOutput >= 0.0F) {
                                /*  (should not happen) */
                                exitg1 = true;
                            } else {
                                /*  Armijo linesearch */
                                step = 1.0F;

                                /*  Returns array x with all values clamped between lower and upper */
                                scale = fmaxf(lower[0], fminf(upper[0], u[0] +
                                               deltaX[0]));
                                z1[0] = scale;
                                uc_idx_0 = scale;
                                absxk = scale * Qu[0];
                                scale = fmaxf(lower[1], fminf(upper[1], u[1] +
                                               deltaX[1]));
                                z1[1] = scale;
                                uc_idx_1 = scale;
                                absxk += scale * Qu[1];
                                scale = fmaxf(lower[2], fminf(upper[2], u[2] +
                                               deltaX[2]));
                                z1[2] = scale;
                                uc_idx_2 = scale;
                                absxk += scale * Qu[2];
                                t = 0.0F;
                                for (i = 0; i < 3; i++) {
                                    t += ((0.5F * z1[0] * Quu[3 * i] + 0.5F *
                                           z1[1] * Quu[3 * i + 1]) + 0.5F *
                                          scale * Quu[3 * i + 2]) * z1[i];
                                }

                                vc = absxk + t;
                                exitg2 = false;
                                while ((!exitg2) && ((vc - value) / (step *
                                         fcnOutput) < 0.1F)) {
                                    step *= 0.6F;

                                    /*  Returns array x with all values clamped between lower and upper */
                                    scale = fmaxf(lower[0], fminf(upper[0], u[0]
                                                   + step * deltaX[0]));
                                    z1[0] = scale;
                                    uc_idx_0 = scale;
                                    absxk = scale * Qu[0];
                                    scale = fmaxf(lower[1], fminf(upper[1], u[1]
                                                   + step * deltaX[1]));
                                    z1[1] = scale;
                                    uc_idx_1 = scale;
                                    absxk += scale * Qu[1];
                                    scale = fmaxf(lower[2], fminf(upper[2], u[2]
                                                   + step * deltaX[2]));
                                    z1[2] = scale;
                                    uc_idx_2 = scale;
                                    absxk += scale * Qu[2];
                                    t = 0.0F;
                                    for (i = 0; i < 3; i++) {
                                        t += ((0.5F * z1[0] * Quu[3 * i] + 0.5F *
                                               z1[1] * Quu[3 * i + 1]) + 0.5F *
                                              scale * Quu[3 * i + 2]) * z1[i];
                                    }

                                    vc = absxk + t;
                                    if (step < 1.0E-20F) {
                                        *result = 2.0F;
                                        exitg2 = true;
                                    }
                                }

                                /*  accept candidate */
                                u[0] = uc_idx_0;
                                u[1] = uc_idx_1;
                                u[2] = uc_idx_2;
                                value = vc;
                                b_iter++;
                            }
                        }
                    }
                }
            }
        }
    }

    if (iter >= 100) {
        *result = 1.0F;
    }

    /*  Results */
    /*  =========================== */
    /*  -1: Hessian is not positive definite */
    /*   0: No descent direction found          (SHOULD NOT OCCUR) */
    /*   1: Maximum main iterations exceeded        */
    /*   2: Maximum line-search iterations exceeded   */
    /*   3: Improvement smaller than tolerance      */
    /*   4: Gradient norm smaller than tolerance     */
    /*   5: All dimensions are clamped  */
}

/*
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
    int jmax;
    int idxAjj;
    float b_A_data[9];
    int n;
    int info;
    int b_info;
    int j;
    bool exitg1;
    int idxA1j;
    float ssq;
    int ix;
    int iy;
    int i;
    int nmj;
    int i1;
    int ia0;
    float c_A_data[9];
    float c;

    /*  Wrapper for MATLAB chol for use with auto coder */
    /*  Inputs: */
    /* =========== */
    /*  A - positive semi-definite matrix */
    A_size_idx_0 = A_size[0];
    jmax = A_size[1];
    idxAjj = A_size[0] * A_size[1];
    if (0 <= idxAjj - 1) {
        memcpy(&b_A_data[0], &A_data[0], idxAjj * sizeof(float));
    }

    n = A_size[1];
    info = 0;
    if (A_size[1] != 0) {
        b_info = -1;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j <= n - 1)) {
            idxA1j = j * n;
            idxAjj = idxA1j + j;
            ssq = 0.0F;
            if (j >= 1) {
                ix = idxA1j;
                iy = idxA1j;
                for (jmax = 0; jmax < j; jmax++) {
                    ssq += b_A_data[ix] * b_A_data[iy];
                    ix++;
                    iy++;
                }
            }

            ssq = b_A_data[idxAjj] - ssq;
            if (ssq > 0.0F) {
                ssq = sqrtf(ssq);
                b_A_data[idxAjj] = ssq;
                if (j + 1 < n) {
                    nmj = (n - j) - 2;
                    ia0 = (idxA1j + n) + 1;
                    idxAjj += n;
                    if ((j != 0) && (nmj + 1 != 0)) {
                        iy = idxAjj;
                        i = ia0 + n * nmj;
                        for (jmax = ia0; n < 0 ? jmax >= i : jmax <= i; jmax +=
                                n) {
                            ix = idxA1j;
                            c = 0.0F;
                            i1 = (jmax + j) - 1;
                            for (info = jmax; info <= i1; info++) {
                                c += b_A_data[info - 1] * b_A_data[ix];
                                ix++;
                            }

                            b_A_data[iy] += -c;
                            iy += n;
                        }
                    }

                    ssq = 1.0F / ssq;
                    i = (idxAjj + n * nmj) + 1;
                    for (jmax = idxAjj + 1; n < 0 ? jmax >= i : jmax <= i; jmax +=
                         n) {
                        b_A_data[jmax - 1] *= ssq;
                    }
                }

                j++;
            } else {
                b_A_data[idxAjj] = ssq;
                b_info = j;
                exitg1 = true;
            }
        }

        info = b_info + 1;
        if (b_info + 1 == 0) {
            jmax = A_size[1];
        } else {
            jmax = b_info;
        }

        for (j = 0; j < jmax; j++) {
            i = j + 2;
            for (idxAjj = i; idxAjj <= jmax; idxAjj++) {
                b_A_data[(idxAjj + A_size_idx_0 * j) - 1] = 0.0F;
            }
        }

        if (1 > jmax) {
            idxAjj = 0;
            jmax = 0;
        } else {
            idxAjj = jmax;
        }

        for (i = 0; i < jmax; i++) {
            for (i1 = 0; i1 < idxAjj; i1++) {
                c_A_data[i1 + idxAjj * i] = b_A_data[i1 + A_size_idx_0 * i];
            }
        }

        A_size_idx_0 = idxAjj;
        idxAjj *= jmax;
        if (0 <= idxAjj - 1) {
            memcpy(&b_A_data[0], &c_A_data[0], idxAjj * sizeof(float));
        }
    }

    *fail = (float)info;
    L_size[0] = A_size_idx_0;
    L_size[1] = jmax;
    idxAjj = A_size_idx_0 * jmax;
    if (0 <= idxAjj - 1) {
        memcpy(&L_data[0], &b_A_data[0], idxAjj * sizeof(float));
    }
}

/*
 * Arguments    : const float x[3500]
 *                const float xg[7]
 *                const float u[1497]
 *                const float K[8982]
 *                const float u_lims[6]
 *                float xnew[3500]
 *                float unew[1497]
 *                float fx[17964]
 *                float fu[8982]
 *                float cx[3000]
 *                float cu[1497]
 *                float cxx[18000]
 *                float cuu[4491]
 *                float *cost
 * Return Type  : void
 */
static void forwardRollout(const float x[3500], const float xg[7], const float
    u[1497], const float K[8982], const float u_lims[6], float xnew[3500], float
    unew[1497], float fx[17964], float fu[8982], float cx[3000], float cu[1497],
    float cxx[18000], float cuu[4491], float *cost)
{
    int i;
    float b_fv[9];
    int k;
    int dx_tmp;
    float dx[6];
    int b_dx_tmp;
    int c_dx_tmp;
    float q_idx_0;
    int q_idx_1_tmp;
    int q_idx_2_tmp;
    int q_idx_3_tmp;
    float xg_tmp;
    float out;
    int b_sign;
    float c[3];
    float q[16];
    int b_k;
    float c_sign[12];
    float f;
    float f1;
    float q_error[4];
    int i1;
    float b_xnew[7];
    int unew_tmp;
    float b_fv1[3];

    /*  Uses an rk method to roll out a trajectory */
    /*  Returns the new trajectory, cost and the derivatives along the trajectory */
    /*  If feed-forward and feed-back controls l, K are non-zero */
    /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
    /*  Sizes */
    /*  Initialize outputs (Don't overwrite x or u) */
    memset(&xnew[0], 0, 3500U * sizeof(float));
    memset(&unew[0], 0, 1497U * sizeof(float));
    *cost = 0.0F;
    for (i = 0; i < 7; i++) {
        xnew[i] = x[i];
    }

    b_fv[0] = 0.0F;
    b_fv[4] = 0.0F;
    b_fv[8] = 0.0F;
    for (k = 0; k < 499; k++) {
        /*  Find the state error vector dx */
        dx_tmp = 7 * k + 4;
        dx[3] = xnew[dx_tmp] - x[dx_tmp];
        b_dx_tmp = 7 * k + 5;
        dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
        c_dx_tmp = 7 * k + 6;
        dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

        /*  Calculate error between qnew and q_nom */
        /*  Defined as conj(q_nom)*qnew */
        q_idx_0 = x[7 * k];
        q_idx_1_tmp = 7 * k + 1;
        q_idx_2_tmp = 7 * k + 2;
        q_idx_3_tmp = 7 * k + 3;

        /*  conjugate */
        /*  Forms the "left matrix" for quaternion multiplication */
        /*  where q*q1 = L(q)q1 */
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
        f = xnew[7 * k];
        for (i = 0; i < 4; i++) {
            f1 = ((q[i] * f + q[i + 4] * xnew[q_idx_1_tmp]) + q[i + 8] *
                  xnew[q_idx_2_tmp]) + q[i + 12] * xnew[q_idx_3_tmp];
            q_error[i] = f1;
            q_idx_0 += f1 * f1;
        }

        q_idx_0 = sqrtf(q_idx_0);
        f = q_error[0] / q_idx_0;
        f1 = q_error[1] / q_idx_0;
        xg_tmp = q_error[2] / q_idx_0;
        q_idx_0 = q_error[3] / q_idx_0;

        /*  re-normalize */
        /*  inverse Cayley Map */
        dx[0] = f1 / f;
        dx[1] = xg_tmp / f;
        dx[2] = q_idx_0 / f;

        /*  Find the new control and ensure it is within the limits */
        for (b_k = 0; b_k < 3; b_k++) {
            f = 0.0F;
            for (i = 0; i < 6; i++) {
                f += K[(b_k + 3 * i) + 18 * k] * dx[i];
            }

            unew_tmp = b_k + 3 * k;
            unew[unew_tmp] = fminf(u_lims[b_k + 3], fmaxf(u_lims[b_k],
                                    u[unew_tmp] - f));
        }

        /*  Step the dynamics forward */
        for (i1 = 0; i1 < 7; i1++) {
            b_xnew[i1] = xnew[i1 + 7 * k];
        }

        satellite_step(b_xnew, *(float (*)[3])&unew[3 * k], *(float (*)[7])&
                       xnew[7 * (k + 1)], *(float (*)[36])&fx[36 * k], *(float (*)
                        [18])&fu[18 * k]);

        /*  Calculate the cost */
        /*  Calculates the cost contribution of a given state and control  */
        /*  Also calculates the 2nd order expansion of the cost function */
        /*  Inputs */
        /* ===================================== */
        /*  x        - [quaternion; omega] (7x1) */
        /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
        /*  terminal - int 0 or 1 */
        /* --------------------------------------------------- */
        /*  control hessian */
        /*  angular velocity hessian */
        /*  Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q)) */
        /*  Also record the sign for use in the backward pass */
        q_idx_0 = xnew[7 * k];
        xg_tmp = ((xg[0] * q_idx_0 + xg[1] * xnew[q_idx_1_tmp]) + xg[2] *
                  xnew[q_idx_2_tmp]) + xg[3] * xnew[q_idx_3_tmp];
        if (xg_tmp + 1.0F < 1.0F - xg_tmp) {
            out = xg_tmp + 1.0F;
            b_sign = 1;
        } else {
            out = 1.0F - xg_tmp;
            b_sign = -1;
        }

        c[0] = xnew[dx_tmp] - xg[4];
        c[1] = xnew[b_dx_tmp] - xg[5];
        c[2] = xnew[c_dx_tmp] - xg[6];

        /*  State cost Hessian */
        for (i = 0; i < 3; i++) {
            b_k = 6 * i + 36 * k;
            cxx[b_k] = (float)(-b_sign * iv[3 * i]) * xg_tmp;
            dx_tmp = 6 * (i + 3) + 36 * k;
            cxx[dx_tmp] = 0.0F;
            cxx[b_k + 1] = (float)(-b_sign * iv[3 * i + 1]) * xg_tmp;
            cxx[dx_tmp + 1] = 0.0F;
            cxx[b_k + 2] = (float)(-b_sign * iv[3 * i + 2]) * xg_tmp;
            cxx[dx_tmp + 2] = 0.0F;
        }

        for (i = 0; i < 6; i++) {
            b_k = 6 * i + 36 * k;
            cxx[b_k + 3] = fv[3 * i];
            cxx[b_k + 4] = fv[3 * i + 1];
            cxx[b_k + 5] = fv[3 * i + 2];
        }

        /*  State cost Jacobian */
        /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
        b_fv[3] = -xnew[q_idx_3_tmp];
        b_fv[6] = xnew[q_idx_2_tmp];
        b_fv[1] = xnew[q_idx_3_tmp];
        b_fv[7] = -xnew[q_idx_1_tmp];
        b_fv[2] = -xnew[q_idx_2_tmp];
        b_fv[5] = xnew[q_idx_1_tmp];
        for (i = 0; i < 3; i++) {
            c_sign[i] = (float)b_sign * -xnew[(i + 7 * k) + 1];
            b_k = 3 * (i + 1);
            c_sign[b_k] = (float)b_sign * (q_idx_0 * (float)iv[i] + b_fv[i]);
            c_sign[b_k + 1] = (float)b_sign * (q_idx_0 * (float)iv[i + 3] +
                b_fv[i + 3]);
            c_sign[b_k + 2] = (float)b_sign * (q_idx_0 * (float)iv[i + 6] +
                b_fv[i + 6]);
        }

        /*  Control cost Hessian & Jacobian */
        f = 0.0F;
        f1 = unew[3 * k];
        i = 3 * k + 1;
        b_dx_tmp = 3 * k + 2;
        for (c_dx_tmp = 0; c_dx_tmp < 3; c_dx_tmp++) {
            b_k = c_dx_tmp + 6 * k;
            cx[b_k] = ((c_sign[c_dx_tmp] * xg[0] + c_sign[c_dx_tmp + 3] * xg[1])
                       + c_sign[c_dx_tmp + 6] * xg[2]) + c_sign[c_dx_tmp + 9] *
                xg[3];
            cx[b_k + 3] = (fv1[c_dx_tmp] * c[0] + fv1[c_dx_tmp + 3] * c[1]) +
                fv1[c_dx_tmp + 6] * c[2];
            b_k = 3 * c_dx_tmp + 9 * k;
            q_idx_0 = fv2[3 * c_dx_tmp];
            cuu[b_k] = q_idx_0;
            dx_tmp = 3 * c_dx_tmp + 1;
            cuu[b_k + 1] = fv2[dx_tmp];
            unew_tmp = 3 * c_dx_tmp + 2;
            cuu[b_k + 2] = fv2[unew_tmp];
            cu[c_dx_tmp + 3 * k] = (fv2[c_dx_tmp] * f1 + fv2[c_dx_tmp + 3] *
                                    unew[i]) + fv2[c_dx_tmp + 6] * unew[b_dx_tmp];
            f += ((0.5F * c[0] * fv1[3 * c_dx_tmp] + 0.5F * c[1] * fv1[dx_tmp])
                  + 0.5F * c[2] * fv1[unew_tmp]) * c[c_dx_tmp];
            b_fv1[c_dx_tmp] = (0.5F * f1 * q_idx_0 + 0.5F * unew[i] * fv2[dx_tmp])
                + 0.5F * unew[b_dx_tmp] * fv2[unew_tmp];
        }

        *cost += (out + f) + ((b_fv1[0] * f1 + b_fv1[1] * unew[3 * k + 1]) +
                              b_fv1[2] * unew[3 * k + 2]);
    }

    /*  Final cost */
    /*  Calculates the cost contribution of a given state and control  */
    /*  Also calculates the 2nd order expansion of the cost function */
    /*  Inputs */
    /* ===================================== */
    /*  x        - [quaternion; omega] (7x1) */
    /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
    /*  terminal - int 0 or 1 */
    /* --------------------------------------------------- */
    /*  control hessian */
    /*  terminal angular velocity hessian  */
    /*  Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q)) */
    /*  Also record the sign for use in the backward pass */
    xg_tmp = ((xg[0] * xnew[3493] + xg[1] * xnew[3494]) + xg[2] * xnew[3495]) +
        xg[3] * xnew[3496];
    if (xg_tmp + 1.0F < 1.0F - xg_tmp) {
        out = xg_tmp + 1.0F;
        b_sign = 1;
    } else {
        out = 1.0F - xg_tmp;
        b_sign = -1;
    }

    c[0] = xnew[3497] - xg[4];
    c[1] = xnew[3498] - xg[5];
    c[2] = xnew[3499] - xg[6];

    /*  State cost Hessian */
    for (i = 0; i < 3; i++) {
        cxx[6 * i + 17964] = (float)(-b_sign * iv[3 * i]) * xg_tmp;
        b_k = 6 * (i + 3);
        cxx[b_k + 17964] = 0.0F;
        cxx[6 * i + 17965] = (float)(-b_sign * iv[3 * i + 1]) * xg_tmp;
        cxx[b_k + 17965] = 0.0F;
        cxx[6 * i + 17966] = (float)(-b_sign * iv[3 * i + 2]) * xg_tmp;
        cxx[b_k + 17966] = 0.0F;
    }

    for (i = 0; i < 6; i++) {
        cxx[6 * i + 17967] = iv1[3 * i];
        cxx[6 * i + 17968] = iv1[3 * i + 1];
        cxx[6 * i + 17969] = iv1[3 * i + 2];
    }

    /*  State cost Jacobian */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    b_fv[0] = 0.0F;
    b_fv[3] = -xnew[3496];
    b_fv[6] = xnew[3495];
    b_fv[1] = xnew[3496];
    b_fv[4] = 0.0F;
    b_fv[7] = -xnew[3494];
    b_fv[2] = -xnew[3495];
    b_fv[5] = xnew[3494];
    b_fv[8] = 0.0F;
    for (i = 0; i < 3; i++) {
        c_sign[i] = (float)b_sign * -xnew[i + 3494];
        b_k = 3 * (i + 1);
        c_sign[b_k] = (float)b_sign * (xnew[3493] * (float)iv[i] + b_fv[i]);
        c_sign[b_k + 1] = (float)b_sign * (xnew[3493] * (float)iv[i + 3] +
            b_fv[i + 3]);
        c_sign[b_k + 2] = (float)b_sign * (xnew[3493] * (float)iv[i + 6] +
            b_fv[i + 6]);
    }

    /*  Control cost Hessian & Jacobian */
    f = 0.0F;
    for (i = 0; i < 3; i++) {
        cx[i + 2994] = ((c_sign[i] * xg[0] + c_sign[i + 3] * xg[1]) + c_sign[i +
                        6] * xg[2]) + c_sign[i + 9] * xg[3];
        cx[i + 2997] = ((float)iv[i] * c[0] + (float)iv[i + 3] * c[1]) + (float)
            iv[i + 6] * c[2];
        f += ((0.5F * c[0] * (float)iv[3 * i] + 0.5F * c[1] * (float)iv[3 * i +
               1]) + 0.5F * c[2] * (float)iv[3 * i + 2]) * c[i];
    }

    *cost += out + f;
}

/*
 * Arguments    : const float x0[7]
 *                const float u0[3]
 *                float x[7]
 *                float fx[36]
 *                float fu[18]
 * Return Type  : void
 */
static void satellite_step(const float x0[7], const float u0[3], float x[7],
    float fx[36], float fu[18])
{
    float qdot_tmp[3];
    float b_qdot_tmp[9];
    int i;
    float wdot_tmp[9];
    float dxdot1_tmp;
    float f;
    float b_wdot_tmp[9];
    float c_wdot_tmp[9];
    float f1;
    float c[3];
    int i1;
    int wdot_tmp_tmp;
    float dxdot1[70];
    int b_dxdot1_tmp;
    float a[9];
    float b_fv[4];
    float b_fv1[12];
    short i2;
    static const short b_a[9] = { -200, 0, 0, 0, -200, 0, 0, 0, -200 };

    short i3;
    int c_dxdot1_tmp;
    float b_c[7];
    float c_a[3];
    static const signed char d_a[9] = { 100, 0, 0, 0, 100, 0, 0, 0, 100 };

    static const signed char b_iv[21] = { 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0,
        100, 0, 0, 0, 0, 0, 0, 0, 100 };

    float dxdot2[70];
    float x1[7];
    float fx_tmp[49];
    float b_fx_tmp[42];
    float b_fv2[49];
    static const signed char b_iv1[49] = { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };

    float b_x0[42];
    float c_fx_tmp[42];
    float fv3[21];

    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Step dynamics */
    /* Explicit midpoint step from x_k to x_{k+1} */
    /*  Calculates the continuous time state derivative and Jacobians */
    /*  kgm^2 */
    /*  Angular velocity */
    /*  Quaternion components */
    /*  Non-linear dynamics */
    qdot_tmp[0] = -x0[1];
    qdot_tmp[1] = -x0[2];
    qdot_tmp[2] = -x0[3];

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    b_qdot_tmp[0] = 0.0F;
    b_qdot_tmp[3] = -x0[3];
    b_qdot_tmp[6] = x0[2];
    b_qdot_tmp[1] = x0[3];
    b_qdot_tmp[4] = 0.0F;
    b_qdot_tmp[7] = -x0[1];
    b_qdot_tmp[2] = -x0[2];
    b_qdot_tmp[5] = x0[1];
    b_qdot_tmp[8] = 0.0F;
    for (i = 0; i < 9; i++) {
        b_qdot_tmp[i] += x0[0] * (float)iv[i];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    wdot_tmp[0] = 0.0F;
    wdot_tmp[3] = -x0[6];
    wdot_tmp[6] = x0[5];
    wdot_tmp[1] = x0[6];
    wdot_tmp[4] = 0.0F;
    wdot_tmp[7] = -x0[4];
    wdot_tmp[2] = -x0[5];
    wdot_tmp[5] = x0[4];
    wdot_tmp[8] = 0.0F;

    /*  Jacobians  */
    for (i = 0; i < 3; i++) {
        dxdot1_tmp = wdot_tmp[i + 3];
        f = wdot_tmp[i + 6];
        f1 = 0.0F;
        for (i1 = 0; i1 < 3; i1++) {
            wdot_tmp_tmp = i + 3 * i1;
            c_wdot_tmp[wdot_tmp_tmp] = (wdot_tmp[i] * fv1[3 * i1] + dxdot1_tmp *
                fv1[3 * i1 + 1]) + f * fv1[3 * i1 + 2];
            f1 += fv1[wdot_tmp_tmp] * x0[i1 + 4];
        }

        c[i] = f1;
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /* !!!!! This changes slightly with magnetorquers !!!! */
    b_wdot_tmp[0] = c_wdot_tmp[0];
    b_wdot_tmp[3] = c_wdot_tmp[3] - (-c[2]);
    b_wdot_tmp[6] = c_wdot_tmp[6] - c[1];
    b_wdot_tmp[1] = c_wdot_tmp[1] - c[2];
    b_wdot_tmp[4] = c_wdot_tmp[4];
    b_wdot_tmp[7] = c_wdot_tmp[7] - (-c[0]);
    b_wdot_tmp[2] = c_wdot_tmp[2] - (-c[1]);
    b_wdot_tmp[5] = c_wdot_tmp[5] - c[0];
    b_wdot_tmp[8] = c_wdot_tmp[8];
    dxdot1[0] = 0.0F;
    for (i = 0; i < 3; i++) {
        dxdot1_tmp = x0[i + 4];
        wdot_tmp_tmp = 7 * (i + 1);
        dxdot1[wdot_tmp_tmp] = 0.5F * -dxdot1_tmp;
        b_dxdot1_tmp = 7 * (i + 4);
        dxdot1[b_dxdot1_tmp] = 0.5F * qdot_tmp[i];
        dxdot1[i + 1] = 0.5F * dxdot1_tmp;
        i2 = b_a[i + 3];
        i3 = b_a[i + 6];
        for (i1 = 0; i1 < 3; i1++) {
            a[i + 3 * i1] = ((float)b_a[i] * b_wdot_tmp[3 * i1] + (float)i2 *
                             b_wdot_tmp[3 * i1 + 1]) + (float)i3 * b_wdot_tmp[3 *
                i1 + 2];
            c_dxdot1_tmp = i1 + 3 * i;
            dxdot1[(i1 + wdot_tmp_tmp) + 1] = 0.5F * -wdot_tmp[c_dxdot1_tmp];
            dxdot1[(i1 + b_dxdot1_tmp) + 1] = 0.5F * b_qdot_tmp[c_dxdot1_tmp];
        }
    }

    for (i = 0; i < 4; i++) {
        dxdot1[7 * i + 4] = 0.0F;
        dxdot1[7 * i + 5] = 0.0F;
        dxdot1[7 * i + 6] = 0.0F;
    }

    for (i = 0; i < 3; i++) {
        wdot_tmp_tmp = 7 * (i + 4);
        dxdot1[wdot_tmp_tmp + 4] = 0.5F * a[3 * i];
        b_dxdot1_tmp = 3 * i + 1;
        dxdot1[wdot_tmp_tmp + 5] = 0.5F * a[b_dxdot1_tmp];
        c_dxdot1_tmp = 3 * i + 2;
        dxdot1[wdot_tmp_tmp + 6] = 0.5F * a[c_dxdot1_tmp];
        for (i1 = 0; i1 < 7; i1++) {
            dxdot1[i1 + 7 * (i + 7)] = b_iv[i1 + 7 * i];
        }

        i1 = i << 2;
        b_fv1[i1] = 0.5F * qdot_tmp[i];
        b_fv1[i1 + 1] = 0.5F * b_qdot_tmp[3 * i];
        b_fv1[i1 + 2] = 0.5F * b_qdot_tmp[b_dxdot1_tmp];
        b_fv1[i1 + 3] = 0.5F * b_qdot_tmp[c_dxdot1_tmp];
    }

    for (i = 0; i < 4; i++) {
        b_fv[i] = (b_fv1[i] * x0[4] + b_fv1[i + 4] * x0[5]) + b_fv1[i + 8] * x0
            [6];
    }

    for (i = 0; i < 3; i++) {
        c[i] = u0[i] - ((c_wdot_tmp[i] * x0[4] + c_wdot_tmp[i + 3] * x0[5]) +
                        c_wdot_tmp[i + 6] * x0[6]);
    }

    for (i = 0; i < 3; i++) {
        c_a[i] = ((float)d_a[i] * c[0] + (float)d_a[i + 3] * c[1]) + (float)
            d_a[i + 6] * c[2];
    }

    b_c[0] = x0[0] + 0.015F * b_fv[0];
    b_c[1] = x0[1] + 0.015F * b_fv[1];
    b_c[2] = x0[2] + 0.015F * b_fv[2];
    b_c[3] = x0[3] + 0.015F * b_fv[3];

    /*  Calculates the continuous time state derivative and Jacobians */
    /*  kgm^2 */
    /*  Angular velocity */
    /*  Quaternion components */
    /*  Non-linear dynamics */
    b_c[4] = x0[4] + 0.015F * c_a[0];
    qdot_tmp[0] = -b_c[1];
    b_c[5] = x0[5] + 0.015F * c_a[1];
    qdot_tmp[1] = -b_c[2];
    b_c[6] = x0[6] + 0.015F * c_a[2];
    qdot_tmp[2] = -b_c[3];

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    b_qdot_tmp[0] = 0.0F;
    b_qdot_tmp[3] = -b_c[3];
    b_qdot_tmp[6] = b_c[2];
    b_qdot_tmp[1] = b_c[3];
    b_qdot_tmp[4] = 0.0F;
    b_qdot_tmp[7] = -b_c[1];
    b_qdot_tmp[2] = -b_c[2];
    b_qdot_tmp[5] = b_c[1];
    b_qdot_tmp[8] = 0.0F;
    for (i = 0; i < 9; i++) {
        b_qdot_tmp[i] += b_c[0] * (float)iv[i];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    wdot_tmp[0] = 0.0F;
    wdot_tmp[3] = -b_c[6];
    wdot_tmp[6] = b_c[5];
    wdot_tmp[1] = b_c[6];
    wdot_tmp[4] = 0.0F;
    wdot_tmp[7] = -b_c[4];
    wdot_tmp[2] = -b_c[5];
    wdot_tmp[5] = b_c[4];
    wdot_tmp[8] = 0.0F;

    /*  Jacobians  */
    for (i = 0; i < 3; i++) {
        dxdot1_tmp = wdot_tmp[i + 3];
        f = wdot_tmp[i + 6];
        f1 = 0.0F;
        for (i1 = 0; i1 < 3; i1++) {
            wdot_tmp_tmp = i + 3 * i1;
            c_wdot_tmp[wdot_tmp_tmp] = (wdot_tmp[i] * fv1[3 * i1] + dxdot1_tmp *
                fv1[3 * i1 + 1]) + f * fv1[3 * i1 + 2];
            f1 += fv1[wdot_tmp_tmp] * b_c[i1 + 4];
        }

        c[i] = f1;
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /* !!!!! This changes slightly with magnetorquers !!!! */
    b_wdot_tmp[0] = c_wdot_tmp[0];
    b_wdot_tmp[3] = c_wdot_tmp[3] - (-c[2]);
    b_wdot_tmp[6] = c_wdot_tmp[6] - c[1];
    b_wdot_tmp[1] = c_wdot_tmp[1] - c[2];
    b_wdot_tmp[4] = c_wdot_tmp[4];
    b_wdot_tmp[7] = c_wdot_tmp[7] - (-c[0]);
    b_wdot_tmp[2] = c_wdot_tmp[2] - (-c[1]);
    b_wdot_tmp[5] = c_wdot_tmp[5] - c[0];
    b_wdot_tmp[8] = c_wdot_tmp[8];
    dxdot2[0] = 0.0F;
    for (i = 0; i < 3; i++) {
        dxdot1_tmp = b_c[i + 4];
        wdot_tmp_tmp = 7 * (i + 1);
        dxdot2[wdot_tmp_tmp] = 0.5F * -dxdot1_tmp;
        b_dxdot1_tmp = 7 * (i + 4);
        dxdot2[b_dxdot1_tmp] = 0.5F * qdot_tmp[i];
        dxdot2[i + 1] = 0.5F * dxdot1_tmp;
        i2 = b_a[i + 3];
        i3 = b_a[i + 6];
        for (i1 = 0; i1 < 3; i1++) {
            a[i + 3 * i1] = ((float)b_a[i] * b_wdot_tmp[3 * i1] + (float)i2 *
                             b_wdot_tmp[3 * i1 + 1]) + (float)i3 * b_wdot_tmp[3 *
                i1 + 2];
            c_dxdot1_tmp = i1 + 3 * i;
            dxdot2[(i1 + wdot_tmp_tmp) + 1] = 0.5F * -wdot_tmp[c_dxdot1_tmp];
            dxdot2[(i1 + b_dxdot1_tmp) + 1] = 0.5F * b_qdot_tmp[c_dxdot1_tmp];
        }
    }

    for (i = 0; i < 4; i++) {
        dxdot2[7 * i + 4] = 0.0F;
        dxdot2[7 * i + 5] = 0.0F;
        dxdot2[7 * i + 6] = 0.0F;
    }

    for (i = 0; i < 3; i++) {
        wdot_tmp_tmp = 7 * (i + 4);
        dxdot2[wdot_tmp_tmp + 4] = 0.5F * a[3 * i];
        b_dxdot1_tmp = 3 * i + 1;
        dxdot2[wdot_tmp_tmp + 5] = 0.5F * a[b_dxdot1_tmp];
        c_dxdot1_tmp = 3 * i + 2;
        dxdot2[wdot_tmp_tmp + 6] = 0.5F * a[c_dxdot1_tmp];
        for (i1 = 0; i1 < 7; i1++) {
            dxdot2[i1 + 7 * (i + 7)] = b_iv[i1 + 7 * i];
        }

        i1 = i << 2;
        b_fv1[i1] = 0.5F * qdot_tmp[i];
        b_fv1[i1 + 1] = 0.5F * b_qdot_tmp[3 * i];
        b_fv1[i1 + 2] = 0.5F * b_qdot_tmp[b_dxdot1_tmp];
        b_fv1[i1 + 3] = 0.5F * b_qdot_tmp[c_dxdot1_tmp];
    }

    for (i = 0; i < 4; i++) {
        b_fv[i] = (b_fv1[i] * b_c[4] + b_fv1[i + 4] * b_c[5]) + b_fv1[i + 8] *
            b_c[6];
    }

    for (i = 0; i < 3; i++) {
        c[i] = u0[i] - ((c_wdot_tmp[i] * b_c[4] + c_wdot_tmp[i + 3] * b_c[5]) +
                        c_wdot_tmp[i + 6] * b_c[6]);
    }

    for (i = 0; i < 3; i++) {
        c_a[i] = ((float)d_a[i] * c[0] + (float)d_a[i + 3] * c[1]) + (float)
            d_a[i + 6] * c[2];
    }

    x1[0] = x0[0] + 0.03F * b_fv[0];
    x1[1] = x0[1] + 0.03F * b_fv[1];
    x1[2] = x0[2] + 0.03F * b_fv[2];
    x1[3] = x0[3] + 0.03F * b_fv[3];
    x1[4] = x0[4] + 0.03F * c_a[0];
    x1[5] = x0[5] + 0.03F * c_a[1];
    x1[6] = x0[6] + 0.03F * c_a[2];

    /*  Normalize the quaternion */
    dxdot1_tmp = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3]
                       * x1[3]);
    x[0] = x1[0] / dxdot1_tmp;
    x[1] = x1[1] / dxdot1_tmp;
    x[2] = x1[2] / dxdot1_tmp;
    x[3] = x1[3] / dxdot1_tmp;
    x[4] = x1[4];
    x[5] = x1[5];
    x[6] = x1[6];

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    for (i = 0; i < 7; i++) {
        for (i1 = 0; i1 < 7; i1++) {
            wdot_tmp_tmp = i1 + 7 * i;
            fx_tmp[wdot_tmp_tmp] = 0.00045F * dxdot2[wdot_tmp_tmp];
        }
    }

    b_qdot_tmp[0] = 0.0F;
    b_qdot_tmp[3] = -x1[3];
    b_qdot_tmp[6] = x1[2];
    b_qdot_tmp[1] = x1[3];
    b_qdot_tmp[4] = 0.0F;
    b_qdot_tmp[7] = -x1[1];
    b_qdot_tmp[2] = -x1[2];
    b_qdot_tmp[5] = x1[1];
    b_qdot_tmp[8] = 0.0F;
    for (i = 0; i < 3; i++) {
        b_fx_tmp[i] = -x1[i + 1];
        b_fx_tmp[i + 3] = 0.0F;
        wdot_tmp_tmp = 6 * (i + 1);
        b_fx_tmp[wdot_tmp_tmp] = x1[0] * (float)iv[i] + b_qdot_tmp[i];
        b_fx_tmp[wdot_tmp_tmp + 3] = 0.0F;
        b_fx_tmp[wdot_tmp_tmp + 1] = x1[0] * (float)iv[i + 3] + b_qdot_tmp[i + 3];
        b_fx_tmp[wdot_tmp_tmp + 4] = 0.0F;
        b_fx_tmp[wdot_tmp_tmp + 2] = x1[0] * (float)iv[i + 6] + b_qdot_tmp[i + 6];
        b_fx_tmp[wdot_tmp_tmp + 5] = 0.0F;
        for (i1 = 0; i1 < 6; i1++) {
            b_fx_tmp[i1 + 6 * (i + 4)] = iv1[i + 3 * i1];
        }
    }

    for (i = 0; i < 7; i++) {
        for (i1 = 0; i1 < 7; i1++) {
            dxdot1_tmp = 0.0F;
            for (wdot_tmp_tmp = 0; wdot_tmp_tmp < 7; wdot_tmp_tmp++) {
                dxdot1_tmp += fx_tmp[i + 7 * wdot_tmp_tmp] * dxdot1[wdot_tmp_tmp
                    + 7 * i1];
            }

            wdot_tmp_tmp = i + 7 * i1;
            b_fv2[wdot_tmp_tmp] = ((float)b_iv1[wdot_tmp_tmp] + 0.03F *
                                   dxdot2[wdot_tmp_tmp]) + dxdot1_tmp;
        }
    }

    b_qdot_tmp[0] = 0.0F;
    b_qdot_tmp[3] = -x0[3];
    b_qdot_tmp[6] = x0[2];
    b_qdot_tmp[1] = x0[3];
    b_qdot_tmp[4] = 0.0F;
    b_qdot_tmp[7] = -x0[1];
    b_qdot_tmp[2] = -x0[2];
    b_qdot_tmp[5] = x0[1];
    b_qdot_tmp[8] = 0.0F;
    for (i = 0; i < 6; i++) {
        for (i1 = 0; i1 < 7; i1++) {
            dxdot1_tmp = 0.0F;
            for (wdot_tmp_tmp = 0; wdot_tmp_tmp < 7; wdot_tmp_tmp++) {
                dxdot1_tmp += b_fx_tmp[i + 6 * wdot_tmp_tmp] *
                    b_fv2[wdot_tmp_tmp + 7 * i1];
            }

            c_fx_tmp[i + 6 * i1] = dxdot1_tmp;
        }
    }

    for (i = 0; i < 3; i++) {
        b_x0[7 * i] = -x0[i + 1];
        wdot_tmp_tmp = 7 * (i + 3);
        b_x0[wdot_tmp_tmp] = 0.0F;
        b_x0[7 * i + 1] = x0[0] * (float)iv[3 * i] + b_qdot_tmp[3 * i];
        b_x0[wdot_tmp_tmp + 1] = 0.0F;
        b_dxdot1_tmp = 3 * i + 1;
        b_x0[7 * i + 2] = x0[0] * (float)iv[b_dxdot1_tmp] +
            b_qdot_tmp[b_dxdot1_tmp];
        b_x0[wdot_tmp_tmp + 2] = 0.0F;
        b_dxdot1_tmp = 3 * i + 2;
        b_x0[7 * i + 3] = x0[0] * (float)iv[b_dxdot1_tmp] +
            b_qdot_tmp[b_dxdot1_tmp];
        b_x0[wdot_tmp_tmp + 3] = 0.0F;
    }

    for (i = 0; i < 6; i++) {
        b_x0[7 * i + 4] = iv1[3 * i];
        b_x0[7 * i + 5] = iv1[3 * i + 1];
        b_x0[7 * i + 6] = iv1[3 * i + 2];
    }

    for (i = 0; i < 6; i++) {
        for (i1 = 0; i1 < 6; i1++) {
            dxdot1_tmp = 0.0F;
            for (wdot_tmp_tmp = 0; wdot_tmp_tmp < 7; wdot_tmp_tmp++) {
                dxdot1_tmp += c_fx_tmp[i + 6 * wdot_tmp_tmp] * b_x0[wdot_tmp_tmp
                    + 7 * i1];
            }

            fx[i + 6 * i1] = dxdot1_tmp;
        }
    }

    for (i = 0; i < 7; i++) {
        for (i1 = 0; i1 < 3; i1++) {
            dxdot1_tmp = 0.0F;
            for (wdot_tmp_tmp = 0; wdot_tmp_tmp < 7; wdot_tmp_tmp++) {
                dxdot1_tmp += fx_tmp[i + 7 * wdot_tmp_tmp] * dxdot1[wdot_tmp_tmp
                    + 7 * (i1 + 7)];
            }

            fv3[i + 7 * i1] = 0.03F * dxdot2[i + 7 * (i1 + 7)] + dxdot1_tmp;
        }
    }

    for (i = 0; i < 6; i++) {
        for (i1 = 0; i1 < 3; i1++) {
            dxdot1_tmp = 0.0F;
            for (wdot_tmp_tmp = 0; wdot_tmp_tmp < 7; wdot_tmp_tmp++) {
                dxdot1_tmp += b_fx_tmp[i + 6 * wdot_tmp_tmp] * fv3[wdot_tmp_tmp
                    + 7 * i1];
            }

            fu[i + 6 * i1] = dxdot1_tmp;
        }
    }
}

/*
 * Arguments    : const float x0[3500]
 *                const float xg[7]
 *                const float u0[1497]
 *                const float u_lims[6]
 *                float x[3500]
 *                float u[1497]
 *                float K[8982]
 *                bool *result
 * Return Type  : void
 */
void milqr(const float x0[3500], const float xg[7], const float u0[1497], const
           float u_lims[6], float x[3500], float u[1497], float K[8982],
           bool *result)
{
    int k;
    float Alphas[11];
    float lambda;
    float dlambda;
    static float x_n[3500];
    static float u_n[1497];
    static float fx_n[17964];
    static float fu_n[8982];
    static float cx_n[3000];
    static float cu_n[1497];
    static float cxx_n[18000];
    static float cuu_n[4491];
    float cost_n;
    static float l[1497];
    float dV[2];
    static float b_fv[8982];
    static float fx[17964];
    static float fu[8982];
    static float cx[3000];
    float cu[1497];
    static float cxx[18000];
    static float cuu[4491];
    float cost;
    int iter;
    int b_iter;
    bool exitg1;
    bool backPassDone;
    int exitg2;
    bool diverge;
    float fcnOutput[1497];
    float b_x;
    float maxval[499];
    float f;
    bool fwdPassDone;
    bool exitg3;
    float expectedChange;
    float z;
    float dcost;

    /*  Solves finite horizon optimal control problem using a */
    /*  multiplicative iterative linear quadratic regulator */
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
    /*  gradient exit criterion */
    /*  lambda criterion for gradient exit */
    /*  minimum accepted cost reduction ratio */
    /*  maximum regularization parameter */
    /*  set lambda = 0 below this value */
    /*  amount to scale dlambda by */
    /*  CONSTANTS */
    for (k = 0; k < 11; k++) {
        Alphas[k] = powf(10.0F, -0.3F * (float)k);
    }

    /*  line search param */
    lambda = 1.0F;
    dlambda = 1.0F;

    /*  error state size (3 param. error representation for attitude) */
    /*  Init matrices for update (otherwise MATLAB coder throws an error) */
    memset(&x_n[0], 0, 3500U * sizeof(float));
    memset(&u_n[0], 0, 1497U * sizeof(float));
    memset(&fx_n[0], 0, 17964U * sizeof(float));
    memset(&fu_n[0], 0, 8982U * sizeof(float));
    memset(&cx_n[0], 0, 3000U * sizeof(float));
    memset(&cu_n[0], 0, 1497U * sizeof(float));
    memset(&cxx_n[0], 0, 18000U * sizeof(float));
    memset(&cuu_n[0], 0, 4491U * sizeof(float));
    cost_n = 0.0F;

    /*  Initial Forward rollout */
    memset(&l[0], 0, 1497U * sizeof(float));
    dV[0] = 0.0F;
    dV[1] = 0.0F;
    memset(&K[0], 0, 8982U * sizeof(float));
    memset(&b_fv[0], 0, 8982U * sizeof(float));
    forwardRollout(x0, xg, u0, b_fv, u_lims, x, u, fx, fu, cx, cu, cxx, cuu,
                   &cost);

    /*  Convergence check params */
    /*  Expected cost change */
    /*  Ratio of cost change to expected cost change */
    *result = false;
    iter = 1;
    b_iter = 0;
    exitg1 = false;
    while ((!exitg1) && (b_iter < 500)) {
        iter = b_iter + 1;

        /*  Backward Pass */
        /* ======================================= */
        backPassDone = false;
        do {
            exitg2 = 0;
            if (!backPassDone) {
                backwardPass(fx, fu, cx, cu, cxx, cuu, lambda, u_lims, u, l, K,
                             dV, &diverge);
                if (diverge) {
                    /*  Increase regularization parameter (lambda) */
                    dlambda = fmaxf(1.6F * dlambda, 1.6F);
                    lambda = fmaxf(lambda * dlambda, 1.0E-6F);
                    if (lambda > 1.0E+10F) {
                        exitg2 = 1;
                    }
                } else {
                    backPassDone = true;
                }
            } else {
                exitg2 = 1;
            }
        } while (exitg2 == 0);

        /*  Check gradient of control, defined as l/u */
        /*  Terminate if sufficiently small (success) */
        for (k = 0; k < 1497; k++) {
            fcnOutput[k] = fabsf(l[k]) / (fabsf(u[k]) + 1.0F);
        }

        for (k = 0; k < 499; k++) {
            b_x = fcnOutput[3 * k];
            f = fcnOutput[3 * k + 1];
            if (b_x < f) {
                b_x = f;
            }

            f = fcnOutput[3 * k + 2];
            if (b_x < f) {
                b_x = f;
            }

            maxval[k] = b_x;
        }

        b_x = maxval[0];
        for (k = 0; k < 498; k++) {
            b_x += maxval[k + 1];
        }

        /*  Avg of max grad at each time step */
        if ((b_x / 499.0F < 0.0001F) && (lambda < 1.0E-5F)) {
            *result = true;
            exitg1 = true;
        } else {
            /*  Forward Line-Search */
            /* =========================================== */
            fwdPassDone = false;
            if (backPassDone) {
                k = 0;
                exitg3 = false;
                while ((!exitg3) && (k < 11)) {
                    b_forwardRollout(x, xg, u, l, K, Alphas[k], u_lims, x_n, u_n,
                                     fx_n, fu_n, cx_n, cu_n, cxx_n, cuu_n,
                                     &cost_n);
                    expectedChange = -Alphas[k] * (dV[0] + Alphas[k] * dV[1]);
                    if (expectedChange > 0.0F) {
                        z = (cost - cost_n) / expectedChange;
                    } else {
                        b_x = cost - cost_n;
                        if (b_x < 0.0F) {
                            b_x = -1.0F;
                        } else {
                            if (b_x > 0.0F) {
                                b_x = 1.0F;
                            }
                        }

                        z = b_x;
                    }

                    if (z > 0.0F) {
                        fwdPassDone = true;
                        exitg3 = true;
                    } else {
                        k++;
                    }
                }
            }

            /*  Parameter Updates */
            /* ============================================= */
            if (fwdPassDone) {
                /*  Decrease Lambda */
                dlambda = fminf(dlambda / 1.6F, 0.625F);
                lambda = lambda * dlambda * (float)(lambda > 1.0E-6F);

                /*  set = 0 if lambda too small */
                dcost = cost - cost_n;

                /*  Update trajectory and controls */
                memcpy(&x[0], &x_n[0], 3500U * sizeof(float));
                memcpy(&u[0], &u_n[0], 1497U * sizeof(float));
                memcpy(&fx[0], &fx_n[0], 17964U * sizeof(float));
                memcpy(&fu[0], &fu_n[0], 8982U * sizeof(float));
                memcpy(&cx[0], &cx_n[0], 3000U * sizeof(float));
                memcpy(&cu[0], &cu_n[0], 1497U * sizeof(float));
                memcpy(&cxx[0], &cxx_n[0], 18000U * sizeof(float));
                memcpy(&cuu[0], &cuu_n[0], 4491U * sizeof(float));
                cost = cost_n;

                /*  Terminate if cost change sufficiently small */
                if (dcost < 1.0E-7F) {
                    *result = true;
                    exitg1 = true;
                } else {
                    b_iter++;
                }
            } else {
                /*  No cost reduction (based on z-value) */
                /*  Increase lambda */
                dlambda = fmaxf(1.6F * dlambda, 1.6F);
                lambda = fmaxf(lambda * dlambda, 1.0E-6F);
                if (lambda > 1.0E+10F) {
                    /*  Lambda too large - solver diverged */
                    exitg1 = true;
                } else {
                    b_iter++;
                }
            }
        }
    }

    if (iter == 500) {
        /*  Ddin't converge completely */
        *result = false;
    }
}

/*
 * File trailer for milqr.c
 *
 * [EOF]
 */
