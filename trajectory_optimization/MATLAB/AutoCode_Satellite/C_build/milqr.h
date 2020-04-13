/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 12-Apr-2020 23:15:53
 */

#ifndef MILQR_H
#define MILQR_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>

/* Function Declarations */
void milqr(float x0[5000][7], const float xg[7], float u0[4999][3], float
           u_lims[2][3], float B_ECI[5000][3], float Js[6][3], const float
           options[10], float x[5000][7], float u[4999][3], float K[4999][6][3],
           bool *result);

#endif

