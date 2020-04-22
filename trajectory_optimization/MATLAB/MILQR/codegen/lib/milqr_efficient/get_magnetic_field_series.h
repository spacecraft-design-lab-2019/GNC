/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: get_magnetic_field_series.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 12:32:15
 */

#ifndef GET_MAGNETIC_FIELD_SERIES_H
#define GET_MAGNETIC_FIELD_SERIES_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "get_magnetic_field_series_types.h"

/* Function Declarations */
void get_magnetic_field_series(const float x0[6], const float t[100], float
  B_eci_vec[297], float X[600]);
void get_magnetic_field_series_initialize(void);
void get_magnetic_field_series_terminate(void);

#endif

/*
 * File trailer for get_magnetic_field_series.h
 *
 * [EOF]
 */
