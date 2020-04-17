/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: get_magnetic_field.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 11-Apr-2020 23:47:25
 */

#ifndef GET_MAGNETIC_FIELD_H
#define GET_MAGNETIC_FIELD_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "get_magnetic_field_types.h"

/* Function Declarations */
void get_magnetic_field(float lat, float lon, float alt, float year, float order,
  float B_NED[3]);
void get_magnetic_field_initialize(void);
void get_magnetic_field_terminate(void);

#endif

/*
 * File trailer for get_magnetic_field.h
 *
 * [EOF]
 */
