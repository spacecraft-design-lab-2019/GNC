/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 12:32:15
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "get_magnetic_field_series.h"
#include "main.h"

/* Function Declarations */
static void argInit_1x100_real32_T(float result[100]);
static void argInit_6x1_real32_T(float result[6]);
static float argInit_real32_T(void);
static void main_get_magnetic_field_series(void);

/* Function Definitions */

/*
 * Arguments    : float result[100]
 * Return Type  : void
 */
static void argInit_1x100_real32_T(float result[100])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 100; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[6]
 * Return Type  : void
 */
static void argInit_6x1_real32_T(float result[6])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 6; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real32_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : float
 */
static float argInit_real32_T(void)
{
  return 0.0F;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_get_magnetic_field_series(void)
{
  float fv1[6];
  float fv2[100];
  float B_eci_vec[297];
  float X[600];

  /* Initialize function 'get_magnetic_field_series' input arguments. */
  /* Initialize function input argument 'x0'. */
  /* Initialize function input argument 't'. */
  /* Call the entry-point 'get_magnetic_field_series'. */
  argInit_6x1_real32_T(fv1);
  argInit_1x100_real32_T(fv2);
  get_magnetic_field_series(fv1, fv2, B_eci_vec, X);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  get_magnetic_field_series_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_get_magnetic_field_series();

  /* Terminate the application.
     You do not need to do this more than one time. */
  get_magnetic_field_series_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
