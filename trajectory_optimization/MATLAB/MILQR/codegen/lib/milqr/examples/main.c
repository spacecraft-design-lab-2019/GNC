/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 16-Apr-2020 13:51:41
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
#include "milqr.h"
#include "main.h"

/* Function Declarations */
static void argInit_3x1000_real32_T(float result[3000]);
static void argInit_3x2_real32_T(float result[6]);
static void argInit_3x999_real32_T(float result[2997]);
static void argInit_7x1000_real32_T(float result[7000]);
static void argInit_7x1_real32_T(float result[7]);
static float argInit_real32_T(void);
static void main_milqr(void);

/* Function Definitions */

/*
 * Arguments    : float result[3000]
 * Return Type  : void
 */
static void argInit_3x1000_real32_T(float result[3000])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 1000; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real32_T();
    }
  }
}

/*
 * Arguments    : float result[6]
 * Return Type  : void
 */
static void argInit_3x2_real32_T(float result[6])
{
  int idx0;
  float result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real32_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 3] = result_tmp;
  }
}

/*
 * Arguments    : float result[2997]
 * Return Type  : void
 */
static void argInit_3x999_real32_T(float result[2997])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 999; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real32_T();
    }
  }
}

/*
 * Arguments    : float result[7000]
 * Return Type  : void
 */
static void argInit_7x1000_real32_T(float result[7000])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 7; idx0++) {
    for (idx1 = 0; idx1 < 1000; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 7 * idx1] = argInit_real32_T();
    }
  }
}

/*
 * Arguments    : float result[7]
 * Return Type  : void
 */
static void argInit_7x1_real32_T(float result[7])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 7; idx0++) {
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
static void main_milqr(void)
{
  static float fv14[7000];
  float fv15[7];
  static float fv16[2997];
  float fv17[6];
  static float fv18[3000];
  static float x[7000];
  static float u[2997];
  static float K[17982];
  bool result;

  /* Initialize function 'milqr' input arguments. */
  /* Initialize function input argument 'x0'. */
  /* Initialize function input argument 'xg'. */
  /* Initialize function input argument 'u0'. */
  /* Initialize function input argument 'u_lims'. */
  /* Initialize function input argument 'B_ECI'. */
  /* Call the entry-point 'milqr'. */
  argInit_7x1000_real32_T(fv14);
  argInit_7x1_real32_T(fv15);
  argInit_3x999_real32_T(fv16);
  argInit_3x2_real32_T(fv17);
  argInit_3x1000_real32_T(fv18);
  milqr(fv14, fv15, fv16, fv17, fv18, x, u, K, &result);
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
  milqr_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_milqr();

  /* Terminate the application.
     You do not need to do this more than one time. */
  milqr_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
