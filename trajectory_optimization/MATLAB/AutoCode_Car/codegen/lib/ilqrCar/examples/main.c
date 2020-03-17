/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 01-Mar-2020 20:37:52
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
#include "main.h"
#include "ilqrCar.h"

/* Function Declarations */
static void argInit_2x2_real_T(double result[4]);
static void argInit_2x500_real_T(double result[1000]);
static void argInit_4x1_real_T(double result[4]);
static void argInit_4x501_real_T(double result[2004]);
static double argInit_real_T(void);
static void main_ilqrCar(void);

/* Function Definitions */

/*
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_2x2_real_T(double result[4])
{
  int idx0;
  double result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 2] = result_tmp;
  }
}

/*
 * Arguments    : double result[1000]
 * Return Type  : void
 */
static void argInit_2x500_real_T(double result[1000])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    for (idx1 = 0; idx1 < 500; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + (idx1 << 1)] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_4x1_real_T(double result[4])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[2004]
 * Return Type  : void
 */
static void argInit_4x501_real_T(double result[2004])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    for (idx1 = 0; idx1 < 501; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + (idx1 << 2)] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_ilqrCar(void)
{
  static double b_dv[2004];
  double dv1[4];
  static double dv2[1000];
  double dv3[4];
  static double x[2004];
  static double u[1000];
  static double K[4000];
  boolean_T result;

  /* Initialize function 'ilqrCar' input arguments. */
  /* Initialize function input argument 'x0'. */
  /* Initialize function input argument 'xg'. */
  /* Initialize function input argument 'u0'. */
  /* Initialize function input argument 'u_lims'. */
  /* Call the entry-point 'ilqrCar'. */
  argInit_4x501_real_T(b_dv);
  argInit_4x1_real_T(dv1);
  argInit_2x500_real_T(dv2);
  argInit_2x2_real_T(dv3);
  ilqrCar(b_dv, dv1, dv2, dv3, x, u, K, &result);
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

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_ilqrCar();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
