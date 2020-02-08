/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 07-Feb-2020 14:37:48
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
#include "iLQRv1.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void argInit_1x399_real_T(real_T result[399]);
static void argInit_2x1_real_T(real_T result[2]);
static void argInit_2x2_real_T(real_T result[4]);
static real_T argInit_real_T(void);
static void main_iLQRv1(void);

/* Function Definitions */

/*
 * Arguments    : real_T result[399]
 * Return Type  : void
 */
static void argInit_1x399_real_T(real_T result[399])
{
  int32_T idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 399; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

/*
 * Arguments    : real_T result[2]
 * Return Type  : void
 */
static void argInit_2x1_real_T(real_T result[2])
{
  real_T result_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result[0] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp;
}

/*
 * Arguments    : real_T result[4]
 * Return Type  : void
 */
static void argInit_2x2_real_T(real_T result[4])
{
  real_T result_tmp_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp_tmp = argInit_real_T();
  result[0] = result_tmp_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[2] = result_tmp_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[3] = argInit_real_T();
}

/*
 * Arguments    : void
 * Return Type  : real_T
 */
static real_T argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_iLQRv1(void)
{
  real_T x0_tmp[2];
  real_T dv[399];
  real_T dv1[4];
  real_T dv2[4];
  real_T xtraj[800];
  real_T utraj[399];
  real_T K[798];

  /* Initialize function 'iLQRv1' input arguments. */
  /* Initialize function input argument 'x0'. */
  argInit_2x1_real_T(x0_tmp);

  /* Initialize function input argument 'xg'. */
  /* Initialize function input argument 'utraj0'. */
  /* Initialize function input argument 'Q'. */
  /* Initialize function input argument 'Qf'. */
  /* Call the entry-point 'iLQRv1'. */
  argInit_1x399_real_T(dv);
  argInit_2x2_real_T(dv1);
  argInit_2x2_real_T(dv2);
  iLQRv1(x0_tmp, x0_tmp, dv, dv1, argInit_real_T(), dv2, argInit_real_T(),
         argInit_real_T(), xtraj, utraj, K);
}

/*
 * Arguments    : int32_T argc
 *                const char * const argv[]
 * Return Type  : int32_T
 */
int32_T main(int32_T argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_iLQRv1();

  /* Terminate the application.
     You do not need to do this more than one time. */
  iLQRv1_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
