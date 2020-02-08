/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 06-Feb-2020 15:05:27
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
static void argInit_1xd400_real_T(double result_data[], int result_size[2]);
static void argInit_2x1_real_T(double result[2]);
static void argInit_2x2_real_T(double result[4]);
static double argInit_real_T(void);
static void main_iLQRv1(void);

/* Function Definitions */

/*
 * Arguments    : double result_data[]
 *                int result_size[2]
 * Return Type  : void
 */
static void argInit_1xd400_real_T(double result_data[], int result_size[2])
{
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result_size[0] = 1;
  result_size[1] = 2;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 2; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_data[idx1] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[2]
 * Return Type  : void
 */
static void argInit_2x1_real_T(double result[2])
{
  double result_tmp;

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
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_2x2_real_T(double result[4])
{
  double result_tmp_tmp;

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
static void main_iLQRv1(void)
{
  double x0_tmp[2];
  double utraj0_data[400];
  int utraj0_size[2];
  double dv[4];
  double dv1[4];
  double xtraj_data[802];
  int xtraj_size[2];
  double utraj_data[400];
  int utraj_size[2];
  double K_data[800];
  int K_size[3];

  /* Initialize function 'iLQRv1' input arguments. */
  /* Initialize function input argument 'x0'. */
  argInit_2x1_real_T(x0_tmp);

  /* Initialize function input argument 'xg'. */
  /* Initialize function input argument 'utraj0'. */
  argInit_1xd400_real_T(utraj0_data, utraj0_size);

  /* Initialize function input argument 'Q'. */
  /* Initialize function input argument 'Qf'. */
  /* Call the entry-point 'iLQRv1'. */
  argInit_2x2_real_T(dv);
  argInit_2x2_real_T(dv1);
  iLQRv1(x0_tmp, x0_tmp, utraj0_data, utraj0_size, dv, argInit_real_T(), dv1,
         argInit_real_T(), argInit_real_T(), xtraj_data, xtraj_size, utraj_data,
         utraj_size, K_data, K_size);
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
