/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 20:32:25
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
#include "iLQRsimple.h"
#include "rt_nonfinite.h"


int main() {

  /* Create inputs to iLQR */
  /* iLQRsimple(const double x0[2], const double xg[2], const double
  utraj0[399], const double Q[4], double R, const double Qf[4], double dt,
  double tol, double xtraj[800], double utraj[399], double K[798]) */
  const int Nx = 2;
  const int Nu = 1;

  const double x0[2] = {0, 0};
  const double xg[2] = {3.14159, 0};
  const double utraj0[399] = { 0 };
  const double Q[4] = {0.01, 0, 0.01, 0};
  const double R = 0.3;
  const double Qf[4] = {100, 0, 100, 0};
  const double dt = 0.01;
  const double tol = 0.001;
  double xtraj[800] = { 0 };
  double utraj[399] = { 0 };
  double K[798] = { 0 };

  // Assign x0 to xtraj
  xtraj[0] = x0[0];
  xtraj[1] = x0[1];


  /* The initialize function is being called automatically from your iLQRsimple function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  // main_iLQRsimple();

  iLQRsimple(x0, xg, utraj0, Q, R, Qf, dt, tol, xtraj, utraj, K);

  /* Terminate the application.
     You do not need to do this more than one time. */
  iLQRsimple_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
