/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 20-Apr-2020 22:26:32
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
#include "milqr_efficient.h"
#include "main.h"
#include "get_magnetic_field_series.h"
#include <stdio.h>
#include "rtwtypes.h"

/* Function Declarations */
static void main_milqr_efficient(bool * retval);

static void main_milqr_efficient(bool * retval)
{
  static int N = 100;
  static int Nx = 7;
  static int Nu = 3;
  static double dt = 1;
  double t[N];
  for (int i = 0;i < N; i++) {
    t[i] = dt*i + 58947.77014;
  }

  double x0[] = {.7071,0,.7071,0,0,0,0};
  double xg[] = {0,0,0,1,0,0,0};

  double r0[] = {5719.8, -755.6, 4121.2, -4.4390, -1.0905, 5.945};

  double B_eci_vec[3*(N-1)];
  double X[6*N];

  get_magnetic_field_series(r0, t, B_eci_vec, X);

  for (int i = 0; i < 3*(N-1); i++) {
    B_eci_vec[i] *= 1e-9;
  }

  printf("The first values of mag_field are: %f \n", B_eci_vec[0]);

  double X0[7*N];
  for (int i = 0; i < N; i++) {
    for (int i0 = 0; i0 < 7; i0++) {
      X0[7*i+i0] = x0[i0];
    }
  }

  printf("The first value of X0 is: %f \n", X0[0]);

  double u_lims[6] = {-100.0, -100.0, -100.0, 100.0, 100.0, 100.0};
  double x[7*N];
  double u[3*(N-1)];
  double u0[3*(N-1)];
  double K[18*(N-1)];

  bool result;

  //printf("%d\n", *result);
  printf("Stepping into MILQR\n");

  printf("INPUTS:\n");
  printf("x0 = ");
  for (int i = 0;i < 7;i++) {
    printf("%f ",x0[i]);
  }
  printf("\n");

  printf("xg = ");
  for (int i = 0;i < 7;i++) {
    printf("%f ",xg[i]);
  }
  printf("\n");

  printf("u0 = ");
  for (int i = 0;i < 3;i++) {
    printf("%f ",u0[i]);
  }
  printf("\n");

  printf("u_lims = ");
  for (int i = 0;i < 6;i++) {
    printf("%f ",u_lims[i]);
  }
  printf("\n");

  printf("B_eci_vec = ");
  for (int i = 0;i < 3;i++) {
    printf("%f ",B_eci_vec[i]);
  }
  printf("\n");  


  milqr_efficient(x0,xg,u0,u_lims,dt,B_eci_vec,x,u,K,&result);
  // printf("%d\n", result);

  printf("The last value of x is: %f %f %f %f %f %f %f \n", 
    x[7*N-6],
    x[7*N-5],
    x[7*N-4],
    x[7*N-3],
    x[7*N-2],
    x[7*N-1],
    x[7*N]);

  *retval = result;

}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  bool retval = 0; 
  main_milqr_efficient(&retval);

  printf("did the algorithm converge?:\n");
  printf(retval ? "true" : "false");

  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
