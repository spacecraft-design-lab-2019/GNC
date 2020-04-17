/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 11-Apr-2020 23:47:25
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
#include "get_magnetic_field.h"
#include "main.h"

/* Function Declarations */
static float argInit_real32_T(void);
static void main_get_magnetic_field(void);

/* Function Definitions */

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
static void main_get_magnetic_field(void)
{
  float B_NED[3];

  /* Initialize function 'get_magnetic_field' input arguments. */
  /* Call the entry-point 'get_magnetic_field'. */
  get_magnetic_field(argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
                     argInit_real32_T(), argInit_real32_T(), B_NED);
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
  get_magnetic_field_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_get_magnetic_field();

  /* Terminate the application.
     You do not need to do this more than one time. */
  get_magnetic_field_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
