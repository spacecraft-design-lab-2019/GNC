/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 19-Mar-2020 19:23:47
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
#include "milqr.h"

/* Function Declarations */
static void argInit_3x2_real32_T(float result[6]);
static void argInit_3x499_real32_T(float result[1497]);
static void argInit_7x1_real32_T(float result[7]);
static void argInit_7x500_real32_T(float result[3500]);
static float argInit_real32_T(void);
static void main_milqr(void);

/* Function Definitions */

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
 * Arguments    : float result[1497]
 * Return Type  : void
 */
static void argInit_3x499_real32_T(float result[1497])
{
    int idx0;
    int idx1;

    /* Loop over the array to initialize each element. */
    for (idx0 = 0; idx0 < 3; idx0++) {
        for (idx1 = 0; idx1 < 499; idx1++) {
            /* Set the value of the array element.
               Change this value to the value that the application requires. */
            result[idx0 + 3 * idx1] = argInit_real32_T();
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
 * Arguments    : float result[3500]
 * Return Type  : void
 */
static void argInit_7x500_real32_T(float result[3500])
{
    int idx0;
    int idx1;

    /* Loop over the array to initialize each element. */
    for (idx0 = 0; idx0 < 7; idx0++) {
        for (idx1 = 0; idx1 < 500; idx1++) {
            /* Set the value of the array element.
               Change this value to the value that the application requires. */
            result[idx0 + 7 * idx1] = argInit_real32_T();
        }
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
    static float b_fv[3500];
    float b_fv1[7];
    static float b_fv2[1497];
    float fv3[6];
    static float x[3500];
    static float u[1497];
    static float K[8982];
    boolean_T result;

    /* Initialize function 'milqr' input arguments. */
    /* Initialize function input argument 'x0'. */
    /* Initialize function input argument 'xg'. */
    /* Initialize function input argument 'u0'. */
    /* Initialize function input argument 'u_lims'. */
    /* Call the entry-point 'milqr'. */
    argInit_7x500_real32_T(b_fv);
    argInit_7x1_real32_T(b_fv1);
    argInit_3x499_real32_T(b_fv2);
    argInit_3x2_real32_T(fv3);
    milqr(b_fv, b_fv1, b_fv2, fv3, x, u, K, &result);
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
    main_milqr();
    return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
