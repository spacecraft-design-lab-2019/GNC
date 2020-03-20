/**
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 17-Mar-2020 19:06:03
 */
/*************************************************************************/


/* Include Files */
#include "main.h"
#include "milqr.h"

/* Function Declarations */
static void initULims(double uLims[6]);
static void initU0(double u0[1497]);
static void initXg(double xg[7]);
static void initX0(double x0[3500]);
static void main_milqr(void);
static void writeToFile(const char* filename, const double x[2004], const int nStates, const int nRows);


/* Initialize the control limits (3x2)*/
static void initULims(double uLims[6])
{
  // The order of these is important (lower first)
  uLims[0] = -0.5;
  uLims[1] = -0.5;
  uLims[2] = -0.5;
  uLims[3] = 0.5;
  uLims[4] = 0.5;
  uLims[5] = 0.5;
}


/* Initialize the controls (3x499)*/
static void initU0(double u0[1497])
{
  /* Loop over the array to initialize each element. */
  for (int idx0 = 0; idx0 < 3; idx0++) {
    for (int idx1 = 0; idx1 < 499; idx1++) {
      u0[idx0 + 3 * idx1] = 0.0;
    }
  }
}


/* Initialize the goal state (7x1)*/
static void initXg(double xg[7])
{
  // Goal state
  double x_goal[7] = {0.0, 0.0, 0.0, 1.0, 0, 0, 0}; // theta_z = pi
  for (int idx = 0; idx < 7; idx++) {
    xg[idx] = x_goal[idx];
  }
}


/* Initialize the trajectory (7x500)*/
static void initX0(double x0[3500])
{
  // Initial state
  double x_init[7] = {0.707, 0.0, 0.0, 0.707, 0.0, 0.0, 0.0}; // theta_z = pi/2

  /* Loop over the array to initialize each element. */
  for (int idx0 = 0; idx0 < 7; idx0++) {
    for (int idx1 = 0; idx1 < 500; idx1++) {
      x0[idx0 + 7 * idx1] = x_init[idx0];
    }
  }
}


/* Set up and run the MILQR optimization */
static void main_milqr(void)
{
  static double x0[3500];
  double xg[7];
  static double u0[1497];
  double uLims[6];
  static double x[3500];
  static double u[1497];
  static double K[8982];
  bool result;

  /* Initialize function 'milqr' input arguments. */
  initX0(x0);
  initXg(xg);
  initU0(u0);
  initULims(uLims);

  /* Call MILQR Optimizer */
  milqr(x0, xg, u0, uLims, x, u, K, &result);
  printf("\nOptimization Done\n");

  const int nStates = 7;
  const int nRows = 500;
  char filename[] = "c_milqr_run1.csv";
  writeToFile(filename, x, nStates, nRows);
}


/**
 * Writes the trajectory optimization result to a csv file
 *
 * filename - char array with the desired filename (no extension)
 * x[]      - the state trajectories from the optimization
 * nStates  - the number of states
 */
static void writeToFile(const char* filename, const double x[3500], const int nStates, const int nRows)
{
  printf("\nWriting data to %s", filename);
 
  FILE* fp = fopen(filename, "w");
  fprintf(fp, "%s", "timestep, x, y, theta, v");

  /* Loop over the array to write each state vector to the file */
  for (int rowIdx = 0; rowIdx < nRows; rowIdx++) {
    fprintf(fp, "\n%d", rowIdx);
    for (int colIdx = 0; colIdx < nStates; colIdx++) {
      fprintf(fp, ",%f", x[(rowIdx * nStates) + colIdx]);
    }
  }
  fclose(fp);
  printf("\nFile successfuly written");
}


int main()
{

  main_milqr();
  return 0;
}

