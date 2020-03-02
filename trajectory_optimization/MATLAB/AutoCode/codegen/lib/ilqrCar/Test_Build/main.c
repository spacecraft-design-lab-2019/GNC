/******************************************************************/

/* Include Files */
#include "main.h"
#include "ilqrCar.h"


/* Function Declarations */
static void initULims(double uLims[4]);
static void initU0(double u0[1000]);
static void initXg(double xg[4]);
static void initX0(double x0[2004]);
static void main_ilqrCar(void);
static void writeToFile(char* filename, double x[2004], const int nStates, const int nRows);


static void initULims(double uLims[4])
{
  // Need to check the order of these is correct
  uLims[0] = -0.5;
  uLims[1] = -2;
  uLims[2] = 0.5;
  uLims[3] = 2;
}


static void initU0(double u0[1000])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    for (idx1 = 0; idx1 < 500; idx1++) {
      u0[idx0 + (idx1 << 1)] = 0.0;
    }
  }
}


static void initXg(double xg[4])
{
  // Set goal state
  xg[0] = 0.0;
  xg[1] = 0.0;
  xg[2] = 0.0;
  xg[3] = 0.0;
}


static void initX0(double x0[2004])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element to zero. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    for (idx1 = 0; idx1 < 501; idx1++) {
      x0[idx0 + (idx1 << 2)] = 0.0;
    }
  }
  // Set the inital state
  x0[0] = 1;
  x0[1] = 1;
  x0[2] = 4.712385;  // 3*pi/2
  x0[3] = 0;
}


static void main_ilqrCar(void)
{
  static double x0[2004];
  double xg[4];
  static double u0[1000];
  double uLims[4];
  static double x[2004];
  static double u[1000];
  static double K[4000];
  bool result;

  // Init x0
  initX0(x0);
  initXg(xg);
  initU0(u0);
  initULims(uLims);
  ilqrCar(x0, xg, u0, uLims, x, u, K, &result);

  const int nStates = 4;
  const int nRows = 501;
  char filename[] = "C_ilqr_run1";
  writeToFile(filename, x, nStates, nRows);
}


/**
 * Writes the trajectory optimization result to a csv file
 *
 * filename - char array with the desired filename (no extension)
 * x[]      - the state trajectories from the optimisation
 * nStates  - the number of states
 */
static void writeToFile(char* filename, double x[2004], const int nStates, const int nRows)
{
  filename = strcat(filename ,".csv");
  printf("Writing data to %s", filename);
 
  FILE* f = fopen(filename, "w+");
  fprintf(f, "timestep, x, y, theta, v");

  /* Loop over the array to write each state vector to the file */
  // for (int rowIdx = 0; rowIdx < nRows; rowIdx++) {
  //   fprintf(f, "\n%d", rowIdx);
  //   for (int colIdx = 0; colIdx < nStates; colIdx++) {
  //     fprintf(f, ",%f", x[(rowIdx * nStates) + colIdx]);
  //   }
  // }
  fclose(f);
  printf("\nFile successfuly written");
}


int main()
{

  main_ilqrCar();
  printf("\nOptimization Done");
  return 0;
}

