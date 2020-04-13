
/* Include Files */
#include "main.h"


/* Function Declarations */
static void initOptions(float result[10]);
static void initULims(float result[2][3]);
static void initU0(float result[4999][3]);
static void initB_ECI(float result[5000][3]);
static void initJ(float result[6][3]);
static void initXg(float result[7]);
static void initX0(float result[5000][7]);
static void main_milqr(void);
static void writeToFile(const char* filename, const float x[5000][7], const int nStates, const int nRows);

/* Function Definitions */

static void initOptions(float options[10])
{
  options[0] = 5000;  // N steps
  options[1] = 0.03;  // dt
  options[2] = 300;   // max iters
  options[3] = 1E-7;  // cost reduction exit tolerance
  options[4] = 1E-4;  // feedforward control change exit criterion
  options[5] = 1E-5;  // max regularizaion param allowed for exit
  options[6] = 0;     // minimum accepted cost reduction ratio
  options[7] = 1E9;   // maximum regularization parameter
  options[8] = 1E-6;  // set lambda = 0 below this value
  options[9] = 1.6;   // amount to scale dlambda by    
}


static void initULims(float ULims[2][3])
{
  for (int idx0 = 0; idx0 < 3; idx0++) {

    ULims[0][idx0] = -100.0;  // lower
    ULims[1][idx0] = 100.0;  //upper
  }
}


static void initU0(float u0[4999][3])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 4999; idx1++) {
      u0[idx1][idx0] = 0.0;
    }
  }
}


static void initB_ECI(float B_ECI[5000][3])
{
  /* Initialize magentic field as sinusoids */
  for (int idx0 = 0; idx0 < 3; idx0++) {
    for (int idx1 = 0; idx1 < 5000; idx1++) {
      float N = 0.04 * idx1; 
      if (idx0 == 0) {
        B_ECI[idx1][idx0] = (40E-6)*sin(N);
      }
      else if (idx0 == 2) {
        B_ECI[idx1][idx0] = (60E-6)*sin(N);
      }
      else {
        B_ECI[idx1][idx0] = (60E-6)*cos(N);
      }
    }
  }
}


static void initJ(float Js[6][3])
{
  // Inertia and inverse [J, Jinv]
  for (int idx0 = 0; idx0 < 3; idx0++) {
    for (int idx1 = 0; idx1 < 6; idx1++) {
      if (idx0 == idx1) {
        Js[idx1][idx0] = 0.01;  // J diagonal
      }
      else if (idx0 == (idx1 - 3)) {
        Js[idx1][idx0] = 100;  // Jinv diagonal
      }
      else {
        Js[idx1][idx0] = 0.0;
      }
    }
  }
}


/* Initialize the goal state (7x1)*/
static void initXg(float xg[7])
{
  // Goal state
  float x_goal[7] = {0.0, 0.0, 0.0, 1.0, 0, 0, 0}; // theta_z = pi
  for (int idx = 0; idx < 7; idx++) {
    xg[idx] = x_goal[idx];
  }
}


static void initX0(float x0[5000][7])
{
  // Initial state
  float x_init[7] = {0.707, 0.0, 0.0, 0.707, 0.0, 0.0, 0.0}; // theta_z = pi/2

  /* Loop over the array to initialize each element. */
  for (int idx0 = 0; idx0 < 7; idx0++) {
    for (int idx1 = 0; idx1 < 5000; idx1++) {
      x0[idx1][idx0] = x_init[idx0];
    }
  }
}


/* Set up and run the MILQR optimization */
static void main_milqr(void)
{
  static float x0[5000][7];
  float xg[7];
  static float u0[4999][3];
  float uLims[2][3];
  float B_ECI[5000][3];
  float Js[6][3];
  float options[10];
  static float x[5000][7];
  static float u[4999][3];
  static float K[4999][6][3];
  bool result;

  /* Initialize function 'milqr' input arguments. */
  initX0(x0);
  initXg(xg);
  initU0(u0);
  initULims(uLims);
  initB_ECI(B_ECI);
  initOptions(options);
  initJ(Js);

  /* Call MILQR Optimizer */
  milqr(x0, xg, u0, uLims, B_ECI, Js, options, x, u, K, &result);
  printf("\nOptimization Done\n");

  const int nStates = 7;
  const int nRows = 5000;
  char filename[] = "c_milqr_run.csv";
  writeToFile(filename, x, nStates, nRows);
}

/**
 * Writes the trajectory optimization result to a csv file
 *
 * filename - char array with the desired filename (no extension)
 * x[]      - the state trajectories from the optimization
 * nStates  - the number of states
 */
static void writeToFile(const char* filename, const float x[5000][7], const int nStates, const int nRows)
{
  printf("\nWriting data to %s", filename);
 
  FILE* fp = fopen(filename, "w");
  fprintf(fp, "%s", "timestep, x, y, theta, v");

  /* Loop over the array to write each state vector to the file */
  for (int rowIdx = 0; rowIdx < nRows; rowIdx++) {
    fprintf(fp, "\n%d", rowIdx);
    for (int colIdx = 0; colIdx < nStates; colIdx++) {
      fprintf(fp, ",%f", x[rowIdx][colIdx]);
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

