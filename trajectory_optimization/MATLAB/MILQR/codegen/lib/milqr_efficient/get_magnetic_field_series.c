/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: get_magnetic_field_series.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 16:22:38
 */

/* Include Files */
#include <math.h>
#include <string.h>
#include "get_magnetic_field_series.h"

/* Function Declarations */
static double b_norm(const double x[3]);
static void get_P_coefficients(double x, double P[49]);
static void get_Pd_coefficients(const double P[49], double x, double Pd[49]);
static void get_magnetic_field(double lat, double lon, double alt, double year,
  double B_NED[3]);
static double rt_roundd(double u);
static void two_body(const double x[6], double xdot[6]);

/* Function Definitions */

/*
 * Arguments    : const double x[3]
 * Return Type  : double
 */
static double b_norm(const double x[3])
{
  double y;
  double scale;
  double absxk;
  double t;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = 1.0 + y * t * t;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = 1.0 + y * t * t;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * sqrt(y);
}

/*
 * Arguments    : double x
 *                double P[49]
 * Return Type  : void
 */
static void get_P_coefficients(double x, double P[49])
{
  double P_tmp;
  int n;
  int m;
  int prev2;
  int b_P_tmp;
  double c_P_tmp;
  memset(&P[0], 0, 49U * sizeof(double));
  P[0] = 1.0;
  P_tmp = sqrt(1.0 - x * x);
  P[8] = P_tmp;
  for (n = 0; n < 6; n++) {
    for (m = 0; m < 7; m++) {
      if ((1 + n != 1) || (m != 1)) {
        prev2 = n - 1;
        if (n - 1 < 0) {
          prev2 = 0;
        }

        if (m < 1 + n) {
          b_P_tmp = m * m;
          c_P_tmp = (1 + n) * (1 + n) - b_P_tmp;
          P[(n + 7 * m) + 1] = (2.0 * (1.0 + (double)n) - 1.0) / sqrt(c_P_tmp) *
            x * P[n + 7 * m] - sqrt((double)(n * n - b_P_tmp) / c_P_tmp) *
            P[prev2 + 7 * m];
        } else {
          P[(n + 7 * m) + 1] = sqrt(1.0 - 1.0 / (2.0 * (double)m)) * P_tmp * P[n
            + 7 * (m - 1)];
        }
      }
    }
  }
}

/*
 * Arguments    : const double P[49]
 *                double x
 *                double Pd[49]
 * Return Type  : void
 */
static void get_Pd_coefficients(const double P[49], double x, double Pd[49])
{
  int n;
  int m;
  int Pd_tmp;
  int prev2;
  double b_Pd_tmp;
  memset(&Pd[0], 0, 49U * sizeof(double));
  Pd[8] = x;
  for (n = 0; n < 6; n++) {
    for (m = 0; m < 7; m++) {
      if ((1 + n != 1) || (m != 1)) {
        if (m < 1 + n) {
          prev2 = n - 1;
          if (n - 1 < 0) {
            prev2 = 0;
          }

          Pd_tmp = m * m;
          b_Pd_tmp = (1 + n) * (1 + n) - Pd_tmp;
          Pd[(n + 7 * m) + 1] = (2.0 * (1.0 + (double)n) - 1.0) / sqrt(b_Pd_tmp)
            * (x * Pd[n + 7 * m] - sqrt(1.0 - x * x) * P[n + 7 * m]) - sqrt
            ((double)(n * n - Pd_tmp) / b_Pd_tmp) * Pd[prev2 + 7 * m];
        } else {
          Pd_tmp = n + 7 * (m - 1);
          Pd[(n + 7 * m) + 1] = sqrt(1.0 - 1.0 / (2.0 * (double)m)) * (sqrt(1.0
            - x * x) * Pd[Pd_tmp] + x * P[Pd_tmp]);
        }
      }
    }
  }
}

/*
 * Arguments    : double lat
 *                double lon
 *                double alt
 *                double year
 *                double B_NED[3]
 * Return Type  : void
 */
static void get_magnetic_field(double lat, double lon, double alt, double year,
  double B_NED[3])
{
  double B_r;
  double B_lat;
  double B_lon;
  double P_tmp;
  double P[49];
  double Pd[49];
  int n;
  int m;
  int i0;
  double coef_tmp;
  double b_coef_tmp;
  double c_coef_tmp;
  static const double g[121] = { 0.0, -29442.0, -2445.1, 1350.7, 907.6, -232.6,
    70.0, 81.6, 24.2, 5.4, -1.9, 0.0, -1501.0, 3012.9, -2352.3, 813.7, 360.1,
    67.7, -76.1, 8.8, 8.8, -6.3, 0.0, 0.0, 1676.7, 1225.6, 120.4, 192.4, 72.7,
    -6.8, -16.9, 3.1, 0.1, 0.0, 0.0, 0.0, 582.0, -334.9, -140.9, -129.9, 51.8,
    -3.2, -3.3, 0.5, 0.0, 0.0, 0.0, 0.0, 70.4, -157.5, -28.9, 15.0, -20.6, 0.7,
    -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 4.1, 13.2, 9.4, 13.4, -13.3, 1.8, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -70.9, -2.8, 11.7, -0.1, -0.7, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.8, -15.9, 8.7, 2.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0,
    9.1, 2.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.5, -1.8, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.6 };

  static const double g_sv[121] = { 0.0, 10.3, -8.7, 3.4, -0.7, -0.2, -0.3, 0.3,
    0.2, 0.0, 0.0, 0.0, 18.1, -3.3, -5.5, 0.2, 0.5, -0.1, -0.2, 0.0, 0.0, 0.0,
    0.0, 0.0, 2.1, -0.7, -9.1, -1.3, -0.7, -0.5, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0,
    -10.1, 4.1, -0.1, 2.1, 1.3, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.3, 1.4,
    -1.2, 0.1, -0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.9, 0.3, -0.6, 0.4, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6, -0.8, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.2, -0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  double d_coef_tmp;
  double e_coef_tmp;
  static const double h[121] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 4797.1, -2845.6, -115.3, 283.3, 47.3, -20.8, -54.1, 10.1,
    -21.6, 3.2, 0.0, 0.0, -641.9, 244.9, -188.7, 197.0, 33.2, -19.5, -18.3, 10.8,
    -0.4, 0.0, 0.0, 0.0, -538.4, 180.9, -119.3, 58.9, 5.7, 13.3, 11.8, 4.6, 0.0,
    0.0, 0.0, 0.0, -329.5, 16.0, -66.7, 24.4, -14.6, -6.8, 4.4, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.2, 7.3, 3.4, 16.2, -6.9, -7.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    62.6, -27.4, 5.7, 7.8, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.2, -9.1,
    1.0, -4.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.1, -4.0, -2.8, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.4, -1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -8.7 };

  static const double h_sv[121] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -26.6, -27.4, 8.2, -1.3, 0.6, 0.0, 0.8, -0.3, 0.0, 0.0, 0.0,
    0.0, -14.1, -0.4, 5.3, 1.7, -2.1, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8,
    2.9, -1.2, -0.7, -0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.2, 3.4, 0.2,
    -0.3, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, -0.6, -0.2, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.1, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  double coef;
  int i1;
  double d3;

  /*      /$ */
  /*      lat is geocentric latitude in degrees */
  /*      lon is longitude in degrees */
  /*      alt is altitude in km */
  /*      year is the fractional year (include months/days essentially) */
  /*   */
  /*      outputs B vector in NED */
  /*      $/ */
  lat = 90.0 - lat;
  lat *= 0.017453292519943295;
  lon *= 0.017453292519943295;

  /*      // radius of earth */
  /*      // year since 2015 for secular variation */
  /*      // magnetic field components */
  B_r = 0.0;
  B_lat = 0.0;
  B_lon = 0.0;

  /*      // IGRF model B field calculation, n/m sets order (5) */
  /*  get coefficients */
  P_tmp = cos(lat);
  get_P_coefficients(P_tmp, P);
  get_Pd_coefficients(P, P_tmp, Pd);
  for (n = 0; n < 6; n++) {
    for (m = 0; m < 7; m++) {
      i0 = (n + 11 * m) + 1;
      coef_tmp = (double)m * lon;
      b_coef_tmp = pow(6371.2 / (6371.2 + alt), (1.0 + (double)n) + 2.0);
      c_coef_tmp = g[i0] + (year - 2015.0) * g_sv[i0];
      d_coef_tmp = sin(coef_tmp);
      e_coef_tmp = h[i0] + (year - 2015.0) * h_sv[i0];
      coef_tmp = cos(coef_tmp);
      coef = b_coef_tmp * (c_coef_tmp * coef_tmp + e_coef_tmp * d_coef_tmp);

      /*              // Radial component */
      i1 = (n + 7 * m) + 1;
      B_r += ((1.0 + (double)n) + 1.0) * coef * P[i1];

      /*              // Colatitudinal component */
      B_lat -= coef * Pd[i1];

      /*              // Address singularity at colatitude of 0 */
      d3 = sin(lat);
      if (d3 == 0.0) {
        /*                  // Longitudinal component */
        B_lon += -P_tmp * b_coef_tmp * (-(g[i0] + (year - 2015.0) * g_sv[i0]) *
          sin((double)m * lon) + (h[i0] + (year - 2015.0) * h_sv[i0]) * cos
          ((double)m * lon)) * Pd[i1];
      } else {
        B_lon += -1.0 / d3 * b_coef_tmp * (double)m * (-c_coef_tmp * d_coef_tmp
          + e_coef_tmp * coef_tmp) * P[i1];
      }
    }
  }

  /*      // NED (North, East, Down) coordinate frame */
  B_NED[0] = -B_lat;
  B_NED[1] = B_lon;
  B_NED[2] = -B_r;
}

/*
 * Arguments    : double u
 * Return Type  : double
 */
static double rt_roundd(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/*
 * Arguments    : const double x[6]
 *                double xdot[6]
 * Return Type  : void
 */
static void two_body(const double x[6], double xdot[6])
{
  double c;
  c = pow(b_norm(*(double (*)[3])&x[0]), 3.0);
  xdot[0] = x[3];
  xdot[3] = -398600.0 * x[0] / c;
  xdot[1] = x[4];
  xdot[4] = -398600.0 * x[1] / c;
  xdot[2] = x[5];
  xdot[5] = -398600.0 * x[2] / c;
}

/*
 * Arguments    : const double x0[6]
 *                const double t[100]
 *                double B_eci_vec[297]
 *                double X[600]
 * Return Type  : void
 */
void get_magnetic_field_series(const double x0[6], const double t[100], double
  B_eci_vec[297], double X[600])
{
  int i;
  int b_i;
  double dt;
  double b_X[6];
  double c_X[6];
  double d0;
  double k1[6];
  double k2[6];
  double k3[6];
  double dv0[9];
  double d1;
  double r_ecef[3];
  double lon;
  double lat;
  double vec_ned[3];
  double d2;
  int B_eci_vec_tmp;

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  This function takes in a series of positions in ECI (in km) as well as a  */
  /*  time series in MJD and outputs a series of ECI magnetic field vectors */
  /*  INPUT: */
  /*  x0 - [r0;v0;q0;w0]; */
  /*  t - time vector in MJD */
  /*  OUTPUT: */
  /*  B_eci_vec - magnetic field in ECI in Teslas */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  This function takes in an initial position and velocity of a satellite along  */
  /*  with a time series and outputs the positions over the time using a simple */
  /*  point mass two body approximation */
  /*  INPUTS: */
  /*  x0 - initial state */
  /*      r0 - initial position in km */
  /*      v0 - initial velocity in km/s */
  /*      q0 - initial quaternion */
  /*      w0 - initial rotation rate */
  /*  t - time series in seconds */
  /*  OUTPUTS: */
  /*  X - array of states with columns: */
  /*      r - position in Earth radii */
  /*      v - velocity in km/s */
  /*  initialize empty array */
  memset(&X[0], 0, 600U * sizeof(double));
  for (i = 0; i < 6; i++) {
    X[i] = x0[i];
  }

  /*  fill it */
  for (b_i = 0; b_i < 99; b_i++) {
    dt = (t[b_i + 1] - t[b_i]) * 86400.0;
    two_body(*(double (*)[6])&X[6 * b_i], b_X);
    for (i = 0; i < 6; i++) {
      d0 = dt * b_X[i];
      k1[i] = d0;
      c_X[i] = X[i + 6 * b_i] + d0 / 2.0;
    }

    two_body(c_X, b_X);
    for (i = 0; i < 6; i++) {
      d0 = dt * b_X[i];
      k2[i] = d0;
      c_X[i] = X[i + 6 * b_i] + d0 / 2.0;
    }

    two_body(c_X, b_X);
    for (i = 0; i < 6; i++) {
      d0 = dt * b_X[i];
      k3[i] = d0;
      c_X[i] = X[i + 6 * b_i] + d0;
    }

    two_body(c_X, b_X);
    for (i = 0; i < 6; i++) {
      b_X[i] = X[i + 6 * b_i] + 0.16666666666666666 * (((k1[i] + 2.0 * k2[i]) +
        2.0 * k3[i]) + dt * b_X[i]);
    }

    for (i = 0; i < 6; i++) {
      X[i + 6 * (b_i + 1)] = b_X[i];
    }
  }

  for (b_i = 0; b_i < 99; b_i++) {
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  This function takes in a position in ECI (in km) and outputs */
    /*  a B vector in teslas that's in the ECI frame */
    /*  INPUTS: */
    /*    r - position in km */
    /*    t - time in MJD */
    /*  first, translate the position into a lat and long */
    /*  this function converts between position in ECI and time to lat, lon and */
    /*  alt */
    /*  INPUTS: */
    /*    r - position in km */
    /*    t - time in MJD */
    /*  OUTPUTS: */
    /*    lat - latitude in deg */
    /*    lon - longitude in deg */
    /*    alt - altitude in km */
    /*  convert to ecef */
    /*  converts a position from eci to ecef */
    /*  find days since 1/1/2000 at 12h */
    /*  converts between MJD and J2000 (days since 1/1/2000 at 12h) */
    /*  MJD epoch: 0h Nov 17, 1858 */
    /*  J2000 epoch: 12h Jan 1, 2000 */
    dt = (280.4606 + 360.9856473 * (t[b_i] - 51544.5)) * 3.1415926535897931 /
      180.0;

    /*  creates a z axis rotation matrix (rad) */
    /*  get the rotation matrix */
    d0 = cos(dt);
    dv0[0] = d0;
    d1 = sin(dt);
    dv0[3] = d1;
    dv0[6] = 0.0;
    dv0[1] = -d1;
    dv0[4] = d0;
    dv0[7] = 0.0;
    dv0[2] = 0.0;
    dv0[5] = 0.0;
    dv0[8] = 1.0;
    for (i = 0; i < 3; i++) {
      r_ecef[i] = 0.0;
      r_ecef[i] = (dv0[i] * X[6 * b_i] + dv0[i + 3] * X[1 + 6 * b_i]) + dv0[i +
        6] * X[2 + 6 * b_i];
    }

    /*  get lat and lon */
    lon = atan2(r_ecef[1], r_ecef[0]) * 180.0 / 3.1415926535897931;
    dt = b_norm(r_ecef);
    lat = asin(r_ecef[2] / dt) * 180.0 / 3.1415926535897931;

    /*  find B_NED */
    /*  gets the year out of MJD */
    get_magnetic_field(lat, lon, dt - 6378.0, rt_roundd((t[b_i] / 365.25 +
      1858.0) + 0.8794520547945206), vec_ned);

    /*  convert to eci */
    lat = lat * 3.1415926535897931 / 180.0;
    lon = lon * 3.1415926535897931 / 180.0;

    /*  this function converts a vector from ned to eci */
    /*  first find ENU */
    d0 = sin(lon);
    dv0[0] = -d0;
    d1 = cos(lon);
    dv0[1] = d1;
    dv0[2] = 0.0;
    d2 = sin(lat);
    dv0[3] = -d2 * d1;
    dv0[4] = -d2 * d0;
    dt = cos(lat);
    dv0[5] = dt;
    dv0[6] = dt * d1;
    dv0[7] = dt * d0;
    dv0[8] = d2;
    for (i = 0; i < 3; i++) {
      B_eci_vec_tmp = i + 3 * b_i;
      B_eci_vec[B_eci_vec_tmp] = 0.0;
      B_eci_vec[B_eci_vec_tmp] = (dv0[i] * vec_ned[0] + dv0[i + 3] * vec_ned[1])
        + dv0[i + 6] * vec_ned[2];
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void get_magnetic_field_series_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void get_magnetic_field_series_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for get_magnetic_field_series.c
 *
 * [EOF]
 */
