//
// Created by Ethan on 10/11/2019.
//

#include "time_functions.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <cassert>
#include <../../pybind11/include/pybind11/pybind11.h>

using namespace std;

// function declaration
double MJD2GMST(double MJD);
double date2MJD(int M, int D, int Y, int HH, int MM, double SS);
bool valid_date(int M, int D, int Y, int HH, int MM, double SS);

int main() {
    // local variable declaration:
    double MJD = 58827.53750000009;
    double ret;

    // calling a function to get max value.
    ret = MJD2GMST(MJD);
    cout << "GMST is : " << ret << endl;
    int M, D, Y, HH, MM;
    double SS;
    M = 1;
    D = 1;
    Y = 2000;
    HH = 12;
    MM = 0;
    SS = 0;
    double ret2 = date2MJD(M, D, Y, HH, MM, SS);
    cout << "MJD is : " << ret2 << endl;
    return 0;
}


double MJD2GMST(double MJD) {
    /*
    Gives the Greenwich Mean Sidereal Time
    Inputs :
    MJD - Modified Julian Day
    Outputs :
    GMST - Greenwich Mean Sidereal Time
    */
    double GMST = 280.4606 + 360.9856473 * (MJD - 51544.5);
    GMST = GMST * M_PI / 180;
    if (GMST >= 2 * M_PI)
    {
        double diff = ceil((GMST - 2 * M_PI) / (2 * M_PI));
        GMST = GMST - 2 * M_PI * diff;
    }

    return GMST;
}

double date2MJD(int M, int D, int Y, int HH, int MM, double SS) {
    /*
     Gives the Modified Julian Date from the date and time using Vallado algorithm
        Inputs:
        M - Month number (January = 1, December = 12)
        D - Day number
        Y - Year
        HH - Hour in 24 hour clock
        MM - Minutes
        SS - Seconds
        Outputs:
        MJD - Modified Julian Date
     */
    int y, m;
    double B, day_frac;

    assert(valid_date(M, D, Y, HH, MM, SS));

    if (M<=2)
    {
        y = Y - 1;
        m = M + 12;
    }
    else
    {
        y = Y;
        m = M;
    }

    if (Y <= 1582 && M <= 10 && D <= 4)
    {
        B = -2 + ((y + 4716) / 4) - 1179;
    }
    else
    {
        B = 2 - floor(y / 100) + floor(floor(y/100)/4);
    }

    day_frac = HH / 24 + MM / 60 / 24 + SS / 60 / 60 / 24;
    double MJD = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + D + B + day_frac - 2401525.0;
    return MJD;
}

bool valid_date(int M, int D, int Y, int HH, int MM, double SS) {
    /*
    Ensures the calendar date is valid
    Inputs:
    M - Month
    Y - Year
    D - Day
    HH - Hours
    MM - Minutes
    SS - Seconds
    Outputs:
    check - true/false
    */
    bool check_1, check_2, check;
    check_1 = false;
    check_2 = false;

    cout << MM << endl;
    if (M <= 12 && D <= 31 && HH <= 24 && MM <= 60 && SS <= 60)
    {
        check_1 = true;
    }
    if (M >= 0 && D >= 0)
    {
        check_2 = true;
    }

    if (check_1 && check_2)
    {
        check = true;
    }
    else
    {
        check = false;
    }
    return check;
}

PYBIND11_MODULE(time_functions_cpp, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("date2MJD", &date2MJD, "A function which returns the sun position");
}