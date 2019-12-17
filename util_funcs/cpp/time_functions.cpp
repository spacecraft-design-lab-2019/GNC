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
namespace py = pybind11;

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
    M = 5;
    D = 10;
    Y = 2020;
    HH = 8;
    MM = 5;
    SS = 3;
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
    double T = (MJD - 51544.5) / 36525.0;
    double gmst = (67310.54841 + (876600.0*3600.0 + 8640184.812866) * T
            + 0.093104 * pow(T,2) - 6.2 * pow(10, -6) * pow(T, 3)) ;
    gmst = fmod(gmst, 86400);
    if (gmst < 0)
    {
        gmst += 86400;
    }
    gmst = gmst * M_PI / 180 / 240;



    return gmst;
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
    double y, m;
    double B, day_frac;

    // assert(valid_date(M, D, Y, HH, MM, SS)); // commenting out for now to avoid the issue of ending when the assert fails

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
        B = -2 + ((y + 4716.0) / 4.0) - 1179.0;
    }
    else
    {
        B = 2 - floor(y / 100.0) + floor(floor(y/100.0)/4.0);
    }

    day_frac = HH / 24.0 + MM / 60.0 / 24.0 + SS / 60.0 / 60.0 / 24.0;
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
    m.doc() = "Time Functions"; // optional module docstring
    m.def("valid_date", &valid_date, "Returns whether the date is valid or not");
    m.def("date2MJD", &date2MJD, "Converts the date to MJD");
    m.def("MJD2GMST", &MJD2GMST, "Converts MJD to GMST");
}
