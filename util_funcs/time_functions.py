import numpy as np
import math
def date2MJD(M, D, Y, HH, MM, SS):
    """
    Gives the Modified Julian Date from the date and time
    Inputs:
    M - Month number (January = 1, December = 12)
    D - Day number
    Y - Year
    HH - Hour in 24 hour clock
    MM - Minutes
    SS - Seconds
    Outputs:
    MJD - Modified Julian Date
    """
    if M <= 2:
        y = Y - 1
        m = M + 12
    else:
        y = Y
        m = M

    if Y <= 1582 and M <= 10 and D <= 4:
        B = -2 + ((y + 4716) / 4) - 1179
    else:
        B = y / 400 - y / 100 + y / 4

    day_frac = HH / 24 + MM / 60 / 24 + SS / 60 / 60 / 24

    MJD = 365 * y - 679004 + np.floor(B) + np.floor(30.6001 * (m + 1)) + D + day_frac;

    return MJD

def MJD2GMST(MJD):
    """
    Gives the Greenwich Mean Sidereal Time
    Inputs:
    MJD - Modified Julian Day
    Outputs:
    GMST - Greenwich Mean Sidereal Time
    """
    GMST = 280.4606 + 360.9856473 * (MJD - 51544.5);
    GMST = GMST * math.pi / 180;
    if GMST >= 2*math.pi:
        diff = np.ceil((GMST-2*math.pi) / (2*math.pi));
        GMST = GMST - 2 * math.pi * diff;

    return GMST