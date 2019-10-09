import math
import numpy as np

def sun_position(MJD):
    """
    This is using the equations given in Motenbruck and Gill's Satellite Orbits book
    Inputs:
    MJD - Modified Julian Day (J2000) as a Real Number
    Outputs:
    r_vec - numpy array with x, y, z of Sun position in ECI at input time
    """
    JD = MJD + 2400000.5
    OplusW = 282.94
    T = (JD - 2451545.0) / 36525

    M = 357.5256 + 35999.049 * T

    long = OplusW + M + 6892 / 3600 * math.sin(M) + 72 / 3600 * math.sin(2*M)
    r_mag = (149.619 - 2.499 * math.cos(M) - 0.021 * math.cos(2*M)) * 10**6

    epsilon = 23.43929111
    r_vec = np.array([r_mag * math.cos(long), r_mag * math.sin(long) * math.cos(epsilon), r_mag * math.sin(long) * math.sin(epsilon)])

    return r_vec

def sat_sun_vect(r, MJD):
    """
    Returns the unit vector from the satellite to the Sun in ECI coordinates
    Inputs:
    r - ECI position of the satellite
    MJD - Julian Day (J2000) as a Real Number
    Outputs:
    r_sat_sun - numpy array giving the unit direction to the Sun from the satellite
    """

    r_sun = sun_position(MJD)
    r_sat_sun = r_sun - r

    r_sat_sun = r_sat_sun / np.linalg.norm(r_sat_sun)


    return r_sat_sun