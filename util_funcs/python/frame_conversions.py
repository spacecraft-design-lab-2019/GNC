import math
from math import cos
from math import sin
import numpy as np

def eci2ecef(GMST):
    """
    Rotation matrix from ECI to ECEF coordinates
    Inputs:
    GMST - Greenwich Mean Sidereal Time
    Outputs:
    R - Rotation matrix from ECI to ECEF
    """
    R = np.array([[cos(GMST), sin(GMST), 0],
                     [ -sin(GMST), cos(GMST), 0],
                     [ 0, 0, 1]])
    return R

def ecef2enu(lat, lon):
    """
    Rotation matrix from ECEF to ENU coordinates
    Inputs:
    lat - Latitude in radians
    lon - Longitude in radians
    Outputs:
    R - Rotation matrix from ECEF to ENU
    """
    Ehat = np.array([-sin(lon), cos(lon), 0])
    Nhat = np.array([-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)])
    Uhat = np.array([cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)])

    R = np.column_stack((Ehat, Nhat, Uhat))
    R = np.transpose(R)
    return R

def ecef2lla(r):
    """
    Gives longitude, latitude, and altitude from ECEF position vector
    Inputs:
    r - Position vector in ECEF
    Outputs:
    lat - geocentric latitude in radians
    lon - longitude in radians
    alt - altitude
    """
    R_earth = 6378.1
    lat = math.asin(r[2] / np.linalg.norm(r))
    lon = math.atan2(r[1], r[0])
    alt = np.linalg.norm(r) - R_earth;

    return lat, lon, alt