import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)

from util_funcs.py_funcs import flight_utils as fut
import numpy as np
import pytest
import math
from util_funcs.py_funcs import frame_conversions as fc
from util_funcs.py_funcs import time_functions as tf
from util_funcs.py_funcs import sun_utils as su

def test_ECI_2_ECEF_1():
    GMST = math.pi
    R_pred = np.array([[-1, 0, 0],
                      [0, -1, 0],
                       [0, 0, 1]])
    np.testing.assert_allclose(np.array(fut.ECI_to_ECEF(GMST)), R_pred, atol=1e-6) # Python test
    np.testing.assert_allclose(fc.eci2ecef(GMST), np.array(fut.ECI_to_ECEF(GMST)),
                               atol=1e-6)  # compare test

def test_ECI_2_ECEF_2():
    GMST = math.pi / 2
    R_pred = np.array([[0, 1, 0],
                       [-1, 0, 0],
                       [0, 0, 1]])
    np.testing.assert_allclose(np.array(fut.ECI_to_ECEF(GMST)), R_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(fc.eci2ecef(GMST), np.array(fut.ECI_to_ECEF(GMST)),
                               atol=1e-6)  # compare test
    
# Test 1: Check function works when month is after March
def test_date2MJD_1():
    M = int(5)
    D = int(10)
    Y = int(2020)
    HH = int(8)
    MM = int(5)
    SS = float(3)
    frac = 8 / 24 + 5 / 24 / 60 + 3 / 24 / 3600
    np.testing.assert_allclose(np.array(fut.date2MJD(M, D, Y, HH, MM, SS)), (58979+frac), atol=1e-6) # Python test
    np.testing.assert_allclose(np.array(fut.date2MJD(M, D, Y, HH, MM, SS)), tf.date2MJD(M, D, Y, HH, MM, SS),
                               atol=1e-6) # compare test

# Test 2: Check function works when month is before March
def test_date2MJD_2():
    M = 2
    D = 16
    Y = 2019
    HH = 15
    MM = 5
    SS = 56
    frac = 15 / 24 + 5 / 24 / 60 + 56 / 24 / 3600
    np.testing.assert_allclose(np.array(fut.date2MJD(M, D, Y, HH, MM, SS)), (58530+frac), atol=1e-6) # Python test
    np.testing.assert_allclose(np.array(fut.date2MJD(M, D, Y, HH, MM, SS)), tf.date2MJD(M, D, Y, HH, MM, SS),
                               atol=1e-6)  # compare test

# Test 1: Checks that MJD is successfully turned into GMST
def test_MJD2GMST_1():
    MJD = 48854.50972222211
    gmst_check = 152.578787810 * math.pi / 180 # example is from Vallado
    np.testing.assert_allclose(np.array(fut.MJD2GMST(MJD)), gmst_check, atol=1e-6) # Python test
    np.testing.assert_allclose(np.array(fut.MJD2GMST(MJD)), tf.MJD2GMST(MJD), atol=1e-6)  # compare test
    
def test_ECEF_to_LLA_1():
    r_ECEF = [7000, 7000, 0]
    rad_Earth = 6378.1
    alt_pred = math.sqrt(2*7000**2) - 6378.1
    lat_pred = 0
    lon_pred = math.pi / 4
    lla_pred = np.array([lat_pred, lon_pred, alt_pred])
    py_lat, py_lon, py_alt = fut.ECEF_to_LLA(r_ECEF, rad_Earth)
    py_lla = np.array([py_lat, py_lon, py_alt])
    
    py_lat2, py_lon2, py_alt2 = fc.ecef2lla(r_ECEF)
    py_lla2 = np.array([py_lat2, py_lon2, py_alt2])

    np.testing.assert_allclose(py_lla, lla_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(py_lla, py_lla2,
                               atol=1e-6)  # compare test
    
def test_ECEF_to_LLA_2():
    r_ECEF = [0, 7000, 7000]
    rad_Earth = 6378.1
    alt_pred = math.sqrt(2*7000**2) - 6378.1
    lat_pred = math.pi / 4
    lon_pred = math.pi / 2
    lla_pred = np.array([lat_pred, lon_pred, alt_pred])
    py_lat, py_lon, py_alt = fut.ECEF_to_LLA(r_ECEF, rad_Earth)
    py_lla = np.array([py_lat, py_lon, py_alt])
    
    py_lat2, py_lon2, py_alt2 = fc.ecef2lla(r_ECEF)
    py_lla2 = np.array([py_lat2, py_lon2, py_alt2])

    np.testing.assert_allclose(py_lla, lla_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(py_lla, py_lla2,
                               atol=1e-6)  # compare test
    
def test_ECEF_to_LLA_3():
    r_ECEF = [0, 0, 7000]
    rad_Earth = 6378.1
    alt_pred = math.sqrt(7000 ** 2) - 6378.1
    lat_pred = math.pi / 2
    lon_pred = 0
    lla_pred = np.array([lat_pred, lon_pred, alt_pred])
    py_lat, py_lon, py_alt = fut.ECEF_to_LLA(r_ECEF, rad_Earth)
    py_lla = np.array([py_lat, py_lon, py_alt])
    
    py_lat2, py_lon2, py_alt2 = fc.ecef2lla(r_ECEF)
    py_lla2 = np.array([py_lat2, py_lon2, py_alt2])
    
    np.testing.assert_allclose(py_lla, lla_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(py_lla, py_lla2,
                               atol=1e-6)  # compare test
    
def test_sun_position_1():
    MJD = 51622 # J2000

    JPL_check = np.array([1.489004447920312E+08, -3.382674902835714E+06, 1.046597880064510E+02]) #given on this MJD
    JPL_check_unit = JPL_check / np.linalg.norm(JPL_check)

    py_unit = fut.sun_position(MJD) / np.linalg.norm(su.sun_position(MJD))
    I_vec = np.array([1, 0, 0])
    ang1 = math.acos(np.dot(py_unit, I_vec))
    ang2 = math.acos(np.dot(JPL_check_unit, I_vec))

    np.testing.assert_allclose(ang1, ang2, atol=5e-5) # Python test, checking that the approximate angle is fairly close to JPL
    np.testing.assert_allclose(su.sun_position(MJD),fut.sun_position(MJD), atol=1e-6) # compare test

def test_sat_sun_vect_1():
    r = [7000, 500, 1000]
    MJD = 51622
    check = su.sun_position(MJD) - r
    check = check / np.linalg.norm(check)

    np.testing.assert_allclose(fut.sat_sun_vect(r, MJD), check, atol=1e-6)
    np.testing.assert_allclose(fut.sat_sun_vect(r, MJD), su.sat_sun_vect(r, MJD), atol=1e-6)
    
def test_ecef2enu_1():
    lat = 0
    lon = math.pi / 6
    test_vec = [1, 0, 0]
    R_pred = np.array([[-1/2, math.sqrt(3)/2 ,0 ],
                       [0, 0, 1],
                       [math.sqrt(3)/2, 1/2, 0]])
    test_rot_vec = np.array([-1/2, 0, math.sqrt(3)/2])
    np.testing.assert_allclose(fut.ecef2enu(lat, lon), R_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(fc.ecef2enu(lat, lon), fut.ecef2enu(lat, lon),
                               atol=1e-6)  # compare test

    np.testing.assert_allclose(np.array(fut.ecef2enu(lat, lon)) @ test_vec, test_rot_vec, atol=1e-6) # checking rotations
    np.testing.assert_allclose(fc.ecef2enu(lat, lon) @ test_vec, np.array(fut.ecef2enu(lat, lon)) @ test_vec,
                               atol=1e-6)  # compare test
