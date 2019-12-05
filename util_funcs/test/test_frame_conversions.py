import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)

from util_funcs.py_funcs import frame_conversions as fc
import numpy as np
import pytest
import math
import frame_conversions_cpp as fccpp


def test_eci2ecef_1():
    GMST = math.pi
    R_pred = np.array([[-1, 0, 0],
                      [0, -1, 0],
                       [0, 0, 1]])
    np.testing.assert_allclose(fc.eci2ecef(GMST), R_pred, atol=1e-6) # Python test
    np.testing.assert_allclose(fccpp.eci2ecef(GMST), R_pred, atol=1e-6) # cpp test
    np.testing.assert_allclose(fc.eci2ecef(GMST), fccpp.eci2ecef(GMST),
                               atol=1e-6)  # compare test

def test_eci2ecef_2():
    GMST = math.pi / 2
    R_pred = np.array([[0, 1, 0],
                       [-1, 0, 0],
                       [0, 0, 1]])
    np.testing.assert_allclose(fc.eci2ecef(GMST), R_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(fccpp.eci2ecef(GMST), R_pred, atol=1e-6)  # cpp test
    np.testing.assert_allclose(fc.eci2ecef(GMST), fccpp.eci2ecef(GMST),
                               atol=1e-6)  # compare test

def test_ecef2enu_1():
    lat = 0
    lon = math.pi / 6
    test_vec = np.array([1, 0, 0])
    R_pred = np.array([[-1/2, math.sqrt(3)/2 ,0 ],
                       [0, 0, 1],
                       [math.sqrt(3)/2, 1/2, 0]])
    test_rot_vec = np.array([-1/2, 0, math.sqrt(3)/2])
    np.testing.assert_allclose(fc.ecef2enu(lat, lon), R_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(fccpp.ecef2enu(lat, lon), R_pred, atol=1e-6)  # cpp test
    np.testing.assert_allclose(fc.ecef2enu(lat, lon), fccpp.ecef2enu(lat, lon),
                               atol=1e-6)  # compare test

    np.testing.assert_allclose(fc.ecef2enu(lat, lon) @ test_vec, test_rot_vec, atol=1e-6) # checking rotations
    np.testing.assert_allclose(fccpp.ecef2enu(lat, lon) @ test_vec, test_rot_vec, atol=1e-6)  # checking rotations
    np.testing.assert_allclose(fc.ecef2enu(lat, lon) @ test_vec, fccpp.ecef2enu(lat, lon) @ test_vec,
                               atol=1e-6)  # compare test

def test_ecef2lla_1():
    r = np.array([7000, 7000, 0])
    alt_pred = math.sqrt(2*7000**2) - 6378.1
    lat_pred = 0
    lon_pred = math.pi / 4
    lla_pred = np.array([lat_pred, lon_pred, alt_pred])
    py_lat, py_lon, py_alt = fc.ecef2lla(r)
    py_lla = np.array([py_lat, py_lon, py_alt])

    np.testing.assert_allclose(py_lla, lla_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(fccpp.ecef2lla(r), lla_pred, atol=1e-6)  # cpp test
    np.testing.assert_allclose(py_lla, fccpp.ecef2lla(r),
                               atol=1e-6)  # compare test

def test_ecef2lla_2():
    r = np.array([0, 7000, 7000])
    alt_pred = math.sqrt(2*7000**2) - 6378.1
    lat_pred = math.pi / 4
    lon_pred = math.pi / 2
    lla_pred = np.array([lat_pred, lon_pred, alt_pred])
    py_lat, py_lon, py_alt = fc.ecef2lla(r)
    py_lla = np.array([py_lat, py_lon, py_alt])

    np.testing.assert_allclose(py_lla, lla_pred, atol=1e-6)  # Python test
    np.testing.assert_allclose(fccpp.ecef2lla(r), lla_pred, atol=1e-6)  # cpp test
    np.testing.assert_allclose(py_lla, fccpp.ecef2lla(r),
                               atol=1e-6)  # compare test

def test_ecef2lla_3():
        r = np.array([0, 0, 7000])
        alt_pred = math.sqrt(7000 ** 2) - 6378.1
        lat_pred = math.pi / 2
        lon_pred = 0
        lla_pred = np.array([lat_pred, lon_pred, alt_pred])
        py_lat, py_lon, py_alt = fc.ecef2lla(r)
        py_lla = np.array([py_lat, py_lon, py_alt])

        np.testing.assert_allclose(py_lla, lla_pred, atol=1e-6)  # Python test
        np.testing.assert_allclose(fccpp.ecef2lla(r), lla_pred, atol=1e-6)  # cpp test
        np.testing.assert_allclose(py_lla, fccpp.ecef2lla(r),
                                   atol=1e-6)  # compare test
        
        
