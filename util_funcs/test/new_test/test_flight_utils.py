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
