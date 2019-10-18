from GNC.util_funcs.py_funcs import time_functions as tf
import numpy as np
import sys
import pytest
from GNC.cmake_build_debug import time_functions_cpp as tfcpp

# Test 1: Check function works when month is after March
def test_date2MJD_1():
    M = int(5)
    D = int(10)
    Y = int(2020)
    HH = int(8)
    MM = int(5)
    SS = float(3)
    frac = 8 / 24 + 5 / 24 / 60 + 3 / 24 / 3600
    np.testing.assert_allclose(tf.date2MJD(M, D, Y, HH, MM, SS), (58979+frac), atol=1e-6) # Python test
    np.testing.assert_allclose(tfcpp.date2MJD(M, D, Y, HH, MM, SS), (58979+frac), atol=1e-6) # cpp test
    np.testing.assert_allclose(tfcpp.date2MJD(M, D, Y, HH, MM, SS), tf.date2MJD(M, D, Y, HH, MM, SS),
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
    np.testing.assert_allclose(tf.date2MJD(M, D, Y, HH, MM, SS), (58530+frac), atol=1e-6) # Python test
    np.testing.assert_allclose(tfcpp.date2MJD(M, D, Y, HH, MM, SS), (58530+frac), atol=1e-6) # cpp test
    np.testing.assert_allclose(tfcpp.date2MJD(M, D, Y, HH, MM, SS), tf.date2MJD(M, D, Y, HH, MM, SS),
                               atol=1e-6)  # compare test

# Test 3: Checks that invalid date is rejected. Valid Date function will have separate unit testing
# def test_date2MJD_3():
#     M = 13
#     D = 14
#     Y = 2019
#     HH = 15
#     MM = 5
#     SS = 56
#     with pytest.raises(Exception):
#         tf.date2MJD(M, D, Y, HH, MM, SS)
#     with pytest.raises(Exception):
#         tfcpp.date2MJD(M, D, Y, HH, MM, SS)

# Test 1: Checks that MJD is successfully turned into GMST
def test_MJD2GMST_1():

    assert True