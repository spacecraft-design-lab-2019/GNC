from util_funcs import time_functions as tf
import numpy as np
import sys
import pytest
# Test 1: Check function works when month is after March
def test_date2MJD_1():
    M = 5
    D = 10
    Y = 2020
    HH = 8
    MM = 5
    SS = 3
    frac = 8 / 24 + 5 / 24 / 60 + 3 / 24 / 3600
    np.testing.assert_allclose(tf.date2MJD(M, D, Y, HH, MM, SS), (58979+frac), atol=1e-6)

# Test 2: Check function works when month is before March
def test_date2MJD_2():
    M = 2
    D = 16
    Y = 2019
    HH = 15
    MM = 5
    SS = 56
    frac = 15 / 24 + 5 / 24 / 60 + 56 / 24 / 3600
    np.testing.assert_allclose(tf.date2MJD(M, D, Y, HH, MM, SS), (58530+frac), atol=1e-6)

# Test 3: Checks that invalid date is rejected. Valid Date function will have separate unit testing
def test_date2MJD_3():
    M = 13
    D = 14
    Y = 2019
    HH = 15
    MM = 5
    SS = 56
    with pytest.raises(Exception):
        tf.date2MJD(M, D, Y, HH, MM, SS)

# Test 1: Checks that MJD is successfully turned into GMST
def test_MJD2GMST_1():

    assert False