
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)


from util_funcs.py_funcs import sun_utils as su
# import sun_utils as su
import numpy as np
import pytest
import math
import sun_utils_cpp as sucpp

# Test 1: Check Sun position
def test_sun_position_1():
    MJD = 51622 # J2000

    JPL_check = np.array([1.489004447920312E+08, -3.382674902835714E+06, 1.046597880064510E+02]) #given on this MJD
    JPL_check_unit = JPL_check / np.linalg.norm(JPL_check)

    py_unit = su.sun_position(MJD) / np.linalg.norm(su.sun_position(MJD))
    I_vec = np.array([1, 0, 0])
    ang1 = math.acos(np.dot(py_unit, I_vec))
    ang2 = math.acos(np.dot(JPL_check_unit, I_vec))

    np.testing.assert_allclose(ang1, ang2, atol=5e-5) # Python test, checking that the approximate angle is fairly close to JPL
    np.testing.assert_allclose(su.sun_position(MJD),sucpp.sun_position(MJD), atol=1e-6) # compare test

def test_sat_sun_vect_1():
    r = np.array([7000, 500, 1000])
    MJD = 51622
    check = su.sun_position(MJD) - r
    check = check / np.linalg.norm(check)

    np.testing.assert_allclose(su.sat_sun_vect(r, MJD), check, atol=1e-6)
    np.testing.assert_allclose(sucpp.sat_sun_vect(r, MJD), check, atol=1e-6)
    np.testing.assert_allclose(sucpp.sat_sun_vect(r, MJD), su.sat_sun_vect(r, MJD), atol=1e-6)

