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