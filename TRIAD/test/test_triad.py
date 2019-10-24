from GNC.util_funcs.py_funcs import time_functions, frame_conversions, sun_utils
from GNC.TRIAD.py_funcs import deterministic_ad
import numpy as np
import pytest
import math
from GNC.cmake_build_debug import sun_utils_cpp, time_functions_cpp, frame_conversions_cpp, triad_cpp
import pyIGRF
from scipy.spatial.transform import Rotation as R


def test_triad_python():
    rot = R.from_quat([0, 0, np.sin(np.pi / 4), np.cos(np.pi / 4)])
    DCM_body2eci = rot.as_dcm()
    DCM_eci2body = np.transpose(DCM_body2eci)
    r_ecef = np.array([7000, 1000, -1000])
    mon = 10
    D = 8
    Y = 2019
    HH = 8
    MM = 10
    SS = 0
    # Converting Date to Times
    MJD = time_functions.date2MJD(mon, D, Y, HH, MM, SS)
    GMST = time_functions.MJD2GMST(MJD)


    # Getting Lat Long and then Magnetic Field
    (lat, lon, alt) = frame_conversions.ecef2lla(r_ecef)
    (_, _, _, B_n, B_e, B_d, _) = pyIGRF.igrf_value(lat, lon, alt, Y)
    B_enu = np.array([B_e, B_n, -B_d])
    B_enu_unit = B_enu / np.linalg.norm(B_enu)

    # Rotation Matrices
    R_ecef2enu = frame_conversions.ecef2enu(lat, lon)
    R_eci2ecef = frame_conversions.eci2ecef(GMST)

    # Converting to ECI and getting unit vector to the Sun
    r_eci = np.transpose(R_eci2ecef) @ r_ecef
    r_sat_sun = sun_utils.sat_sun_vect(r_eci, MJD)
    r_sat_sun_unit = r_sat_sun / np.linalg.norm(r_sat_sun)

    # Getting unit vector for Magnetic Field
    B_eci_unit = np.transpose(R_ecef2enu) @ np.transpose(R_eci2ecef) @ B_enu_unit
    V = np.column_stack((r_sat_sun_unit, B_eci_unit))

    # Spoofing measurements
    M = DCM_eci2body @ V
    # numpy.random.normal(0, 0.01, size=1)
    M[0, 0] = M[0, 0] + np.random.normal(0, 0.001, size=1)
    M[1, 0] = M[1, 0] + np.random.normal(0, 0.001, size=1)
    # M[:, 0] = M[:, 0] / np.linalg.norm(M[:,0])


    R_eci2body_py = deterministic_ad.triad_ad(M, V)
    R_eci2body_cpp = triad_cpp.triad_ad(M, V)
    np.testing.assert_allclose(R_eci2body_py, DCM_eci2body, atol=5e-3)  # Python test
    np.testing.assert_allclose(R_eci2body_cpp, DCM_eci2body, atol=5e-3)  # cpp test
    np.testing.assert_allclose(R_eci2body_cpp, R_eci2body_py, atol=5e-3)  # compare test
