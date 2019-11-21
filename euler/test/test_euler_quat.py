import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)



import numpy as np
import pytest
import math
import euler_cpp as ecpp
from math import cos, sin , pi, sqrt

def test_quat_rot1():
	theta = pi / 4
	q = np.array([cos(theta /2), 0, 0, sin(theta/2)]) # 45 degree rotation about the z axis body to inertial

	xb = np.array([1, 0, 0])
	zb = np.array([0, 0, 1])
	x_sol = np.array([cos(theta), sin(theta), 0])
	z_sol = np.array([0, 0, 1])

	np.testing.assert_allclose(ecpp.rotate_vec(xb, q), x_sol, atol=1e-15)
	np.testing.assert_allclose(ecpp.rotate_vec(zb, q), z_sol, atol=1e-15)


def test_quat_rot2():
	theta = pi / 2
	q = np.array([cos(theta /2), 0, 0, sin(theta/2)]) # 90 degree rotation about the z axis body to inertial

	xb = np.array([1, 0, 0])
	x_sol = np.array([cos(theta), sin(theta), 0])

	np.testing.assert_allclose(ecpp.rotate_vec(xb, q), x_sol, atol=1e-15)


def test_quat_rot3():
	theta = pi / 2
	q = np.array([cos(theta /2), sin(theta/2), 0, 0]) # 90 degree rotation about the x axis body to inertial

	xb = np.array([1, 0, 0])
	x_sol = np.array([1, 0, 0])

	np.testing.assert_allclose(ecpp.rotate_vec(xb, q), x_sol, atol=1e-15)

def test_quat_rot4():
	theta = pi / 3
	q = np.array([cos(theta /2), sin(theta/2), 0, 0]) # 60 degree rotation about the x axis body to inertial

	xb = np.array([0, 1, 0])
	x_sol = np.array([0, cos(theta), sin(theta)])

	np.testing.assert_allclose(ecpp.rotate_vec(xb, q), x_sol, atol=1e-15)

def test_quat_rot5():
	theta = pi / 4
	q = np.array([cos(theta /2), 0, 0, sin(theta/2)]) # 45 degree rotation about the z axis body to inertial

	zb = np.array([0, 0, 1])
	z_sol = np.array([0, 0, 1])

	np.testing.assert_allclose(ecpp.rotate_vec(zb, q), z_sol, atol=1e-15)

def test_quat_rot6():
	theta = pi / 4
	q = np.array([cos(theta /2), 0, 0, sin(theta/2)]) # 45 degree rotation about the z axis body to inertial

	yb = np.array([0, 1, 0])
	y_sol = np.array([-sin(theta), cos(theta), 0])

	np.testing.assert_allclose(ecpp.rotate_vec(yb, q), y_sol, atol=1e-15)

def test_quat_inverse():
	vec = np.array([.33,.66,.99])
	q = np.array([sqrt(.25),sqrt(.25),sqrt(.25),sqrt(.25)])

	vec_rotated = ecpp.rotate_vec(vec, q)
	q_inv = ecpp.get_inverse_quaternion(q)
	vec_rot_back = ecpp.rotate_vec(vec_rotated, q_inv)

	np.testing.assert_allclose(vec,vec_rot_back,atol = 1e-15)