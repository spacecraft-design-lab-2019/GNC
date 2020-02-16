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


#NOTE: All Quat rot functions test Lq and Rq as well as they are used. This also tests the skew matrix (hat) function. 
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

def test_quat_deriv():
	# Test using rotation and rotation rate solely about z axis
	theta = pi / 4.0
	d_theta = .0000001
	dt = .00001

	# Numerical derivative
	q_1 = np.array([cos(theta / 2), 0.0, 0.0, sin(theta / 2)])  # 45 degree rotation about the z axis body to inertial
	q_2 = np.array([cos((theta + d_theta)/2), 0.0, 0.0, sin((theta + d_theta)/2)])
	q_dot_num = (q_2-q_1)/dt

	# Analytical derivative
	w = np.array([0.0, 0.0, d_theta/dt])
	q_dot_anal = ecpp.get_q_dot(q_2, w)

	np.testing.assert_allclose(q_dot_num, q_dot_anal, atol=1e-5)

def test_get_w_dot():
	# Test checking Euler Equations
	I = np.array([[5, 0, 0], [0, 10, 0], [0, 0, 100]])
	M = np.array([2, -1, 5])
	w = np.array([2, 5, -10])

	wdot_pred = np.array([4502/5, -1901/10, -45/100])

	np.testing.assert_allclose(ecpp.get_w_dot(w, M, I), wdot_pred, 1e-15)

def test_get_q_dot():
	# Test checking Kinetmatic Equations
	q = np.array([cos(pi/4), 0, 0, sin(pi/4)]) # 90 degree rotation about z from body to ECI
	w = np.array([2, 5, -10])

	qdot_pred = 0.5 * np.array([10*sqrt(2)/2, sqrt(2) - 5 * sqrt(2)/2, sqrt(2) + 5*sqrt(2)/2, -10*sqrt(2)/2])

	np.testing.assert_allclose(ecpp.get_q_dot(q, w), qdot_pred, 1e-15)

def test_get_attitude_derivative():
	# Test checking Kinetmatic Equations

	I = np.array([[5, 0, 0], [0, 10, 0], [0, 0, 100]])
	M = np.array([2, -1, 5])
	x = np.array([cos(pi/4), 0, 0, sin(pi/4), 2, 5, -10])
	time = 0		# Needs time in order to work with scipy ODE integrator
	x_pred = np.array([10*sqrt(2)/4, sqrt(2)/2 - 5 * sqrt(2)/4, sqrt(2)/2 + 5*sqrt(2)/4, -10*sqrt(2)/4, 4502/5, -1901/10, -45/100])
	np.testing.assert_allclose(ecpp.get_attitude_derivative(time, x, M, I), x_pred, 1e-15)
