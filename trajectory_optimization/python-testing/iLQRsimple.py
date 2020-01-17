# Author: Nick Goodson
# 16th Jan 2020


import numpy as np


def iLQRsimple(x0, utraj0, xg, Q, R, Qf, dt, tol):
	"""
	Simple iLQR for testing C++ functions
	"""
	Nx = x0.shape(0)
	Nu = utraj0.shape(0)
	N = utraj0.shape(1) + 1

	A = np.zeros(Nx, Nx, N)
	B = np.zeros(Nx, Nu, N)

	J = 0
	xtraj = np.zeros(Nx, N)
	xtraj[:, 0] = x0

	# Forward simulate using baseline controls
	for k in range(0, N):
		J = J + 0.5 * (xtraj[:, k] - xg).T @ Q @ (xtraj[:, k] - xg) + 0.5 * (utraj[:, k]).T @ R @ utraj[:, k]
		xtraj[:, k+1], A[: ,: , k+1], B[:, :, k+1] = rkstep(xtraj[:, k], utraj[:, k], dt)

