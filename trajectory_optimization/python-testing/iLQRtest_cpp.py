import os
import sys
import inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
currentdir = os.path.dirname(os.path.abspath(filename))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)

sys.path.insert(0, parentdir)
sys.path.insert(0, gncdir)

# import CPP module
import iLQRsimple_cpp as ilqr
import numpy as np
import matplotlib.pyplot as plt


# Test the algorithm
def main():

	# sizes
	N = 250
	Nx = 2
	Nu = 1

	# Cost matrices
	Qf = np.eye(Nx) * 30
	Q = np.eye(Nx) * 0.01
	R = np.array([[0.3]])

	xtraj = np.zeros((Nx, N))
	utraj = np.zeros((Nu, N-1))
	K = np.zeros((Nu, Nx*(N-1)))
	Jhist = []

	# Initial and final conditions
	x0 = np.zeros((2))
	xg = np.array([np.pi, 0])
	xtraj[:, 0] = x0

	dt = 0.01
	tol = 0.001

	# Call the C++ version of iLQRsimple
	result = ilqr.iLQRsimple(xg, Q, R, Qf, dt, tol, xtraj, utraj, K, Jhist)


	# xtraj, utraj, K, Jhist = iLQRsimple_py(x0, xg,utraj0, Q, R, Qf, dt, tol) # python version

	# Plot results
	# fig, ax = plt.subplots(2, 2, 1)
	# fig.title("iLQR results")
	# ax[0, 0].plot(xtraj[0, :])
	# ax[0, 0].set_title("angle")
	# ax[0, 1].plot(xtraj[1, :])
	# ax[0, 1].set_title("angular rate")
	# ax[1, 0].plot(utraj)
	# ax[1, 0].set_title("Control")
	# ax[1, 1].plot(Jhist)


if __name__ == "__main__":
	main()