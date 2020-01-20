"""
 * Author: Nick Goodson 
 * Jan 18th 2020
 *
 * This code is a template, different sections can be pulled out and replaced
 * with calls to C++ functions to test the C++ functionality
 *
"""

import pdb
import os
import sys
import inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
currentdir = os.path.dirname(os.path.abspath(filename))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)

sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)

# import CPP modules here
import numpy as np
import matplotlib.pyplot as plt


# Test the algorithm
def main():

	N = 250
	Nx = 2
	Nu = 1

	Qf = np.eye(Nx) * 30
	Q = np.eye(Nx) * 0.01
	R = 0.3

	x0 = np.zeros((2))
	xg = np.array([np.pi, 0])
	utraj0 = np.zeros((Nu, N))

	dt = 0.01
	tol = 0.001

	# TODO: Replace this with call to C++ version
	xtraj, utraj, K, Jhist = iLQRsimple_py(x0, xg,utraj0, Q, R, Qf, dt, tol) 

	# Plot results
	fig, ax = plt.subplots(2, 2, 1)
	fig.title("iLQR results")
	ax[0, 0].plot(xtraj[0, :])
	ax[0, 0].set_title("angle")
	ax[0, 1].plot(xtraj[1, :])
	ax[0, 1].set_title("angular rate")
	ax[1, 0].plot(utraj)
	ax[1, 0].set_title("Control")
	ax[1, 1].plot(Jhist)


def iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol):
	"""
	Simple iLQR for testing C++ functions
	"""
	Nx = x0.shape[0]
	Nu = utraj0.shape[0]
	N = utraj0.shape[1] - 1

	A = np.zeros((Nx, Nx, N))
	B = np.zeros((Nx, Nu, N))

	J = 0
	Jhist = []
	xtraj = np.zeros((Nx, N+1))
	xtraj[:, 0] = x0
	utraj = utraj0


	# Forward simulate using initial controls
	"""TODO: Test with C++ 'rkstep()' function"""
	for k in range(0, N-1):
		J = J + 0.5 * (xtraj[:, k] - xg) @ Q * (xtraj[:, k] - xg) + 0.5 * R*utraj[:, k].T @ utraj[:, k]
		xtraj[:, k+1], A[:, :, k+1], B[:, :, k+1] = rkstep_py(xtraj[:, k], utraj[:, k], dt)

	J = J + 0.5 * (xtraj[:, N] - xg).T @ Qf * (xtraj[:, N] - xg)
	Jhist.append(J)

	S = np.zeros((Nx, Nx, N+1))
	s = np.zeros((Nx, N+1))
	K = np.zeros((Nu, Nx, N))
	l = (tol+1)*np.ones((Nu, N))

	# Backward pass
	"""TODO: Test replacement of entire while-loop with C++"""
	count = 0
	while (np.max(np.abs(l)) > tol):

		count += 1

		S[:, :, N] = Qf
		s[:, N] = Qf @ (xtraj[:, N] - xg)

		for k in range(N-1, -1, -1):
		
			# Calculate cost gradients for this time step
			q = Q @ (xtraj[:, k] - xg)
			r = R * utraj[:, k]
		
			# Calculate l and K
			LH = (R + B[:,:,k].T @ S[:,:,k+1] @ B[:,:,k])
			l[:,k] = np.linalg.solve(LH, (r + B[:,:,k].T @ s[:,k+1]))
			K[:,:,k] = np.linalg.solve(LH, B[:,:,k].T @ S[:,:,k+1] @ A[:,:,k])
		
			# Calculate new S and s        
			S[:,:,k] = Q + R*K[:,:,k].T @ K[:,:,k] + (A[:,:,k]-B[:,:,k] @ K[:,:,k]).T@S[:,:,k+1]@(A[:,:,k]-B[:,:,k]@K[:,:,k])
			s[:,k] = q - K[:,:,k].T*r + R*K[:,:,k].T@l[:,k] + (A[:,:,k]-B[:,:,k]@K[:,:,k]).T@(s[:,k+1] - S[:,:,k+1]@B[:,:,k]@l[:,k])

		# Forward pass line search with new l and K
		unew = np.zeros((Nu,N))
		xnew = np.zeros((Nx,N+1))
		xnew[:,1] = xtraj[:,1]
		alpha = 1.0
		Jnew = J+1
		while (Jnew > J):
			Jnew = 0
			for k in range(0, N-1):
				unew[:,k] = utraj[:,k] - alpha*l[:,k] - K[:,:,k]@(xnew[:,k]-xtraj[:,k])
				xnew[:,k+1], A[:,:,k], B[:,:,k] = rkstep_py(xnew[:, k], unew[:, k], dt)
				Jnew = Jnew + 0.5*(xnew[:,k]-xg).T@Q@(xnew[:,k]-xg) + 0.5*R*unew[:,k].T@unew[:,k]

			Jnew = Jnew + 0.5*(xnew[:,N]-xg).T@Qf@(xnew[:,N]-xg)
			alpha = 0.5*alpha   

		xtraj = xnew
		utraj = unew
		J = Jnew
		Jhist.append(J)

		print("Final l = {}".format(np.max(np.abs(l))), "alpha = {}".format(2*alpha))
		print("count = {}".format(count))

	return xtraj, utraj, K, Jhist


def rkstep_py(x0, u0, dt):

	# Define constants
	Nx = 2
	Nu = 1

	xdot1, dxdot1 = pendulumDynamics_py(0, x0, u0)
	xdot2, dxdot2 = pendulumDynamics_py(0, x0 + 0.5*xdot1*dt, u0)

	x1 = x0 + dt * xdot2

	A1 = dxdot1[:, 0:Nx]
	A2 = dxdot2[:, 0:Nx]
	B1 = dxdot1[:, Nx:]
	B2 = dxdot2[:, Nx:]

	A = np.eye(2) + dt*A2 + 0.5*dt*dt*A2@A1
	B = dt*B2 + 0.5*dt*dt*A2@B1

	return x1, A, B


def pendulumDynamics_py(t, x, u):

	Nx = 2

	m = 1.0
	l = .5
	b = 0.1
	lc = 0.5
	I = 0.25
	g = 9.81

	xdot = np.zeros((Nx))
	xdot[0] = x[1]
	xdot[1] = (u - m*g*lc*np.sin(x[0]) - b*x[1])/I

	A = np.array([[0, 1], [-m*g*lc*np.cos(x[0])/I, -b/I]])
	B = np.array([[0], [1/I]])

	dxdot = np.hstack((A, B))

	return xdot, dxdot


if __name__ == "__main__":
	main()

