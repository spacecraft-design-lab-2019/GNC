import numpy as np
import math
from predict import predict
from numpy import linalg as LA
import sys
sys.path.append('../TRIAD/py_funcs/')
import deterministic_ad as triad

# initial conditions for MEKF
xk = np.array([0,0,0,1,0.1,0.1,0.1])
w = np.array([0.15,0.5,0.1])
dt = 0.01
W = 3.04617419786709e-10 * np.eye(6)
# 10deg and 10deg/s 1 sigma uncertainty 
P = (10*math.pi/180)**2*np.eye(6) 

# run predict step
'''
xn - predicted state (mu_k+1|k)
P  - predicted covariance (sigma_k+1|k)
A  - linearized state transition matrix
'''
xn,A = predict(xk,w,dt)
P = A*P*np.transpose(A)+W

# find rotation matrix using TRIAD
M = np.transpose(np.array([[1,0,0],[0,1,0]]))
V = np.transpose(np.array([[0,1,0],[1,0,0]]))
R = triad.triad_ad(M,V)

Q = np.zeros([6,6])
Q[0:3,0:3] = R
Q[3:6,3:6] = R

# change for real measurement at next step
rB1 = np.array([1,0,0])
rB2 = np.array([0,1,0])
rB = np.append(rB1,rB2)

rN1 = rB1 = np.array([1,0,0])
rN2 = np.array([0,1,0])
rN = np.append(rN1,rN2)

# innovation step
z = rB - np.matmul(Q,rN)

C = np.array([rB,np.zeros([7,1])])
