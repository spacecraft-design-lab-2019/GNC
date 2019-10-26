import numpy as np
import math
from predict import predict
from measurement import measurement
from innovation import innovation
from numpy import linalg as LA
import sys
sys.path.append('../TRIAD/py_funcs/')
import deterministic_ad as triad

# initial conditions for MEKF
xk = np.array([0,0,0,1,0.1,0.1,0.1])
w = np.array([0.15,0.5,0.1])
dt = 0.01
W = 3.04617419786709e-10 * np.eye(6)
V = 3.04617419786709e-10 * np.eye(6)
# 10deg and 10deg/s 1 sigma uncertainty 
P = (10*math.pi/180)**2*np.eye(6) 
rN1 = np.array([1,0,0]) # for now...
rN2 = np.array([0,1,0]) # for now...
rN  = np.append(rN1,rN2)
rB1 = np.array([1,0,0]) # for now...
rB2 = np.array([0,1,0]) # for now...
rB  = np.append(rB1,rB2)
q   = np.array([1,0,0,0]) # for now...


# run predict step
'''
xn - predicted state (mu_k+1|k)
P  - predicted covariance (sigma_k+1|k)
A  - linearized state transition matrix
'''
xn,A = predict(xk,w,dt)
P = A*P*np.transpose(A)+W

# run measurement step
'''
'''
y,R,C = measurement(q,rN)

# run innovation step to find z
z,S = innovation(R, rN, rB, P, V, C)

# Kalman Gain
L = P @ C.T @ LA.inv(S)



