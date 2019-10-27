import numpy as np
import math
from predict import predict
from measurement import measurement
from innovation import innovation
from update import update

from numpy import linalg as LA
import sys
sys.path.append('../../TRIAD/py_funcs/')
import deterministic_ad
sys.path.append('../../util_funcs/py_funcs/')
import DCM2q 

# import data from sample mat file
from scipy.io import loadmat
MEKF_inputs = loadmat('mekf_inputs.mat')

rN1 = MEKF_inputs['rN1']
rN2 = MEKF_inputs['rN2']
rN = np.vstack((rN1,rN2))

whist = MEKF_inputs['whist']

rB2hist = MEKF_inputs['rB2hist']
rB1hist = MEKF_inputs['rB1hist']
# initial conditions:
rB1 = np.reshape(rB1hist[:,0],(3,1))
rB2 = np.reshape(rB2hist[:,0],(3,1))

# 10deg and 10deg/s 1 sigma uncertainty 
W = MEKF_inputs['W']
V = MEKF_inputs['V']
dt = MEKF_inputs['dt']


Pk = (10*math.pi/180)**2*np.eye(6) 

R0   = deterministic_ad.triad_ad(np.hstack((rB1,rB2)),np.hstack((rN1,rN2)))
q0    = DCM2q.DCM2quat(R0)  
q0 = np.reshape(q0,(4,1))

# initial conditions for MEKF
Beta0 = np.reshape(np.array([0.1,0.1,0.1]),(3,1))
xk = np.vstack((q0,Beta0))


# start sim
for ii in range(len(rB1hist.T)):
    # run predict step
    '''
    xn - predicted state (mu_k+1|k)
    Pn - predicted covariance (sigma_k+1|k)
    A  - linearized state transition matrix
    W  - Process noise covariance
    V  - Measurement noise covariance
    rN - Vector measurement in newtonian frame
    rB - Vector measurement in body frame
    '''
    # get measurement vector
    rB1 = np.reshape(rB1hist[:,ii],(3,1))
    rB2 = np.reshape(rB2hist[:,ii],(3,1))
    rB = np.vstack((rB1,rB2))
    
    xn,A = predict(xk,whist[:,ii],dt)
    Pn = A*Pk*np.transpose(A)+W
    
    # run measurement step
    y,R,C = measurement(xn[0:4],rN)
    
    # run innovation step to find z
    z,S = innovation(R, rN, rB, Pn, V, C)
    
    # Kalman Gain
    L = Pn @ np.transpose(C) @ LA.inv(S)
    
    # update step
    xk, Pk = update(L, z, xn, Pn, V, C)
    
    
    

