import numpy as np
from numpy import linalg as LA
import math
def predict(xk,w,dt):
    '''
    Predict step of MEKF - state x = [q, Beta]
    qk*[cos(theta),rsin(theta)] uses Hamilton's quaternion mult.
    ** scalar first quaternion ** 
    '''
    q = xk[0:4]
    b = xk[4:7]
    
    theta = LA.norm(w-b)*dt
    r = (w-b)/LA.norm(w-b)
    
    xn = np.zeros(7)
    
    xn[0] = q[0]*math.cos(theta)-q[1]*r[0]*math.sin(theta)-\
            q[2]*r[1]*math.sin(theta)-q[3]*r[2]*math.sin(theta)
    xn[1] = q[0]*r[0]*math.sin(theta)+q[1]*math.cos(theta)\
            +q[2]*r[2]*math.sin(theta)-q[3]*r[1]*math.sin(theta)
    xn[2] = q[0]*r[1]*math.sin(theta)-q[1]*r[2]*math.sin(theta)\
            +q[2]*math.cos(theta)+q[3]*r[0]*math.sin(theta)
    xn[3] = q[0]*r[2]*math.sin(theta)+q[1]*r[1]*math.sin(theta)\
            -q[2]*r[0]*math.sin(theta)+q[3]*math.cos(theta)
    
    xn[4:7] = b[4:7]
    
    A = np.zeros(6,6)
    # fill out A matrix
    for ii in range(3):
        print(ii)
        A[ii,ii] = math.exp(-r[ii]*dt)
        A[ii,ii+3] = -dt
        A[ii+3,ii+3] = 1
        
    return xn, A


