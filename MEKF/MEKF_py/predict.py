import numpy as np
from numpy import linalg as LA
import math

def SkewSymmetric(a):
    a1 = a[0]
    a2 = a[1]
    a3 = a[2]
    a_skew = np.array([[0,-a3,a2],[a3,0,-a1],[-a2,a1,0]])
    return a_skew

def predict(xk,w,dt):
    '''
    Predict step of MEKF - state x = [q, Beta]
    qk*[cos(theta),rsin(theta)] uses Hamilton's quaternion mult.
    ** scalar first quaternion ** 


    '''
    q = xk[0:4]
    b = xk[4:7]
    
    u = w + b
    theta = LA.norm(u-b)*dt
    r = (u-b)/LA.norm(u-b)
    
    xn = np.zeros(7)
    
    # q2 == s
    s = np.array([math.cos(theta/2),r*math.sin(theta)])
    
    # propagate dynamics forward
    xn[0] = q[0]*s[0]-np.transpose(q[1:4])*s[1:4]
    xn[1:4] = q[0]*s[1:4]+s[0]*q[1:4]+np.cross(q[1:4],s[1:4])
    xn[4:7] = b
    
    # isolate vector component of quaternion
    V = np.array([[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
    # left quaternion multiply
    Ls = np.zeros(4,4)
    Ls[0,0] = s[0]
    Ls[0,1:4] = -np.transpose(s[1:4])
    Ls[1:4,0] = s[1:4]
    Ls[1:4,1:4] = s[0]*np.identity(3) + SkewSymmetric(s[1:4])
    
    # right quaternion multiply
    Rs = np.zeros(4,4)
    Rs[0,0] = s[0]
    Rs[0,1:4] = -np.transpose(s[1:4])
    Rs[1:4,0] = s[1:4]
    Rs[1:4,1:4] = s[0]*np.identity(3) - SkewSymmetric(s[1:4])
    
    
    A = np.zeros([6,6])
    # fill out A matrix
    for ii in range(3):
        A[ii,ii] = V*np.transpose(Ls)*Rs*np.transpose(V)
        A[ii,ii+3] = -0.5*dt
        A[ii+3,ii+3] = 1
        
    return xn, A



