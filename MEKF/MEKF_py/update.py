import numpy as np
import math


def quatMult(q1,q2):
    qn = np.zeros(4)
    # propagate dynamics forward
    qn[0] = q1[0]*q2[0]-np.transpose(q1[1:4])*q2[1:4]
    qn[1:4] = q1[0]*q2[1:4]+q2[0]*q1[1:4]+np.cross(q1[1:4],q2[1:4])
    return qn

def update(L, z, xn, Pn, V, C):
    
    dx = L*z
    dphi = dx[0:3]
    dbeta = dx[3:6]
    temp_q = np.array([[math.sqrt(1-np.transpose(dphi)*dphi)],[dphi]])
    
    q = quatMult(xn[0:4],temp_q)
    b = xn[4:7] + dbeta
    
    x = np.array([[q],[b]])
    P = (np.identity(3)-L*C)*Pn*np.transpose(np.identity(3)-L*C)+L*V*np.transpose(L)
    
    return x, P