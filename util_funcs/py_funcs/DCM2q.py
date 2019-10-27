import math
import numpy as np
from numpy import linalg as LA


def trace(A):
    trace = 0
    for ii in range(len(A)):
        trace += A[ii,ii]
    return trace

def DCM2quat(A):
    ''' 
    A should be a 3x3 array
    output quaternion first
    '''
    B4 = math.sqrt(1/4*(1+trace(A)))
    B1 = math.sqrt(1/4*(1+2*A[0,0]-trace(A)));
    B2 = math.sqrt(1/4*(1+2*A[1,1]-trace(A)));
    B3 = math.sqrt(1/4*(1+2*A[2,2]-trace(A)));
    B = np.array([B1,B2,B3,B4])
    if B1 == max(B):
        q1 = B1;
        q2 = (A[0,1]+A[1,0])/(4*B1);
        q3 = (A[0,2]+A[2,0])/(4*B1);
        qs = (A[1,2]-A[2,1])/(4*B1);
    elif B2 == max(B):
        q1 = (A[0,1]+A[1,0])/(4*B2);
        q2 = B2;
        q3 = (A[1,2]+A[2,1])/(4*B2);
        qs = (A[2,0]-A[0,2])/(4*B2);
    elif B3 == max(B):
        q1 = (A[2,0]+A[0,2])/(4*B3);
        q2 = (A[2,1]+A[1,2])/(4*B3);
        q3 = B3;
        qs = (A[0,1]-A[1,0])/(4*B3);
    else:
        q1 = (A[1,2]-A[2,1])/(4*B4);
        q2 = (A[2,0]-A[0,2])/(4*B4);
        q3 = (A[0,1]-A[1,0])/(4*B4);
        qs = B4;

    q = np.array([qs,q1,q2,q3])
    q = q/LA.norm(q)

    return q
