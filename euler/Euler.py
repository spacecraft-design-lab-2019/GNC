import sys, math, numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

'''
TODO:
    - write quaternion multiplication function
    - write quaternion difference function
    - enforce quaternion normalization --> enforced now in q_dot function-Paul
'''

# rate of change of angular rate
def get_w_dot(w,M,I):
    '''
    Takes in angular rate, net torque (3x1, principal frame), and the
    principal moment of inertia matrix (3x3, principal frame). Returns the rates of change
    of the angular rates as a 3x1.
    Inputs:
        w - angular velocity vector, 3x1, principal frame, rad/s
        M - vector of moments, 3x1, principal frame, N-m
        I - principal moment of inertia matrix, 3x3, kg-m^2
    Outputs:
        w_dot - angular acceleration vector, 3x1, principal frame, rad/s^2
    '''
    wx_dot = (I[1][1]-I[2][2])/I[0][0]*w[1]*w[2]+M[0]/I[0][0]
    wy_dot = (I[2][2]-I[0][0])/I[1][1]*w[0]*w[2]+M[1]/I[1][1]
    wz_dot = (I[0][0]-I[1][1])/I[2][2]*w[0]*w[1]+M[2]/I[2][2]
    return np.squeeze(np.array([[wx_dot],[wy_dot],[wz_dot]]))

def get_q_dot(q,w):
    '''
    Takes in a quaternion and the rotation rate vector,
    and returns the time derivative of the quaternion 
    Inputs:
        q -  4x1, scalar last, normalized
        w -  3x1, rad/s
    Outputs:
        q_dot - 4x1, scalar last, 1/sec
    '''
    
    # Enforce normalized quaternion before finding derivative
    q_hat = q/np.linalg.norm(q)
    
    q1_dot =  w[0]*q_hat[3]/2-w[1]*q_hat[2]/2+w[2]*q_hat[1]/2
    q2_dot =  w[0]*q_hat[2]/2+w[1]*q_hat[3]/2-w[2]*q_hat[0]/2
    q3_dot =  w[1]*q_hat[0]/2-w[0]*q_hat[1]/2+w[2]*q_hat[3]/2
    q4_dot = -w[0]*q_hat[0]/2-w[1]*q_hat[1]/2-w[2]*q_hat[2]/2
    return np.squeeze(np.array([[q1_dot],[q2_dot],[q3_dot],[q4_dot]]))
    
# state dot equations
def get_attitude_derivative(t,x,M,I):  
    '''
    Takes in an attitude state,
    a set of moments, and the spacecraft moment of inertia matrix and outputs the derivative
    of that attitude state
    Inputs:
        x - spacecraft attitude state, [q;w], 7x1, principal frame
        M - vector of moments, 3x1, principal frame, N-m
        I - principal moment of inertia matrix, 3x3, kg-m^2
    Outputs:
        x_dot - [q_dot;w_dot], 7x1
    '''
    # rotational rate derivative
    w_dot = get_w_dot(x[4:7],M,I)
    
    # quaternion derivative
    q_dot = get_q_dot(x[0:4],x[4:7])

    return np.concatenate((q_dot,w_dot))

def quat2DCM(quat):
    '''
    Takes in a quaternion converts to a direction cosine matrix
    Inputs:
        q - quaternion, scalar first, (4x1)
        !!!!!!!! Why is this scalar first, but q_dot function uses scalar last!!!!!!-Paul
    Outputs:
        DCM - 3x3 direction cosine matrix
    TODO:
        - clarify quaternion convention
        - clarify quaternion rotation direction convention (eci2attitude, attitude2eci?)
        - clarify DCM rotation direction convention
    '''
    q1 = quat[0]
    q2 = quat[1]
    q3 = quat[2]
    q4 = quat[3]
    q_vec = np.array([[q1],[q2],[q3]])
    # calculate skew symmetric matrix Q
    Q = np.array([[0,-q3,q2],[q3,0,-q1],[-q2,q1,0]])
    # calculate DCM
    DCM = (q4**2-np.transpose(q_vec)*q_vec)*np.eye(3)-2*q4*Q+2*q_vec*np.transpose(q_vec)
    return DCM

