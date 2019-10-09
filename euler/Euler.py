import sys, math, numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

pi = math.pi

# inertia properties (add real later)
Ixx = 75
Iyy = 100
Izz = 125
I = np.array([[Ixx, 0.0, 0.0],[0.0, Iyy, 0.0], [0.0, 0.0, Izz]])

# initial conditions, radians & rad/s
q_0 = np.array([[0.0],[0.0],[0.0],[1.0]])                     # initial quaternion, scalar last
w_0 = np.array([[0.1*pi/180.],[0.1*pi/180.],[pi/180.]])  # initial rotation rate

# moments, Nm
M_0 = np.array([[0.0],[0.0],[0.0]])

# initial state: quaternion, rotation rate
x_0 = np.squeeze(np.concatenate((q_0,w_0)))

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
    q1_dot =  w[0]*q[3]/2-w[1]*q[2]/2+w[2]*q[1]/2
    q2_dot =  w[0]*q[2]/2+w[1]*q[3]/2-w[2]*q[0]/2
    q3_dot =  w[1]*q[0]/2-w[0]*q[1]/2+w[2]*q[3]/2
    q4_dot = -w[0]*q[0]/2-w[1]*q[1]/2-w[2]*q[2]/2
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



    
    

#def euler_integrate(x,dt,tf):
#    n = int(tf/dt)
#    t = np.array([0])
#    for i in range(n):
#        # euler integrator
#        y = np.reshape(x[:,-1],(10,1)) + np.reshape(dt*x_dot(x[:,-1],I),(10,1))
#        # normalize the quaternion
#        y[6:10] = y[6:10]/np.linalg.norm(y[6:10])
#        x = np.append(x,y,axis=1)
#        t = np.append(t,dt*i)
#    return t,x
t0 = 0
tf = 10000
dt = 1
n = int(tf/dt)
times = np.linspace(t0,tf,n)

#[t,y] = euler_integrate(x, dt, tf)
# y = [w, M, q]

integrated_solution = integrate.odeint(get_attitude_derivative,x_0,times,(M_0,I),tfirst=True)

w_out = y[0:3,:]
M_out = y[3:6,:]
q_out = y[6:10,:]

def quat2DCM(quat):
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

DCM_out = quat2DCM(q_out[:,0])

#plt.plot(t,y[0,:])
