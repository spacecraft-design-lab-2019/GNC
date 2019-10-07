import sys, math, numpy as np
import matplotlib.pyplot as plt

pi = math.pi

# inertia properties (add real later)
Ixx = 75
Iyy = 100
Izz = 125
I = np.array([[Ixx, 0, 0],[0, Iyy, 0], [0, 0, Izz]])

# deg --> rad, initial conditions
wx0 = 0.1*pi/180
wy0 = 0.1*pi/180
wz0 = pi/180

# moments, Nm
Mx0 = 0
My0 = 0
Mz0 = 0

# initial quaternion: {0,0,0,1} scalar last
x = []
x = np.array([[wx0],[wy0],[wz0],[Mx0],[My0],[Mz0],[0],[0],[0],[1]])

# state dot equations
def x_dot(x,I):
    x_dot = []
    wx_dot = (I[1][1]-I[2][2])/I[0][0]*x[1]*x[2]+x[3]/I[0][0]
#    print(x_dot)
    wy_dot = (I[2][2]-I[0][0])/I[1][1]*x[0]*x[2]+x[4]/I[1][1]
    wz_dot = (I[0][0]-I[1][1])/I[2][2]*x[0]*x[1]+x[5]/I[2][2]
    # moment derivative
    Mx_dot = np.array([0])
    My_dot = np.array([0])
    Mz_dot = np.array([0])
    # quaternion derivative
    q1_dot = x[0]*x[9]/2-x[1]*x[8]/2+x[2]*x[7]/2
    q2_dot = x[0]*x[8]/2+x[1]*x[9]/2-x[2]*x[6]/2
    q3_dot = x[1]*x[6]/2-x[0]*x[7]/2+x[2]*x[9]/2
    q4_dot = -x[0]*x[6]/2-x[1]*x[7]/2-x[2]*x[8]/2
    x_dot = np.array([wx_dot,wy_dot,wz_dot,Mx_dot,My_dot,Mz_dot,q1_dot,q2_dot,q3_dot,q4_dot])
    return x_dot

def euler_integrate(x,dt,tf):
    n = int(tf/dt)
    t = np.array([0])
    for i in range(n):
        # euler integrator
        y = np.reshape(x[:,-1],(10,1)) + np.reshape(dt*x_dot(x[:,-1],I),(10,1))
        # normalize the quaternion
        y[6:10] = y[6:10]/np.linalg.norm(y[6:10])
        x = np.append(x,y,axis=1)
        t = np.append(t,dt*i)
    return t,x

tf = 10000
dt = 1
n = int(tf/dt)

[t,y] = euler_integrate(x, dt, tf)
# y = [w, M, q]

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
