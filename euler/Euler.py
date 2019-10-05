import sys, math, numpy as np
from scipy.integrate import odeint
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
        y = np.reshape(x[:,-1],(10,1)) + np.reshape(dt*x_dot(x[:,-1],I),(10,1))
        x = np.append(x,y,axis=1)
        t = np.append(t,dt*i)
    return t,x

tf = 10000
dt = 1
n = int(tf/dt)
[t,y] = euler_integrate(x, dt, tf)

#plt.plot(t,y[0,:])
