# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:53:17 2019

@author: Paul DeTrempe

@description: example script for calling and testing dynamics/kinematics functions
"""

from euler import quat2DCM, get_attitude_derivative, get_q_dot, get_w_dot
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# clear figures
plt.close('all')
        
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

t0 = 0
tf = 10000
dt = 1
n = int(tf/dt)
times = np.linspace(t0,tf,n)
# need tfirst = true for t,y ordered inputs. Include parameters/extra arguments as tuple.
y = integrate.odeint(get_attitude_derivative,x_0,times,(M_0,I),tfirst=True)


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(y[:,4],y[:,5],y[:,6])
ax.set_title('herpolode')
ax.set_xlabel('omega_x')
ax.set_ylabel('omega_y')
ax.set_zlabel('omega_z')
