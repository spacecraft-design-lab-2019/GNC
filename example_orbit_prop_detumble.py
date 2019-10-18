'''
Script integrating detumble with orbit/magnetic field knowledge
'''

from euler import quat2DCM, get_attitude_derivative, get_q_dot, get_w_dot
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import numpy as np
import scipy.integrate as integrate


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


