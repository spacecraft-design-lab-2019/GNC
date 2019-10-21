'''
Script integrating detumble with orbit/magnetic field knowledge
'''

from euler import quat2DCM, get_attitude_derivative, get_q_dot, get_w_dot
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import numpy as np
import scipy.integrate as integrate
from orbit_propagation import get_orbit_pos

# clear figures
plt.close('all')
        
pi = math.pi

# inertia properties (add real later)
Ixx = 75
Iyy = 100
Izz = 125
I = np.array([[Ixx, 0.0, 0.0],[0.0, Iyy, 0.0], [0.0, 0.0, Izz]])

# initial attitude conditions, radians & rad/s
q_0 = np.array([[0.0],[0.0],[0.0],[1.0]])                     # initial quaternion, scalar last
w_0 = np.array([[0.1*pi/180.],[0.1*pi/180.],[pi/180.]])  # initial rotation rate

# initial orbit state conditions, TLE+epoch
epoch = '2013-12-14T14:18:37.00'
line1 = ('1 25635U 99008B   13348.59627062  .00000650  00000-0  16622-3 0  9860')
line2 = ('2 25635 096.4421 173.2395 0141189 010.0389 029.8678 14.46831495780970')
TLE = {'line1': line1, 'line2': line2}

mean_motion = 14.46/(24*3600)*2*math.pi # mean motion, radians/second
period = 2*pi/mean_motion                      # Period, seconds

# feed in a vector of times and plot orbit
t0 = 0
tf = period
n = 1000
times = np.linspace(t0,period,n)

# preallocate position storage matrix
positions = np.zeros((n,3))

for i in range(len(times)):
    positions[i,:] = get_orbit_pos(TLE, epoch, times[i])

