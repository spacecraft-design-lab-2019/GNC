'''
Script integrating detumble with orbit/magnetic field knowledge
'''

from euler import quat2DCM, get_attitude_derivative, get_q_dot, get_w_dot
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import numpy as np
import scipy.integrate as integrate
from orbit_propagation import get_orbit_pos, get_B_field_at_point
from GNC.cmake_build_debug import SGP4_cpp as SGP4
from util_funcs.py_funcs.frame_conversions import eci2ecef

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
line2 = ('2 25635  96.4421 173.2395 0141189  10.0389  29.8678 14.46831495780970')
TLE = {'line1': line1, 'line2': line2}

mean_motion = 14.46/(24*3600)*2*math.pi # mean motion, radians/second
period = 2*pi/mean_motion                      # Period, seconds

# feed in a vector of times and plot orbit
t0 = 0
tf = 6000
n = 101
times = np.linspace(t0,tf,n)

# preallocate position storage matrix
positions_ECI = np.zeros((n,3))
positions_ECEF = np.zeros((n,3))
B_fields = np.zeros((n,3))

# extract position info at all times (from Python)
for i in range(len(times)):
    positions_ECI[i,:] = get_orbit_pos(TLE, epoch, times[i])
    # get magnetic field at all these positions
    # convert to ECEF
    B_fields[i,:] = get_B_field_at_point(positions[i,:])
    # North, East, Down

# extract position info from C++
# typerun     - type of run                    verification 'v', catalog 'c', manual 'm'                         
# typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
# opsmode     - mode of operation afspc or improved 'a', 'i'
#	whichconst  - which set of constants to use  72, 84

# get gravity constants first
# wgs72 = SGP4.get_gravconsttype(72)
#satrec = SGP4.twoline2rv_wrapper(line1, line2, 72)
#satrec_ptr = SGP4.get_new_satrec()

# plot trajectory
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(positions_ECI[:,0],positions_ECI[:,1],positions_ECI[:,2])
ax.set_title('Orsted orbit')

# Let's 
