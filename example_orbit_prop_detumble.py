'''
Script integrating detumble with orbit/magnetic field knowledge
'''

from euler import quat2DCM, get_attitude_derivative, get_q_dot, get_w_dot
from detumble.py_funcs import detumble_B_cross,detumble_B_dot,get_B_dot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import scipy.integrate as integrate
from orbit_propagation import get_orbit_pos, get_B_field_at_point
# from GNC.cmake_build_debug import SGP4_cpp as SGP4
# from util_funcs.py_funcs.frame_conversions import eci2ecef
import time_functions_cpp as tfcpp
import frame_conversions_cpp as fccpp
import time

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
w_0 = np.array([[.3],[.3],[.3]])  # initial rotation rate, rad/s
# initial state: quaternion, rotation rate
x_0 = np.squeeze(np.concatenate((q_0,w_0)))

# initial orbit state conditions, TLE+epoch
epoch = '2013-12-14T14:18:37.00'
line1 = ('1 25635U 99008B   13348.59627062  .00000650  00000-0  16622-3 0  9860')
line2 = ('2 25635  96.4421 173.2395 0141189  10.0389  29.8678 14.46831495780970')
TLE = {'line1': line1, 'line2': line2}

# initial orbit time (fair warning, this is the time for PyCubed, not Orsted)
MJD = 58827.53750000009
GMST_0 = tfcpp.MJD2GMST(MJD)

mean_motion = 14.46/(24*3600)*2*math.pi # mean motion, radians/second
period = 2*pi/mean_motion                      # Period, seconds

# feed in a vector of times and plot orbit
t0 = 0.0
tf = 600.0
tstep = .1
times = np.arange(t0,tf,tstep)
n = len(times)

# preallocate position storage matrix
positions_ECI = np.zeros((n-1,3))
positions_ECEF = np.zeros((n-1,3))
B_field_body = np.zeros((n-1,3))
w_vec = np.zeros((n-1,3))


# Define function for calculating full state derivative
t = time.time()

# extract position info at all times (from Python)
x = x_0;
for i in range(len(times)-1):
    positions_ECI[i,:] = get_orbit_pos(TLE, epoch, times[i])

    # convert to ECEF
    R_ECI2ECEF = fccpp.eci2ecef(GMST_0 + times[i])

    positions_ECEF[i,:] = np.transpose(np.matmul(R_ECI2ECEF, np.transpose(positions_ECI[i,:])))
    lat, lon, alt = fccpp.ecef2lla(np.transpose(positions_ECEF[i,:]))
    R_ECEF2ENU = fccpp.ecef2enu(lat, lon)

    # get magnetic field at position
    B_field_ENU = get_B_field_at_point(positions_ECEF[i,:]) # North, East, Down


    # get magnetic field in body frame (for detumble algorithm)
    R_body2ECI = quat2DCM(np.transpose(x[0:4]))
    B_field_ECEF = np.matmul(np.transpose(R_ECEF2ENU),B_field_ENU)
    B_field_ECI = np.matmul(np.transpose(R_ECI2ECEF),B_field_ECEF)
    B_field_body[i,:] = np.transpose(np.matmul(np.transpose(R_body2ECI),B_field_ECI))

    # Get B_dot based on previous measurement
    if i>0:
        B_dot = get_B_dot(np.transpose(B_field_body[i-1,:]),np.transpose(B_field_body[i,:]),tstep)
        Moment = detumble_B_dot(np.transpose(B_field_body[i,:]),B_dot,k=1.0)
        # Normalize and scale Moment
        k = 10.0
        M = k * Moment / np.linalg.norm(Moment)
    else:
        M = np.zeros((3,1))



    # Use magnetic field and angular rate to find commanded momentum
    # M = detumble_B_cross(np.transpose(x[3:6]),np.transpose(B_field_body),k=1)

    # Propagate dynamics/kinematics forward using commanded moment
    y = integrate.odeint(get_attitude_derivative, x, (times[i],times[i+1]), (M, I), tfirst=True)

    # Store angular velocity
    w_vec[i,:] = y[-1,4:7]

    # Update full attitude state
    x = y[-1,:]

elapsed = time.time() - t
print(elapsed)

# need tfirst = true for t,y ordered inputs. Include parameters/extra arguments as tuple.
#

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
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D(positions_ECI[:,0],positions_ECI[:,1],positions_ECI[:,2])
# ax.set_title('Orsted orbit')
# # Esoteric plotting function
# with plt.rc_context(rc={'interactive': False}):
#     plt.show()
# plt.show(block = True)

# plot angular velocity over time
fig2 = plt.figure()
plt.plot(w_vec[:,0])
plt.plot(w_vec[:,1])
plt.plot(w_vec[:,2])

# plot B field components as a function of time



# Let's 
