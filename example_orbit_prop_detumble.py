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
import detumble_cpp as dcpp
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
q_0 = np.array([[1.0],[0.0],[0.0],[0.0]])                     # initial quaternion, scalar last
w_0 = np.array([[.05],[.05],[.05]])  # initial rotation rate, rad/s
# initial state: quaternion, rotation rate
x_0 = np.squeeze(np.concatenate((q_0,w_0)))

# initial orbit state conditions, TLE+epoch
epoch = '2020-05-10T08:05:03.00'
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
tf = period
tstep = 2
times = np.arange(t0,tf,tstep)
n = len(times)

# preallocate position storage matrix
positions_ECI = np.zeros((n-1,3))
positions_ECEF = np.zeros((n-1,3))
B_field_body = np.zeros((n-1,3))
B_field_ECI_vec = np.zeros((n-1,3))
B_field_NED_vec = np.zeros((n-1,3))
B_dot_body = np.zeros((n-1,3))
w_vec = np.zeros((n-1,3))
q_vec = np.zeros((n-1,4))
M_vec = np.zeros((n-1,3))


# Define function for calculating full state derivative
t = time.time()

# extract position info at all times (from Python)
x = x_0;
for i in range(len(times)-1):
    # Get GMST at this time
    GMST = tfcpp.MJD2GMST(MJD + times[i] / 60.0 / 60.0 / 24.0)

    positions_ECI[i,:] = get_orbit_pos(TLE, epoch, times[i])

    # convert to ECEF
    R_ECI2ECEF = fccpp.eci2ecef(GMST)

    positions_ECEF[i,:] = np.transpose(R_ECI2ECEF @ np.transpose(positions_ECI[i,:]))
    lat, lon, alt = fccpp.ecef2lla(np.transpose(positions_ECEF[i,:]))
    R_ECEF2ENU = fccpp.ecef2enu(lat, lon)

    # get magnetic field at position
    B_field_NED = get_B_field_at_point(positions_ECEF[i,:]) # North, East, Down
    B_field_NED_vec[i,:] = np.transpose(B_field_NED)        # store for later analysis
    B_field_ENU = np.array([[B_field_NED[1]],[B_field_NED[0]],[-B_field_NED[2]]])    # north, east, down to east, north, up

    # get magnetic field in body frame (for detumble algorithm)
    R_body2ECI = quat2DCM(np.transpose(x[0:4]))
    B_field_ECEF = np.transpose(R_ECEF2ENU) @ B_field_ENU
    B_field_ECI = np.transpose(R_ECI2ECEF) @ B_field_ECEF

    # Correct to max expected value of Earth's magnetic field if pyIGRF throws a huge value
    # if i > 0:
        # if np.linalg.norm(B_field_ECI) > 7e-08:
        #     B_field_ECI = B_field_ECI/np.linalg.norm(B_field_ECI)*np.linalg.norm(B_field_ECI_vec[i-1,:])#* 7e-08


    B_field_ECI_vec[i,:] = np.transpose(B_field_ECI)
    B_field_body[i,:] = np.transpose( np.transpose(R_body2ECI) @ B_field_ECI)

    # Get B_dot based on previous measurement
    if i>0:
        B_dot = get_B_dot(np.transpose(B_field_body[i-1,:]),np.transpose(B_field_body[i,:]),tstep)
        B_dot_body[i,:] = np.transpose(B_dot)

        # # Validate B_dot algorithm
        # k_B_dot = -4.0*math.pi/period*(2)*Ixx*9.0e7
        # Moment = dcpp.detumble_B_dot(np.transpose(B_field_body[i,:]),B_dot, k_B_dot)
        # # Torque on spacecraft
        # M = np.cross(Moment, np.transpose(B_field_body[i, :]))

        # Validate B_cross_c++
        k_B_cross = 4.0*math.pi/period*(2)*Ixx*2e-1
        M = dcpp.detumble_B_cross(np.transpose(x[4:7]),np.transpose(B_field_body[i,:]),k_B_cross)

        # # Validate B_cross_python
        # k_B_cross = 4.0*math.pi/period*(2)*Ixx*1.0e-1
        # M = detumble_B_cross(np.transpose(x[4:7]), np.transpose(B_field_body[i, :]), k_B_cross)

    else:
        M = np.zeros((3,1))

    # store for plotting
    M_vec[i,:] = np.transpose(M)

    # Propagate dynamics/kinematics forward using commanded moment
    y = integrate.odeint(get_attitude_derivative, x, (times[i],times[i+1]), (M, I), tfirst=True)

    # Store angular velocity
    w_vec[i,:] = y[-1,4:7]
    q_vec[i,:] = y[-1,0:4]

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

# # plot trajectory
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D(positions_ECI[:,0],positions_ECI[:,1],positions_ECI[:,2])
# ax.set_title('Orbit, ECI')
# # # Esoteric plotting function
# with plt.rc_context(rc={'interactive': False}):
#     plt.show()
# #     plt.show(block = True)

# plot angular velocity over time
fig2 = plt.figure()
plt.plot(w_vec[:,0])
plt.plot(w_vec[:,1])
plt.plot(w_vec[:,2])
plt.title('Angular velocity components')


# # Plot trajectory in ECEF
# fig3 = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D(positions_ECEF[:,0],positions_ECEF[:,1],positions_ECEF[:,2])
# ax.set_title('Orbit, ECEF')
# # # Esoteric plotting function
# with plt.rc_context(rc={'interactive': False}):
#     plt.show()

# plot B field components as a function of time
fig4 = plt.figure()
plt.plot(B_field_body[:,0])
plt.plot(B_field_body[:,1])
plt.plot(B_field_body[:,2])
plt.title('Magnetic field components, B')
plt.show()

# plot B dotcomponents as a function of time
fig5 = plt.figure()
plt.plot(B_dot_body[:,0])
plt.plot(B_dot_body[:,1])
plt.plot(B_dot_body[:,2])
plt.title('Magnetic field rate of change, B_dot')
plt.show()

# plot B dot components as a function of time
fig5 = plt.figure()
plt.plot(B_field_ECI_vec[:,0])
plt.plot(B_field_ECI_vec[:,1])
plt.plot(B_field_ECI_vec[:,2])
plt.title('Magnetic field ECI')
plt.show()

# plot quaternion over time
fig6 = plt.figure()
plt.plot(q_vec[:,0])
plt.plot(q_vec[:,1])
plt.plot(q_vec[:,2])
plt.plot(q_vec[:,3])
plt.title('quaternion components')
plt.show()
# plot Moment over time
fig7 = plt.figure()
plt.plot(M_vec[:,0])
plt.plot(M_vec[:,1])
plt.plot(M_vec[:,2])
plt.title('Moment components')
plt.show()
# plot norm of velocity vector over time
fig8 = plt.figure()
plt.plot(times[0:n-1]/period,np.linalg.norm(w_vec,axis=1))
plt.title('B_dot convergence')
plt.xlabel('Period')
plt.ylabel('Norm of angular rate, [rad/s]')
plt.show()
# Plot North, East, Down (directly from IGRF) to see if singularities coming from pyIGRF or a coordinate transformation
fig9 = plt.figure()
plt.plot(B_field_NED_vec[:,0])
plt.plot(B_field_NED_vec[:,1])
plt.plot(B_field_NED_vec[:,2])
plt.title('Components of B field in NED (from pyIGRF)')
plt.show()