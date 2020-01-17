'''
Script integrating detumble with orbit/magnetic field knowledge
'''

# from detumble.py_funcs import detumble_B_cross,detumble_B_dot,get_B_dot, detumble_B_dot_bang_bang
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
import euler_cpp as ecpp

# clear figures
plt.close('all')
        
pi = math.pi


# inertia properties (add real later)
# Ixx = 0.34375
# Iyy = 0.34375
# Izz = 0.5
I = np.array([[17,0,0],[0,18,0],[0,0,22]])
max_dipoles = np.array([[8.8e-3], [1.373e-2], [8.2e-3]])

# initial attitude conditions, radians & rad/s
q_0 = np.array([[math.sqrt(4.0)/4.0],[math.sqrt(4.0)/4.0],[math.sqrt(4.0)/4.0],[math.sqrt(4.0)/4.0]])                     # initial quaternion, scalar last
w_0 = np.array([[.03],[.03],[.03]])  # initial rotation rate, rad/s
# initial state: quaternion, rotation rate
x_0 = np.squeeze(np.concatenate((q_0,w_0)))

# initial orbit state conditions, TLE+epoch
epoch = '2019-12-30T00:00:00.00'
line1 = ('1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991')
line2 = ('2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102')
TLE = {'line1': line1, 'line2': line2}

# initial orbit time (fair warning, this is the time for PyCubed, not Orsted)
MJD = 58847.0
GMST_0 = tfcpp.MJD2GMST(MJD)

mean_motion = 14.46/(24*3600)*2*math.pi # mean motion, radians/second
period = 2*pi/mean_motion                      # Period, seconds

# feed in a vector of times and plot orbit
t0 = 0.0
tf = 60000
tstep = .1
times = np.arange(t0,tf,tstep)
n = len(times)

# preallocate position storage matrix
positions_ECI = np.zeros((n-1,3))
positions_ECEF = np.zeros((n-1,3))  
B_field_body = np.zeros((n-1,3))
B_field_ECI_vec = np.zeros((n-1,3))
B_field_NED_vec = np.zeros((n-1,3))
B_dot_body = np.zeros((n-1,3))
w_vec_b_dot = np.zeros((n-1,3))
w_vec_b_cross_bang_bang = np.zeros((n-1,3))
w_vec_b_cross_directional = np.zeros((n-1,3))
q_vec = np.zeros((n-1,4))
M_vec = np.zeros((n-1,3))
dipole_vec = np.zeros((n-1,3))

# Define function for calculating full state derivative
t = time.time()


#--------------------------------------------------------B_dot--------------------------------------------------------------
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
    # R_body2ECI = quat2DCM(np.transpose(x[0:4]))
    B_field_ECEF = np.transpose(R_ECEF2ENU) @ B_field_ENU
    B_field_ECI = np.transpose(R_ECI2ECEF) @ B_field_ECEF

    B_field_ECI_vec[i,:] = np.transpose(B_field_ECI)
    q_ECI2body = ecpp.get_inverse_quaternion(x[0:4])
    B_field_body[i,:] = ecpp.rotate_vec(B_field_ECI, q_ECI2body)

    # Get B_dot based on previous measurement
    if i>0:
        B_1 = np.transpose(B_field_body[i-1,:])
        B_2 = np.transpose(B_field_body[i,:])
        B_dot = dcpp.get_B_dot(B_1,B_2,tstep)
        B_dot_body[i,:] = np.transpose(B_dot)

        # Validate B_dot bang bang control law:
        # include 1e-9 factor to get nanoTesla back into SI units
        dipole = dcpp.detumble_B_dot_bang_bang(np.transpose(B_dot_body[i,:]),max_dipoles)
        dipole_vec[i,:] = np.transpose(dipole)
        bang_bang_gain = 1e-9 
        M = np.cross(np.squeeze(dipole), bang_bang_gain*np.transpose(B_field_body[i, :]))
    else:
        M = np.zeros((3,1))

    # store for plotting
    M_vec[i,:] = np.transpose(M)

    # Propagate dynamics/kinematics forward using commanded moment
    y = integrate.odeint(ecpp.get_attitude_derivative, x, (times[i],times[i+1]), args = (M, I), tfirst=True)

    # Store angular velocity
    w_vec_b_dot[i,:] = y[-1,4:7]
    q_vec[i,:] = y[-1,0:4]

    # Update full attitude state
    x = np.transpose(y[-1,:])

elapsed = time.time() - t
print(elapsed)

# #---------------------------------------------B_cross_bang_bang----------------------------------------------------
# x = x_0;
# for i in range(len(times)-1):
#     # Get GMST at this time
#     GMST = tfcpp.MJD2GMST(MJD + times[i] / 60.0 / 60.0 / 24.0)

#     positions_ECI[i,:] = get_orbit_pos(TLE, epoch, times[i])

#     # convert to ECEF
#     R_ECI2ECEF = fccpp.eci2ecef(GMST)

#     positions_ECEF[i,:] = np.transpose(R_ECI2ECEF @ np.transpose(positions_ECI[i,:]))
#     lat, lon, alt = fccpp.ecef2lla(np.transpose(positions_ECEF[i,:]))
#     R_ECEF2ENU = fccpp.ecef2enu(lat, lon)

#     # get magnetic field at position
#     B_field_NED = get_B_field_at_point(positions_ECEF[i,:]) # North, East, Down
#     B_field_NED_vec[i,:] = np.transpose(B_field_NED)        # store for later analysis
#     B_field_ENU = np.array([[B_field_NED[1]],[B_field_NED[0]],[-B_field_NED[2]]])    # north, east, down to east, north, up

#     # get magnetic field in body frame (for detumble algorithm)
#     # R_body2ECI = quat2DCM(np.transpose(x[0:4]))
#     B_field_ECEF = np.transpose(R_ECEF2ENU) @ B_field_ENU
#     B_field_ECI = np.transpose(R_ECI2ECEF) @ B_field_ECEF

#     B_field_ECI_vec[i,:] = np.transpose(B_field_ECI)
#     q_ECI2body = ecpp.get_inverse_quaternion(x[0:4])
#     B_field_body[i,:] = ecpp.rotate_vec(B_field_ECI, q_ECI2body)

#     omega = x[4:7]
#     gain = 1
#     dipole = dcpp.detumble_B_cross_bang_bang(np.transpose(omega), np.transpose(B_field_body[i,:]), gain, max_dipoles)
#     bang_bang_gain = 1e-9 
#     M = np.cross(np.squeeze(dipole), bang_bang_gain*np.transpose(B_field_body[i, :]))


#     # store for plotting
#     M_vec[i,:] = np.transpose(M)

#     # Propagate dynamics/kinematics forward using commanded moment
#     y = integrate.odeint(ecpp.get_attitude_derivative, x, (times[i],times[i+1]), (M, I), tfirst=True)

#     # Store angular velocity
#     w_vec_b_cross_bang_bang[i,:] = y[-1,4:7]
#     q_vec[i,:] = y[-1,0:4]

#     # Update full attitude state
#     x = y[-1,:]

# elapsed = time.time() - t
# print(elapsed)

# #-----------------------------------------------B_cross_directional----------------------------------------------------------
# x = x_0;
# for i in range(len(times)-1):
#     # Get GMST at this time
#     GMST = tfcpp.MJD2GMST(MJD + times[i] / 60.0 / 60.0 / 24.0)

#     positions_ECI[i,:] = get_orbit_pos(TLE, epoch, times[i])

#     # convert to ECEF
#     R_ECI2ECEF = fccpp.eci2ecef(GMST)

#     positions_ECEF[i,:] = np.transpose(R_ECI2ECEF @ np.transpose(positions_ECI[i,:]))
#     lat, lon, alt = fccpp.ecef2lla(np.transpose(positions_ECEF[i,:]))
#     R_ECEF2ENU = fccpp.ecef2enu(lat, lon)

#     # get magnetic field at position
#     B_field_NED = get_B_field_at_point(positions_ECEF[i,:]) # North, East, Down
#     B_field_NED_vec[i,:] = np.transpose(B_field_NED)        # store for later analysis
#     B_field_ENU = np.array([[B_field_NED[1]],[B_field_NED[0]],[-B_field_NED[2]]])    # north, east, down to east, north, up

#     # get magnetic field in body frame (for detumble algorithm)
#     # R_body2ECI = quat2DCM(np.transpose(x[0:4]))
#     B_field_ECEF = np.transpose(R_ECEF2ENU) @ B_field_ENU
#     B_field_ECI = np.transpose(R_ECI2ECEF) @ B_field_ECEF

#     B_field_ECI_vec[i,:] = np.transpose(B_field_ECI)
#     q_ECI2body = ecpp.get_inverse_quaternion(x[0:4])
#     B_field_body[i,:] = ecpp.rotate_vec(B_field_ECI, q_ECI2body)

#     omega = x[4:7]
#     gain = 1000
#     dipole = dcpp.detumble_B_cross_directional(np.transpose(omega), np.transpose(B_field_body[i,:]), gain, max_dipoles)
#     bang_bang_gain = 1e-9 
#     M = np.cross(np.squeeze(dipole), bang_bang_gain*np.transpose(B_field_body[i, :]))


#     # store for plotting
#     M_vec[i,:] = np.transpose(M)

#     # Propagate dynamics/kinematics forward using commanded moment
#     y = integrate.odeint(ecpp.get_attitude_derivative, x, (times[i],times[i+1]), (M, I), tfirst=True)

#     # Store angular velocity
#     w_vec_b_cross_directional[i,:] = y[-1,4:7]
#     q_vec[i,:] = y[-1,0:4]

#     # Update full attitude state
#     x = y[-1,:]

# elapsed = time.time() - t
# print(elapsed)


#----------------------------------------------Plotting--------------------------------------------------------------------
# plot norm of velocity vector over time
fig1 = plt.figure()
plt.plot(times[0:n-1]/3600.0,np.linalg.norm(w_vec_b_dot,axis=1))
plt.legend(('B_dot','B_cross_bang_bang','B_cross_directional'))
plt.title('Detumble algorithm convergence')
plt.xlabel('Period')
plt.ylabel('Norm of angular rate, [rad/s]')

# plot angular velocity over time
fig2 = plt.figure()
plt.plot(times[0:n-1]/3600.0,w_vec_b_dot[:,0])
plt.plot(times[0:n-1]/3600.0,w_vec_b_dot[:,1])
plt.plot(times[0:n-1]/3600.0,w_vec_b_dot[:,2])
plt.title('Angular velocity components')

# plot B field components as a function of time
fig3 = plt.figure()
plt.plot(times[0:n-1]/3600.0,B_field_body[:,0])
plt.plot(times[0:n-1]/3600.0,B_field_body[:,1])
plt.plot(times[0:n-1]/3600.0,B_field_body[:,2])
plt.title('Body Magnetic field components, B')

# plot B dotcomponents as a function of time
fig4 = plt.figure()
plt.plot(times[0:n-1]/3600.0,B_dot_body[:,0])
plt.plot(times[0:n-1]/3600.0,B_dot_body[:,1])
plt.plot(times[0:n-1]/3600.0,B_dot_body[:,2])
plt.title('Magnetic field rate of change, B_dot')
plt.show()
# plot B dot components as a function of time
fig5 = plt.figure()
plt.plot(times[0:n-1]/3600.0,B_field_ECI_vec[:,0])
plt.plot(times[0:n-1]/3600.0,B_field_ECI_vec[:,1])
plt.plot(times[0:n-1]/3600.0,B_field_ECI_vec[:,2])
plt.title('Magnetic field ECI')


# plot quaternion over time
fig6 = plt.figure()
plt.plot(times[0:n-1]/3600.0,q_vec[:,0])
plt.plot(times[0:n-1]/3600.0,q_vec[:,1])
plt.plot(times[0:n-1]/3600.0,q_vec[:,2])
plt.plot(times[0:n-1]/3600.0,q_vec[:,3])
plt.title('quaternion components')

# plot Moment over time
fig7 = plt.figure()
plt.plot(times[0:n-1]/3600.0,M_vec[:,0])
plt.plot(times[0:n-1]/3600.0,M_vec[:,1])
plt.plot(times[0:n-1]/3600.0,M_vec[:,2])
plt.title('Moment components')

# plot dipole over time
fig8 = plt.figure()
plt.plot(times[0:n-1]/3600.0,dipole_vec[:,0])
plt.plot(times[0:n-1]/3600.0,dipole_vec[:,1])
plt.plot(times[0:n-1]/3600.0,dipole_vec[:,2])
plt.title('dipole components')
plt.show()



