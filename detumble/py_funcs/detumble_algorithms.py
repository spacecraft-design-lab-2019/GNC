# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 19:36:55 2019

@author: Paul DeTrempe

@description: Implementation of magnetorquer detumbling laws (See Wertz 7.5)
"""
import numpy as np

def detumble_B_cross(omega,B,k):
    '''
    Takes in the angular rate and magnetic field vectors (in principal frame)
    (as 3x1 np arrays) and returns a 3x1 vector of control torque to detumble.
    
    - based on the paper by AVANZINI AND GIULIETTI
    TODO: find optimal k for our system
    '''
    # try adjusting k to the velocity
    # k = -2*omega
    b = B/np.linalg.norm(B) # normalize magnetic field vector
    L = -k* (np.identity(3) - b @ np.transpose(b)) @ omega
    return L
    
def detumble_B_dot(B,B_dot,k=1.0):
    '''
    Takes in magnetic field, magnetic field rate (as 3x1 vectors, in principal frame), and control gain (scalar)
    and returns a 3x1 control moment
    
    TODO: find optimal k for our system
    '''
    m = -k*B_dot/np.linalg.norm(B)
    return m
    
def get_B_dot(B1,B2,dt):
    '''
    Takes in two magnetic field measurements (3x1) and the timestep between them and returns
    a first order approximation to the rate of change of the magnetic field (3x1)
    '''
    B_dot = (B2-B1)/dt
    return B_dot
    