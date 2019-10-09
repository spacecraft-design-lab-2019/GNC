# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 19:36:55 2019

@author: Paul DeTrempe

@description: Implementation of magnetorquer detumbling laws (See Wertz 7.5)
"""
import numpy as np

def detumble_avanzini(omega,B,k=1):
    '''
    Takes in the angular rate and magnetic field vectors (in principal frame)
    (as 3x1 np arrays) and returns a 3x1 vector of control torque to detumble.

    TODO: find optimal k for our system
    '''
    b = B/np.linalg.norm(B) # normalize magnetic field vector
    L = -k*np.matmul((np.identity(3) - np.matmul(b,np.transpose(b))),omega)
    return L
    
def detumble_B_dot(B,B_dot,k=1):
    '''
    Takes in magnetic field, magnetic field rate (as 3x1 vectors, in principal frame), and control gain (scalar)
    and returns a 3x1 control moment
    
    TODO: find optimal k for our system
    '''
    m = -k*B_dot/np.linalg.norm(B)
    return m
    
    