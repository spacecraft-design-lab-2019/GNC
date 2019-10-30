# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 20:08:34 2019

@author: pdetr
"""
import numpy as np
from detumble.py_funcs import detumble_B_cross,detumble_B_dot,get_B_dot
from math import sqrt

# Example call of -omega detumble
omega = np.array([0.0, 0.0, 1.0])
B_1 = np.array([0.0, 1.0, 0.0])
B_2 = np.array([-sqrt(2), sqrt(2), 0.0])
dt = .1 # seconds

L_1 = detumble_B_cross(omega,B_2)
B_dot = get_B_dot(B_1,B_2,dt)
L_2 = detumble_B_dot(B_2,B_dot)

print(L_1)
print(L_2)
