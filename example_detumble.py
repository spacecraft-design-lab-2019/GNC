# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 20:08:34 2019

@author: pdetr
"""
import numpy as np
from detumble import detumble_avanzini


# Example call of -omega detumble
omega = np.array([0.0,0.0,1.0])
B = np.array([0.0,1.0,0.0])

L = detumble_avanzini(omega,B)

print(L)