import sys
import numpy as np

sys.path.append('../../euler/') 
import Euler

def measurement(q,rN):
    ''' 
    R = Rotation from N frame to B frame
    '''
    R = Euler.quat2DCM(q)
    rB1 = (R@rN[0:3])
    rB2 = (R@rN[3:6])
    y = np.append(rB1,rB2)
    return y,R

