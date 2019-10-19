import sys
import numpy as np

sys.path.append('../../euler/') 
import Euler

def measurement(q,rN):
    ''' 
    Q = Rotation from N frame to B frame
    '''
    Q = Euler.quat2DCM(q)
    rB = (Q@rN)
    y = rB
    C = np.array([rB,np.zeros(len(rB))]).T
    return C,y

