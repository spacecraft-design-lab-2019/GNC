import numpy as np
def innovation(R,rN, rB, P,V):
    Q = np.zeros([6,6])
    Q[0:3,0:3] = R
    Q[3:6,3:6] = R

    
    z = rB - Q@rN
    

    C = np.array([rB,np.zeros(len(rB)), np.zeros(len(rB)),np.zeros(len(rB)),\
                  np.zeros(len(rB)),np.zeros(len(rB))])

    S = C@P@C.T+V
    
    return z, C, S
