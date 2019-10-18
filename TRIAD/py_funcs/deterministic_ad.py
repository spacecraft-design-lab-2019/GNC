import numpy as np
def triad_ad(M, V):
    """
    Computes rotation matrix from inertial to body frame using measurement vectors and modelled measurement vectors
    Inputs:
    M - Matrix where each column is a measurement vector in the body frame - Note, most accurate measurement should be in the first column
    V - Matrix where each column is a modeled vector of the corresponding measurement vector from M, in the inertial frame
    Outputs:
    R - rotation matrix from inertial to body frame
    """
    assert M.shape[0] == V.shape[0] == 3, "Either M or V vectors aren't length 3"

    if M.shape[1] == 2 and V.shape[1] == 2:
        m2 = np.cross(M[:, 0], M[:, 1])
        m3 = np.cross(M[:, 0], m2)
        v2 = np.cross(V[:, 0], V[:, 1])
        v3 = np.cross(V[:, 0], v2)

        R = np.column_stack((M[:, 0], m2, m3)) @ np.linalg.inv(np.column_stack((V[:, 0], v2, v3)))
    else:
        R = M @ np.linalg.inv(V)

    return R