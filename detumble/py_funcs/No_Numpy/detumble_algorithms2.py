import flight_utils as ut
import math

def get_B_dot(B1,B2,dt):
    '''
    Takes in two magnetic field measurements (3x1) and the timestep between them and returns
    a first order approximation to the rate of change of the magnetic field (3x1)
    '''
    B_dot = [x*(1/dt) for x in sub(B2,B1)]
    return B_dot

def sub(vec1, vec2):
    """
    Inputs: 2 vectors
    Outputs: vector1 - vector2 in vector form
    """
    return [x - y for x, y in zip(vec1, vec2)]

def detumble_B_dot_bang_bang(B_dot, max_dipoles = [[8.8e-3],[1.373e-2],[8.2e-3]]):
    '''
    :param B_dot: 3x1 vector of the rate of change of Earth's magnetic field, [nanoTesla/s], as measured in the spacecraft body
                    frame.
    :param max_dipoles: 3x1 vector of max dipole values [Amp-m^2] in the x, y, and z body (principal) axis directions.
    :return: m, a 3x1 vector of the dipole commanded by the bang-bang b_dot control law.
    references: Wertz, equation (7.54)
    '''
    m = [0,0,0]
    m[0] = max_dipoles[0] * -math.copysign(1, ut.dot([1.0, 0.0, 0.0] , B_dot))
    m[1] = max_dipoles[1] * -math.copysign(1, ut.dot([0.0, 1.0, 0.0] , B_dot))
    m[2] = max_dipoles[2] * -math.copysign(1, ut.dot([0.0, 0.0, 1.0], B_dot))
    return m