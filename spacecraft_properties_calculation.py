'''
Script for calculating spacecraft properties as designs change. All units should be SI or SI derived (kg, m, Tesla, etc.)
'''

#-------------------A priori properties----------------
rho = 1100000       # Density of Windform, [kg/m^3]
l = .05             # Side length of a Pocket-Qube, [m]
m_x = 8.8e-3        # Moments of magnetorquer triad, [A-m^2]
m_y = 1.373-2
m_z = 8.2e-3

#-------------------Derived properties-----------------
mass = rho*l**3
Ix = Iy = Iz = mass*l**2
print(Ix)
m_min = min([m_x,m_y,m_z])





