import sys, math, numpy as np

pi = math.pi

# inertia properties (add real later)
Ixx = 75
Iyy = 100
Izz = 125
I = np.array([[Ixx, 0, 0],[0, Iyy, 0], [0, 0, Izz]])

# deg --> rad, initial conditions
wx0 = 0.1*pi/180
wy0 = 0.1*pi/180
wz0 = pi/180
w = []
w.append([wx0, wy0, wz0])
print(w[0])
# quaternion, scalar last
q = []
q.append([0, 0, 0, 1])
print(q[0])

# moments, Nm
Mx0 = 0
My0 = 0
Mz0 = 0
M = []
M.append([Mx0, My0, Mz0])



print(I)
