import numpy as np
import math
from predict import predict

xk = np.array([0,0,0,1,0.1,0.1,0.1])
w = np.array([0.15,0.5,0.1])
dt = 0.01

xn,A = predict(xk,w,dt)
