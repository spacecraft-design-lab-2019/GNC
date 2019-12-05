import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

import math
import pyIGRF

import magnetic_field_cpp


### Setting up parameters to run sweep
alt = 400
year = 2015
min_lat = -90
max_lat = 90
min_lon = -180
max_lon = 180
delta = 2

lat_range = np.arange(min_lat, max_lat, delta)
lon_range = np.arange(min_lon, max_lon, delta)



mag_field_strength_cpp = np.zeros((np.size(lat_range), np.size(lon_range)))
mag_field_strength_igrf = np.zeros((np.size(lat_range), np.size(lon_range)))
err = np.zeros((np.size(lat_range), np.size(lon_range)))
for i, lat in enumerate(lat_range):
	for j, lon in enumerate(lon_range):
		
		B_vec = magnetic_field_cpp.get_magnetic_field(lat, lon, alt, year, 10)
		mag_field_strength_cpp[i, j] = np.linalg.norm(B_vec)

		D, I, H, Bx, By, Bz, F = pyIGRF.igrf_value(lat, lon, alt, year)
		mag_field_strength_igrf[i, j] = F
		err[i, j] = abs(mag_field_strength_cpp[i, j] - mag_field_strength_igrf[i, j])

X, Y = np.meshgrid(lon_range, lat_range)


fig1, ax = plt.subplots()
CS = ax.contour(X, Y, mag_field_strength_cpp, 30)
CB = fig1.colorbar(CS, shrink=0.8, extend='both')
ax.set_title('5th Order IGRF (Our Model)')
# plt.show()

fig2, ax2 = plt.subplots()
CS = ax2.contour(X, Y, mag_field_strength_igrf, 30)
CB = fig2.colorbar(CS, shrink=0.8, extend='both')
ax2.set_title('Full IGRF (True Model)')
plt.show()

fig2, ax2 = plt.subplots()
CS = ax2.contour(X, Y, err, 30)
CB = fig2.colorbar(CS, shrink=0.8, extend='both')
ax2.set_title('Error')
plt.show()


