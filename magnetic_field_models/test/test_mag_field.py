import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)

import numpy as np
import pytest
import math
import magnetic_field_cpp as mfcpp

# COMPARISONS DONE USING A MATLAB FUNCTION THAT MATCHES ONLINE CALCULATORS IT CAN BE FOUND HERE:
# https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model


def test_mag_field_1():
	alt = 400
	lat = 45
	lon = 30
	year = 2019
	order = 10
	# MATLAB CALL: igrf('1-Jan-2019', 45, 30, 6371.2 + 400, 'geocentric')' with nmax on line 343 as the original line
	mag_field_pred = np.array([18769.8277724085, 1819.07911050589, 36004.3487237974]) # from internet, FULL IGRF so expect error. 
	#This is found using http://www.geomag.bgs.ac.uk/cgi-bin/igrfsynth
	np.testing.assert_allclose(mfcpp.get_magnetic_field(lat, lon, alt, year, order), mag_field_pred, atol=25) # cpp test
	# I put an allowable error of 25 nT because of the order difference

def test_mag_field_2():
	# Testing against the truncated matlab function 
	# Note: Matlab function matches online calculators to double precision when using full model, so truncating will give exact answer to that order
	# MATLAB CALL: igrf('1-Jan-2019', 45, 30, 6371.2 + 400, 'geocentric')' with nmax on line 343 changed to 5
	alt = 400
	lat = 45
	lon = 30
	year = 2019
	order = 5
	mag_field_pred = np.array([18806.8694251561, 2023.78102588801, 36242.8613842709])
	np.testing.assert_allclose(mfcpp.get_magnetic_field(lat, lon, alt, year, order), mag_field_pred, atol=1e-15) # cpp test




