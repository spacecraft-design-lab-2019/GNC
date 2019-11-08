'''

Author: Tunde Demuren
Date: 10-10-2019
Description: Orbit propagation functions to derive orbit state (position and magnetic field)

'''
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
gncdir = os.path.dirname(parentdir)
docdir = os.path.dirname(gncdir)
sys.path.insert(0,parentdir)
sys.path.insert(0, gncdir)
sys.path.insert(0, docdir)

from util_funcs.py_funcs import frame_conversions as fc
import numpy as np
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.io import twoline2rv
import pyIGRF
import math
# import pyproj

def get_orbit_pos(TLE, epoch, sec_since_epoch, wgs=wgs84):
    '''
    determine position in ECI from TLE, epoch (as a string), and seconds past epoch 
    (as a float)
    '''
    # parse TLE
    line1 = TLE['line1']
    line2 = TLE['line2']
    # parse epoch
    epoch = epoch.split('-')
    year, month = int(epoch[0]), int(epoch[1])
    day_time = epoch[2].split('T')
    day = int(day_time[0])
    time = day_time[1].split(':')
    hour, minute = int(time[0]), int(time[1])
    second = int(time[2].split('.')[0])
    # load in spacecraft object
    satellite = twoline2rv(line1, line2, wgs)
    # propagate epoch + time
    r, v = satellite.propagate(year, month, day, hour, minute, second + sec_since_epoch)

    return r

def get_orbit_state(TLE, epoch, sec_since_epoch, wgs=wgs84):
    '''
    determine state in ECI from TLE, epoch (as a string), and seconds past epoch 
    (as a float)
    '''
    # parse TLE
    line1 = TLE['line1']
    line2 = TLE['line2']
    # parse epoch
    epoch = epoch.split('-')
    year, month = int(epoch[0]), int(epoch[1])
    day_time = epoch[2].split('T')
    day = int(day_time[0])
    time = day_time[1].split(':')
    hour, minute = int(time[0]), int(time[1])
    second = int(time[2].split('.')[0])
    # load in spacecraft object
    satellite = twoline2rv(line1, line2, wgs)
    # propagate epoch + time
    r, v = satellite.propagate(year, month, day, hour, minute, second + sec_since_epoch)
    state = np.concatenate([np.array(r),np.array(v)])
    return state


def get_orbit_magnetic(TLE, epoch, sec_past_epoch, wgs=wgs84):
    '''
    determine magnetic field properties from orbit state, in ECEF!!!!
    '''
    r = get_orbit_pos(TLE, epoch, sec_past_epoch, wgs)
    # parse orbit state
    x = r[0]
    y = r[1]
    z = r[2]
    epoch = epoch.split('-')
    year = int(epoch[0])
    # Convert ECEF position to geodetic
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
    lat, lon, alt = fc.ecef2lla(r)
    # Calculate magnetic field properties
    lat = lat * 180 / math.pi
    lon = lon * 180 / math.pi
    D,I,H,Bx,By,Bz,F = pyIGRF.igrf_value(lat, lon, alt, year)

    return Bx, By, Bz

def get_B_field_at_point(r, year=2019):
    '''
    Find the North-East-Down B-field (in nT) in ECEF coordinates 
    '''
    # parse position
    x = r[0]
    y = r[1]
    z = r[2]
    # Convert ECEF position to geodetic
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
    lat, lon, alt = fc.ecef2lla(r)
    lat = lat * 180 / math.pi
    lon = lon * 180 / math.pi
    # Calculate magnetic field properties
    D, I, H, Bx, By, Bz, F = pyIGRF.igrf_value(lat, lon, alt, year)

    return Bx, By, Bz

