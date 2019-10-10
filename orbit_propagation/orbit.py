'''

Author: Tunde Demuren
Date: 10-10-2019
Description: Orbit propagation functions to derive orbit state (position and magnetic field)

'''

import numpy as np
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.io import twoline2rv
import pyIGRF
import pyproj

def get_orbit_pos(TLE, epoch, wgs=wgs84):
    '''
    determine position in ecef from orbital elements (TLE) and time (epoch)
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
    # propagate epoch
    r, v = satellite.propagate(year, month, day, hour, minute, second)

    return r


def get_orbit_magnetic(TLE, epoch, wgs=wgs84):
    '''
    determine magnetic field properties from orbit state
    '''
    r = get_orbit_pos(TLE, epoch, wgs)
    # parse orbit state
    x = r[0]
    y = r[1]
    z = r[2]
    epoch = epoch.split('-')
    year = int(epoch[0])
    # Convert ECEF position to geodetic
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
    # Calculate magnetic field properties
    D,I,H,Bx,By,Bz,F = pyIGRF.igrf_variation(lat, lon, alt, year)

    return Bx, By, Bz



