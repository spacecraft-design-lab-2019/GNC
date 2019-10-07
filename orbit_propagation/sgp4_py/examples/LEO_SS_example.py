from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import pyIGRF
import pyproj
'''
Function that loads in satellite object and returns propagated orbit
'''

# Orsted satellite TLE (https://www.calsky.com/observer/tle.cgi?satid=99008B&tdt=2456641.33063657)
line1 = ('1 25635U 99008B   13348.59627062  .00000650  00000-0  16622-3 0  9860')
line2 = ('2 25635 096.4421 173.2395 0141189 010.0389 029.8678 14.46831495780970')

# load in spacecraft object
satellite = twoline2rv(line1, line2, wgs72)
# propagte orbit based on epoch
year = 2013
day = 14
month = 12
hour = 14
minute = 18
second = 37
r, v = satellite.propagate(year, month, day, hour, minute, second)
x = r[0]
y = r[1]
z = r[2]
# Convert ECEF position to geodetic
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
# Calculate magnetic field properties
(D,I,H,Bx,By,Bz,F) = pyIGRF.igrf_variation(lat, lon, alt, year)

    
