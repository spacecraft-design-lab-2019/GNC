from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

'''
Function that loads in satellite object and returns orbit
'''

# Orsted satellite TLE (https://www.calsky.com/observer/tle.cgi?satid=99008B&tdt=2456641.33063657)
line1 = ('1 25635U 99008B   13348.59627062  .00000650  00000-0  16622-3 0  9860')
line2 = ('2 25635 096.4421 173.2395 0141189 010.0389 029.8678 14.46831495780970')

def orbit_prop(TLE, epoch):
    # parse TLE
    line1 = TLE[0]
    line2 = TLE[1]
    # load in spacecraft object
    satellite = twoline2rv(line1, line2, wgs72)
    # propagte orbit based on epoch
    year = epoch[0]
    day = epoch[1]
    month = epoch[2]
    hour = epoch[3]
    minute = epoch[4]
    second = epoch[5]
    x, v = satellite.propagate(year, month, day, hour, minute, second)

    return x,v
