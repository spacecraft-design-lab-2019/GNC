import math

def sun_position(MJD):
    """
    This is using the equations given in Motenbruck and Gill's Satellite Orbits book
    Inputs:
    MJD - Modified Julian Day (J2000) as a Real Number
    Outputs:
    r_vec - list with x, y, z of Sun position in ECI at input time
    """
    JD = MJD + 2400000.5
    OplusW = 282.94
    T = (JD - 2451545.0) / 36525

    M = math.radians(357.5256 + 35999.049 * T)

    long = math.radians(OplusW + math.degrees(M) + 6892 / 3600 * math.sin(M) + 72 / 3600 * math.sin(2*M))
    r_mag = (149.619 - 2.499 * math.cos(M) - 0.021 * math.cos(2*M)) * 10**6

    epsilon = math.radians(23.43929111)
    r_vec = [r_mag * math.cos(long), r_mag * math.sin(long) * math.cos(epsilon), r_mag * math.sin(long) * math.sin(epsilon)]

    return r_vec

def sat_sun_vect(r, MJD):
    """
    Returns the unit vector from the satellite to the Sun in ECI coordinates
    Inputs:
    r - ECI position of the satellite
    MJD - Julian Day (J2000) as a Real Number
    Outputs:
    r_sat_sun - numpy array giving the unit direction to the Sun from the satellite
    """

    r_sun = sun_position(MJD)
    r_sat_sun = [x - y for x, y in zip(r_sun,r)]

    r_sat_sun = [x/norm2(r_sat_sun) for x in r_sat_sun]

    return r_sat_sun

def norm2(v):
    """
    Caluclates the 2-norm of a vector
    """
    # return (v.T @ v) ** (0.5)
    return math.sqrt(sum(x*x for x in v))

def ECEF_to_LLA(r_ECEF, rad_Earth):
    """
    Function: ECEF_to_LLA
        Converts position vector in ECEF to geocentric coordinates.
    Inputs:
        r_ECEF: position vector in Earth Centered Earth Fixed (ECEF)
        rad_Earth: radius of the Earth [km]
    Outputs:
        lat:    latitude    [rad]
        long:   longitude   [rad]
        alt:    altitude    [rad]
    """

    lat = math.asin(r_ECEF[2] / norm2(r_ECEF))
    lon = math.atan2(r_ECEF[1], r_ECEF[0])

    alt = norm2(r_ECEF) - rad_Earth

    return lat, lon, alt

def dot(v1, v2):
    return sum(x*y for x,y in zip(v1,v2))

def transpose(M):
    I = range(len(M))
    J = range(len(M[0]))
    return [[M[i][j] for i in I] for j in J]

def ECI_to_ECEF(GMST):
    """
    Rotation matrix from ECI to ECEF coordinates
    Inputs:
    GMST - Greenwich Mean Sidereal Time
    Outputs:
    R - Rotation matrix from ECI to ECEF
    """

    rotation = [[math.cos(GMST), math.sin(GMST), 0],
                [-math.sin(GMST), math.cos(GMST), 0],
                [0, 0, 1]]

    return rotation

def date2MJD(M, D, Y, HH, MM, SS):
    """
    Gives the Modified Julian Date from the date and time using Vallado algorithm
    Inputs:
    M - Month number (January = 1, December = 12)
    D - Day number
    Y - Year
    HH - Hour in 24 hour clock
    MM - Minutes
    SS - Seconds
    Outputs:
    MJD - Modified Julian Date
    """
    assert valid_date(M, D, Y, HH, MM, SS)
    if M <= 2:
        y = Y - 1
        m = M + 12
    else:
        y = Y
        m = M

    if Y <= 1582 and M <= 10 and D <= 4:
        B = -2 + ((y + 4716) / 4) - 1179
    else:
        # B = y / 400 - y / 100 + y / 4
        B = 2 - math.floor(y/100) + math.floor(math.floor(y/100) / 4)

    day_frac = HH / 24 + MM / 60 / 24 + SS / 60 / 60 / 24

    # MJD = 365 * y - 679004 + np.floor(B) + np.floor(30.6001 * (m + 1)) + D + day_frac;
    MJD = math.floor(365.25*(y + 4716)) + math.floor(30.6001*(m+1))+ D + B + day_frac - 2401525.0
    return MJD

def valid_date(M, D, Y, HH, MM, SS):
    """
    Ensures the calendar date is valid
    Inputs:
    M - Month
    Y - Year
    D - Day
    HH - Hours
    MM - Minutes
    SS - Seconds
    Outputs:
    check - true/false
    """
    check1 = False
    check2 = False
    check3 = False
    if isinstance(M, int) and isinstance(Y, int) and isinstance(D, int):
        check1 = True
    if M <= 12 and D <= 31 and HH <= 24 and MM <= 60 and SS <= 60:
        check2 = True
    if M >= 0 and D >= 0:
        check3 = True
    if check1 and check2 and check3:
        check = True
    else:
        check = False

    return check

def MJD2GMST(mjd):
    """
    Function: mjd_2_GMST
        Calculates Greenwich Mean Sidereal Time
    Inputs:
        mjd - modified julian day
    Outputs:
        GMST - Greenwich Mean Sidereal Time [rad]
    Reference: Vallado
    """
    d = (mjd - 51544.5) / 36525.0
    GMST = 67310.54841 + (876600*3600 + 8640184.812866) * d + 0.093104 * d**2 - 6.2 * 10**(-6)*d**3
    return (GMST % 86400) / 240 * (math.pi / 180)

def ecef2enu(lat, lon):
    """
    Rotation matrix from ECEF to ENU coordinates
    Inputs:
    lat - Latitude in radians
    lon - Longitude in radians
    Outputs:
    R - Rotation matrix from ECEF to ENU
    """
    Ehat = [-math.sin(lon), math.cos(lon), 0]
    Nhat = [-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)]
    Uhat = [math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)]
#
#    R = np.column_stack((Ehat, Nhat, Uhat))
    R = transpose(R)
    return R

