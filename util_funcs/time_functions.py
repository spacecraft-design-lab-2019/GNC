import numpy as np

def date2MJD(M, D, Y, HH, MM, SS):
    if M <= 2:
        y = Y - 1
        m = M + 12
    else:
        y = Y
        m = M

    if Y <= 1582 and M <= 10 and D <= 4:
        B = -2 + ((y + 4716) / 4) - 1179
    else:
        B = y / 400 - y / 100 + y / 4

    day_frac = HH / 24 + MM / 60 / 24 + SS / 60 / 60 / 24

    MJD = 365 * y - 679004 + np.floor(B) + np.floor(30.6001 * (m + 1)) + D + day_frac;

    return MJD

