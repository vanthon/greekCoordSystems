# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 18:14:28 2019

@author: anthon
"""

from math import sqrt, cos, sin, tan
from meridians import meridianLength, latFromMeridianLength


#def directTM(lat, lon, projName):
def directTM(lat, lon, ellipsoid, m0, lat0, lon0, E0, N0):
    """
    geodetic values in radians
    ellipsoid of class Ellipsoid
    m0 scale factor
    E0, N0 false Easting, Northing in meters
    returns projected coordinates of a point in meters
    """
    t = tan(lat)
    hta2 = ellipsoid.ePrimeSquared * (cos(lat)**2)
    dlon = lon-lon0

    T1 = m0 * meridianLength(lat, lat0, ellipsoid)
    T2 = (sin(lat) * cos(lat)) / 2
    T3 = ((sin(lat) * (cos(lat)**3)) / 24) * (5 - t**2 + (9 * hta2) + (4 *
                                              hta2**2))
    T4 = ((sin(lat) * (cos(lat)**5)) / 720) *\
         (61 - (58 * (t**2)) + t**4 + (270 * hta2) - (330*(t**2)*hta2) +
          (445 * hta2**2) + (324 * hta2**3) - (680*(t**2)*hta2**2) +
          (88 * hta2**4) - (660*(t**2)*hta2**3) - (192 * (t**2) * hta2**4))
    T5 = ((sin(lat) * (cos(lat)**7)) / 40320) *\
         (1385 - (3111 * (t**2)) + (543 * (t**4)) - t**6)
    T6 = cos(lat)
    T7 = ((cos(lat)**3) / 6) * (1 - t**2 + hta2)
    T8 = ((cos(lat)**5) / 120) *\
        (5 - (18 * (t**2)) + t**4 + (14 * hta2) - (58 * (t**2) * hta2) +
         (13 * (hta2**2)) + (4 * (hta2**3)) - (64 * (t**2) * (hta2**2)) -
         (24 * (t**2) * (hta2**3)))
    T9 = ((cos(lat)**7) / 5040) * (61 - (479 * (t**2)) + (179 * (t**4)) - t**6)

    W = sqrt(1 - ellipsoid.eSquared * (sin(lat))**2)
    Nc = ellipsoid.a / W

    E = E0 + (m0 * Nc) * ((T6 * dlon) + (T7 * (dlon**3)) + (T8 * (dlon**5)) +
                          (T9 * (dlon**7)))
    N = N0 + T1 + (m0 * Nc) * ((T2*dlon**2) + (T3*dlon**4) + (T4*dlon**6) +
                               (T5*dlon**8))
    return E, N


#def invTM(E, N, projName):
def invTM(E, N, ellipsoid, m0, lat0, lon0, E0, N0):
    """
    geodetic values in radians
    ellipsoid of class Ellipsoid
    m0 scale factor
    E0, N0 false Easting, Northing in meters
    returns projected coordinates of a point in meters
    """
    S = N / m0
    lat = latFromMeridianLength(lat0, S, ellipsoid)

    t = tan(lat)
    hta2 = ellipsoid.ePrimeSquared * (cos(lat)**2)
    E = E - E0
    W = sqrt(1 - ellipsoid.eSquared * (sin(lat))**2)
    N = ellipsoid.a / W
    Q = E / (m0 * N)
    M = ellipsoid.a * (1 - ellipsoid.eSquared) / (W**3)

    T10 = t / (2 * m0**2 * M * N)
    T11 = (t / (24 * m0**4 * M * N**3)) *\
          (5 + (3 * t**2) + hta2 - (4 * hta2**2) - (9 * t**2 * hta2))
    T12 = (t / (720 * m0**6 * M * N**5)) *\
          (61 + (90 * t**2) + (45 * t**4) + (46 * hta2) - (252 * t**2 * hta2) -
           (3 * hta2**2) + (100 * hta2**3) - (66 * t**2 * hta2**2) -
           (90 * t**4 * hta2) + (88 * hta2**4) + (225 * t**4 * hta2**2) +
           (84 * t**2 * hta2**3) - (192 * t**2 * hta2**4))
    T13 = (t / (40320 * m0**8 * M * N**7)) *\
          (1385 + (3633 * t**2) + (4095 * t**4) + (1575 * t**6))
    T14 = 1 / cos(lat)
    T15 = (1 + (2 * t**2) + hta2) / (6 * cos(lat))
    T16 = (5 + (6 * hta2) + (28 * t**2) - (3 * hta2**2) + (8 * t**2 * hta2) +
           (24 * t**4) - (4 * hta2**3) + (4 * t**2 * hta2**2) +
           (24 * t**2 * hta2**3)) / (120 * cos(lat))
    T17 = (61 + (66 * t**2) + (1320 * t**4) + (720 * t**6)) / 5040

    lat = lat - (T10 * E**2) + (T11 * E**4) - (T12 * E**6) + (T13 * E**8)
    lon = lon0 + (T14*Q) - (T15 * Q**3) + (T16 * Q**5) - (T17 * Q**7)
    return lat, lon
