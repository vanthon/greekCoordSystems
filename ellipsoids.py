# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 18:21:50 2019

@author: anthon
"""
from math import sqrt, cos, sin, atan2
from loadFiles import loadEllipsoidParams 
        
class Ellipsoid():
    """
    Ellipsoid class
    """
    def __init__(self, ellName):
        """
        ellName: string, 'WGS84', 'GRS80', 'BESSEL', 'HAYFORD'
        """
        self.name, self.a, self.b, self.f, self.eSquared, self.ePrimeSquared =\
        loadEllipsoidParams(ellName)
    
    def __repr__(self):
        return "EllipsoidID: " + self.name


def geodeticToECEF(lat, lon, h, ellipsoid):
    """
    lat, lon and h: in radians and meters/None
    ellipsoid: class ellipsoid
    returns ECEF Coords X, Y, Z in meters
    """
    if h is None: h = 0

    W = sqrt(1 - ellipsoid.eSquared * (sin(lat))**2)
    N = ellipsoid.a / W
    X = (N + h) * cos(lat) * cos(lon)
    Y = (N + h) * cos(lat) * sin(lon)
    Z = ((1 - ellipsoid.eSquared) * N + h) * sin(lat)
    return X, Y, Z


def ECEFtoGeodetic(x, y, z, ellipsoid):
    """
    x, y, z in meters
    ellipsoid: class ellipsoid
    returns Geodetic lat, lon in radians (between -pi & pi) and h in meters
    """    
    lon = atan2(y, x)
    P = sqrt(x**2 + y**2)
    lat0 = atan2((z * (1 + ellipsoid.ePrimeSquared)), P)
    dlat = 1
    while dlat > 10**-12:
        W = sqrt(1 - ellipsoid.eSquared * (sin(lat0))**2)
        N = ellipsoid.a / W
        lat = atan2((z + ellipsoid.eSquared * N * sin(lat0)), P)
        dlat = abs(lat - lat0)
        lat0 = lat
    h = (z / sin(lat)) - ((1 - ellipsoid.eSquared) * N)
    return lat, lon, h


def meridianLength(lat0, lat, ellipsoid):
    """
    geodetic values in radians
    ellipsoid: class ellipsoid
    returns meridian length in meters
    """
    e2 = ellipsoid.eSquared
    A0 = 1 - ((1/4)*e2) - ((3/64)*(e2**2)) - ((5/256)*(e2**3)) -\
             ((175/16384)*(e2**4))
    A2 = (3/8) * e2 * (1 + ((1/4)*e2) + ((15/128)*(e2**2)) +
                       ((35/512)*(e2**3)))
    A4 = (15/256) * (e2**2) * (1 + ((3/4)*e2) + ((35/64)*(e2**2)))
    A6 = (35/3072) * (e2**3) * (1 + ((5/4)*e2))
    A8 = (315/131072) * (e2**4)

    K0 = A0 * (lat - lat0)
    K2 = A2 * (sin(2*lat) - sin(2*lat0))
    K4 = A4 * (sin(4*lat) - sin(4*lat0))
    K6 = A6 * (sin(6*lat) - sin(6*lat0))
    K8 = A8 * (sin(8*lat) - sin(8*lat0))

    return abs(ellipsoid.a * (K0 - K2 + K4 - K6 + K8))


def latFromMeridianLength(lat0, S, ellipsoid):
    """
    lat0 in radians
    ellipsoid: class ellipsoid
    returns lat in radians
    """
    e2 = ellipsoid.eSquared
    A0 = 1 - ((1/4)*e2) - ((3/64)*(e2**2)) - ((5/256)*(e2**3)) -\
             ((175/16384)*(e2**4))
    A2 = (3/8) * e2 * (1 + ((1/4)*e2) + ((15/128)*(e2**2)) +
                       ((35/512)*(e2**3)))
    A4 = (15/256) * (e2**2) * (1 + ((3/4)*e2) + ((35/64)*(e2**2)))
    A6 = (35/3072) * (e2**3) * (1 + ((5/4)*e2))
    A8 = (315/131072) * (e2**4)

    dlat0 = S / (ellipsoid.a * A0)

    diff = 1
    while diff > 10**-12:
        latAvg = lat0 + (dlat0 / 2)
        K0 = S / (ellipsoid.a * A0)
        K2 = (2 * A2 / A0) * sin(dlat0) * cos(2 * latAvg)
        K4 = (2 * A4 / A0) * sin(2 * dlat0) * cos(4 * latAvg)
        K6 = (2 * A6 / A0) * sin(3 * dlat0) * cos(6 * latAvg)
        K8 = (2 * A8 / A0) * sin(4 * dlat0) * cos(8 * latAvg)
        dlat = K0 + K2 - K4 + K6 - K8
        diff = abs(dlat - dlat0)
        dlat0 = dlat

    return lat0 + dlat
