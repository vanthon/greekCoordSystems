# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:27:10 2019

@author: Antonios Vatalis
"""
from math import sqrt, cos, sin, atan2
from angles import DMStoRads, radToDMS


class GeoPoint():
    def __init__(self, x, y, z, coordType, ellipsoid, geoidHeight = None):
        """
        This method initializes an object of class GeoPoint

        x, y, z: geocentric or geodetic coordinates of point
                 coordType: 'ECEF" -> x, y, z: in meters
                 coordType: 'latlonh-(r/dd/dm/dms)"
                             x: lar in (radians/dd/dm/dms)
                             y: lon in (radians/dd/dm/dms)
                             z: h in meters
        coordType: one of the following strings:
                        'latlonh-r', 'latlonh-dd', 'latlonh-dm', 'latlonh-dms',
                        'ECEF'
        ellipsoid: object of class Ellipsoid
        geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
        """
        self.x = x
        self.y = y
        self.z = z
        self.coordType = coordType
        self.ellipsoid = ellipsoid
        self.geoidHeight = geoidHeight
        # default precs
        self.precs = {'ECEF': [3, 3],
                 'latlonh-r': [6, 3],
                 'latlonh-dd': [6, 3],
                 'latlonh-dm': [6, 3],
                 'latlonh-dms': [5, 3]}

    def radians(self):
        """
        This method transforms attributes x, y in radians
        (not if attribute coordType is 'ECEF')
        """
        if self.coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
            self.coordType = "latlonh-r"
            self.x = DMStoRads(self.x)
            self.y = DMStoRads(self.y)
            
    def latlonh(self, newCoordType):
        """
        This method transforms attributes x, y, z in geodetic coordinates: 
            lat, lon, h
        The type of lat, lon is defined by newCoordType parameter

        newCoordType: one of the following strings: 'r', 'dd', 'dm', 'dms'
        """
        # latlonh to latlonh
        if self.coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm", 
                              "latlonh-dms"]:
            self.radians()
            if newCoordType != 'r':
                self.x = radToDMS(self.x, newCoordType)
                self.y = radToDMS(self.y, newCoordType)
        # ECEF to latlonh
        elif self.coordType == "ECEF":
            geodeticPoint = ECEFtoGeodetic(self, newCoordType)
            self.x, self.y, self.z =\
                geodeticPoint.x, geodeticPoint.y, geodeticPoint.z
        
        self.coordType = "latlonh-" + newCoordType
        
    def ECEF(self):
        """
        This method transforms attributes x, y, z in ECEF coordinates
        """
        # latlonh to ECEF
        if self.coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm",
                              "latlonh-dms"]:
            self.radians()
            ecefPoint = geodeticToECEF(self)
            self.x, self.y, self.z = ecefPoint.x, ecefPoint.y, ecefPoint.z
        self.coordType = "ECEF"

    def setPrec(self, precXY, precZ):
        """
        This method specifies the precision of attributes strX, strY, strZ

        precXY: int, precision of attributes strX, strY
        precZ:  int, precision of attribute strZ

        e.g. a = GeoPoint(x, y, z, coordType, ellipsoid)
             a.strPrec(precXY, precZ)
             print(a) # prints a with its coords in desired precisions
        """
        # set precision of x and y
        if type(self.x) is str:
            self.strX = radToDMS(DMStoRads(self.x),
                                 self.coordType.split('-')[1], precXY)
            self.strY = radToDMS(DMStoRads(self.y),
                                 self.coordType.split('-')[1], precXY)
        else:
            self.strX = str(round(self.x, precXY))
            self.strY = str(round(self.y, precXY))
        # set precision of z
        if self.z is not None:
            self.strZ = str(round(self.z, precZ))
        else:
            self.strZ = str(self.z)
        if self.geoidHeight is not None:
            self.strGeoidHeight = str(round(self.geoidHeight, precZ))
        else:
            self.strGeoidHeight = str(self.geoidHeight)

        self.precs[self.coordType] = [precXY, precZ]

    def __repr__(self):
        # default precision of strX, strY, strZ
        self.setPrec(self.precs[self.coordType][0],
                     self.precs[self.coordType][1])

        return "(" + ", ".join([self.strX, self.strY, self.strZ,
                                self.coordType, self.ellipsoid.name,
                                self.strGeoidHeight]) + ")"


class MapPoint():
    def __init__(self, E, N, U, projection, geoidHeight = None):
        """
        This method initializes an object of class MapPoint

        E, N, U: map Eastings and Northings of point in meters
        U:       orthometric height of point in meters
        projection: object of class Projection where the MapPoint belongs to
        lat0, lon0: in radians. They matter only if projection belongs to 
                                    HATT family.
        geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
        """
        self.E = E
        self.N = N
        self.U = U
        self.projection = projection
        self.ellipsoid = self.projection.ellipsoid
        # geoid height
        self.geoidHeight = geoidHeight
        # default prec
        self.prec = 3

    def setPrec(self, prec):
        """
        This method specifies the precision of attributes strEm, strN, strU

        prec: int

        e.g. a = MapPoint(E, N, U, projection)
             a.strPrec(prec)
             print(a) # prints a with its coords in desired precision
        """
        # set precision of E, N
        self.strE = str(round(self.E, prec))
        self.strN = str(round(self.N, prec))
        # set precision of U
        if self.U is not None:
            self.strU = str(round(self.U, prec))
        else:
            self.strU = str(self.U)
        if self.geoidHeight is not None:
            self.strGeoidHeight = str(round(self.geoidHeight, prec))
        else:
            self.strGeoidHeight = str(self.geoidHeight)

        self.prec = prec

    def __repr__(self):
        # set default precision
        self.setPrec(self.prec)

        return "(" + ", ".join([self.strE, self.strN, self.strU, "ENU",
                                self.projection.name,
                                self.strGeoidHeight]) + ")"


def geodeticToECEF(point):
    """
    lat, lon and h: in radians and meters/None
    ellipsoid: class ellipsoid
    returns ECEF Coords X, Y, Z in meters
    """
    point.latlonh('r')
    ellipsoid = point.ellipsoid
    geoidHeight = point.geoidHeight
    lat = point.x
    lon = point.y
    h = point.z
    if h is None:
        h = 0

    W = sqrt(1 - ellipsoid.eSquared * (sin(lat))**2)
    N = ellipsoid.a / W
    X = (N + h) * cos(lat) * cos(lon)
    Y = (N + h) * cos(lat) * sin(lon)
    Z = ((1 - ellipsoid.eSquared) * N + h) * sin(lat)
    return GeoPoint(X, Y, Z, 'ECEF', ellipsoid, geoidHeight)


def ECEFtoGeodetic(point, coordType):
    """
    x, y, z in meters
    ellipsoid: class ellipsoid
    returns Geodetic lat, lon in radians (between -pi & pi) and h in meters
    """
    x = point.x
    y = point.y
    z = point.z
    ellipsoid = point.ellipsoid
    geoidHeight = point.geoidHeight

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

    gPt = GeoPoint(lat, lon, h, 'latlonh-r', ellipsoid, geoidHeight)
    gPt.latlonh(coordType)
    return gPt


def project(point, projection):
    """
    point: class GeoPoint
    projection: class Projection
    """
    try:
        if point.ellipsoid.name != projection.ellipsoid.name:
            raise Exception()
    except Exception as error:
        print("Error: ")
        print("Your projection type ellipsoid differs from the")
        print("original one")
    point.latlonh('r')
    lat, lon, h, geoidHeight = point.x, point.y, point.z, point.geoidHeight
    E, N = projection.directEqns(lat, lon, *projection.params[2:])
    try:
        U = h - geoidHeight
    except TypeError:
        U = None
    
    return MapPoint(E, N, U, projection, geoidHeight)


def unproject(point, ellipsoid, coordType):
    """
    point: class MspPoint
    ellipsoid: class Ellipsoid
    coordType: string - latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms
    """
    try:
        if point.projection.ellipsoid.name != ellipsoid.name:
            raise Exception()
    except Exception as error:
        print("Error: ")
        print("Your projection type ellipsoid differs from the")
        print("original one")
    E, N, U, geoidHeight = point.E, point.N, point.U, point.geoidHeight
    lat, lon = point.projection.invEqns(E, N, *point.projection.params[2:])
    try:
        h = U + geoidHeight
    except TypeError:
        h = None
    
    gPt = GeoPoint(lat, lon, h, 'latlonh-r', ellipsoid, geoidHeight)
    if coordType == 'ECEF':
        gPt.ECEF()
    elif coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
        gPt.latlonh(coordType.split('-')[-1])

    return gPt
