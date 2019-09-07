# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:27:10 2019

@author: Antonios Vatalis
"""
from math import sqrt, cos, sin, atan2
from angles import DMStoRads, radToDMS
from loadFiles import loadCoordsFile
from modelers import Ellipsoid, Projection
from copy import deepcopy


class GeoPoint():
    def __init__(self, ptID, x, y, z, coordType, ellipsoid, geoidHeight = None):
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
        self.id = str(ptID)
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
        self.coords = [self.x, self.y, self.z, self.geoidHeight]

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

    def __sub__(self, other):
        _coordType = self.coordType
        self.ECEF()
        other.ECEF()

        coords = []
        for x, y in zip(self.coords, other.coords):
            try:
                coords.append(x-y)
            except TypeError:
                coords.append(None)
        
        gPt = GeoPoint(self.id, coords[0], coords[1], coords[2],
                       self.coordType, self.projection, coords[3])
        if _coordType != 'ECEF':
            gPt.latlonh(_coordType.split('-')[1])
        return gPt
    
    def __repr__(self):
        # default precision of strX, strY, strZ
        self.setPrec(self.precs[self.coordType][0],
                     self.precs[self.coordType][1])

        return "(" + ", ".join([self.id, self.strX, self.strY, self.strZ,
                                self.coordType, self.ellipsoid.name,
                                self.strGeoidHeight]) + ")"


class MapPoint():
    def __init__(self, ptID, E, N, U, projection, geoidHeight = None):
        """
        This method initializes an object of class MapPoint

        E, N, U: map Eastings and Northings of point in meters
        U:       orthometric height of point in meters
        projection: object of class Projection where the MapPoint belongs to
        lat0, lon0: in radians. They matter only if projection belongs to 
                                    HATT family.
        geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
        """
        self.id = str(ptID)
        self.x = E
        self.y = N
        self.z = U
        self.projection = projection
        self.ellipsoid = self.projection.ellipsoid
        # geoid height
        self.geoidHeight = geoidHeight
        # default prec
        self.prec = 3
        self.coords = [self.x, self.y, self.z, self.geoidHeight]

    def setPrec(self, prec):
        """
        This method specifies the precision of attributes strEm, strN, strU

        prec: int

        e.g. a = MapPoint(E, N, U, projection)
             a.strPrec(prec)
             print(a) # prints a with its coords in desired precision
        """
        # set precision of E, N
        self.strX = str(round(self.x, prec))
        self.strY = str(round(self.y, prec))
        # set precision of U
        if self.z:
            self.strZ = str(round(self.z, prec))
        else:
            self.strZ = str(self.z)
        if self.geoidHeight:
            self.strGeoidHeight = str(round(self.geoidHeight, prec))
        else:
            self.strGeoidHeight = str(self.geoidHeight)

        self.prec = prec

    def __sub__(self, other):
        coords = []
        for x, y in zip(self.coords, other.coords):
            try:
                coords.append(x-y)
            except TypeError:
                coords.append(None)
        
        return MapPoint(self.id, coords[0], coords[1], coords[2], 
                        self.projection, coords[3])

    def __repr__(self):
        # set default precision
        self.setPrec(self.prec)

        return "(" + ", ".join([self.id, self.strX, self.strY, self.strZ, "ENU",
                                self.projection.name,
                                self.strGeoidHeight]) + ")"


def _geodeticToECEF(point):
    """
    lat, lon and h: in radians and meters/None
    ellipsoid: class ellipsoid
    returns ECEF Coords X, Y, Z in meters
    """
    _point = deepcopy(point)
    if point.coordType != 'ECEF':
        ptID = _point.id
        _point.latlonh('r')
        ellipsoid = _point.ellipsoid
        geoidHeight = _point.geoidHeight
        lat = _point.x
        lon = _point.y
        h = _point.z
        if h is None:
            h = 0
    
        W = sqrt(1 - ellipsoid.eSquared * (sin(lat))**2)
        N = ellipsoid.a / W
        X = (N + h) * cos(lat) * cos(lon)
        Y = (N + h) * cos(lat) * sin(lon)
        Z = ((1 - ellipsoid.eSquared) * N + h) * sin(lat)
        return GeoPoint(ptID, X, Y, Z, 'ECEF', ellipsoid, geoidHeight)
    else:
        return _point


def geodeticToECEF(points):
    if isinstance(points, list):
        return [_geodeticToECEF(point) for point in points]
    else:
        return _geodeticToECEF(points)


def _ECEFtoGeodetic(point, coordType):
    """
    x, y, z in meters
    ellipsoid: class ellipsoid
    returns Geodetic lat, lon in radians (between -pi & pi) and h in meters
    """
    _point = deepcopy(point)
    if point.coordType == 'ECEF':
        ptID = _point.id
        x = _point.x
        y = _point.y
        z = _point.z
        ellipsoid = _point.ellipsoid
        geoidHeight = _point.geoidHeight
    
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
    
        _point = GeoPoint(ptID, lat, lon, h, 'latlonh-r', ellipsoid, geoidHeight)
        _point.latlonh(coordType)
        return _point
    else:
        _point.latlonh(coordType)
        return _point


def ECEFtoGeodetic(points, coordType):
    if isinstance(points, list):
        return [_ECEFtoGeodetic(point, coordType) for point in points]
    else:
        return _ECEFtoGeodetic(points, coordType)


def _project(point, projection):
    """
    point: class GeoPoint
    projection: class Projection
    """
    _point = deepcopy(point)
    if isinstance(point, GeoPoint):
        try:
            if point.ellipsoid.name != projection.ellipsoid.name:
                raise Exception()
        except Exception as error:
            raise Exception("Error: Your projection type ellipsoid differs"
                            " from the original one")
        _point.latlonh('r')
        lat, lon, h, geoidHeight = _point.x, _point.y, _point.z, _point.geoidHeight
        E, N = projection.directEqns(lat, lon)
        try:
            U = h - geoidHeight
        except TypeError:
            U = 0
        
        return MapPoint(point.id, E, N, U, projection, geoidHeight)
    else:
        return _point


def project(points, projection):
    """
    point2: list of GeoPoints or single GeoPoint
    projection: class Projection
    """
    if isinstance(points, list):
        return [_project(point, projection) for point in points]
    else:
        return _project(points, projection)
            
            
def _unproject(point, coordType):
    """
    point: class MapPoint
    coordType: string - ECEF, latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms
    """
    if isinstance(point, MapPoint):
        E, N, U = point.x, point.y, point.z
        ellipsoid, geoidHeight = point.ellipsoid, point.geoidHeight
        lat, lon = point.projection.invEqns(E, N)
        try:
            h = U + geoidHeight
        except TypeError:
            h = 0
        _point = GeoPoint(point.id, lat, lon, h, 'latlonh-r', ellipsoid, geoidHeight)
    else:
        _point = deepcopy(point)

    if coordType == 'ECEF':
        _point.ECEF()
    elif coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm", "latlonh-dms"]:
        _point.latlonh(coordType.split('-')[-1])
    return _point


def unproject(points, coordType):
    """
    point: class MapPoint
    coordType: string - ECEF, latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms
    """
    if isinstance(points, list):
        return [_unproject(point, coordType) for point in points]
    else:
        return _unproject(points, coordType)


def fileDataToPoints(fileName, model, coordType = None, lastClmIsGeoid = True):
    """
     coordType: ECEF, latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms, ENU
    """
    if coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
        latlonIsStr = True
    else:
        latlonIsStr = False
    data = loadCoordsFile(fileName)  # list of lists
    points = []
    for sublist in data:
        if latlonIsStr:
            sublist = sublist[:3] + [float(x) if x else None for x in sublist[3:]]
        else:
            sublist = [sublist[0]] + [float(x) if x else None for x in sublist[1:]]

        if isinstance(model, Ellipsoid):
            args = sublist[:-1] + [coordType, model] + [sublist[-1]]
            points.append(GeoPoint(*args))
        elif isinstance(model, Projection):
            args = sublist[:-1] + [model] + [sublist[-1]]
            points.append(MapPoint(*args))

    if not lastClmIsGeoid:
        for point in points:
            point.geoidHeight = point.z - point.geoidHeight

    return points
