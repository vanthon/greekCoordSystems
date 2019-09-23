# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:27:10 2019

@author: Antonios Vatalis
"""
from math import sqrt, cos, sin, atan2
from angles import DMStoRads, radToDMS, toDMS, toRads, Angle
from loadFiles import loadCoordsFile
from modelers import Ellipsoid, Projection
from copy import deepcopy

# class GeoPoint:
#     def __init__(self, ptID, x, y, z, coordType, ellipsoid,
#                  geoidHeight=None, GMHeight=None):
#         """
#         This method initializes an object of class GeoPoint
#
#         x, y, z: geocentric or geodetic coordinates of point
#                  coordType: 'ECEF" -> x, y, z: in meters
#                  coordType: 'latlonh-(r/dd/dm/dms)"
#                              x: lar in (radians/dd/dm/dms)
#                              y: lon in (radians/dd/dm/dms)
#                              z: h in meters
#         coordType: one of the following strings:
#                         'latlonh-r', 'latlonh-dd', 'latlonh-dm', 'latlonh-dms',
#                         'ECEF'
#         ellipsoid: object of class Ellipsoid
#         geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
#         """
#         self.id = str(ptID)
#         self.x = x
#         self.y = y
#         self.z = z
#         self.coordType = coordType
#         self.ellipsoid = ellipsoid
#         self.geoidHeight = geoidHeight
#         self.GMHeight = GMHeight
#         # default precs
#         self.precs = {'ECEF': [3, 3],
#                  'latlonh-r': [6, 3],
#                  'latlonh-dd': [6, 3],
#                  'latlonh-dm': [6, 3],
#                  'latlonh-dms': [5, 3]}
#         self.coords = [self.x, self.y, self.z, self.geoidHeight, self.GMHeight]
#
#     def radians(self):
#         """
#         This method transforms attributes x, y in radians
#         (not if attribute coordType is 'ECEF')
#         """
#         if self.coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
#             self.coordType = "latlonh-r"
#             self.x = DMStoRads(self.x)
#             self.y = DMStoRads(self.y)
#
#     def latlonh(self, newCoordType):
#         """
#         This method transforms attributes x, y, z in geodetic coordinates:
#             lat, lon, h
#         The type of lat, lon is defined by newCoordType parameter
#
#         newCoordType: one of the following strings: 'r', 'dd', 'dm', 'dms'
#         """
#         # latlonh to latlonh
#         if self.coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm",
#                               "latlonh-dms"]:
#             self.radians()
#             if newCoordType != 'r':
#                 self.x = radToDMS(self.x, newCoordType)
#                 self.y = radToDMS(self.y, newCoordType)
#         # ECEF to latlonh
#         elif self.coordType == "ECEF":
#             geodeticPoint = ECEFtoGeodetic(self, newCoordType)
#             self.x, self.y, self.z =\
#                 geodeticPoint.x, geodeticPoint.y, geodeticPoint.z
#
#         self.coordType = "latlonh-" + newCoordType
#
#     def ECEF(self):
#         """
#         This method transforms attributes x, y, z in ECEF coordinates
#         """
#         # latlonh to ECEF
#         if self.coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm",
#                               "latlonh-dms"]:
#             self.radians()
#             ecefPoint = geodeticToECEF(self)
#             self.x, self.y, self.z = ecefPoint.x, ecefPoint.y, ecefPoint.z
#         self.coordType = "ECEF"
#
#     def setPrec(self, precXY, precZ):
#         """
#         This method specifies the precision of attributes strX, strY, strZ
#
#         precXY: int, precision of attributes strX, strY
#         precZ:  int, precision of attribute strZ
#
#         e.g. a = GeoPoint(x, y, z, coordType, ellipsoid)
#              a.strPrec(precXY, precZ)
#              print(a) # prints a with its coords in desired precisions
#         """
#         # set precision of x and y
#         if type(self.x) is str:
#             self.strX = radToDMS(DMStoRads(self.x),
#                                  self.coordType.split('-')[1], precXY)
#             self.strY = radToDMS(DMStoRads(self.y),
#                                  self.coordType.split('-')[1], precXY)
#         else:
#             self.strX = str(round(self.x, precXY))
#             self.strY = str(round(self.y, precXY))
#         # set precision of z
#         if self.z is not None:
#             self.strZ = str(round(self.z, precZ))
#         else:
#             self.strZ = str(self.z)
#         if self.geoidHeight:
#             self.strGeoidHeight = str(round(self.geoidHeight, precZ))
#         else:
#             self.strGeoidHeight = str(self.geoidHeight)
#         if self.GMHeight:
#             self.strGMHeight = str(round(self.GMHeight, precZ))
#         else:
#             self.strGMHeight = str(self.GMHeight)
#
#         self.precs[self.coordType] = [precXY, precZ]
#
#     def __sub__(self, other):
#         _coordType = self.coordType
#         self.ECEF()
#         other.ECEF()
#
#         coords = []
#         for x, y in zip(self.coords, other.coords):
#             try:
#                 coords.append(x-y)
#             except TypeError:
#                 coords.append(None)
#
#         gPt = GeoPoint(self.id, coords[0], coords[1], coords[2],
#                        self.coordType, self.projection, coords[3], coords[4])
#         if _coordType != 'ECEF':
#             gPt.latlonh(_coordType.split('-')[1])
#         return gPt
#
#     def __repr__(self):
#         # default precision of strX, strY, strZ
#         self.setPrec(self.precs[self.coordType][0],
#                      self.precs[self.coordType][1])
#
#         return "(" + ", ".join([self.id, self.strX, self.strY, self.strZ,
#                                 self.coordType, self.ellipsoid.name,
#                                 self.strGeoidHeight, self.strGMHeight]) + ")"


# class MapPoint:
#     def __init__(self, ptID, E, N, U, projection,
#                  geoidHeight=None, GMHeight=None):
#         """
#         This method initializes an object of class MapPoint
#
#         E, N, U: map Eastings and Northings of point in meters
#         U:       orthometric height of point in meters
#         projection: object of class Projection where the MapPoint belongs to
#         lat0, lon0: in radians. They matter only if projection belongs to
#                                     HATT family.
#         geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
#         """
#         self.id = str(ptID)
#         self.x = E
#         self.y = N
#         self.z = U
#         self.projection = projection
#         self.ellipsoid = self.projection.ellipsoid
#         self.geoidHeight = geoidHeight
#         self.GMHeight = GMHeight
#         # default prec
#         self.prec = 3
#         self.coords = [self.x, self.y, self.z, self.geoidHeight, self.GMHeight]
#
#     def setPrec(self, prec):
#         """
#         This method specifies the precision of attributes strEm, strN, strU
#
#         prec: int
#
#         e.g. a = MapPoint(E, N, U, projection)
#              a.strPrec(prec)
#              print(a) # prints a with its coords in desired precision
#         """
#         # set precision of E, N
#         self.strX = str(round(self.x, prec))
#         self.strY = str(round(self.y, prec))
#         # set precision of U
#         if self.z:
#             self.strZ = str(round(self.z, prec))
#         else:
#             self.strZ = str(self.z)
#         if self.geoidHeight:
#             self.strGeoidHeight = str(round(self.geoidHeight, prec))
#         else:
#             self.strGeoidHeight = str(self.geoidHeight)
#         if self.GMHeight:
#             self.strGMHeight = str(round(self.GMHeight, prec))
#         else:
#             self.strGMHeight = str(self.GMHeight)
#
#         self.prec = prec
#
#     def __sub__(self, other):
#         coords = []
#         for x, y in zip(self.coords, other.coords):
#             try:
#                 coords.append(x-y)
#             except TypeError:
#                 coords.append(None)
#
#         return MapPoint(self.id, coords[0], coords[1], coords[2],
#                         self.projection, coords[3], coords[4])
#
#     def __repr__(self):
#         # set default precision
#         self.setPrec(self.prec)
#
#         return "(" + ", ".join([self.id, self.strX, self.strY, self.strZ,
#                                 "ENU", self.projection.name,
#                                 self.strGeoidHeight, self.strGMHeight]) + ")"


# def _geodeticToECEF(point):
#     """
#     lat, lon and h: in radians and meters/None
#     ellipsoid: class ellipsoid
#     returns ECEF Coords X, Y, Z in meters
#     """
#     _point = deepcopy(point)
#     if point.coordType != 'ECEF':
#         ptID = _point.id
#         _point.latlonh('r')
#         ellipsoid = _point.ellipsoid
#         geoidHeight = _point.geoidHeight
#         lat = _point.x
#         lon = _point.y
#         h = _point.z
#         if h is None:
#             h = 0
#
#         W = sqrt(1 - ellipsoid.eSquared * (sin(lat))**2)
#         N = ellipsoid.a / W
#         X = (N + h) * cos(lat) * cos(lon)
#         Y = (N + h) * cos(lat) * sin(lon)
#         Z = ((1 - ellipsoid.eSquared) * N + h) * sin(lat)
#         return GeoPoint(ptID, X, Y, Z, 'ECEF', ellipsoid, geoidHeight)
#     else:
#         return _point


# def geodeticToECEF(points):
#     if isinstance(points, list):
#         return [_geodeticToECEF(point) for point in points]
#     else:
#         return _geodeticToECEF(points)


# def _ECEFtoGeodetic(point, coordType):
#     """
#     x, y, z in meters
#     ellipsoid: class ellipsoid
#     returns Geodetic lat, lon in radians (between -pi & pi) and h in meters
#     """
#     _point = deepcopy(point)
#     if point.coordType == 'ECEF':
#         ptID = _point.id
#         x = _point.x
#         y = _point.y
#         z = _point.z
#         ellipsoid = _point.ellipsoid
#         geoidHeight = _point.geoidHeight
#
#         lon = atan2(y, x)
#         P = sqrt(x**2 + y**2)
#         lat0 = atan2((z * (1 + ellipsoid.ePrimeSquared)), P)
#         dlat = 1
#         while dlat > 10**-12:
#             W = sqrt(1 - ellipsoid.eSquared * (sin(lat0))**2)
#             N = ellipsoid.a / W
#             lat = atan2((z + ellipsoid.eSquared * N * sin(lat0)), P)
#             dlat = abs(lat - lat0)
#             lat0 = lat
#         h = (z / sin(lat)) - ((1 - ellipsoid.eSquared) * N)
#
#         _point = GeoPoint(ptID, lat, lon, h, 'latlonh-r',
#                           ellipsoid, geoidHeight)
#         _point.latlonh(coordType)
#         return _point
#     else:
#         _point.latlonh(coordType)
#         return _point


# def ECEFtoGeodetic(points, coordType):
#     if isinstance(points, list):
#         return [_ECEFtoGeodetic(point, coordType) for point in points]
#     else:
#         return _ECEFtoGeodetic(points, coordType)


# def _project(point, projection):
#     """
#     point: class GeoPoint
#     projection: class Projection
#     """
#     _point = deepcopy(point)
#     if isinstance(point, GeoPoint):
#         try:
#             if point.ellipsoid.name != projection.ellipsoid.name:
#                 raise Exception()
#         except Exception as error:
#             raise Exception("Error: Your projection type ellipsoid differs"
#                             " from the original one")
#         _point.latlonh('r')
#         lat, lon, h = _point.x, _point.y, _point.z
#         geoidHeight, GMHeight = _point.geoidHeight, _point.GMHeight
#         E, N = projection.directEqns(lat, lon)
#         try:
#             U = h - geoidHeight
#         except TypeError:
#             U = 0
#
#         return MapPoint(point.id, E, N, U, projection, geoidHeight, GMHeight)
#     else:
#         return _point


# def project(points, projection):
#     """
#     point2: list of GeoPoints or single GeoPoint
#     projection: class Projection
#     """
#     if isinstance(points, list):
#         return [_project(point, projection) for point in points]
#     else:
#         return _project(points, projection)
#
#
# def _unproject(point, coordType):
#     """
#     point: class MapPoint
#     coordType: string - ECEF, latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms
#     """
#     if isinstance(point, MapPoint):
#         E, N, U = point.x, point.y, point.z
#         ellipsoid, geoidHeight = point.ellipsoid, point.geoidHeight
#         GMHeight = point.GMHeight
#         lat, lon = point.projection.invEqns(E, N)
#         try:
#             h = U + geoidHeight
#         except TypeError:
#             h = 0
#         _point = GeoPoint(point.id, lat, lon, h, 'latlonh-r', ellipsoid,
#                           geoidHeight, GMHeight)
#     else:
#         _point = deepcopy(point)
#
#     if coordType == 'ECEF':
#         _point.ECEF()
#     elif coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm", "latlonh-dms"]:
#         _point.latlonh(coordType.split('-')[-1])
#     return _point
#
#
# def unproject(points, coordType):
#     """
#     point: class MapPoint
#     coordType: string - ECEF, latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms
#     """
#     if isinstance(points, list):
#         return [_unproject(point, coordType) for point in points]
#     else:
#         return _unproject(points, coordType)


#########################################################################################

# class ECEFpointS:
#     def __init__(self, points):
#         if not all(isinstance(x, ECEFpoint) for x in points):
#             raise Exception('Your points list should consist of ECEFpoints')
#         self.points = points  # list of points
#         self.coordSystem = points[0].coordSystem
#
#     def geodetic(self):
#         return GeoPointS([x.geodetic() for x in self.points])
#
#     def project(self, projection):
#         return MapPointS([x.project(projection) for x in self.points])
#
#     def __repr__(self):
#         return '\n'.join(map(str, self.points))
#
#
# class GeoPointS:
#     def __init__(self, points):
#         if not all(isinstance(x, GeoPoint) for x in points):
#             raise Exception('Your points list should consist of GeoPoints')
#         self.points = points  # list of points
#         self.coordSystem = points[0].coordSystem
#
#     def ECEF(self):
#         return ECEFpointS([x.ECEF() for x in self.points])
#
#     def project(self, projection):
#         return MapPointS([x.project(projection) for x in self.points])
#
#     def __repr__(self):
#         return '\n'.join(map(str, self.points))
#
#
# class MapPointS:
#     def __init__(self, points):
#         if not all(isinstance(x, MapPoint) for x in points):
#             raise Exception('Your points list should consist of MapPoints')
#         self.points = points  # list of points
#         self.coordSystem = points[0].coordSystem
#
#     def geodetic(self):
#         return GeoPointS([x.geodetic() for x in self.points])
#
#     def ECEF(self):
#         return ECEFpointS([x.ECEF() for x in self.points])
#
#     def __repr__(self):
#         return '\n'.join(map(str, self.points))


def _geodeticToECEF(lat, lon, h, ellipsoid):
    """
    lat, lon and h: in Angles and meters/None
    ellipsoid: class ellipsoid
    returns ECEF Coords X, Y, Z in meters
    """
    lat = lat.radians
    lon = lon.radians
    if h is None:
        h = 0

    W = sqrt(1 - ellipsoid.eSquared * (sin(lat)) ** 2)
    N = ellipsoid.a / W
    X = (N + h) * cos(lat) * cos(lon)
    Y = (N + h) * cos(lat) * sin(lon)
    Z = ((1 - ellipsoid.eSquared) * N + h) * sin(lat)
    return X, Y, Z


def _ECEFtoGeodetic(X, Y, Z, ellipsoid):
    """
    X, Y, Z in meters
    ellipsoid: class ellipsoid
    returns Geodetic lat, lon in Angle and h in meters
    """
    lon = atan2(Y, X)
    P = sqrt(X ** 2 + Y ** 2)
    lat0 = atan2((Z * (1 + ellipsoid.ePrimeSquared)), P)
    dlat = 1
    while dlat > 10 ** -12:
        W = sqrt(1 - ellipsoid.eSquared * (sin(lat0)) ** 2)
        N = ellipsoid.a / W
        lat = atan2((Z + ellipsoid.eSquared * N * sin(lat0)), P)
        dlat = abs(lat - lat0)
        lat0 = lat
    h = (Z / sin(lat)) - ((1 - ellipsoid.eSquared) * N)

    lat = Angle(lat)
    lat.angleType = 'dms'
    lon = Angle(lon)
    lon.angleType = 'dms'
    return lat, lon, h


class ECEFpoint:
    def __init__(self, ptID=None, X=None, Y=None, Z=None, ellipsoid=None,
                 geoidHeight=None, GMHeight=None):
        self.id = ptID
        self.X = X
        self.Y = Y
        self.Z = Z
        self.geoidHeight = geoidHeight
        self.GMHeight = GMHeight
        self.coordSystem = ellipsoid
        self.prec = 3

    @property
    def U(self):
        point = self.geodetic()
        try:
            return point.h - self.geoidHeight
        except TypeError:
            return None

    @U.setter
    def U(self, value):
        point = self.geodetic()
        try:
            self.geoidHeight = point.h - value
        except TypeError:
            self.geoidHeight = None

    @property
    def texts(self):
        texts = {}
        keys = ['id', 'X', 'Y', 'Z', 'coordSystem', 'geoidHeight', 'GmHeight']
        values = [self.id, self.X, self.Y, self.Z, self.coordSystem,
                  self.geoidHeight, self.GMHeight]
        for key, value in zip(keys, values):
            try:
                texts[key] = str(round(value, self.prec))
            except TypeError:
                texts[key] = str(value)
        return texts

    def geodetic(self):
        lat, lon, h = _ECEFtoGeodetic(self.X, self.Y, self.Z,
                                      self.coordSystem)
        pt = GeoPoint(self.id, lat, lon, h, self.coordSystem,
                      self.geoidHeight, self.GMHeight)
        pt.U = self.U
        return pt

    def project(self, projection):
        if projection.ellipsoid.name != self.coordSystem.name:
            raise Exception('Define a projection with same Ellipsoid as'
                            'your point uses')
        point = self.geodetic()
        E, N = projection.directEqns(point.lat, point.lon)
        pt = MapPoint(self.id, E, N, self.U, projection,
                      self.geoidHeight, self.GMHeight)
        pt.h = point.h
        return pt

    def transform(self, transformation):
        return transformation.transform(self)

    def __str__(self):
        return ", ".join(self.texts.values())

    def __repr__(self):
        return "(" + str(self) + ")"


class GeoPoint:
    def __init__(self, ptID=None, lat=None, lon=None, h=None, ellipsoid=None,
                 geoidHeight=None, GMHeight=None):
        """
        lat, lon of type Angle
        """
        self.id = ptID
        self.lat = lat
        self.lon = lon
        self.h = h
        self.geoidHeight = geoidHeight
        self.GMHeight = GMHeight
        self.coordSystem = ellipsoid
        self.precLatLon = 5
        self.precHeights = 3

    @property
    def latLonType(self):
        return self.lat.angleType

    @latLonType.setter
    def latLonType(self, newType):
        self.lat.angleType = newType
        self.lon.angleType = newType

    @property
    def U(self):
        try:
            return self.h - self.geoidHeight
        except TypeError:
            return None

    @U.setter
    def U(self, value):
        try:
            self.geoidHeight = self.h - value
        except TypeError:
            self.geoidHeight = None

    @property
    def texts(self):
        self.lat.prec = self.lon.prec = self.precLatLon

        texts = {}
        keys = ['id', 'lat', 'lon', 'h', 'coordSystem', 'geoidHeight', 'GmHeight']
        values = [self.id, self.lat, self.lon, self.h, self.coordSystem,
                  self.geoidHeight, self.GMHeight]
        for key, value in zip(keys, values):
            try:
                texts[key] = str(round(value, self.precHeights))
            except TypeError:
                texts[key] = str(value)
        return texts

    def ECEF(self):
        X, Y, Z = _geodeticToECEF(self.lat, self.lon,
                                  self.h, self.coordSystem)
        pt = ECEFpoint(self.id, X, Y, Z, self.coordSystem,
                       self.geoidHeight, self.GMHeight)
        pt.U = self.U
        return pt

    def project(self, projection):
        if projection.ellipsoid.name != self.coordSystem.name:
            raise Exception('Define a projection with same Ellipsoid as'
                            'your point uses')
        E, N = projection.directEqns(self.lat, self.lon)
        pt = MapPoint(self.id, E, N, self.U, projection,
                      self.geoidHeight, self.GMHeight)
        pt.h = self.h
        return pt

    def transform(self, transformation):
        return transformation.transform(self)

    def __str__(self):
        return ", ".join(self.texts.values())

    def __repr__(self):
        return "(" + str(self) + ")"


class MapPoint:
    def __init__(self, ptID=None, E=None, N=None, U=None, projection=None,
                 geoidHeight=None, GMHeight=None):
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
        self.E = E
        self.N = N
        self.U = U
        self.coordSystem = projection
        # self.ellipsoid = self.coordSystem.ellipsoid
        self.geoidHeight = geoidHeight
        self.GMHeight = GMHeight
        self.prec = 3

    @property
    def h(self):
        try:
            return self.U + self.geoidHeight
        except TypeError:
            return None

    @h.setter
    def h(self, value):
        try:
            self.geoidHeight = value - self.U
        except TypeError:
            self.geoidHeight = None

    @property
    def texts(self):
        texts = {}
        keys = ['id', 'E', 'N', 'U', 'coordSystem', 'geoidHeight', 'GmHeight']
        values = [self.id, self.E, self.N, self.U, self.coordSystem,
                  self.geoidHeight, self.GMHeight]
        for key, value in zip(keys, values):
            try:
                texts[key] = str(round(value, self.prec))
            except TypeError:
                texts[key] = str(value)
        return texts

    def geodetic(self):
        lat, lon = self.coordSystem.invEqns(self.E, self.N)
        pt = GeoPoint(self.id, lat, lon, self.h, self.coordSystem.ellipsoid,
                      self.geoidHeight, self.GMHeight)
        pt.U = self.U
        return pt

    def ECEF(self):
        pt = self.geodetic()
        pt = pt.ECEF()
        pt.U = self.U
        return pt

    def transform(self, transformation):
        return transformation.transform(self)

    def __str__(self):
        return ", ".join(self.texts.values())

    def __repr__(self):
        return "(" + str(self) + ")"


class Points:
    """
    Points: class of GeoPoints, ECEFpoints or MapPoints
    returns its corresponding list of ECEFpoints
    """
    def __init__(self, points):
        self.points = points  # list of points
        if not any([self.all(GeoPoint), self.all(ECEFpoint), self.all(MapPoint)]):
            raise Exception('Your points list should consist of '
                            'GeoPoints, ECEFpoints or MapPoints')
        if not all([point.coordSystem is points[0].coordSystem for point in points]):
            raise Exception('All points within list should be on the same '
                            'coordinate system')

        self.coordSystem = points[0].coordSystem

    def all(self, pointType):
        return all([isinstance(x, pointType) for x in self.points])

    @property
    def elmsType(self):
        if self.all(GeoPoint):
            return GeoPoint
        elif self.all(ECEFpoint):
            return ECEFpoint
        elif self.all(MapPoint):
            return MapPoint
        else:
            raise Exception('Your points list should consist of '
                            'GeoPoints, ECEFpoints or MapPoints')

    @property
    def geoidHeight(self):
        return [x.geoidHeight for x in self]

    @geoidHeight.setter
    def geoidHeight(self, geoids: list):
        for a, b in zip(self, geoids):
            a.geoidHeight = b

    @property
    def GMHeight(self):
        return [x.GMHeight for x in self]

    @GMHeight.setter
    def GMHeight(self, heights):
        for a, b in zip(self, heights):
            a.GMHeight = b

    def ECEF(self):
        if self.elmsType is not ECEFpoint:
            return Points([x.ECEF() for x in self])
        else:
            return deepcopy(self)

    def geodetic(self):
        if self.elmsType is not GeoPoint:
            return Points([x.geodetic() for x in self])
        else:
            return deepcopy(self)

    def project(self, projection):
        if self.elmsType is not MapPoint:
            return Points([x.project(projection) for x in self])
        else:
            return deepcopy(self)

    def transform(self, transformation):
        return Points([x.transform(transformation) for x in self])

    def __add__(self, other):
        return Points(self.points + other.points)

    def __len__(self):
        return len(self.points)

    def __iter__(self):
        for point in self.points:
            yield point

    def __repr__(self):
        return '\n'.join(map(str, self))


def importPoints(fileName, clmsSpecification, pointType, coordSystem, dlm=','):
    """"""
    if {'h', 'U', 'geoidHeight'}.issubset(set(clmsSpecification)):
        raise Exception('Two heights among [h, U, geoidHeight] could be present '
                        'within your columns specification')

    data = loadCoordsFile(fileName, clmsSpecification, dlm=dlm)

    points = []
    for i in range(len(list(data.values())[0])):
        if pointType is ECEFpoint:
            point = ECEFpoint()
        elif pointType is GeoPoint:
            point = GeoPoint()
        elif pointType is MapPoint:
            point = MapPoint()
        else:
            raise Exception('Define properly the type of imported points')
        for key in data.keys():
            if key in ['lat-r', 'lon-r', 'lat-dd', 'lon-dd', 'lat-dm',
                       'lon-dm', 'lat-dms', 'lon-dms']:
                setattr(point, key.split('-')[0], data[key][i])
            else:
                setattr(point, key, data[key][i])
            setattr(point, 'coordSystem', coordSystem)
        points.append(point)

    return Points(points)
