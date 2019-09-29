# -*- coding: utf-8 -*-
"""
@author: Antonios Vatalis
"""
from math import sqrt, cos, sin, atan2
from angles import Angle
from loadFiles import loadCoordsFile
from copy import deepcopy
from geoid import GeoidHeight


def geodeticToECEF(lat, lon, h, ellipsoid):
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


def ECEFtoGeodetic(X, Y, Z, ellipsoid):
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
        # point = self.geodetic()
        lat, lon, h = ECEFtoGeodetic(self.X, self.Y, self.Z,
                                     self.coordSystem)
        try:
            return h - self.geoidHeight
        except TypeError:
            return None

    @U.setter
    def U(self, value):
        # point = self.geodetic()
        lat, lon, h = ECEFtoGeodetic(self.X, self.Y, self.Z,
                                     self.coordSystem)
        try:
            self.geoidHeight = h - value
        except TypeError:
            self.geoidHeight = None

    def setGMHeight(self):
        if self.coordSystem.name != 'WGS84':
            raise Exception('The ellipsoid should be WGS84')
        pt = self.geodetic()
        lat = float(pt.lat.dd)
        lon = float(pt.lon.dd)
        gmHeight = GeoidHeight()
        self.GMHeight = gmHeight.get(lat, lon)

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
        lat, lon, h = ECEFtoGeodetic(self.X, self.Y, self.Z,
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

    def setGMHeight(self):
        if self.coordSystem.name != 'WGS84':
            raise Exception('The ellipsoid should be WGS84')
        lat = float(self.lat.dd)
        lon = float(self.lon.dd)
        gmHeight = GeoidHeight()
        self.GMHeight = gmHeight.get(lat, lon)

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
        X, Y, Z = geodeticToECEF(self.lat, self.lon,
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

    def setGMHeight(self):
        if self.coordSystem.ellipsoid.name != 'WGS84':
            raise Exception('The ellipsoid should WGS84')
        pt = self.geodetic()
        lat = float(pt.lat.dd)
        lon = float(pt.lon.dd)
        gmHeight = GeoidHeight()
        self.GMHeight = gmHeight.get(lat, lon)

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
