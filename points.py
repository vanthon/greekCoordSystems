# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:27:10 2019

@author: Antonios Vatalis
"""
from angles import DMStoRads, radToDMS
from ellipsoids import Ellipsoid, geodeticToECEF, ECEFtoGeodetic
from projections import Projection


class GeoPoint():
    def __init__(self, x, y, z, coordType, ellipsoid, geoidHeight = None):
        """
        This method initializes an object of class GeoPoint

        x, y, z: geocentric or geodetic coordinates of point
                 coordType: 'ECEF" -> x, y, z: in meters
                 coordType: 'latlonh-(r/dd/dm/dms)" -> x: lar in (radians/dd/dm/dms)
                                                       y: lon in (radians/dd/dm/dms)
                                                       z: h in meters
        coordType: one of the following strings:
                        'latlonh-r', 'latlonh-dd', 'latlonh-dm', 'latlonh-dms', 'ECEF'
        ellipsoid: object of class Ellipsoid
        geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
        """
        self.x = x
        self.y = y
        self.z = z
        self.coordType = coordType
        self.ellipsoid = ellipsoid
        self.geoidHeight = geoidHeight
        # strX, strY, strZ for __repr__
        self.strX = str(x)
        self.strY = str(y)
        self.strZ = str(z)
        # default precs
        if self.coordType == 'ECEF':
            self.precXY = 3
        elif self.coordType == 'latlonh-dms':
            self.precXY = 5
        else:
            self.precXY = 6
        self.precZ = 3

    def strXYZ(self, strx, stry, strz = None):
        """
        This method updates attributes strX, strY, strZ

        strx, stry, strz: strings
        """
        self.strX = strx
        self.strY = stry
        if strz is not None:
            self.strZ = strz

    def updPrecs(self, precXY, precZ = None):
        """
        This method updates attributes precXY, precZ

        precXY, precZ: ints
        """
        self.precXY = precXY
        if precZ is not None:
            self.precZ = precZ

    def radians(self):
        """
        This method transforms attributes x, y in radians
        (not if attribute coordType is 'ECEF')
        """
        if self.coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
            self.coordType = "latlonh-r"
            self.x = DMStoRads(self.x)
            self.y = DMStoRads(self.y)

        self.strXYZ(str(self.x), str(self.y))
        self.updPrecs(6)
            
    def latlonh(self, newCoordType):
        """
        This method transforms attributes x, y, z in geodetic coordinates: lat, lon, h
        The type of lat, lon is defined by newCoordType parameter

        newCoordType: one of the following strings: 'r', 'dd', 'dm', 'dms'
        """
        # latlonh to latlonh
        if self.coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
            self.radians()
        # ECEF to latlonh
        elif self.coordType == "ECEF":
            self.x, self.y, self.z =\
            ECEFtoGeodetic(self.x, self.y, self.z, self.ellipsoid) # radians
        
        self.coordType = "latlonh-" + newCoordType
        if newCoordType in ["dd", "dm", "dms"]:
            self.x = radToDMS(self.x, newCoordType)
            self.y = radToDMS(self.y, newCoordType)

        self.strXYZ(str(self.x), str(self.y), str(self.z))
        if self.coordType == 'latlonh-dms':
            self.updPrecs(5, 3)
        else:
            self.updPrecs(6, 3)
        
    def ECEF(self):
        """
        This method transforms attributes x, y, z in ECEF coordinates
        """
        # latlonh to ECEF
        if self.coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm",
                              "latlonh-dms"]:
            self.radians()
            self.x, self.y, self.z = \
                geodeticToECEF(self.x, self.y, self.z, self.ellipsoid)
        self.coordType = "ECEF"

        self.strXYZ(str(self.x), str(self.y), str(self.z))
        self.updPrecs(3, 3)
    
    def mapPoint(self, projection):
        """
        This method maps the GeoPoint to a specified projection returning
        a MapPoint

        projection: object of class Projection that specifies the desired
                    projection where the GeoPoint has to be mapped to
        geoidHeight: Geoid height of point in meters
        """
        try:
            if projection.ellipsoid.name != self.ellipsoid.name:
                raise Exception()
        except Exception as error:
            print("Error: ")
            print("Your projection type ellipsoid differs from the")
            print("original one")
        else:
            self.latlonh("r")
            E, N = projection.toMap(self.x, self.y)
            if self.geoidHeight is not None and self.z is not None:
                U = self.z - self.geoidHeight
            else:
                U = None
            return MapPoint(E, N, U, projection, self.geoidHeight)

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
            strX = radToDMS(DMStoRads(self.x), self.coordType[8:], precXY)
            strY = radToDMS(DMStoRads(self.y), self.coordType[8:], precXY)
        else:
            strX = str(round(self.x, precXY))
            strY = str(round(self.y, precXY))
        # set precision of z
        if self.z is not None:
            strZ = str(round(self.z, precZ))
        else:
            strZ = str(self.z)

        self.strXYZ(strX, strY, strZ)
        self.precXY = precXY
        self.precZ = precZ

    def __repr__(self):
        # default precision of strX, strY, strZ
        self.setPrec(self.precXY, self.precZ)

        return "(" + ", ".join([self.strX, self.strY, self.strZ, self.coordType,
                                self.ellipsoid.name]) + ")"


class MapPoint():
    def __init__(self, E, N, U, projection, geoidHeight = None):
        """
        This method initializes an object of class MapPoint

        E, N, U: map Eastings and Northings of point in meters
        U:       orthometric height of point in meters
        projection: object of class Projection where the MapPoint belongs to
        lat0, lon0: in radians. They matter only if projection belongs to HATT family.
        geoidHeight: geoid Height of Point in meters reffering to its ellipsoid
        """
        self.E = E
        self.N = N
        self.U = U
        self.projection = projection
        self.ellipsoid = self.projection.ellipsoid
        # geoid height
        self.geoidHeight = geoidHeight
        # strE, strN, strU for __repr__ method
        self.strE = str(E)
        self.strN = str(N)
        self.strU = str(U)
        # default prec
        self.precENU = 3

    def strENU(self, strE, strN, strU = None):
        """
        This method updates attributes strE, strN, strU

        strE, strN, strU: strings
        """
        # strX, strY, strZ for __repr__
        self.strE = strE
        self.strN = strN
        if strU is not None:
            self.strU = strU
    
    def geoPoint(self, coordType):
        """
        This method transforms the MapPoint to a GeoPoint and returns
        the latter

        coordType: Specifies coordinates format of returned GeoPoint.
                   It could be one o the following strings:
                    "ECEF", "latlonh-r", "latlonh-dd", "latlonh-dm", "latlonh-dms"
        geoidHeight: Geoid height of point in meters
        """
        lat, lon = self.projection.fromMap(self.E, self.N)
        if self.geoidHeight is not None and self.U is not None:
            h = self.U + self.geoidHeight
        else:
            h = None

        gPoint = GeoPoint(lat, lon, h, "latlonh-r", self.ellipsoid, self.geoidHeight)
        if coordType == "ECEF":
            gPoint.ECEF()
        elif coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm", 
                           "latlonh-dms"]:
            gPoint.latlonh(coordType[8:])
        return gPoint

    def setPrec(self, precENU):
        """
        This method specifies the precision of attributes strEm, strN, strU

        precENU: int

        e.g. a = MapPoint(E, N, U, projection)
             a.strPrec(precENU)
             print(a) # prints a with its coords in desired precision
        """
        # set precision of E, N
        strE = str(round(self.E, precENU))
        strN = str(round(self.N, precENU))
        # set precision of U
        if self.U is not None:
            strU = str(round(self.U, precENU))
        else:
            strU = str(self.U)

        self.strENU(strE, strN, strU)
        self.precENU = precENU

    def __repr__(self):
        # set default precision
        self.setPrec(self.precENU)

        return "(" + ", ".join([self.strE, self.strN, self.strU, "ENU",
                                self.projection.name]) + ")"
