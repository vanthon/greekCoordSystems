# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:27:10 2019

@author: anthon
"""
from angles import DMStoRads, radToDMS
from ellipsoids import Ellipsoid, geodeticToECEF, ECEFtoGeodetic
from projections import Projection


class GeoPoint():
    """
    x, y, z in meters/radians/d/dm/dms
    coordType: latlonh-r, latlonh-dd, latlonh-dm, latlonh-dms, ECEF (meters)
    ellipsoid of class Ellipsoid
    """

    def __init__(self, x, y, z, coordType, ellipsoid):
        self.x = x
        self.y = y
        self.z = z
        self.coordType = coordType
        self.ellipsoid = ellipsoid

    def radians(self):
        if self.coordType in ["latlonh-dd", "latlonh-dm", "latlonh-dms"]:
            self.coordType = "latlonh-r"
            self.x = DMStoRads(self.x)
            self.y = DMStoRads(self.y)
            
    def latlonh(self, newCoordType):
        """
        newCoordType: r, dd, dm, dms
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
        
    def ECEF(self):
        # latlonh to ECEF
        if self.coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm",
                              "latlonh-dms"]:
            self.radians()
            self.x, self.y, self.z = \
                geodeticToECEF(self.x, self.y, self.z, self.ellipsoid)
        self.coordType = "ECEF"
    
    def mapPoint(self, projection, geoidHeight = None):
        """
        geoidHeight in meters
        projection of class Projection
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
            if geoidHeight != None and self.z != None:
                U = self.z - geoidHeight
            else:
                U = None
            return MapPoint(E, N, U, projection)

    def __repr__(self, precXY = 3, precZ = 3):
        if type(self.x) is not str and type(self.y) is not str:
            strX = str(round(self.x, precXY))
            strY = str(round(self.y, precXY))
        else:
            strX = self.x
            strY = self.y
            
        if self.z is not None: strZ = str(round(self.z, precZ))
        else: strZ = str(self.z)
        
        return ", ".join([strX, strY, strZ, self.coordType,
                          self.ellipsoid.name])


class MapPoint():
    """
    E, N, U in meters
    coordType: ENU
    projection of class Projection
    """
    
    def __init__(self, E, N, U, projection, lat0 = None, lon0 = None):
        self.E = E
        self.N = N
        self.U = U
        self.projection = projection
        self.ellipsoid = self.projection.ellipsoid
        # HATT make a try block
        self.lat0 = lat0
        self.lon0 = lon0
    
    def geoPoint(self, coordType, geoidHeight = None):
        """
        geoidHeight in meters
        returns Geopoint of coordType latlonh-r with height = 0
        """
        lat, lon = self.projection.fromMap(self.E, self.N)
        if geoidHeight != None and self.U != None:
            h = self.U + geoidHeight
        else:
            h = None
        
        # make a try block here
        gP = GeoPoint(lat, lon, h, "latlonh-r", self.ellipsoid)
        if coordType == "ECEF":
            gP.ECEF()
        elif coordType in ["latlonh-r", "latlonh-dd", "latlonh-dm", 
                           "latlonh-dms"]:
            gP.latlonh(coordType[8:])
        return gP
        
    def __repr__(self, precXY = 3, precZ = 3):
        strE = str(round(self.E, precXY)); strN = str(round(self.N, precXY))
            
        if self.U is not None: strU = str(round(self.U, precZ))
        else: strU = str(self.U)
        
        return ", ".join([strE, strN, strU, "ENU", self.projection.name])

        
        























