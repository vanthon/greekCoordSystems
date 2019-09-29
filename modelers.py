# -*- coding: utf-8 -*-
"""
@author: Antonios Vatalis
"""
from importlib import import_module
from loadFiles import loadEllipsoidParams, loadProjParams
from csv import reader


class Ellipsoid:
    """
    Creates Ellipsoid objects
    """
    def __init__(self, ellName):
        """
        Arguments
            ellName: string that would be used as argument in order
                     to load refEllipsoids.txt file and filter its lines.
                     If refEllipsoids.txt file contains some ellipsoid
                     where its id is identical with ellname, function
                     creates Ellipsoid object, otherwise an error occures.
        """
        params = loadEllipsoidParams(ellName)
        self.name = params['ellipsID']
        self.a = params['semiMajorAxis']
        self.b = params['semiMinorAxis']
        self.f = params['flattening']
        self.eSquared = params['eccentricitySquared']
        self.ePrimeSquared = params['eccentricityPrimeSquared']

    def __str__(self):
        return "Ellipsoid: " + self.name

    def __repr__(self):
        return str(self)


class Projection:
    """
    Creates Projection objects
    """
    def __init__(self, projName, lat0=None, lon0=None):
        """
        Arguments
            projName: string that would be used as argument in order
                      to load projections.txt file and filter its lines.
                      If projections.txt file contains some projection
                      where its id is identical with projName, function
                      creates Ellipsoid object, otherwise an error occures.
            lat0, lon0: numbers or strings (They should be both of same type).
                        If they are numbers their values treated as radians.
                        If they are strings their values treated as degrees
                        depending their format, i.e. dd or dm or dms, e.g.
                        -0:3422332, -0:59.324432, 12:43:09.78349 respectively.
                        They used in order to overload lat0, lon0 within
                        projections.txt file.
        """     
        params = loadProjParams(projName)
        self.name = params['projID']
        self.family = params['projFamily']
        self.ellipsoid = Ellipsoid(params['ellipsoid'])
        self.m0 = params['scale']     
        self.E0 = params['falseEasting']
        self.N0 = params['falseNorthing']
        self.lat0 = params['lat0']
        self.lon0 = params['lon0']
        # overload lat0, lon0 (concerns HATT projection)
        if all([lat0, lon0]):
            self.lat0 = lat0
            self.lon0 = lon0
        # load projection direct and inverse eqn's
        eqnsFamily = import_module(self.family + "eqns")
        self._directEqns = getattr(eqnsFamily, "direct" + self.family)
        self._invEqns = getattr(eqnsFamily, "inv" + self.family)
    
    def directEqns(self, lat, lon):
        """
        This method projects lat, lon geodetic coordinates to map coordinates.
        Arguments
            lat, lon: radians or strings (both of same type)
                      if number: angle in radians
                      if string: angle in degrees of the following formats
                                 decimal degrees (dd), e.g. -22.0912233
                                 degrees:minutes (dm), e.g. 34:23.45902
                                 degrees:minutes:seconds (dms), e.g. 0:1:9.082
        Returns
            E, N coordinates in meters
        """
        return self._directEqns(lat, lon, self)
    
    def invEqns(self, E, N):
        """
        This method converts E, N map coordinates to geodetic coordinates.
        Arguments
            E, N: numbers as map coordinates in meters
            latlonType: string. One of the followings:
                            r: for radians
                            dd: for decimal degrees, e.g. -0.342423
                            dm: for degrees:minutes, e.g. -0:3.321
                            dms: for degrees:minutes:seconds, e.g. 1:9:32.4324
        Returns
            lat, lon tuple in radinas or degrees depending latlonType
        """
        lat, lon = self._invEqns(E, N, self)  # Angles
        return lat, lon

    def __str__(self):
        return ", ".join(["Projection: " + self.name,
                          "Ellipsoid: " + self.ellipsoid.name])

    def __repr__(self):
        return str(self)


def filterProj(ellipsoid):
    """
    This function loads projections.txt file, filters it based on given
    Ellipsoid object and returns its first corresponding Projection object.
    Arguments
        ellipsoid: object of Ellipsoid class. ellipsoid.name represents
                   the third element of every row (header excluded) in
                   projections.txt file.
    Returns
        object of type Projection if ellipsoid.name exists in projections.txt
        file, otherwise Exception raises
    """
    with open("projections.txt") as fObj:
        readCSV = reader(fObj, delimiter=',')
        for line in readCSV:
            if line[2].strip() == ellipsoid.name:
                return Projection(line[0].strip())
        raise Exception("The ellipsoid you defined does not exist")
