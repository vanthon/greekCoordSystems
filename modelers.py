# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:25:03 2019

@author: anthon
"""
from importlib import import_module
from loadFiles import loadEllipsoidParams, loadProjParams
from angles import DMStoRads
from csv import reader


class Ellipsoid():
    """
    Ellipsoid class
    """
    def __init__(self, ellName):
        """
        ellName: string, 'WGS84', 'GRS80', 'BESSEL', 'HAYFORD'
        """
        params = loadEllipsoidParams(ellName)
        self.name = params['ellipsID']
        self.a = float(params['semiMajorAxis'])
        self.b = float(params['semiMinorAxis'])
        self.f = float(params['flattening'])
        self.eSquared = float(params['eccentricitySquared'])
        self.ePrimeSquared = float(params['eccentricityPrimeSquared'])

    def __repr__(self):
        return "EllipsoidID: " + self.name


class Projection():
    def __init__(self, projName, lat0 = None, lon0 = None):
        """
        projName string: "HATT", "TM3-WEST", "TM3-CENTRAL", "TM3-EAST", "TM87",
                         "UTM-34N", "UTM-35N"
        lat0, lon0: numeral for radians
                    string for dd, dm, dms. e.g. '23.33', '3:34.2', '-0:3:4.33'
        """     
        params = loadProjParams(projName)
        self.name = params['projID']
        self.family = params['projFamily']
        self.ellipsoid = Ellipsoid(params['ellipsoid'])
        self.m0 = params['scale']     
        self.E0 = params['falseEasting']
        self.N0 = params['falseNorthing']
        try:
            self.m0 = float(self.m0)       
            self.E0 = float(self.E0)
            self.N0 = float(self.N0)
        except TypeError:
            pass
        
        self.lat0 = params['lat0']
        self.lon0 = params['lon0']
        if all([self.lat0, self.lon0]):
            self.lat0 = DMStoRads(self.lat0)
            self.lon0 = DMStoRads(self.lon0)
        
        if all([lat0, lon0]):
            if all([isinstance(lat0, str), isinstance(lon0, str)]):
                lat0 = DMStoRads(lat0)
                lon0 = DMStoRads(lon0)
            self.lat0 = lat0
            self.lon0 = lon0

        eqnsFamily = import_module(self.family + "eqns")
        self._directEqns = getattr(eqnsFamily, "direct" + self.family)
        self._invEqns = getattr(eqnsFamily, "inv" + self.family)
    
    def directEqns(self, lat, lon):
        return self._directEqns(lat, lon, self)  # returns E, N in meters
    
    def invEqns(self, E, N):
        return self._invEqns(E, N, self)  #returns lat, lon in radians
    
    def __repr__(self):
        return ", ".join(["projID: " + self.name, "projFamily: " +
                          self.family, "EllipsoidID: " + self.ellipsoid.name])


def projWithCommonEll(ellipsoid):
    """"""
    with open("projections.txt") as fObj:
        readCSV = reader(fObj, delimiter = ',')
        for line in readCSV:
            if line[2].strip() == ellipsoid.name:
                return Projection(line[0].strip())
        raise Exception("The ellipsoid you defined does not exist in the system")


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    