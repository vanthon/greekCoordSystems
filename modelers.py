# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:25:03 2019

@author: anthon
"""
from math import radians
from importlib import import_module
from loadFiles import loadEllipsoidParams, loadProjParams
from angles import DMStoRads


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


class Projection():
    def __init__(self, projName, lat0 = None, lon0 = None):
        """
        projName string: "HATT", "TM3-WEST", "TM3-CENTRAL", "TM3-EAST", "TM87",
                         "UTM-34N", "UTM-35N"
        lat0, lon0: numeral for radians
                    string for dd, dm, dms. e.g. '23.33', '3:34.2', '-0:3:4.33'
        """
        params = loadProjParams(projName)
        self.family = params[1]
        params[2] = Ellipsoid(params[2]) # ellName to Ellipsoid class
        try:
            params[4] = radians(params[4]) # lat0 from degrees to radians
            params[5] = radians(params[5]) # lon0 from degrees to radians
        except IndexError:
            pass
            
        if all([type(lat0) is str, type(lon0) is str]):
            lat0 = DMStoRads(lat0)
            lon0 = DMStoRads(lon0)
            params = params + [lat0, lon0]
        else:
            if all([lat0 is not None, lon0 is not None]):
                params = params + [lat0, lon0]
        
        self.name = params[0]
        self.ellipsoid = params[2]
        self.params = params
        
        eqnsFamily = import_module(self.family + "eqns")
        self.directEqns = getattr(eqnsFamily, "direct" + self.family)
        self.invEqns = getattr(eqnsFamily, "inv" + self.family)

    def __repr__(self):
        return ", ".join(["projID: " + self.name, "projFamily: " +
                          self.family, "EllipsoidID: " + self.ellipsoid.name])
