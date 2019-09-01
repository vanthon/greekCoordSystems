# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 20:42:35 2019

@author: Antonios Vatalis
"""
from math import radians
from loadFiles import loadProjParams
from ellipsoids import Ellipsoid
from importlib import import_module


class Projection():
    def __init__(self, projName, lat0 = None, lon0 = None):
        """
        projName string: "HATT", "TM3-WEST", "TM3-CENTRAL", "TM3-EAST", "TM87",
                         "UTM-34N", "UTM-35N"
        lat0, lon0 in radians
        """
        params = loadProjParams(projName)
        self.family = params[1]
        params[2] = Ellipsoid(params[2]) # ellName to Ellipsoid class
        if self.family == 'TM':  
            params[4] = radians(params[4]) # lat0 from degrees to radians
            params[5] = radians(params[5]) # lon0 from degrees to radians  
        if lat0 != None and lon0 != None:
            params = params + [lat0, lon0]
        
        self.name = params[0]
        self.ellipsoid = params[2]
        self.params = params

    def toMap(self, lat, lon):
        """
        lat, lon in radians
        return E, N in meters
        """
        eqnsFamily = import_module(self.family + "eqns")
        directEqns = getattr(eqnsFamily, "direct" + self.family)
        E, N = directEqns(lat, lon, *self.params[2:])
        return E, N

    def fromMap(self, E, N):
        """
        lat, lon in radians
        return lat, lon in radians
        """
        eqnsFamily = import_module(self.family + "eqns")
        invEqns = getattr(eqnsFamily, "inv" + self.family)
        lat, lon = invEqns(E, N, *self.params[2:])
        return lat, lon

    def __repr__(self):
        return ", ".join(["projID: " + self.name, "projFamily: " +
                          self.family, "EllipsoidID: " + self.ellipsoid.name])

