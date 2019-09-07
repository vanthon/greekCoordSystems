# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 11:23:43 2019

@author: anthon
"""
from csv import reader


def loadEllipsoidParams(ellName):
    """
    ellName: string, 'WGS84', 'GRS80', 'BESSEL', 'HAYFORD'
    """
    with open("refEllipsoids.txt") as fObj:
        readCSV = reader(fObj, delimiter = ',')
        header = [elm.strip() for elm in next(readCSV)]
        values = len(header) * [None]
        for line in readCSV:
            if line[0].strip() == ellName:
                for indx, elm in enumerate(line):
                    values[indx] = elm.strip()
                return dict(zip(header, values))
        raise Exception("The ellipsoid you defined does not exist " +\
                        "in the system")


def loadProjParams(projName):
    """
    projName string: "TM3-WEST", "TM3-CENTRAL", "TM3-EAST", "TM87", 
                     "UTM-34N", "UTM-35N"
    """  
    with open("projections.txt") as fObj:
        readCSV = reader(fObj, delimiter = ',')
        header = [elm.strip() for elm in next(readCSV)]
        values = len(header) * [None]
        for line in readCSV:
            if line[0].strip() == projName:
                for indx, elm in enumerate(line):
                    values[indx] = elm.strip()
                return dict(zip(header, values))
        raise Exception("The projection you defined does not exist " +\
                    "in the system")


def loadCoordsFile(fileName):
    """"""
    with open(fileName) as fObj:
        readCSV = reader(fObj, delimiter = ',')
        data = []
        for line in readCSV:
            if line:
                PtIDxyzGeoid = 5 * [None]  # [PtID, x, y, z, GeoidHeight]
                for count, lineEl in enumerate(line):
                    PtIDxyzGeoid[count] = lineEl.strip()
                data.append(PtIDxyzGeoid)
        return data
