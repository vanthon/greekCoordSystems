# -*- coding: utf-8 -*-
"""
@author: Antonios Vatalis
"""
from csv import reader
from angles import DMStoRads, Angle
from collections import defaultdict


def loadEllipsoidParams(ellName):
    """
    This function loads refEllipsoids.txt file, which includes different
    ellipsoids with their geometric characteristics. A sample row of this
    file has the following format:
    WGS84, 6378137.000, 6356752.3142, 0.0033528106, 0.0066943799, 0.0067394964
    where each comma delimited value corresponds to: ellipsID, semiMajorAxis,
    semiMinorAxis, flattening, eccentricitySquared, eccentricityPrimeSquared
    geometric quantities respectively.
    Arguments
        ellName: string that specifies ellipsoid's ID,
                 could be one of the followings:
                    'WGS84', 'GRS80', 'BESSEL', 'HAYFORD'

                 * If user desires to add another ellipsoid has to modify
                   refEllipsoids.txt file following its header's instructions
    Returns
        dictionary of ellipsoid's parameters, i.e.
        {ellipsID: string, semiMajorAxis: number, semiMinorAxis: number,
        flattening: number, eccentricitySquared: number,
        eccentricityPrimeSquared: number}
    """
    with open("refEllipsoids.txt") as fObj:
        readCSV = reader(fObj, delimiter=',')
        header = [value.strip() for value in next(readCSV)]
        values = []
        for line in readCSV:
            if line[0].strip() == ellName:
                for value in line:
                    values.append(value.strip())
                values[1:] = [float(x) for x in values[1:]]
                return dict(zip(header, values))
        raise Exception("The ellipsoid you defined does not exist "
                        "in refEllipsoids.txt file")


def loadProjParams(projName):
    """
    This function loads projections.txt file, which includes different
    projections with their correspoding parameters. A sample raw of this
    file has the following format:
    TM3-EAST, TM, BESSEL, 0.9999, 34, 26:42:58.815, 200000, 0
    where each comma delimited value corresponds to:
    projID, projFamily, ellipsoid, scale, lat0, lon0, falseEasting,
    falseNorthing respectively.
    It should be stated here that lat0, lon0 are defined in degrees
    (dms format, e.g. -0:23:9.01239573) within file.
    Arguments
        projName: string that specifies projection's ID,
                  could be one of the followings:
                    'TM3-WEST', 'TM3-CENTRAL', 'TM3-EAST', 'TM87',
                     'UTM-34N', 'UTM-35N', 'HATT', 'TM87-WGS84'

                  * If user desires to add another projection has to modify
                    projections.txt file following its header's instructions
    Returns
        dictionary of projection's parameters, i.e.
        {'projID': string, 'projFamily': string, 'ellipsoid': string,
        'scale': number, 'lat0': number, 'lon0': number,
        'falseEasting': number, 'falseNorthing': number}

        * values of dictionary could be None if they are missing from file
        * CAUTION. lat0, lon0 inserted as degrees in dms format but returned
          in radians
    """  
    with open("projections.txt") as fObj:
        readCSV = reader(fObj, delimiter=',')
        header = [value.strip() for value in next(readCSV)]
        projParams = dict.fromkeys(header)
        for line in readCSV:
            if line[0].strip() == projName:
                for key, value in zip(header, line):
                    projParams[key] = value.strip()
                for key in ['scale', 'falseEasting', 'falseNorthing']:
                    if projParams[key]:
                        projParams[key] = float(projParams[key])
                for key in ['lat0', 'lon0']:
                    if projParams[key]:
                        # projParams[key] = DMStoRads(projParams[key])
                        projParams[key] = Angle(projParams[key])
                return projParams
        raise Exception("The projection you defined does not exist "
                        "in projections.txt file")


def loadCoordsFile(fileName, clmsSpecification, dlm=','):
    """
    This function loads a file that contains points along with their
    respective coordinates. A sample row of this file may have the
    following format:
    p4,4478329.405,1921444.508,4101378.606, 41.464, 42.651
    where each delimited value corresponds to:
    id, E/X/lat, N/Y/lon, U/Z/h, geoidHeight, GMHeight respectively.
    and
    id: point id
    ENU: map coordinates
    XYZ: ECEF coordinates
    latlonh: geodetic coordinates
    geoidHeight: geoid undulation, i.e. geoidHeight = h - U
    GMHeight: geoid undulation coming from a geopotential model as EGM08

    * All coordinates values expressed in meters, except lat, lon where their
      values could be radians or degrees in the following formats:
      decimal degrees-dd, e.g. -0.1927322, degrees:minutes-dm, e.g. 11:2.3474,
      degrees:minutes:seconds-dms, e.g. -0:59:59.0001123

    Arguments
        fileName: string that specifies the name of file to be loaded
        clmsSpecification: list of strings, that specifies what each column
                           represent and in what order the columns are layout.
                           possible list any subset of:
                           ['id', 'X', 'Y', 'Z', lat-r/dd/dm/dms,
                           lon-r/dd/dm/dms, h, E, N, U, 'geoidHeight',
                           'GMHeight']
                           e.g. ['id', 'X', 'Y', 'Z'] or
                                ['id', 'Z', 'Y', 'X', 'geoidHeight']

                            lat-r/dd/dm/dms means that it could be
                            lat-r or lat-dd or lat-dm or lat-dms. The same
                            holds for lon-r/dd/dm/dms
        dlm: one character string that specifies file's delimiter
             e.g. ' ', ',', '\t'
    Returns
        any possible subset (concerning keys) of the follwoing dictionary:
        {'id': [strings], 'X': [numbers], 'Y': [numbers], 'Z': [numbers],
        'lat-r/dd/dm/dms': [numbers/strings],
        'lon-r/dd/dm/dms': [numbers/strings], 'h': [numbers], 'E': [numbers],
        'N': [numbers], 'U': [numbers]
        'geoidHeight': [numbers], 'GMHeight': [numbers]}
    """
    if len(set(clmsSpecification)) != len(clmsSpecification):
        raise Exception("Your columns specification contain duplicates")
    allClms = {'id', 'X', 'Y', 'Z', 'lat-r', 'lon-r', 'lat-dd', 'lon-dd',
               'lat-dm', 'lon-dm', 'lat-dms', 'lon-dms', 'h', 'E', 'N', 'U',
               'geoidHeight', 'GMHeight'}
    if not set(clmsSpecification).issubset(allClms):
        raise Exception("Define columns specification properly")

    with open(fileName) as fObj:
        readCSV = reader(fObj, delimiter=dlm)
        data = defaultdict(list)
        for line in readCSV:
            if line:
                for key, value in zip(clmsSpecification, line):
                    data[key].append(value.strip())
        latlonInDeg = ['lat-dd', 'lon-dd', 'lat-dm', 'lon-dm',
                       'lat-dms', 'lon-dms']
        for key in set(data.keys()) - {'id'}:
            # if key not in latlonInDeg:
            #     data[key] = [float(x) for x in data[key]]
            if key in ['lat-r', 'lon-r']:
                data[key] = [Angle(float(x)) for x in data[key]]
            elif key in latlonInDeg:
                data[key] = [Angle(x) for x in data[key]]
            else:
                data[key] = [float(x) for x in data[key]]
        return data

