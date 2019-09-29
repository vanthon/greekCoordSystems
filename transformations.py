# -*- coding: utf-8 -*-
"""
@author: Antonios Vatalis
"""
import numpy as np
from points import MapPoint, GeoPoint, ECEFpoint, Points
from modelers import filterProj
from copy import deepcopy
from math import cos, sin


def sim2D(source: Points, target: Points, scale=True, rotation=True, shift=True):
    """

    :param source: list of lists, [[E1, N1], [E2, N2], ...]
    :param target: list of lists, [[E'1, N'1], [E'2, N'2], ...]
    :param shift:
    :param rotation:
    :param scale:
    :return: dictionary of parameters, c, d, tx, ty
    """
    if not (source.all(MapPoint) and target.all(MapPoint)):
        raise Exception('Points should be of MapPoint type')

    A = []  # array_like
    for point in source:
        E, N = point.E, point.N
        scaleRotArray = shiftArray = [[], []]
        if rotation:
            scaleRotArray = [[E, N], [N, -E]]
        else:
            if scale:
                scaleRotArray = [[E], [N]]
        if shift:
            shiftArray = [[1, 0], [0, 1]]
        for x, y in zip(scaleRotArray, shiftArray):
            A.append(x + y)

    b = [coord for point in target for coord in [point.E, point.N]]
    paramsValues = np.linalg.lstsq(A, b, rcond=None)
    paramsValues = paramsValues[0]

    paramsKeys = ['c', 'd', 'tx', 'ty']
    if (not rotation) and (not scale):
        paramsKeys.pop(0)
        paramsKeys.pop(1)
    if not rotation:
        paramsKeys.pop(1)
    if not shift:
        paramsKeys.pop(2)
        paramsKeys.pop(3)

    params = {'c': 1, 'd': 0, 'tx': 0, 'ty': 0}
    params.update(zip(paramsKeys, paramsValues))
    params['coordSystem'] = target.coordSystem

    return params


def applySim2D(point: MapPoint, params: dict):
    """

    :param point: list of E, N coordinates, i.e [E, N]
    :param params: dictionary of paramters c, d, tx, ty
    :return: E', N' coordinates
    """
    if not isinstance(point, MapPoint):
        raise Exception('point should be of MapPoint')
    E, N = point.E, point.N
    Epr = params['c'] * E + params['d'] * N + params['tx']
    Npr = -params['d'] * E + params['c'] * N + params['ty']
    return MapPoint(point.id, Epr, Npr, point.U, params['coordSystem'],
                    point.geoidHeight, point.GMHeight)


def sim3D(source: Points, target: Points, scale=True, rotation=True, shift=True):
    """

    :param source:
    :param target:
    :param scale:
    :param rotation:
    :param shift:
    :return:
    """
    if not (source.all(ECEFpoint) and target.all(ECEFpoint)):
        raise Exception('Points should be of ECEFpoint type')

    A = []  # array_like
    for point in source:
        X, Y, Z = point.X, point.Y, point.Z
        scaleArray = rotArray = shiftArray = [[], [], []]
        if scale:
            scaleArray = [[X], [Y], [Z]]
        if rotation:
            rotArray = [[0, -Z, Y], [Z, 0, -X], [-Y, X, 0]]
        if shift:
            shiftArray = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        for x, y, z in zip(scaleArray, rotArray, shiftArray):
            A.append(x+y+z)

    b = [x - y for a, b in zip(target, source) for x, y in
         zip([a.X, a.Y, a.Z], [b.X, b.Y, b.Z])]

    paramsValues = np.linalg.lstsq(A, b, rcond=None)
    paramsValues = paramsValues[0]

    paramsKeys = [x if y else [] for x, y in
                  zip([['dm'], ['ex', 'ey', 'ez'], ['tx', 'ty', 'tz']],
                      [scale, rotation, shift])]
    paramsKeys = [x for sublist in paramsKeys for x in sublist]

    params = {'dm': 0, 'ex': 0, 'ey': 0, 'ez': 0, 'tx': 0, 'ty': 0, 'tz': 0}
    params.update(zip(paramsKeys, paramsValues))
    params['coordSystem'] = target.coordSystem

    return params


def applySim3D(point: ECEFpoint, params: dict):
    if not isinstance(point, ECEFpoint):
        raise Exception('point should be of type ECEFpoint')

    X, Y, Z = point.X, point.Y, point.Z
    Xpr = (params['dm'] * X + params['ez'] * Y - params['ey'] * Z +
           params['tx'] + X)
    Ypr = (-params['ez'] * X + params['dm'] * Y + params['ex'] * Z +
           params['ty'] + Y)
    Zpr = (params['ey'] * X - params['ex'] * Y + params['dm'] * Z +
           params['tz'] + Z)
    return ECEFpoint(point.id, Xpr, Ypr, Zpr, params['coordSystem'],
                     point.geoidHeight, point.GMHeight)


def plane(source: Points, GM=False):
    if not source.all(MapPoint):
        raise Exception('points should be of type MapPoint')

    if GM:
        geoids = [point.geoidHeight - point.GMHeight for point in source]
    else:
        geoids = [point.geoidHeight for point in source]
    designArray = [[1, point.E, point.N] for point in source]
    paramsValues = np.linalg.lstsq(designArray, geoids, rcond=None)
    paramsValues = paramsValues[0]
    paramsKeys = ['a0', 'a1', 'a2']
    params = dict(zip(paramsKeys, paramsValues))
    params['GM'] = GM
    return params


def applyPlane(point, params):
    if not isinstance(point, MapPoint):
        raise Exception('point should be of MapPoint type')

    if params['GM']:
        if point.GMHeight:
            GMHeight = point.GMHeight
        else:
            raise Exception('You apply a transformation that includes GMHeight'
                            'while your point does not have a valid GMHeight')
    else:
        GMHeight = 0

    geoidHeight = params['a0'] + params['a1']*point.E + params['a2']*point.N + GMHeight
    mapPoint = deepcopy(point)
    mapPoint.geoidHeight = geoidHeight
    return mapPoint


def poly2d(source: Points, GM=False):
    if not source.all(MapPoint):
        raise Exception('points should be of type MapPoint')

    if GM:
        geoids = [point.geoidHeight - point.GMHeight for point in source]
    else:
        geoids = [point.geoidHeight for point in source]

    designArray = [[1, point.E, point.N, point.E * point.N, point.E ** 2, point.N ** 2] for point in source]
    paramsValues = np.linalg.lstsq(designArray, geoids, rcond=None)
    paramsValues = paramsValues[0]
    paramsKeys = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5']
    params = dict(zip(paramsKeys, paramsValues))
    params['GM'] = GM
    return params


def applyPoly2d(point: MapPoint, params):
    if not isinstance(point, MapPoint):
        raise Exception('point should be of MapPoint type')

    if params['GM']:
        if point.GMHeight:
            GMHeight = point.GMHeight
        else:
            raise Exception('You apply a transformation that includes GMHeight'
                            'while your point does not have a valid GMHeight')
    else:
        GMHeight = 0

    geoidHeight = (params['a0'] + params['a1'] * point.E + params['a2'] * point.N +
                   params['a3'] * point.E * point.N + params['a4'] * point.E ** 2 +
                   params['a5'] * point.N ** 2) + GMHeight
    mapPoint = deepcopy(point)
    mapPoint.geoidHeight = geoidHeight
    return mapPoint


def sphere(source: Points, GM=False):
    if not source.all(GeoPoint):
        raise Exception('points should be of type GeoPoint')

    if GM:
        geoids = [point.geoidHeight - point.GMHeight for point in source]
    else:
        geoids = [point.geoidHeight for point in source]

    designArray = [[1, cos(point.lat.radians) * cos(point.lon.radians), cos(point.lat.radians) * sin(point.lon.radians),
                    sin(point.lat.radians)] for point in source]
    paramsValues = np.linalg.lstsq(designArray, geoids, rcond=None)
    paramsValues = paramsValues[0]
    paramsKeys = ['a0', 'a1', 'a2', 'a3']
    params = dict(zip(paramsKeys, paramsValues))
    params['GM'] = GM
    return params


def applySphere(point: GeoPoint, params):
    if not isinstance(point, GeoPoint):
        raise Exception('point should be of MapPoint type')

    if params['GM']:
        if point.GMHeight:
            GMHeight = point.GMHeight
        else:
            raise Exception('You apply a transformation that includes GMHeight'
                            'while your point does not have a valid GMHeight')
    else:
        GMHeight = 0

    geoidHeight = (params['a0'] +
                   params['a1']*cos(point.lat.radians)*cos(point.lon.radians) +
                   params['a2']*cos(point.lat.radians)*sin(point.lon.radians) +
                   params['a3']*sin(point.lat.radians)) + GMHeight
    geoPoint = deepcopy(point)
    geoPoint.geoidHeight = geoidHeight
    return geoPoint


class Transformation:
    def __init__(self, transName: str, source: Points, target=None,
                 scale=True, rotation=True, shift=True, GM=False):
        """
        source: class of Points
        target: class of Points
        system: type Ellipsoid or Projection
        """
        self.source = source
        self.target = target
        self.coordSystem = None
        self.scale = scale
        self.rotation = rotation
        self.shift = shift
        self.GM = GM
        self.params = None

        self.transName = transName

        # if self.GM:
        #     GMheights = egm08(self.source)  # list of geoid undulations of source point
        #     for point, height in zip(self.source, GMheights):
        #         point.GMHeight = height

        self.transType = {'sim2D': self.applySim2D, 'sim3D': self.applySim3D,
                          'plane': self.applyPlane, 'poly2d': self.applyPoly2D,
                          'sphere': self.applySphere}

    @property
    def sim2D(self):
        """
        returns: list o transformation parameters
                 [c, d, tx, ty] if source & target list of MapPoints
        """
        self.coordSystem = self.target.coordSystem
        return sim2D(self.source, self.target, self.scale, self.rotation,
                     self.shift)

    def applySim2D(self, point):
        return applySim2D(point, self.sim2D)

    @property
    def sim3D(self):
        self.coordSystem = self.target.coordSystem
        return sim3D(self.source, self.target, self.scale, self.rotation,
                     self.shift)

    def applySim3D(self, point):
        return applySim3D(point, self.sim3D)

    @property
    def plane(self):
        """"""
        return plane(self.source, self.GM)

    def applyPlane(self, point):
        return applyPlane(point, self.plane)

    @property
    def poly2d(self):
        """"""
        return poly2d(self.source, self.GM)

    def applyPoly2D(self, point):
        return applyPoly2d(point, self.poly2d)

    @property
    def sphere(self):
        """"""
        return sphere(self.source, self.GM)

    def applySphere(self, point):
        return applySphere(point, self.sphere)

    def affine(self):
        """"""

    def transform(self, point):
        return self.transType[self.transName](point)


def transform(transName, source: Points, target: Points, sourceRest: Points, GM=False):
    """
    source: Points
    target: Points
    sourceRest: list of GeoPoints or MapPoints
    """
    if GM and (transName in ['plane', 'poly2d', 'sphere']):
        if not all(source.GMHeight):
            for point in source:
                point.setGMHeight()
        if not all(sourceRest.GMHeight):
            for point in sourceRest:
                point.setGMHeight()

    if transName == 'sim3D':
        # Make all ECEFpoints
        ECEFsource = source.ECEF()
        ECEFtarget = target.ECEF()
        ECEFsourceRest = sourceRest.ECEF()
        #
        transformation = Transformation('sim3D', ECEFsource, ECEFtarget)
        sourceTR = ECEFsource.transform(transformation)
        sourceRestTR = ECEFsourceRest.transform(transformation)
        if target.elmsType is GeoPoint:
            sourceTR = sourceTR.geodetic()
            sourceRestTR = sourceRestTR.geodetic()
        elif target.elmsType is MapPoint:
            sourceTR = sourceTR.project(target.coordSystem)
            sourceRestTR = sourceRestTR.project(target.coordSystem)
        return sourceTR, sourceRestTR

    elif transName == 'sim2D':
        if target.elmsType is not MapPoint:
            raise Exception('')
        # Make all ECEFpoints
        ECEFsource = source.ECEF()
        ECEFtarget = target.ECEF()
        ECEFsourceRest = sourceRest.ECEF()
        #
        transformation = Transformation('sim3D', ECEFsource, ECEFtarget,
                                        scale=False, rotation=False, shift=True)
        sourceTR = ECEFsource.transform(transformation)
        sourceRestTR = ECEFsourceRest.transform(transformation)
        #
        sourceTR = sourceTR.project(target.coordSystem)
        sourceRestTR = sourceRestTR.project(target.coordSystem)
        transformation = Transformation('sim2D', sourceTR, target)
        sourceTR = sourceTR.transform(transformation)
        sourceRestTR = sourceRestTR.transform(transformation)
        return sourceTR, sourceRestTR

    elif transName in ['plane', 'poly2d']:
        sourceRestTR = deepcopy(sourceRest)
        if source.elmsType is not MapPoint:
            projection = filterProj(source.coordSystem)
            source = source.project(projection)
            sourceRest = sourceRest.project(projection)
        transformation = Transformation(transName, source, GM=GM)
        mapSourceRestTR = sourceRest.transform(transformation)
        sourceRestTR.geoidHeight = mapSourceRestTR.geoidHeight
        return sourceRestTR

    elif transName in 'sphere':
        sourceRestTR = deepcopy(sourceRest)
        if source.elmsType is not GeoPoint:
            source = source.geodetic()
            sourceRest = sourceRest.geodetic()
        transformation = Transformation(transName, source, GM=GM)
        mapSourceRestTR = sourceRest.transform(transformation)
        sourceRestTR.geoidHeight = mapSourceRestTR.geoidHeight
        return sourceRestTR
