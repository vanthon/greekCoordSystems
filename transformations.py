# -*- coding: utf-8 -*-
import numpy as np
from points import MapPoint, GeoPoint, project, unproject, geodeticToECEF
from modelers import projWithCommonEll
from copy import copy, deepcopy
from math import cos, sin


class Transformation:
    def __init__(self, transName, source, target = None,
                 scale=True, rotation=True, shift=True):
        """
        source: list of MapPoints or GeoPoints
        target: list of MapPoints or GeoPoints
        system: type Ellipsoid or Projection
        """
        self.name = transName
        self.system = None
        self.source = source
        self.target = target
        self.scale = scale
        self.rotation = rotation
        self.shift = shift
        self.params = None

        transType = {'similarity': self.similarity, 'affine': self.affine, 'plane': self.plane,
                     'poly2d': self.poly2d, 'spheare': self.spheare}
        transType[transName]()

    def similarity(self):
        """
        returns: list o transformation parameters
                 [c, d, tx, ty] if source & target list of MapPoints
        """
        if all(isinstance(x, MapPoint) for x in (self.source + self.target)):
            self.system = self.target[0].projection

            designArray = []  # array_like
            for point in self.source:
                scaleRotArray = shiftArray = [[], []]
                if self.rotation:
                    scaleRotArray = [[point.x, point.y], [point.y, -point.x]]
                else:
                    if self.scale:
                        scaleRotArray = [[point.x], [point.y]]
                if self.shift:
                    shiftArray = [[1, 0], [0, 1]]
                for x, y in zip(scaleRotArray, shiftArray):
                    designArray.append(x+y)

            target = [coord for point in self.target for coord in [point.x, point.y]]
            paramsValues = np.linalg.lstsq(designArray, target, rcond=None)
            paramsValues = paramsValues[0]

            paramsKeys = ['c', 'd', 'tx', 'ty']
            if not self.rotation and self.scale:
                paramsKeys.pop('c')
                paramsKeys.pop('d')
            if not self.rotation:
                paramsKeys.pop('d')
            if not self.shift:
                paramsKeys.pop('tx')
                paramsKeys.pop('ty')

            params = {'c': 1, 'd': 0, 'tx': 0, 'ty': 0}
            params.update(zip(paramsKeys, paramsValues))

            self.params = params

        elif all(isinstance(x, GeoPoint) for x in (self.source + self.target)):
            self.system = self.target[0].ellipsoid

            source = geodeticToECEF(self.source)
            target = geodeticToECEF(self.target)
            designArray = []  # array_like
            for point in source:
                scaleArray = rotArray = shiftArray = [[], [], []]
                if self.scale:
                    scaleArray = [[point.x], [point.y], [point.z]]
                if self.rotation:
                    rotArray = [[0, -point.z, point.y], [point.z, 0, -point.x],
                                [-point.y, point.x, 0]]
                if self.shift:
                    shiftArray = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

                for x, y, z in zip(scaleArray, rotArray, shiftArray):
                    designArray.append(x+y+z)

            target = [coord for a, b in zip(target, source)
                      for coord in [a.x - b.x, a.y - b.y, a.z - b.z]]

            paramsValues = np.linalg.lstsq(designArray, target, rcond=None)
            paramsValues = paramsValues[0]

            paramsKeys = [x if y else [] for x, y in
                          zip([['dm'], ['ex', 'ey', 'ez'], ['tx', 'ty', 'tz']],
                              [self.scale, self.rotation, self.shift])]
            paramsKeys = [x for sublist in paramsKeys for x in sublist]

            params = {'dm': 0, 'ex': 0, 'ey': 0, 'ez': 0, 'tx': 0, 'ty': 0, 'tz': 0}
            params.update(zip(paramsKeys, paramsValues))

            self.params = params
        else:
            raise Exception('Points should be of MapPoint or GeoPoint class')

    def plane(self):
        """"""
        if isinstance(self.source[0], GeoPoint):
            projection = projWithCommonEll(self.source[0].ellipsoid)
        else:
            projection = self.source[0].projection
        self.system = projection
        source = project(self.source, projection)

        geoids = [point.geoidHeight for point in source]
        designArray = [[1, point.x, point.y] for point in source]
        paramsValues = np.linalg.lstsq(designArray, geoids, rcond=None)
        paramsValues = paramsValues[0]
        paramsKeys = ['a0', 'a1', 'a2']
        self.params = dict(zip(paramsKeys, paramsValues))

    def poly2d(self):
        """"""
        if isinstance(self.source[0], GeoPoint):
            projection = projWithCommonEll(self.source[0].ellipsoid)
        else:
            projection = self.source[0].projection
        self.system = projection
        source = project(self.source, projection)

        geoids = [point.geoidHeight for point in source]
        designArray = [[1, point.x, point.y, point.x*point.y, point.x**2, point.y**2] for point in source]
        paramsValues = np.linalg.lstsq(designArray, geoids, rcond=None)
        paramsValues = paramsValues[0]
        paramsKeys = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5']
        self.params = dict(zip(paramsKeys, paramsValues))

    def spheare(self):
        """"""
        source = unproject(self.source, 'latlonh-r')
        geoids = [point.geoidHeight for point in source]
        designArray = [[1, cos(point.x)*cos(point.y), cos(point.x)*sin(point.y),
                        sin(point.x)] for point in source]
        paramsValues = np.linalg.lstsq(designArray, geoids, rcond=None)
        paramsValues = paramsValues[0]
        paramsKeys = ['a0', 'a1', 'a2', 'a3']
        self.params = dict(zip(paramsKeys, paramsValues))

    def affine(self):
        """"""


def _transform(point, transformation):
    params = transformation.params
    if transformation.name == 'similarity':
        if isinstance(point, MapPoint):
            E = params['c']*point.x + params['d']*point.y + params['tx']
            N = -params['d']*point.x + params['c']*point.y + params['ty']
            return MapPoint(point.id, E, N, point.z, transformation.system, None)
        elif isinstance(point, GeoPoint):
            coordType = point.coordType
            _point = geodeticToECEF(point)
            X = params['dm']*_point.x + params['ez']*_point.y -\
                params['ey']*_point.z + params['tx'] + _point.x
            Y = -params['ez']*_point.x + params['dm']*_point.y +\
                params['ex']*_point.z + params['ty'] + _point.y
            Z = params['ey']*_point.x - params['ex']*_point.y +\
                params['dm']*_point.z + params['tz'] + _point.z
            _point = GeoPoint(point.id, X, Y, Z, 'ECEF', transformation.system, None)
            if coordType != 'ECEF':
                _point.latlonh(point.coordType.split('-')[1])
            return _point
    elif transformation.name == 'plane':
        _point = project(point, transformation.system)
        geoidHeight = params['a0'] + params['a1']*_point.x + params['a2']*_point.y
        _point = deepcopy(point)
        _point.geoidHeight = geoidHeight
        return _point
    elif transformation.name == 'poly2d':
        _point = project(point, transformation.system)
        geoidHeight = (params['a0'] + params['a1']*_point.x + params['a2']*_point.y +
                       params['a3']*_point.x*_point.y + params['a4']*_point.x**2 +
                       params['a5']*_point.y**2)
        _point = deepcopy(point)
        _point.geoidHeight = geoidHeight
        return _point
    elif transformation.name == 'spheare':
        _point = unproject(point, 'latlonh-r')
        geoidHeight = (params['a0'] + params['a1']*cos(_point.x)*cos(_point.y) +
                       params['a2']*cos(_point.x)*sin(_point.y) + params['a3']*sin(_point.x))
        _point = deepcopy(point)
        _point.geoidHeight = geoidHeight
        return _point
    else:
        _point = deepcopy(point)
        return point


def transform(points, transformation):
    if isinstance(points, list):
        return [_transform(point, transformation) for point in points]
    else:
        return _transform(points, transformation)


def transformPointSets(transName, source, target, sourceRest = None, ifMap2Dor3D = '2D'):
    """
    """
    # Make them all GeoPoints  
    _source = unproject(source, 'ECEF')
    _target = unproject(target, 'ECEF')     
    _sourceRest = unproject(sourceRest, 'ECEF')
    #
    if isinstance(target[0], GeoPoint):
        transformation = Transformation(transName, _source, _target)
        _source = transform(_source, transformation)
        _sourceRest = transform(_sourceRest, transformation)
    elif isinstance(target[0], MapPoint):
        if ifMap2Dor3D == '2D':
            transformation = Transformation(transName, _source, _target,
                                            False, False, True)
            _source = transform(_source, transformation)
            _sourceRest = transform(_sourceRest, transformation)
            
            _source = project(_source, target[0].projection)
            _sourceRest = project(_sourceRest, target[0].projection)

            transformation = Transformation(transName, _source, target)
            _source = transform(_source, transformation)
            _sourceRest = transform(_sourceRest, transformation)
        elif ifMap2Dor3D == '3D':
            transformation = Transformation(transName, _source, _target)
            _source = transform(_source, transformation)
            _sourceRest = transform(_sourceRest, transformation)
            
            _source = project(_source, target[0].projection)
            _sourceRest = project(_sourceRest, target[0].projection)
    
    return _source, _sourceRest
            
            
        
    
    
    
    




































