# -*- coding: utf-8 -*-
import numpy as np
from points import MapPoint, GeoPoint


def transform(point, transformation):
    params = transformation.params
    if transformation.name == 'similarity':
        if isinstance(point, MapPoint):
            E = params['c']*point.E + params['d']*point.N + params['tx']
            N = -params['d']*point.E + params['c']*point.N + params['ty']
            return MapPoint(point.id, E, N, None, point.projection, None)
        elif isinstance(point, GeoPoint):
            coordType = point.coordType
            point.ECEF()
            X = params['dm']*point.x + params['ez']*point.y -\
                params['ey']*point.z + params['tx'] + point.x
            Y = -params['ez']*point.x + params['dm']*point.y +\
                params['ex']*point.z + params['ty'] + point.y
            Z = params['ey']*point.x - params['ex']*point.y +\
                params['dm']*point.z + params['tz'] + point.z
            gPt = GeoPoint(point.id, X, Y, Z, 'ECEF', point.ellipsoid, None)
            if coordType != 'ECEF':
                gPt.latlonh(point.coordType.split('-')[1])
            return gPt


class Transformation:
    def __init__(self, transName, source, target, scale=True, rotation=True, shift=True):
        """
        source: list of MapPoints or GeoPoints
        target: list of MapPoints or GeoPoints
        """
        self.name = transName
        self.source = source
        self.target = target
        self.scale = scale
        self.rotation = rotation
        self.shift = shift
        self.params = None

        transType = {'similarity': self.similarity, 'affine': self.affine}
        transType[transName]()

    def shift(self, tx, ty, tz=None):
        """Takes both kind of points, returns both"""
        if tz:
            self.params = [tx, ty, tz]
        else:
            self.params = [tx, ty]

    def scale(self, scale):
        """Takes both kind of points, returns both"""
        self.params = [scale]

    def rotate(self, alpha, beta = None, gama = None):
        """1 arg for MapPoints, 3 for GeoPoints"""
        self.params = []

    def similarity(self):
        """
        returns: list o transformation parameters
                 [c, d, tx, ty] if source & target list of MapPoints
        """
        if all(isinstance(x, MapPoint) for x in (self.source + self.target)):
            designArray = []  # array_like
            for point in self.source:
                scaleRotArray = shiftArray = [[], []]
                if self.rotation:
                    scaleRotArray = [[point.E, point.N], [point.N, -point.E]]
                else:
                    if self.scale:
                        scaleRotArray = [[point.E], [point.N]]
                if self.shift:
                    shiftArray = [[1, 0], [0, 1]]
                for x, y in zip(scaleRotArray, shiftArray):
                    designArray.append(x+y)

            target = [coord for point in self.target for coord in [point.E, point.N]]
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
            source = [(point.ECEF(), point)[1] for point in self.source]
            target = [(point.ECEF(), point)[1] for point in self.target]
            designArray = []  # array_like
            for point in self.source:
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


    def affine(self):
        """"""

    def geoidSurface(self):
        """Takes GeoPoint returns GeoPoint augmented with its geoid height"""
