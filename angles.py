# -*- coding: utf-8 -*-
"""
@author: Antonios Vatalis
"""
from math import degrees, radians, pi
latAthensBessel = '23:42:58.815'  # Athens benchmark latitude in dms (Bessel ellipsoid)


def radToDMS(rads, degFormat, prec=None):
    """
    This functions returns radians to degrees in various formats.
    Arguments
        rads: number in radians
              valid value: 0 <= rads <= 2*pi
                           -2pi <= rads <= 0
        degFormat: string that specifies the desired format of returned degrees
                   "dd" for decimal degrees, e.g. -22.0912233
                   "dm" for degrees:minutes, e.g. 34:23.45902
                   "dms" for degrees:minutes:seconds, e.g. -0:29:59.0387822
        prec: integer that specifies the derided precision on decimal part of returned degrees
    Returns
        string as degrees on the format that argument degFormat specifies
    """
    if not any([0 <= rads <= 2*pi, -2*pi <= rads <= 0]):
        raise ValueError('rads should be within [0, 2pi], [0, -2pi]')

    if rads < 0:
        dSign = "-"
        rads = abs(rads)
    else:
        dSign = ""
    dms = degrees(rads)
    
    if degFormat == "dms":
        m, s = divmod(dms * 3600, 60)
        if prec is not None:
            s = round(s, prec)
        d, m = divmod(m, 60)
        d = int(d)
        m = int(m)
        dms = [str(x) for x in [d, m, s]]
        dms = ":".join(dms)
    elif degFormat == "dm":
        d, m = divmod(dms * 60, 60)
        d = int(d)
        if prec is not None:
            m = round(m, prec)
        dms = [str(x) for x in [d, m]]
        dms = ":".join(dms)
    elif degFormat == "dd":
        if prec is not None:
            dms = round(dms, prec)
        dms = str(dms)
    else:
        raise Exception("The returned angle should have one of the "
                        "following formats: 'dd' or 'dm' or 'dms'")
    return dSign + dms


def DMStoRads(dms, prec=None):
    """
    This function returns degrees of various formats to radians
    Arguments
        dms: string as degrees on the following formats
             decimal degrees of type dd, e.g. -22.0912233
             degrees:minutes of type dm, e.g. 34:23.45902
             degrees:minutes:seconds of type dms, e.g. -0:29:59.0387822
    Returns
        number in radians
    """
    dms = dms.strip()
    
    if dms.startswith('-'):
        dSign = -1
    else:
        dSign = 1
    
    dms = dms.split(':')

    dms = [abs(float(x)) for x in dms]
    if len(dms) == 3:
        dms[2] = dms[2] / 3600
        dms[1] = dms[1] / 60
    elif len(dms) == 2:
        dms[1] = dms[1] / 60
    dms = sum(dms)
    rads = radians(dSign * dms)

    if not any([0 <= rads <= 2 * pi, -2 * pi <= rads <= 0]):
        raise ValueError('degrees should be within [0, 360], [0, -360]')

    if prec:
        return round(rads, prec)
    else:
        return rads


def toRads(theta, prec=None):
    if isinstance(theta, str):
        return DMStoRads(theta, prec)
    else:
        return theta


def toDMS(theta, degFormat, prec=None):
    if isinstance(theta, str):
        theta = toRads(theta)
    return radToDMS(theta, degFormat, prec)


def shiftAngle(angle, dangle, shAngleFormat='dms'):
    """
    This function shifts an angle by a specified offset value
    Arguments
        angle: number or string as angle to be shifted
               number: angle in radians
               string: angle in degrees of the following formats
                        decimal degrees (dd), e.g. -22.0912233
                        degrees:minutes (dm), e.g. 34:23.45902
                        degrees:minutes:seconds (dms), e.g. -0:19:9.0382
        dangle: number or string as offset value with same type as lon
        shAngleFormat: string that specifies the format of returned shifted angle
                         "r" for radians, e.g. 0.0121345
                            where the returned value is a number,
                         "dd" for decimal degrees, e.g. -22.0912233
                            where the returned value is a string,
                         "dm" for degrees:minutes, e.g. 34:23.45902
                            where the returned value is a string,
                         "dms" for degrees:minutes:seconds, e.g. -0:29:59.0387822
                            where the returned value is a string
                       If not specified returns the shifted angle in dms format
    Returns
        number or string as shifted angle of lonFormat format
    """
    if isinstance(angle, str):
        angle = DMStoRads(angle)
    if isinstance(dangle, str):
        dangle = DMStoRads(dangle)
    angle = angle + dangle
    if shAngleFormat in ['dd', 'dm', 'dms']:
        return radToDMS(angle, shAngleFormat)
    elif shAngleFormat == 'r':
        return angle
    else:
        raise Exception('The shifted angle should be on the following formats: '
                        'r, dd, dm, dms')


class Angle:
    def __init__(self, value=None):
        # check if value within appropriate boundaries
        # moreover, check if value is of correct type of string
        if value:
            toRads(value)
            toDMS(value, 'dms')
        # end of check
        self.value = value
        self.prec = 5

    @property
    def angleType(self):
        if self.value:
            degFormat = {0: 'dd', 1: 'dm', 2: 'dms'}
            try:
                return degFormat[self.value.count(':')]
            except AttributeError:
                return 'r'
        else:
            return None

    @angleType.setter
    def angleType(self, newType):
        if self.value:
            if newType in ['dd', 'dm', 'dms']:
                self.value = toDMS(self.value, newType)
            elif newType == 'r':
                self.value = toRads(self.value)
            else:
                raise Exception('new lat, lon type should be r, dd, dm, dms')
        else:
            raise Exception('You cannot set new angle type for Angle class '
                            'that has initialized with None value')

    @property
    def radians(self):
        if self.value:
            return toRads(self.value)
        else:
            return None

    @property
    def dd(self):
        if self.radians:
            return toDMS(self.radians, 'dd')
        else:
            return None

    @property
    def dm(self):
        if self.radians:
            return toDMS(self.radians, 'dm')
        else:
            return None

    @property
    def dms(self):
        if self.radians:
            return toDMS(self.radians, 'dms')
        else:
            return None

    def shift(self, other=None):
        try:
            if other is not None:
                theta = self.radians + other.radians
            else:
                theta = self.radians + Angle(latAthensBessel).radians
        except TypeError:
            return None
        theta = Angle(theta)
        theta.angleType = 'dms'
        return theta

    def unshift(self, other=None):
        try:
            if other is not None:
                theta = self.radians - other.radians
            else:
                theta = self.radians - Angle(latAthensBessel).radians
        except TypeError:
            return None
        theta = Angle(theta)
        theta.angleType = 'dms'
        return theta

    def __bool__(self):
        if self.value is not None:
            return True
        else:
            return False

    def __str__(self):
        if self.angleType in ['dd', 'dm', 'dms']:
            return toDMS(self.value, self.angleType, self.prec)
        elif self.angleType == 'r':
            return str(round(self.value, self.prec))
        else:
            return str(None)

    def __repr__(self):
        return str(self.value)
