# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 18:16:44 2019

@author: anthon
"""
from math import sqrt, cos, sin, tan

#def directHATT(lat, lon, ellipsoid, lat0, lon0):
def directHATT(lat, lon, projection):
    """
    lat, lon, lat0, lon0 in radians
    ellipsoid of class Ellipsoid
    returns x, y in meters
    """
    ellipsoid = projection.ellipsoid
    lat0 = projection.lat0
    lon0 = projection.lon0
    
    hta2o = ellipsoid.ePrimeSquared * ((cos(lat0))**2)
    to = tan(lat0)
    Wo = sqrt(1 - ellipsoid.eSquared * ((sin(lat0))**2))
    No = ellipsoid.a / Wo
    Mo = (ellipsoid.a * (1 - ellipsoid.eSquared)) / (Wo**3)

    dlat = lat-lat0
    dlon = lon-lon0

    x = No*(cos(lat0))*dlon -\
        Mo*(sin(lat0))*dlon*dlat -\
        (Mo*(cos(lat0))*(2+9*hta2o*(to**2))*dlon*(dlat**2))/6 -\
        (No*(cos(lat0))*((sin(lat0))**2)*(dlon**3))/6 -\
        (No*(sin(lat0))*(1-2*((sin(lat0))**2))*(dlon**3)*dlat)/6

    y = (Mo*dlat) +\
        ((No*(cos(lat0))*(sin(lat0))*(dlon**2))/2) +\
        ((3*(Mo**2)*ellipsoid.ePrimeSquared*(cos(lat0))*(sin(lat0)) *
         (dlat**2))/(2*No)) +\
        ((Mo*(1-4*((sin(lat0))**2)+hta2o*((cos(lat0))**2))*dlat *
         (dlon**2))/6) +\
        (((Mo*ellipsoid.ePrimeSquared*(1-2*((sin(lat0))**2))*(to**2)) *
         (dlat**3))/2) +\
        ((No*(sin(lat0))*(cos(lat0))*(1-2*((sin(lat0))**2)) *
         (dlon**4))/24) -\
        ((No*(sin(lat0))*(cos(lat0))*(dlat**2)*(dlon**2))/3)

    return x, y

#def invHATT(x, y, ellipsoid, lat0, lon0):
def invHATT(x, y, projection):
    """
    x, y in meters
    lat0, lon0 in radians
    returns lat, lon in radians
    """
    ellipsoid = projection.ellipsoid
    lat0 = projection.lat0
    lon0 = projection.lon0
    
    hta2o = ellipsoid.ePrimeSquared * ((cos(lat0))**2)
    to = tan(lat0)
    Wo = sqrt(1 - ellipsoid.eSquared * ((sin(lat0))**2))
    No = ellipsoid.a / Wo
    Mo = (ellipsoid.a*(1 - ellipsoid.eSquared)) / (Wo**3)

    dlat = (y/Mo) -\
           (to*(x**2))/(2*Mo*No) -\
           (3*ellipsoid.ePrimeSquared*(sin(2*lat0))*(y**2)) / (4*Mo*No) -\
           ((1+3*(to**2)-ellipsoid.ePrimeSquared*(1-10*((sin(lat0))**2))) *
            (x**2)*y) / (6*Mo*(No**2)) -\
           (ellipsoid.ePrimeSquared*(1-2*((sin(lat0)**2)))*(y**3)) /\
           (2*(Mo**2)*No) +\
           (to*(1+3*(to**2))*(x**4))/(24*(Mo**2)*(No**2)) -\
           (to*(2+3*(to**2))*(x**2)*(y**2))/(6*(Mo**2)*(No**2))

    dlon = (x/(No*(cos(lat0)))) +\
           (to*x*y)/((No**2)*(cos(lat0))) +\
           ((1+3*(to**2)+hta2o)*x*(y**2))/(3*(No**3)*(cos(lat0))) -\
           ((to**2)*(x**3))/(3*(No**3)*(cos(lat0))) +\
           (to*(2+3*(to**2))*x*(y**3))/(3*(No**4)*(cos(lat0))) -\
           (to*(1+3*(to**2))*(x**3)*y)/(3*(No**4)*(cos(lat0)))

    lat = dlat + lat0  # in radians
    lon = dlon + lon0  # in radians
    return lat, lon
