# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 18:18:09 2019

@author: Antonios Vatalis
"""

from math import degrees, radians


def radToDMS(rads, degFormat, prec = None):
    """
    rads: radians
    degFormat: "dd", "dm", "dms" of string type
    prec: desired precision of decimal part, type integer
    returns d, dm, dms of string type
            e.g. "-34.4784326" or "+34:29.04784" or "-0:29:19.03847"
    """
    if rads < 0:
        dSign = "-"; rads = abs(rads)
    else: dSign = ""
    dms = degrees(rads)
    
    if degFormat == "dms":
        m, s = divmod(dms * 3600, 60)
        if prec is not None:
            s = round(s, prec)
        d, m = divmod(m, 60)
        d = int(d); m = int(m)
        dms = [str(x) for x in [d, m, s]]; dms = ":".join(dms)
    elif degFormat == "dm":
        d, m = divmod(dms * 60, 60)
        d = int(d)
        if prec is not None:
            m = round(m, prec)
        dms = [str(x) for x in [d, m]]; dms = ":".join(dms)
    elif degFormat == "dd":
        if prec is not None:
            dms = round(dms, prec)
        dms = str(dms)
    else:
        raise Exception("degFormat one of the strings: 'dd' or 'dm' or 'dms'")
    return dSign + dms


def DMStoRads(dms):
    """
    dms: dd or dm or dms of string type
         e.g. "-0.4784326" or "+34:29.04784" or "-34:29:19.03847" 
    returns radians
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
    return radians(dSign * dms)
