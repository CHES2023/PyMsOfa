# -*- coding: utf-8 -*-

"""
Created on Sat Aug  9 04:15:10 2025

@author: PMO
Description: SOFA constants (the SOFA library header file sofam.h)
SOFA release: 2023-10-11
"""
import numpy as np

# Basic mathematical constants
DPI   = 3.141592653589793238462643           # π
D2PI  = 6.283185307179586476925287           # 2π
DR2D  = 57.29577951308232087679815           # Radians -> degrees
DD2R  = 1.745329251994329576923691e-2        # Degrees -> radians
DR2AS = 206264.8062470963551564734           # Radians -> arcseconds
DAS2R = 4.848136811095359935899141e-6        # Arcseconds -> radians
DS2R  = 7.272205216643039903848712e-5        # Time seconds -> radians
TURNAS = 1296000.0                           # Arcseconds in one full circle
DMAS2R = DAS2R / 1e3                         # Milliarcseconds -> radians


# Time and astronomical calendar constants
DTY   = 365.242198781                        # Length of B1900 tropical year (days)
DAYSEC = 86400.0                             # Seconds per day
DJY   = 365.25                               # Days per Julian year
DJC   = 36525.0                              # Days per Julian century
DJM   = 365250.0                             # Days per Julian millennium
DJ00  = 2451545.0                            # Julian Date of J2000.0 epoch
DJM0  = 2400000.5                            # MJD zero point (JD)
DJM00 = 51544.5                              # MJD of J2000.0 epoch
DJM77 = 43144.0                              # MJD of 1977-01-01.0
TTMTAI = 32.184                              # TT - TAI (seconds)


# Astronomical physical constants
DAU   = 149597870.7e3                        # Astronomical Unit (m, IAU 2012)
CMPS  = 299792458.0                          # Speed of light (m/s)
AULT  = DAU / CMPS                           # Light travel time for 1 au (seconds)
DC    = DAYSEC / AULT                        # Speed of light (au/day)
ELG   = 6.969290134e-10                      # L_G
ELB   = 1.550519768e-8                       # L_B
TDB0  = -6.55e-5                             # TDB offset (seconds)
SRS   = 1.97412574336e-8                     # Solar Schwarzschild radius (au)
 

# -------------------------------
# Auxiliary math functions using numpy
# -------------------------------

def dint(a):
    """Truncate toward zero."""
    return np.ceil(a) if a < 0.0 else np.floor(a)

def dnint(a):
    """Round to nearest integer."""
    if np.abs(a) < 0.5:
        return 0.0
    return np.ceil(a - 0.5) if a < 0.0 else np.floor(a + 0.5)

def dsign(a, b):
    """Return |a| with the sign of b."""
    return -np.abs(a) if b < 0.0 else np.abs(a)

def gmax(a, b):
    """Return maximum of a and b."""
    return a if a > b else b

def gmin(a, b):
    """Return minimum of a and b."""
    return a if a < b else b


# Reference ellipsoid    
WGS84 = 1
GRS80 = 2
WGS72 = 3
