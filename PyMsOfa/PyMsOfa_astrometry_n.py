"""
Created on Sat Aug  9  04:15:10 2023
Done    on Mon Oct  26 18:35:26 2025
@author: Dr. Jianghui JI  (jijh@pmo.ac.cn)
Description: SOFA Astrometry
"""

import numpy as np
try:                                  # installed as part of the PyMsOfa package
    from .sofa_const import *
    from .PyMsOfa_basic_n import *
    from .PyMsOfa_time_n import (pymEpj, pymEpb2jd,
                                 pymTaitt, pymUtctai, pymUtcut1)
    from .PyMsOfa_earth_attitude_n import (pymPnm06a, pymFw2m, pymBpn2xy, pymS06,
                                           pymC2ixys, pymPom00, pymEors, pymEra00,
                                           pymObl06, pymPmat06, pymSp00)
except ImportError:                   # flat layout: all modules in one folder
    from sofa_const import *
    from PyMsOfa_basic_n import *
    from PyMsOfa_time_n import (pymEpj, pymEpb2jd,
                                pymTaitt, pymUtctai, pymUtcut1)
    from PyMsOfa_earth_attitude_n import (pymPnm06a, pymFw2m, pymBpn2xy, pymS06,
                                          pymC2ixys, pymPom00, pymEors, pymEra00,
                                          pymObl06, pymPmat06, pymSp00)
         
      
#2025-09-28
def pymAb(pnat, v, s, bm1):
    """
    Apply aberration to transform natural direction into proper direction.

    Parameters
    ----------
    pnat : array_like, shape (3,)
        Natural direction to the source (unit vector).
    v : array_like, shape (3,)
        Observer barycentric velocity in units of c.
    s : float
        Distance between the Sun and the observer (au).
    bm1 : float
        sqrt(1 - |v|^2): reciprocal of Lorentz factor.

    Returns
    -------
    ppr : ndarray, shape (3,)
        Proper direction to source (unit vector).
    """
    pnat = np.asarray(pnat, dtype=float)
    v    = np.asarray(v, dtype=float)

    # Scalar product of pnat and v
    pdv = pymPdp(pnat, v)

    w1 = 1.0 + pdv / (1.0 + bm1)
    w2 = SRS / s

    # Compute temporary vector p
    p = pnat * bm1 + w1 * v + w2 * (v - pdv * pnat)

    # Normalize result
    r = np.linalg.norm(p)
    if r == 0.0:
        raise ValueError("Resulting vector has zero norm, invalid input.")
    ppr = p / r
        
    return ppr


def pymAe2hd(az, el, phi):
    """
    Horizon to equatorial coordinates: transform azimuth and altitude
    to hour angle and declination.

    Parameters
    ----------
    az : float
        Azimuth (radians), north = 0, east = +pi/2.
    el : float
        Altitude (radians), informally elevation.
    phi : float
        Site latitude (radians).

    Returns
    -------
    ha : float
        Hour angle (radians), in range +/- pi.
    dec : float
        Declination (radians), in range +/- pi/2.
    """

    # Useful trig functions
    sa = np.sin(az)
    ca = np.cos(az)
    se = np.sin(el)
    ce = np.cos(el)
    sp = np.sin(phi)
    cp = np.cos(phi)

    # HA, Dec unit vector
    x = - ca * ce * sp + se * cp
    y = - sa * ce
    z =   ca * ce * cp + se * sp

    # To spherical
    r  = np.sqrt(x * x + y * y)
    ha = np.arctan2(y, x) if r != 0.0 else 0.0
    dec = np.arctan2(z, r)

    return ha, dec


def pymApcg(date1, date2, ebpv, ehp):
    """
    For a geocentric observer, prepare star-independent astrometry
    parameters for transformations between ICRS and GCRS coordinates.
    The Earth ephemeris is supplied by the caller.

    The parameters produced by this function are required in the
    parallax, light deflection, and aberration parts of the astrometric
    transformation chain.

    Parameters
    ----------
    date1, date2 : float
        TDB as 2-part Julian Date
    ebpv : ndarray shape (2,3)
        Earth barycentric position (au) and velocity (au/day)
    ehp : ndarray shape (3,)
        Earth heliocentric position (au)

    Returns
    -------
    astrom : pymASTROM
        Object holding star-independent astrometry parameters
    """
 
    # Geocentric observer: position and velocity zero
    pv = np.zeros((2, 3))

    # Call pymApcs with observer at geocenter
    astrom = pymApcs(date1, date2, pv, ebpv, ehp)

    return astrom
 

def pymApcg13(date1, date2):
    """
    For a geocentric observer, prepare star-independent astrometry parameters
    for transformations between ICRS and GCRS coordinates. The Earth ephemeris
    is obtained from SOFA models (no ephemeris input required).

    The parameters produced by this function are required in the
    parallax, light deflection, and aberration parts of the astrometric
    transformation chain.

    Parameters
    ----------
    date1, date2 : float
        TDB as 2-part Julian Date

    Returns
    -------
    astrom : pymASTROM
        Object holding star-independent astrometry parameters
    """

    # Get Earth ephemeris from SOFA models
    # ehpv: heliocentric & barycentric Earth position/velocity (au, au/day)
    ehpv, ebpv = pymEpv00(date1, date2)  # returns ehpv[2,3], ebpv[2,3]

    # Compute astrometry parameters at geocenter
    astrom = pymApcg(date1, date2, ebpv, ehpv[0])

    return astrom



def pymApci(date1, date2, ebpv, ehp, x, y, s):
    """
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and geocentric CIRS
    coordinates. Earth ephemeris and CIP/CIO are supplied by the caller.

    Parameters
    ----------
    date1, date2 : float
        TDB as 2-part Julian Date
    ebpv : ndarray shape (2,3)
        Earth barycentric position/velocity (au, au/day)
    ehp : ndarray shape (3,)
        Earth heliocentric position (au)
    x, y : float
        CIP X,Y coordinates
    s : float
        CIO locator s (radians)

    Returns
    -------
    astrom : pymASTROM
        Star-independent astrometry parameters
    """
    #astrom = pymASTROM()

    # Compute star-independent astrometry parameters for geocenter
    astrom = pymApcg(date1, date2, ebpv, ehp)

    # CIO-based bias-precession-nutation matrix
    astrom.bpn = pymC2ixys(x, y, s)

    return astrom


def pymApci13(date1, date2):
    """
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and geocentric CIRS
    coordinates. SOFA models are used to predict Earth ephemeris and CIP/CIO.

    Parameters
    ----------
    date1, date2 : float
        TDB as 2-part Julian Date

    Returns
    -------
    astrom : pymASTROM
        Star-independent astrometry parameters
    eo : float
        Equation of the origins (ERA-GST, radians)
    """
   
    #  Earth barycentric & heliocentric position/velocity
    ehpv, ebpv = pymEpv00(date1, date2)  # returns (ehpv[2,3], ebpv[2,3])

    #  Classical NPB matrix, IAU 2006/2000A
    r = pymPnm06a(date1, date2)  # 3x3 matrix

    #  Extract CIP X,Y from NPB matrix
    x, y = pymBpn2xy(r)

    #  CIO locator s
    s = pymS06(date1, date2, x, y)

    #  Compute star-independent astrometry parameters
    astrom = pymApci(date1, date2, ebpv, ehpv[0], x, y, s)

    #  Equation of the origins
    eo = pymEors(r, s)

    return astrom, eo

 
#2025-10-05
def pymApco(date1, date2, ebpv, ehp,
             x, y, s, theta, elong, phi, hm,
             xp, yp, sp, refa, refb):

    """
    Prepare star-independent astrometry parameters for a terrestrial observer.

    It computes the parameters required for transformations
    between ICRS and observed coordinates for a ground-based station, given
    the Earth ephemeris, Earth-rotation information, site coordinates,
    polar motion and refraction constants.

    Parameters
    ----------
    date1, date2 : float
        TDB as a 2-part Julian Date (not otherwise used in this function,
        but passed on to pymApcs).
    ebpv : ndarray (2,3)
        Earth barycentric position & velocity (au, au/day).
    ehp : ndarray (3,)
        Earth heliocentric position (au).
    x, y : float
        CIP X, Y coordinates (unit vector components).
    s : float
        CIO locator s (radians).
    theta : float
        Earth rotation angle (radians).
    elong : float
        Longitude (radians, east positive).
    phi : float
        Geodetic latitude (radians).
    hm : float
        Height above ellipsoid (m).
    xp, yp : float
        Polar motion coordinates (radians).
    sp : float
        TIO locator s' (radians).
    refa, refb : float
        Refraction constants A,B (radians).
    astrom : pymASTROM or None
        Optional existing astrometry object to fill. If None, a new one is created.

    Returns
    -------
    astrom : pymASTROM
        Populated astrometry parameters structure.
    """

    # Create astrom
    astrom = pymASTROM()

    # CIO-based BPN matrix (C2I(x,y,s)) 
    r_bpn = pymC2ixys(x, y, s)

    # Observer's geocentric position & velocity in CIRS (m, m/s)
    pvc = pymPvtob(elong, phi, hm, xp, yp, sp, theta)   # shape (2,3)

    # Rotate observer PV into GCRS: pv = r_bpn * pvc  (pvtob produced CIRS)
    pv = pymTrxpv(r_bpn, pvc)

    # Compute ICRS <-> GCRS parameters  
    # pymApcs should accept and return the astrom structure 
    astrom = pymApcs(date1, date2, pv, ebpv, ehp)
     
    # Form the rotation matrix r that maps CIRS -> apparent(local HA,Dec)
    # using the same sequence as the C code: Rz(theta+sp) * Ry(-xp) * Rx(-yp) * Rz(elong)
    '''
    ri = np.eye(3)
    r = pymRz(theta + sp, ri)      # Rz(theta+sp) @ I
    r = pymRy(-xp, r)              # Ry(-xp) @ previous
    r = pymRx(-yp, r)              # Rx(-yp) @ previous
    r = pymRz(elong, r)            # Rz(elong) @ previous
    '''
    #ri = pymIr()
    #r  = pymRz(elong, ri) @ pymRx(-yp, ri) @ pymRy(-xp, ri) @ pymRz(theta + sp, ri)
    
    r  = pymRz(elong) @ pymRx(-yp) @ pymRy(-xp) @ pymRz(theta + sp)
    
    # Local Earth rotation angle (eral)
    a = r[0, 0]
    b = r[0, 1]
    astrom.eral = np.arctan2(b, a) if (a != 0.0 or b != 0.0) else 0.0
    
    # Polar motion relative to local meridian: xpl, ypl
    astrom.xpl = np.arctan2(r[0, 2], np.hypot(r[0, 0], r[0, 1]))
    # ypl = -atan2( r[1,2], r[2,2] ) with safe check
    a = r[1, 2]
    b = r[2, 2]
    astrom.ypl = -np.arctan2(a, b) if (a != 0.0 or b != 0.0) else 0.0
    
    # Adjusted longitude (normalize to -pi..+pi)
    astrom.along = pymAnpm(astrom.eral - theta)
    
    # Latitude functions
    astrom.sphi = np.sin(phi)
    astrom.cphi = np.cos(phi)

    # Refraction constants
    astrom.refa = refa
    astrom.refb = refb

    # Disable diurnal aberration (redundant here)
    astrom.diurab = 0.0

    # Store the CIO-based BPN matrix in the astrom record
    astrom.bpn = np.array(r_bpn, copy=True)
   
    return astrom


def pymPvtob(elong, phi, hm, xp, yp, sp, theta):
    """
    Compute the geocentric position and velocity of a terrestrial observing station.
    
    Parameters
    ----------
    elong : float
        Longitude (radians, east positive)
    phi : float
        Geodetic latitude (radians)
    hm : float
        Height above reference ellipsoid (meters)
    xp, yp : float
        Polar motion coordinates (radians)
    sp : float
        TIO locator s' (radians)
    theta : float
        Earth rotation angle (radians)
        
    Returns
    -------
    pv : np.ndarray, shape (2,3)
        Position (meters) and velocity (meters per UT1 second) in CIRS frame.
    """
  
    OM = 1.00273781191135448 * D2PI / DAYSEC  # Earth rotation rate (rad/UT1s)

    #  Geodetic to geocentric (WGS84)
    xyz = pymGd2gc(1, elong, phi, hm)  # (3,)

    #  Apply polar motion & TIO locator
    rpm = pymPom00(xp, yp, sp)            # 3x3
    
    xyz = rpm.T @ xyz                

    #  Rotate by Earth rotation angle
    s, c = np.sin(theta), np.cos(theta)
    pos  = np.array([c * xyz[0] - s * xyz[1],
                     s * xyz[0] + c * xyz[1],
                                      xyz[2]])

    #  Velocity from Earth's rotation
    vel = OM * np.array([-s * xyz[0] - c * xyz[1],
                          c * xyz[0] - s * xyz[1],
                                            0.0])

    return np.vstack((pos, vel))


def pymGd2gce(a, f, elong, phi, height):
    """
    Transform geodetic coordinates to geocentric for a reference ellipsoid
    of specified form.

    Parameters
    ----------
    a : float
        Equatorial radius of the ellipsoid (e.g. WGS84: 6378137.0 meters).
    f : float
        Flattening of the ellipsoid (e.g. WGS84: 1/298.257223563).
    elong : float
        Geodetic longitude (radians, east positive).
    phi : float
        Geodetic latitude (radians).
    height : float
        Height above the ellipsoid (same units as a).

    Returns
    -------
    xyz : ndarray, shape (3,)
        Geocentric Cartesian coordinates (same units as a).

    Raises
    ------
    ValueError
        If the input geometry is illegal (denominator d <= 0).
    """

    sp = np.sin(phi)
    cp = np.cos(phi)

    # Flattening parameter squared
    w = (1.0 - f) ** 2

    # Denominator
    d = cp * cp + w * sp * sp
    if d <= 0.0:
        raise ValueError("illegal ellipsoid geometry (denominator <= 0)")

    # Auxiliary values
    ac = a / np.sqrt(d)      # radius in prime vertical
    at = w * ac              # radius in polar direction

    # Compute geocentric position
    r = (ac + height) * cp
    x = r * np.cos(elong)
    y = r * np.sin(elong)
    z = (at + height) * sp

    xyz = np.array([x, y, z])

    return xyz

 

def pymEform(n):
    """
    Return reference ellipsoid parameters.

    Parameters
    ----------
    n : int
        Ellipsoid identifier:
            1 : WGS84
            2 : GRS80
            3 : WGS72

    Returns
    -------
    a : float
        Equatorial radius [m]
    f : float
        Flattening

    Raises
    ------
    ValueError
        If the ellipsoid identifier is invalid.
    """
    # Reference ellipsoids dictionary
    ellipsoids = {
        1: (6378137.0, 1.0 / 298.257223563),  # WGS84
        2: (6378137.0, 1.0 / 298.257222101),  # GRS80
        3: (6378135.0, 1.0 / 298.26),         # WGS72
    }

    # Lookup ellipsoids dictionary
    params = ellipsoids.get(n)
    if params is None:
        raise ValueError("invalid ellipsoid identifier: %r (valid: 1, 2, 3)" % (n,))

    a, f = params

    return a, f



def pymGd2gc(n, elong, phi, height):
    """
    Transform geodetic coordinates to geocentric using a specified reference ellipsoid.
    
    Parameters
    ----------
    n : int
        Ellipsoid identifier (1=WGS84, 2=GRS80, 3=WGS72)
    elong : float
        Longitude (radians, east +ve)
    phi : float
        Geodetic latitude (radians)
    height : float
        Height above ellipsoid [m]
    
    Returns
    -------
    xyz : ndarray, shape (3,)
        Geocentric Cartesian coordinates [m]

    Raises
    ------
    ValueError
        If the ellipsoid identifier is invalid, or the input geometry is illegal.
    """
    # Get ellipsoid parameters (raises ValueError on a bad identifier)
    a, f = pymEform(n)

    # Transform to geocentric (raises ValueError on illegal geometry)
    xyz = pymGd2gce(a, f, elong, phi, height)

    return xyz

 

#2025-10-06
def pymApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
              phpa, tc, rh, wl):
    """
    Compute star-independent astrometry parameters for terrestrial observers.

    Parameters
    ----------
    utc1, utc2 : float
        UTC as 2-part quasi Julian Date.
    dut1 : float
        UT1 - UTC in seconds.
    elong, phi : float
        Geodetic longitude and latitude (radians).
    hm : float
        Height above ellipsoid (m).
    xp, yp : float
        Polar motion coordinates (radians).
    phpa, tc, rh, wl : float
        Pressure (hPa), temperature (C), relative humidity (0-1), wavelength (micrometers).
    astrom : pymASTROM or None
        Optional preallocated astrometry object.

    Returns
    -------
    astrom : pymASTROM
        Updated astrometry parameters.
    eo : float
        Equation of the origins (radians).
   
    """
    astrom = pymASTROM()

    # UTC -> TAI -> TT, UTC -> UT1  
    tai1, tai2 = pymUtctai(utc1, utc2)
    tt1,  tt2  = pymTaitt(tai1, tai2)
    ut11, ut12 = pymUtcut1(utc1, utc2, dut1)

    # Earth ephemeris ===
    ehpv, ebpv = pymEpv00(tt1, tt2)

    # Precession-nutation  
    r    = pymPnm06a(tt1, tt2)
    x, y = pymBpn2xy(r)
    s    = pymS06(tt1, tt2, x, y)

    # Earth rotation  
    theta = pymEra00(ut11, ut12)
    sp    = pymSp00(tt1, tt2)

    # Refraction constants  
    refa, refb = pymRefco(phpa, tc, rh, wl)

    # Star-independent astrometry parameters ===
    astrom = pymApco(tt1, tt2, ebpv, ehpv[0], x, y, s, theta,
                     elong, phi, hm, xp, yp, sp, refa, refb)

    # Equation of the origins  
    eo = pymEors(r, s)

    return astrom, eo



def pymRefco(phpa, tc, rh, wl):
    """
    Determine the constants A and B in the atmospheric refraction model
    dZ = A tan Z + B tan^3 Z.
    
    Parameters
    ----------
    phpa : float
        Pressure at the observer (hPa)
    tc : float
        Ambient temperature at the observer (deg C)
    rh : float
        Relative humidity at the observer (0-1)
    wl : float
        Wavelength (micrometers)
    
    Returns
    -------
    refa : float
        Coefficient of tan(Z) in radians
    refb : float
        Coefficient of tan^3(Z) in radians
    """
    # Optical/IR or radio case: switch at 100 microns
    optic = wl <= 100.0

    # Restrict input parameters to safe ranges using numpy.clip
    t = np.clip(tc, -150.0, 200.0)
    p = np.clip(phpa, 0.0, 10000.0)
    r = np.clip(rh, 0.0, 1.0)
    w = np.clip(wl, 0.1, 1e6)

    # Water vapor pressure at the observer
    if p > 0.0:
        ps = 10 ** ((0.7859 + 0.03477*t) / (1.0 + 0.00412*t)) * (1.0 + p*(4.5e-6 + 6e-10*t**2))
        pw = r * ps / (1.0 - (1.0 - r)*ps/p)
    else:
        pw = 0.0

    # Temperature in Kelvin
    tk = t + 273.15

    # Refractive index minus 1 at the observer
    if optic:
        wlsq = w**2
        gamma = ((77.53484e-6 + (4.39108e-7 + 3.666e-9/wlsq)/wlsq) * p - 11.2684e-6*pw) / tk
    else:
        gamma = (77.6890e-6*p - (6.3938e-6 - 0.375463/tk)*pw) / tk

    # Beta factor from Stone (with empirical adjustment for radio)
    beta = 4.4474e-6 * tk
    if not optic:
        beta -= 0.0074 * pw * beta

    # Refraction constants from Green
    refa =  gamma * (1.0 - beta)
    refb = -gamma * (beta - gamma/2.0)

    return refa, refb


# Astrometry ASTROM class
class pymASTROM:
    """
    Python version of iauASTROM structure.
    Holds star-independent astrometry parameters.
    """
    def __init__(self):
        self.pmt = 0.0
        self.eb = np.zeros(3)     # SSB to observer (au)
        self.eh = np.zeros(3)     # Sun to observer (unit vector)
        self.em = 0.0             # Distance from Sun to observer (au)
        self.v  = np.zeros(3)     # Barycentric velocity (c units)
        self.bm1 = 1.0            # Reciprocal of Lorentz factor
        self.bpn = np.identity(3) # Bias-precession-nutation matrix

        # Other unchanged parameters (placeholders)
        self.along = 0.0
        self.phi   = 0.0
        self.xpl   = 0.0
        self.ypl   = 0.0
        self.sphi  = 0.0
        self.cphi  = 1.0
        self.diurab= 0.0
        self.eral  = 0.0
        self.refa  = 0.0
        self.refb  = 0.0
 
    def print(self):
       
        eb_c = tuple([x for x in self.eb])
        eh_c = tuple([x for x in self.eh])
        v_c  = tuple([x for x in self.v])
        
        bpn_c  = np.array([x for x in self.bpn])
        print("pmt:",    self.pmt)
        print("eb:",     eb_c)
        print("eh:",     eh_c)
        print("em:",     self.em) 
        print("v:",      v_c)
        print("bm1:",    self.bm1) 
        print("bpn:",    bpn_c)
        print("along:",  self.along)   
        print("phi:",    self.phi)     
        print("xpl:",    self.xpl)    
        print("ypl:",    self.ypl)    
        print("sphi:",   self.sphi)    
        print("cphi:",   self.cphi)    
        print("diurab:", self.diurab) 
        print("eral:",   self.eral)
        print("refa:",   self.refa)  
        print("refb:",   self.refb) 
        
    def __str__(self):
         lines = [
            f"pmt: {self.pmt}",
            f"eb: {self.eb}",
            f"eh: {self.eh}",
            f"em: {self.em}",
            f"v: {self.v}",
            f"bm1: {self.bm1}",
            f"bpn:\n{self.bpn}",
            f"along: {self.along}",
            f"phi: {self.phi}",
            f"xpl: {self.xpl}",
            f"ypl: {self.ypl}",
            f"sphi: {self.sphi}",
            f"cphi: {self.cphi}",
            f"diurab: {self.diurab}",
            f"eral: {self.eral}",
            f"refa: {self.refa}",
            f"refb: {self.refb}"
         ]
         return "\n".join(lines)     

def pymApcs(date1, date2, pv, ebpv, ehp):
    """
    For an observer whose geocentric position and velocity are known,
    prepare star-independent astrometry parameters for transformations
    between ICRS and GCRS.  The Earth ephemeris is supplied by the
    caller.
    
    The parameters produced by this function are required in the space
    motion, parallax, light deflection and aberration parts of the
    astrometric transformation chain.
    
    Parameters
    ----------
    date1, date2 : float
        TDB as 2-part Julian Date
    pv : ndarray shape (2,3)
        Observer's geocentric position (m) and velocity (m/s)
    ebpv : ndarray shape (2,3)
        Earth barycentric position (au) and velocity (au/day)
    ehp : ndarray shape (3,)
        Earth heliocentric position (au)

    Returns
    -------
    astrom : pymASTROM
        Object holding star-independent astrometry parameters
        - pmt   : float       PM time interval (SSB, Julian years)
        - eb    : ndarray(3,) SSB to observer (vector, au)
        - eh    : ndarray(3,) Sun to observer (unit vector)
        - em    : float       Distance from Sun to observer (au)
        - v     : ndarray(3,) Barycentric observer velocity (vector, fraction of c)
        - bm1   : float       sqrt(1 - |v|^2): reciprocal Lorenz factor
        - bpn   : ndarray(3,3) Bias-precession-nutation matrix
        - along, xpl, ypl, sphi, cphi, diurab, eral, refa, refb
                       unchanged placeholders
    """

    AUDMS = DAU / DAYSEC     # au/day -> m/s
    CR = AULT / DAYSEC       # Light time for 1 au (days)

    pv   = np.asarray(pv, dtype=float).reshape(2, 3)
    ebpv = np.asarray(ebpv, dtype=float).reshape(2, 3)
    ehp  = np.asarray(ehp, dtype=float).reshape(3)

    astrom = pymASTROM()

    # Proper motion time interval (SSB, Julian years)
    astrom.pmt = ((date1 - DJ00) + date2) / DJY

    # Adjust Earth ephemeris to observer
    dp = pv[0] / DAU        # m -> au
    dv = pv[1] / AUDMS      # m/s -> au/day
    pb = ebpv[0] + dp
    vb = ebpv[1] + dv
    ph = ehp + dp

    # Barycentric position of observer (au)
    astrom.eb = pb

    # Heliocentric direction & distance
    astrom.em = np.linalg.norm(ph)
    astrom.eh = ph / astrom.em

    # Barycentric velocity in units of c
    astrom.v = vb * CR
    v2 = np.dot(astrom.v, astrom.v)
    astrom.bm1 = np.sqrt(1.0 - v2)

    # Reset NPB matrix = identity
    astrom.bpn = np.identity(3)

    return astrom

 
#2025-10-04
def pymApcs13(date1, date2, pv):
    """
    For an observer whose geocentric position and velocity are known,
    prepare star-independent astrometry parameters for transformations
    between ICRS and GCRS.  The Earth ephemeris is computed internally
    using the SOFA ephemeris model (pymEpv00).

    The parameters produced by this function are required in the space
    motion, parallax, light deflection and aberration parts of the
    astrometric transformation chain.

    Parameters
    ----------
    date1, date2 : float
        TDB as 2-part Julian Date
    pv : ndarray shape (2,3)
        Observer's geocentric position (m) and velocity (m/s)

    Returns
    -------
    astrom : pymASTROM
        Object holding star-independent astrometry parameters
        
        - pmt   : float       PM time interval (SSB, Julian years)
        - eb    : ndarray(3,) SSB to observer (vector, au)
        - eh    : ndarray(3,) Sun to observer (unit vector)
        - em    : float       Distance from Sun to observer (au)
        - v     : ndarray(3,) Barycentric observer velocity (vector, fraction of c)
        - bm1   : float       sqrt(1 - |v|^2): reciprocal Lorenz factor
        - bpn   : ndarray(3,3) Bias-precession-nutation matrix
        - along, xpl, ypl, sphi, cphi, diurab, eral, refa, refb
                       unchanged placeholders

    """

    # Earth barycentric & heliocentric position/velocity (au, au/day)
    ehpv, ebpv = pymEpv00(date1, date2)

    # Compute star-independent astrometry parameters
    astrom = pymApcs(date1, date2, pv, ebpv, ehpv[0])

    return astrom


def pymAper(theta, astrom):
    """
    Update only the Earth rotation angle in the astrometry parameters.

    Parameters
    ----------
    theta : float or array_like
        Earth rotation angle in radians
    astrom : ASTROM
        Star-independent astrometry parameters with 'along' attribute

    Returns
    -------
    None
        Updates astrom.eral in place
    """
    # Use numpy for potential array input
    astrom.eral = np.asarray(theta) + astrom.along
    
    return astrom


def pymAper13(ut11, ut12, astrom):
    """
    Update only the Earth rotation angle in the astrometry parameters
    using UT1 as input.

    Parameters
    ----------
    ut11, ut12 : float or array_like
        UT1 as a 2-part Julian Date
    astrom : ASTROM
        Star-independent astrometry parameters with 'along' attribute

    Returns
    -------
    None
        Updates astrom.eral in place
    """
    # Compute Earth Rotation Angle (ERA)
    theta  = pymEra00(ut11, ut12)
    # Update astrom structure
    astrom = pymAper(theta, astrom)
    
    return astrom


def pymApio(sp, theta, elong, phi, hm, xp, yp, refa, refb):
    """
    Prepare star-independent astrometry parameters for a terrestrial observer
    for transformations between CIRS and observed coordinates. 
    The caller provides the Earth orientation information, site coordinates,
    and refraction constants.

    Parameters
    ----------
    sp : float
        TIO locator s' (radians)
    theta : float
        Earth rotation angle (radians)
    elong : float
        Observer’s longitude (radians, east positive)
    phi : float
        Observer’s geodetic latitude (radians)
    hm : float
        Observer’s height above ellipsoid (meters)
    xp : float
        Polar motion coordinate xp (radians)
    yp : float
        Polar motion coordinate yp (radians)
    refa : float
        Atmospheric refraction constant A (radians)
    refb : float
        Atmospheric refraction constant B (radians)
   
    Returns
    -------
    astrom : pymASTROM
        Updated astrometry parameters including:
        - eral : local Earth rotation angle (radians)
        - along : adjusted longitude + s' (radians)
        - xpl, ypl : polar motion with respect to local meridian (radians)
        - sphi, cphi : sine and cosine of geodetic latitude
        - diurab : magnitude of diurnal aberration vector
        - refa, refb : atmospheric refraction constants (radians)
    """
    astrom = pymASTROM()
    
    # Rotation matrix from CIRS to apparent [HA,Dec]
    '''
    r = pymIr()              # identity 3x3
    pymRz(theta + sp, r)     # rotate Z
    pymRy(-xp, r)            # rotate Y
    pymRx(-yp, r)            # rotate X
    pymRz(elong, r)          # rotate Z
    '''
    #ri = pymIr()
    #r  = pymRz(elong, ri) @ pymRx(-yp, ri) @ pymRy(-xp, ri) @ pymRz(theta + sp, ri)
    r  = pymRz(elong) @ pymRx(-yp) @ pymRy(-xp) @ pymRz(theta + sp)
    
    # Local Earth rotation angle
    a, b = r[0,0], r[0,1]
    astrom.eral = np.arctan2(b, a) if (a != 0.0 or b != 0.0) else 0.0

    # Polar motion wrt local meridian
    c = r[0,2]
    astrom.xpl = np.arctan2(c, np.sqrt(a*a + b*b))
    a, b = r[1,2], r[2,2]
    astrom.ypl = -np.arctan2(a, b) if (a != 0.0 or b != 0.0) else 0.0

    # Adjusted longitude
    astrom.along = pymAnpm(astrom.eral - theta)

    # Latitude functions
    astrom.sphi = np.sin(phi)
    astrom.cphi = np.cos(phi)

    # Observer geocentric position/velocity
    pv = pymPvtob(elong, phi, hm, xp, yp, sp, theta)
    astrom.diurab = np.linalg.norm(pv[1, :2]) / CMPS

    # Refraction constants
    astrom.refa = refa
    astrom.refb = refb
    
    return astrom


def pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
              phpa, tc, rh, wl):
    """
    Prepare star-independent astrometry parameters for a terrestrial observer
    for transformations between CIRS and observed coordinates. 
    The caller provides UTC, site coordinates, ambient air conditions, and observing wavelength.

    Parameters
    ----------
    utc1 : float
        First part of UTC as a 2-part Julian Date (e.g., integer part)
    utc2 : float
        Second part of UTC as a 2-part Julian Date (e.g., fractional day)
    dut1 : float
        UT1-UTC in seconds
    elong : float
        Observer’s longitude (radians, east positive)
    phi : float
        Observer’s geodetic latitude (radians)
    hm : float
        Observer’s height above ellipsoid (meters)
    xp : float
        Polar motion coordinate xp (radians)
    yp : float
        Polar motion coordinate yp (radians)
    phpa : float
        Atmospheric pressure at observer (hPa)
    tc : float
        Ambient temperature at observer (°C)
    rh : float
        Relative humidity at observer (0–1)
    wl : float
        Observing wavelength (micrometers)
  
    Returns
    -------
    astrom : pymASTROM
        Updated astrometry parameters including:
        - eral : local Earth rotation angle (radians)
        - along : adjusted longitude + s' (radians)
        - xpl, ypl : polar motion with respect to local meridian (radians)
        - sphi, cphi : sine and cosine of geodetic latitude
        - diurab : magnitude of diurnal aberration vector
        - refa, refb : atmospheric refraction constants (radians)
    
    """
    #astrom = pymASTROM() 
    # UTC -> TAI
    tai1, tai2 = pymUtctai(utc1, utc2)
    
    # TAI -> TT
    tt1, tt2   = pymTaitt(tai1, tai2)
    
    # UTC -> UT1
    ut11, ut12 = pymUtcut1(utc1, utc2, dut1)
    
    # TIO locator s'
    sp = pymSp00(tt1, tt2)
    
    # Earth rotation angle
    theta = pymEra00(ut11, ut12)
    
    # Refraction constants
    refa, refb = pymRefco(phpa, tc, rh, wl)
    
    # Fill astrometry parameters
    astrom = pymApio(sp, theta, elong, phi,  hm,  xp,  yp,  refa,  refb)
    
    return astrom


def pymAtcc13(rc, dc, pr, pd, px, rv, date1, date2):
    """
    Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS
    astrometric place.

    Parameters:
    -----------
    rc, dc : float
        ICRS RA,Dec at J2000.0 (radians)
    pr : float
        RA proper motion (radians/year)
    pd : float
        Dec proper motion (radians/year)
    px : float
        parallax (arcsec)
    rv : float
        radial velocity (km/s, +ve if receding)
    date1, date2 : float
        TDB as a 2-part Julian Date

    Returns:
    --------
    tuple : (ra, da)
        ra : float - ICRS astrometric RA (radians)
        da : float - ICRS astrometric Dec (radians)
    """
    
    # Star-independent astrometry parameters
    astrom, eo = pymApci13(date1, date2)
    
    # Catalog ICRS (epoch J2000.0) to astrometric
    ra, da = pymAtccq(rc, dc, pr, pd, px, rv, astrom)
    
    return ra, da

def pymAtccq(rc, dc, pr, pd, px, rv, astrom):
    """
    Transformation of a star's ICRS catalog entry (epoch J2000.0)
    into ICRS astrometric place using precomputed star-independent
    astrometry parameters.

    Parameters
    ----------
    rc : float
        ICRS RA at J2000.0 (radians)
    dc : float
        ICRS Dec at J2000.0 (radians)
    pr : float
        RA proper motion (radians/year)
    pd : float
        Dec proper motion (radians/year)
    px : float
        Parallax (arcseconds)
    rv : float
        Radial velocity (km/s, positive if receding)
    astrom : pymASTROM
        Star-independent astrometry parameters.

    Returns
    -------
    ra : float
        ICRS astrometric RA (radians)
    da : float
        ICRS astrometric Dec (radians)

    """
    # BCRS coordinate vector including proper motion and parallax
    p = pymPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)

    # Convert vector to spherical coordinates (RA, Dec)
    w, da = pymC2s(p)

    # Normalize RA to 0..2pi
    ra = pymAnp(w)

    return ra, da


def pymPmpx(rc, dc, pr, pd, px, rv, pmt, pob):
    """
    Proper motion and parallax  

    Parameters
    ----------
    rc, dc : float
        ICRS right ascension and declination at catalog epoch [radians]
    pr, pd : float
        Proper motion in RA and Dec [radians/year]
    px : float
        Parallax [arcsec]
    rv : float
        Radial velocity [km/s, +ve if receding]
    pmt : float
        Proper motion time interval (Julian years)
    pob : ndarray, shape (3,)
        SSB to observer vector [au]

    Returns
    -------
    pco : ndarray, shape (3,)
        Coordinate direction (BCRS unit vector)
    """
    # Constants
    VF = DAYSEC * DJM / DAU            # Conversion factor: km/s to au/year
    AULTY = AULT / DAYSEC / DJY        # Light travel time for 1 AU in Julian years
    pxr = px * DAS2R                   # Parallax in radians
    pob = np.asarray(pob, dtype=float).reshape(3)

    # Convert spherical coordinates to Cartesian unit vector
    sr, cr = np.sin(rc), np.cos(rc)
    sd, cd = np.sin(dc), np.cos(dc)
    p = np.array([cr * cd, sr * cd, sd])

    # Proper motion time interval including approximate Roemer effect
    dt = pmt + np.dot(p, pob) * AULTY

    # Space motion components (radians/year)
    w = VF * rv * pxr
    pdz = pd * p[2]
    pm = np.array([
            -pr * p[1] - pdz * cr + w * p[0],
             pr * p[0] - pdz * sr + w * p[1],
             pd * cd + w * p[2]
    ])

    # Corrected coordinate direction including proper motion and parallax
    pco = p + dt * pm - pxr * pob

    # Normalize the resulting vector
    _, pco = pymPn(pco)
    
    return pco


def pymAtci13(rc, dc, pr, pd, px, rv, date1, date2):
    """
    Transform ICRS star data, epoch J2000.0, to CIRS coordinates.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)
    dc : float
        ICRS declination at J2000.0 (radians)
    pr : float
        Proper motion in RA (radians/year), dRA/dt (not multiplied by cos(Dec))
    pd : float
        Proper motion in Dec (radians/year)
    px : float
        Parallax (arcseconds)
    rv : float
        Radial velocity (km/s, positive if receding)
    date1 : float
        TDB as a 2-part Julian Date (part 1)
    date2 : float
        TDB as a 2-part Julian Date (part 2)

    Returns
    -------
    ri : float
        CIRS right ascension (radians)
    di : float
        CIRS declination (radians)
    eo : float
        Equation of the origins (ERA - GST, radians)

    """

    # Compute astrometry parameters for ICRS to CIRS transformation (2013 model)
    astrom, eo = pymApci13(date1, date2)

    # Apply the ICRS -> CIRS transformation
    ri, di = pymAtciq(rc, dc, pr, pd, px, rv, astrom)

    return ri, di, eo


def pymAtciq(rc, dc, pr, pd, px, rv, astrom):
    """
    Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    star-independent astrometry parameters.

    Parameters:
    -----------
    rc, dc : float
        ICRS RA,Dec at J2000.0 (radians)
    pr : float
        RA proper motion (radians/year)
    pd : float
        Dec proper motion (radians/year)
    px : float
        parallax (arcsec)
    rv : float
        radial velocity (km/s, +ve if receding)
    astrom : pymASTROM
        star-independent astrometry parameters

    Returns:
    --------
    tuple : (ri, di)
        ri : float - CIRS RA (radians)
        di : float - CIRS Dec (radians)
    """
    
    # Proper motion and parallax, giving BCRS coordinate direction
    pco = pymPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)
    
    # Light deflection by the Sun, giving BCRS natural direction
    pnat = pymLdsun(pco, astrom.eh, astrom.em)
    
    # Aberration, giving GCRS proper direction
    ppr = pymAb(pnat, astrom.v, astrom.em, astrom.bm1)
    
    # Bias-precession-nutation, giving CIRS proper direction
    pi = pymRxp(astrom.bpn, ppr)
    
    # CIRS RA,Dec
    w, di = pymC2s(pi)
    ri    = pymAnp(w)
    
    return ri, di


def pymLd(bm, p, q, e, em, dlim):
    """
    Apply light deflection by a solar-system body (NumPy version of iauLd).

    Parameters
    ----------
    bm : float
        Mass of the gravitating body (solar masses).
    p : ndarray, shape (3,)
        Direction from observer to source (unit vector).
    q : ndarray, shape (3,)
        Direction from body to source (unit vector).
    e : ndarray, shape (3,)
        Direction from body to observer (unit vector).
    em : float
        Distance from body to observer (au).
    dlim : float
        Deflection limiter (phi^2 / 2, where phi is angular separation in radians).

    Returns
    -------
    p1 : ndarray, shape (3,)
        Deflected direction from observer to source (not normalized).
    
    Notes
    -----
    1) Based on Eq. (70) in Klioner (2003) and Eq. (7.63) in the
       Explanatory Supplement (Urban & Seidelmann 2013).
    2) The returned vector is not normalized, but the deviation from
       unit magnitude is negligible.
    3) Vectors must be of unit magnitude; no validation performed.
    """

    # Compute q · (q + e)
    p = np.asarray(p, dtype=float).reshape(3)
    q = np.asarray(q, dtype=float).reshape(3)
    e = np.asarray(e, dtype=float).reshape(3)

    qpe   = q + e
    qdqpe = pymPdp(q, qpe)
    

    # Compute weight term: 2 * G * M / (em * c^2 * (q·(q + e)))
    w = bm * SRS / em / max(qdqpe, dlim)

    # Compute p × (e × q)
    eq  = pymPxp(e, q)
    peq = pymPxp(p, eq)

    # Apply light deflection correction
    p1 = p + w * peq

    return p1

def pymLdsun(p, e, em):
    """
    Deflection of starlight by the Sun.  (Python version of iauLdsun)

    Parameters
    ----------
    p : ndarray, shape (3,)
        Direction from observer to star (unit vector).
    e : ndarray, shape (3,)
        Direction from Sun to observer (unit vector).
    em : float
        Distance from Sun to observer (au).

    Returns
    -------
    p1 : ndarray, shape (3,)
        Observer to deflected star (unit vector, not normalized).

    """
    # Compute the square of the Sun-observer distance
    em2 = em * em
    if em2 < 1.0:
        em2 = 1.0

    # Deflection limiter (smaller for distant observers)
    dlim = 1e-6 / (em2 if em2 > 1.0 else 1.0)

    # Apply the solar light deflection
    p1 = pymLd(1.0, p, p, e, em, dlim)

    return p1

 

#2025-10-09
class pymLDBODY:
    """
    Equivalent to IAU SOFA structure 'iauLDBODY'.

    Attributes
    ----------
    bm : float
        Mass of the body (solar masses)
    dl : float
        Deflection limiter (radians^2 / 2)
    pv : np.ndarray, shape (2,3)
        Barycentric position and velocity of the body:
        pv[0] = position (au)
        pv[1] = velocity (au/day)
    """
    def __init__(self, bm, dl, pv):
        self.bm = float(bm)
        self.dl = float(dl)
        self.pv = np.array(pv, dtype=float).reshape(2, 3)


def pymLdn(n, bodies, ob, sc):
    """
    Apply light deflection caused by multiple solar-system bodies
    to the observed direction of a star.

    Parameters
    ----------
    n : int
        Number of deflecting bodies.
    bodies : list[pymLDBODY]
        List of body objects with mass, deflection limiter, and PV data.
    ob : array-like, shape (3,)
        Observer barycentric position (au).
    sc : array-like, shape (3,)
        Star direction unit vector (observer to star).

    Returns
    -------
    np.ndarray, shape (3,)
        Deflected star direction (not normalized, matching iauLdn).
    """
    # Physical constants 
    
    CR = AULT / DAYSEC    # Light time for 1 au (days)

    ob = np.asarray(ob, dtype=float).reshape(3)
    sn = np.asarray(sc, dtype=float).reshape(3)
     
    if n == 0:
        return sn.copy()

    # Process each body 
    for body in bodies:
        # Body barycentric position and velocity
        rb = body.pv[0]
        vb = body.pv[1]

        # Body-to-observer vector
        v = ob - rb

        # Time correction
        dt = np.dot(sn, v) * CR
        dt = min(dt, 0.0)

        # Backtrack body to time light passed
        ev = v - dt * vb

        # Normalize body-observer vector
        em = np.linalg.norm(ev)
        if em == 0.0:
            continue
        e = ev / em

        # Apply deflection
        sn = pymLd(body.bm, sn, sn, e, em, body.dl)
        

    return sn

 
def pymAtciqn(rc, dc, pr, pd, px, rv, astrom, n, bodies):
    """
    Quick ICRS (J2000.0) to CIRS transformation with multiple light-deflecting bodies.

    Parameters
    ----------
    rc, dc : float
        ICRS RA, Dec at J2000.0 (radians)
    pr, pd : float
        Proper motion in RA and Dec (radians/year)
    px : float
        Parallax (arcsec)
    rv : float
        Radial velocity (km/s, positive if receding)
    astrom : object
        Star-independent astrometry parameters (with attributes:
        pmt, eb, v, em, bm1, bpn)
    bodies : list of pymLDBODY
        Light-deflecting bodies

    Returns
    -------
    ri, di : float
        CIRS RA, Dec (radians)
    """
    #  Proper motion & parallax -> BCRS coordinate direction
    pco = pymPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)

    #  Light deflection -> BCRS natural direction
    pnat = pymLdn(n, bodies, astrom.eb, pco)

    #  Aberration -> GCRS proper direction
    ppr = pymAb(pnat, astrom.v, astrom.em, astrom.bm1)

    #  Bias-precession-nutation -> CIRS proper direction
    pi = pymRxp(astrom.bpn, ppr)

    #  Convert to spherical coordinates (RA, Dec)
    w, di = pymC2s(pi)

    #  Normalize RA to 0–2π
    ri = pymAnp(w)

    return ri, di


#2025-10-09
def pymAtciqz(rc, dc, astrom):
    """
    Quick ICRS -> CIRS transformation (zero parallax, zero proper motion).
    
    Parameters
    ----------
    rc, dc : float
        ICRS right ascension and declination (radians)
    astrom : object
        Astrometry parameters, with attributes:
            eh   : ndarray(3,)   Unit vector Sun->observer
            em   : float         Distance Sun->observer (au)
            v    : ndarray(3,)   Barycentric observer velocity (in units of c)
            bm1  : float         sqrt(1 - |v|^2)
            bpn  : ndarray(3,3)  Bias-precession-nutation matrix

    Returns
    -------
    ri, di : float
        CIRS right ascension and declination (radians)
    """

    #  ICRS coordinate direction (unit vector, BCRS)
    pco = pymS2c(rc, dc)

    #  Light deflection due to the Sun
    pnat = pymLdsun(pco, astrom.eh, astrom.em)

    #  Apply stellar aberration
    ppr = pymAb(pnat, astrom.v, astrom.em, astrom.bm1)

    #  Apply bias-precession-nutation matrix
    pi = pymRxp(astrom.bpn, ppr)

    #  Convert to spherical coordinates
    ri, di = pymC2s(pi)

    #  Normalize RA to [0, 2π)
    ri = pymAnp(ri)

    return ri, di


def pymAtco13(rc, dc, pr, pd, px, rv, 
              utc1, utc2, dut1,
              elong, phi, hm, xp, yp,
              phpa, tc, rh, wl):
    """
    Convert ICRS RA, Dec to observed coordinates using UTC, site data,
    atmospheric conditions, and observing wavelength.

    Parameters
    ----------
    rc, dc : float
        ICRS right ascension and declination (radians, epoch J2000.0)
    pr, pd : float
        Proper motion in RA and Dec (radians/year)
    px : float
        Parallax (arcseconds)
    rv : float
        Radial velocity (km/s, positive if receding)
    utc1, utc2 : float
        UTC as 2-part Julian Date
    dut1 : float
        UT1 - UTC (seconds)
    elong, phi : float
        Site longitude (east +ve, radians) and geodetic latitude (radians)
    hm : float
        Site height above ellipsoid (meters)
    xp, yp : float
        Polar motion coordinates (radians)
    phpa : float
        Atmospheric pressure at observer (hPa)
    tc : float
        Temperature (°C)
    rh : float
        Relative humidity (0–1)
    wl : float
        Observing wavelength (micrometers)

    Returns
    -------
    aob : float
        Observed azimuth (radians: N=0, E=90)
    zob : float
        Observed zenith distance (radians)
    hob : float
        Observed hour angle (radians)
    dob : float
        Observed declination (radians)
    rob : float
        Observed right ascension (CIO-based, radians)
    eo : float
        Equation of the origins (ERA - GST, radians)
    """

    # Compute star-independent astrometric parameters
    astrom, eo = pymApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                           phpa, tc, rh, wl)

    # Convert ICRS -> CIRS
    ri, di = pymAtciq(rc, dc, pr, pd, px, rv, astrom)

    # Convert CIRS -> observed
    aob, zob, hob, dob, rob = pymAtioq(ri, di, astrom)

    return aob, zob, hob, dob, rob, eo


def pymAtic13(ri, di, date1, date2):
    """
    Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

    Parameters
    ----------
    ri, di : float
        CIRS geocentric RA, Dec (radians)
    date1, date2 : float
        TDB as 2-part Julian Date (any convenient split)

    Returns
    -------
    rc, dc : float
        ICRS astrometric RA, Dec (radians)
    eo : float
        Equation of the origins (ERA - GST, radians)
    """

    # Compute star-independent astrometry parameters
    astrom, eo = pymApci13(date1, date2)

    # Quick CIRS to ICRS astrometric transformation
    rc, dc = pymAticq(ri, di, astrom)

    return rc, dc, eo



def pymAticq(ri, di, astrom):
    """
    Quick CIRS RA,Dec to ICRS astrometric place using star-independent astrometry.

    Parameters
    ----------
    ri, di : float
        CIRS right ascension and declination (radians)
    astrom : object
        Star-independent astrometry parameters, must contain:
            bpn  : ndarray(3,3), bias-precession-nutation matrix
            v    : ndarray(3,), observer barycentric velocity (unit c)
            em   : float, distance from Sun to observer (au)
            bm1  : float, reciprocal Lorentz factor sqrt(1-|v|^2)
            eh   : ndarray(3,), Sun to observer unit vector

    Returns
    -------
    rc, dc : float
        ICRS astrometric right ascension and declination (radians)
    """

    #  CIRS RA,Dec to Cartesian  
    pi = pymS2c(ri, di)

    #  Bias-precession-nutation, GCRS proper direction  
    ppr = astrom.bpn.T @ pi  # transpose of BPN times vector

    #  Aberration correction (iterative)  
    d = np.zeros(3)
    for _ in range(2):
        before = ppr - d
        before /= np.linalg.norm(before)
        after = pymAb(before, astrom.v, astrom.em, astrom.bm1)
        d = after - before
    pnat = (ppr - d)
    pnat /= np.linalg.norm(pnat)

    #  Light deflection by Sun (iterative)  
    d = np.zeros(3)
    for _ in range(5):
        before = pnat - d
        before /= np.linalg.norm(before)
        after = pymLdsun(before, astrom.eh, astrom.em)
        d = after - before
    pco = pnat - d
    pco /= np.linalg.norm(pco)

    #  Cartesian to ICRS RA,Dec 
    rc_tmp, dc = pymC2s(pco)
    rc = np.mod(rc_tmp, 2*np.pi)

    return rc, dc


def pymAticqn(ri, di, astrom, n, bodies):
    """
    Quick CIRS to ICRS astrometric place transformation with n light-deflecting bodies.

    Parameters
    ----------
    ri, di : float
        CIRS RA, Dec (radians)
    astrom : object
        Star-independent astrometry parameters
    n : int
        Number of light-deflecting bodies
    bodies : list of LDBODY
        List of light-deflecting bodies

    Returns
    -------
    rc, dc : float
        ICRS astrometric RA, Dec (radians)
    """
    
    #  CIRS RA,Dec to Cartesian
    pi = pymS2c(ri, di)

    #  Apply bias-precession-nutation to get GCRS proper direction
    ppr = pymTrxp(astrom.bpn, pi)

    #  Stellar aberration (iterative, 2 iterations)
    d = np.zeros(3)
    for _ in range(2):
        before = ppr - d
        before /= np.linalg.norm(before)
        after = pymAb(before, astrom.v, astrom.em, astrom.bm1)
        d = after - before
        pnat = ppr - d
        pnat /= np.linalg.norm(pnat)

    #  Light deflection by n bodies (iterative, 5 iterations)
    d[:] = 0.0
    for _ in range(5):
        before = pnat - d
        before /= np.linalg.norm(before)
        after = pymLdn(n, bodies, astrom.eb, before)
        d = after - before
        pco = pnat - d
        pco /= np.linalg.norm(pco)

    #  Cartesian to spherical
    w, dc = pymC2s(pco)
    rc = pymAnp(w)

    return rc, dc


#2025-10-11
def pymAtio13(ri, di, utc1, utc2, dut1,
              elong, phi, hm, xp, yp,
              phpa, tc, rh, wl):
    """
    CIRS RA, Dec to observed place.  The caller supplies UTC, site
    coordinates, ambient air conditions and observing wavelength.

    Parameters
    ----------
    ri, di : float
        CIRS right ascension and declination (radians)
    utc1, utc2 : float
        UTC as a two-part Julian Date
    dut1 : float
        UT1 - UTC (seconds)
    elong : float
        Observer’s east-positive longitude (radians)
    phi : float
        Geodetic latitude (radians)
    hm : float
        Height above WGS84 ellipsoid (m)
    xp, yp : float
        Polar motion coordinates (radians)
    phpa : float
        Atmospheric pressure (hPa)
    tc : float
        Ambient temperature (°C)
    rh : float
        Relative humidity (0–1)
    wl : float
        Wavelength (micrometers)

    Returns
    -------
    aob : float
        Observed azimuth (radians, N=0,E=90)
    zob : float
        Observed zenith distance (radians)
    hob : float
        Observed hour angle (radians)
    dob : float
        Observed declination (radians)
    rob : float
        Observed right ascension (radians)
    """

    #  Compute CIRS => observed astrometry parameters
    astrom = pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                       phpa, tc, rh, wl)

    #  Transform CIRS (ri, di)  => observed (az, zd, ha, dec, ra)
    aob, zob, hob, dob, rob = pymAtioq(ri, di, astrom)

    #  Return computed observed coordinates
    return aob, zob, hob, dob, rob


def pymAtioq(ri, di, astrom):
    """
    Quick CIRS -> observed place transformation using pymS2c and pymC2s

    Parameters
    ----------
    ri, di : float
        CIRS right ascension and declination (radians)
    astrom : object
        Star-independent astrometry parameters, must contain:
            eral   : float, local Earth rotation angle
            xpl    : float, polar motion x (radians)
            ypl    : float, polar motion y (radians)
            sphi   : float, sine of geodetic latitude
            cphi   : float, cosine of geodetic latitude
            diurab : float, magnitude of diurnal aberration vector
            refa   : float, refraction constant A (radians)
            refb   : float, refraction constant B (radians)

    Returns
    -------
    aob, zob, hob, dob, rob : float
        Observed azimuth, zenith distance, hour angle, declination,
        and CIO-based right ascension (radians)
    """

    # Cartesian -HA,Dec using pymS2c 
    ha = ri - astrom.eral
    v = pymS2c(ha, di)
    x, y, z = v

    # Polar motion 
    sx, cx = np.sin(astrom.xpl), np.cos(astrom.xpl)
    sy, cy = np.sin(astrom.ypl), np.cos(astrom.ypl)

    xhd = cx * x + sx * z
    yhd = sx * sy * x + cy * y - cx * sy * z
    zhd = -sx * cy * x + sy * y + cx * cy * z

    #  Diurnal aberration  
    f = 1.0 - astrom.diurab * yhd
    xhdt = f * xhd
    yhdt = f * (yhd + astrom.diurab)
    zhdt = f * zhd

    #  Rotate to Azimuth/Altitude coordinates  
    xaet = astrom.sphi * xhdt - astrom.cphi * zhdt
    yaet = yhdt
    zaet = astrom.cphi * xhdt + astrom.sphi * zhdt

    #  Azimuth (N=0, E=90)  
    aob = np.arctan2(yaet, -xaet) if (xaet != 0 or yaet != 0) else 0.0

    #  Refraction  
    CELMIN, SELMIN = 1e-6, 0.05
    r = max(np.hypot(xaet, yaet), CELMIN)
    z = max(zaet, SELMIN)

    tz = r / z
    w = astrom.refb * tz**2
    del_ref = (astrom.refa + w) * tz / (1 + (astrom.refa + 3*w) / z**2)

    cosdel = 1.0 - 0.5 * del_ref**2
    f = cosdel - del_ref * z / r
    xaeo = xaet * f
    yaeo = yaet * f
    zaeo = cosdel * zaet + del_ref * r

    #  Observed HA,Dec vector  
    v_obs = np.array([
        astrom.sphi * xaeo + astrom.cphi * zaeo,
        yaeo,
        -astrom.cphi * xaeo + astrom.sphi * zaeo
    ])

    # Use pymC2s to get spherical coordinates (HA, Dec)
    hmobs, dcobs = pymC2s(v_obs)

    # Observed hour angle and declination
    hob = -hmobs
    dob = dcobs

    # Observed RA (CIO-based)
    rob = np.mod(astrom.eral + hmobs, 2*np.pi)

    # Observed zenith distance
    zob = np.arctan2(np.hypot(xaeo, yaeo), zaeo)

    # Normalize azimuth
    aob = np.mod(aob, 2*np.pi)

    return aob, zob, hob, dob, rob


def pymAtoc13(coord_type, ob1, ob2, utc1, utc2, dut1,
              elong, phi, hm, xp, yp, phpa, tc, rh, wl):
    """
    Converts observed coordinates (Az/ZD, HA/Dec, or RA/Dec)
    to ICRS astrometric coordinates.

    Parameters
    ----------
    coord_type : str
        Type of observed coordinates ('A', 'H', or 'R').
    ob1, ob2 : float
        Observed coordinates in radians (Az/ZD, HA/Dec, or RA/Dec).
    utc1, utc2 : float
        UTC as 2-part Julian Date.
    dut1 : float
        UT1 - UTC (seconds).
    elong, phi : float
        Observer's geodetic longitude and latitude (radians).
    hm : float
        Height above ellipsoid (m).
    xp, yp : float
        Polar motion coordinates (radians).
    phpa : float
        Pressure (hPa).
    tc : float
        Temperature (Celsius).
    rh : float
        Relative humidity (0–1).
    wl : float
        Wavelength (micrometers).

    Returns
    -------
    rc, dc : tuple of float
        ICRS right ascension and declination (radians).
    """

    #  Compute astrometry parameters
    astrom, eo = pymApco13(utc1, utc2, dut1, elong, phi, hm,
                           xp, yp, phpa, tc, rh, wl)

    #  Transform observed -> CIRS
    ri, di = pymAtoiq(coord_type, ob1, ob2, astrom)

    #  Transform CIRS -> ICRS
    rc, dc = pymAticq(ri, di, astrom)

    return rc, dc


def pymAtoi13(coord_type, ob1, ob2,  utc1, utc2, dut1,
              elong, phi, hm, xp, yp, phpa, tc, rh, wl) :
    """
    Observed place to CIRS. The caller supplies UTC, site coordinates,
    ambient air conditions and observing wavelength.

    Parameters:
    -----------
    coord_type : str
        Type of coordinates - "R", "H" or "A"
        - "R": observed right ascension and declination
        - "H": observed hour angle (west +ve) and declination
        - "A": observed azimuth (north zero, east 90 deg) and zenith distance
    ob1, ob2 : float
        Observed coordinates in radians (Az/ZD, HA/Dec, or RA/Dec).
    utc1, utc2 : float
        UTC as 2-part Julian Date.
    dut1 : float
        UT1 - UTC (seconds).
    elong, phi : float
        Observer's geodetic longitude and latitude (radians).
    hm : float
        Height above ellipsoid (m).
    xp, yp : float
        Polar motion coordinates (radians).
    phpa : float
        Pressure (hPa).
    tc : float
        Temperature (Celsius).
    rh : float
        Relative humidity (0–1).
    wl : float
        Wavelength (micrometers).
        
    Returns:
    --------
    ri 
        CIRS right ascension (CIO-based, radians)
    di 
        CIRS declination (radians)
    """
    # Star-independent astrometry parameters for CIRS->observed
    astrom = pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                              phpa, tc, rh, wl)
 
    # Transform observed to CIRS
    ri, di = pymAtoiq(coord_type, ob1, ob2, astrom)
 
    return  ri, di 
 

def pymAtoiq(coord_type, ob1, ob2, astrom):
    """
    Quick observed place to CIRS (ICRS -> CIRS inverse).
    Equivalent to SOFA iauAtoiq, but written in Python with NumPy.

    Parameters
    ----------
    coord_type : str
        'R', 'H', or 'A'
    ob1, ob2 : float
        Observed coordinates (radians)
    astrom : pymASTROM
        Star-independent astrometry parameters.

    Returns
    -------
    ri, di : float
        CIRS right ascension and declination (radians)
    """

    SELMIN = 0.05  # Minimum sin(alt) for refraction

    sphi = astrom.sphi
    cphi = astrom.cphi
    c = coord_type[0].upper()
    if c not in ['R', 'H', 'A']:
        c = 'A'

    # Case A: Az, ZD
    if c == 'A':
        ce = np.sin(ob2)
        xaeo = -np.cos(ob1) * ce
        yaeo = np.sin(ob1) * ce
        zaeo = np.cos(ob2)
    else:
        # If RA,Dec -> HA,Dec
        c1 = astrom.eral - ob1 if c == 'R' else ob1

        # Convert to Cartesian (-HA,Dec)
        v = pymS2c(-c1, ob2)
        xmhdo, ymhdo, zmhdo = v

        # Convert to Cartesian Az,El (S=0,E=90)
        xaeo = sphi * xmhdo - cphi * zmhdo
        yaeo = ymhdo
        zaeo = cphi * xmhdo + sphi * zmhdo

    # Observed azimuth and zenith distance
    az = np.arctan2(yaeo, xaeo) if (xaeo != 0 or yaeo != 0) else 0.0
    sz = np.hypot(xaeo, yaeo)
    zdo = np.arctan2(sz, zaeo)

    # Refraction correction
    refa = astrom.refa
    refb = astrom.refb
    tz = sz / (zaeo if zaeo > SELMIN else SELMIN)
    dref = (refa + refb * tz**2) * tz
    zdt = zdo + dref

    # To Cartesian Az,ZD
    ce = np.sin(zdt)
    xaet = np.cos(az) * ce
    yaet = np.sin(az) * ce
    zaet = np.cos(zdt)

    # Cartesian Az,ZD -> Cartesian -HA,Dec
    xmhda = sphi * xaet + cphi * zaet
    ymhda = yaet
    zmhda = -cphi * xaet + sphi * zaet

    # Diurnal aberration
    f = 1.0 + astrom.diurab * ymhda
    xhd = f * xmhda
    yhd = f * (ymhda - astrom.diurab)
    zhd = f * zmhda

    # Polar motion
    sx, cx = np.sin(astrom.xpl), np.cos(astrom.xpl)
    sy, cy = np.sin(astrom.ypl), np.cos(astrom.ypl)

    v = np.zeros(3)
    v[0] = cx * xhd + sx * sy * yhd - sx * cy * zhd
    v[1] = cy * yhd + sy * zhd
    v[2] = sx * xhd - cx * sy * yhd + cx * cy * zhd

    # Convert to spherical coordinates (-HA, Dec)
    hma, di = pymC2s(v)

    # Right ascension (CIRS)
    #ri = (astrom.eral + hma) % (2 * np.pi)
    ri = pymAnp(astrom.eral + hma)

    return ri, di



def pymHd2ae(ha, dec, phi):
    """
    Convert equatorial coordinates (hour angle, declination)
    to horizon coordinates (azimuth, altitude).

    Python translation of the IAU SOFA function iauHd2ae.

    Parameters
    ----------
    ha : float
        Local hour angle [radians]
    dec : float
        Declination [radians]
    phi : float
        Observer’s geodetic latitude [radians]

    Returns
    -------
    az : float
        Azimuth [radians], range [0, 2π), north=0, east=π/2
    el : float
        Altitude (elevation) [radians], range [-π/2, +π/2]

    """

    # Useful trig functions
    sh, ch = np.sin(ha),  np.cos(ha)
    sd, cd = np.sin(dec), np.cos(dec)
    sp, cp = np.sin(phi), np.cos(phi)

    # Cartesian horizon coordinates
    x = -ch * cd * sp + sd * cp
    y = -sh * cd
    z = ch * cd * cp + sd * sp

    # Convert to spherical (azimuth, elevation)
    r = np.hypot(x, y)  # sqrt(x^2 + y^2)
    a = np.arctan2(y, x) if r != 0.0 else 0.0

    # Normalize azimuth to [0, 2π)
    #az = a if a >= 0.0 else a + 2.0 * np.pi
    az = a if a >= 0.0 else a +  D2PI
    el = np.arctan2(z, r)

    return az, el



def pymHd2pa(ha, dec, phi):
    """
    Parallactic angle for a given hour angle and declination.

    Python translation of the IAU SOFA function iauHd2pa.

    Parameters
    ----------
    ha : float
        Hour angle [radians]
    dec : float
        Declination [radians]
    phi : float
        Observer's geodetic latitude [radians]

    Returns
    -------
    pa : float
        Parallactic angle [radians], range (-π, +π)

    Notes
    -----
    The parallactic angle is the position angle of the vertical,
    i.e. the angle between the directions to the north celestial pole
    and to the zenith.
    """

    cp   = np.cos(phi)
    sqsz = cp * np.sin(ha)
    cqsz = np.sin(phi) * np.cos(dec) - cp * np.sin(dec) * np.cos(ha)

    # Handle the pole case: if both components are zero, return 0.0
    pa = np.arctan2(sqsz, cqsz) if (sqsz != 0.0 or cqsz != 0.0) else 0.0

    return pa

 
#2025-10-20
pym_pm_safe_msg = {
    -1: "System error",
     1: "Distance overridden",
     2: "Excessive velocity",
     4: "Solution didn't converge"
}

def pymPmsafe(ra1, dec1, pmr1, pmd1, px1, rv1,
              ep1a, ep1b, ep2a, ep2b):
    """
    Update star catalog data for space motion, with special handling for zero parallax.

    Parameters
    ----------
    ra1 : float
        Right ascension at the initial epoch (radians)
    dec1 : float
        Declination at the initial epoch (radians)
    pmr1 : float
        Proper motion in RA (radians/year)
    pmd1 : float
        Proper motion in Dec (radians/year)
    px1 : float
        Parallax (arcseconds)
    rv1 : float
        Radial velocity (km/s, +ve = receding)
    ep1a, ep1b : float
        "Before" epoch, split into two parts (Julian Dates)
    ep2a, ep2b : float
        "After" epoch, split into two parts (Julian Dates)

    Returns
    -------
    ra2, dec2, pmr2, pmd2, px2, rv2 : float
        Updated star parameters at the new epoch

    Raises
    ------
    ValueError
        If distance is overridden, velocity is excessive, or solution didn't converge.
    """

    # Minimum allowed parallax (arcsec)
    PXMIN = 5e-7

    # Factor giving maximum allowed transverse speed ~0.1 c
    F = 326.0

    # Proper motion in one year (radians)
    pm = pymSeps(ra1, dec1, ra1 + pmr1, dec1 + pmd1)

    # Adjust parallax to avoid warnings
    px_adj = max(px1, pm * F, PXMIN)
    status = 0
    if px_adj != px1:
        status |= 1  # Distance overridden

    # Compute updated star parameters
    #try:
    ra2, dec2, pmr2, pmd2, px2, rv2 = pymStarpm(
            ra1, dec1, pmr1, pmd1, px_adj, rv1,
            ep1a, ep1b, ep2a, ep2b
        )
    #except Exception as e:
    #    raise ValueError(pym_pm_safe_msg[-1]) from e

    # Raise warnings as ValueError if any
     
    if status != 0:
        messages = [msg for idx, msg in pym_pm_safe_msg.items() if idx > 0 and 
                    (status & idx)]
        raise ValueError("; ".join(messages))

    return ra2, dec2, pmr2, pmd2, px2, rv2



def pymTpors(xi, eta, a, b):
    """
    Tangent plane projection: given star rectangular coordinates (xi, eta)
    and its spherical coordinates (a, b), determine the spherical coordinates
    of the tangent point.

    Parameters
    ----------
    xi, eta : float
        Rectangular coordinates of star image (radians).
    a, b : float
        Star's spherical coordinates (radians).

    Returns
    -------
    a01, b01, a02, b02 : float or None
        Tangent point's spherical coordinates, solution 1 and 2.
        If no solution exists, returns None for all.
    """

    xi2 = xi**2
    r = np.sqrt(1.0 + xi2 + eta**2)
    sb = np.sin(b)
    cb = np.cos(b)
    rsb = r * sb
    rcb = r * cb
    w2 = rcb**2 - xi2

    if w2 < 0.0:
        # No solution
        return None, None, None, None

    w = np.sqrt(w2)
    s = rsb - eta * w
    c = rsb * eta + w
    if xi == 0.0 and w == 0.0:
        w = 1.0
    a01 = pymAnp(a - np.arctan2(xi, w))
    b01 = np.arctan2(s, c)

    # Second solution
    w = -w
    s = rsb - eta * w
    c = rsb * eta + w
    a02 = pymAnp(a - np.arctan2(xi, w))
    b02 = np.arctan2(s, c)

    return a01, b01, a02, b02


def pymTporv(xi, eta, v):
    """
    Tangent plane projection: given star rectangular coordinates (xi, eta)
    and its direction cosines, determine the direction cosines of the tangent point.

    Parameters
    ----------
    xi, eta : float
        Rectangular coordinates of star image (radians)
    v : array_like, shape (3,)
        Star's direction cosines (unit vector)

    Returns
    -------
    v01, v02 : ndarray, shape (3,) or None
        Tangent point's direction cosines, solution 1 and 2.
        If no solution exists, returns (None, None)
    """

    x, y, z = v
    rxy = np.hypot(x, y)
    xi2 = xi**2
    eta2p1 = eta**2 + 1.0
    r = np.sqrt(xi2 + eta2p1)
    rsb = r * z
    rcb = r * rxy
    w2 = rcb**2 - xi2

    if w2 <= 0.0:
        return None, None

    sqrt_w2 = np.sqrt(w2)
    denom = eta2p1 * np.sqrt(rxy**2 * (w2 + xi2))

    # First solution
    w = sqrt_w2
    c = (rsb*eta + w) / denom
    v01 = c * np.array([x*w + y*xi, y*w - x*xi, 0.0])
    v01[2] = (rsb - eta*w) / eta2p1

    # Second solution
    w = -sqrt_w2
    c = (rsb*eta + w) / denom
    v02 = c * np.array([x*w + y*xi, y*w - x*xi, 0.0])
    v02[2] = (rsb - eta*w) / eta2p1

    return v01, v02


#2025-10-21
def pymTpsts(xi, eta, a0, b0):
    """
    Tangent plane projection (gnomonic): given the star's rectangular coordinates
    (xi, eta) and the spherical coordinates of the tangent point (a0, b0),
    compute the star's spherical coordinates (a, b).

    Parameters
    ----------
    xi, eta : float
        Rectangular coordinates of star image (radians).
    a0, b0 : float
        Spherical coordinates (radians) of the tangent point.

    Returns
    -------
    a, b : float
        Spherical coordinates (radians) of the star.
    """
    sb0, cb0 = np.sin(b0), np.cos(b0)

    # Compute intermediate quantities
    d = cb0 - eta * sb0
    r = np.sqrt(xi**2 + d**2)

    # Compute spherical coordinates
    a = np.arctan2(xi, d) + a0
    b = np.arctan2(sb0 + eta * cb0, r)

    # Normalize longitude a into [0, 2π)
    #a = np.mod(a, 2.0 * np.pi)
    a = np.mod(a, D2PI)
    
    return a, b


def pymTpstv(xi, eta, v0):
    """
    Tangent plane projection (gnomonic): given the star's rectangular coordinates (xi, eta)
    and the direction cosines of the tangent point v0, compute the direction cosines of the star.

    Parameters
    ----------
    xi, eta : float
        Rectangular coordinates of star image (radians).
    v0 : array_like, shape (3,)
        Tangent point's direction cosines (unit vector).

    Returns
    -------
    v : ndarray, shape (3,)
        Star's direction cosines (unit vector).
    """
    x, y, z = v0
    r = np.hypot(x, y)

    # Handle polar case
    if r == 0.0:
        r = 1e-20
        x = r

    # Star vector length in tangent plane
    f = np.sqrt(1.0 + xi**2 + eta**2)

    # Apply the transformation
    v = np.empty(3)
    v[0] = (x - (xi*y + eta*x*z) / r) / f
    v[1] = (y + (xi*x - eta*y*z) / r) / f
    v[2] = (z + eta*r) / f

    return v
 

pym_tpxes_msg = {
    0: 'OK',
    1: 'star too far from axis',
    2: 'antistar on tangent plane',
    3: 'antistar too far from axis'
}

def pymTpxes(a, b, a0, b0):
    """
    Tangent plane (gnomonic) projection: compute rectangular coordinates
    of a star given its spherical coordinates and the tangent point.

    Parameters
    ----------
    a, b : float
        Star's spherical coordinates (radians).
    a0, b0 : float
        Tangent point's spherical coordinates (radians).

    Returns
    -------
    xi, eta : float
        Rectangular coordinates of the star in the tangent plane (radians).

    Raises
    ------
    ValueError
        If the status is 1, 2, or 3; the message describes the problem.
    """
    TINY = 1e-6

    # Precompute trigonometric functions
    sb,  cb  = np.sin(b),  np.cos(b)
    sb0, cb0 = np.sin(b0), np.cos(b0)
    sda, cda = np.sin(a - a0), np.cos(a - a0)

    # Reciprocal of star vector length to tangent plane
    d = sb*sb0 + cb*cb0*cda

    # Determine status code
    if d > TINY:
        j = 0
    elif d >= 0.0:
        j = 1
        d = TINY
    elif d > -TINY:
        j = 2
        d = -TINY
    else:
        j = 3

    # Compute tangent plane coordinates
    xi  = cb * sda / d
    eta = (sb*cb0 - cb*sb0*cda) / d

    # Raise ValueError if status is warning/error
    if j != 0:
        raise ValueError(f"{pym_tpxes_msg[j]} (status={j})")

    return xi, eta


pym_tpxev_msg = pym_tpxes_msg 
 
def pymTpxev(v, v0):
    """
    Tangent plane (gnomonic) projection: compute rectangular coordinates
    of a star given its direction cosines and the tangent point.

    Parameters
    ----------
    v : array_like, shape (3,)
        Direction cosines of the star.
    v0 : array_like, shape (3,)
        Direction cosines of the tangent point.

    Returns
    -------
    xi, eta : float
        Rectangular coordinates of the star in the tangent plane.

    Raises
    ------
    ValueError
        If the star is too far or antistar cases occur.
    """
    TINY = 1e-6
    v  = np.asarray(v,  dtype=float)
    v0 = np.asarray(v0, dtype=float)

    # Tangent point vector length in xy-plane
    r = np.linalg.norm(v0[:2])
    if r == 0.0:
        r = 1e-20
        v0[0] = r

    # Reciprocal of star vector length to tangent plane
    w = np.dot(v[:2], v0[:2])
    d = w + v[2]*v0[2]

    # Determine status
    if d > TINY:
        j = 0
    elif d >= 0.0:
        j = 1
        d = TINY
    elif d > -TINY:
        j = 2
        d = -TINY
    else:
        j = 3

    # Compute tangent plane coordinates
    d *= r
    xi = (v[1]*v0[0] - v[0]*v0[1]) / d
    eta = (v[2]*r**2 - v0[2]*w) / d

    if j != 0:
        raise ValueError(f"{pym_tpxev_msg[j]} (status={j})")

    return xi, eta


def pymEceq06(date1, date2, dl, db):
    """
    Transform ecliptic coordinates (mean equinox and ecliptic of date)
    to ICRS RA, Dec using IAU 2006 precession model, using pym* utilities.

    Parameters
    ----------
    date1, date2 : float
        TT as 2-part Julian date
    dl, db : float
        Ecliptic longitude and latitude in radians

    Returns
    -------
    dr, dd : float
        ICRS right ascension and declination in radians
    """
    #  Spherical to Cartesian using pymS2c
    v1 = pymS2c(dl, db)

    #  Rotation matrix: ICRS equatorial to ecliptic using pymEcm06
    rm = pymEcm06(date1, date2)

    #  Transformation from ecliptic to ICRS (transpose multiplication)
    v2 = rm.T @ v1

    #  Cartesian to spherical using pymC2s
    a, b = pymC2s(v2)

    #  Normalize angles using pymAnp and pymAnpm
    dr = pymAnp(a)   # 0 to 2pi
    dd = pymAnpm(b)  # -pi to +pi

    return dr, dd


def pymEcm06(date1, date2):
    """
    ICRS equatorial to ecliptic rotation matrix, IAU 2006.

    Parameters
    ----------
    date1, date2 : float
        TT as a 2-part Julian Date.

    Returns
    -------
    rm : ndarray, shape (3,3)
        Rotation matrix from ICRS to ecliptic coordinates of the date.
    """
    # Mean obliquity of the ecliptic, IAU 2006
    ob = pymObl06(date1, date2)

    # Precession-bias matrix, IAU 2006
    bp = pymPmat06(date1, date2)

    # Equatorial of date to ecliptic rotation (about x-axis)
    e = np.eye(3)
    e = pymRx(ob, e)
    

    # ICRS to ecliptic: matrix multiplication
    rm = e @ bp

    return rm



#2025-10-01
def pymEpv00(date1, date2):
    """
    Earth position and velocity, heliocentric and barycentric, with
    respect to the Barycentric Celestial Reference System.

    Parameters:
    -----------
    date1, date2 : float
        TDB date as Julian Date

    Returns:
    --------
    pvh : numpy.ndarray
        (2, 3) array of heliocentric Earth position/velocity [au, au/d]
    pvb : numpy.ndarray  
        (2, 3) array of barycentric Earth position/velocity [au, au/d]

    Raises:
    -------
    ValueError
        If date outside the range 1900-2100 AD
    """

    # Check date range
    t = ((date1 - DJ00) + date2) / DJY
    if abs(t) > 100.0:
        raise ValueError("date outside the range 1900-2100 AD")

    # Matrix elements for orienting the analytical model to DE405
    ori_mat = np.array([
        [ 1.0,             0.000000211284, -0.000000091603],
        [-0.000000230286,  0.917482137087, -0.397776982902],
        [ 0.0,             0.397776982902,  0.917482137087]
    ])
    
    
    # Optimized harmonic function for position and velocity computation
    def pv_harmonic_terms(coeff_sets, t, t2):
        """
        Compute position and velocity using harmonic terms with optimized vectorization
        
        Parameters:
        -----------
        coeff_sets : list of ndarray
            List of coefficient arrays for T^0, T^1, T^2 terms
        t : float
            Time in Julian years
        t2 : float
            Time squared
            
        Returns:
        --------
        position : float
            Position component
        velocity_d : float
            Velocity derivative component (before division by DJY)
        """
        xyz_t  = 0.0
        xyzd_t = 0.0
        
       
        # Process each coefficient set (T^0, T^1, T^2)
        for power, coeffs in enumerate(coeff_sets):
            if coeffs.size == 0:
                continue
                
            a = coeffs[:, 0]  # amplitudes
            b = coeffs[:, 1]  # phases  
            c = coeffs[:, 2]  # frequencies
            
            # Vectorized computation of trigonometric terms
            ct    = c * t
            phase = b + ct
            cos_p = np.cos(phase)
            sin_p = np.sin(phase)
            
            # Common terms for all power levels
            a_cos   = a * cos_p
            a_c_sin = a * c * sin_p
            
            if power == 0:    # T^0 terms
                xyz_t  += np.sum(a_cos)
                xyzd_t -= np.sum(a_c_sin)
                
            elif power == 1:  # T^1 terms
                xyz_t  += np.sum(a_cos * t)
                xyzd_t += np.sum(a_cos - a_c_sin * t)
                
            elif power == 2:  # T^2 terms
                xyz_t  += np.sum(a_cos * t2)
                xyzd_t += np.sum(2.0 * a_cos * t - a_c_sin * t2)
        
        return xyz_t, xyzd_t

    # Sun-to-Earth, T^0, X
    e0x = np.array([
        [0.9998292878132e+00, 0.1753485171504e+01, 0.6283075850446e+01],
        [0.8352579567414e-02, 0.1710344404582e+01, 0.1256615170089e+02],
        [0.5611445335148e-02, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.1046664295572e-03, 0.1667225416770e+01, 0.1884922755134e+02],
        [0.3110842534677e-04, 0.6687513390251e+00, 0.8399684731857e+02],
        [0.2552413503550e-04, 0.5830637358413e+00, 0.5296909721118e+00],
        [0.2137207845781e-04, 0.1092330954011e+01, 0.1577343543434e+01],
        [0.1680240182951e-04, 0.4955366134987e+00, 0.6279552690824e+01],
        [0.1679012370795e-04, 0.6153014091901e+01, 0.6286599010068e+01],
        [0.1445526946777e-04, 0.3472744100492e+01, 0.2352866153506e+01],
                                                                        
        [0.1091038246184e-04, 0.3689845786119e+01, 0.5223693906222e+01],
        [0.9344399733932e-05, 0.6073934645672e+01, 0.1203646072878e+02],
        [0.8993182910652e-05, 0.3175705249069e+01, 0.1021328554739e+02],
        [0.5665546034116e-05, 0.2152484672246e+01, 0.1059381944224e+01],
        [0.6844146703035e-05, 0.1306964099750e+01, 0.5753384878334e+01],
        [0.7346610905565e-05, 0.4354980070466e+01, 0.3981490189893e+00],
        [0.6815396474414e-05, 0.2218229211267e+01, 0.4705732307012e+01],
        [0.6112787253053e-05, 0.5384788425458e+01, 0.6812766822558e+01],
        [0.4518120711239e-05, 0.6087604012291e+01, 0.5884926831456e+01],
        [0.4521963430706e-05, 0.1279424524906e+01, 0.6256777527156e+01],
                                                                        
        [0.4497426764085e-05, 0.5369129144266e+01, 0.6309374173736e+01],
        [0.4062190566959e-05, 0.5436473303367e+00, 0.6681224869435e+01],
        [0.5412193480192e-05, 0.7867838528395e+00, 0.7755226100720e+00],
        [0.5469839049386e-05, 0.1461440311134e+01, 0.1414349524433e+02],
        [0.5205264083477e-05, 0.4432944696116e+01, 0.7860419393880e+01],
        [0.2149759935455e-05, 0.4502237496846e+01, 0.1150676975667e+02],
        [0.2279109618501e-05, 0.1239441308815e+01, 0.7058598460518e+01],
        [0.2259282939683e-05, 0.3272430985331e+01, 0.4694002934110e+01],
        [0.2558950271319e-05, 0.2265471086404e+01, 0.1216800268190e+02],
        [0.2561581447555e-05, 0.1454740653245e+01, 0.7099330490126e+00],
                                                                        
        [0.1781441115440e-05, 0.2962068630206e+01, 0.7962980379786e+00],
        [0.1612005874644e-05, 0.1473255041006e+01, 0.5486777812467e+01],
        [0.1818630667105e-05, 0.3743903293447e+00, 0.6283008715021e+01],
        [0.1818601377529e-05, 0.6274174354554e+01, 0.6283142985870e+01],
        [0.1554475925257e-05, 0.1624110906816e+01, 0.2513230340178e+02],
        [0.2090948029241e-05, 0.5852052276256e+01, 0.1179062909082e+02],
        [0.2000176345460e-05, 0.4072093298513e+01, 0.1778984560711e+02],
        [0.1289535917759e-05, 0.5217019331069e+01, 0.7079373888424e+01],
        [0.1281135307881e-05, 0.4802054538934e+01, 0.3738761453707e+01],
        [0.1518229005692e-05, 0.8691914742502e+00, 0.2132990797783e+00],
                                                                        
        [0.9450128579027e-06, 0.4601859529950e+01, 0.1097707878456e+02],
        [0.7781119494996e-06, 0.1844352816694e+01, 0.8827390247185e+01],
        [0.7733407759912e-06, 0.3582790154750e+01, 0.5507553240374e+01],
        [0.7350644318120e-06, 0.2695277788230e+01, 0.1589072916335e+01],
        [0.6535928827023e-06, 0.3651327986142e+01, 0.1176985366291e+02],
        [0.6324624183656e-06, 0.2241302375862e+01, 0.6262300422539e+01],
        [0.6298565300557e-06, 0.4407122406081e+01, 0.6303851278352e+01],
        [0.8587037089179e-06, 0.3024307223119e+01, 0.1672837615881e+03],
        [0.8299954491035e-06, 0.6192539428237e+01, 0.3340612434717e+01],
        [0.6311263503401e-06, 0.2014758795416e+01, 0.7113454667900e-02],
                                                                        
        [0.6005646745452e-06, 0.3399500503397e+01, 0.4136910472696e+01],
        [0.7917715109929e-06, 0.2493386877837e+01, 0.6069776770667e+01],
        [0.7556958099685e-06, 0.4159491740143e+01, 0.6496374930224e+01],
        [0.6773228244949e-06, 0.4034162934230e+01, 0.9437762937313e+01],
        [0.5370708577847e-06, 0.1562219163734e+01, 0.1194447056968e+01],
        [0.5710804266203e-06, 0.2662730803386e+01, 0.6282095334605e+01],
        [0.5709824583726e-06, 0.3985828430833e+01, 0.6284056366286e+01],
        [0.5143950896447e-06, 0.1308144688689e+01, 0.6290189305114e+01],
        [0.5088010604546e-06, 0.5352817214804e+01, 0.6275962395778e+01],
        [0.4960369085172e-06, 0.2644267922349e+01, 0.6127655567643e+01],
                                                                        
        [0.4803137891183e-06, 0.4008844192080e+01, 0.6438496133249e+01],
        [0.5731747768225e-06, 0.3794550174597e+01, 0.3154687086868e+01],
        [0.4735947960579e-06, 0.6107118308982e+01, 0.3128388763578e+01],
        [0.4808348796625e-06, 0.4771458618163e+01, 0.8018209333619e+00],
        [0.4115073743137e-06, 0.3327111335159e+01, 0.8429241228195e+01],
        [0.5230575889287e-06, 0.5305708551694e+01, 0.1336797263425e+02],
        [0.5133977889215e-06, 0.5784230738814e+01, 0.1235285262111e+02],
        [0.5065815825327e-06, 0.2052064793679e+01, 0.1185621865188e+02],
        [0.4339831593868e-06, 0.3644994195830e+01, 0.1726015463500e+02],
        [0.3952928638953e-06, 0.4930376436758e+01, 0.5481254917084e+01],
                                                                        
        [0.4898498111942e-06, 0.4542084219731e+00, 0.9225539266174e+01],
        [0.4757490209328e-06, 0.3161126388878e+01, 0.5856477690889e+01],
        [0.4727701669749e-06, 0.6214993845446e+00, 0.2544314396739e+01],
        [0.3800966681863e-06, 0.3040132339297e+01, 0.4265981595566e+00],
        [0.3257301077939e-06, 0.8064977360087e+00, 0.3930209696940e+01],
        [0.3255810528674e-06, 0.1974147981034e+01, 0.2146165377750e+01],
        [0.3252029748187e-06, 0.2845924913135e+01, 0.4164311961999e+01],
        [0.3255505635308e-06, 0.3017900824120e+01, 0.5088628793478e+01],
        [0.2801345211990e-06, 0.6109717793179e+01, 0.1256967486051e+02],
        [0.3688987740970e-06, 0.2911550235289e+01, 0.1807370494127e+02],
                                                                        
        [0.2475153429458e-06, 0.2179146025856e+01, 0.2629832328990e-01],
        [0.3033457749150e-06, 0.1994161050744e+01, 0.4535059491685e+01],
        [0.2186743763110e-06, 0.5125687237936e+01, 0.1137170464392e+02],
        [0.2764777032774e-06, 0.4822646860252e+00, 0.1256262854127e+02],
        [0.2199028768592e-06, 0.4637633293831e+01, 0.1255903824622e+02],
        [0.2046482824760e-06, 0.1467038733093e+01, 0.7084896783808e+01],
        [0.2611209147507e-06, 0.3044718783485e+00, 0.7143069561767e+02],
        [0.2286079656818e-06, 0.4764220356805e+01, 0.8031092209206e+01],
        [0.1855071202587e-06, 0.3383637774428e+01, 0.1748016358760e+01],
        [0.2324669506784e-06, 0.6189088449251e+01, 0.1831953657923e+02],
                                                                        
        [0.1709528015688e-06, 0.5874966729774e+00, 0.4933208510675e+01],
        [0.2168156875828e-06, 0.4302994009132e+01, 0.1044738781244e+02],
        [0.2106675556535e-06, 0.3800475419891e+01, 0.7477522907414e+01],
        [0.1430213830465e-06, 0.1294660846502e+01, 0.2942463415728e+01],
        [0.1388396901944e-06, 0.4594797202114e+01, 0.8635942003952e+01],
        [0.1922258844190e-06, 0.4943044543591e+00, 0.1729818233119e+02],
        [0.1888460058292e-06, 0.2426943912028e+01, 0.1561374759853e+03],
        [0.1789449386107e-06, 0.1582973303499e+00, 0.1592596075957e+01],
        [0.1360803685374e-06, 0.5197240440504e+01, 0.1309584267300e+02],
        [0.1504038014709e-06, 0.3120360916217e+01, 0.1649636139783e+02],
                                                                        
        [0.1382769533389e-06, 0.6164702888205e+01, 0.7632943190217e+01],
        [0.1438059769079e-06, 0.1437423770979e+01, 0.2042657109477e+02],
        [0.1326303260037e-06, 0.3609688799679e+01, 0.1213955354133e+02],
        [0.1159244950540e-06, 0.5463018167225e+01, 0.5331357529664e+01],
        [0.1433118149136e-06, 0.6028909912097e+01, 0.7342457794669e+01],
        [0.1234623148594e-06, 0.3109645574997e+01, 0.6279485555400e+01],
        [0.1233949875344e-06, 0.3539359332866e+01, 0.6286666145492e+01],
        [0.9927196061299e-07, 0.1259321569772e+01, 0.7234794171227e+01],
        [0.1242302191316e-06, 0.1065949392609e+01, 0.1511046609763e+02],
        [0.1098402195201e-06, 0.2192508743837e+01, 0.1098880815746e+02],
                                                                        
        [0.1158191395315e-06, 0.4054411278650e+01, 0.5729506548653e+01],
        [0.9048475596241e-07, 0.5429764748518e+01, 0.9623688285163e+01],
        [0.8889853269023e-07, 0.5046586206575e+01, 0.6148010737701e+01],
        [0.1048694242164e-06, 0.2628858030806e+01, 0.6836645152238e+01],
        [0.1112308378646e-06, 0.4177292719907e+01, 0.1572083878776e+02],
        [0.8631729709901e-07, 0.1601345232557e+01, 0.6418140963190e+01],
        [0.8527816951664e-07, 0.2463888997513e+01, 0.1471231707864e+02],
        [0.7892139456991e-07, 0.3154022088718e+01, 0.2118763888447e+01],
        [0.1051782905236e-06, 0.4795035816088e+01, 0.1349867339771e+01],
        [0.1048219943164e-06, 0.2952983395230e+01, 0.5999216516294e+01],
                                                                        
        [0.7435760775143e-07, 0.5420547991464e+01, 0.6040347114260e+01],
        [0.9869574106949e-07, 0.3695646753667e+01, 0.6566935184597e+01],
        [0.9156886364226e-07, 0.3922675306609e+01, 0.5643178611111e+01],
        [0.7006834356188e-07, 0.1233968624861e+01, 0.6525804586632e+01],
        [0.9806170182601e-07, 0.1919542280684e+01, 0.2122839202813e+02],
        [0.9052289673607e-07, 0.4615902724369e+01, 0.4690479774488e+01],
        [0.7554200867893e-07, 0.1236863719072e+01, 0.1253985337760e+02],
        [0.8215741286498e-07, 0.3286800101559e+00, 0.1097355562493e+02],
        [0.7185178575397e-07, 0.5880942158367e+01, 0.6245048154254e+01],
        [0.7130726476180e-07, 0.7674871987661e+00, 0.6321103546637e+01],
                                                                        
        [0.6650894461162e-07, 0.6987129150116e+00, 0.5327476111629e+01],
        [0.7396888823688e-07, 0.3576824794443e+01, 0.5368044267797e+00],
        [0.7420588884775e-07, 0.5033615245369e+01, 0.2354323048545e+02],
        [0.6141181642908e-07, 0.9449927045673e+00, 0.1296430071988e+02],
        [0.6373557924058e-07, 0.6206342280341e+01, 0.9517183207817e+00],
        [0.6359474329261e-07, 0.5036079095757e+01, 0.1990745094947e+01],
        [0.5740173582646e-07, 0.6105106371350e+01, 0.9555997388169e+00],
        [0.7019864084602e-07, 0.7237747359018e+00, 0.5225775174439e+00],
        [0.6398054487042e-07, 0.3976367969666e+01, 0.2407292145756e+02],
        [0.7797092650498e-07, 0.4305423910623e+01, 0.2200391463820e+02],
                                                                        
        [0.6466760000900e-07, 0.3500136825200e+01, 0.5230807360890e+01],
        [0.7529417043890e-07, 0.3514779246100e+01, 0.1842262939178e+02],
        [0.6924571140892e-07, 0.2743457928679e+01, 0.1554202828031e+00],
        [0.6220798650222e-07, 0.2242598118209e+01, 0.1845107853235e+02],
        [0.5870209391853e-07, 0.2332832707527e+01, 0.6398972393349e+00],
        [0.6263953473888e-07, 0.2191105358956e+01, 0.6277552955062e+01],
        [0.6257781390012e-07, 0.4457559396698e+01, 0.6288598745829e+01],
        [0.5697304945123e-07, 0.3499234761404e+01, 0.1551045220144e+01],
        [0.6335438746791e-07, 0.6441691079251e+00, 0.5216580451554e+01],
        [0.6377258441152e-07, 0.2252599151092e+01, 0.5650292065779e+01],
                                                                        
        [0.6484841818165e-07, 0.1992812417646e+01, 0.1030928125552e+00],
        [0.4735551485250e-07, 0.3744672082942e+01, 0.1431416805965e+02],
        [0.4628595996170e-07, 0.1334226211745e+01, 0.5535693017924e+00],
        [0.6258152336933e-07, 0.4395836159154e+01, 0.2608790314060e+02],
        [0.6196171366594e-07, 0.2587043007997e+01, 0.8467247584405e+02],
        [0.6159556952126e-07, 0.4782499769128e+01, 0.2394243902548e+03],
        [0.4987741172394e-07, 0.7312257619924e+00, 0.7771377146812e+02],
        [0.5459280703142e-07, 0.3001376372532e+01, 0.6179983037890e+01],
        [0.4863461189999e-07, 0.3767222128541e+01, 0.9027992316901e+02],
        [0.5349912093158e-07, 0.3663594450273e+01, 0.6386168663001e+01],
                                                                        
        [0.5673725607806e-07, 0.4331187919049e+01, 0.6915859635113e+01],
        [0.4745485060512e-07, 0.5816195745518e+01, 0.6282970628506e+01],
        [0.4745379005326e-07, 0.8323672435672e+00, 0.6283181072386e+01],
        [0.4049002796321e-07, 0.3785023976293e+01, 0.6254626709878e+01],
        [0.4247084014515e-07, 0.2378220728783e+01, 0.7875671926403e+01],
        [0.4026912363055e-07, 0.2864103423269e+01, 0.6311524991013e+01],
        [0.4062935011774e-07, 0.2415408595975e+01, 0.3634620989887e+01],
        [0.5347771048509e-07, 0.3343479309801e+01, 0.2515860172507e+02],
        [0.4829494136505e-07, 0.2821742398262e+01, 0.5760498333002e+01],
        [0.4342554404599e-07, 0.5624662458712e+01, 0.7238675589263e+01],
                                                                        
        [0.4021599184361e-07, 0.5557250275009e+00, 0.1101510648075e+02],
        [0.4104900474558e-07, 0.3296691780005e+01, 0.6709674010002e+01],
        [0.4376532905131e-07, 0.3814443999443e+01, 0.6805653367890e+01],
        [0.3314590480650e-07, 0.3560229189250e+01, 0.1259245002418e+02],
        [0.3232421839643e-07, 0.5185389180568e+01, 0.1066495398892e+01],
        [0.3541176318876e-07, 0.3921381909679e+01, 0.9917696840332e+01],
        [0.3689831242681e-07, 0.4190658955386e+01, 0.1192625446156e+02],
        [0.3890605376774e-07, 0.5546023371097e+01, 0.7478166569050e-01],
        [0.3038559339780e-07, 0.6231032794494e+01, 0.1256621883632e+02],
        [0.3137083969782e-07, 0.6207063419190e+01, 0.4292330755499e+01],
                                                                        
        [0.4024004081854e-07, 0.1195257375713e+01, 0.1334167431096e+02],
        [0.3300234879283e-07, 0.1804694240998e+01, 0.1057540660594e+02],
        [0.3635399155575e-07, 0.5597811343500e+01, 0.6208294184755e+01],
        [0.3032668691356e-07, 0.3191059366530e+01, 0.1805292951336e+02],
        [0.2809652069058e-07, 0.4094348032570e+01, 0.3523159621801e-02],
        [0.3696955383823e-07, 0.5219282738794e+01, 0.5966683958112e+01],
        [0.3562894142503e-07, 0.1037247544554e+01, 0.6357857516136e+01],
        [0.3510598524148e-07, 0.1430020816116e+01, 0.6599467742779e+01],
        [0.3617736142953e-07, 0.3002911403677e+01, 0.6019991944201e+01],
        [0.2624524910730e-07, 0.2437046757292e+01, 0.6702560555334e+01],
                                                                        
        [0.2535824204490e-07, 0.1581594689647e+01, 0.3141537925223e+02],
        [0.3519787226257e-07, 0.5379863121521e+01, 0.2505706758577e+03],
        [0.2578406709982e-07, 0.4904222639329e+01, 0.1673046366289e+02],
        [0.3423887981473e-07, 0.3646448997315e+01, 0.6546159756691e+01],
        [0.2776083886467e-07, 0.3307829300144e+01, 0.1272157198369e+02],
        [0.3379592818379e-07, 0.1747541251125e+01, 0.1494531617769e+02],
        [0.3050255426284e-07, 0.1784689432607e-01, 0.4732030630302e+01],
        [0.2652378350236e-07, 0.4420055276260e+01, 0.5863591145557e+01],
        [0.2374498173768e-07, 0.3629773929208e+01, 0.2388894113936e+01],
        [0.2716451255140e-07, 0.3079623706780e+01, 0.1202934727411e+02],
                                                                        
        [0.3038583699229e-07, 0.3312487903507e+00, 0.1256608456547e+02],
        [0.2220681228760e-07, 0.5265520401774e+01, 0.1336244973887e+02],
        [0.3044156540912e-07, 0.4766664081250e+01, 0.2908881142201e+02],
        [0.2731859923561e-07, 0.5069146530691e+01, 0.1391601904066e+02],
        [0.2285603018171e-07, 0.5954935112271e+01, 0.6076890225335e+01],
        [0.2025006454555e-07, 0.4061789589267e+01, 0.4701116388778e+01],
        [0.2012597519804e-07, 0.2485047705241e+01, 0.6262720680387e+01],
        [0.2003406962258e-07, 0.4163779209320e+01, 0.6303431020504e+01],
        [0.2207863441371e-07, 0.6923839133828e+00, 0.6489261475556e+01],
        [0.2481374305624e-07, 0.5944173595676e+01, 0.1204357418345e+02],
                                                                        
        [0.2130923288870e-07, 0.4641013671967e+01, 0.5746271423666e+01],
        [0.2446370543391e-07, 0.6125796518757e+01, 0.1495633313810e+00],
        [0.1932492759052e-07, 0.2234572324504e+00, 0.1352175143971e+02],
        [0.2600122568049e-07, 0.4281012405440e+01, 0.4590910121555e+01],
        [0.2431754047488e-07, 0.1429943874870e+00, 0.1162474756779e+01],
        [0.1875902869209e-07, 0.9781803816948e+00, 0.6279194432410e+01],
        [0.1874381139426e-07, 0.5670368130173e+01, 0.6286957268481e+01],
        [0.2156696047173e-07, 0.2008985006833e+01, 0.1813929450232e+02],
        [0.1965076182484e-07, 0.2566186202453e+00, 0.4686889479442e+01],
        [0.2334816372359e-07, 0.4408121891493e+01, 0.1002183730415e+02],
                                                                        
        [0.1869937408802e-07, 0.5272745038656e+01, 0.2427287361862e+00],
        [0.2436236460883e-07, 0.4407720479029e+01, 0.9514313292143e+02],
        [0.1761365216611e-07, 0.1943892315074e+00, 0.1351787002167e+02],
        [0.2156289480503e-07, 0.1418570924545e+01, 0.6037244212485e+01],
        [0.2164748979255e-07, 0.4724603439430e+01, 0.2301353951334e+02],
        [0.2222286670853e-07, 0.2400266874598e+01, 0.1266924451345e+02],
        [0.2070901414929e-07, 0.5230348028732e+01, 0.6528907488406e+01],
        [0.1792745177020e-07, 0.2099190328945e+01, 0.6819880277225e+01],
        [0.1841802068445e-07, 0.3467527844848e+00, 0.6514761976723e+02],
        [0.1578401631718e-07, 0.7098642356340e+00, 0.2077542790660e-01],
                                                                        
        [0.1561690152531e-07, 0.5943349620372e+01, 0.6272439236156e+01],
        [0.1558591045463e-07, 0.7040653478980e+00, 0.6293712464735e+01],
        [0.1737356469576e-07, 0.4487064760345e+01, 0.1765478049437e+02],
        [0.1434755619991e-07, 0.2993391570995e+01, 0.1102062672231e+00],
        [0.1482187806654e-07, 0.2278049198251e+01, 0.1052268489556e+01],
        [0.1424812827089e-07, 0.1682114725827e+01, 0.1311972100268e+02],
        [0.1380282448623e-07, 0.3262668602579e+01, 0.1017725758696e+02],
        [0.1811481244566e-07, 0.3187771221777e+01, 0.1887552587463e+02],
        [0.1504446185696e-07, 0.5650162308647e+01, 0.7626583626240e-01],
        [0.1740776154137e-07, 0.5487068607507e+01, 0.1965104848470e+02],
                                                                        
        [0.1374339536251e-07, 0.5745688172201e+01, 0.6016468784579e+01],
        [0.1761377477704e-07, 0.5748060203659e+01, 0.2593412433514e+02],
        [0.1535138225795e-07, 0.6226848505790e+01, 0.9411464614024e+01],
        [0.1788140543676e-07, 0.6189318878563e+01, 0.3301902111895e+02],
        [0.1375002807996e-07, 0.5371812884394e+01, 0.6327837846670e+00],
        [0.1242115758632e-07, 0.1471687569712e+01, 0.3894181736510e+01],
        [0.1450977333938e-07, 0.4143836662127e+01, 0.1277945078067e+02],
        [0.1297579575023e-07, 0.9003477661957e+00, 0.6549682916313e+01],
        [0.1462667934821e-07, 0.5760505536428e+01, 0.1863592847156e+02],
        [0.1381774374799e-07, 0.1085471729463e+01, 0.2379164476796e+01],
                                                                        
        [0.1682333169307e-07, 0.5409870870133e+01, 0.1620077269078e+02],
        [0.1190812918837e-07, 0.1397205174601e+01, 0.1149965630200e+02],
        [0.1221434762106e-07, 0.9001804809095e+00, 0.1257326515556e+02],
        [0.1549934644860e-07, 0.4262528275544e+01, 0.1820933031200e+02],
        [0.1252138953050e-07, 0.1411642012027e+01, 0.6993008899458e+01],
        [0.1237078905387e-07, 0.2844472403615e+01, 0.2435678079171e+02],
        [0.1446953389615e-07, 0.5295835522223e+01, 0.3813291813120e-01],
        [0.1388446457170e-07, 0.4969428135497e+01, 0.2458316379602e+00],
        [0.1019339179228e-07, 0.2491369561806e+01, 0.6112403035119e+01],
        [0.1258880815343e-07, 0.4679426248976e+01, 0.5429879531333e+01],
                                                                        
        [0.1297768238261e-07, 0.1074509953328e+01, 0.1249137003520e+02],
        [0.9913505718094e-08, 0.4735097918224e+01, 0.6247047890016e+01],
        [0.9830453155969e-08, 0.4158649187338e+01, 0.6453748665772e+01],
        [0.1192615865309e-07, 0.3438208613699e+01, 0.6290122169689e+01],
        [0.9835874798277e-08, 0.1913300781229e+01, 0.6319103810876e+01],
        [0.9639087569277e-08, 0.9487683644125e+00, 0.8273820945392e+01],
        [0.1175716107001e-07, 0.3228141664287e+01, 0.6276029531202e+01],
        [0.1018926508678e-07, 0.2216607854300e+01, 0.1254537627298e+02],
        [0.9500087869225e-08, 0.2625116459733e+01, 0.1256517118505e+02],
        [0.9664192916575e-08, 0.5860562449214e+01, 0.6259197520765e+01],
                                                                        
        [0.9612858712203e-08, 0.7885682917381e+00, 0.6306954180126e+01],
        [0.1117645675413e-07, 0.3932148831189e+01, 0.1779695906178e+02],
        [0.1158864052160e-07, 0.9995605521691e+00, 0.1778273215245e+02],
        [0.9021043467028e-08, 0.5263769742673e+01, 0.6172869583223e+01],
        [0.8836134773563e-08, 0.1496843220365e+01, 0.1692165728891e+01],
        [0.1045872200691e-07, 0.7009039517214e+00, 0.2204125344462e+00],
        [0.1211463487798e-07, 0.4041544938511e+01, 0.8257698122054e+02],
        [0.8541990804094e-08, 0.1447586692316e+01, 0.6393282117669e+01],
        [0.1038720703636e-07, 0.4594249718112e+00, 0.1550861511662e+02],
        [0.1126722351445e-07, 0.3925550579036e+01, 0.2061856251104e+00],
                                                                        
        [0.8697373859631e-08, 0.4411341856037e+01, 0.9491756770005e+00],
        [0.8869380028441e-08, 0.2402659724813e+01, 0.3903911373650e+01],
        [0.9247014693258e-08, 0.1401579743423e+01, 0.6267823317922e+01],
        [0.9205062930950e-08, 0.5245978000814e+01, 0.6298328382969e+01],
        [0.8000745038049e-08, 0.3590803356945e+01, 0.2648454860559e+01],
        [0.9168973650819e-08, 0.2470150501679e+01, 0.1498544001348e+03],
        [0.1075444949238e-07, 0.1328606161230e+01, 0.3694923081589e+02],
        [0.7817298525817e-08, 0.6162256225998e+01, 0.4804209201333e+01],
        [0.9541469226356e-08, 0.3942568967039e+01, 0.1256713221673e+02],
        [0.9821910122027e-08, 0.2360246287233e+00, 0.1140367694411e+02],
                                                                        
        [0.9897822023777e-08, 0.4619805634280e+01, 0.2280573557157e+02],
        [0.7737289283765e-08, 0.3784727847451e+01, 0.7834121070590e+01],
        [0.9260204034710e-08, 0.2223352487601e+01, 0.2787043132925e+01],
        [0.7320252888486e-08, 0.1288694636874e+01, 0.6282655592598e+01],
        [0.7319785780946e-08, 0.5359869567774e+01, 0.6283496108294e+01],
        [0.7147219933778e-08, 0.5516616675856e+01, 0.1725663147538e+02],
        [0.7946502829878e-08, 0.2630459984567e+01, 0.1241073141809e+02],
        [0.9001711808932e-08, 0.2849815827227e+01, 0.6281591679874e+01],
        [0.8994041507257e-08, 0.3795244450750e+01, 0.6284560021018e+01],
        [0.8298582787358e-08, 0.5236413127363e+00, 0.1241658836951e+02],
                                                                        
        [0.8526596520710e-08, 0.4794605424426e+01, 0.1098419223922e+02],
        [0.8209822103197e-08, 0.1578752370328e+01, 0.1096996532989e+02],
        [0.6357049861094e-08, 0.5708926113761e+01, 0.1596186371003e+01],
        [0.7370473179049e-08, 0.3842402530241e+01, 0.4061219149443e+01],
        [0.7232154664726e-08, 0.3067548981535e+01, 0.1610006857377e+03],
        [0.6328765494903e-08, 0.1313930030069e+01, 0.1193336791622e+02],
        [0.8030064908595e-08, 0.3488500408886e+01, 0.8460828644453e+00],
        [0.6275464259232e-08, 0.1532061626198e+01, 0.8531963191132e+00],
        [0.7051897446325e-08, 0.3285859929993e+01, 0.5849364236221e+01],
        [0.6161593705428e-08, 0.1477341999464e+01, 0.5573142801433e+01],
                                                                        
        [0.7754683957278e-08, 0.1586118663096e+01, 0.8662240327241e+01],
        [0.5889928990701e-08, 0.1304887868803e+01, 0.1232342296471e+02],
        [0.5705756047075e-08, 0.4555333589350e+01, 0.1258692712880e+02],
        [0.5964178808332e-08, 0.3001762842062e+01, 0.5333900173445e+01],
        [0.6712446027467e-08, 0.4886780007595e+01, 0.1171295538178e+02],
        [0.5941809275464e-08, 0.4701509603824e+01, 0.9779108567966e+01],
        [0.5466993627395e-08, 0.4588357817278e+01, 0.1884211409667e+02],
        [0.6340512090980e-08, 0.1164543038893e+01, 0.5217580628120e+02],
        [0.6325505710045e-08, 0.3919171259645e+01, 0.1041998632314e+02],
        [0.6164789509685e-08, 0.2143828253542e+01, 0.6151533897323e+01],
                                                                        
        [0.5263330812430e-08, 0.6066564434241e+01, 0.1885275071096e+02],
        [0.5597087780221e-08, 0.2926316429472e+01, 0.4337116142245e+00],
        [0.5396556236817e-08, 0.3244303591505e+01, 0.6286362197481e+01],
        [0.5396615148223e-08, 0.3404304703662e+01, 0.6279789503410e+01],
        [0.7091832443341e-08, 0.8532377803192e+00, 0.4907302013889e+01],
        [0.6572352589782e-08, 0.4901966774419e+01, 0.1176433076753e+02],
        [0.5960236060795e-08, 0.1874672315797e+01, 0.1422690933580e-01],
        [0.5125480043511e-08, 0.3735726064334e+01, 0.1245594543367e+02],
        [0.5928241866410e-08, 0.4502033899935e+01, 0.6414617803568e+01],
        [0.5249600357424e-08, 0.4372334799878e+01, 0.1151388321134e+02],
                                                                        
        [0.6059171276087e-08, 0.2581617302908e+01, 0.6062663316000e+01],
        [0.5295235081662e-08, 0.2974811513158e+01, 0.3496032717521e+01],
        [0.5820561875933e-08, 0.1796073748244e+00, 0.2838593341516e+00],
        [0.4754696606440e-08, 0.1981998136973e+01, 0.3104930017775e+01],
        [0.6385053548955e-08, 0.2559174171605e+00, 0.6133512519065e+01],
        [0.6589828273941e-08, 0.2750967106776e+01, 0.4087944051283e+02],
        [0.5383376567189e-08, 0.6325947523578e+00, 0.2248384854122e+02],
        [0.5928941683538e-08, 0.1672304519067e+01, 0.1581959461667e+01],
        [0.4816060709794e-08, 0.3512566172575e+01, 0.9388005868221e+01],
        [0.6003381586512e-08, 0.5610932219189e+01, 0.5326786718777e+01],
                                                                        
        [0.5504225393105e-08, 0.4037501131256e+01, 0.6503488384892e+01],
        [0.5353772620129e-08, 0.6122774968240e+01, 0.1735668374386e+03],
        [0.5786253768544e-08, 0.5527984999515e+01, 0.1350651127443e+00],
        [0.5065706702002e-08, 0.9980765573624e+00, 0.1248988586463e+02],
        [0.5972838885276e-08, 0.6044489493203e+01, 0.2673594526851e+02],
        [0.5323585877961e-08, 0.3924265998147e+01, 0.4171425416666e+01],
        [0.5210772682858e-08, 0.6220111376901e+01, 0.2460261242967e+02],
        [0.4726549040535e-08, 0.3716043206862e+01, 0.7232251527446e+01],
        [0.6029425105059e-08, 0.8548704071116e+00, 0.3227113045244e+03],
        [0.4481542826513e-08, 0.1426925072829e+01, 0.5547199253223e+01],
                                                                        
        [0.5836024505068e-08, 0.7135651752625e-01, 0.7285056171570e+02],
        [0.4137046613272e-08, 0.5330767643283e+01, 0.1087398597200e+02],
        [0.5171977473924e-08, 0.4494262335353e+00, 0.1884570439172e+02],
        [0.5694429833732e-08, 0.2952369582215e+01, 0.9723862754494e+02],
        [0.4009158925298e-08, 0.3500003416535e+01, 0.6244942932314e+01],
        [0.4784939596873e-08, 0.6196709413181e+01, 0.2929661536378e+02],
        [0.3983725022610e-08, 0.5103690031897e+01, 0.4274518229222e+01],
        [0.3870535232462e-08, 0.3187569587401e+01, 0.6321208768577e+01],
        [0.5140501213951e-08, 0.1668924357457e+01, 0.1232032006293e+02],
        [0.3849034819355e-08, 0.4445722510309e+01, 0.1726726808967e+02],
                                                                        
        [0.4002383075060e-08, 0.5226224152423e+01, 0.7018952447668e+01],
        [0.3890719543549e-08, 0.4371166550274e+01, 0.1491901785440e+02],
        [0.4887084607881e-08, 0.5973556689693e+01, 0.1478866649112e+01],
        [0.3739939287592e-08, 0.2089084714600e+01, 0.6922973089781e+01],
        [0.5031925918209e-08, 0.4658371936827e+01, 0.1715706182245e+02],
        [0.4387748764954e-08, 0.4825580552819e+01, 0.2331413144044e+03],
        [0.4147398098865e-08, 0.3739003524998e+01, 0.1376059875786e+02],
        [0.3719089993586e-08, 0.1148941386536e+01, 0.6297302759782e+01],
        [0.3934238461056e-08, 0.1559893008343e+01, 0.7872148766781e+01],
        [0.3672471375622e-08, 0.5516145383612e+01, 0.6268848941110e+01],
                                                                        
        [0.3768911277583e-08, 0.6116053700563e+01, 0.4157198507331e+01],
        [0.4033388417295e-08, 0.5076821746017e+01, 0.1567108171867e+02],
        [0.3764194617832e-08, 0.8164676232075e+00, 0.3185192151914e+01],
        [0.4840628226284e-08, 0.1360479453671e+01, 0.1252801878276e+02],
        [0.4949443923785e-08, 0.2725622229926e+01, 0.1617106187867e+03],
        [0.4117393089971e-08, 0.6054459628492e+00, 0.5642198095270e+01],
        [0.3925754020428e-08, 0.8570462135210e+00, 0.2139354194808e+02],
        [0.3630551757923e-08, 0.3552067338279e+01, 0.6294805223347e+01],
        [0.3627274802357e-08, 0.3096565085313e+01, 0.6271346477544e+01],
        [0.3806143885093e-08, 0.6367751709777e+00, 0.1725304118033e+02],
                                                                        
        [0.4433254641565e-08, 0.4848461503937e+01, 0.7445550607224e+01],
        [0.3712319846576e-08, 0.1331950643655e+01, 0.4194847048887e+00],
        [0.3849847534783e-08, 0.4958368297746e+00, 0.9562891316684e+00],
        [0.3483955430165e-08, 0.2237215515707e+01, 0.1161697602389e+02],
        [0.3961912730982e-08, 0.3332402188575e+01, 0.2277943724828e+02],
        [0.3419978244481e-08, 0.5785600576016e+01, 0.1362553364512e+02],
        [0.3329417758177e-08, 0.9812676559709e-01, 0.1685848245639e+02],
        [0.4207206893193e-08, 0.9494780468236e+00, 0.2986433403208e+02],
        [0.3268548976410e-08, 0.1739332095686e+00, 0.5749861718712e+01],
        [0.3321880082685e-08, 0.1423354800666e+01, 0.6279143387820e+01],
                                                                        
        [0.4503173010852e-08, 0.2314972675293e+00, 0.1385561574497e+01],
        [0.4316599090954e-08, 0.1012646782616e+00, 0.4176041334900e+01],
        [0.3283493323850e-08, 0.5233306881265e+01, 0.6287008313071e+01],
        [0.3164033542343e-08, 0.4005597257511e+01, 0.2099539292909e+02],
        [0.4159720956725e-08, 0.5365676242020e+01, 0.5905702259363e+01],
        [0.3565176892217e-08, 0.4284440620612e+01, 0.3932462625300e-02],
        [0.3514440950221e-08, 0.4270562636575e+01, 0.7335344340001e+01],
        [0.3540596871909e-08, 0.5953553201060e+01, 0.1234573916645e+02],
        [0.2960769905118e-08, 0.1115180417718e+01, 0.2670964694522e+02],
        [0.2962213739684e-08, 0.3863811918186e+01, 0.6408777551755e+00],
                                                                        
        [0.3883556700251e-08, 0.1268617928302e+01, 0.6660449441528e+01],
        [0.2919225516346e-08, 0.4908605223265e+01, 0.1375773836557e+01],
        [0.3115158863370e-08, 0.3744519976885e+01, 0.3802769619140e-01],
        [0.4099438144212e-08, 0.4173244670532e+01, 0.4480965020977e+02],
        [0.2899531858964e-08, 0.5910601428850e+01, 0.2059724391010e+02],
        [0.3289733429855e-08, 0.2488050078239e+01, 0.1081813534213e+02],
        [0.3933075612875e-08, 0.1122363652883e+01, 0.3773735910827e+00],
        [0.3021403764467e-08, 0.4951973724904e+01, 0.2982630633589e+02],
        [0.2798598949757e-08, 0.5117057845513e+01, 0.1937891852345e+02],
        [0.3397421302707e-08, 0.6104159180476e+01, 0.6923953605621e+01],
                                                                        
        [0.3720398002179e-08, 0.1184933429829e+01, 0.3066615496545e+02],
        [0.3598484186267e-08, 0.3505282086105e+01, 0.6147450479709e+01],
        [0.3694594027310e-08, 0.2286651088141e+01, 0.2636725487657e+01],
        [0.2680444152969e-08, 0.1871816775482e+00, 0.6816289982179e+01],
        [0.3497574865641e-08, 0.3143251755431e+01, 0.6418701221183e+01],
        [0.3130274129494e-08, 0.2462167316018e+01, 0.1235996607578e+02],
        [0.3241119069551e-08, 0.4256374004686e+01, 0.1652265972112e+02],
        [0.2601960842061e-08, 0.4970362941425e+01, 0.1045450126711e+02],
        [0.2690601527504e-08, 0.2372657824898e+01, 0.3163918923335e+00],
        [0.2908688152664e-08, 0.4232652627721e+01, 0.2828699048865e+02],
                                                                        
        [0.3120456131875e-08, 0.3925747001137e+00, 0.2195415756911e+02],
        [0.3148855423384e-08, 0.3093478330445e+01, 0.1172006883645e+02],
        [0.3051044261017e-08, 0.5560948248212e+01, 0.6055599646783e+01],
        [0.2826006876660e-08, 0.5072790310072e+01, 0.5120601093667e+01],
        [0.3100034191711e-08, 0.4998530231096e+01, 0.1799603123222e+02],
        [0.2398771640101e-08, 0.2561739802176e+01, 0.6255674361143e+01],
        [0.2384002842728e-08, 0.4087420284111e+01, 0.6310477339748e+01],
        [0.2842146517568e-08, 0.2515048217955e+01, 0.5469525544182e+01],
        [0.2847674371340e-08, 0.5235326497443e+01, 0.1034429499989e+02],
        [0.2903722140764e-08, 0.1088200795797e+01, 0.6510552054109e+01],
                                                                        
        [0.3187610710605e-08, 0.4710624424816e+01, 0.1693792562116e+03],
        [0.3048869992813e-08, 0.2857975896445e+00, 0.8390110365991e+01],
        [0.2860216950984e-08, 0.2241619020815e+01, 0.2243449970715e+00],
        [0.2701117683113e-08, 0.6651573305272e-01, 0.6129297044991e+01],
        [0.2509891590152e-08, 0.1285135324585e+01, 0.1044027435778e+02],
        [0.2623200252223e-08, 0.2981229834530e+00, 0.6436854655901e+01],
        [0.2622541669202e-08, 0.6122470726189e+01, 0.9380959548977e+01],
        [0.2818435667099e-08, 0.4251087148947e+01, 0.5934151399930e+01],
        [0.2365196797465e-08, 0.3465070460790e+01, 0.2470570524223e+02],
        [0.2358704646143e-08, 0.5791603815350e+01, 0.8671969964381e+01],
                                                                        
        [0.2388299481390e-08, 0.4142483772941e+01, 0.7096626156709e+01],
        [0.1996041217224e-08, 0.2101901889496e+01, 0.1727188400790e+02],
        [0.2687593060336e-08, 0.1526689456959e+01, 0.7075506709219e+02],
        [0.2618913670810e-08, 0.2397684236095e+01, 0.6632000300961e+01],
        [0.2571523050364e-08, 0.5751929456787e+00, 0.6206810014183e+01],
        [0.2582135006946e-08, 0.5595464352926e+01, 0.4873985990671e+02],
        [0.2372530190361e-08, 0.5092689490655e+01, 0.1590676413561e+02],
        [0.2357178484712e-08, 0.4444363527851e+01, 0.3097883698531e+01],
        [0.2451590394723e-08, 0.3108251687661e+01, 0.6612329252343e+00],
        [0.2370045949608e-08, 0.2608133861079e+01, 0.3459636466239e+02],
                                                                        
        [0.2268997267358e-08, 0.3639717753384e+01, 0.2844914056730e-01],
        [0.1731432137906e-08, 0.1741898445707e+00, 0.2019909489111e+02],
        [0.1629869741622e-08, 0.3902225646724e+01, 0.3035599730800e+02],
        [0.2206215801974e-08, 0.4971131250731e+01, 0.6281667977667e+01],
        [0.2205469554680e-08, 0.1677462357110e+01, 0.6284483723224e+01],
        [0.2148792362509e-08, 0.4236259604006e+01, 0.1980482729015e+02],
        [0.1873733657847e-08, 0.5926814998687e+01, 0.2876692439167e+02],
        [0.2026573758959e-08, 0.4349643351962e+01, 0.2449240616245e+02],
        [0.1807770325110e-08, 0.5700940482701e+01, 0.2045286941806e+02],
        [0.1881174408581e-08, 0.6601286363430e+00, 0.2358125818164e+02],
                                                                        
        [0.1368023671690e-08, 0.2211098592752e+01, 0.2473415438279e+02],
        [0.1720017916280e-08, 0.4942488551129e+01, 0.1679593901136e+03],
        [0.1702427665131e-08, 0.1452233856386e+01, 0.3338575901272e+03],
        [0.1414032510054e-08, 0.5525357721439e+01, 0.1624205518357e+03],
        [0.1652626045364e-08, 0.4108794283624e+01, 0.8956999012000e+02],
        [0.1642957769686e-08, 0.7344335209984e+00, 0.5267006960365e+02],
        [0.1614952403624e-08, 0.3541213951363e+01, 0.3332657872986e+02],
        [0.1535988291188e-08, 0.4031094072151e+01, 0.3852657435933e+02],
        [0.1593193738177e-08, 0.4185136203609e+01, 0.2282781046519e+03],
        [0.1074569126382e-08, 0.1720485636868e+01, 0.8397383534231e+02],
                                                                        
        [0.1074408214509e-08, 0.2758613420318e+01, 0.8401985929482e+02],
        [0.9700199670465e-09, 0.4216686842097e+01, 0.7826370942180e+02],
        [0.1258433517061e-08, 0.2575068876639e+00, 0.3115650189215e+03],
        [0.1240303229539e-08, 0.4800844956756e+00, 0.1784300471910e+03],
        [0.9018345948127e-09, 0.3896756361552e+00, 0.5886454391678e+02],
        [0.1135301432805e-08, 0.3700805023550e+00, 0.7842370451713e+02],
        [0.9215887951370e-09, 0.4364579276638e+01, 0.1014262087719e+03],
        [0.1055401054147e-08, 0.2156564222111e+01, 0.5660027930059e+02],
        [0.1008725979831e-08, 0.5454015785234e+01, 0.4245678405627e+02],
        [0.7217398104321e-09, 0.1597772562175e+01, 0.2457074661053e+03],
                                                                        
        [0.6912033134447e-09, 0.5824090621461e+01, 0.1679936946371e+03],
        [0.6833881523549e-09, 0.3578778482835e+01, 0.6053048899753e+02],
        [0.4887304205142e-09, 0.3724362812423e+01, 0.9656299901946e+02],
        [0.5173709754788e-09, 0.5422427507933e+01, 0.2442876000072e+03],
        [0.4671353097145e-09, 0.2396106924439e+01, 0.1435713242844e+03],
        [0.5652608439480e-09, 0.2804028838685e+01, 0.8365903305582e+02],
        [0.5604061331253e-09, 0.1638816006247e+01, 0.8433466158131e+02],
        [0.4712723365400e-09, 0.8979003224474e+00, 0.3164282286739e+03],
        [0.4909967465112e-09, 0.3210426725516e+01, 0.4059982187939e+03],
        [0.4771358267658e-09, 0.5308027211629e+01, 0.1805255418145e+03],
                                                                        
        [0.3943451445989e-09, 0.2195145341074e+01, 0.2568537517081e+03],
        [0.3952109120244e-09, 0.5081189491586e+01, 0.2449975330562e+03],
        [0.3788134594789e-09, 0.4345171264441e+01, 0.1568131045107e+03],
        [0.3738330190479e-09, 0.2613062847997e+01, 0.3948519331910e+03],
        [0.3099866678136e-09, 0.2846760817689e+01, 0.1547176098872e+03],
        [0.2002962716768e-09, 0.4921360989412e+01, 0.2268582385539e+03],
        [0.2198291338754e-09, 0.1130360117454e+00, 0.1658638954901e+03],
        [0.1491958330784e-09, 0.4228195232278e+01, 0.2219950288015e+03],
        [0.1475384076173e-09, 0.3005721811604e+00, 0.3052819430710e+03],
        [0.1661626624624e-09, 0.7830125621203e+00, 0.2526661704812e+03],
                                                                        
        [0.9015823460025e-10, 0.3807792942715e+01, 0.4171445043968e+03] 
        
    ])

    # Sun-to-Earth, T^0, Y
    e0y = np.array([
        [0.9998921098898e+00, 0.1826583913846e+00, 0.6283075850446e+01],  
        [-0.2442700893735e-01,0.0000000000000e+00, 0.0000000000000e+00],  
        [0.8352929742915e-02, 0.1395277998680e+00, 0.1256615170089e+02],  
        [0.1046697300177e-03, 0.9641423109763e-01, 0.1884922755134e+02],  
        [0.3110841876663e-04, 0.5381140401712e+01, 0.8399684731857e+02],  
        [0.2570269094593e-04, 0.5301016407128e+01, 0.5296909721118e+00],  
        [0.2147389623610e-04, 0.2662510869850e+01, 0.1577343543434e+01],  
        [0.1680344384050e-04, 0.5207904119704e+01, 0.6279552690824e+01],  
        [0.1679117312193e-04, 0.4582187486968e+01, 0.6286599010068e+01],  
        [0.1440512068440e-04, 0.1900688517726e+01, 0.2352866153506e+01],  
          
        [0.1135139664999e-04, 0.5273108538556e+01, 0.5223693906222e+01],  
        [0.9345482571018e-05, 0.4503047687738e+01, 0.1203646072878e+02],  
        [0.9007418719568e-05, 0.1605621059637e+01, 0.1021328554739e+02],  
        [0.5671536712314e-05, 0.5812849070861e+00, 0.1059381944224e+01],  
        [0.7451401861666e-05, 0.2807346794836e+01, 0.3981490189893e+00],  
        [0.6393470057114e-05, 0.6029224133855e+01, 0.5753384878334e+01],  
        [0.6814275881697e-05, 0.6472990145974e+00, 0.4705732307012e+01],  
        [0.6113705628887e-05, 0.3813843419700e+01, 0.6812766822558e+01],  
        [0.4503851367273e-05, 0.4527804370996e+01, 0.5884926831456e+01],  
        [0.4522249141926e-05, 0.5991783029224e+01, 0.6256777527156e+01],  
          
        [0.4501794307018e-05, 0.3798703844397e+01, 0.6309374173736e+01],  
        [0.5514927480180e-05, 0.3961257833388e+01, 0.5507553240374e+01],  
        [0.4062862799995e-05, 0.5256247296369e+01, 0.6681224869435e+01],  
        [0.5414900429712e-05, 0.5499032014097e+01, 0.7755226100720e+00],  
        [0.5463153987424e-05, 0.6173092454097e+01, 0.1414349524433e+02],  
        [0.5071611859329e-05, 0.2870244247651e+01, 0.7860419393880e+01],  
        [0.2195112094455e-05, 0.2952338617201e+01, 0.1150676975667e+02],  
        [0.2279139233919e-05, 0.5951775132933e+01, 0.7058598460518e+01],  
        [0.2278386100876e-05, 0.4845456398785e+01, 0.4694002934110e+01],  
        [0.2559088003308e-05, 0.6945321117311e+00, 0.1216800268190e+02],  
          
        [0.2561079286856e-05, 0.6167224608301e+01, 0.7099330490126e+00],  
        [0.1792755796387e-05, 0.1400122509632e+01, 0.7962980379786e+00],  
        [0.1818715656502e-05, 0.4703347611830e+01, 0.6283142985870e+01],  
        [0.1818744924791e-05, 0.5086748900237e+01, 0.6283008715021e+01],  
        [0.1554518791390e-05, 0.5331008042713e-01, 0.2513230340178e+02],  
        [0.2063265737239e-05, 0.4283680484178e+01, 0.1179062909082e+02],  
        [0.1497613520041e-05, 0.6074207826073e+01, 0.5486777812467e+01],  
        [0.2000617940427e-05, 0.2501426281450e+01, 0.1778984560711e+02],  
        [0.1289731195580e-05, 0.3646340599536e+01, 0.7079373888424e+01],  
        [0.1282657998934e-05, 0.3232864804902e+01, 0.3738761453707e+01],  
          
        [0.1528915968658e-05, 0.5581433416669e+01, 0.2132990797783e+00],  
        [0.1187304098432e-05, 0.5453576453694e+01, 0.9437762937313e+01],  
        [0.7842782928118e-06, 0.2823953922273e+00, 0.8827390247185e+01],  
        [0.7352892280868e-06, 0.1124369580175e+01, 0.1589072916335e+01],  
        [0.6570189360797e-06, 0.2089154042840e+01, 0.1176985366291e+02],  
        [0.6324967590410e-06, 0.6704855581230e+00, 0.6262300422539e+01],  
        [0.6298289872283e-06, 0.2836414855840e+01, 0.6303851278352e+01],  
        [0.6476686465855e-06, 0.4852433866467e+00, 0.7113454667900e-02],  
        [0.8587034651234e-06, 0.1453511005668e+01, 0.1672837615881e+03],  
        [0.8068948788113e-06, 0.9224087798609e+00, 0.6069776770667e+01],  
          
        [0.8353786011661e-06, 0.4631707184895e+01, 0.3340612434717e+01],  
        [0.6009324532132e-06, 0.1829498827726e+01, 0.4136910472696e+01],  
        [0.7558158559566e-06, 0.2588596800317e+01, 0.6496374930224e+01],  
        [0.5809279504503e-06, 0.5516818853476e+00, 0.1097707878456e+02],  
        [0.5374131950254e-06, 0.6275674734960e+01, 0.1194447056968e+01],  
        [0.5711160507326e-06, 0.1091905956872e+01, 0.6282095334605e+01],  
        [0.5710183170746e-06, 0.2415001635090e+01, 0.6284056366286e+01],  
        [0.5144373590610e-06, 0.6020336443438e+01, 0.6290189305114e+01],  
        [0.5103108927267e-06, 0.3775634564605e+01, 0.6275962395778e+01],  
        [0.4960654697891e-06, 0.1073450946756e+01, 0.6127655567643e+01],  
          
        [0.4786385689280e-06, 0.2431178012310e+01, 0.6438496133249e+01],  
        [0.6109911263665e-06, 0.5343356157914e+01, 0.3154687086868e+01],  
        [0.4839898944024e-06, 0.5830833594047e-01, 0.8018209333619e+00],  
        [0.4734822623919e-06, 0.4536080134821e+01, 0.3128388763578e+01],  
        [0.4834741473290e-06, 0.2585090489754e+00, 0.7084896783808e+01],  
        [0.5134858581156e-06, 0.4213317172603e+01, 0.1235285262111e+02],  
        [0.5064004264978e-06, 0.4814418806478e+00, 0.1185621865188e+02],  
        [0.3753476772761e-06, 0.1599953399788e+01, 0.8429241228195e+01],  
        [0.4935264014283e-06, 0.2157417556873e+01, 0.2544314396739e+01],  
        [0.3950929600897e-06, 0.3359394184254e+01, 0.5481254917084e+01],  
          
        [0.4895849789777e-06, 0.5165704376558e+01, 0.9225539266174e+01],  
        [0.4215241688886e-06, 0.2065368800993e+01, 0.1726015463500e+02],  
        [0.3796773731132e-06, 0.1468606346612e+01, 0.4265981595566e+00],  
        [0.3114178142515e-06, 0.3615638079474e+01, 0.2146165377750e+01],  
        [0.3260664220838e-06, 0.4417134922435e+01, 0.4164311961999e+01],  
        [0.3976996123008e-06, 0.4700866883004e+01, 0.5856477690889e+01],  
        [0.2801459672924e-06, 0.4538902060922e+01, 0.1256967486051e+02],  
        [0.3638931868861e-06, 0.1334197991475e+01, 0.1807370494127e+02],  
        [0.2487013269476e-06, 0.3749275558275e+01, 0.2629832328990e-01],  
        [0.3034165481994e-06, 0.4236622030873e+00, 0.4535059491685e+01],  
          
        [0.2676278825586e-06, 0.5970848007811e+01, 0.3930209696940e+01],  
        [0.2764903818918e-06, 0.5194636754501e+01, 0.1256262854127e+02],  
        [0.2485149930507e-06, 0.1002434207846e+01, 0.5088628793478e+01],  
        [0.2199305540941e-06, 0.3066773098403e+01, 0.1255903824622e+02],  
        [0.2571106500435e-06, 0.7588312459063e+00, 0.1336797263425e+02],  
        [0.2049751817158e-06, 0.3444977434856e+01, 0.1137170464392e+02],  
        [0.2599707296297e-06, 0.1873128542205e+01, 0.7143069561767e+02],  
        [0.1785018072217e-06, 0.5015891306615e+01, 0.1748016358760e+01],  
        [0.2324833891115e-06, 0.4618271239730e+01, 0.1831953657923e+02],  
        [0.1709711119545e-06, 0.5300003455669e+01, 0.4933208510675e+01],  
          
        [0.2107159351716e-06, 0.2229819815115e+01, 0.7477522907414e+01],  
        [0.1750333080295e-06, 0.6161485880008e+01, 0.1044738781244e+02],  
        [0.2000598210339e-06, 0.2967357299999e+01, 0.8031092209206e+01],  
        [0.1380920248681e-06, 0.3027007923917e+01, 0.8635942003952e+01],  
        [0.1412460470299e-06, 0.6037597163798e+01, 0.2942463415728e+01],  
        [0.1888459803001e-06, 0.8561476243374e+00, 0.1561374759853e+03],  
        [0.1788370542585e-06, 0.4869736290209e+01, 0.1592596075957e+01],  
        [0.1360893296167e-06, 0.3626411886436e+01, 0.1309584267300e+02],  
        [0.1506846530160e-06, 0.1550975377427e+01, 0.1649636139783e+02],  
        [0.1800913376176e-06, 0.2075826033190e+01, 0.1729818233119e+02],  
          
        [0.1436261390649e-06, 0.6148876420255e+01, 0.2042657109477e+02],  
        [0.1220227114151e-06, 0.4382583879906e+01, 0.7632943190217e+01],  
        [0.1337883603592e-06, 0.2036644327361e+01, 0.1213955354133e+02],  
        [0.1159326650738e-06, 0.3892276994687e+01, 0.5331357529664e+01],  
        [0.1352853128569e-06, 0.1447950649744e+01, 0.1673046366289e+02],  
        [0.1433408296083e-06, 0.4457854692961e+01, 0.7342457794669e+01],  
        [0.1234701666518e-06, 0.1538818147151e+01, 0.6279485555400e+01],  
        [0.1234027192007e-06, 0.1968523220760e+01, 0.6286666145492e+01],  
        [0.1244024091797e-06, 0.5779803499985e+01, 0.1511046609763e+02],  
        [0.1097934945516e-06, 0.6210975221388e+00, 0.1098880815746e+02],  
          
        [0.1254611329856e-06, 0.2591963807998e+01, 0.1572083878776e+02],  
        [0.1158247286784e-06, 0.2483612812670e+01, 0.5729506548653e+01],  
        [0.9039078252960e-07, 0.3857554579796e+01, 0.9623688285163e+01],  
        [0.9108024978836e-07, 0.5826368512984e+01, 0.7234794171227e+01],  
        [0.8887068108436e-07, 0.3475694573987e+01, 0.6148010737701e+01],  
        [0.8632374035438e-07, 0.3059070488983e-01, 0.6418140963190e+01],  
        [0.7893186992967e-07, 0.1583194837728e+01, 0.2118763888447e+01],  
        [0.8297650201172e-07, 0.8519770534637e+00, 0.1471231707864e+02],  
        [0.1019759578988e-06, 0.1319598738732e+00, 0.1349867339771e+01],  
        [0.1010037696236e-06, 0.9937860115618e+00, 0.6836645152238e+01],  
          
        [0.1047727548266e-06, 0.1382138405399e+01, 0.5999216516294e+01],  
        [0.7351993881086e-07, 0.3833397851735e+01, 0.6040347114260e+01],  
        [0.9868771092341e-07, 0.2124913814390e+01, 0.6566935184597e+01],  
        [0.7007321959390e-07, 0.5946305343763e+01, 0.6525804586632e+01],  
        [0.6861411679709e-07, 0.4574654977089e+01, 0.7238675589263e+01],  
        [0.7554519809614e-07, 0.5949232686844e+01, 0.1253985337760e+02],  
        [0.9541880448335e-07, 0.3495242990564e+01, 0.2122839202813e+02],  
        [0.7185606722155e-07, 0.4310113471661e+01, 0.6245048154254e+01],  
        [0.7131360871710e-07, 0.5480309323650e+01, 0.6321103546637e+01],  
        [0.6651142021039e-07, 0.5411097713654e+01, 0.5327476111629e+01],  
          
        [0.8538618213667e-07, 0.1827849973951e+01, 0.1101510648075e+02],  
        [0.8634954288044e-07, 0.5443584943349e+01, 0.5643178611111e+01],  
        [0.7449415051484e-07, 0.2011535459060e+01, 0.5368044267797e+00],  
        [0.7421047599169e-07, 0.3464562529249e+01, 0.2354323048545e+02],  
        [0.6140694354424e-07, 0.5657556228815e+01, 0.1296430071988e+02],  
        [0.6353525143033e-07, 0.3463816593821e+01, 0.1990745094947e+01],  
        [0.6221964013447e-07, 0.1532259498697e+01, 0.9517183207817e+00],  
        [0.5852480257244e-07, 0.1375396598875e+01, 0.9555997388169e+00],  
        [0.6398637498911e-07, 0.2405645801972e+01, 0.2407292145756e+02],  
        [0.7039744069878e-07, 0.5397541799027e+01, 0.5225775174439e+00],  
          
        [0.6977997694382e-07, 0.4762347105419e+01, 0.1097355562493e+02],  
        [0.7460629558396e-07, 0.2711944692164e+01, 0.2200391463820e+02],  
        [0.5376577536101e-07, 0.2352980430239e+01, 0.1431416805965e+02],  
        [0.7530607893556e-07, 0.1943940180699e+01, 0.1842262939178e+02],  
        [0.6822928971605e-07, 0.4337651846959e+01, 0.1554202828031e+00],  
        [0.6220772380094e-07, 0.6716871369278e+00, 0.1845107853235e+02],  
        [0.6586950799043e-07, 0.2229714460505e+01, 0.5216580451554e+01],  
        [0.5873800565771e-07, 0.7627013920580e+00, 0.6398972393349e+00],  
        [0.6264346929745e-07, 0.6202785478961e+00, 0.6277552955062e+01],  
        [0.6257929115669e-07, 0.2886775596668e+01, 0.6288598745829e+01],  
          
        [0.5343536033409e-07, 0.1977241012051e+01, 0.4690479774488e+01],  
        [0.5587849781714e-07, 0.1922923484825e+01, 0.1551045220144e+01],  
        [0.6905100845603e-07, 0.3570757164631e+01, 0.1030928125552e+00],  
        [0.6178957066649e-07, 0.5197558947765e+01, 0.5230807360890e+01],  
        [0.6187270224331e-07, 0.8193497368922e+00, 0.5650292065779e+01],  
        [0.5385664291426e-07, 0.5406336665586e+01, 0.7771377146812e+02],  
        [0.6329363917926e-07, 0.2837760654536e+01, 0.2608790314060e+02],  
        [0.4546018761604e-07, 0.2933580297050e+01, 0.5535693017924e+00],  
        [0.6196091049375e-07, 0.4157871494377e+01, 0.8467247584405e+02],  
        [0.6159555108218e-07, 0.3211703561703e+01, 0.2394243902548e+03],  
          
        [0.4995340539317e-07, 0.1459098102922e+01, 0.4732030630302e+01],  
        [0.5457031243572e-07, 0.1430457676136e+01, 0.6179983037890e+01],  
        [0.4863461418397e-07, 0.2196425916730e+01, 0.9027992316901e+02],  
        [0.5342947626870e-07, 0.2086612890268e+01, 0.6386168663001e+01],  
        [0.5674296648439e-07, 0.2760204966535e+01, 0.6915859635113e+01],  
        [0.4745783120161e-07, 0.4245368971862e+01, 0.6282970628506e+01],  
        [0.4745676961198e-07, 0.5544725787016e+01, 0.6283181072386e+01],  
        [0.4049796869973e-07, 0.2213984363586e+01, 0.6254626709878e+01],  
        [0.4248333596940e-07, 0.8075781952896e+00, 0.7875671926403e+01],  
        [0.4027178070205e-07, 0.1293268540378e+01, 0.6311524991013e+01],  
          
        [0.4066543943476e-07, 0.3986141175804e+01, 0.3634620989887e+01],  
        [0.4858863787880e-07, 0.1276112738231e+01, 0.5760498333002e+01],  
        [0.5277398263530e-07, 0.4916111741527e+01, 0.2515860172507e+02],  
        [0.4105635656559e-07, 0.1725805864426e+01, 0.6709674010002e+01],  
        [0.4376781925772e-07, 0.2243642442106e+01, 0.6805653367890e+01],  
        [0.3235827894693e-07, 0.3614135118271e+01, 0.1066495398892e+01],  
        [0.3073244740308e-07, 0.2460873393460e+01, 0.5863591145557e+01],  
        [0.3088609271373e-07, 0.5678431771790e+01, 0.9917696840332e+01],  
        [0.3393022279836e-07, 0.3814017477291e+01, 0.1391601904066e+02],  
        [0.3038686508802e-07, 0.4660216229171e+01, 0.1256621883632e+02],  
          
        [0.4019677752497e-07, 0.5906906243735e+01, 0.1334167431096e+02],  
        [0.3288834998232e-07, 0.9536146445882e+00, 0.1620077269078e+02],  
        [0.3889973794631e-07, 0.3942205097644e+01, 0.7478166569050e-01],  
        [0.3050438987141e-07, 0.1624810271286e+01, 0.1805292951336e+02],  
        [0.3601142564638e-07, 0.4030467142575e+01, 0.6208294184755e+01],  
        [0.3689015557141e-07, 0.3648878818694e+01, 0.5966683958112e+01],  
        [0.3563471893565e-07, 0.5749584017096e+01, 0.6357857516136e+01],  
        [0.2776183170667e-07, 0.2630124187070e+01, 0.3523159621801e-02],  
        [0.2922350530341e-07, 0.1790346403629e+01, 0.1272157198369e+02],  
        [0.3511076917302e-07, 0.6142198301611e+01, 0.6599467742779e+01],  
          
        [0.3619351007632e-07, 0.1432421386492e+01, 0.6019991944201e+01],  
        [0.2561254711098e-07, 0.2302822475792e+01, 0.1259245002418e+02],  
        [0.2626903942920e-07, 0.8660470994571e+00, 0.6702560555334e+01],  
        [0.2550187397083e-07, 0.6069721995383e+01, 0.1057540660594e+02],  
        [0.2535873526138e-07, 0.1079020331795e-01, 0.3141537925223e+02],  
        [0.3519786153847e-07, 0.3809066902283e+01, 0.2505706758577e+03],  
        [0.3424651492873e-07, 0.2075435114417e+01, 0.6546159756691e+01],  
        [0.2372676630861e-07, 0.2057803120154e+01, 0.2388894113936e+01],  
        [0.2710980779541e-07, 0.1510068488010e+01, 0.1202934727411e+02],  
        [0.3038710889704e-07, 0.5043617528901e+01, 0.1256608456547e+02],  
          
        [0.2220364130585e-07, 0.3694793218205e+01, 0.1336244973887e+02],  
        [0.3025880825460e-07, 0.5450618999049e-01, 0.2908881142201e+02],  
        [0.2784493486864e-07, 0.3381164084502e+01, 0.1494531617769e+02],  
        [0.2294414142438e-07, 0.4382309025210e+01, 0.6076890225335e+01],  
        [0.2012723294724e-07, 0.9142212256518e+00, 0.6262720680387e+01],  
        [0.2036357831958e-07, 0.5676172293154e+01, 0.4701116388778e+01],  
        [0.2003474823288e-07, 0.2592767977625e+01, 0.6303431020504e+01],  
        [0.2207144900109e-07, 0.5404976271180e+01, 0.6489261475556e+01],  
        [0.2481664905135e-07, 0.4373284587027e+01, 0.1204357418345e+02],  
        [0.2674949182295e-07, 0.5859182188482e+01, 0.4590910121555e+01],  
          
        [0.2450554720322e-07, 0.4555381557451e+01, 0.1495633313810e+00],  
        [0.2601975986457e-07, 0.3933165584959e+01, 0.1965104848470e+02],  
        [0.2199860022848e-07, 0.5227977189087e+01, 0.1351787002167e+02],  
        [0.2448121172316e-07, 0.4858060353949e+01, 0.1162474756779e+01],  
        [0.1876014864049e-07, 0.5690546553605e+01, 0.6279194432410e+01],  
        [0.1874513219396e-07, 0.4099539297446e+01, 0.6286957268481e+01],  
        [0.2156380842559e-07, 0.4382594769913e+00, 0.1813929450232e+02],  
        [0.1981691240061e-07, 0.1829784152444e+01, 0.4686889479442e+01],  
        [0.2329992648539e-07, 0.2836254278973e+01, 0.1002183730415e+02],  
        [0.1765184135302e-07, 0.2803494925833e+01, 0.4292330755499e+01],  
          
        [0.2436368366085e-07, 0.2836897959677e+01, 0.9514313292143e+02],  
        [0.2164089203889e-07, 0.6127522446024e+01, 0.6037244212485e+01],  
        [0.1847755034221e-07, 0.3683163635008e+01, 0.2427287361862e+00],  
        [0.1674798769966e-07, 0.3316993867246e+00, 0.1311972100268e+02],  
        [0.2222542124356e-07, 0.8294097805480e+00, 0.1266924451345e+02],  
        [0.2071074505925e-07, 0.3659492220261e+01, 0.6528907488406e+01],  
        [0.1608224471835e-07, 0.4774492067182e+01, 0.1352175143971e+02],  
        [0.1857583439071e-07, 0.2873120597682e+01, 0.8662240327241e+01],  
        [0.1793018836159e-07, 0.5282441177929e+00, 0.6819880277225e+01],  
        [0.1575391221692e-07, 0.1320789654258e+01, 0.1102062672231e+00],  
          
        [0.1840132009557e-07, 0.1917110916256e+01, 0.6514761976723e+02],  
        [0.1760917288281e-07, 0.2972635937132e+01, 0.5746271423666e+01],  
        [0.1561779518516e-07, 0.4372569261981e+01, 0.6272439236156e+01],  
        [0.1558687885205e-07, 0.5416424926425e+01, 0.6293712464735e+01],  
        [0.1951359382579e-07, 0.3094448898752e+01, 0.2301353951334e+02],  
        [0.1569144275614e-07, 0.2802103689808e+01, 0.1765478049437e+02],  
        [0.1479130389462e-07, 0.2136435020467e+01, 0.2077542790660e-01],  
        [0.1467828510764e-07, 0.7072627435674e+00, 0.1052268489556e+01],  
        [0.1627627337440e-07, 0.3947607143237e+01, 0.6327837846670e+00],  
        [0.1503498479758e-07, 0.4079248909190e+01, 0.7626583626240e-01],  
          
        [0.1297967708237e-07, 0.6269637122840e+01, 0.1149965630200e+02],  
        [0.1374416896634e-07, 0.4175657970702e+01, 0.6016468784579e+01],  
        [0.1783812325219e-07, 0.1476540547560e+01, 0.3301902111895e+02],  
        [0.1525884228756e-07, 0.4653477715241e+01, 0.9411464614024e+01],  
        [0.1451067396763e-07, 0.2573001128225e+01, 0.1277945078067e+02],  
        [0.1297713111950e-07, 0.5612799618771e+01, 0.6549682916313e+01],  
        [0.1462784012820e-07, 0.4189661623870e+01, 0.1863592847156e+02],  
        [0.1384185980007e-07, 0.2656915472196e+01, 0.2379164476796e+01],  
        [0.1221497599801e-07, 0.5612515760138e+01, 0.1257326515556e+02],  
        [0.1560574525896e-07, 0.4783414317919e+01, 0.1887552587463e+02],  
          
        [0.1544598372036e-07, 0.2694431138063e+01, 0.1820933031200e+02],  
        [0.1531678928696e-07, 0.4105103489666e+01, 0.2593412433514e+02],  
        [0.1349321503795e-07, 0.3082437194015e+00, 0.5120601093667e+01],  
        [0.1252030290917e-07, 0.6124072334087e+01, 0.6993008899458e+01],  
        [0.1459243816687e-07, 0.3733103981697e+01, 0.3813291813120e-01],  
        [0.1226103625262e-07, 0.1267127706817e+01, 0.2435678079171e+02],  
        [0.1019449641504e-07, 0.4367790112269e+01, 0.1725663147538e+02],  
        [0.1380789433607e-07, 0.3387201768700e+01, 0.2458316379602e+00],  
        [0.1019453421658e-07, 0.9204143073737e+00, 0.6112403035119e+01],  
        [0.1297929434405e-07, 0.5786874896426e+01, 0.1249137003520e+02],  
          
        [0.9912677786097e-08, 0.3164232870746e+01, 0.6247047890016e+01],  
        [0.9829386098599e-08, 0.2586762413351e+01, 0.6453748665772e+01],  
        [0.1226807746104e-07, 0.6239068436607e+01, 0.5429879531333e+01],  
        [0.1192691755997e-07, 0.1867380051424e+01, 0.6290122169689e+01],  
        [0.9836499227081e-08, 0.3424716293727e+00, 0.6319103810876e+01],  
        [0.9642862564285e-08, 0.5661372990657e+01, 0.8273820945392e+01],  
        [0.1165184404862e-07, 0.5768367239093e+01, 0.1778273215245e+02],  
        [0.1175794418818e-07, 0.1657351222943e+01, 0.6276029531202e+01],  
        [0.1018948635601e-07, 0.6458292350865e+00, 0.1254537627298e+02],  
        [0.9500383606676e-08, 0.1054306140741e+01, 0.1256517118505e+02],  
          
        [0.1227512202906e-07, 0.2505278379114e+01, 0.2248384854122e+02],  
        [0.9664792009993e-08, 0.4289737277000e+01, 0.6259197520765e+01],  
        [0.9613285666331e-08, 0.5500597673141e+01, 0.6306954180126e+01],  
        [0.1117906736211e-07, 0.2361405953468e+01, 0.1779695906178e+02],  
        [0.9611378640782e-08, 0.2851310576269e+01, 0.2061856251104e+00],  
        [0.8845354852370e-08, 0.6208777705343e+01, 0.1692165728891e+01],  
        [0.1054046966600e-07, 0.5413091423934e+01, 0.2204125344462e+00],  
        [0.1215539124483e-07, 0.5613969479755e+01, 0.8257698122054e+02],  
        [0.9932460955209e-08, 0.1106124877015e+01, 0.1017725758696e+02],  
        [0.8785804715043e-08, 0.2869224476477e+01, 0.9491756770005e+00],  
          
        [0.8538084097562e-08, 0.6159640899344e+01, 0.6393282117669e+01],  
        [0.8648994369529e-08, 0.1374901198784e+01, 0.4804209201333e+01],  
        [0.1039063219067e-07, 0.5171080641327e+01, 0.1550861511662e+02],  
        [0.8867983926439e-08, 0.8317320304902e+00, 0.3903911373650e+01],  
        [0.8327495955244e-08, 0.3605591969180e+01, 0.6172869583223e+01],  
        [0.9243088356133e-08, 0.6114299196843e+01, 0.6267823317922e+01],  
        [0.9205657357835e-08, 0.3675153683737e+01, 0.6298328382969e+01],  
        [0.1033269714606e-07, 0.3313328813024e+01, 0.5573142801433e+01],  
        [0.8001706275552e-08, 0.2019980960053e+01, 0.2648454860559e+01],  
        [0.9171858254191e-08, 0.8992015524177e+00, 0.1498544001348e+03],  
          
        [0.1075327150242e-07, 0.2898669963648e+01, 0.3694923081589e+02],  
        [0.9884866689828e-08, 0.4946715904478e+01, 0.1140367694411e+02],  
        [0.9541835576677e-08, 0.2371787888469e+01, 0.1256713221673e+02],  
        [0.7739903376237e-08, 0.2213775190612e+01, 0.7834121070590e+01],  
        [0.7311962684106e-08, 0.3429378787739e+01, 0.1192625446156e+02],  
        [0.9724904869624e-08, 0.6195878564404e+01, 0.2280573557157e+02],  
        [0.9251628983612e-08, 0.6511509527390e+00, 0.2787043132925e+01],  
        [0.7320763787842e-08, 0.6001083639421e+01, 0.6282655592598e+01],  
        [0.7320296650962e-08, 0.3789073265087e+01, 0.6283496108294e+01],  
        [0.7947032271039e-08, 0.1059659582204e+01, 0.1241073141809e+02],  
          
        [0.9005277053115e-08, 0.1280315624361e+01, 0.6281591679874e+01],  
        [0.8995601652048e-08, 0.2224439106766e+01, 0.6284560021018e+01],  
        [0.8288040568796e-08, 0.5234914433867e+01, 0.1241658836951e+02],  
        [0.6359381347255e-08, 0.4137989441490e+01, 0.1596186371003e+01],  
        [0.8699572228626e-08, 0.1758411009497e+01, 0.6133512519065e+01],  
        [0.6456797542736e-08, 0.5919285089994e+01, 0.1685848245639e+02],  
        [0.7424573475452e-08, 0.5414616938827e+01, 0.4061219149443e+01],  
        [0.7235671196168e-08, 0.1496516557134e+01, 0.1610006857377e+03],  
        [0.8104015182733e-08, 0.1919918242764e+01, 0.8460828644453e+00],  
        [0.8098576535937e-08, 0.3819615855458e+01, 0.3894181736510e+01],  
          
        [0.6275292346625e-08, 0.6244264115141e+01, 0.8531963191132e+00],  
        [0.6052432989112e-08, 0.5037731872610e+00, 0.1567108171867e+02],  
        [0.5705651535817e-08, 0.2984557271995e+01, 0.1258692712880e+02],  
        [0.5789650115138e-08, 0.6087038140697e+01, 0.1193336791622e+02],  
        [0.5512132153377e-08, 0.5855668994076e+01, 0.1232342296471e+02],  
        [0.7388890819102e-08, 0.2443128574740e+01, 0.4907302013889e+01],  
        [0.5467593991798e-08, 0.3017561234194e+01, 0.1884211409667e+02],  
        [0.6388519802999e-08, 0.5887386712935e+01, 0.5217580628120e+02],  
        [0.6106777149944e-08, 0.3483461059895e+00, 0.1422690933580e-01],  
        [0.7383420275489e-08, 0.5417387056707e+01, 0.2358125818164e+02],  
          
        [0.5505208141738e-08, 0.2848193644783e+01, 0.1151388321134e+02],  
        [0.6310757462877e-08, 0.2349882520828e+01, 0.1041998632314e+02],  
        [0.6166904929691e-08, 0.5728575944077e+00, 0.6151533897323e+01],  
        [0.5263442042754e-08, 0.4495796125937e+01, 0.1885275071096e+02],  
        [0.5591828082629e-08, 0.1355441967677e+01, 0.4337116142245e+00],  
        [0.5397051680497e-08, 0.1673422864307e+01, 0.6286362197481e+01],  
        [0.5396992745159e-08, 0.1833502206373e+01, 0.6279789503410e+01],  
        [0.6572913000726e-08, 0.3331122065824e+01, 0.1176433076753e+02],  
        [0.5123421866413e-08, 0.2165327142679e+01, 0.1245594543367e+02],  
        [0.5930495725999e-08, 0.2931146089284e+01, 0.6414617803568e+01],  
          
        [0.6431797403933e-08, 0.4134407994088e+01, 0.1350651127443e+00],  
        [0.5003182207604e-08, 0.3805420303749e+01, 0.1096996532989e+02],  
        [0.5587731032504e-08, 0.1082469260599e+01, 0.6062663316000e+01],  
        [0.5935263407816e-08, 0.8384333678401e+00, 0.5326786718777e+01],  
        [0.4756019827760e-08, 0.3552588749309e+01, 0.3104930017775e+01],  
        [0.6599951172637e-08, 0.4320826409528e+01, 0.4087944051283e+02],  
        [0.5902606868464e-08, 0.4811879454445e+01, 0.5849364236221e+01],  
        [0.5921147809031e-08, 0.9942628922396e-01, 0.1581959461667e+01],  
        [0.5505382581266e-08, 0.2466557607764e+01, 0.6503488384892e+01],  
        [0.5353771071862e-08, 0.4551978748683e+01, 0.1735668374386e+03],  
          
        [0.5063282210946e-08, 0.5710812312425e+01, 0.1248988586463e+02],  
        [0.5926120403383e-08, 0.1333998428358e+01, 0.2673594526851e+02],  
        [0.5211016176149e-08, 0.4649315360760e+01, 0.2460261242967e+02],  
        [0.5347075084894e-08, 0.5512754081205e+01, 0.4171425416666e+01],  
        [0.4872609773574e-08, 0.1308025299938e+01, 0.5333900173445e+01],  
        [0.4727711321420e-08, 0.2144908368062e+01, 0.7232251527446e+01],  
        [0.6029426018652e-08, 0.5567259412084e+01, 0.3227113045244e+03],  
        [0.4321485284369e-08, 0.5230667156451e+01, 0.9388005868221e+01],  
        [0.4476406760553e-08, 0.6134081115303e+01, 0.5547199253223e+01],  
        [0.5835268277420e-08, 0.4783808492071e+01, 0.7285056171570e+02],  
          
        [0.5172183602748e-08, 0.5161817911099e+01, 0.1884570439172e+02],  
        [0.5693571465184e-08, 0.1381646203111e+01, 0.9723862754494e+02],  
        [0.4060634965349e-08, 0.3876705259495e+00, 0.4274518229222e+01],  
        [0.3967398770473e-08, 0.5029491776223e+01, 0.3496032717521e+01],  
        [0.3943754005255e-08, 0.1923162955490e+01, 0.6244942932314e+01],  
        [0.4781323427824e-08, 0.4633332586423e+01, 0.2929661536378e+02],  
        [0.3871483781204e-08, 0.1616650009743e+01, 0.6321208768577e+01],  
        [0.5141741733997e-08, 0.9817316704659e-01, 0.1232032006293e+02],  
        [0.4002385978497e-08, 0.3656161212139e+01, 0.7018952447668e+01],  
        [0.4901092604097e-08, 0.4404098713092e+01, 0.1478866649112e+01],  
          
        [0.3740932630345e-08, 0.5181188732639e+00, 0.6922973089781e+01],  
        [0.4387283718538e-08, 0.3254859566869e+01, 0.2331413144044e+03],  
        [0.5019197802033e-08, 0.3086773224677e+01, 0.1715706182245e+02],  
        [0.3834931695175e-08, 0.2797882673542e+01, 0.1491901785440e+02],  
        [0.3760413942497e-08, 0.2892676280217e+01, 0.1726726808967e+02],  
        [0.3719717204628e-08, 0.5861046025739e+01, 0.6297302759782e+01],  
        [0.4145623530149e-08, 0.2168239627033e+01, 0.1376059875786e+02],  
        [0.3932788425380e-08, 0.6271811124181e+01, 0.7872148766781e+01],  
        [0.3686377476857e-08, 0.3936853151404e+01, 0.6268848941110e+01],  
        [0.3779077950339e-08, 0.1404148734043e+01, 0.4157198507331e+01],  
          
        [0.4091334550598e-08, 0.2452436180854e+01, 0.9779108567966e+01],  
        [0.3926694536146e-08, 0.6102292739040e+01, 0.1098419223922e+02],  
        [0.4841000253289e-08, 0.6072760457276e+01, 0.1252801878276e+02],  
        [0.4949340130240e-08, 0.1154832815171e+01, 0.1617106187867e+03],  
        [0.3761557737360e-08, 0.5527545321897e+01, 0.3185192151914e+01],  
        [0.3647396268188e-08, 0.1525035688629e+01, 0.6271346477544e+01],  
        [0.3932405074189e-08, 0.5570681040569e+01, 0.2139354194808e+02],  
        [0.3631322501141e-08, 0.1981240601160e+01, 0.6294805223347e+01],  
        [0.4130007425139e-08, 0.2050060880201e+01, 0.2195415756911e+02],  
        [0.4433905965176e-08, 0.3277477970321e+01, 0.7445550607224e+01],  
          
        [0.3851814176947e-08, 0.5210690074886e+01, 0.9562891316684e+00],  
        [0.3485807052785e-08, 0.6653274904611e+00, 0.1161697602389e+02],  
        [0.3979772816991e-08, 0.1767941436148e+01, 0.2277943724828e+02],  
        [0.3402607460500e-08, 0.3421746306465e+01, 0.1087398597200e+02],  
        [0.4049993000926e-08, 0.1127144787547e+01, 0.3163918923335e+00],  
        [0.3420511182382e-08, 0.4214794779161e+01, 0.1362553364512e+02],  
        [0.3640772365012e-08, 0.5324905497687e+01, 0.1725304118033e+02],  
        [0.3323037987501e-08, 0.6135761838271e+01, 0.6279143387820e+01],  
        [0.4503141663637e-08, 0.1802305450666e+01, 0.1385561574497e+01],  
        [0.4314560055588e-08, 0.4812299731574e+01, 0.4176041334900e+01],  
          
        [0.3294226949110e-08, 0.3657547059723e+01, 0.6287008313071e+01],  
        [0.3215657197281e-08, 0.4866676894425e+01, 0.5749861718712e+01],  
        [0.4129362656266e-08, 0.3809342558906e+01, 0.5905702259363e+01],  
        [0.3137762976388e-08, 0.2494635174443e+01, 0.2099539292909e+02],  
        [0.3514010952384e-08, 0.2699961831678e+01, 0.7335344340001e+01],  
        [0.3327607571530e-08, 0.3318457714816e+01, 0.5436992986000e+01],  
        [0.3541066946675e-08, 0.4382703582466e+01, 0.1234573916645e+02],  
        [0.3216179847052e-08, 0.5271066317054e+01, 0.3802769619140e-01],  
        [0.2959045059570e-08, 0.5819591585302e+01, 0.2670964694522e+02],  
        [0.3884040326665e-08, 0.5980934960428e+01, 0.6660449441528e+01],  
          
        [0.2922027539886e-08, 0.3337290282483e+01, 0.1375773836557e+01],  
        [0.4110846382042e-08, 0.5742978187327e+01, 0.4480965020977e+02],  
        [0.2934508411032e-08, 0.2278075804200e+01, 0.6408777551755e+00],  
        [0.3966896193000e-08, 0.5835747858477e+01, 0.3773735910827e+00],  
        [0.3286695827610e-08, 0.5838898193902e+01, 0.3932462625300e-02],  
        [0.3720643094196e-08, 0.1122212337858e+01, 0.1646033343740e+02],  
        [0.3285508906174e-08, 0.9182250996416e+00, 0.1081813534213e+02],  
        [0.3753880575973e-08, 0.5174761973266e+01, 0.5642198095270e+01],  
        [0.3022129385587e-08, 0.3381611020639e+01, 0.2982630633589e+02],  
        [0.2798569205621e-08, 0.3546193723922e+01, 0.1937891852345e+02],  
          
        [0.3397872070505e-08, 0.4533203197934e+01, 0.6923953605621e+01],  
        [0.3708099772977e-08, 0.2756168198616e+01, 0.3066615496545e+02],  
        [0.3599283541510e-08, 0.1934395469918e+01, 0.6147450479709e+01],  
        [0.3688702753059e-08, 0.7149920971109e+00, 0.2636725487657e+01],  
        [0.2681084724003e-08, 0.4899819493154e+01, 0.6816289982179e+01],  
        [0.3495993460759e-08, 0.1572418915115e+01, 0.6418701221183e+01],  
        [0.3130770324995e-08, 0.8912190180489e+00, 0.1235996607578e+02],  
        [0.2744353821941e-08, 0.3800821940055e+01, 0.2059724391010e+02],  
        [0.2842732906341e-08, 0.2644717440029e+01, 0.2828699048865e+02],  
        [0.3046882682154e-08, 0.3987793020179e+01, 0.6055599646783e+01],  
          
        [0.2399072455143e-08, 0.9908826440764e+00, 0.6255674361143e+01],  
        [0.2384306274204e-08, 0.2516149752220e+01, 0.6310477339748e+01],  
        [0.2977324500559e-08, 0.5849195642118e+01, 0.1652265972112e+02],  
        [0.3062835258972e-08, 0.1681660100162e+01, 0.1172006883645e+02],  
        [0.3109682589231e-08, 0.5804143987737e+00, 0.2751146787858e+02],  
        [0.2903920355299e-08, 0.5800768280123e+01, 0.6510552054109e+01],  
        [0.2823221989212e-08, 0.9241118370216e+00, 0.5469525544182e+01],  
        [0.3187949696649e-08, 0.3139776445735e+01, 0.1693792562116e+03],  
        [0.2922559771655e-08, 0.3549440782984e+01, 0.2630839062450e+00],  
        [0.2436302066603e-08, 0.4735540696319e+01, 0.3946258593675e+00],  
          
        [0.3049473043606e-08, 0.4998289124561e+01, 0.8390110365991e+01],  
        [0.2863682575784e-08, 0.6709515671102e+00, 0.2243449970715e+00],  
        [0.2641750517966e-08, 0.5410978257284e+01, 0.2986433403208e+02],  
        [0.2704093466243e-08, 0.4778317207821e+01, 0.6129297044991e+01],  
        [0.2445522177011e-08, 0.6009020662222e+01, 0.1171295538178e+02],  
        [0.2623608810230e-08, 0.5010449777147e+01, 0.6436854655901e+01],  
        [0.2079259704053e-08, 0.5980943768809e+01, 0.2019909489111e+02],  
        [0.2820225596771e-08, 0.2679965110468e+01, 0.5934151399930e+01],  
        [0.2365221950927e-08, 0.1894231148810e+01, 0.2470570524223e+02],  
        [0.2359682077149e-08, 0.4220752950780e+01, 0.8671969964381e+01],  
          
        [0.2387577137206e-08, 0.2571783940617e+01, 0.7096626156709e+01],  
        [0.1982102089816e-08, 0.5169765997119e+00, 0.1727188400790e+02],  
        [0.2687502389925e-08, 0.6239078264579e+01, 0.7075506709219e+02],  
        [0.2207751669135e-08, 0.2031184412677e+01, 0.4377611041777e+01],  
        [0.2618370214274e-08, 0.8266079985979e+00, 0.6632000300961e+01],  
        [0.2591951887361e-08, 0.8819350522008e+00, 0.4873985990671e+02],  
        [0.2375055656248e-08, 0.3520944177789e+01, 0.1590676413561e+02],  
        [0.2472019978911e-08, 0.1551431908671e+01, 0.6612329252343e+00],  
        [0.2368157127199e-08, 0.4178610147412e+01, 0.3459636466239e+02],  
        [0.1764846605693e-08, 0.1506764000157e+01, 0.1980094587212e+02],  
          
        [0.2291769608798e-08, 0.2118250611782e+01, 0.2844914056730e-01],  
        [0.2209997316943e-08, 0.3363255261678e+01, 0.2666070658668e+00],  
        [0.2292699097923e-08, 0.4200423956460e+00, 0.1484170571900e-02],  
        [0.1629683015329e-08, 0.2331362582487e+01, 0.3035599730800e+02],  
        [0.2206492862426e-08, 0.3400274026992e+01, 0.6281667977667e+01],  
        [0.2205746568257e-08, 0.1066051230724e+00, 0.6284483723224e+01],  
        [0.2026310767991e-08, 0.2779066487979e+01, 0.2449240616245e+02],  
        [0.1762977622163e-08, 0.9951450691840e+00, 0.2045286941806e+02],  
        [0.1368535049606e-08, 0.6402447365817e+00, 0.2473415438279e+02],  
        [0.1720598775450e-08, 0.2303524214705e+00, 0.1679593901136e+03],  
          
        [0.1702429015449e-08, 0.6164622655048e+01, 0.3338575901272e+03],  
        [0.1414033197685e-08, 0.3954561185580e+01, 0.1624205518357e+03],  
        [0.1573768958043e-08, 0.2028286308984e+01, 0.3144167757552e+02],  
        [0.1650705184447e-08, 0.2304040666128e+01, 0.5267006960365e+02],  
        [0.1651087618855e-08, 0.2538461057280e+01, 0.8956999012000e+02],  
        [0.1616409518983e-08, 0.5111054348152e+01, 0.3332657872986e+02],  
        [0.1537175173581e-08, 0.5601130666603e+01, 0.3852657435933e+02],  
        [0.1593191980553e-08, 0.2614340453411e+01, 0.2282781046519e+03],  
        [0.1499480170643e-08, 0.3624721577264e+01, 0.2823723341956e+02],  
        [0.1493807843235e-08, 0.4214569879008e+01, 0.2876692439167e+02],  
          
        [0.1074571199328e-08, 0.1496911744704e+00, 0.8397383534231e+02],  
        [0.1074406983417e-08, 0.1187817671922e+01, 0.8401985929482e+02],  
        [0.9757576855851e-09, 0.2655703035858e+01, 0.7826370942180e+02],  
        [0.1258432887565e-08, 0.4969896184844e+01, 0.3115650189215e+03],  
        [0.1240336343282e-08, 0.5192460776926e+01, 0.1784300471910e+03],  
        [0.9016107005164e-09, 0.1960356923057e+01, 0.5886454391678e+02],  
        [0.1135392360918e-08, 0.5082427809068e+01, 0.7842370451713e+02],  
        [0.9216046089565e-09, 0.2793775037273e+01, 0.1014262087719e+03],  
        [0.1061276615030e-08, 0.3726144311409e+01, 0.5660027930059e+02],  
        [0.1010110596263e-08, 0.7404080708937e+00, 0.4245678405627e+02],  
          
        [0.7217424756199e-09, 0.2697449980577e-01, 0.2457074661053e+03],  
        [0.6912003846756e-09, 0.4253296276335e+01, 0.1679936946371e+03],  
        [0.6871814664847e-09, 0.5148072412354e+01, 0.6053048899753e+02],  
        [0.4887158016343e-09, 0.2153581148294e+01, 0.9656299901946e+02],  
        [0.5161802866314e-09, 0.3852750634351e+01, 0.2442876000072e+03],  
        [0.5652599559057e-09, 0.1233233356270e+01, 0.8365903305582e+02],  
        [0.4710812608586e-09, 0.5610486976767e+01, 0.3164282286739e+03],  
        [0.4909977500324e-09, 0.1639629524123e+01, 0.4059982187939e+03],  
        [0.4772641839378e-09, 0.3737100368583e+01, 0.1805255418145e+03],  
        [0.4487562567153e-09, 0.1158417054478e+00, 0.8433466158131e+02],  
          
        [0.3943441230497e-09, 0.6243502862796e+00, 0.2568537517081e+03],  
        [0.3952236913598e-09, 0.3510377382385e+01, 0.2449975330562e+03],  
        [0.3788898363417e-09, 0.5916128302299e+01, 0.1568131045107e+03],  
        [0.3738329328831e-09, 0.1042266763456e+01, 0.3948519331910e+03],  
        [0.2451199165151e-09, 0.1166788435700e+01, 0.1435713242844e+03],  
        [0.2436734402904e-09, 0.3254726114901e+01, 0.2268582385539e+03],  
        [0.2213605274325e-09, 0.1687210598530e+01, 0.1658638954901e+03],  
        [0.1491521204829e-09, 0.2657541786794e+01, 0.2219950288015e+03],  
        [0.1474995329744e-09, 0.5013089805819e+01, 0.3052819430710e+03],  
        [0.1661939475656e-09, 0.5495315428418e+01, 0.2526661704812e+03],  
          
        [0.9015946748003e-10, 0.2236989966505e+01, 0.4171445043968e+03] 
         
    ])

    # Sun-to-Earth, T^0, Z
    e0z = np.array([
        [0.2796207639075e-05, 0.3198701560209e+01, 0.8433466158131e+02],
        [0.1016042198142e-05, 0.5422360395913e+01, 0.5507553240374e+01],
        [0.8044305033647e-06, 0.3880222866652e+01, 0.5223693906222e+01],
        [0.4385347909274e-06, 0.3704369937468e+01, 0.2352866153506e+01],
        [0.3186156414906e-06, 0.3999639363235e+01, 0.1577343543434e+01],
        [0.2272412285792e-06, 0.3984738315952e+01, 0.1047747311755e+01],
        [0.1645620103007e-06, 0.3565412516841e+01, 0.5856477690889e+01],
        [0.1815836921166e-06, 0.4984507059020e+01, 0.6283075850446e+01],
        [0.1447461676364e-06, 0.3702753570108e+01, 0.9437762937313e+01],
        [0.1430760876382e-06, 0.3409658712357e+01, 0.1021328554739e+02],
        
        [0.1120445753226e-06, 0.4829561570246e+01, 0.1414349524433e+02],
        [0.1090232840797e-06, 0.2080729178066e+01, 0.6812766822558e+01],
        [0.9715727346551e-07, 0.3476295881948e+01, 0.4694002934110e+01],
        [0.1036267136217e-06, 0.4056639536648e+01, 0.7109288135493e+02],
        [0.8752665271340e-07, 0.4448159519911e+01, 0.5753384878334e+01],
        [0.8331864956004e-07, 0.4991704044208e+01, 0.7084896783808e+01],
        [0.6901658670245e-07, 0.4325358994219e+01, 0.6275962395778e+01],
        [0.9144536848998e-07, 0.1141826375363e+01, 0.6620890113188e+01],
        [0.7205085037435e-07, 0.3624344170143e+01, 0.5296909721118e+00],
        [0.7697874654176e-07, 0.5554257458998e+01, 0.1676215758509e+03],
        
        [0.5197545738384e-07, 0.6251760961735e+01, 0.1807370494127e+02],
        [0.5031345378608e-07, 0.2497341091913e+01, 0.4705732307012e+01],
        [0.4527110205840e-07, 0.2335079920992e+01, 0.6309374173736e+01],
        [0.4753355798089e-07, 0.7094148987474e+00, 0.5884926831456e+01],
        [0.4296951977516e-07, 0.1101916352091e+01, 0.6681224869435e+01],
        [0.3855341568387e-07, 0.1825495405486e+01, 0.5486777812467e+01],
        [0.5253930970990e-07, 0.4424740687208e+01, 0.7860419393880e+01],
        [0.4024630496471e-07, 0.5120498157053e+01, 0.1336797263425e+02],
        [0.4061069791453e-07, 0.6029771435451e+01, 0.3930209696940e+01],
        [0.3797883804205e-07, 0.4435193600836e+00, 0.3154687086868e+01],
        
        [0.2933033225587e-07, 0.5124157356507e+01, 0.1059381944224e+01],
        [0.3503000930426e-07, 0.5421830162065e+01, 0.6069776770667e+01],
        [0.3670096214050e-07, 0.4582101667297e+01, 0.1219403291462e+02],
        [0.2905609437008e-07, 0.1926566420072e+01, 0.1097707878456e+02],
        [0.2466827821713e-07, 0.6090174539834e+00, 0.6496374930224e+01],
        [0.2691647295332e-07, 0.1393432595077e+01, 0.2200391463820e+02],
        [0.2150554667946e-07, 0.4308671715951e+01, 0.5643178611111e+01],
        [0.2237481922680e-07, 0.8133968269414e+00, 0.8635942003952e+01],
        [0.1817741038157e-07, 0.3755205127454e+01, 0.3340612434717e+01],
        [0.2227820762132e-07, 0.2759558596664e+01, 0.1203646072878e+02],
        
        [0.1944713772307e-07, 0.5699645869121e+01, 0.1179062909082e+02],
        [0.1527340520662e-07, 0.1986749091746e+01, 0.3981490189893e+00],
        [0.1577282574914e-07, 0.3205017217983e+01, 0.5088628793478e+01],
        [0.1424738825424e-07, 0.6256747903666e+01, 0.2544314396739e+01],
        [0.1616563121701e-07, 0.2601671259394e+00, 0.1729818233119e+02],
        [0.1401210391692e-07, 0.4686939173506e+01, 0.7058598460518e+01],
        [0.1488726974214e-07, 0.2815862451372e+01, 0.2593412433514e+02],
        [0.1692626442388e-07, 0.4956894109797e+01, 0.1564752902480e+03],
        [0.1123571582910e-07, 0.2381192697696e+01, 0.3738761453707e+01],
        [0.9903308606317e-08, 0.4294851657684e+01, 0.9225539266174e+01],
        
        [0.9174533187191e-08, 0.3075171510642e+01, 0.4164311961999e+01],
        [0.8645985631457e-08, 0.5477534821633e+00, 0.8429241228195e+01],
        [-0.1085876492688e-07,0.0000000000000e+00, 0.0000000000000e+00],
        [0.9264309077815e-08, 0.5968571670097e+01, 0.7079373888424e+01],
        [0.8243116984954e-08, 0.1489098777643e+01, 0.1044738781244e+02],
        [0.8268102113708e-08, 0.3512977691983e+01, 0.1150676975667e+02],
        [0.9043613988227e-08, 0.1290704408221e+00, 0.1101510648075e+02],
        [0.7432912038789e-08, 0.1991086893337e+01, 0.2608790314060e+02],
        [0.8586233727285e-08, 0.4238357924414e+01, 0.2986433403208e+02],
        [0.7612230060131e-08, 0.2911090150166e+01, 0.4732030630302e+01],
        
        [0.7097787751408e-08, 0.1908938392390e+01, 0.8031092209206e+01],
        [0.7640237040175e-08, 0.6129219000168e+00, 0.7962980379786e+00],
        [0.7070445688081e-08, 0.1380417036651e+01, 0.2146165377750e+01],
        [0.7690770957702e-08, 0.1680504249084e+01, 0.2122839202813e+02],
        [0.8051292542594e-08, 0.5127423484511e+01, 0.2942463415728e+01],
        [0.5902709104515e-08, 0.2020274190917e+01, 0.7755226100720e+00],
        [0.5134567496462e-08, 0.2606778676418e+01, 0.1256615170089e+02],
        [0.5525802046102e-08, 0.1613011769663e+01, 0.8018209333619e+00],
        [0.5880724784221e-08, 0.4604483417236e+01, 0.4690479774488e+01],
        [0.5211699081370e-08, 0.5718964114193e+01, 0.8827390247185e+01],
        
        [0.4891849573562e-08, 0.3689658932196e+01, 0.2132990797783e+00],
        [0.5150246069997e-08, 0.4099769855122e+01, 0.6480980550449e+02],
        [0.5102434319633e-08, 0.5660834602509e+01, 0.3379454372902e+02],
        [0.5083405254252e-08, 0.9842221218974e+00, 0.4136910472696e+01],
        [0.4206562585682e-08, 0.1341363634163e+00, 0.3128388763578e+01],
        [0.4663249683579e-08, 0.8130132735866e+00, 0.5216580451554e+01],
        [0.4099474416530e-08, 0.5791497770644e+01, 0.4265981595566e+00],
        [0.4628251220767e-08, 0.1249802769331e+01, 0.1572083878776e+02],
        [0.5024068728142e-08, 0.4795684802743e+01, 0.6290189305114e+01],
        [0.5120234327758e-08, 0.3810420387208e+01, 0.5230807360890e+01],
        
        [0.5524029815280e-08, 0.1029264714351e+01, 0.2397622045175e+03],
        [0.4757415718860e-08, 0.3528044781779e+01, 0.1649636139783e+02],
        [0.3915786131127e-08, 0.5593889282646e+01, 0.1589072916335e+01],
        [0.4869053149991e-08, 0.3299636454433e+01, 0.7632943190217e+01],
        [0.3649365703729e-08, 0.1286049002584e+01, 0.6206810014183e+01],
        [0.3992493949002e-08, 0.3100307589464e+01, 0.2515860172507e+02],
        [0.3320247477418e-08, 0.6212683940807e+01, 0.1216800268190e+02],
        [0.3287123739696e-08, 0.4699118445928e+01, 0.7234794171227e+01],
        [0.3472776811103e-08, 0.2630507142004e+01, 0.7342457794669e+01],
        [0.3423253294767e-08, 0.2946432844305e+01, 0.9623688285163e+01],
        
        [0.3896173898244e-08, 0.1224834179264e+01, 0.6438496133249e+01],
        [0.3388455337924e-08, 0.1543807616351e+01, 0.1494531617769e+02],
        [0.3062704716523e-08, 0.1191777572310e+01, 0.8662240327241e+01],
        [0.3270075600400e-08, 0.5483498767737e+01, 0.1194447056968e+01],
        [0.3101209215259e-08, 0.8000833804348e+00, 0.3772475342596e+02],
        [0.2780883347311e-08, 0.4077980721888e+00, 0.5863591145557e+01],
        [0.2903605931824e-08, 0.2617490302147e+01, 0.1965104848470e+02],
        [0.2682014743119e-08, 0.2634703158290e+01, 0.7238675589263e+01],
        [0.2534360108492e-08, 0.6102446114873e+01, 0.6836645152238e+01],
        [0.2392564882509e-08, 0.3681820208691e+01, 0.5849364236221e+01],
        
        [0.2656667254856e-08, 0.6216045388886e+01, 0.6133512519065e+01],
        [0.2331242096773e-08, 0.5864949777744e+01, 0.4535059491685e+01],
        [0.2287898363668e-08, 0.4566628532802e+01, 0.7477522907414e+01],
        [0.2336944521306e-08, 0.2442722126930e+01, 0.1137170464392e+02],
        [0.3156632236269e-08, 0.1626628050682e+01, 0.2509084901204e+03],
        [0.2982612402766e-08, 0.2803604512609e+01, 0.1748016358760e+01],
        [0.2774031674807e-08, 0.4654002897158e+01, 0.8223916695780e+02],
        [0.2295236548638e-08, 0.4326518333253e+01, 0.3378142627421e+00],
        [0.2190714699873e-08, 0.4519614578328e+01, 0.2908881142201e+02],
        [0.2191495845045e-08, 0.3012626912549e+01, 0.1673046366289e+02],
        
        [0.2492901628386e-08, 0.1290101424052e+00, 0.1543797956245e+03],
        [0.1993778064319e-08, 0.3864046799414e+01, 0.1778984560711e+02],
        [0.1898146479022e-08, 0.5053777235891e+01, 0.2042657109477e+02],
        [0.1918280127634e-08, 0.2222470192548e+01, 0.4165496312290e+02],
        [0.1916351061607e-08, 0.8719067257774e+00, 0.7737595720538e+02],
        [0.1834720181466e-08, 0.4031491098040e+01, 0.2358125818164e+02],
        [0.1249201523806e-08, 0.5938379466835e+01, 0.3301902111895e+02],
        [0.1477304050539e-08, 0.6544722606797e+00, 0.9548094718417e+02],
        [0.1264316431249e-08, 0.2059072853236e+01, 0.8399684731857e+02],
        [0.1203526495039e-08, 0.3644813532605e+01, 0.4558517281984e+02],
        
        [0.9221681059831e-09, 0.3241815055602e+01, 0.7805158573086e+02],
        [0.7849278367646e-09, 0.5043812342457e+01, 0.5217580628120e+02],
        [0.7983392077387e-09, 0.5000024502753e+01, 0.1501922143975e+03],
        [0.7925395431654e-09, 0.1398734871821e-01, 0.9061773743175e+02],
        [0.7640473285886e-09, 0.5067111723130e+01, 0.4951538251678e+02],
        [0.5398937754482e-09, 0.5597382200075e+01, 0.1613385000004e+03],
        [0.5626247550193e-09, 0.2601338209422e+01, 0.7318837597844e+02],
        [0.5525197197855e-09, 0.5814832109256e+01, 0.1432335100216e+03],
        [0.5407629837898e-09, 0.3384820609076e+01, 0.3230491187871e+03],
        [0.3856739119801e-09, 0.1072391840473e+01, 0.2334791286671e+03],
        
        [0.3856425239987e-09, 0.2369540393327e+01, 0.1739046517013e+03],
        [0.4350867755983e-09, 0.5255575751082e+01, 0.1620484330494e+03],
        [0.3844113924996e-09, 0.5482356246182e+01, 0.9757644180768e+02],
        [0.2854869155431e-09, 0.9573634763143e+00, 0.1697170704744e+03],
        [0.1719227671416e-09, 0.1887203025202e+01, 0.2265204242912e+03],
        [0.1527846879755e-09, 0.3982183931157e+01, 0.3341954043900e+03],
        [0.1128229264847e-09, 0.2787457156298e+01, 0.3119028331842e+03]  
        
    ])

    # Sun-to-Earth, T^1, X
    e1x = np.array([
        [0.1234046326004e-05, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.5150068824701e-06, 0.6002664557501e+01, 0.1256615170089e+02],
        [0.1290743923245e-07, 0.5959437664199e+01, 0.1884922755134e+02],
        [0.1068615564952e-07, 0.2015529654209e+01, 0.6283075850446e+01],
        [0.2079619142538e-08, 0.1732960531432e+01, 0.6279552690824e+01],
        [0.2078009243969e-08, 0.4915604476996e+01, 0.6286599010068e+01],
        [0.6206330058856e-09, 0.3616457953824e+00, 0.4705732307012e+01],
        [0.5989335313746e-09, 0.3802607304474e+01, 0.6256777527156e+01],
        [0.5958495663840e-09, 0.2845866560031e+01, 0.6309374173736e+01],
        [0.4866923261539e-09, 0.5213203771824e+01, 0.7755226100720e+00],
        
        [0.4267785823142e-09, 0.4368189727818e+00, 0.1059381944224e+01],
        [0.4610675141648e-09, 0.1837249181372e-01, 0.7860419393880e+01],
        [0.3626989993973e-09, 0.2161590545326e+01, 0.5753384878334e+01],
        [0.3563071194389e-09, 0.1452631954746e+01, 0.5884926831456e+01],
        [0.3557015642807e-09, 0.4470593393054e+01, 0.6812766822558e+01],
        [0.3210412089122e-09, 0.5195926078314e+01, 0.6681224869435e+01],
        [0.2875473577986e-09, 0.5916256610193e+01, 0.2513230340178e+02],
        [0.2842913681629e-09, 0.1149902426047e+01, 0.6127655567643e+01],
        [0.2751248215916e-09, 0.5502088574662e+01, 0.6438496133249e+01],
        [0.2481432881127e-09, 0.2921989846637e+01, 0.5486777812467e+01],
        
        [0.2059885976560e-09, 0.3718070376585e+01, 0.7079373888424e+01],
        [0.2015522342591e-09, 0.5979395259740e+01, 0.6290189305114e+01],
        [0.1995364084253e-09, 0.6772087985494e+00, 0.6275962395778e+01],
        [0.1957436436943e-09, 0.2899210654665e+01, 0.5507553240374e+01],
        [0.1651609818948e-09, 0.6228206482192e+01, 0.1150676975667e+02],
        [0.1822980550699e-09, 0.1469348746179e+01, 0.1179062909082e+02],
        [0.1675223159760e-09, 0.3813910555688e+01, 0.7058598460518e+01],
        [0.1706491764745e-09, 0.3004380506684e+00, 0.7113454667900e-02],
        [0.1392952362615e-09, 0.1440393973406e+01, 0.7962980379786e+00],
        [0.1209868266342e-09, 0.4150425791727e+01, 0.4694002934110e+01],
        
        [0.1009827202611e-09, 0.3290040429843e+01, 0.3738761453707e+01],
        [0.1047261388602e-09, 0.4229590090227e+01, 0.6282095334605e+01],
        [0.1047006652004e-09, 0.2418967680575e+01, 0.6284056366286e+01],
        [0.9609993143095e-10, 0.4627943659201e+01, 0.6069776770667e+01],
        [0.9590900593873e-10, 0.1894393939924e+01, 0.4136910472696e+01],
        [0.9146249188071e-10, 0.2010647519562e+01, 0.6496374930224e+01],
        [0.8545274480290e-10, 0.5529846956226e-01, 0.1194447056968e+01],
        [0.8224377881194e-10, 0.1254304102174e+01, 0.1589072916335e+01],
        [0.6183529510410e-10, 0.3360862168815e+01, 0.8827390247185e+01],
        [0.6259255147141e-10, 0.4755628243179e+01, 0.8429241228195e+01],
        
        [0.5539291694151e-10, 0.5371746955142e+01, 0.4933208510675e+01],
        [0.7328259466314e-10, 0.4927699613906e+00, 0.4535059491685e+01],
        [0.6017835843560e-10, 0.5776682001734e-01, 0.1255903824622e+02],
        [0.7079827775243e-10, 0.4395059432251e+01, 0.5088628793478e+01],
        [0.5170358878213e-10, 0.5154062619954e+01, 0.1176985366291e+02],
        [0.4872301838682e-10, 0.6289611648973e+00, 0.6040347114260e+01],
        [0.5249869411058e-10, 0.5617272046949e+01, 0.3154687086868e+01],
        [0.4716172354411e-10, 0.3965901800877e+01, 0.5331357529664e+01],
        [0.4871214940964e-10, 0.4627507050093e+01, 0.1256967486051e+02],
        [0.4598076850751e-10, 0.6023631226459e+01, 0.6525804586632e+01],
        
        [0.4562196089485e-10, 0.4138562084068e+01, 0.3930209696940e+01],
        [0.4325493872224e-10, 0.1330845906564e+01, 0.7632943190217e+01],
        [0.5673781176748e-10, 0.2558752615657e+01, 0.5729506548653e+01],
        [0.3961436642503e-10, 0.2728071734630e+01, 0.7234794171227e+01],
        [0.5101868209058e-10, 0.4113444965144e+01, 0.6836645152238e+01],
        [0.5257043167676e-10, 0.6195089830590e+01, 0.8031092209206e+01],
        [0.5076613989393e-10, 0.2305124132918e+01, 0.7477522907414e+01],
        [0.3342169352778e-10, 0.5415998155071e+01, 0.1097707878456e+02],
        [0.3545881983591e-10, 0.3727160564574e+01, 0.4164311961999e+01],
        [0.3364063738599e-10, 0.2901121049204e+00, 0.1137170464392e+02],
        
        [0.3357039670776e-10, 0.1652229354331e+01, 0.5223693906222e+01],
        [0.4307412268687e-10, 0.4938909587445e+01, 0.1592596075957e+01],
        [0.3405769115435e-10, 0.2408890766511e+01, 0.3128388763578e+01],
        [0.3001926198480e-10, 0.4862239006386e+01, 0.1748016358760e+01],
        [0.2778264787325e-10, 0.5241168661353e+01, 0.7342457794669e+01],
        [0.2676159480666e-10, 0.3423593942199e+01, 0.2146165377750e+01],
        [0.2954273399939e-10, 0.1881721265406e+01, 0.5368044267797e+00],
        [0.3309362888795e-10, 0.1931525677349e+01, 0.8018209333619e+00],
        [0.2810283608438e-10, 0.2414659495050e+01, 0.5225775174439e+00],
        [0.3378045637764e-10, 0.4238019163430e+01, 0.1554202828031e+00],
        
        [0.2558134979840e-10, 0.1828225235805e+01, 0.5230807360890e+01],
        [0.2273755578447e-10, 0.5858184283998e+01, 0.7084896783808e+01],
        [0.2294176037690e-10, 0.4514589779057e+01, 0.1726015463500e+02],
        [0.2533506099435e-10, 0.2355717851551e+01, 0.5216580451554e+01],
        [0.2716685375812e-10, 0.2221003625100e+01, 0.8635942003952e+01],
        [0.2419043435198e-10, 0.5955704951635e+01, 0.4690479774488e+01],
        [0.2521232544812e-10, 0.1395676848521e+01, 0.5481254917084e+01],
        [0.2630195021491e-10, 0.5727468918743e+01, 0.2629832328990e-01],
        [0.2548395840944e-10, 0.2628351859400e-03, 0.1349867339771e+01] 
    ])

    # Sun-to-Earth, T^1, Y
    e1y = np.array([
        [0.9304690546528e-06, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.5150715570663e-06, 0.4431807116294e+01, 0.1256615170089e+02],
        [0.1290825411056e-07, 0.4388610039678e+01, 0.1884922755134e+02],
        [0.4645466665386e-08, 0.5827263376034e+01, 0.6283075850446e+01],
        [0.2079625310718e-08, 0.1621698662282e+00, 0.6279552690824e+01],
        [0.2078189850907e-08, 0.3344713435140e+01, 0.6286599010068e+01],
        [0.6207190138027e-09, 0.5074049319576e+01, 0.4705732307012e+01],
        [0.5989826532569e-09, 0.2231842216620e+01, 0.6256777527156e+01],
        [0.5961360812618e-09, 0.1274975769045e+01, 0.6309374173736e+01],
        [0.4874165471016e-09, 0.3642277426779e+01, 0.7755226100720e+00],
        
        [0.4283834034360e-09, 0.5148765510106e+01, 0.1059381944224e+01],
        [0.4652389287529e-09, 0.4715794792175e+01, 0.7860419393880e+01],
        [0.3751707476401e-09, 0.6617207370325e+00, 0.5753384878334e+01],
        [0.3559998806198e-09, 0.6155548875404e+01, 0.5884926831456e+01],
        [0.3558447558857e-09, 0.2898827297664e+01, 0.6812766822558e+01],
        [0.3211116927106e-09, 0.3625813502509e+01, 0.6681224869435e+01],
        [0.2875609914672e-09, 0.4345435813134e+01, 0.2513230340178e+02],
        [0.2843109704069e-09, 0.5862263940038e+01, 0.6127655567643e+01],
        [0.2744676468427e-09, 0.3926419475089e+01, 0.6438496133249e+01],
        [0.2481285237789e-09, 0.1351976572828e+01, 0.5486777812467e+01],
        
        [0.2060338481033e-09, 0.2147556998591e+01, 0.7079373888424e+01],
        [0.2015822358331e-09, 0.4408358972216e+01, 0.6290189305114e+01],
        [0.2001195944195e-09, 0.5385829822531e+01, 0.6275962395778e+01],
        [0.1953667642377e-09, 0.1304933746120e+01, 0.5507553240374e+01],
        [0.1839744078713e-09, 0.6173567228835e+01, 0.1179062909082e+02],
        [0.1643334294845e-09, 0.4635942997523e+01, 0.1150676975667e+02],
        [0.1768051018652e-09, 0.5086283558874e+01, 0.7113454667900e-02],
        [0.1674874205489e-09, 0.2243332137241e+01, 0.7058598460518e+01],
        [0.1421445397609e-09, 0.6186899771515e+01, 0.7962980379786e+00],
        [0.1255163958267e-09, 0.5730238465658e+01, 0.4694002934110e+01],
        
        [0.1013945281961e-09, 0.1726055228402e+01, 0.3738761453707e+01],
        [0.1047294335852e-09, 0.2658801228129e+01, 0.6282095334605e+01],
        [0.1047103879392e-09, 0.8481047835035e+00, 0.6284056366286e+01],
        [0.9530343962826e-10, 0.3079267149859e+01, 0.6069776770667e+01],
        [0.9604637611690e-10, 0.3258679792918e+00, 0.4136910472696e+01],
        [0.9153518537177e-10, 0.4398599886584e+00, 0.6496374930224e+01],
        [0.8562458214922e-10, 0.4772686794145e+01, 0.1194447056968e+01],
        [0.8232525360654e-10, 0.5966220721679e+01, 0.1589072916335e+01],
        [0.6150223411438e-10, 0.1780985591923e+01, 0.8827390247185e+01],
        [0.6272087858000e-10, 0.3184305429012e+01, 0.8429241228195e+01],
        
        [0.5540476311040e-10, 0.3801260595433e+01, 0.4933208510675e+01],
        [0.7331901699361e-10, 0.5205948591865e+01, 0.4535059491685e+01],
        [0.6018528702791e-10, 0.4770139083623e+01, 0.1255903824622e+02],
        [0.5150530724804e-10, 0.3574796899585e+01, 0.1176985366291e+02],
        [0.6471933741811e-10, 0.2679787266521e+01, 0.5088628793478e+01],
        [0.5317460644174e-10, 0.9528763345494e+00, 0.3154687086868e+01],
        [0.4832187748783e-10, 0.5329322498232e+01, 0.6040347114260e+01],
        [0.4716763555110e-10, 0.2395235316466e+01, 0.5331357529664e+01],
        [0.4871509139861e-10, 0.3056663648823e+01, 0.1256967486051e+02],
        [0.4598417696768e-10, 0.4452762609019e+01, 0.6525804586632e+01],
        
        [0.5674189533175e-10, 0.9879680872193e+00, 0.5729506548653e+01],
        [0.4073560328195e-10, 0.5939127696986e+01, 0.7632943190217e+01],
        [0.5040994945359e-10, 0.4549875824510e+01, 0.8031092209206e+01],
        [0.5078185134679e-10, 0.7346659893982e+00, 0.7477522907414e+01],
        [0.3769343537061e-10, 0.1071317188367e+01, 0.7234794171227e+01],
        [0.4980331365299e-10, 0.2500345341784e+01, 0.6836645152238e+01],
        [0.3458236594757e-10, 0.3825159450711e+01, 0.1097707878456e+02],
        [0.3578859493602e-10, 0.5299664791549e+01, 0.4164311961999e+01],
        [0.3370504646419e-10, 0.5002316301593e+01, 0.1137170464392e+02],
        [0.3299873338428e-10, 0.2526123275282e+01, 0.3930209696940e+01],
        
        [0.4304917318409e-10, 0.3368078557132e+01, 0.1592596075957e+01],
        [0.3402418753455e-10, 0.8385495425800e+00, 0.3128388763578e+01],
        [0.2778460572146e-10, 0.3669905203240e+01, 0.7342457794669e+01],
        [0.2782710128902e-10, 0.2691664812170e+00, 0.1748016358760e+01],
        [0.2711725179646e-10, 0.4707487217718e+01, 0.5296909721118e+00],
        [0.2981760946340e-10, 0.3190260867816e+00, 0.5368044267797e+00],
        [0.2811672977772e-10, 0.3196532315372e+01, 0.7084896783808e+01],
        [0.2863454474467e-10, 0.2263240324780e+00, 0.5223693906222e+01],
        [0.3333464634051e-10, 0.3498451685065e+01, 0.8018209333619e+00],
        [0.3312991747609e-10, 0.5839154477412e+01, 0.1554202828031e+00],
        
        [0.2813255564006e-10, 0.8268044346621e+00, 0.5225775174439e+00],
        [0.2665098083966e-10, 0.3934021725360e+01, 0.5216580451554e+01],
        [0.2349795705216e-10, 0.5197620913779e+01, 0.2146165377750e+01],
        [0.2330352293961e-10, 0.2984999231807e+01, 0.1726015463500e+02],
        [0.2728001683419e-10, 0.6521679638544e+00, 0.8635942003952e+01],
        [0.2484061007669e-10, 0.3468955561097e+01, 0.5230807360890e+01],
        [0.2646328768427e-10, 0.1013724533516e+01, 0.2629832328990e-01],
        [0.2518630264831e-10, 0.6108081057122e+01, 0.5481254917084e+01],
        [0.2421901455384e-10, 0.1651097776260e+01, 0.1349867339771e+01],
        [0.6348533267831e-11, 0.3220226560321e+01, 0.8433466158131e+02] 
    ])

    # Sun-to-Earth, T^1, Z
    e1z = np.array([
        [0.2278290449966e-05, 0.3413716033863e+01, 0.6283075850446e+01],
        [0.5429458209830e-07, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.1903240492525e-07, 0.3370592358297e+01, 0.1256615170089e+02],
        [0.2385409276743e-09, 0.3327914718416e+01, 0.1884922755134e+02],
        [0.8676928342573e-10, 0.1824006811264e+01, 0.5223693906222e+01],
        [0.7765442593544e-10, 0.3888564279247e+01, 0.5507553240374e+01],
        [0.7066158332715e-10, 0.5194267231944e+01, 0.2352866153506e+01],
        [0.7092175288657e-10, 0.2333246960021e+01, 0.8399684731857e+02],
        [0.5357582213535e-10, 0.2224031176619e+01, 0.5296909721118e+00],
        [0.3828035865021e-10, 0.2156710933584e+01, 0.6279552690824e+01],
        
        [0.3824857220427e-10, 0.1529755219915e+01, 0.6286599010068e+01],
        [0.3286995181628e-10, 0.4879512900483e+01, 0.1021328554739e+02] 
    ])

    # Sun-to-Earth, T^2, X
    e2x = np.array([
        [-0.4143818297913e-10,0.0000000000000e+00, 0.0000000000000e+00],
        [0.2171497694435e-10, 0.4398225628264e+01, 0.1256615170089e+02],
        [0.9845398442516e-11, 0.2079720838384e+00, 0.6283075850446e+01],
        [0.9256833552682e-12, 0.4191264694361e+01, 0.1884922755134e+02],
        [0.1022049384115e-12, 0.5381133195658e+01, 0.8399684731857e+02]
    ])

    # Sun-to-Earth, T^2, Y
    e2y = np.array([
        [0.5063375872532e-10, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.2173815785980e-10, 0.2827805833053e+01, 0.1256615170089e+02],
        [0.1010231999920e-10, 0.4634612377133e+01, 0.6283075850446e+01],
        [0.9259745317636e-12, 0.2620612076189e+01, 0.1884922755134e+02],
        [0.1022202095812e-12, 0.3809562326066e+01, 0.8399684731857e+02]
    ])

    # Sun-to-Earth, T^2, Z
    e2z = np.array([
        [0.9722666114891e-10, 0.5152219582658e+01, 0.6283075850446e+01],
        [-0.3494819171909e-11,0.0000000000000e+00, 0.0000000000000e+00],
        [0.6713034376076e-12, 0.6440188750495e+00, 0.1256615170089e+02]
    ])

    # SSB-to-Sun, T^0, X
    s0x = np.array([
        [0.4956757536410e-02, 0.3741073751789e+01, 0.5296909721118e+00],
        [0.2718490072522e-02, 0.4016011511425e+01, 0.2132990797783e+00],
        [0.1546493974344e-02, 0.2170528330642e+01, 0.3813291813120e-01],
        [0.8366855276341e-03, 0.2339614075294e+01, 0.7478166569050e-01],
        [0.2936777942117e-03, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.1201317439469e-03, 0.4090736353305e+01, 0.1059381944224e+01],
        [0.7578550887230e-04, 0.3241518088140e+01, 0.4265981595566e+00],
        [0.1941787367773e-04, 0.1012202064330e+01, 0.2061856251104e+00],
        [0.1889227765991e-04, 0.3892520416440e+01, 0.2204125344462e+00],
        [0.1937896968613e-04, 0.4797779441161e+01, 0.1495633313810e+00],
        
        [0.1434506110873e-04, 0.3868960697933e+01, 0.5225775174439e+00],
        [0.1406659911580e-04, 0.4759766557397e+00, 0.5368044267797e+00],
        [0.1179022300202e-04, 0.7774961520598e+00, 0.7626583626240e-01],
        [0.8085864460959e-05, 0.3254654471465e+01, 0.3664874755930e-01],
        [0.7622752967615e-05, 0.4227633103489e+01, 0.3961708870310e-01],
        [0.6209171139066e-05, 0.2791828325711e+00, 0.7329749511860e-01],
        [0.4366435633970e-05, 0.4440454875925e+01, 0.1589072916335e+01],
        [0.3792124889348e-05, 0.5156393842356e+01, 0.7113454667900e-02],
        [0.3154548963402e-05, 0.6157005730093e+01, 0.4194847048887e+00],
        [0.3088359882942e-05, 0.2494567553163e+01, 0.6398972393349e+00],
        
        [0.2788440902136e-05, 0.4934318747989e+01, 0.1102062672231e+00],
        [0.3039928456376e-05, 0.4895077702640e+01, 0.6283075850446e+01],
        [0.2272258457679e-05, 0.5278394064764e+01, 0.1030928125552e+00],
        [0.2162007057957e-05, 0.5802978019099e+01, 0.3163918923335e+00],
        [0.1767632855737e-05, 0.3415346595193e-01, 0.1021328554739e+02],
        [0.1349413459362e-05, 0.2001643230755e+01, 0.1484170571900e-02],
        [0.1170141900476e-05, 0.2424750491620e+01, 0.6327837846670e+00],
        [0.1054355266820e-05, 0.3123311487576e+01, 0.4337116142245e+00],
        [0.9800822461610e-06, 0.3026258088130e+01, 0.1052268489556e+01],
        [0.1091203749931e-05, 0.3157811670347e+01, 0.1162474756779e+01],
        
        [0.6960236715913e-06, 0.8219570542313e+00, 0.1066495398892e+01],
        [0.5689257296909e-06, 0.1323052375236e+01, 0.9491756770005e+00],
        [0.6613172135802e-06, 0.2765348881598e+00, 0.8460828644453e+00],
        [0.6277702517571e-06, 0.5794064466382e+01, 0.1480791608091e+00],
        [0.6304884066699e-06, 0.7323555380787e+00, 0.2243449970715e+00],
        [0.4897850467382e-06, 0.3062464235399e+01, 0.3340612434717e+01],
        [0.3759148598786e-06, 0.4588290469664e+01, 0.3516457698740e-01],
        [0.3110520548195e-06, 0.1374299536572e+01, 0.6373574839730e-01],
        [0.3064708359780e-06, 0.4222267485047e+01, 0.1104591729320e-01],
        [0.2856347168241e-06, 0.3714202944973e+01, 0.1510475019529e+00],
        
        [0.2840945514288e-06, 0.2847972875882e+01, 0.4110125927500e-01],
        [0.2378951599405e-06, 0.3762072563388e+01, 0.2275259891141e+00],
        [0.2714229481417e-06, 0.1036049980031e+01, 0.2535050500000e-01],
        [0.2323551717307e-06, 0.4682388599076e+00, 0.8582758298370e-01],
        [0.1881790512219e-06, 0.4790565425418e+01, 0.2118763888447e+01],
        [0.2261353968371e-06, 0.1669144912212e+01, 0.7181332454670e-01],
        [0.2214546389848e-06, 0.3937717281614e+01, 0.2968341143800e-02],
        [0.2184915594933e-06, 0.1129169845099e+00, 0.7775000683430e-01],
        [0.2000164937936e-06, 0.4030009638488e+01, 0.2093666171530e+00],
        [0.1966105136719e-06, 0.8745955786834e+00, 0.2172315424036e+00],
        
        [0.1904742332624e-06, 0.5919743598964e+01, 0.2022531624851e+00],
        [0.1657399705031e-06, 0.2549141484884e+01, 0.7358765972222e+00],
        [0.1574070533987e-06, 0.5277533020230e+01, 0.7429900518901e+00],
        [0.1832261651039e-06, 0.3064688127777e+01, 0.3235053470014e+00],
        [0.1733615346569e-06, 0.3011432799094e+01, 0.1385174140878e+00],
        [0.1549124014496e-06, 0.4005569132359e+01, 0.5154640627760e+00],
        [0.1637044713838e-06, 0.1831375966632e+01, 0.8531963191132e+00],
        [0.1123420082383e-06, 0.1180270407578e+01, 0.1990721704425e+00],
        [0.1083754165740e-06, 0.3414101320863e+00, 0.5439178814476e+00],
        [0.1156638012655e-06, 0.6130479452594e+00, 0.5257585094865e+00],
        
        [0.1142548785134e-06, 0.3724761948846e+01, 0.5336234347371e+00],
        [0.7921463895965e-07, 0.2435425589361e+01, 0.1478866649112e+01],
        [0.7428600285231e-07, 0.3542144398753e+01, 0.2164800718209e+00],
        [0.8323211246747e-07, 0.3525058072354e+01, 0.1692165728891e+01],
        [0.7257595116312e-07, 0.1364299431982e+01, 0.2101180877357e+00],
        [0.7111185833236e-07, 0.2460478875808e+01, 0.4155522422634e+00],
        [0.6868090383716e-07, 0.4397327670704e+01, 0.1173197218910e+00],
        [0.7226419974175e-07, 0.4042647308905e+01, 0.1265567569334e+01],
        [0.6955642383177e-07, 0.2865047906085e+01, 0.9562891316684e+00],
        [0.7492139296331e-07, 0.5014278994215e+01, 0.1422690933580e-01],
        
        [0.6598363128857e-07, 0.2376730020492e+01, 0.6470106940028e+00],
        [0.7381147293385e-07, 0.3272990384244e+01, 0.1581959461667e+01],
        [0.6402909624032e-07, 0.5302290955138e+01, 0.9597935788730e-01],
        [0.6237454263857e-07, 0.5444144425332e+01, 0.7084920306520e-01],
        [0.5241198544016e-07, 0.4215359579205e+01, 0.5265099800692e+00],
        [0.5144463853918e-07, 0.1218916689916e+00, 0.5328719641544e+00],
        [0.5868164772299e-07, 0.2369402002213e+01, 0.7871412831580e-01],
        [0.6233195669151e-07, 0.1254922242403e+01, 0.2608790314060e+02],
        [0.6068463791422e-07, 0.5679713760431e+01, 0.1114304132498e+00],
        [0.4359361135065e-07, 0.6097219641646e+00, 0.1375773836557e+01],
        
        [0.4686510366826e-07, 0.4786231041431e+01, 0.1143987543936e+00],
        [0.3758977287225e-07, 0.1167368068139e+01, 0.1596186371003e+01],
        [0.4282051974778e-07, 0.1519471064319e+01, 0.2770348281756e+00],
        [0.5153765386113e-07, 0.1860532322984e+01, 0.2228608264996e+00],
        [0.4575129387188e-07, 0.7632857887158e+00, 0.1465949902372e+00],
        [0.3326844933286e-07, 0.1298219485285e+01, 0.5070101000000e-01],
        [0.3748617450984e-07, 0.1046510321062e+01, 0.4903339079539e+00],
        [0.2816756661499e-07, 0.3434522346190e+01, 0.2991266627620e+00],
        [0.3412750405039e-07, 0.2523766270318e+01, 0.3518164938661e+00],
        [0.2655796761776e-07, 0.2904422260194e+01, 0.6256703299991e+00],
        
        [0.2963597929458e-07, 0.5923900431149e+00, 0.1099462426779e+00],
        [0.2539523734781e-07, 0.4851947722567e+01, 0.1256615170089e+02],
        [0.2283087914139e-07, 0.3400498595496e+01, 0.6681224869435e+01],
        [0.2321309799331e-07, 0.5789099148673e+01, 0.3368040641550e-01],
        [0.2549657649750e-07, 0.3991856479792e-01, 0.1169588211447e+01],
        [0.2290462303977e-07, 0.2788567577052e+01, 0.1045155034888e+01],
        [0.1945398522914e-07, 0.3290896998176e+01, 0.1155361302111e+01],
        [0.1849171512638e-07, 0.2698060129367e+01, 0.4452511715700e-02],
        [0.1647199834254e-07, 0.3016735644085e+01, 0.4408250688924e+00],
        [0.1529530765273e-07, 0.5573043116178e+01, 0.6521991896920e-01],
        
        [0.1433199339978e-07, 0.1481192356147e+01, 0.9420622223326e+00],
        [0.1729134193602e-07, 0.1422817538933e+01, 0.2108507877249e+00],
        [0.1716463931346e-07, 0.3469468901855e+01, 0.2157473718317e+00],
        [0.1391206061378e-07, 0.6122436220547e+01, 0.4123712502208e+00],
        [0.1404746661924e-07, 0.1647765641936e+01, 0.4258542984690e-01],
        [0.1410452399455e-07, 0.5989729161964e+01, 0.2258291676434e+00],
        [0.1089828772168e-07, 0.2833705509371e+01, 0.4226656969313e+00],
        [0.1047374564948e-07, 0.5090690007331e+00, 0.3092784376656e+00],
        [0.1358279126532e-07, 0.5128990262836e+01, 0.7923417740620e-01],
        [0.1020456476148e-07, 0.9632772880808e+00, 0.1456308687557e+00],
        
        [0.1033428735328e-07, 0.3223779318418e+01, 0.1795258541446e+01],
        [0.1412435841540e-07, 0.2410271572721e+01, 0.1525316725248e+00],
        [0.9722759371574e-08, 0.2333531395690e+01, 0.8434341241180e-01],
        [0.9657334084704e-08, 0.6199270974168e+01, 0.1272681024002e+01],
        [0.1083641148690e-07, 0.2864222292929e+01, 0.7032915397480e-01],
        [0.1067318403838e-07, 0.5833458866568e+00, 0.2123349582968e+00],
        [0.1062366201976e-07, 0.4307753989494e+01, 0.2142632012598e+00],
        [0.1236364149266e-07, 0.2873917870593e+01, 0.1847279083684e+00],
        [0.1092759489593e-07, 0.2959887266733e+01, 0.1370332435159e+00],
        [0.8912069362899e-08, 0.5141213702562e+01, 0.2648454860559e+01],
        
        [0.9656467707970e-08, 0.4532182462323e+01, 0.4376440768498e+00],
        [0.8098386150135e-08, 0.2268906338379e+01, 0.2880807454688e+00],
        [0.7857714675000e-08, 0.4055544260745e+01, 0.2037373330570e+00],
        [0.7288455940646e-08, 0.5357901655142e+01, 0.1129145838217e+00],
        [0.9450595950552e-08, 0.4264926963939e+01, 0.5272426800584e+00],
        [0.9381718247537e-08, 0.7489366976576e-01, 0.5321392641652e+00],
        [0.7079052646038e-08, 0.1923311052874e+01, 0.6288513220417e+00],
        [0.9259004415344e-08, 0.2970256853438e+01, 0.1606092486742e+00],
        [0.8259801499742e-08, 0.3327056314697e+01, 0.8389694097774e+00],
        [0.6476334355779e-08, 0.2954925505727e+01, 0.2008557621224e+01],
        
        [0.5984021492007e-08, 0.9138753105829e+00, 0.2042657109477e+02],
        [0.5989546863181e-08, 0.3244464082031e+01, 0.2111650433779e+01],
        [0.6233108606023e-08, 0.4995232638403e+00, 0.4305306221819e+00],
        [0.6877299149965e-08, 0.2834987233449e+01, 0.9561746721300e-02],
        [0.8311234227190e-08, 0.2202951835758e+01, 0.3801276407308e+00],
        [0.6599472832414e-08, 0.4478581462618e+01, 0.1063314406849e+01],
        [0.6160491096549e-08, 0.5145858696411e+01, 0.1368660381889e+01],
        [0.6164772043891e-08, 0.3762976697911e+00, 0.4234171675140e+00],
        [0.6363248684450e-08, 0.3162246718685e+01, 0.1253008786510e-01],
        [0.6448587520999e-08, 0.3442693302119e+01, 0.5287268506303e+00],
        
        [0.6431662283977e-08, 0.8977549136606e+00, 0.5306550935933e+00],
        [0.6351223158474e-08, 0.4306447410369e+01, 0.5217580628120e+02],
        [0.5476721393451e-08, 0.3888529177855e+01, 0.2221856701002e+01],
        [0.5341772572619e-08, 0.2655560662512e+01, 0.7466759693650e-01],
        [0.5337055758302e-08, 0.5164990735946e+01, 0.7489573444450e-01],
        [0.5373120816787e-08, 0.6041214553456e+01, 0.1274714967946e+00],
        [0.5392351705426e-08, 0.9177763485932e+00, 0.1055449481598e+01],
        [0.6688495850205e-08, 0.3089608126937e+01, 0.2213766559277e+00],
        [0.5072003660362e-08, 0.4311316541553e+01, 0.2132517061319e+00],
        [0.5070726650455e-08, 0.5790675464444e+00, 0.2133464534247e+00],
        
        [0.5658012950032e-08, 0.2703945510675e+01, 0.7287631425543e+00],
        [0.4835509924854e-08, 0.2975422976065e+01, 0.7160067364790e-01],
        [0.6479821978012e-08, 0.1324168733114e+01, 0.2209183458640e-01],
        [0.6230636494980e-08, 0.2860103632836e+01, 0.3306188016693e+00],
        [0.4649239516213e-08, 0.4832259763403e+01, 0.7796265773310e-01],
        [0.6487325792700e-08, 0.2726165825042e+01, 0.3884652414254e+00],
        [0.4682823682770e-08, 0.6966602455408e+00, 0.1073608853559e+01],
        [0.5704230804976e-08, 0.5669634104606e+01, 0.8731175355560e-01],
        [0.6125413585489e-08, 0.1513386538915e+01, 0.7605151500000e-01],
        [0.6035825038187e-08, 0.1983509168227e+01, 0.9846002785331e+00],
        
        [0.4331123462303e-08, 0.2782892992807e+01, 0.4297791515992e+00],
        [0.4681107685143e-08, 0.5337232886836e+01, 0.2127790306879e+00],
        [0.4669105829655e-08, 0.5837133792160e+01, 0.2138191288687e+00],
        [0.5138823602365e-08, 0.3080560200507e+01, 0.7233337363710e-01],
        [0.4615856664534e-08, 0.1661747897471e+01, 0.8603097737811e+00],
        [0.4496916702197e-08, 0.2112508027068e+01, 0.7381754420900e-01],
        [0.4278479042945e-08, 0.5716528462627e+01, 0.7574578717200e-01],
        [0.3840525503932e-08, 0.6424172726492e+00, 0.3407705765729e+00],
        [0.4866636509685e-08, 0.4919244697715e+01, 0.7722995774390e-01],
        [0.3526100639296e-08, 0.2550821052734e+01, 0.6225157782540e-01],
        
        [0.3939558488075e-08, 0.3939331491710e+01, 0.5268983110410e-01],
        [0.4041268772576e-08, 0.2275337571218e+01, 0.3503323232942e+00],
        [0.3948761842853e-08, 0.1999324200790e+01, 0.1451108196653e+00],
        [0.3258394550029e-08, 0.9121001378200e+00, 0.5296435984654e+00],
        [0.3257897048761e-08, 0.3428428660869e+01, 0.5297383457582e+00],
        [0.3842559031298e-08, 0.6132927720035e+01, 0.9098186128426e+00],
        [0.3109920095448e-08, 0.7693650193003e+00, 0.3932462625300e-02],
        [0.3132237775119e-08, 0.3621293854908e+01, 0.2346394437820e+00],
        [0.3942189421510e-08, 0.4841863659733e+01, 0.3180992042600e-02],
        [0.3796972285340e-08, 0.1814174994268e+01, 0.1862120789403e+00],
        
        [0.3995640233688e-08, 0.1386990406091e+01, 0.4549093064213e+00],
        [0.2875013727414e-08, 0.9178318587177e+00, 0.1905464808669e+01],
        [0.3073719932844e-08, 0.2688923811835e+01, 0.3628624111593e+00],
        [0.2731016580075e-08, 0.1188259127584e+01, 0.2131850110243e+00],
        [0.2729549896546e-08, 0.3702160634273e+01, 0.2134131485323e+00],
        [0.3339372892449e-08, 0.7199163960331e+00, 0.2007689919132e+00],
        [0.2898833764204e-08, 0.1916709364999e+01, 0.5291709230214e+00],
        [0.2894536549362e-08, 0.2424043195547e+01, 0.5302110212022e+00],
        [0.3096872473843e-08, 0.4445894977497e+01, 0.2976424921901e+00],
        [0.2635672326810e-08, 0.3814366984117e+01, 0.1485980103780e+01],
        
        [0.3649302697001e-08, 0.2924200596084e+01, 0.6044726378023e+00],
        [0.3127954585895e-08, 0.1842251648327e+01, 0.1084620721060e+00],
        [0.2616040173947e-08, 0.4155841921984e+01, 0.1258454114666e+01],
        [0.2597395859860e-08, 0.1158045978874e+00, 0.2103781122809e+00],
        [0.2593286172210e-08, 0.4771850408691e+01, 0.2162200472757e+00],
        [0.2481823585747e-08, 0.4608842558889e+00, 0.1062562936266e+01],
        [0.2742219550725e-08, 0.1538781127028e+01, 0.5651155736444e+00],
        [0.3199558469610e-08, 0.3226647822878e+00, 0.7036329877322e+00],
        [0.2666088542957e-08, 0.1967991731219e+00, 0.1400015846597e+00],
        [0.2397067430580e-08, 0.3707036669873e+01, 0.2125476091956e+00],
        
        [0.2376570772738e-08, 0.1182086628042e+01, 0.2140505503610e+00],
        [0.2547228007887e-08, 0.4906256820629e+01, 0.1534957940063e+00],
        [0.2265575594114e-08, 0.3414949866857e+01, 0.2235935264888e+00],
        [0.2464381430585e-08, 0.4599122275378e+01, 0.2091065926078e+00],
        [0.2433408527044e-08, 0.2830751145445e+00, 0.2174915669488e+00],
        [0.2443605509076e-08, 0.4212046432538e+01, 0.1739420156204e+00],
        [0.2319779262465e-08, 0.9881978408630e+00, 0.7530171478090e-01],
        [0.2284622835465e-08, 0.5565347331588e+00, 0.7426161660010e-01],
        [0.2467268750783e-08, 0.5655708150766e+00, 0.2526561439362e+00],
        [0.2808513492782e-08, 0.1418405053408e+01, 0.5636314030725e+00],
        
        [0.2329528932532e-08, 0.4069557545675e+01, 0.1056200952181e+01],
        [0.9698639532817e-09, 0.1074134313634e+01, 0.7826370942180e+02] 
        
    ])

    # SSB-to-Sun, T^0, Y
    s0y = np.array([
        [0.4955392320126e-02, 0.2170467313679e+01, 0.5296909721118e+00],
        [0.2722325167392e-02, 0.2444433682196e+01, 0.2132990797783e+00],
        [0.1546579925346e-02, 0.5992779281546e+00, 0.3813291813120e-01],
        [0.8363140252966e-03, 0.7687356310801e+00, 0.7478166569050e-01],
        [0.3385792683603e-03, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.1201192221613e-03, 0.2520035601514e+01, 0.1059381944224e+01],
        [0.7587125720554e-04, 0.1669954006449e+01, 0.4265981595566e+00],
        [0.1964155361250e-04, 0.5707743963343e+01, 0.2061856251104e+00],
        [0.1891900364909e-04, 0.2320960679937e+01, 0.2204125344462e+00],
        [0.1937373433356e-04, 0.3226940689555e+01, 0.1495633313810e+00],
        
        [0.1437139941351e-04, 0.2301626908096e+01, 0.5225775174439e+00],
        [0.1406267683099e-04, 0.5188579265542e+01, 0.5368044267797e+00],
        [0.1178703080346e-04, 0.5489483248476e+01, 0.7626583626240e-01],
        [0.8079835186041e-05, 0.1683751835264e+01, 0.3664874755930e-01],
        [0.7623253594652e-05, 0.2656400462961e+01, 0.3961708870310e-01],
        [0.6248667483971e-05, 0.4992775362055e+01, 0.7329749511860e-01],
        [0.4366353695038e-05, 0.2869706279678e+01, 0.1589072916335e+01],
        [0.3829101568895e-05, 0.3572131359950e+01, 0.7113454667900e-02],
        [0.3175733773908e-05, 0.4535372530045e+01, 0.4194847048887e+00],
        [0.3092437902159e-05, 0.9230153317909e+00, 0.6398972393349e+00],
        
        [0.2874168812154e-05, 0.3363143761101e+01, 0.1102062672231e+00],
        [0.3040119321826e-05, 0.3324250895675e+01, 0.6283075850446e+01],
        [0.2699723308006e-05, 0.2917882441928e+00, 0.1030928125552e+00],
        [0.2134832683534e-05, 0.4220997202487e+01, 0.3163918923335e+00],
        [0.1770412139433e-05, 0.4747318496462e+01, 0.1021328554739e+02],
        [0.1377264209373e-05, 0.4305058462401e+00, 0.1484170571900e-02],
        [0.1127814538960e-05, 0.8538177240740e+00, 0.6327837846670e+00],
        [0.1055608090130e-05, 0.1551800742580e+01, 0.4337116142245e+00],
        [0.9802673861420e-06, 0.1459646735377e+01, 0.1052268489556e+01],
        [0.1090329461951e-05, 0.1587351228711e+01, 0.1162474756779e+01],
        
        [0.6959590025090e-06, 0.5534442628766e+01, 0.1066495398892e+01],
        [0.5664914529542e-06, 0.6030673003297e+01, 0.9491756770005e+00],
        [0.6607787763599e-06, 0.4989507233927e+01, 0.8460828644453e+00],
        [0.6269725742838e-06, 0.4222951804572e+01, 0.1480791608091e+00],
        [0.6301889697863e-06, 0.5444316669126e+01, 0.2243449970715e+00],
        [0.4891042662861e-06, 0.1490552839784e+01, 0.3340612434717e+01],
        [0.3457083123290e-06, 0.3030475486049e+01, 0.3516457698740e-01],
        [0.3032559967314e-06, 0.2652038793632e+01, 0.1104591729320e-01],
        [0.2841133988903e-06, 0.1276744786829e+01, 0.4110125927500e-01],
        [0.2855564444432e-06, 0.2143368674733e+01, 0.1510475019529e+00],
        
        [0.2765157135038e-06, 0.5444186109077e+01, 0.6373574839730e-01],
        [0.2382312465034e-06, 0.2190521137593e+01, 0.2275259891141e+00],
        [0.2808060365077e-06, 0.5735195064841e+01, 0.2535050500000e-01],
        [0.2332175234405e-06, 0.9481985524859e-01, 0.7181332454670e-01],
        [0.2322488199659e-06, 0.5180499361533e+01, 0.8582758298370e-01],
        [0.1881850258423e-06, 0.3219788273885e+01, 0.2118763888447e+01],
        [0.2196111392808e-06, 0.2366941159761e+01, 0.2968341143800e-02],
        [0.2183810335519e-06, 0.4825445110915e+01, 0.7775000683430e-01],
        [0.2002733093326e-06, 0.2457148995307e+01, 0.2093666171530e+00],
        [0.1967111767229e-06, 0.5586291545459e+01, 0.2172315424036e+00],
        
        [0.1568473250543e-06, 0.3708003123320e+01, 0.7429900518901e+00],
        [0.1852528314300e-06, 0.4310638151560e+01, 0.2022531624851e+00],
        [0.1832111226447e-06, 0.1494665322656e+01, 0.3235053470014e+00],
        [0.1746805502310e-06, 0.1451378500784e+01, 0.1385174140878e+00],
        [0.1555730966650e-06, 0.1068040418198e+01, 0.7358765972222e+00],
        [0.1554883462559e-06, 0.2442579035461e+01, 0.5154640627760e+00],
        [0.1638380568746e-06, 0.2597913420625e+00, 0.8531963191132e+00],
        [0.1159938593640e-06, 0.5834512021280e+01, 0.1990721704425e+00],
        [0.1083427965695e-06, 0.5054033177950e+01, 0.5439178814476e+00],
        [0.1156480369431e-06, 0.5325677432457e+01, 0.5257585094865e+00],
        
        [0.1141308860095e-06, 0.2153403923857e+01, 0.5336234347371e+00],
        [0.7913146470946e-07, 0.8642846847027e+00, 0.1478866649112e+01],
        [0.7439752463733e-07, 0.1970628496213e+01, 0.2164800718209e+00],
        [0.7280277104079e-07, 0.6073307250609e+01, 0.2101180877357e+00],
        [0.8319567719136e-07, 0.1954371928334e+01, 0.1692165728891e+01],
        [0.7137705549290e-07, 0.8904989440909e+00, 0.4155522422634e+00],
        [0.6900825396225e-07, 0.2825717714977e+01, 0.1173197218910e+00],
        [0.7245757216635e-07, 0.2481677513331e+01, 0.1265567569334e+01],
        [0.6961165696255e-07, 0.1292955312978e+01, 0.9562891316684e+00],
        [0.7571804456890e-07, 0.3427517575069e+01, 0.1422690933580e-01],
        
        [0.6605425721904e-07, 0.8052192701492e+00, 0.6470106940028e+00],
        [0.7375477357248e-07, 0.1705076390088e+01, 0.1581959461667e+01],
        [0.7041664951470e-07, 0.4848356967891e+00, 0.9597935788730e-01],
        [0.6322199535763e-07, 0.3878069473909e+01, 0.7084920306520e-01],
        [0.5244380279191e-07, 0.2645560544125e+01, 0.5265099800692e+00],
        [0.5143125704988e-07, 0.4834486101370e+01, 0.5328719641544e+00],
        [0.5871866319373e-07, 0.7981472548900e+00, 0.7871412831580e-01],
        [0.6300822573871e-07, 0.5979398788281e+01, 0.2608790314060e+02],
        [0.6062154271548e-07, 0.4108655402756e+01, 0.1114304132498e+00],
        [0.4361912339976e-07, 0.5322624319280e+01, 0.1375773836557e+01],
        
        [0.4417005920067e-07, 0.6240817359284e+01, 0.2770348281756e+00],
        [0.4686806749936e-07, 0.3214977301156e+01, 0.1143987543936e+00],
        [0.3758892132305e-07, 0.5879809634765e+01, 0.1596186371003e+01],
        [0.5151351332319e-07, 0.2893377688007e+00, 0.2228608264996e+00],
        [0.4554683578572e-07, 0.5475427144122e+01, 0.1465949902372e+00],
        [0.3442381385338e-07, 0.5992034796640e+01, 0.5070101000000e-01],
        [0.2831093954933e-07, 0.5367350273914e+01, 0.3092784376656e+00],
        [0.3756267090084e-07, 0.5758171285420e+01, 0.4903339079539e+00],
        [0.2816374679892e-07, 0.1863718700923e+01, 0.2991266627620e+00],
        [0.3419307025569e-07, 0.9524347534130e+00, 0.3518164938661e+00],
        
        [0.2904250494239e-07, 0.5304471615602e+01, 0.1099462426779e+00],
        [0.2471734511206e-07, 0.1297069793530e+01, 0.6256703299991e+00],
        [0.2539620831872e-07, 0.3281126083375e+01, 0.1256615170089e+02],
        [0.2281017868007e-07, 0.1829122133165e+01, 0.6681224869435e+01],
        [0.2275319473335e-07, 0.5797198160181e+01, 0.3932462625300e-02],
        [0.2547755368442e-07, 0.4752697708330e+01, 0.1169588211447e+01],
        [0.2285979669317e-07, 0.1223205292886e+01, 0.1045155034888e+01],
        [0.1913386560994e-07, 0.1757532993389e+01, 0.1155361302111e+01],
        [0.1809020525147e-07, 0.4246116108791e+01, 0.3368040641550e-01],
        [0.1649213300201e-07, 0.1445162890627e+01, 0.4408250688924e+00],
        
        [0.1834972793932e-07, 0.1126917567225e+01, 0.4452511715700e-02],
        [0.1439550648138e-07, 0.6160756834764e+01, 0.9420622223326e+00],
        [0.1487645457041e-07, 0.4358761931792e+01, 0.4123712502208e+00],
        [0.1731729516660e-07, 0.6134456753344e+01, 0.2108507877249e+00],
        [0.1717747163567e-07, 0.1898186084455e+01, 0.2157473718317e+00],
        [0.1418190430374e-07, 0.4180286741266e+01, 0.6521991896920e-01],
        [0.1404844134873e-07, 0.7654053565412e-01, 0.4258542984690e-01],
        [0.1409842846538e-07, 0.4418612420312e+01, 0.2258291676434e+00],
        [0.1090948346291e-07, 0.1260615686131e+01, 0.4226656969313e+00],
        [0.1357577323612e-07, 0.3558248818690e+01, 0.7923417740620e-01],
        
        [0.1018154061960e-07, 0.5676087241256e+01, 0.1456308687557e+00],
        [0.1412073972109e-07, 0.8394392632422e+00, 0.1525316725248e+00],
        [0.1030938326496e-07, 0.1653593274064e+01, 0.1795258541446e+01],
        [0.1180081567104e-07, 0.1285802592036e+01, 0.7032915397480e-01],
        [0.9708510575650e-08, 0.7631889488106e+00, 0.8434341241180e-01],
        [0.9637689663447e-08, 0.4630642649176e+01, 0.1272681024002e+01],
        [0.1068910429389e-07, 0.5294934032165e+01, 0.2123349582968e+00],
        [0.1063716179336e-07, 0.2736266800832e+01, 0.2142632012598e+00],
        [0.1234858713814e-07, 0.1302891146570e+01, 0.1847279083684e+00],
        [0.8912631189738e-08, 0.3570415993621e+01, 0.2648454860559e+01],
        
        [0.1036378285534e-07, 0.4236693440949e+01, 0.1370332435159e+00],
        [0.9667798501561e-08, 0.2960768892398e+01, 0.4376440768498e+00],
        [0.8108314201902e-08, 0.6987781646841e+00, 0.2880807454688e+00],
        [0.7648364324628e-08, 0.2499017863863e+01, 0.2037373330570e+00],
        [0.7286136828406e-08, 0.3787426951665e+01, 0.1129145838217e+00],
        [0.9448237743913e-08, 0.2694354332983e+01, 0.5272426800584e+00],
        [0.9374276106428e-08, 0.4787121277064e+01, 0.5321392641652e+00],
        [0.7100226287462e-08, 0.3530238792101e+00, 0.6288513220417e+00],
        [0.9253056659571e-08, 0.1399478925664e+01, 0.1606092486742e+00],
        [0.6636432145504e-08, 0.3479575438447e+01, 0.1368660381889e+01],
        
        [0.6469975312932e-08, 0.1383669964800e+01, 0.2008557621224e+01],
        [0.7335849729765e-08, 0.1243698166898e+01, 0.9561746721300e-02],
        [0.8743421205855e-08, 0.3776164289301e+01, 0.3801276407308e+00],
        [0.5993635744494e-08, 0.5627122113596e+01, 0.2042657109477e+02],
        [0.5981008479693e-08, 0.1674336636752e+01, 0.2111650433779e+01],
        [0.6188535145838e-08, 0.5214925208672e+01, 0.4305306221819e+00],
        [0.6596074017566e-08, 0.2907653268124e+01, 0.1063314406849e+01],
        [0.6630815126226e-08, 0.2127643669658e+01, 0.8389694097774e+00],
        [0.6156772830040e-08, 0.5082160803295e+01, 0.4234171675140e+00],
        [0.6446960563014e-08, 0.1872100916905e+01, 0.5287268506303e+00],
        
        [0.6429324424668e-08, 0.5610276103577e+01, 0.5306550935933e+00],
        [0.6302232396465e-08, 0.1592152049607e+01, 0.1253008786510e-01],
        [0.6399244436159e-08, 0.2746214421532e+01, 0.5217580628120e+02],
        [0.5474965172558e-08, 0.2317666374383e+01, 0.2221856701002e+01],
        [0.5339293190692e-08, 0.1084724961156e+01, 0.7466759693650e-01],
        [0.5334733683389e-08, 0.3594106067745e+01, 0.7489573444450e-01],
        [0.5392665782110e-08, 0.5630254365606e+01, 0.1055449481598e+01],
        [0.6682075673789e-08, 0.1518480041732e+01, 0.2213766559277e+00],
        [0.5079130495960e-08, 0.2739765115711e+01, 0.2132517061319e+00],
        [0.5077759793261e-08, 0.5290711290094e+01, 0.2133464534247e+00],
        
        [0.4832037368310e-08, 0.1404473217200e+01, 0.7160067364790e-01],
        [0.6463279674802e-08, 0.6038381695210e+01, 0.2209183458640e-01],
        [0.6240592771560e-08, 0.1290170653666e+01, 0.3306188016693e+00],
        [0.4672013521493e-08, 0.3261895939677e+01, 0.7796265773310e-01],
        [0.6500650750348e-08, 0.1154522312095e+01, 0.3884652414254e+00],
        [0.6344161389053e-08, 0.6206111545062e+01, 0.7605151500000e-01],
        [0.4682518370646e-08, 0.5409118796685e+01, 0.1073608853559e+01],
        [0.5329460015591e-08, 0.1202985784864e+01, 0.7287631425543e+00],
        [0.5701588675898e-08, 0.4098715257064e+01, 0.8731175355560e-01],
        [0.6030690867211e-08, 0.4132033218460e+00, 0.9846002785331e+00],
        
        [0.4336256312655e-08, 0.1211415991827e+01, 0.4297791515992e+00],
        [0.4688498808975e-08, 0.3765479072409e+01, 0.2127790306879e+00],
        [0.4675578609335e-08, 0.4265540037226e+01, 0.2138191288687e+00],
        [0.4225578112158e-08, 0.5237566010676e+01, 0.3407705765729e+00],
        [0.5139422230028e-08, 0.1507173079513e+01, 0.7233337363710e-01],
        [0.4619995093571e-08, 0.9023957449848e-01, 0.8603097737811e+00],
        [0.4494776255461e-08, 0.5414930552139e+00, 0.7381754420900e-01],
        [0.4274026276788e-08, 0.4145735303659e+01, 0.7574578717200e-01],
        [0.5018141789353e-08, 0.3344408829055e+01, 0.3180992042600e-02],
        [0.4866163952181e-08, 0.3348534657607e+01, 0.7722995774390e-01],
        
        [0.4111986020501e-08, 0.4198823597220e+00, 0.1451108196653e+00],
        [0.3356142784950e-08, 0.5609144747180e+01, 0.1274714967946e+00],
        [0.4070575554551e-08, 0.7028411059224e+00, 0.3503323232942e+00],
        [0.3257451857278e-08, 0.5624697983086e+01, 0.5296435984654e+00],
        [0.3256973703026e-08, 0.1857842076707e+01, 0.5297383457582e+00],
        [0.3830771508640e-08, 0.4562887279931e+01, 0.9098186128426e+00],
        [0.3725024005962e-08, 0.2358058692652e+00, 0.1084620721060e+00],
        [0.3136763921756e-08, 0.2049731526845e+01, 0.2346394437820e+00],
        [0.3795147256194e-08, 0.2432356296933e+00, 0.1862120789403e+00],
        [0.2877342229911e-08, 0.5631101279387e+01, 0.1905464808669e+01],
        
        [0.3076931798805e-08, 0.1117615737392e+01, 0.3628624111593e+00],
        [0.2734765945273e-08, 0.5899826516955e+01, 0.2131850110243e+00],
        [0.2733405296885e-08, 0.2130562964070e+01, 0.2134131485323e+00],
        [0.2898552353410e-08, 0.3462387048225e+00, 0.5291709230214e+00],
        [0.2893736103681e-08, 0.8534352781543e+00, 0.5302110212022e+00],
        [0.3095717734137e-08, 0.2875061429041e+01, 0.2976424921901e+00],
        [0.2636190425832e-08, 0.2242512846659e+01, 0.1485980103780e+01],
        [0.3645512095537e-08, 0.1354016903958e+01, 0.6044726378023e+00],
        [0.2808173547723e-08, 0.6705114365631e-01, 0.6225157782540e-01],
        [0.2625012866888e-08, 0.4775705748482e+01, 0.5268983110410e-01],
        
        [0.2572233995651e-08, 0.2638924216139e+01, 0.1258454114666e+01],
        [0.2604238824792e-08, 0.4826358927373e+01, 0.2103781122809e+00],
        [0.2596886385239e-08, 0.3200388483118e+01, 0.2162200472757e+00],
        [0.3228057304264e-08, 0.5384848409563e+01, 0.2007689919132e+00],
        [0.2481601798252e-08, 0.5173373487744e+01, 0.1062562936266e+01],
        [0.2745977498864e-08, 0.6250966149853e+01, 0.5651155736444e+00],
        [0.2669878833811e-08, 0.4906001352499e+01, 0.1400015846597e+00],
        [0.3203986611711e-08, 0.5034333010005e+01, 0.7036329877322e+00],
        [0.3354961227212e-08, 0.6108262423137e+01, 0.4549093064213e+00],
        [0.2400407324558e-08, 0.2135399294955e+01, 0.2125476091956e+00],
        
        [0.2379905859802e-08, 0.5893721933961e+01, 0.2140505503610e+00],
        [0.2550844302187e-08, 0.3331940762063e+01, 0.1534957940063e+00],
        [0.2268824211001e-08, 0.1843418461035e+01, 0.2235935264888e+00],
        [0.2464700891204e-08, 0.3029548547230e+01, 0.2091065926078e+00],
        [0.2436814726024e-08, 0.4994717970364e+01, 0.2174915669488e+00],
        [0.2443623894745e-08, 0.2645102591375e+01, 0.1739420156204e+00],
        [0.2318701783838e-08, 0.5700547397897e+01, 0.7530171478090e-01],
        [0.2284448700256e-08, 0.5268898905872e+01, 0.7426161660010e-01],
        [0.2468848123510e-08, 0.5276280575078e+01, 0.2526561439362e+00],
        [0.2814052350303e-08, 0.6130168623475e+01, 0.5636314030725e+00],
        
        [0.2243662755220e-08, 0.6631692457995e+00, 0.8886590321940e-01],
        [0.2330795855941e-08, 0.2499435487702e+01, 0.1056200952181e+01],
        [0.9757679038404e-09, 0.5796846023126e+01, 0.7826370942180e+02] 
    ])

    # SSB-to-Sun, T^0, Z
    s0z = np.array([
        [0.1181255122986e-03, 0.4607918989164e+00, 0.2132990797783e+00],
        [0.1127777651095e-03, 0.4169146331296e+00, 0.5296909721118e+00],
        [0.4777754401806e-04, 0.4582657007130e+01, 0.3813291813120e-01],
        [0.1129354285772e-04, 0.5758735142480e+01, 0.7478166569050e-01],
        [-0.1149543637123e-04,0.0000000000000e+00, 0.0000000000000e+00],
        [0.3298730512306e-05, 0.5978801994625e+01, 0.4265981595566e+00],
        [0.2733376706079e-05, 0.7665413691040e+00, 0.1059381944224e+01],
        [0.9426389657270e-06, 0.3710201265838e+01, 0.2061856251104e+00],
        [0.8187517749552e-06, 0.3390675605802e+00, 0.2204125344462e+00],
        [0.4080447871819e-06, 0.4552296640088e+00, 0.5225775174439e+00],
        
        [0.3169973017028e-06, 0.3445455899321e+01, 0.5368044267797e+00],
        [0.2438098615549e-06, 0.5664675150648e+01, 0.3664874755930e-01],
        [0.2601897517235e-06, 0.1931894095697e+01, 0.1495633313810e+00],
        [0.2314558080079e-06, 0.3666319115574e+00, 0.3961708870310e-01],
        [0.1962549548002e-06, 0.3167411699020e+01, 0.7626583626240e-01],
        [0.2180518287925e-06, 0.1544420746580e+01, 0.7113454667900e-02],
        [0.1451382442868e-06, 0.1583756740070e+01, 0.1102062672231e+00],
        [0.1358439007389e-06, 0.5239941758280e+01, 0.6398972393349e+00],
        [0.1050585898028e-06, 0.2266958352859e+01, 0.3163918923335e+00],
        [0.1050029870186e-06, 0.2711495250354e+01, 0.4194847048887e+00],
        
        [0.9934920679800e-07, 0.1116208151396e+01, 0.1589072916335e+01],
        [0.1048395331560e-06, 0.3408619600206e+01, 0.1021328554739e+02],
        [0.8370147196668e-07, 0.3810459401087e+01, 0.2535050500000e-01],
        [0.7989856510998e-07, 0.3769910473647e+01, 0.7329749511860e-01],
        [0.5441221655233e-07, 0.2416994903374e+01, 0.1030928125552e+00],
        [0.4610812906784e-07, 0.5858503336994e+01, 0.4337116142245e+00],
        [0.3923022803444e-07, 0.3354170010125e+00, 0.1484170571900e-02],
        [0.2610725582128e-07, 0.5410600646324e+01, 0.6327837846670e+00],
        [0.2455279767721e-07, 0.6120216681403e+01, 0.1162474756779e+01],
        [0.2375530706525e-07, 0.6055443426143e+01, 0.1052268489556e+01],
        
        [0.1782967577553e-07, 0.3146108708004e+01, 0.8460828644453e+00],
        [0.1581687095238e-07, 0.6255496089819e+00, 0.3340612434717e+01],
        [0.1594657672461e-07, 0.3782604300261e+01, 0.1066495398892e+01],
        [0.1563448615040e-07, 0.1997775733196e+01, 0.2022531624851e+00],
        [0.1463624258525e-07, 0.1736316792088e+00, 0.3516457698740e-01],
        [0.1331585056673e-07, 0.4331941830747e+01, 0.9491756770005e+00],
        [0.1130634557637e-07, 0.6152017751825e+01, 0.2968341143800e-02],
        [0.1028949607145e-07, 0.2101792614637e+00, 0.2275259891141e+00],
        [0.1024074971618e-07, 0.4071833211074e+01, 0.5070101000000e-01],
        [0.8826956060303e-08, 0.4861633688145e+00, 0.2093666171530e+00],
        
        [0.8572230171541e-08, 0.5268190724302e+01, 0.4110125927500e-01],
        [0.7649332643544e-08, 0.5134543417106e+01, 0.2608790314060e+02],
        [0.8581673291033e-08, 0.2920218146681e+01, 0.1480791608091e+00],
        [0.8430589300938e-08, 0.3604576619108e+01, 0.2172315424036e+00],
        [0.7776165501012e-08, 0.3772942249792e+01, 0.6373574839730e-01],
        [0.8311070234408e-08, 0.6200412329888e+01, 0.3235053470014e+00],
        [0.6927365212582e-08, 0.4543353113437e+01, 0.8531963191132e+00],
        [0.6791574208598e-08, 0.2882188406238e+01, 0.7181332454670e-01],
        [0.5593100811839e-08, 0.1776646892780e+01, 0.7429900518901e+00],
        [0.4553381853021e-08, 0.3949617611240e+01, 0.7775000683430e-01],
        
        [0.5758000450068e-08, 0.3859251775075e+01, 0.1990721704425e+00],
        [0.4281283457133e-08, 0.1466294631206e+01, 0.2118763888447e+01],
        [0.4206935661097e-08, 0.5421776011706e+01, 0.1104591729320e-01],
        [0.4213751641837e-08, 0.3412048993322e+01, 0.2243449970715e+00],
        [0.5310506239878e-08, 0.5421641370995e+00, 0.5154640627760e+00],
        [0.3827450341320e-08, 0.8887314524995e+00, 0.1510475019529e+00],
        [0.4292435241187e-08, 0.1405043757194e+01, 0.1422690933580e-01],
        [0.3189780702289e-08, 0.1060049293445e+01, 0.1173197218910e+00],
        [0.3226611928069e-08, 0.6270858897442e+01, 0.2164800718209e+00],
        [0.2893897608830e-08, 0.5117563223301e+01, 0.6470106940028e+00],
        
        [0.3239852024578e-08, 0.4079092237983e+01, 0.2101180877357e+00],
        [0.2956892222200e-08, 0.1594917021704e+01, 0.3092784376656e+00],
        [0.2980177912437e-08, 0.5258787667564e+01, 0.4155522422634e+00],
        [0.3163725690776e-08, 0.3854589225479e+01, 0.8582758298370e-01],
        [0.2662262399118e-08, 0.3561326430187e+01, 0.5257585094865e+00],
        [0.2766689135729e-08, 0.3180732086830e+00, 0.1385174140878e+00],
        [0.2411600278464e-08, 0.3324798335058e+01, 0.5439178814476e+00],
        [0.2483527695131e-08, 0.4169069291947e+00, 0.5336234347371e+00],
        [0.7788777276590e-09, 0.1900569908215e+01, 0.5217580628120e+02]  
    ])

    # SSB-to-Sun, T^1, X
    s1x = np.array([
        [-0.1296310361520e-07, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.8975769009438e-08, 0.1128891609250e+01, 0.4265981595566e+00],
        [0.7771113441307e-08, 0.2706039877077e+01, 0.2061856251104e+00],
        [0.7538303866642e-08, 0.2191281289498e+01, 0.2204125344462e+00],
        [0.6061384579336e-08, 0.3248167319958e+01, 0.1059381944224e+01],
        [0.5726994235594e-08, 0.5569981398610e+01, 0.5225775174439e+00],
        [0.5616492836424e-08, 0.5057386614909e+01, 0.5368044267797e+00],
        [0.1010881584769e-08, 0.3473577116095e+01, 0.7113454667900e-02],
        [0.7259606157626e-09, 0.3651858593665e+00, 0.6398972393349e+00],
        [0.8755095026935e-09, 0.1662835408338e+01, 0.4194847048887e+00],
        
        [0.5370491182812e-09, 0.1327673878077e+01, 0.4337116142245e+00],
        [0.5743773887665e-09, 0.4250200846687e+01, 0.2132990797783e+00],
        [0.4408103140300e-09, 0.3598752574277e+01, 0.1589072916335e+01],
        [0.3101892374445e-09, 0.4887822983319e+01, 0.1052268489556e+01],
        [0.3209453713578e-09, 0.9702272295114e+00, 0.5296909721118e+00],
        [0.3017228286064e-09, 0.5484462275949e+01, 0.1066495398892e+01],
        [0.3200700038601e-09, 0.2846613338643e+01, 0.1495633313810e+00],
        [0.2137637279911e-09, 0.5692163292729e+00, 0.3163918923335e+00],
        [0.1899686386727e-09, 0.2061077157189e+01, 0.2275259891141e+00],
        [0.1401994545308e-09, 0.4177771136967e+01, 0.1102062672231e+00],
        
        [0.1578057810499e-09, 0.5782460597335e+01, 0.7626583626240e-01],
        [0.1237713253351e-09, 0.5705900866881e+01, 0.5154640627760e+00],
        [0.1313076837395e-09, 0.5163438179576e+01, 0.3664874755930e-01],
        [0.1184963304860e-09, 0.3054804427242e+01, 0.6327837846670e+00],
        [0.1238130878565e-09, 0.2317292575962e+01, 0.3961708870310e-01],
        [0.1015959527736e-09, 0.2194643645526e+01, 0.7329749511860e-01],
        [0.9017954423714e-10, 0.2868603545435e+01, 0.1990721704425e+00],
        [0.8668024955603e-10, 0.4923849675082e+01, 0.5439178814476e+00],
        [0.7756083930103e-10, 0.3014334135200e+01, 0.9491756770005e+00],
        [0.7536503401741e-10, 0.2704886279769e+01, 0.1030928125552e+00],
        
        [0.5483308679332e-10, 0.6010983673799e+01, 0.8531963191132e+00],
        [0.5184339620428e-10, 0.1952704573291e+01, 0.2093666171530e+00],
        [0.5108658712030e-10, 0.2958575786649e+01, 0.2172315424036e+00],
        [0.5019424524650e-10, 0.1736317621318e+01, 0.2164800718209e+00],
        [0.4909312625978e-10, 0.3167216416257e+01, 0.2101180877357e+00],
        [0.4456638901107e-10, 0.7697579923471e+00, 0.3235053470014e+00],
        [0.4227030350925e-10, 0.3490910137928e+01, 0.6373574839730e-01],
        [0.4095456040093e-10, 0.5178888984491e+00, 0.6470106940028e+00],
        [0.4990537041422e-10, 0.3323887668974e+01, 0.1422690933580e-01],
        [0.4321170010845e-10, 0.4288484987118e+01, 0.7358765972222e+00],
        
        [0.3544072091802e-10, 0.6021051579251e+01, 0.5265099800692e+00],
        [0.3480198638687e-10, 0.4600027054714e+01, 0.5328719641544e+00],
        [0.3440287244435e-10, 0.4349525970742e+01, 0.8582758298370e-01],
        [0.3330628322713e-10, 0.2347391505082e+01, 0.1104591729320e-01],
        [0.2973060707184e-10, 0.4789409286400e+01, 0.5257585094865e+00],
        [0.2932606766089e-10, 0.5831693799927e+01, 0.5336234347371e+00],
        [0.2876972310953e-10, 0.2692638514771e+01, 0.1173197218910e+00],
        [0.2827488278556e-10, 0.2056052487960e+01, 0.2022531624851e+00],
        [0.2515028239756e-10, 0.7411863262449e+00, 0.9597935788730e-01],
        [0.2853033744415e-10, 0.3948481024894e+01, 0.2118763888447e+01]  
        
    ])

    # SSB-to-Sun, T^1, Y
    s1y = np.array([
        [0.8989047573576e-08, 0.5840593672122e+01, 0.4265981595566e+00],
        [0.7815938401048e-08, 0.1129664707133e+01, 0.2061856251104e+00],
        [0.7550926713280e-08, 0.6196589104845e+00, 0.2204125344462e+00],
        [0.6056556925895e-08, 0.1677494667846e+01, 0.1059381944224e+01],
        [0.5734142698204e-08, 0.4000920852962e+01, 0.5225775174439e+00],
        [0.5614341822459e-08, 0.3486722577328e+01, 0.5368044267797e+00],
        [0.1028678147656e-08, 0.1877141024787e+01, 0.7113454667900e-02],
        [0.7270792075266e-09, 0.5077167301739e+01, 0.6398972393349e+00],
        [0.8734141726040e-09, 0.9069550282609e-01, 0.4194847048887e+00],
        [0.5377371402113e-09, 0.6039381844671e+01, 0.4337116142245e+00],
        
        [0.4729719431571e-09, 0.2153086311760e+01, 0.2132990797783e+00],
        [0.4458052820973e-09, 0.5059830025565e+01, 0.5296909721118e+00],
        [0.4406855467908e-09, 0.2027971692630e+01, 0.1589072916335e+01],
        [0.3101659310977e-09, 0.3317677981860e+01, 0.1052268489556e+01],
        [0.3016749232545e-09, 0.3913703482532e+01, 0.1066495398892e+01],
        [0.3198541352656e-09, 0.1275513098525e+01, 0.1495633313810e+00],
        [0.2142065389871e-09, 0.5301351614597e+01, 0.3163918923335e+00],
        [0.1902615247592e-09, 0.4894943352736e+00, 0.2275259891141e+00],
        [0.1613410990871e-09, 0.2449891130437e+01, 0.1102062672231e+00],
        [0.1576992165097e-09, 0.4211421447633e+01, 0.7626583626240e-01],
        
        [0.1241637259894e-09, 0.4140803368133e+01, 0.5154640627760e+00],
        [0.1313974830355e-09, 0.3591920305503e+01, 0.3664874755930e-01],
        [0.1181697118258e-09, 0.1506314382788e+01, 0.6327837846670e+00],
        [0.1238239742779e-09, 0.7461405378404e+00, 0.3961708870310e-01],
        [0.1010107068241e-09, 0.6271010795475e+00, 0.7329749511860e-01],
        [0.9226316616509e-10, 0.1259158839583e+01, 0.1990721704425e+00],
        [0.8664946419555e-10, 0.3353244696934e+01, 0.5439178814476e+00],
        [0.7757230468978e-10, 0.1447677295196e+01, 0.9491756770005e+00],
        [0.7693168628139e-10, 0.1120509896721e+01, 0.1030928125552e+00],
        [0.5487897454612e-10, 0.4439380426795e+01, 0.8531963191132e+00],
        
        [0.5196118677218e-10, 0.3788856619137e+00, 0.2093666171530e+00],
        [0.5110853339935e-10, 0.1386879372016e+01, 0.2172315424036e+00],
        [0.5027804534813e-10, 0.1647881805466e+00, 0.2164800718209e+00],
        [0.4922485922674e-10, 0.1594315079862e+01, 0.2101180877357e+00],
        [0.6155599524400e-10, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.4447147832161e-10, 0.5480720918976e+01, 0.3235053470014e+00],
        [0.4144691276422e-10, 0.1931371033660e+01, 0.6373574839730e-01],
        [0.4099950625452e-10, 0.5229611294335e+01, 0.6470106940028e+00],
        [0.5060541682953e-10, 0.1731112486298e+01, 0.1422690933580e-01],
        [0.4293615946300e-10, 0.2714571038925e+01, 0.7358765972222e+00],
        
        [0.3545659845763e-10, 0.4451041444634e+01, 0.5265099800692e+00],
        [0.3479112041196e-10, 0.3029385448081e+01, 0.5328719641544e+00],
        [0.3438516493570e-10, 0.2778507143731e+01, 0.8582758298370e-01],
        [0.3297341285033e-10, 0.7898709807584e+00, 0.1104591729320e-01],
        [0.2972585818015e-10, 0.3218785316973e+01, 0.5257585094865e+00],
        [0.2931707295017e-10, 0.4260731012098e+01, 0.5336234347371e+00],
        [0.2897198149403e-10, 0.1120753978101e+01, 0.1173197218910e+00],
        [0.2832293240878e-10, 0.4597682717827e+00, 0.2022531624851e+00],
        [0.2864348326612e-10, 0.2169939928448e+01, 0.9597935788730e-01],
        [0.2852714675471e-10, 0.2377659870578e+01, 0.2118763888447e+01]  
        
    ])

    # SSB-to-Sun, T^1, Z
    s1z = np.array([
        [0.5444220475678e-08, 0.1803825509310e+01, 0.2132990797783e+00],
        [0.3883412695596e-08, 0.4668616389392e+01, 0.5296909721118e+00],
        [0.1334341434551e-08, 0.0000000000000e+00, 0.0000000000000e+00],
        [0.3730001266883e-09, 0.5401405918943e+01, 0.2061856251104e+00],
        [0.2894929197956e-09, 0.4932415609852e+01, 0.2204125344462e+00],
        [0.2857950357701e-09, 0.3154625362131e+01, 0.7478166569050e-01],
        [0.2499226432292e-09, 0.3657486128988e+01, 0.4265981595566e+00],
        [0.1937705443593e-09, 0.5740434679002e+01, 0.1059381944224e+01],
        [0.1374894396320e-09, 0.1712857366891e+01, 0.5368044267797e+00],
        [0.1217248678408e-09, 0.2312090870932e+01, 0.5225775174439e+00],
        
        [0.7961052740870e-10, 0.5283368554163e+01, 0.3813291813120e-01],
        [0.4979225949689e-10, 0.4298290471860e+01, 0.4194847048887e+00],
        [0.4388552286597e-10, 0.6145515047406e+01, 0.7113454667900e-02],
        [0.2586835212560e-10, 0.3019448001809e+01, 0.6398972393349e+00]
    ])

    # SSB-to-Sun, T^2, X
    s2x = np.array([
        [0.1603551636587e-11, 0.4404109410481e+01, 0.2061856251104e+00],
        [0.1556935889384e-11, 0.4818040873603e+00, 0.2204125344462e+00],
        [0.1182594414915e-11, 0.9935762734472e+00, 0.5225775174439e+00],
        [0.1158794583180e-11, 0.3353180966450e+01, 0.5368044267797e+00],
        [0.9597358943932e-12, 0.5567045358298e+01, 0.2132990797783e+00],
        [0.6511516579605e-12, 0.5630872420788e+01, 0.4265981595566e+00],
        [0.7419792747688e-12, 0.2156188581957e+01, 0.5296909721118e+00],
        [0.3951972655848e-12, 0.1981022541805e+01, 0.1059381944224e+01],
        [0.4478223877045e-12, 0.0000000000000e+00, 0.0000000000000e+00]
    ])

    # SSB-to-Sun, T^2, Y
    s2y = np.array([
        [0.1609114495091e-11, 0.2831096993481e+01, 0.2061856251104e+00],
        [0.1560330784946e-11, 0.5193058213906e+01, 0.2204125344462e+00],
        [0.1183535479202e-11, 0.5707003443890e+01, 0.5225775174439e+00],
        [0.1158183066182e-11, 0.1782400404928e+01, 0.5368044267797e+00],
        [0.1032868027407e-11, 0.4036925452011e+01, 0.2132990797783e+00],
        [0.6540142847741e-12, 0.4058241056717e+01, 0.4265981595566e+00],
        [0.7305236491596e-12, 0.6175401942957e+00, 0.5296909721118e+00],
        [-0.5580725052968e-12,0.0000000000000e+00, 0.0000000000000e+00],
        [0.3946122651015e-12, 0.4108265279171e+00, 0.1059381944224e+01]
    ])

    # SSB-to-Sun, T^2, Z
    s2z = np.array([
        [0.3749920358054e-12, 0.3230285558668e+01, 0.2132990797783e+00],
        [0.2735037220939e-12, 0.6154322683046e+01, 0.5296909721118e+00]
    ])

    # Group coefficients
    ce0 = [e0x, e0y, e0z]
    ce1 = [e1x, e1y, e1z]
    ce2 = [e2x, e2y, e2z]
    cs0 = [s0x, s0y, s0z]
    cs1 = [s1x, s1y, s1z]
    cs2 = [s2x, s2y, s2z]

    t2 = t * t

    # Pre-allocate arrays
    ph = np.zeros(3)  # heliocentric position
    vh = np.zeros(3)  # heliocentric velocity  
    pb = np.zeros(3)  # barycentric position
    vb = np.zeros(3)  # barycentric velocity

    # Process each coordinate component using optimized harmonic function
    for i in range(3):
        # Sun to Earth components
        sun_coeffs = [ce0[i], ce1[i], ce2[i]]
        ph[i], vh_d = pv_harmonic_terms(sun_coeffs, t, t2)
        vh[i] = vh_d / DJY

        # SSB to Sun components
        ssb_coeffs = [cs0[i], cs1[i], cs2[i]]
        pb_a, vb_d_a = pv_harmonic_terms(ssb_coeffs, t, t2)
        pb[i] = ph[i] + pb_a
        vb[i] = vh[i] + vb_d_a / DJY

    # Apply orientation matrix transformation
    pvh = np.zeros((2, 3))
    pvb = np.zeros((2, 3))

    # For position and velocity
    pvh[0] = ori_mat @ ph
    pvh[1] = ori_mat @ vh

    # For barycentric
    pvb[0] = ori_mat @ pb
    pvb[1] = ori_mat @ vb

    return pvh, pvb


def pymEqec06(date1, date2, dr, dd):
    """
    Transform ICRS equatorial coordinates to ecliptic coordinates
    (mean equinox and ecliptic of date) using IAU 2006 precession model,
    with pym* utilities.

    Parameters
    ----------
    date1, date2 : float
        TT as 2-part Julian date
    dr, dd : float
        ICRS right ascension and declination in radians

    Returns
    -------
    dl, db : float
        Ecliptic longitude and latitude in radians
    """
    #  Spherical to Cartesian
    v1 = pymS2c(dr, dd)

    #  Rotation matrix: ICRS equatorial to ecliptic
    rm = pymEcm06(date1, date2)

    #  Transformation from ICRS to ecliptic (matrix-vector product)
    v2 = rm @ v1

    #  Cartesian to spherical
    a, b = pymC2s(v2)

    #  Normalize angles
    dl = pymAnp(a)   # 0 to 2pi
    db = pymAnpm(b)  # -pi to +pi

    return dl, db


def pymFk5hip():
    """
    FK5 to Hipparcos rotation and spin vectors.

    Returns
    -------
    r5h : ndarray, shape (3,3)
        Rotation matrix: transforms FK5 vectors to Hipparcos frame
    s5h : ndarray, shape (3,)
        Spin vector: FK5 spin with respect to Hipparcos [radians/year]

    Notes
    -----
    Reference: Mignard & Froeschlé, A&A 354, 732–739 (2000)
    """
    # FK5 wrt Hipparcos orientation (radians)
    epx = -19.9e-3 * DAS2R
    epy =  -9.1e-3 * DAS2R
    epz =  22.9e-3 * DAS2R

    # FK5 wrt Hipparcos spin (radians per Julian year)
    omx = -0.30e-3 * DAS2R
    omy =  0.60e-3 * DAS2R
    omz =  0.70e-3 * DAS2R

    # Orientation vector
    v = np.array([epx, epy, epz])

    # Convert rotation vector -> matrix
    r5h = pymRv2m(v)

    # Spin vector
    s5h = np.array([omx, omy, omz])

    return r5h, s5h



def pymFk5hz(r5, d5, date1, date2):
    """
    Transform an FK5 (J2000.0) star position into the Hipparcos system,
    assuming zero Hipparcos proper motion.

    Parameters
    ----------
    r5 : float
        FK5 right ascension [radians], equinox J2000.0, at date.
    d5 : float
        FK5 declination [radians], equinox J2000.0, at date.
    date1, date2 : float
        TDB date (Julian Date, split in any convenient way).

    Returns
    -------
    rh : float
        Hipparcos right ascension [radians].
    dh : float
        Hipparcos declination [radians].
    """

    # Interval from given date to fundamental epoch J2000.0 (Julian years)
    t = -((date1 - DJ00) + date2) / DJY

    # FK5 barycentric position vector
    p5e = pymS2c(r5, d5)

    # FK5 to Hipparcos orientation matrix and spin vector
    r5h, s5h = pymFk5hip()

    # Accumulated Hipparcos wrt FK5 spin over that interval
    vst = pymSxp(t, s5h)

    # Express the accumulated spin as a rotation matrix
    rst = pymRv2m(vst)

    # Derotate the FK5 vector back to the given date
    p5 = pymTrxp(rst, p5e)

    # Rotate the vector into the Hipparcos system
    ph = pymRxp(r5h, p5)

    # Convert Hipparcos vector to spherical coordinates
    w, dh = pymC2s(ph)
    rh    = pymAnp(w)

    return rh, dh



def pymFk45z(r1950, d1950, bepoch):
    """
    Convert a B1950.0 FK4 star position to J2000.0 FK5,
    assuming zero FK5 proper motion.

    Parameters
    ----------
    r1950 : float
        B1950 right ascension [radians].
    d1950 : float
        B1950 declination [radians].
    bepoch : float
        Besselian epoch of the observation.

    Returns
    -------
    r2000 : float
        J2000 FK5 right ascension [radians].
    d2000 : float
        J2000 FK5 declination [radians].
    """
    PMF = 100.0 * DR2AS  # rad/year to arcsec/century

    # Canonical constants (Seidelmann 1992)
    a_vec = np.array([-1.62557e-6, -0.31919e-6, -0.13843e-6])
    ad_vec = np.array([1.245e-3, -1.580e-3, -0.659e-3])

    # 3x3x2 matrix for transformation
    em = np.array([
        [[+0.9999256782, -0.0111820611, -0.0048579477],
         [+0.0111820610, +0.9999374784, -0.0000271765],
         [+0.0048579479, -0.0000271474, +0.9999881997]],
        [[-0.000551, -0.238565, +0.435739],
         [+0.238514, -0.002667, -0.008541],
         [-0.435623, +0.012254, +0.002117]]
    ])

    # Convert spherical coordinates to position vector
    r0 = pymS2c(r1950, d1950)

    # Adjust vector A to give zero proper motion in FK5
    w = (bepoch - 1950) / PMF
    p = pymPpsp(a_vec, w, ad_vec)

    # Remove E-terms
    p = pymPpsp(p, -pymPdp(r0, p), r0)
    p = pymPmp(r0, p)

    # Vectorized conversion to Fricke system pv-vector
    pv = np.tensordot(em, p, axes=([2],[0]))  # shape (2,3)

    # Allow for fictitious proper motion
    djm0, djm = pymEpb2jd(bepoch)
    w = (pymEpj(djm0, djm) - 2000.0) / PMF
    pv = pymPvu(w, pv)

    # Convert back to spherical coordinates
    r2000, d2000 = pymC2s(pv[0])
    r2000 = pymAnp(r2000)

    return r2000, d2000


def pymFk52h(r5, d5, dr5, dd5, px5, rv5):
    """
    Transform FK5 (J2000.0) star data into the Hipparcos system.

    Parameters
    ----------
    r5, d5 : float
        FK5 RA and Dec (radians, J2000.0)
    dr5, dd5 : float
        Proper motion in RA, Dec (rad/year)
    px5 : float
        Parallax (arcsec)
    rv5 : float
        Radial velocity (km/s, positive = receding)

    Returns
    -------
    rh, dh : float
        Hipparcos RA, Dec (radians)
    drh, ddh : float
        Proper motion in RA, Dec (rad/year)
    pxh : float
        Parallax (arcsec)
    rvh : float
        Radial velocity (km/s)
    """
    # Convert star catalog data to space motion pv-vector 
    pv5 = pymStarpv(r5, d5, dr5, dd5, px5, rv5)  # shape (2,3)

    # FK5 to Hipparcos orientation matrix and spin vector
    r5h, s5h = pymFk5hip()  # r5h: 3x3, s5h: 3-vector

    # Convert spin from per year to per day
    s5h = s5h / 365.25

    # Rotate FK5 position into Hipparcos system
    pvh    = np.zeros((2, 3))
    pvh[0] = pymRxp(r5h, pv5[0])

    # Spin contribution to space motion using pymPxp
    wxp = pymPxp(pv5[0], s5h)

    # Add spin contribution to FK5 space motion
    #vv = pymPpp(wxp, pv5[1])
    vv = wxp +  pv5[1] 

    # Rotate FK5 space motion into Hipparcos system
    pvh[1] = pymRxp(r5h, vv)

    # Convert Hipparcos pv-vector to star catalog data
    rh, dh, drh, ddh, pxh, rvh = pymPvstar(pvh)

    return rh, dh, drh, ddh, pxh, rvh

 

def pymFk54z(r2000, d2000, bepoch):
    """
    Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
    proper motion in FK5 and parallax.

    Parameters
    ----------
    r2000 : float
        J2000.0 RA (rad)
    d2000 : float
        J2000.0 Dec (rad)
    bepoch : float
        Besselian epoch (e.g., 1950.0)

    Returns
    -------
    r1950 : float
        B1950.0 RA (rad)
    d1950 : float
        B1950.0 Dec (rad)
    dr1950 : float
        B1950.0 fictitious proper motion in RA (rad/trop.yr)
    dd1950 : float
        B1950.0 fictitious proper motion in Dec (rad/trop.yr)
    """
    # --- FK5 to FK4 conversion assuming zero pm, parallax, radial velocity ---
    r, d, pr, pd, _, _ = pymFk524(r2000, d2000, 0.0, 0.0, 0.0, 0.0)

    # --- Convert to Cartesian (unit vector) ---
    
    p = pymS2c(r, d);
    # --- Fictitious proper motion in Cartesian ---
    v = np.zeros(3)
    v[0] = - pr * p[1] - pd * np.cos(r) * np.sin(d)
    v[1] =   pr * p[0] - pd * np.sin(r) * np.sin(d)
    v[2] =               pd * np.cos(d)
 

    # --- Apply the motion for epoch difference ---
    dt = bepoch - 1950.0
    p += dt * v

   #  Convert back to spherical coordinates 
    r1950, d1950 = pymC2s(p)
    r1950 = pymAnp(r1950)  # normalize RA to [0, 2π)

    #  Return B1950.0 coordinates and proper motion 
    dr1950, dd1950 = pr, pd
    
    return r1950, d1950, dr1950, dd1950



def pymFk425(r1950, d1950, dr1950, dd1950, p1950, v1950):
    """
    Convert B1950.0 FK4 star catalog data to J2000.0 FK5 system.

    Parameters
    ----------
    r1950 : float
        B1950.0 RA (rad)
    d1950 : float
        B1950.0 Dec (rad)
    dr1950 : float
        Proper motion in RA (rad/tropical year)
    dd1950 : float
        Proper motion in Dec (rad/tropical year)
    p1950 : float
        Parallax (arcsec)
    v1950 : float
        Radial velocity (km/s, +ve = receding)

    Returns
    -------
    r2000 : float
        J2000.0 RA (rad)
    d2000 : float
        J2000.0 Dec (rad)
    dr2000 : float
        Proper motion in RA (rad/Julian year)
    dd2000 : float
        Proper motion in Dec (rad/Julian year)
    p2000 : float
        Parallax (arcsec)
    v2000 : float
        Radial velocity (km/s)
    """
 
    PMF = 100.0 * DR2AS          # rad/yr -> arcsec/century
    VF  = 21.095                 # km/s to au per tropical century
    TINY = 1e-30

    # Constant pv-vectors from Seidelmann (1992) 
    a = np.array([
        [-1.62557e-6, -0.31919e-6, -0.13843e-6],
        [ 1.245e-3,   -1.580e-3,   -0.659e-3]
    ])

    em = np.array([
        [
            [[+0.9999256782, -0.0111820611, -0.0048579477],
             [+0.00000242395018, -0.00000002710663, -0.00000001177656]],
            [[+0.0111820610, +0.9999374784, -0.0000271765],
             [+0.00000002710663, +0.00000242397878, -0.00000000006587]],
            [[+0.0048579479, -0.0000271474, +0.9999881997],
             [+0.00000001177656, -0.00000000006582, +0.00000242410173]]
        ],
        [
            [[-0.000551, -0.238565, +0.435739],
             [+0.99994704, -0.01118251, -0.00485767]],
            [[+0.238514, -0.002667, -0.008541],
             [+0.01118251, +0.99995883, -0.00002718]],
            [[-0.435623, +0.012254, +0.002117],
             [+0.00485767, -0.00002714, +1.00000956]]
        ]
    ])

    #  Convert input to FK4 pv-vector 
    ur = dr1950 * PMF
    ud = dd1950 * PMF
    px = p1950
    rv = v1950
    pxvf = px * VF
    w = rv * pxvf

    # spherical -> pv (FK4, B1950)
    r0 = pymS2pv(r1950, d1950, 1.0, ur, ud, w)

    #  Apply E-terms of aberration 
    pv1 = r0 - a
    pv1 += np.array([
        np.dot(r0[0], a[0]) * r0[0],
        np.dot(r0[0], a[1]) * r0[0]
    ])

    #  Transform to FK5 using einsum 
    pv2 = np.einsum('kl,ijkl->ij', pv1, em)

    #  Convert back to spherical catalog form 
    r, d, w, ur, ud, rd = pymPv2s(pv2)

    if px > TINY:
        rv = rd / pxvf
        px = px / w

     
    r2000 = pymAnp(r)
    d2000 = d
    dr2000 = ur / PMF
    dd2000 = ud / PMF
    v2000 = rv
    p2000 = px

    return r2000, d2000, dr2000, dd2000, p2000, v2000


def pymFk524(r2000, d2000, dr2000, dd2000, p2000, v2000):
    """
    Convert J2000.0 FK5 star catalog data to B1950.0 FK4 

    Parameters
    ----------
    r2000, d2000 : float
        J2000.0 right ascension and declination (radians, FK5)
    dr2000, dd2000 : float
        J2000.0 proper motions (radians per Julian year, FK5)
    p2000 : float
        Parallax (arcseconds)
    v2000 : float
        Radial velocity (km/s, positive = receding)

    Returns
    -------
    r1950, d1950 : float
        B1950.0 FK4 right ascension and declination (radians)
    dr1950, dd1950 : float
        B1950.0 FK4 proper motions (radians per tropical year)
    p1950 : float
        Parallax (arcseconds)
    v1950 : float
        Radial velocity (km/s)
    """
    #  Constants  
    PMF = 100.0 * DR2AS  # rad/year -> arcsec/century
    TINY = 1e-30
    VF = 21.095  # km/s -> AU/tropical century

    # A and Adot vectors
    a = np.array([
        [-1.62557e-6, -0.31919e-6, -0.13843e-6],
        [ 1.245e-3,   -1.580e-3,   -0.659e-3]
    ])

    # em tensor (3x2x3x2) from Seidelmann 1992
    em = np.array([
        [
            [[+0.9999256795,     +0.0111814828,     +0.0048590039],
             [-0.00000242389840, -0.00000002710544, -0.00000001177742]],
            [[-0.0111814828,     +0.9999374849,     -0.0000271771],
             [+0.00000002710544, -0.00000242392702, +0.00000000006585]],
            [[-0.0048590040,     -0.0000271557,     +0.9999881946],
             [+0.00000001177742, +0.00000000006585, -0.00000242404995]]
        ],
        [
            [[-0.000551,         +0.238509,         -0.435614],
             [+0.99990432,       +0.01118145,       +0.00485852]],
            [[-0.238560,         -0.002667,         +0.012254],
             [-0.01118145,       +0.99991613,       -0.00002717]],
            [[+0.435730,         -0.008541,         +0.002117],
             [-0.00485852,       -0.00002716,       +0.99996684]]
        ]
    ])

    # FK5 input  
    r  = r2000
    d  = d2000
    ur = dr2000 * PMF
    ud = dd2000 * PMF
    px = p2000
    rv = v2000
    pxvf = px * VF
    radial = rv * pxvf

    # Convert to pv-vector  
    r0 = pymS2pv(r, d, 1.0, ur, ud, radial)  # shape (2,3)

    # Convert pv-vector to Bessel-Newcomb system  
    r1 = np.einsum('ijkl,kl->ij', em, r0)
    
    w  = np.linalg.norm(r1[0])
    p1 = r1[0] + (w * a[0] - np.dot(r1[0], a[0]) * r1[0])

    # Apply E-terms using numpy  
    w   = np.linalg.norm(p1)
    pv0 = r1[0] + (w * a[0] - np.dot(r1[0], a[0]) * r1[0])
    pv1 = r1[1] + (w * a[1] - np.dot(r1[0], a[1]) * pv0)
    pv  = np.vstack([pv0, pv1])

    #  Convert pv -> spherical coordinates 
    r, d, w, ur, ud, rd = pymPv2s(pv)

    #  Adjust parallax and velocity 
    if px > TINY:
        rv = rd / pxvf
        px /= w

    #  Return FK4 results 
    r1950 = pymAnp(r)
    #  r1950,   d1950, dr1950, dd1950,  p1950, v1950
    return r1950, d, ur / PMF, ud / PMF, px, rv


def pymG2icrs(dl, db):
    """
    Transform Galactic coordinates (l, b) to ICRS equatorial coordinates (RA, Dec).

    Parameters
    ----------
    dl : float
        Galactic longitude in radians
    db : float
        Galactic latitude in radians

    Returns
    -------
    dr : float
        ICRS right ascension in radians
    dd : float
        ICRS declination in radians
    """
    # ICRS to Galactic rotation matrix (Hipparcos canonical)
    r = np.array([
        [-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116],
        [ 0.494109427875583673525222371358, -0.444829629960011178146614061616,  0.746982244497218890527388004556],
        [-0.867666149019004701181616534570, -0.198076373431201528180486091412,  0.455983776175066922272100478348]
    ])

    # Spherical to Cartesian
    v1 = pymS2c(dl, db)

    # Galactic to ICRS (transpose of rotation matrix)
    v2 = r.T @ v1   

    # Cartesian to spherical
    dr, dd = pymC2s(v2)

    # Normalize angles
    dr = pymAnp(dr)   # 0 to 2pi
    dd = pymAnpm(dd)  # -pi to +pi

    return dr, dd



def pymGc2gd(n, xyz):
    """
    Transform geocentric coordinates to geodetic using the specified reference ellipsoid.

    Parameters
    ----------
    n : int
        Ellipsoid identifier:
            1 : WGS84
            2 : GRS80
            3 : WGS72
    xyz : array-like of float, shape (3,)
        Geocentric coordinates [x, y, z] in meters.

    Returns
    -------
    elong : float
        Longitude (radians, east positive)
    phi : float
        Geodetic latitude (radians)
    height : float
        Height above ellipsoid (meters)

    Raises
    ------
    ValueError
        If the ellipsoid identifier is invalid, or the coordinates are illegal.
    """

    a, f = pymEform(n)

    #  pymGc2gde  (longitude, latitude, height); it raises on illegal input
    elong, phi, height = pymGc2gde(a, f, xyz)

    return elong, phi, height


def pymGc2gde(a, f, xyz):
    """
    Transform geocentric coordinates to geodetic coordinates
    for a reference ellipsoid of given equatorial radius and flattening.

    Parameters
    ----------
    a : float
        Equatorial radius of ellipsoid [m]
    f : float
        Flattening
    xyz : array-like, shape (3,)
        Geocentric coordinates [x, y, z] in meters

    Returns
    -------
    lon : float
        Longitude (radians, east positive)
    lat : float
        Geodetic latitude (radians)
    h : float
        Height above ellipsoid [m]

    Raises
    ------
    ValueError
        If the ellipsoid parameters are illegal (flattening out of range,
        non-positive equatorial radius, or invalid eccentricity).
    """

    xyz = np.asarray(xyz, dtype=float).reshape(3)

    # Validate ellipsoid parameters
    if f < 0.0 or f >= 1.0:
        raise ValueError("illegal flattening f = %r (must satisfy 0 <= f < 1)" % (f,))
    if a <= 0.0:
        raise ValueError("illegal equatorial radius a = %r (must be > 0)" % (a,))

    e2 = (2.0 - f) * f
    e4t = 1.5 * e2 * e2
    ec2 = 1.0 - e2
    if ec2 <= 0.0:
        raise ValueError("illegal ellipsoid: second eccentricity <= 0")
    ec = np.sqrt(ec2)
    b = a * ec

    x, y, z = xyz
    p2 = x*x + y*y

    # Longitude
    lon = np.arctan2(y, x) if p2 > 0.0 else 0.0
    absz = np.abs(z)

    aeps2 = a*a * 1e-32  # small number to avoid division by zero

    if p2 > aeps2:
        p = np.sqrt(p2)
        s0 = absz / a
        pn = p / a
        zc = ec * s0

        # Halley correction factors
        c0 = ec * pn
        c02 = c0 * c0
        c03 = c02 * c0
        s02 = s0 * s0
        s03 = s02 * s0
        a02 = c02 + s02
        a0 = np.sqrt(a02)
        a03 = a02 * a0
        d0 = zc * a03 + e2 * s03
        f0 = pn * a03 - e2 * c03
        b0 = e4t * s02 * c02 * pn * (a0 - ec)
        s1 = d0 * f0 - b0 * s0
        cc = ec * (f0*f0 - b0 * c0)

        lat = np.arctan(s1 / cc)
        s12 = s1 * s1
        cc2 = cc * cc
        h = (p * cc + absz * s1 - a * np.sqrt(ec2*s12 + cc2)) / np.sqrt(s12 + cc2)
    else:
        # Polar case
        lat = DPI / 2.0
        h = absz - b

    if z < 0:
        lat = -lat

    return lon, lat, h


def pymH2fk5(rh, dh, drh, ddh, pxh, rvh):
    """
    Transform Hipparcos star data into FK5 (J2000.0) system.

    Parameters
    ----------
    rh : float
        Hipparcos right ascension [radians]
    dh : float
        Hipparcos declination [radians]
    drh : float
        Proper motion in RA (rad/year)
    ddh : float
        Proper motion in Dec (rad/year)
    pxh : float
        Parallax [arcsec]
    rvh : float
        Radial velocity [km/s, positive = receding]

    Returns
    -------
    r5, d5 : float
        FK5 right ascension, declination [radians]
    dr5, dd5 : float
        Proper motions (rad/year)
    px5 : float
        Parallax [arcsec]
    rv5 : float
        Radial velocity [km/s]
    """

    #  Hipparcos barycentric pv-vector 
    pvh = pymStarpv(rh, dh, drh, ddh, pxh, rvh)  # shape (2,3)

    #  FK5 -> Hipparcos rotation matrix and spin vector 
    r5h, s5h = pymFk5hip()  # rotation matrix (3x3), spin vector (3,)

    #  Convert spin units from per year to per day 
    s5h /= 365.25

    #  Orient spin into Hipparcos system 
    sh = pymRxp(r5h, s5h)

    #  De-orient Hipparcos position into FK5 system 
    pv5_pos = r5h.T @ pvh[0]

    #  Compute extra space motion from spin × position 
    wxp = np.cross(pvh[0], sh)
    #wxp = pymPxp(pvh[0], sh)  

    #  Subtract spin contribution from Hipparcos motion 
    vv = pvh[1] - wxp

    #  De-orient Hipparcos space motion into FK5 system 
    pv5_vel = r5h.T @ vv

    #  Combine position + velocity 
    pv5 = np.array([pv5_pos, pv5_vel])

    #  Convert back to catalog form (FK5 spherical) 
    r5, d5, dr5, dd5, px5, rv5 = pymPvstar(pv5)

    return r5, d5, dr5, dd5, px5, rv5


def pymHfk5z(rh, dh, date1, date2):
    """
    Transform a Hipparcos star position into FK5 (J2000.0),
    assuming zero Hipparcos proper motion.

    Parameters
    ----------
    rh, dh : float
        Hipparcos right ascension, declination [radians]
    date1, date2 : float
        TDB Julian date (any split, e.g., JD or J2000 format)

    Returns
    -------
    r5, d5 : float
        FK5 right ascension, declination [radians]
    dr5, dd5 : float
        FK5 fictitious proper motions [radians/year]
    """
 
    #  Time interval from J2000.0 to given date (Julian years)
    t = ((date1 - DJ00) + date2) / DJY

    #  Hipparcos unit position vector
    ph = pymS2c(rh, dh)  # shape (3,)

    #  FK5 -> Hipparcos rotation matrix & spin vector
    r5h, s5h = pymFk5hip()

    #  Rotate spin into Hipparcos system
    sh = pymRxp(r5h, s5h)

    #  Accumulated spin (per JY)
    vst = t * s5h

    #  Express accumulated spin as rotation matrix
    rst = pymRv2m(vst)

    #  Combined rotation: accumulated spin then FK5 -> Hipparcos
    r5ht = r5h @ rst

    #  De-orient + de-spin Hipparcos position into FK5 J2000.0
    pv5e_pos = r5ht.T @ ph

    #  Apply spin to position to get space motion
    vv = np.cross(sh, ph)

    #  De-orient + de-spin space motion into FK5 J2000.0
    pv5e_vel = r5ht.T @ vv

    #  FK5 pv-vector -> spherical
    pv5e = np.array([pv5e_pos, pv5e_vel])
 
    w, d5, r, dr5, dd5, v  =  pymPv2s(pv5e)    
    
    r5 = pymAnp(w)

    return r5, d5, dr5, dd5


def pymIcrs2g(dr, dd):
    """
    Transform ICRS equatorial coordinates (RA, Dec) to Galactic coordinates (l, b).

    Parameters
    ----------
    dr : float
        ICRS right ascension in radians
    dd : float
        ICRS declination in radians

    Returns
    -------
    dl : float
        Galactic longitude in radians
    db : float
        Galactic latitude in radians
    """
    # ICRS to Galactic rotation matrix (Hipparcos canonical, full precision)
    r = np.array([
        [-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116],
        [ 0.494109427875583673525222371358, -0.444829629960011178146614061616,  0.746982244497218890527388004556],
        [-0.867666149019004701181616534570, -0.198076373431201528180486091412,  0.455983776175066922272100478348]
    ])

    # Spherical to Cartesian
    v1 = pymS2c(dr, dd)

    # ICRS to Galactic (matrix-vector product)
    v2 = r @ v1   

    # Cartesian to spherical
    dl, db = pymC2s(v2)

    # Normalize angles
    dl = pymAnp(dl)   # 0 to 2pi
    db = pymAnpm(db)  # -pi to +pi

    return dl, db



def pymLteceq(epj, dl, db):
    """
    Transform ecliptic coordinates (mean equinox and ecliptic of date)
    to ICRS RA,Dec using long-term precession model (Vondrak et al.).

    Parameters
    ----------
    epj : float
        Julian epoch (TT)
    dl, db : float
        Ecliptic longitude and latitude in radians

    Returns
    -------
    dr, dd : float
        ICRS right ascension and declination in radians
    """
    #  Spherical to Cartesian
    v1 = pymS2c(dl, db)

    #  Rotation matrix, long-term precession
    rm = pymLtecm(epj)

    #  Transformation from ecliptic to ICRS (transpose matrix)
    v2 = rm.T @ v1   

    #  Cartesian to spherical
    a, b = pymC2s(v2)

    #  Normalize angles
    dr = pymAnp(a)   # 0 to 2pi
    dd = pymAnpm(b)  # -pi to +pi

    return dr, dd


def pymLtecm(epj):
    """
    Long-term precession: rotation matrix from ICRS equatorial to ecliptic frame.
    (Vondrak et al. 2011, 2012; IAU SOFA standard)

    Parameters
    ----------
    epj : float
        Julian epoch (TT)

    Returns
    -------
    rm : ndarray, shape (3, 3)
        ICRS-to-ecliptic rotation matrix at epoch `epj`.

    Notes
    -----
    The matrix transforms in the sense:
        E_ep = rm @ P_ICRS
    where P_ICRS is a direction vector in the ICRS frame.

    The model agrees with IAU 2006 precession at J2000.0,
    remaining within 100 µas over ±200 years and within a few
    arcseconds over ±200,000 years.
 
    """

    # Frame bias (IERS Conventions 2010, Eq. 5.21 / 5.33)
    dx, de, dr = -0.016617 * DAS2R, -0.0068192 * DAS2R, -0.0146 * DAS2R
 

    # Equator pole unit vector (Z-equator)
    p = pymLtpequ(epj)

    # Ecliptic pole unit vector (Z-ecliptic)
    z = pymLtpecl(epj)
 
 
    # Equinox direction (x-axis)
    w    = pymPxp(p, z)
    _, x = pymPn(w)

    # Middle row (y-axis)
    y = pymPxp(z, x)

    #  Construct rotation matrix 
    rm = np.stack((x, y, z), axis=0)
    
    #  Apply frame bias correction matrix 
    bias = np.array([
        [1.0,  -dr,   dx],
        [dr,   1.0,   de],
        [-dx,  -de,  1.0]
    ])
    
    # Combine with frame bias
    rm = rm @ bias.T  
    
    return rm



def pymLtpecl(epj):
    """
    Long-term precession of the ecliptic (Vondrak et al. 2011, 2012 model).
   
    Parameters
    ----------
    epj : float
        Julian epoch (TT).

    Returns
    -------
    vec : ndarray, shape (3,)
        Ecliptic pole unit vector, referred to the J2000.0 mean equator and equinox.

    Notes
    -----
    - Based on the IAU SOFA routine 'iauLtpecl' (release 2023-10-11).
    - The model remains accurate within 100 microarcseconds over the 20th–21st centuries,
      and within a few arcseconds over ±200,000 years.
    """

    # Constants
    eps0 = 84381.406 * DAS2R  # Obliquity at J2000.0 (radians)

    # Polynomial coefficients (aligned)
    pqpol = np.array([
        [ 5851.607687, -0.1189000, -0.00028913,  0.000000101],
        [-1600.886300,  1.1689818, -0.00000020, -0.000000437]
    ])

    # Periodic coefficients (aligned)
    pqper = np.array([
        [ 708.15, -5486.751211, -684.661560,  667.666730, -5523.863691],
        [2309.00,   -17.127623, 2446.283880, -2354.886252,  -549.747450],
        [1620.00,  -617.517403,  399.671049,  -428.152441,  -310.998056],
        [ 492.20,   413.442940, -356.652376,   376.202861,   421.535876],
        [1183.00,    78.614193, -186.387003,   184.778874,   -36.776172],
        [ 622.00,  -180.732815, -316.800070,   335.321713,  -145.278396],
        [ 882.00,   -87.676083,  198.296701,  -185.138669,   -34.744450],
        [ 547.00,    46.140315,  101.135679,  -120.972830,    22.885731]
    ])

    # Centuries since J2000.0
    t = (epj - 2000.0) / 100.0
    w = D2PI * t

    # Vectorized periodic terms
    a = w / pqper[:, 0]
    s, c = np.sin(a), np.cos(a)
    p = np.sum(c * pqper[:, 1] + s * pqper[:, 3])
    q = np.sum(c * pqper[:, 2] + s * pqper[:, 4])

    # Polynomial terms (Horner’s rule)
    T = np.array([1.0, t, t**2, t**3])
    p += np.dot(pqpol[0], T)
    q += np.dot(pqpol[1], T)

    # Convert to radians
    p *= DAS2R
    q *= DAS2R

    # Compute ecliptic pole vector
    w = 1.0 - p * p - q * q
    w = np.sqrt(np.maximum(w, 0.0))
    s0, c0 = np.sin(eps0), np.cos(eps0)

    vec = np.array([
        p,
        -q * c0 - w * s0,
        -q * s0 + w * c0
    ])

    return vec

#2025-10-17
def pymLtpequ(epj):
    """
    Compute the long-term precession of the equator.
  
    Parameters
    ----------
    epj : float or ndarray
        Julian epoch (TT). Can be scalar or array.

    Returns
    -------
    veq : ndarray, shape (..., 3)
        Equator pole unit vector(s) in J2000.0 mean equator and equinox frame.

    Notes
    -----
    Implements the Vondrak et al. (2011, 2012) long-term precession model.
    Vectorized for batch input using NumPy broadcasting.

    References
    ----------
    Vondrak, J., Capitaine, N. & Wallace, P. (2011), A&A 534, A22  
    Vondrak, J., Capitaine, N. & Wallace, P. (2012), A&A 541, C1
    """

    
    # Polynomial coefficients (X, Y)
    xypol = np.array([
        [   5453.282155,    0.4252841,  -0.00037173,  -0.000000152 ],
        [ -73750.930350,   -0.7675452,  -0.00018725,   0.000000231 ]
    ])  # shape (2,4)

    
    # Periodic coefficients
    xyper = np.array([
        [  256.75,   -819.940624,   75004.344875,   81491.287984,   1558.515853 ],
        [  708.15,  -8444.676815,     624.033993,     787.163481,   7774.939698 ],
        [  274.20,   2600.009459,    1251.136893,    1251.296102,  -2219.534038 ],
        [  241.45,   2755.175630,   -1102.212834,   -1257.950837,  -2523.969396 ],
        [ 2309.00,   -167.659835,   -2660.664980,   -2966.799730,    247.850422 ],
        [  492.20,    871.855056,     699.291817,     639.744522,   -846.485643 ],
        [  396.10,     44.769698,     153.167220,     131.600209,  -1393.124055 ],
        [  288.90,   -512.313065,    -950.865637,    -445.040117,    368.526116 ],
        [  231.10,   -819.415595,     499.754645,     584.522874,    749.045012 ],
        [ 1610.00,   -538.071099,    -145.188210,     -89.756563,    444.704518 ],
        [  620.00,   -189.793622,     558.116553,     524.429630,    235.934465 ],
        [  157.87,   -402.922932,     -23.923029,     -13.549067,    374.049623 ],
        [  220.30,    179.516345,    -165.405086,    -210.157124,   -171.330180 ],
        [ 1200.00,     -9.814756,       9.344131,     -44.919798,    -22.899655 ]
    ])  # shape (14,5)

    # Convert epoch(s) to NumPy array
    epj = np.atleast_1d(epj)
    t   = (epj - 2000.0) / 100.0          # centuries since J2000
    w   = D2PI * t                        # angular frequency term

    
    # Periodic terms (vectorized)
    a = w[:, None] / xyper[:, 0]          # (n, 14)
    s = np.sin(a)
    c = np.cos(a)

    x_per = np.sum(c * xyper[:, 1] + s * xyper[:, 3], axis=1)
    y_per = np.sum(c * xyper[:, 2] + s * xyper[:, 4], axis=1)

    
    # Polynomial terms
    powers = np.vstack([t**i for i in range(4)])  # shape (4, n)
    x_poly = np.dot(xypol[0], powers).T
    y_poly = np.dot(xypol[1], powers).T

    
    # Combine and convert to radians
    x = (x_per + x_poly) * DAS2R
    y = (y_per + y_poly) * DAS2R

    
    # Compute z component
    z_sq = 1.0 - x**2 - y**2
    z = np.where(z_sq > 0.0, np.sqrt(z_sq), 0.0)

    veq = np.stack((x, y, z), axis=-1)
    
    
    return veq[0] if veq.shape[0] == 1 else veq
 

def pymLteqec(epj, ra, dec):
    """
    Transform ICRS right ascension and declination to ecliptic coordinates
    (mean equinox and ecliptic of date) using the long-term precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)
    ra : float or ndarray
        ICRS right ascension in radians
    dec : float or ndarray
        ICRS declination in radians

    Returns
    -------
    lam : float or ndarray
        Ecliptic longitude in radians (range: 0 to 2π)
    beta : float or ndarray
        Ecliptic latitude in radians (range: −π to +π)

    Notes
    -----
    The transformation converts (RA, Dec) in the ICRS frame to
    (λ, β) referred to the mean ecliptic and equinox of the given epoch.

    """

    #  Convert spherical to Cartesian  
    v_eq = pymS2c(ra, dec)

    #  Get long-term precession rotation matrix  
    rm = pymLtecm(epj)

    #  Apply rotation: ICRS -> ecliptic  
    v_ec = np.dot(rm, v_eq)

    #  Convert Cartesian back to spherical  
    lam, beta = pymC2s(v_ec)

    #  Normalize longitude and latitude  
    lam  = pymAnp(lam)     # 0 ≤ λ < 2π
    beta = pymAnpm(beta)   # −π < β ≤ +π

    return lam, beta



def pymPvstar(pv):
    """
    Convert star position+velocity vector to catalog coordinates.

    Parameters
    ----------
    pv : array_like, shape (2,3)
        Position and velocity vector [au, au/day]

    Returns
    -------
    ra, dec : float
        Right ascension and declination (radians)
    pmr, pmd : float
        Proper motions in RA and Dec (radians/year)
    px : float
        Parallax (arcseconds)
    rv : float
        Radial velocity (km/s, positive = receding)
    """

    #  Decompose position vector into modulus and direction
    r = np.linalg.norm(pv[0])
    if r == 0.0:
        raise ValueError("Null position vector")

    pu = pv[0] / r

    #  Isolate radial velocity component
    vr = np.dot(pu, pv[1])
    ur = vr * pu

    #  Isolate transverse velocity component
    ut = pv[1] - ur
    vt = np.linalg.norm(ut)

    #  Special-relativity parameters
    bett = vt / DC
    betr = vr / DC

    d = 1.0 + betr
    w = betr**2 + bett**2
    if d == 0.0 or w >= 1.0:
        raise ValueError("Superluminal speed")

    delta = - w / (np.sqrt(1.0 - w) + 1.0)

    #  Scale tangential velocity into observed frame
    ust = ut / d

    #  Compute observed radial velocity vector
    usr = pu * DC * (betr - delta) / d

    #  Combine radial and tangential vectors
    pv[1] = usr + ust

    #  Convert Cartesian to spherical
    ra, dec, r_obs, rad, decd, rd = pymPv2s(pv)
    if r_obs == 0.0:
        raise ValueError("Null position vector")

    #  Normalize RA to 0..2pi
    #ra = ra % (2*np.pi)
    ra = pymAnp(ra)

    #  Proper motions (radians/year)
    pmr = rad  * DJY
    pmd = decd * DJY

    #  Parallax (arcsec)
    px = DR2AS / r_obs

    #  Radial velocity (km/s)
    rv = 1e-3 * rd * DAU / DAYSEC

    return ra, dec, pmr, pmd, px, rv


def pymStarpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b):
    """
    Update star catalog data for space motion from one epoch to another.

    Parameters
    ----------
    ra1, dec1 : float
        Right ascension and declination at "before" epoch (radians)
    pmr1, pmd1 : float
        Proper motions in RA and Dec (radians/year)
    px1 : float
        Parallax (arcseconds)
    rv1 : float
        Radial velocity (km/s, positive = receding)
    ep1a, ep1b : float
        "Before" epoch (Julian date split)
    ep2a, ep2b : float
        "After" epoch (Julian date split)

    Returns
    -------
    ra2, dec2, pmr2, pmd2, px2, rv2 : float
        Updated star catalog parameters at the "after" epoch
    """

    #  Convert initial star data to pv-vector
    pv1 = pymStarpv(ra1, dec1, pmr1, pmd1, px1, rv1)

    #  Light-time at observation (days)
    tl1 = np.linalg.norm(pv1[0]) / DC

    #  Time interval from "before" to "after" (days)
    dt = (ep2a - ep1a) + (ep2b - ep1b)

    #  Move star along its track to geometric position at "after" epoch
    pv = pymPvu(dt + tl1, pv1)

    #  Compute observed light time at "after" epoch
    r2  = np.dot(pv[0], pv[0])
    rdv = np.dot(pv[0], pv[1])
    v2  = np.dot(pv[1], pv[1])
    c2mv2 = DC**2 - v2
    if c2mv2 <= 0.0:
        raise ValueError("Velocity exceeds speed of light")

    tl2 = (-rdv + np.sqrt(rdv**2 + c2mv2*r2)) / c2mv2

    #  Move from observed place at "before" to observed place at "after"
    pv2 = pymPvu(dt + (tl1 - tl2), pv1)

    #  Convert pv-vector back to catalog parameters at "after" epoch
    ra2, dec2, pmr2, pmd2, px2, rv2 = pymPvstar(pv2)

    return ra2, dec2, pmr2, pmd2, px2, rv2


def pymStarpv(ra, dec, pmr, pmd, px, rv):
    """
    Convert star catalog coordinates to position+velocity vector.

    
    Parameters:
    -----------
    ra : float
        Right ascension (radians)
    dec : float
        Declination (radians)
    pmr : float
        RA proper motion (radians/year)
    pmd : float
        Dec proper motion (radians/year)
    px : float
        Parallax (arcseconds)
    rv : float
        Radial velocity (km/s, positive = receding)

    Returns:
    --------
    pv : np.ndarray
        pv-vector (au, au/day), shape (2, 3)

    Raises:
    -------
    ValueError
        If relativistic solution doesn't converge
    """
    # Constants for starpv calculations
    PXMIN = 1e-7  # Smallest allowed parallax
    VMAX  = 0.5  # Largest allowed speed (fraction of c)
    IMAX  = 100  # Maximum iterations for relativistic solution
    
    # Distance (au)
    if px >= PXMIN:
        parallax_used = px
    else:
        parallax_used = PXMIN
    
    r = DR2AS / parallax_used

    # Radial speed (au/day)
    rd = DAYSEC * rv * 1e3 / DAU

    # Proper motion (radian/day)
    rad  = pmr / DJY
    decd = pmd / DJY

    # To pv-vector (au, au/day) using pymS2pv
    pv = pymS2pv(ra, dec, r, rad, decd, rd)

    # If excessive velocity, arbitrarily set it to zero
    v = np.linalg.norm(pv[1])
    if v / DC > VMAX:
        pv[1] = np.zeros(3)

    # Apply relativistic correction
    # Isolate radial component
    r_pos = np.linalg.norm(pv[0])
    if r_pos == 0.0:
        pu = np.zeros(3)
    else:
        pu = pv[0] / r_pos
        
    vsr = np.dot(pu, pv[1])
    usr = vsr * pu
    
    # Isolate transverse component
    ust = pv[1] - usr
    vst = np.linalg.norm(ust)
    
    # Special-relativity dimensionless parameters
    betsr = vsr / DC
    betst = vst / DC
    
    # Determine observed-to-inertial correction terms iteratively
    bett = betst
    betr = betsr
    
    # Initialize convergence tracking variables
    odd, oddel, od, odel = 0.0, 0.0, 0.0, 0.0
    conv = False
    
    for i in range(IMAX):
        d = 1.0 + betr
        w = betr * betr + bett * bett
        
        # Check for invalid sqrt argument
        if w >= 1.0:
            raise ValueError("Relativistic correction failed: velocity too high")
            
        del_val = -w / (np.sqrt(1.0 - w) + 1.0)
        
        new_betr = d * betsr + del_val
        new_bett = d * betst
        
        if i > 0:
            dd = abs(d - od)
            ddel = abs(del_val - odel)
            if (i > 1) and (dd >= odd) and (ddel >= oddel):
                conv = True
                break
            odd, oddel = dd, ddel
        
        od, odel = d, del_val
        betr, bett = new_betr, new_bett
    
    if not conv:
        raise ValueError("Relativistic correction failed to converge")

    # Scale observed tangential velocity vector into inertial (au/d)
    ut = d * ust

    # Compute inertial radial velocity vector (au/d)
    ur = DC * (d * betsr + del_val) * pu

    # Combine the two to obtain the inertial space velocity vector
    pv[1] = ur + ut
    
    return pv