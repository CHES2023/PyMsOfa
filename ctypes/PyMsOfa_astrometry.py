# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 02:44:08 2023
Done    on Mon Jun  5 02:11:37 2023
@author: Dr. Jianghui JI  (jijh@pmo.ac.cn)
"""
#from   numba  import jit
#from __future__ import print_function
import ctypes
from   ctypes import *
import numpy  as  np
import numpy.ctypeslib  as  nt
import warnings  as  ws 
import os
import sys
import platform as pf

#dll = ctypes.cdll.LoadLibrary
#lib = dll('./libsofa_c.so')

#compile sofa_c.c to create libsofa_c.so shared libary
#gcc -shared -fPIC -o libsofa_c.so sofa_a.c   

#path = os.getcwd()
#lib  = CDLL(os.path.join(path, 'libsofa_c.so'))


'''
   Descr: return file name/path information
#      Input: 
          filepath: full path of a file

#      Return values:
          dirname: directory name of a file
          fullname: full name of a file
          filename: file name of a file
          extension: file extension
'''
def getFileName(filepath):

    __file__  = filepath
    
    dirname   = os.path.dirname(__file__)
    fullname  = os.path.split(__file__)[-1]
    filename  = os.path.split(__file__)[-1].split('.')[0]
    extension = os.path.split(__file__)[-1].split('.')[1]
    
    return dirname, fullname, filename, extension

     
def getFilePath(filepath):

    __file__  = filepath
    
    dirname   = os.path.dirname(__file__)
    fullname  = os.path.split(__file__)[-1]
    filename  = os.path.split(__file__)[-1].split('.')[0]
    extension = os.path.split(__file__)[-1].split('.')[1]
    
    return dirname 

#place libsofa_c.so in the same directory with PyMsOfa.py
cur_dir    = sys.argv[0]
#dirname, fullname, filename, extension = getFileName(cur_dir)
dirname    = getFilePath(cur_dir)
libs_file  = 'libsofa_c.so'
libs_path  = os.path.join(dirname, libs_file) 
#lib       = CDLL(libs_path)


#'''
#for windows
if   pf.system().lower() == 'windows':
           
          #user_path  = './'
          #libs_file  = 'libsofa_c.so'
          #libs_path  = user_path + libs_file
          #print(libs_path)
          #lib  = CDLL(libs_path)
           
          #print(libs_path)
          lib  = CDLL(libs_path)
           
#for linux 
elif pf.system().lower() == 'linux':  
          dirname = os.getcwd()
          lib     = CDLL(os.path.join(dirname, libs_file)) 
          #lib = CDLL('./libsofa_c.so')
          #print(libs_path)
          #lib  = CDLL(libs_path)
 
#for macOS
elif pf.system().lower() == 'darwin':  

          lib  = CDLL(libs_path)

#'''

c_int4 = c_int * 4

array_1d_double = nt.ndpointer(shape=(1,3), dtype=np.double, flags='C') 
vector_double   = nt.ndpointer(shape=(1,3), dtype=np.double, flags='C') 
array_double_2  = nt.ndpointer(shape=(1,2), dtype=np.double, flags='C') 
c_double_p = POINTER(c_double)   

#void iauAb(double pnat[3], double v[3], double s, double bm1,
#           double ppr[3])    
def pymAb(pnat, v, s, bm1):
    '''
    Apply aberration to transform natural direction into proper direction.

    Parameters
    ----------
    pnat : numpy.matrix(1,3)
        natural direction to the source (unit vector)
    v : numpy.matrix(1,3)
        observer barycentric velocity in units of c
    s : float
        distance between the Sun and the observer (au)
    bm1 : float
        sqrt(1-|v|^2): reciprocal of Lorenz factor

    Returns
    -------
    ppr : numpy.matrix(1,3)
        proper direction to source (unit vector)

    '''
    lib.iauAb.argtypes = [vector_double, vector_double, c_double, c_double,
                          vector_double]
    lib.iauAb.restype  = None
    ppr = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    lib.iauAb(pnat, v, s, bm1, ppr)
    
    return ppr

#void iauAe2hd (double az, double el, double phi,
#               double *ha, double *dec)

def pymAe2hd(az, el, phi):
    '''
    Horizon to equatorial coordinates:  transform azimuth and altitude
    to hour angle and declination.

    Parameters
    ----------
    az : float
        azimuth
    el : float
        altitude (informally, elevation)
    phi : float
        site latitude

    Returns
    -------
    ha : float
        hour angle (local)
    dec : float
        declination

    '''
    lib.iauAe2hd.argtypes = [c_double, c_double, c_double]
    lib.iauAe2hd.restype  = None
    ha  = c_double()
    dec = c_double()
    lib.iauAe2hd(az, el, phi, byref(ha),  byref(dec))
    
    return ha.value, dec.value 

'''
/* Star-independent astrometry parameters */
typedef struct {
   double pmt;        /* PM time interval (SSB, Julian years) */
   double eb[3];      /* SSB to observer (vector, au) */
   double eh[3];      /* Sun to observer (unit vector) */
   double em;         /* distance from Sun to observer (au) */
   double v[3];       /* barycentric observer velocity (vector, c) */
   double bm1;        /* sqrt(1-|v|^2): reciprocal of Lorenz factor */
   double bpn[3][3];  /* bias-precession-nutation matrix */
   double along;      /* longitude + s' + dERA(DUT) (radians) */
   double phi;        /* geodetic latitude (radians) */
   double xpl;        /* polar motion xp wrt local meridian (radians) */
   double ypl;        /* polar motion yp wrt local meridian (radians) */
   double sphi;       /* sine of geodetic latitude */
   double cphi;       /* cosine of geodetic latitude */
   double diurab;     /* magnitude of diurnal aberration vector */
   double eral;       /* "local" Earth rotation angle (radians) */
   double refa;       /* refraction constant A (radians) */
   double refb;       /* refraction constant B (radians) */
} pymASTROM;
/* (Vectors eb, eh, em and v are all with respect to BCRS axes.) */

class Point(ctypes.Structure):
    _fields_ = [("x", ctypes.c_int),
                ("y", ctypes.c_int)]
'''

#fortran version for array [3,2]
vector_double2f = nt.ndpointer(shape=(3,2), dtype=np.double, flags='C')
#C  version for array [2,3]
vector_double2  = nt.ndpointer(shape=(2,3), dtype=np.double, flags='C')

vector_double3  = nt.ndpointer(shape=(3,3), dtype=np.double, flags='C')

c_dbl3          = c_double*3
c_dbl33         = c_dbl3*3
c_dbl32         = c_dbl3*2

class pymASTROM(ctypes.Structure):
    _fields_ = [("pmt",    c_double),
                ("eb",     c_dbl3),
                ("eh",     c_dbl3),
                ("em",     c_double),
                ("v",      c_dbl3),
                ("bm1",    c_double),
                ("bpn",    c_dbl33),
                ("along",  c_double),
                ("phi",    c_double),
                ("xpl",    c_double),
                ("ypl",    c_double),
                ("sphi",   c_double),
                ("cphi",   c_double),
                ("diurab", c_double),
                ("eral",   c_double),
                ("refa",   c_double),
                ("refb",   c_double),
                ]
    
   
    def print(pymASTROM):
        #eb_c = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
        #eb_c = (c_dbl3)() 
        #eh_c = (c_dbl3)() 
        #v_c  = (c_dbl3)()
        #bpn_c= np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
        
        eb_c = tuple([x for x in pymASTROM.eb])
        eh_c = tuple([x for x in pymASTROM.eh])
        v_c  = tuple([x for x in pymASTROM.v])
        #v_c  = np.array([x for x in pymASTROM.v])
        #bpn_c= tuple([x for x in pymASTROM.bpn])
        bpn_c  = np.array([x for x in pymASTROM.bpn])
        print("pmt:",    pymASTROM.pmt)
        print("eb:",     eb_c)
        print("eh:",     eh_c)
        print("em:",     pymASTROM.em) 
        print("v:",      v_c)
        print("bm1:",    pymASTROM.bm1) 
        print("bpn:",    bpn_c)
        print("along:",  pymASTROM.along)   
        print("phi:",    pymASTROM.phi)     
        print("xpl:",    pymASTROM.xpl)    
        print("ypl:",    pymASTROM.ypl)    
        print("sphi:",   pymASTROM.sphi)    
        print("cphi:",   pymASTROM.cphi)    
        print("diurab:", pymASTROM.diurab) 
        print("eral:",   pymASTROM.eral)
        print("refa:",   pymASTROM.refa)  
        print("refb:",   pymASTROM.refb)
         

'''
class pymASTROM(ctypes.Structure):
    _fields_ = [("pmt",    c_double),
                ("eb",     c_dbl3),
                ("eh",     c_dbl3),
                ("em",     c_double),
                ("v",      c_dbl3),
                ("bm1",    c_double),
                ("bpn",    c_dbl33),
                ("along",  c_double),
                ("phi",    c_double),
                ("xpl",    c_double),
                ("ypl",    c_double),
                ("sphi",   c_double),
                ("cphi",   c_double),
                ("diurab", c_double),
                ("eral",   c_double),
                ("refa",   c_double),
                ("refb",   c_double),
                ]
    
    #def __init__(self, pmt, eb, eh, em, v, bm1, bpn, along, phi,
    #                   xpl, ypl, sphi, cphi, diurab, eral, refa, refb):
    def __init__(self, pymASTROM):
    
        self.pmt   = pymASTROM.pmt      
        self.eb    = pymASTROM.eb       
        self.eh    = pymASTROM.eh       
        self.em    = pymASTROM.em       
        self.v     = pymASTROM.v        
        self.bm1   = pymASTROM.bm1      
        self.bpn   = pymASTROM.bpn      
        self.along = pymASTROM.along    
        self.phi   = pymASTROM.phi      
        self.xpl   = pymASTROM.xpl      
        self.ypl   = pymASTROM.ypl      
        self.sphi  = pymASTROM.sphi     
        self.cphi  = pymASTROM.cphi     
        self.diurab= pymASTROM.diurab   
        self.eral  = pymASTROM.eral     
        self.refa  = pymASTROM.refa     
        self.refb  = pymASTROM.refb     

        
    def print(self):    
    #def print(pymASTROM):
        #eb_c = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
        #eb_c = (c_dbl3)() 
        #eh_c = (c_dbl3)() 
        #v_c  = (c_dbl3)()
        #bpn_c= np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
        
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
'''

c_star_p =  POINTER(pymASTROM)

#void iauApcg(double date1, double date2,
#             double ebpv[2][3], double ehp[3],
#             pymASTROM *astrom)    

def pymApcg(date1, date2, ebpv,  ehp):
    '''
    For a geocentric observer, prepare star-independent astrometry
    parameters for transformations between ICRS and GCRS coordinates.
    The Earth ephemeris is supplied by the caller.
    
    The parameters produced by this function are required in the
    parallax, light deflection and aberration parts of the astrometric
    transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    
    ebpv : numpy.matrix(2,3)
        Earth barycentric pos/vel (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged
    '''
    lib.iauApcg.argtypes = [c_double, c_double, vector_double2, vector_double,
                            c_star_p]
    lib.iauApcg.restype  = None

#Note: no pointer here, use byref to return structure/class    
    astrom  = pymASTROM()

    lib.iauApcg(date1, date2, ebpv,  ehp,  astrom)
    return  astrom

#void iauApcg13(double date1, double date2, pymASTROM *astrom)
def pymApcg13(date1, date2):
    '''
    For a geocentric observer, prepare star-independent astrometry
    parameters for transformations between ICRS and GCRS coordinates.
    The caller supplies the date, and SOFA models are used to predict
    the Earth ephemeris.
    
    The parameters produced by this function are required in the
    parallax, light deflection and aberration parts of the astrometric
    transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauApcg13.argtypes = [c_double, c_double,   c_star_p]
    lib.iauApcg13.restype  = None
  
    astrom  = pymASTROM()
   
    lib.iauApcg13(date1, date2, astrom)
    return  astrom


#void iauApci(double date1, double date2,
#             double ebpv[2][3], double ehp[3],
#             double x, double y, double s,
#             pymASTROM *astrom)
def pymApci(date1, date2, ebpv,  ehp,  x, y, s): 
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and geocentric CIRS
    coordinates.  The Earth ephemeris and CIP/CIO are supplied by the
    caller.
   
    The parameters produced by this function are required in the
    parallax, light deflection, aberration, and bias-precession-nutation
    parts of the astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date     
    ebpv : numpy.matrix(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)    
    x : flaot
        CIP X,Y (components of unit vector)    
    y : float
        CIP X,Y (components of unit vector)    
    s : float
        the CIO locator s (radians)

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauApci.argtypes = [c_double, c_double, vector_double2, vector_double,
                            c_double, c_double, c_double, 
                            c_star_p]
    lib.iauApci.restype  = None
 
    astrom  = pymASTROM()
    lib.iauApci(date1, date2, ebpv,  ehp,  x,  y,  s,  astrom)
    return  astrom

#void iauApci13(double date1, double date2,
#               pymASTROM *astrom, double *eo)    

def pymApci13(date1, date2):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and geocentric CIRS
    coordinates.  The caller supplies the date, and SOFA models are used
    to predict the Earth ephemeris and CIP/CIO.
    
    The parameters produced by this function are required in the
    parallax, light deflection, aberration, and bias-precession-nutation
    parts of the astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date     

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
    eo : float
        equation of the origins (ERA-GST)
        astrom:
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged
    
    '''
    lib.iauApci13.argtypes = [c_double, c_double, 
                              c_star_p, c_double_p]
    lib.iauApci13.restype  = None
 
    astrom  = pymASTROM()
    eo      = c_double()
    lib.iauApci13(date1, date2, astrom, byref(eo))
    return  astrom, eo.value 

#void iauApco(double date1, double date2,
#             double ebpv[2][3], double ehp[3],
#             double x, double y, double s, double theta,
#             double elong, double phi, double hm,
#             double xp, double yp, double sp,
#             double refa, double refb,
#             pymASTROM *astrom)   
def pymApco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, 
            xp, yp, sp, refa, refb): 
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and observed
    coordinates.  The caller supplies the Earth ephemeris, the Earth
    rotation information and the refraction constants as well as the
    site coordinates.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date     
    ebpv : numpy.matrix(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)    
    x : flaot
        CIP X,Y (components of unit vector)    
    y : float
        CIP X,Y (components of unit vector)    
    s : float
        the CIO locator s (radians)    
    theta : float
        Earth rotation angle (radians)
    elong : float
        longitude (radians, east +ve)
    phi : float
        latitude (geodetic, radians)
    hm : float
        height above ellipsoid (m, geodetic)
    xp : float
        polar motion coordinates (radians)
    yp : float
        polar motion coordinates (radians)
    sp : float
        the TIO locator s' (radians)
    refa : float
        refraction constant A (radians)
    refb : float
        refraction constant B (radians)

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    lib.iauApco.argtypes = [c_double, c_double, vector_double2, vector_double,
                            c_double, c_double, c_double, c_double,  
                            c_double, c_double, c_double,
                            c_double, c_double, c_double,
                            c_double, c_double, 
                            c_star_p]
    lib.iauApco.restype  = None
 
    astrom  = pymASTROM()
    
    lib.iauApco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, 
            xp, yp, sp, refa, refb, astrom)
    return  astrom

#int iauApco13(double utc1, double utc2, double dut1,
#              double elong, double phi, double hm, double xp, double yp,
#              double phpa, double tc, double rh, double wl,
#              pymASTROM *astrom, double *eo)
def pymApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and observed
    coordinates.  The caller supplies UTC, site coordinates, ambient air
    conditions and observing wavelength, and SOFA models are used to
    obtain the Earth ephemeris, CIP/CIO and refraction constants.
   
    The parameters produced by this function are required in the
    parallax, light deflection, aberration, and bias-precession-nutation
    parts of the ICRS/CIRS transformations.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date    
    utc2 : float
        UTC as a 2-part quasi Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float 
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    
        
    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    
    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
    eo : float
        equation of the origins (ERA-GST)
        astrom:
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    lib.iauApco13.argtypes = [c_double, c_double, c_double,
                              c_double, c_double, c_double, c_double, c_double, 
                              c_double, c_double, c_double, c_double,
                              c_star_p, c_double_p]
    lib.iauApco13.restype  = c_int
 
    astrom  = pymASTROM()
    eo      = c_double()
    
    j = lib.iauApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl, astrom, byref(eo))
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
    return  astrom, eo.value 

#void iauApcs(double date1, double date2, double pv[2][3],
#             double ebpv[2][3], double ehp[3],
#             pymASTROM *astrom)

def pymApcs(date1, date2, pv, ebpv, ehp): 
    '''
    For an observer whose geocentric position and velocity are known,
    prepare star-independent astrometry parameters for transformations
    between ICRS and GCRS.  The Earth ephemeris is supplied by the
    caller.
    
    The parameters produced by this function are required in the space
    motion, parallax, light deflection and aberration parts of the
    astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    
    pv : numpy.matrix(2,3)
        observer's geocentric pos/vel (m, m/s)    
    ebpv : numpy.matrix(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)    

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauApcs.argtypes = [c_double, c_double, vector_double2, 
                            vector_double2, vector_double,
                            c_star_p]
    lib.iauApcs.restype  = None
 
    astrom  = pymASTROM()
    
    lib.iauApcs(date1, date2, pv, ebpv, ehp, astrom)
    return  astrom  


#void iauApcs13(double date1, double date2, double pv[2][3],
#               pymASTROM *astrom)  
def pymApcs13(date1, date2, pv): 
    '''
    For an observer whose geocentric position and velocity are known,
    prepare star-independent astrometry parameters for transformations
    between ICRS and GCRS.  The Earth ephemeris is from SOFA models.
    
    The parameters produced by this function are required in the space
    motion, parallax, light deflection and aberration parts of the
    astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    
    pv : numpy.matrix(2,3)
        observer's geocentric pos/vel (m, m/s)    

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauApcs13.argtypes = [c_double, c_double, vector_double2, 
                              c_star_p]
    lib.iauApcs13.restype  = None
 
    astrom  = pymASTROM()
    
    lib.iauApcs13(date1, date2, pv, astrom)
    return  astrom 

#void iauAper(double theta, pymASTROM *astrom)  
def pymAper(theta, astrom):
    '''
    In the star-independent astrometry parameters, update only the
    Earth rotation angle, supplied by the caller explicitly.

    Parameters
    ----------
    theta : float
        Earth rotation angle (radians)    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : not used
        >eb : numpy.matrix(1,3) : not used
        >eh : numpy.matrix(1,3) : not used
        >em : float : not used
        >v : numpy.matrix(1,3) : not used
        >bm1 : float : not used
        >bpn : numpy.matrix(3,3) : not used
        >along : float : longitude + s' (radians)
        >xpl : float : not used
        >ypl : float : not used
        >sphi : float : not used
        >cphi : float : not used
        >diurab : float : not used
        >eral : float : not used
        >refa : float : not used
        >refb : float : not used

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauAper.argtypes = [c_double, c_star_p]
    lib.iauAper.restype  = None

    lib.iauAper(theta, astrom)
    
    return   

#void iauAper13(double ut11, double ut12, pymASTROM *astrom)  
def pymAper13(ut11, ut12, astrom):
    '''
    In the star-independent astrometry parameters, update only the
    Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date    
    ut12 : float
        UT1 as a 2-part Julian Date    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : not used
        >eb : numpy.matrix(1,3) : not used
        >eh : numpy.matrix(1,3) : not used
        >em : float : not used
        >v : numpy.matrix(1,3) : not used
        >bm1 : float : not used
        >bpn : numpy.matrix(3,3) : not used
        >along : float : longitude + s' (radians)
        >xpl : float : not used
        >ypl : float : not used
        >sphi : float : not used
        >cphi : float : not used
        >diurab : float : not used
        >eral : float : not used
        >refa : float : not used
        >refb : float : not used

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauAper13.argtypes = [c_double, c_double, c_star_p]
    lib.iauAper13.restype  = None

    lib.iauAper13(ut11, ut12, astrom)
    
    return  astrom

#void iauApio(double sp, double theta,
#             double elong, double phi, double hm, double xp, double yp,
#             double refa, double refb,
#             pymASTROM *astrom)    

def pymApio(sp, theta, elong, phi,  hm,  xp,  yp,  refa,  refb):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between CIRS and observed
    coordinates.  The caller supplies the Earth orientation information
    and the refraction constants as well as the site coordinates.

    Parameters
    ----------
    sp : float
        the TIO locator s' (radians)    
    theta : float
        Earth rotation angle (radians)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    refa : float
        refraction constant A (radians)    
    refb : float
        refraction constant B (radians)

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    lib.iauApio.argtypes = [c_double, c_double, c_double, c_double, c_double,
                            c_double, c_double, c_double, c_double,
                            c_star_p]
    lib.iauApio.restype  = None
    
    astrom  = pymASTROM()
    lib.iauApio(sp, theta, elong, phi, hm, xp,  yp,  refa,  refb, astrom)
    
    return   astrom  

#int iauApio13(double utc1, double utc2, double dut1,
#              double elong, double phi, double hm, double xp, double yp,
#              double phpa, double tc, double rh, double wl,
#              pymASTROM *astrom)
    
def pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between CIRS and observed
    coordinates.  The caller supplies UTC, site coordinates, ambient air
    conditions and observing wavelength.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    lib.iauApio13.argtypes = [c_double, c_double, c_double,
                              c_double, c_double, c_double, c_double, c_double, 
                              c_double, c_double, c_double, c_double,
                              c_star_p]
    lib.iauApio13.restype  =  c_int
 
    astrom  = pymASTROM()
      
    j = lib.iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                     phpa, tc, rh,  wl, astrom)
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
    return  astrom

#void iauAtcc13(double rc, double dc,
#               double pr, double pd, double px, double rv,
#               double date1, double date2,
#               double *ra, double *da)

def pymAtcc13(rc,  dc,  pr,  pd,  px,  rv, date1, date2):
    '''
    Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS
    astrometric place.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    ra : float
        ICRS astrometric RA (radians)    
    da : float
        ICRS astrometric Dec (radians)

    '''
    lib.iauAtcc13.argtypes = [c_double, c_double, 
                              c_double, c_double, c_double, c_double,
                              c_double, c_double, 
                              c_double_p, c_double_p]
    lib.iauAtcc13.restype  = None
    
    ra  = c_double()
    da  = c_double()
    
    lib.iauAtcc13(rc,  dc,  pr,  pd,  px,  rv, date1, date2, 
                  byref(ra), byref(da))
    
    return  ra.value, da.value   

#void iauAtccq(double rc, double dc,
#              double pr, double pd, double px, double rv,
#              pymASTROM *astrom, double *ra, double *da) 
def pymAtccq(rc,  dc,  pr,  pd,  px,  rv,  astrom): 
    '''
    Quick transformation of a star's ICRS catalog entry (epoch J2000.0)
    into ICRS astrometric place, given precomputed star-independent
    astrometry parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    If the parallax and proper motions are zero the transformation has
    no effect.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ra : float
        ICRS astrometric RA (radians)    
    da : float
        ICRS astrometric Dec (radians)

    '''
    lib.iauAtccq.argtypes = [c_double, c_double,  
                             c_double, c_double, c_double, c_double,
                             c_star_p, c_double_p, c_double_p]
    lib.iauAtccq.restype  =  None
 
    ra  = c_double()
    da  = c_double()
      
    lib.iauAtccq(rc,  dc,  pr,  pd,  px,  rv,  astrom, 
                     byref(ra), byref(da)) 
    return  ra.value, da.value 

#void iauAtci13(double rc, double dc,
#               double pr, double pd, double px, double rv,
#               double date1, double date2,
#               double *ri, double *di, double *eo)   
def pymAtci13(rc,  dc,  pr,  pd,  px,  rv,  date1, date2):
    '''
    Transform ICRS star data, epoch J2000.0, to CIRS.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    
    eo : float
        equation of the origins (ERA-GST)

    '''
    lib.iauAtci13.argtypes = [c_double, c_double, 
                              c_double, c_double, c_double, c_double,
                              c_double, c_double, 
                              c_double_p, c_double_p, c_double_p]
    lib.iauAtci13.restype  = None
    
    ri  = c_double()
    di  = c_double()
    eo  = c_double()
    lib.iauAtci13(rc,  dc,  pr,  pd,  px,  rv,  date1, date2, 
                  byref(ri), byref(di), byref(eo))
    
    return  ri.value, di.value, eo.value    

#void iauAtciq(double rc, double dc,
#              double pr, double pd, double px, double rv,
#              pymASTROM *astrom, double *ri, double *di)
def pymAtciq(rc,  dc,  pr,  pd,  px,  rv,  astrom):
    '''
    Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    star-independent astrometry parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    If the parallax and proper motions are zero the iauAtciqz function
    can be used instead.
    
    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    

    '''
    lib.iauAtciq.argtypes = [c_double, c_double, 
                             c_double, c_double, c_double, c_double,
                             c_star_p, c_double_p, c_double_p]
    lib.iauAtciq.restype  = None
    
    ri  = c_double()
    di  = c_double()
     
    lib.iauAtciq(rc,  dc,  pr,  pd,  px,  rv,  astrom, 
                 byref(ri), byref(di)) 
                 
    return  ri.value, di.value 


'''
                                                                       
/* Body parameters for light deflection */                             
typedef struct {                                                       
   double bm;         /* mass of the body (solar masses) */            
   double dl;         /* deflection limiter (radians^2/2) */           
   double pv[2][3];   /* barycentric PV of the body (au, au/day) */    
} iauLDBODY;                                                           
                                                                    
'''
#2023-03-12, iauLDBODY class from C structure 
class iauLDBODY(ctypes.Structure):
    _fields_ = [("bm:",    c_double),
                ("dl:",    c_double),
                ("pv:",    c_dbl32),
                ]
    
    def __init__(self, bm, dl, pv):
   
        self.bm  = bm
        self.dl  = dl
        self.pv  = pv 
        
    #@property     
    def print(iauLDBODY):
         pv_c  = np.array([x for x in iauLDBODY.pv])
         print("bm:", iauLDBODY.bm)
         print("dl:", iauLDBODY.dl)
         print("pv:", pv_c)
    
  
#void iauAtciqn(double rc, double dc, double pr, double pd,
#               double px, double rv, pymASTROM *astrom,
#               int n, iauLDBODY b[], double *ri, double *di)

#using 2-dim array to replace iauLDBody structure, it works
c_2d_double = nt.ndpointer(dtype=np.double, ndim=2, flags="C")

def pymAtciqn(rc,  dc,  pr,  pd,  px,  rv,  astrom,  n,  b):
    '''
    Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    star-independent astrometry parameters plus a list of light-
    deflecting bodies.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    If the only light-deflecting body to be taken into account is the
    Sun, the iauAtciq function can be used instead.  If in addition the
    parallax and proper motions are zero, the iauAtciqz function can be
    used.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)
    n : int
        number of bodies    
    b : iauLDBODY class
        data for each of the n bodies    
        >bm : float : mass of the body (solar masses)    
        >dl : float : deflection limiter
        >pv : numpy.matrix(2,3) : barycentric PV of the body (au, au/day)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  

    '''
    lib.iauAtciqn.argtypes = [c_double,   c_double,  c_double,  c_double,
                              c_double,   c_double,  c_star_p, 
                              c_int,      c_2d_double, 
                              c_double_p, c_double_p]
    lib.iauAtciqn.restype  = None
    
    ri  = c_double()
    di  = c_double()
     
         
    lib.iauAtciqn(rc,  dc,  pr,  pd,  px,  rv,  astrom,  n,  b,
                  byref(ri), byref(di))
    
    return  ri.value, di.value

#void iauAtciqz(double rc, double dc, pymASTROM *astrom,
#               double *ri, double *di)
#2023-03-13    
def pymAtciqz(rc,  dc,  astrom): 
    '''
    Quick ICRS to CIRS transformation, given precomputed star-
    independent astrometry parameters, and assuming zero parallax and
    proper motion.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    The corresponding function for the case of non-zero parallax and
    proper motion is iauAtciq.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  

    '''
    lib.iauAtciqz.argtypes = [c_double, c_double,  
                              c_star_p, c_double_p, c_double_p]
    lib.iauAtciqz.restype  =  None
 
    ri  = c_double()
    di  = c_double()
      
    lib.iauAtciqz(rc,  dc,  astrom,  byref(ri), byref(di)) 
    return  ri.value, di.value 

#int iauAtco13(double rc, double dc,
#              double pr, double pd, double px, double rv,
#              double utc1, double utc2, double dut1,
#              double elong, double phi, double hm, double xp, double yp,
#              double phpa, double tc, double rh, double wl,
#              double *aob, double *zob, double *hob,
#              double *dob, double *rob, double *eo) 

def pymAtco13(rc,  dc,  pr,  pd,  px,  rv,
              utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
    '''
    ICRS RA,Dec to observed place.  The caller supplies UTC, site
    coordinates, ambient air conditions and observing wavelength.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    aob : float
        observed azimuth (radians: N=0,E=90)     
    zob : float
        observed zenith distance (radians)    
    hob : float
        observed hour angle (radians)    
    dob : float
        observed declination (radians)    
    rob : float
        observed right ascension (CIO-based, radians)    
    eo : float
        equation of the origins (ERA-GST)

    '''
    lib.iauAtco13.argtypes = [c_double, c_double, 
                              c_double, c_double, c_double, c_double,
                              c_double, c_double, c_double, 
                              c_double, c_double, c_double, c_double, c_double,
                              c_double, c_double, c_double, c_double,
                              c_double_p, c_double_p, c_double_p, 
                              c_double_p, c_double_p, c_double_p ]
    lib.iauAtco13.restype  = c_int
 
    aob     = c_double()
    zob     = c_double()
    hob     = c_double()
    dob     = c_double()
    rob     = c_double()
    eo      = c_double()
   
    j = lib.iauAtco13(rc,  dc,  pr,  pd,  px,  rv,
                      utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      byref(aob), byref(zob), byref(hob), 
                      byref(dob), byref(rob),byref(eo))
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
    return  aob.value, zob.value, hob.value, dob.value, rob.value, eo.value   

#2023-03-15
#void iauAtic13(double ri, double di, double date1, double date2,
#               double *rc, double *dc, double *eo)  
def pymAtic13(ri,  di,  date1,  date2):
    '''
    Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    eo : float
        equation of the origins (ERA-GST)

    '''
    lib.iauAtic13.argtypes = [c_double, c_double, c_double, c_double,
                              c_double_p, c_double_p, c_double_p]
    lib.iauAtic13.restype  = None
    
    rc  = c_double()
    dc  = c_double()
    eo  = c_double()
    lib.iauAtic13(ri,  di,  date1,  date2,  byref(rc), byref(dc), byref(eo))
    
    return  rc.value, dc.value, eo.value 


#void iauAticq(double ri, double di, pymASTROM *astrom,
#              double *rc, double *dc)      
def pymAticq(ri,  di,  astrom): 
    '''
    Quick CIRS RA,Dec to ICRS astrometric place, given the star-
    independent astrometry parameters.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    lib.iauAticq.argtypes = [c_double, c_double,  c_star_p,
                             c_double_p, c_double_p]
    lib.iauAticq.restype  =  None
 
    rc  = c_double()
    dc  = c_double()
      
    lib.iauAticq(ri,  di,  astrom,  byref(rc), byref(dc)) 
    return  rc.value, dc.value    

#2023-04-16
#it now work well for 2-dim of b rather than the structure iauLDBody      
def pymAticqn(ri,  di,  astrom,  n,  b):
    '''
    Quick CIRS to ICRS astrometric place transformation, given the star-
    independent astrometry parameters plus a list of light-deflecting
    bodies.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
    or iauApcs[13].
    
    If the only light-deflecting body to be taken into account is the
    Sun, the iauAticq function can be used instead.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)
    n : int
        number of bodies    
    b : iauLDBODY class
        data for each of the n bodies    
        >bm : float : mass of the body (solar masses)    
        >dl : float : deflection limiter
        >pv : numpy.matrix(2,3) : barycentric PV of the body (au, au/day)

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    lib.iauAticqn.argtypes = [c_double,   c_double,   c_star_p,
                              c_int,      c_2d_double, 
                              c_double_p, c_double_p]
    lib.iauAticqn.restype  = None
    
    rc  = c_double()
    dc  = c_double()
     
    lib.iauAticqn(ri,  di,  astrom,  n,  b,
                  byref(rc), byref(dc))
    
    return  rc.value, dc.value  

def pymAtio13(ri,  di, 
              utc1, utc2, dut1, 
              elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
    '''
    CIRS RA,Dec to observed place.  The caller supplies UTC, site
    coordinates, ambient air conditions and observing wavelength.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    aob : float
        observed azimuth (radians: N=0,E=90)     
    zob : float
        observed zenith distance (radians)    
    hob : float
        observed hour angle (radians)    
    dob : float
        observed declination (radians)    
    rob : float
        observed right ascension (CIO-based, radians)   
    '''
    
    lib.iauAtio13.argtypes = [c_double, c_double, 
                              c_double, c_double, c_double, 
                              c_double, c_double, c_double, c_double, c_double,
                              c_double, c_double, c_double, c_double,
                              c_double_p, c_double_p, c_double_p, 
                              c_double_p, c_double_p ]
    lib.iauAtio13.restype  = c_int
 
    aob     = c_double()
    zob     = c_double()
    hob     = c_double()
    dob     = c_double()
    rob     = c_double()
   
    j = lib.iauAtio13(ri,  di,   
                      utc1, utc2, dut1, 
                      elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      byref(aob), byref(zob), byref(hob), 
                      byref(dob), byref(rob) )
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
    return  aob.value, zob.value, hob.value, dob.value, rob.value

def pymAtioq(ri,  di, astrom):
    '''
    Quick CIRS to observed place transformation.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling iauApio[13] or iauApco[13].

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    aob : float
        observed azimuth (radians: N=0,E=90)     
    zob : float
        observed zenith distance (radians)    
    hob : float
        observed hour angle (radians)    
    dob : float
        observed declination (radians)    
    rob : float
        observed right ascension (CIO-based, radians)   

    '''
    lib.iauAtioq.argtypes = [ c_double, c_double, c_star_p,
                              c_double_p, c_double_p, c_double_p, 
                              c_double_p, c_double_p ]
    lib.iauAtioq.restype  = None
 
    aob     = c_double()
    zob     = c_double()
    hob     = c_double()
    dob     = c_double()
    rob     = c_double()
   
    j = lib.iauAtioq( ri,  di,   astrom,
                      byref(aob), byref(zob), byref(hob), 
                      byref(dob), byref(rob) )
 
    return  aob.value, zob.value, hob.value, dob.value, rob.value

def pymAtoc13(stype,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl): 
    '''
    Observed place at a groundbased site to to ICRS astrometric RA,Dec.
    The caller supplies UTC, site coordinates, ambient air conditions
    and observing wavelength.

    Parameters
    ----------
    stype : ctypes.c_char_p
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)  

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    lib.iauAtoc13.argtypes = [c_char_p, c_double, c_double, 
                              c_double, c_double, c_double, 
                              c_double, c_double, c_double, c_double, c_double,
                              c_double, c_double, c_double, c_double,
                              c_double_p, c_double_p ]
    lib.iauAtoc13.restype  = c_int
 
    rc  = c_double()
    dc  = c_double()
   
    j = lib.iauAtoc13(stype,   ob1,   ob2,
                      utc1,  utc2,  dut1, 
                      elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      byref(rc), byref(dc))
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
        
    return  rc.value, dc.value  

def pymAtoi13(stype,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl): 
    '''
    Observed place to CIRS.  The caller supplies UTC, site coordinates,
    ambient air conditions and observing wavelength.

    Parameters
    ----------
    stype : ctypes.c_char_p
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)  
    
    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    
    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     

    '''
    lib.iauAtoi13.argtypes = [c_char_p, c_double, c_double, 
                              c_double, c_double, c_double, 
                              c_double, c_double, c_double, c_double, c_double,
                              c_double, c_double, c_double, c_double,
                              c_double_p, c_double_p ]
    lib.iauAtoi13.restype  = c_int
 
    ri  = c_double()
    di  = c_double()
   
    j = lib.iauAtoi13(stype,   ob1,   ob2,
                      utc1,  utc2,  dut1, 
                      elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      byref(ri), byref(di))
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
        
    return  ri.value, di.value  

def pymAtoiq(stype,   ob1,   ob2,  astrom):
    '''
    Quick observed place to CIRS, given the star-independent astrometry
    parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling iauApio[13] or iauApco[13].

    Parameters
    ----------
    stype : ctypes.c_char_p
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     

    '''
    lib.iauAtoiq.argtypes = [ c_char_p, c_double, c_double, c_star_p,
                              c_double_p, c_double_p ]
    lib.iauAtoiq.restype  = None
 
    ri  = c_double()
    di  = c_double()
   
    lib.iauAtoiq(stype,   ob1,   ob2,  astrom,
                     byref(ri), byref(di))
        
    return  ri.value, di.value  

#2023-04-17 
def pymHd2ae(ha, dec, phi):
    '''
    Equatorial to horizon coordinates:  transform hour angle and
    declination to azimuth and altitude.

    Parameters
    ----------
    ha : float
        hour angle (local)    
    dec : float
        declination    
    phi : float
        site latitude

    Returns
    -------
    az : float
        azimuth    
    el : float
        altitude (informally, elevation)

    '''
    lib.iauHd2ae.argtypes = [ c_double, c_double, c_double,
                              c_double_p, c_double_p ]
    lib.iauHd2ae.restype  = None
 
    az  = c_double()
    el  = c_double()
   
    lib.iauHd2ae(ha, dec, phi, 
                 byref(az), byref(el))
        
    return  az.value, el.value  

def pymHd2pa(ha, dec, phi):
    '''
    Parallactic angle for a given hour angle and declination.

    Parameters
    ----------
    ha : float
        hour angle    
    dec : float
        declination    
    phi : float
        site latitude

    Returns
    -------
    function value : flaot
        parallactic angle

    '''
    lib.iauHd2pa.argtypes = [c_double, c_double, c_double]
    lib.iauHd2pa.restype  = c_double

    return  lib.iauHd2pa(ha, dec, phi)

#2023-04-18
def pymLd(bm, p, q, e, em, dlim):
    '''
    Apply light deflection by a solar-system body, as part of
    transforming coordinate direction into natural direction.

    Parameters
    ----------
    bm : float
        mass of the gravitating body (solar masses)    
    p : numpy.matrix(1,3)
        direction from observer to source (unit vector)    
    q : numpy.matrix(1,3)
        direction from body to source (unit vector)    
    e : numpy.matrix(1,3)
        direction from body to observer (unit vector)    
    em : flaot
        distance from body to observer (au)    
    dlim : flaot
        deflection limiter

    Returns
    -------
    p1 : numpy.matrix(1,3)
        observer to deflected source (unit vector)

    '''
    lib.iauLd.argtypes =[c_double, vector_double, vector_double, vector_double,
                         c_double, c_double, vector_double]
    
    lib.iauLd.restype  = None
    
    p1 = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauLd(bm, p, q, e, em, dlim, p1) 

    return  p1

#c_2d_double = nt.ndpointer(dtype=np.double, ndim=2, flags="C")

def pymLdn(n, b, ob, sc):
    '''
    For a star, apply light deflection by multiple solar-system bodies,
    as part of transforming coordinate direction into natural direction.

    Parameters
    ----------
    n : int
        number of bodies    
    b : iauLDBODY class
        data for each of the n bodies    
        >bm : float : mass of the body (solar masses)    
        >dl : float : deflection limiter
        >pv : numpy.matrix(2,3) : barycentric PV of the body (au, au/day)
    ob : numpy.matrix(1,3)
        barycentric position of the observer (au)    
    sc : numpy.matrix(1,3)
        observer to star coord direction (unit vector)

    Returns
    -------
    sn : numpy.matrix(1,3)
        observer to deflected star (unit vector)

    '''                                  
    lib.iauLdn.argtypes =[c_int, c_2d_double, vector_double, vector_double,
                          vector_double]
    
    lib.iauLdn.restype  = None
    
    sn = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauLdn(n, b, ob, sc, sn) 

    return  sn

def pymLdsun(p, e, em):
    '''
    Deflection of starlight by the Sun.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        direction from observer to star (unit vector)    
    e : numpy.matrix(1,3)
        direction from Sun to observer (unit vector)    
    em : float
        distance from Sun to observer (au)

    Returns
    -------
    p1 : numpy.matrix(1,3)
        observer to deflected star (unit vector)

    '''
    lib.iauLdsun.argtypes =[vector_double, vector_double, c_double,
                            vector_double]
    
    lib.iauLdsun.restype  = None
    
    p1 = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauLdsun(p, e, em, p1) 

    return  p1

#2023-05-01
def pymPmpx(rc, dc, pr, pd, px, rv, pmt, pob):
    '''
    Proper motion and parallax.

    Parameters
    ----------
    rc : float
        ICRS RA,Dec at catalog epoch (radians)    
    dc : float
        ICRS RA,Dec at catalog epoch (radians)    
    pr : float
        RA proper motion (radians/year)    
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)    
    pmt : float
        proper motion time interval (SSB, Julian years)    
    pob : numpy.matrix(1,3)
        SSB to observer vector (au)    

    Returns
    -------
    pco : numpy.matrix(1,3)
        coordinate direction (BCRS unit vector)

    '''
    lib.iauPmpx.argtypes =[c_double, c_double, c_double, c_double,
                           c_double, c_double, c_double, vector_double,
                           vector_double]
    
    lib.iauPmpx.restype  = None
    
    pco = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPmpx(rc, dc, pr, pd, px, rv, pmt, pob, pco) 

    return  pco

pym_pmsafe_msg = {
                  4: 'solution not converge',
                  2: 'excessive velocity',
                  1: 'distance overridden',
                  0: 'no warnings or errors',
                 -1: 'system error',
                  }  
 
def pymPmsafe(ra1,  dec1, pmr1, pmd1,
              px1,  rv1,
              ep1a, ep1b, ep2a, ep2b):
    '''
    Star proper motion:  update star catalog data for space motion, with
    special handling to handle the zero parallax case.

    Parameters
    ----------
    ra1 : float
        right ascension (radians), before    
    dec1 : float
        declination (radians), before    
    pmr1 : float
        RA proper motion (radians/year), before     
    pmd1 : float
        Dec proper motion (radians/year), before     
    px1 : float
        parallax (arcseconds), before    
    rv1 : float
        radial velocity (km/s, +ve = receding), before    
    ep1a : float
        "before" epoch, part A    
    ep1b : float
        "before" epoch, part B    
    ep2a : float
        "after" epoch, part A     
    ep2b : float
        "after" epoch, part B    

    Raises
    ------
    ValueError
        4: 'solution not converge',
        2: 'excessive velocity',
        1: 'distance overridden',
        0: 'no warnings or errors',
       -1: 'system error',
        
    Returns
    -------
    ra2 : float
        right ascension (radians), after    
    dec2 : float
        declination (radians), after    
    pmr2 : float
        RA proper motion (radians/year), after      
    pmd2 : float
        Dec proper motion (radians/year), after       
    px2 : float
        parallax (arcseconds), after    
    rv2 : float
        radial velocity (km/s, +ve = receding), after    

    '''
    lib.iauPmsafe.argtypes = [c_double, c_double, c_double, c_double,
                              c_double, c_double,   
                              c_double, c_double, c_double, c_double, 
                              c_double_p, c_double_p, c_double_p, c_double_p,
                              c_double_p, c_double_p]
    
    lib.iauPmsafe.restype  = c_int
 
    ra2   = c_double()
    dec2  = c_double()
    pmr2  = c_double()
    pmd2  = c_double()
    
    px2   = c_double()
    rv2   = c_double()
   
    j = lib.iauPmsafe(ra1,  dec1, pmr1, pmd1,
                      px1,  rv1,
                      ep1a, ep1b,  ep2a, ep2b, 
                      byref(ra2),  byref(dec2),  byref(pmr2),  byref(pmd2),
                      byref(px2),  byref(rv2))
   
#debug message should check 
    if   j < 0:
        raise ValueError(pym_pmsafe_msg[j])
    elif j > 0:
         ws.warn(pym_pmsafe_msg[j], UserWarning, 2)
        
    return ra2.value, dec2.value, pmr2.value, pmd2.value, px2.value, rv2.value  

def pymPvtob(elong,  phi,  hm,  xp,  yp,  sp,  theta):
    '''
    Position and velocity of a terrestrial observing station.

    Parameters
    ----------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    hm : float
        height above ref. ellipsoid (geodetic, m)    
    xp : float
        coordinates of the pole (radians)    
    yp : float
        coordinates of the pole (radians)     
    sp : float
        the TIO locator s' (radians)    
    theta : float
        Earth rotation angle (radians)

    Returns
    -------
    pv : numpy.matrix(2,3)
        position/velocity vector (m, m/s, CIRS)

    '''
    lib.iauPvtob.argtypes =[c_double, c_double, c_double, 
                            c_double, c_double, c_double, c_double,
                            vector_double2]
    
    lib.iauPvtob.restype  = None
    
    pv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauPvtob(elong,  phi,  hm,
                 xp,  yp,  sp,  theta,  pv) 

    return  pv

def pymRefco(phpa,   tc,   rh,   wl):
    '''
    Determine the constants A and B in the atmospheric refraction model
    dZ = A tan Z + B tan^3 Z.

    Parameters
    ----------
    phpa : float
        pressure at the observer (hPa = millibar)    
    tc : float
        ambient temperature at the observer (deg C)     
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)

    Returns
    -------
    refa : float
        tan Z coefficient (radians)     
    refb : float
        tan^3 Z coefficient (radians)    

    '''
    lib.iauRefco.argtypes =[c_double, c_double, c_double, c_double, 
                            c_double_p, c_double_p ]
    
    lib.iauRefco.restype  = None
    
    refa   = c_double()
    refb   = c_double()
    
    lib.iauRefco(phpa,   tc,   rh,   wl,
                 byref(refa),   byref(refb)) 

    return  refa.value, refb.value

def pymTpors(xi,   eta,   a,   b):
    '''
    In the tangent plane projection, given the rectangular coordinates
    of a star and its spherical coordinates, determine the spherical
    coordinates of the tangent point.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image    
    a : float
        star's spherical coordinates    
    b : float
        star's spherical coordinates

    Returns
    -------
    a01 : float
        tangent point's spherical coordinates, Soln. 1    
    b01 : float
        tangent point's spherical coordinates, Soln. 1    
    a02 : float
        tangent point's spherical coordinates, Soln. 2    
    b02 : float
        tangent point's spherical coordinates, Soln. 2

    '''
    lib.iauTpors.argtypes =[c_double, c_double, c_double, c_double, 
                            c_double_p, c_double_p, c_double_p, c_double_p]
    
    lib.iauTpors.restype  = c_int
    
    a01   = c_double()
    b01   = c_double()
    a02   = c_double()
    b02   = c_double()
    
    lib.iauTpors(xi,   eta,   a,   b,
                 byref(a01),  byref(b01), byref(a02), byref(b02)) 

    return  a01.value, b01.value, a02.value, b02.value 

#2023-05-02
pym_tporv_msg = {
                  0: 'no solutions returned',
                  1: 'only the first solution is useful',
                  2: 'both solutions are useful',
                 }  
def pymTporv(xi,   eta,   v):
    '''
    In the tangent plane projection, given the rectangular coordinates
    of a star and its direction cosines, determine the direction
    cosines of the tangent point.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image     
    v : numpy.matrix(3,3)
        star's direction cosines

    Raises
    ------
    ValueError
        0: 'no solutions returned',
        1: 'only the first solution is useful',
        2: 'both solutions are useful',

    Returns
    -------
    v01 : numpy.matrix(3,3)
        tangent point's direction cosines, Solution 1    
    v02 : numpy.matrix(3,3)
        tangent point's direction cosines, Solution 2    

    '''
    lib.iauTporv.argtypes =[c_double, c_double, vector_double,  
                            vector_double, vector_double]
    
    lib.iauTporv.restype  = c_int
    
    v01   = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    v02   = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
     
    
    j = lib.iauTporv(xi,   eta,   v,  v01,  v02)
    
    #if   j == 0:
    #     raise ValueError(pym_tporv_msg[j])
    #elif j > 0:
    #     ws.warn(pym_tporv_msg[j], UserWarning, 2)    
     
    return  v01, v02

#2023-05-03
def pymTpsts(xi,   eta,   a0,   b0):
    '''
    In the tangent plane projection, given the star's rectangular
    coordinates and the spherical coordinates of the tangent point,
    solve for the spherical coordinates of the star.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image        
    a0 : float
        tangent point's spherical coordinates    
    b0 : float
        tangent point's spherical coordinates

    Returns
    -------
    a : float
        star's spherical coordinates    
    b : float
        star's spherical coordinates

    '''
    lib.iauTpsts.argtypes =[c_double, c_double, c_double, c_double, 
                            c_double_p, c_double_p]
    
    lib.iauTpsts.restype  = None
    
    a   = c_double()
    b   = c_double()
   
    
    lib.iauTpsts(xi,   eta,   a0,   b0,
                 byref(a),  byref(b)) 

    return  a.value, b.value

def pymTpstv(xi,   eta,   v0):
    '''
    In the tangent plane projection, given the star's rectangular
    coordinates and the direction cosines of the tangent point, solve
    for the direction cosines of the star.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image     
    v0 : numpy.matrix(3,3)
        tangent point's direction cosines

    Returns
    -------
    v : numpy.matrix(3,3)
        star's direction cosines

    '''
    lib.iauTpstv.argtypes =[c_double, c_double, vector_double,  
                            vector_double]
    
    lib.iauTpstv.restype  = None
    
    v = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauTpstv(xi,   eta,   v0,   v)
     
    return  v

pym_tpxes_msg = {
                  0: 'OK',
                  1: 'star too far from axis',
                  2: 'antistar on tangent plane',
                  3: 'antistar too far from axis',
                 }  

def pymTpxes(a,   b,   a0,   b0):
    '''
    In the tangent plane projection, given celestial spherical
    coordinates for a star and the tangent point, solve for the star's
    rectangular coordinates in the tangent plane.

    Parameters
    ----------
    a : float
        star's spherical coordinates    
    b : float
        star's spherical coordinates    
    a0 : float
        tangent point's spherical coordinates    
    b0 : float
        tangent point's spherical coordinates 

    Raises
    ------
    ValueError
        0: 'OK',
        1: 'star too far from axis',
        2: 'antistar on tangent plane',
        3: 'antistar too far from axis',

    Returns
    -------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image

    '''
    lib.iauTpxes.argtypes =[c_double, c_double, c_double, c_double,
                            c_double_p, c_double_p]
    
    lib.iauTpxes.restype  = c_int
    
    xi   = c_double()
    eta  = c_double()
    
    j = lib.iauTpxes(a,   b,   a0,   b0, 
                     byref(xi),   byref(eta))
    
    if  j > 0:
        ws.warn(pym_tpxes_msg[j], UserWarning, 2) 
        
    return  xi.value, eta.value



def pymTpxev(v,   v0):
    '''
    In the tangent plane projection, given celestial direction cosines
    for a star and the tangent point, solve for the star's rectangular
    coordinates in the tangent plane.

    Parameters
    ----------
    v : numpy.matrix(3,3)
        direction cosines of star    
    v0 : numpy.matrix(3,3)
        direction cosines of tangent point

    Returns
    -------
    xi : float
        tangent plane coordinates of star    
    eta : float
        tangent plane coordinates of star

    '''
    lib.iauTpxev.argtypes =[vector_double, vector_double,
                            c_double_p, c_double_p]
    
    lib.iauTpxev.restype  = c_int
    
    xi   = c_double()
    eta  = c_double()
    
    j = lib.iauTpxev(v,  v0,   byref(xi),   byref(eta))
    
    if  j > 0:
        ws.warn(pym_tpxes_msg[j], UserWarning, 2)   
     
    return  xi.value, eta.value


def pymEceq06(date1,   date2,  dl,  db):
    '''
    Transformation from ecliptic coordinates (mean equinox and ecliptic
    of date) to ICRS RA,Dec, using the IAU 2006 precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date     
    dl : float
        ecliptic longitude (radians)
    db : float
        ecliptic latitude (radians)

    Returns
    -------
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    '''
    lib.iauEceq06.argtypes = [c_double,    c_double,   c_double,   c_double,
                              c_double_p,  c_double_p]
    lib.iauEceq06.restype  =  None
    
    dr     = c_double()
    dd     = c_double()
  
    
    lib.iauEceq06(date1,   date2,   dl,   db,  byref(dr),  byref(dd))
    
    return  dr.value,  dd.value 


def pymEcm06(date1, date2):
    '''
    ICRS equatorial to ecliptic rotation matrix, IAU 2006.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date     

    Returns
    -------
    rm : numpy.matrix(3,3)
        ICRS to ecliptic rotation matrix

    '''
    lib.iauEcm06.argtypes = [c_double,  c_double,  vector_double3]
    lib.iauEcm06.restype  = None
    
    rm = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauEcm06(date1, date2, rm)
    
    return  rm 

def pymEform(n):
    '''
    Earth reference ellipsoids.

    Parameters
    ----------
    n : int
        ellipsoid identifier
        
    n    ellipsoid    
    1     WGS84    
    2     GRS80    
    3     WGS72    
    
    Returns
    -------
    a : float
        equatorial radius (meters)    
    f : float
        flattening

    '''
    lib.iauEform.argtypes = [c_int,   c_double_p,  c_double_p]
    lib.iauEform.restype  = c_int
    
    a     = c_double()
    f     = c_double()
    
    lib.iauEform(n,  byref(a),  byref(f))
    
    return  a.value,  f.value

    
def pymEpv00(date1, date2): 
    '''
    Earth position and velocity, heliocentric and barycentric, with
    respect to the Barycentric Celestial Reference System.

    Parameters
    ----------
    date1 : float
        TDB date    
    date2 : float
        TDB date

    Returns
    -------
    pvh : numpy.matrix(2,3)
        heliocentric Earth position/velocity    
    pvb : numpy.matrix(2,3)
        barycentric Earth position/velocity

    '''
    lib.iauEpv00.argtypes = [c_double, c_double, vector_double2, vector_double2]
    lib.iauEpv00.restype  = c_int
 
    pvh = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    pvb = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauEpv00(date1, date2, pvh, pvb)
    
    return  pvh, pvb


def pymEqec06(date1,   date2,  dr,  dd):
    '''
    Transformation from ICRS equatorial coordinates to ecliptic
    coordinates (mean equinox and ecliptic of date) using IAU 2006
    precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian date    
    date2 : float
        TT as a 2-part Julian date    
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    Returns
    -------
    dl : flaot
        ecliptic longitude (radians)    
    db : float
        ecliptic latitude (radians)

    '''
    lib.iauEqec06.argtypes = [c_double,    c_double,   c_double,   c_double,
                              c_double_p,  c_double_p]
    lib.iauEqec06.restype  =  None
    
    dl     = c_double()
    db     = c_double()
  
    
    lib.iauEqec06(date1,   date2,  dr,  dd,  byref(dl),  byref(db))
    
    return  dl.value,  db.value 


def pymFk5hip():
    '''
    FK5 to Hipparcos rotation and spin.

    Returns
    -------
    r5h : numpy.matrix(3,3)
        r-matrix: FK5 rotation wrt Hipparcos    
    s5h : numpy.matrix(1,3)
        r-vector: FK5 spin wrt Hipparcos

    '''
    lib.iauFk5hip.argtypes = [vector_double3,   vector_double]
    
    lib.iauFk5hip.restype  = None
    
    r5h = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    s5h = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauFk5hip(r5h, s5h)
    
    return  r5h, s5h  


def pymFk5hz(r5,   d5,   date1,   date2):
    '''
    Transform an FK5 (J2000.0) star position into the system of the
    Hipparcos catalogue, assuming zero Hipparcos proper motion.

    Parameters
    ----------
    r5 : flaot
        FK5 RA (radians), equinox J2000.0, at date    
    d5 : flaot
        FK5 Dec (radians), equinox J2000.0, at date     
    date1 : flaot
        TDB date    
    date2 : flaot
        TDB date    

    Returns
    -------
    rh : float
        Hipparcos RA (radians)    
    dh : float
        Hipparcos Dec (radians)

    '''
    lib.iauFk5hz.argtypes = [c_double,    c_double,   c_double,   c_double,
                             c_double_p,  c_double_p]
    lib.iauFk5hz.restype  =  None
    
    rh     = c_double()
    dh     = c_double()
    
    lib.iauFk5hz(r5,   d5,   date1,   date2,  byref(rh),  byref(dh))
    
    return  rh.value,  dh.value   

def pymFk45z(r1950,  d1950,  bepoch):
    '''
    Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
    proper motion in the FK5 system.

    Parameters
    ----------
    r1950 : float
        B1950.0 FK4 RA at epoch (rad)    
    d1950 : float
        B1950.0 FK4 Dec at epoch (rad)    
    bepoch : float
        Besselian epoch (e.g. 1979.3)

    Returns
    -------
    r2000 : float
        J2000.0 FK5 RA (rad)     
    d2000 : float
        J2000.0 FK5 Dec (rad)

    '''
    lib.iauFk45z.argtypes = [c_double,    c_double,   c_double,   
                             c_double_p,  c_double_p]
    lib.iauFk45z.restype  =  None
    
    r2000    = c_double()
    d2000    = c_double()
    
    lib.iauFk45z(r1950,  d1950,  bepoch,  byref(r2000),  byref(d2000))
    
    return  r2000.value,  d2000.value

     
def pymFk52h(r5,  d5,  dr5,  dd5,  px5,  rv5):
    '''
    Transform FK5 (J2000.0) star data into the Hipparcos system.

    Parameters
    ----------
    all FK5, equinox J2000.0, epoch J2000.0     
    r5 : float
        RA (radians)    
    d5 : float
        Dec (radians)    
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    px5 : float
        parallax (arcsec)    
    rv5 : float
        radial velocity (km/s, positive = receding)

    Returns
    -------
    all Hipparcos, epoch J2000.0    
    rh : float
        RA (radians)    
    dh : float
        Dec (radians)    
    drh : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    ddh : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    pxh : float
        parallax (arcsec)    
    rvh : float
        radial velocity (km/s, positive = receding)

    '''
    lib.iauFk52h.argtypes = [c_double,    c_double,   
                             c_double,    c_double,   c_double,    c_double,
                             c_double_p,  c_double_p,
                             c_double_p,  c_double_p, c_double_p,  c_double_p]
    lib.iauFk52h.restype  =  None
    
    rh     = c_double()
    dh     = c_double()
    drh    = c_double()
    ddh    = c_double()
    pxh    = c_double()
    rvh    = c_double()
    
    
    lib.iauFk52h(r5,   d5,   dr5,  dd5,  px5,  rv5,  
                 byref(rh),  byref(dh),
                 byref(drh),  byref(ddh), byref(pxh),  byref(rvh))
    
    return  rh.value,  dh.value,  drh.value,  ddh.value,  pxh.value,  rvh.value 


def pymFk54z(r2000,  d2000,  bepoch):
    '''
    Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
    proper motion in FK5 and parallax.

    Parameters
    ----------
    r2000 : float
        J2000.0 FK5 RA (rad)    
    d2000 : float
        J2000.0 FK5 Dec (rad)    
    bepoch : float
        Besselian epoch (e.g. 1950.0)

    Returns
    -------
    r1950 : float
        B1950.0 FK4 RA (rad) at epoch BEPOCH    
    d1950 : float
        B1950.0 FK4 Dec (rad) at epoch BEPOCH    
    dr1950 : float
        B1950.0 FK4 proper motions (rad/trop.yr)    
    dd1950 : float
        B1950.0 FK4 proper motions (rad/trop.yr)    

    '''    
    lib.iauFk54z.argtypes = [c_double,    c_double,   c_double,   
                             c_double_p,  c_double_p,
                             c_double_p,  c_double_p]
    lib.iauFk54z.restype  =  None
    
    r1950    = c_double()
    d1950    = c_double()
    dr1950   = c_double()
    dd1950   = c_double()
    
    lib.iauFk54z(r2000,  d2000,  bepoch,  byref(r1950),  byref(d1950),
                 byref(dr1950),  byref(dd1950))
    
    return  r1950.value,  d1950.value,  dr1950.value,  dd1950.value 


def pymFk425(r1950,  d1950,  dr1950,  dd1950,  p1950,  v1950):
    '''
    Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

    Parameters
    ----------
    all B1950.0, FK4    
    r1950 : float
        B1950.0 RA (rad)    
    d1950 : float
        B1950.0 Dec (rad)    
    dr1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    dd1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    p1950 : float
        parallax (arcsec)    
    v1950 : float
        radial velocity (km/s, +ve = moving away)

    Returns
    -------
    all J2000.0, FK5    
    r2000 : float
        J2000.0 RA (rad)    
    d2000 : float
        J2000.0 Dec (rad)    
    dr2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    dd2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    p2000 : float
        parallax (arcsec)    
    v2000 : float
        radial velocity (km/s, +ve = moving away)

    '''
    lib.iauFk425.argtypes = [c_double,    c_double,   c_double, 
                             c_double,    c_double,   c_double, 
                             c_double_p,  c_double_p, c_double_p,
                             c_double_p,  c_double_p, c_double_p]
    lib.iauFk425.restype  =  None
    
    r2000    = c_double()
    d2000    = c_double()
    
    dr2000   = c_double()
    dd2000   = c_double()
    
    p2000    = c_double()
    v2000    = c_double()
    
    
    lib.iauFk425(r1950,  d1950,  dr1950,  dd1950,  p1950,  v1950,
                 byref(r2000),   byref(d2000),
                 byref(dr2000),  byref(dd2000),
                 byref(p2000),   byref(v2000))
    
    return  r2000.value, d2000.value, dr2000.value, dd2000.value, p2000.value, v2000.value 

    
def pymFk524(r2000, d2000, dr2000, dd2000, p2000, v2000):
    '''
    Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

    Parameters
    ----------
    all J2000.0, FK5    
    r2000 : float
        J2000.0 RA (rad)    
    d2000 : float
        J2000.0 Dec (rad)    
    dr2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    dd2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    p2000 : float
        parallax (arcsec)    
    v2000 : float
        radial velocity (km/s, +ve = moving away)
        
    Returns
    -------
    all B1950.0, FK4    
    r1950 : float
        B1950.0 RA (rad)    
    d1950 : float
        B1950.0 Dec (rad)    
    dr1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    dd1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    p1950 : float
        parallax (arcsec)    
    v1950 : float
        radial velocity (km/s, +ve = moving away)

    '''
    lib.iauFk524.argtypes = [c_double,    c_double,   c_double, 
                             c_double,    c_double,   c_double, 
                             c_double_p,  c_double_p, c_double_p,
                             c_double_p,  c_double_p, c_double_p]
    lib.iauFk524.restype  =  None
    
    r1950    = c_double()
    d1950    = c_double()
    
    dr1950   = c_double()
    dd1950   = c_double()
    
    p1950    = c_double()
    v1950    = c_double()
    
    
    lib.iauFk524(r2000, d2000, dr2000, dd2000, p2000, v2000,
                 byref(r1950),   byref(d1950),
                 byref(dr1950),  byref(dd1950),
                 byref(p1950),   byref(v1950))
    
    return  r1950.value, d1950.value, dr1950.value, dd1950.value, p1950.value, v1950.value


def pymG2icrs(dl, db):
    '''
    Transformation from Galactic Coordinates to ICRS.

    Parameters
    ----------
    dl : float
        galactic longitude (radians)    
    db : float
        galactic latitude (radians)

    Returns
    -------
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    '''
    lib.iauG2icrs.argtypes = [c_double,    c_double,   
                              c_double_p,  c_double_p]
    lib.iauG2icrs.restype  =  None
    
    dr     = c_double()
    dd     = c_double()
    
    lib.iauG2icrs(dl, db,  byref(dr),  byref(dd))
    
    return  dr.value,  dd.value

    
def pymGc2gd(n, xyz):
    '''
    Transform geocentric coordinates to geodetic using the specified
    reference ellipsoid.

    Parameters
    ----------
    n : int
        ellipsoid identifier    
    xyz : numpy.matrix(1,3)
        geocentric vector
        
    n    ellipsoid    
    1     WGS84    
    2     GRS80    
    3     WGS72    
    
    Returns
    -------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    '''
    lib.iauGc2gd.argtypes = [c_int,   vector_double,  
                             c_double_p,  c_double_p, c_double_p]
    lib.iauGc2gd.restype  = c_int
    
    elong     = c_double()
    phi       = c_double()
    height    = c_double()
    
    lib.iauGc2gd(n, xyz,  byref(elong),  byref(phi), byref(height))
    
    return  elong.value,  phi.value,  height.value


def pymGc2gde(a, f, xyz):
    '''
    Transform geocentric coordinates to geodetic for a reference
    ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius    
    f : float
        flattening    
    xyz : numpy.matrix(1,3)
        geocentric vector

    Returns
    -------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    '''
    lib.iauGc2gde.argtypes = [c_double,    c_double,  vector_double,  
                              c_double_p,  c_double_p, c_double_p]
    lib.iauGc2gde.restype  = c_int
    
    elong     = c_double()
    phi       = c_double()
    height    = c_double()
    
    lib.iauGc2gde(a, f, xyz,  byref(elong),  byref(phi),  byref(height))
    
    return  elong.value,  phi.value,  height.value 

def pymGd2gc(n, elong,  phi,  height):
    '''
    Transform geodetic coordinates to geocentric using the specified
    reference ellipsoid.

    Parameters
    ----------
    n : int
        ellipsoid identifier
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    n    ellipsoid    
    1     WGS84    
    2     GRS80    
    3     WGS72    
    
    Returns
    -------
    xyz : numpy.matrix(1,3)
        geocentric vector

    '''
    lib.iauGd2gc.argtypes = [c_int,   c_double,  c_double, c_double, 
                             vector_double]
    lib.iauGd2gc.restype  = c_int
    
    xyz  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C')) 
    
    lib.iauGd2gc(n,  elong,  phi,  height,  xyz)
    
    return  xyz

def pymGd2gce(a,  f,  elong,  phi,  height):
    '''
    Transform geodetic coordinates to geocentric for a reference
    ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius    
    f : float
        flattening     
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    Returns
    -------
    xyz : numpy.matrix(1,3)
        geocentric vector

    '''
    lib.iauGd2gce.argtypes = [c_double,   c_double,  c_double,  
                              c_double, c_double,   vector_double]
    lib.iauGd2gce.restype  = c_int
    
    xyz  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C')) 
    
    lib.iauGd2gce(a,  f,   elong,  phi,  height,  xyz)
    
    return  xyz   

def pymH2fk5(rh,  dh,  drh,  ddh,  pxh,  rvh):
    '''
    Transform Hipparcos star data into the FK5 (J2000.0) system.

    Parameters
    ----------
    all Hipparcos, epoch J2000.0    
    rh : float
        RA (radians)    
    dh : float
        Dec (radians)    
    drh : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    ddh : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    pxh : float
        parallax (arcsec)    
    rvh : float
        radial velocity (km/s, positive = receding)
        
    Returns
    -------
    all FK5, equinox J2000.0, epoch J2000.0     
    r5 : float
        RA (radians)    
    d5 : float
        Dec (radians)    
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    px5 : float
        parallax (arcsec)    
    rv5 : float
        radial velocity (km/s, positive = receding)

    '''
    lib.iauH2fk5.argtypes = [c_double,    c_double,   
                             c_double,    c_double,   c_double,    c_double,
                             c_double_p,  c_double_p,
                             c_double_p,  c_double_p, c_double_p,  c_double_p]
    lib.iauH2fk5.restype  =  None
    
    r5     = c_double()
    d5     = c_double()
    dr5    = c_double()
    dd5    = c_double()
    px5    = c_double()
    rv5    = c_double()
    
    
    lib.iauH2fk5(rh,  dh,  drh,  ddh,  pxh,  rvh,  
                 byref(r5),  byref(d5),
                 byref(dr5),  byref(dd5), byref(px5), byref(rv5))
    
    return  r5.value,  d5.value,  dr5.value,  dd5.value,  px5.value,  rv5.value


def pymHfk5z(rh,   dh,   date1,   date2):
    '''
    Transform a Hipparcos star position into FK5 J2000.0, assuming
    =zero Hipparcos proper motion.

    Parameters
    ----------
    rh : float
        Hipparcos RA (radians)    
    dh : float
        Hipparcos Dec (radians)    
    date1 : float
        TDB date    
    date2 : float
        TDB date

    Returns
    -------
    all FK5, equinox J2000.0, date date1+date2      
    r5 : float
        RA (radians)    
    d5 : float
        Dec (radians)    
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear)    

    '''
    lib.iauHfk5z.argtypes = [c_double,    c_double,   c_double,   c_double,
                             c_double_p,  c_double_p, c_double_p,  c_double_p]
    lib.iauHfk5z.restype  =  None
    
    r5     = c_double()
    d5     = c_double()
    dr5    = c_double()
    dd5    = c_double()
    
    
    lib.iauHfk5z(rh,   dh,   date1,   date2,  
                 byref(r5),  byref(d5), byref(dr5),  byref(dd5))
    
    return  r5.value,  d5.value, dr5.value,  dd5.value   


def pymIcrs2g(dr, dd):
    '''
    Transformation from ICRS to Galactic Coordinates.

    Parameters
    ----------
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    Returns
    -------
    dl : float
        galactic longitude (radians)    
    db : float
        galactic latitude (radians)

    '''
    lib.iauIcrs2g.argtypes = [c_double,    c_double,   
                              c_double_p,  c_double_p]
    lib.iauIcrs2g.restype  =  None
    
    dl    = c_double()
    db    = c_double()
    
    lib.iauIcrs2g(dr,  dd,  byref(dl),  byref(db))
    
    return  dl.value,  db.value   


def pymLteceq(epj, dl, db):
    '''
    Transformation from ecliptic coordinates (mean equinox and ecliptic
    of date) to ICRS RA,Dec, using a long-term precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)
    dl : float
        ecliptic longitude (radians)    
    db : float
        ecliptic latitude (radians)

    Returns
    -------
    dr : float
        ICRS right ascension (radians)     
    dd : float
        ICRS right declination (radians)

    '''
    lib.iauLteceq.argtypes = [c_double,    c_double, c_double,   
                              c_double_p,  c_double_p]
    lib.iauLteceq.restype  =  None
    
    dr     = c_double()
    dd     = c_double()
    
    lib.iauLteceq(epj, dl, db,  byref(dr),  byref(dd))
    
    return  dr.value,  dd.value 


def pymLtecm(epj):
    '''
    ICRS equatorial to ecliptic rotation matrix, long-term.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)

    Returns
    -------
    rm : numpy.matrix(3,3)
        ICRS to ecliptic rotation matrix

    '''
    lib.iauLtecm.argtypes = [c_double, vector_double3]
    lib.iauLtecm.restype  =  None
    
    rm = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauLtecm(epj, rm)
    
    return rm     


def pymLteqec(epj, dr, dd):
    '''
    Transformation from ICRS equatorial coordinates to ecliptic
    coordinates (mean equinox and ecliptic of date) using a long-term
    precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    Returns
    -------
    dl : float
        ecliptic longitude (radians)    
    db : float
        ecliptic latitude (radians)

    '''
    lib.iauLteqec.argtypes = [c_double,    c_double,  c_double,  
                              c_double_p,  c_double_p]
    lib.iauLteqec.restype  =  None
    
    dl    = c_double()
    db    = c_double()
    
    lib.iauLteqec(epj, dr,  dd,  byref(dl),  byref(db))
    
    return  dl.value,  db.value   

pym_pvstar_msg = {
                  0: 'OK',
                 -1: 'superluminal speed',
                 -2: 'null position vector' 
                  }  

def pymPvstar(pv):
    '''
    Convert star position+velocity vector to catalog coordinates.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector (au, au/day)

    Raises
    ------
    ValueError
        0: 'OK',
       -1: 'superluminal speed',
       -2: 'null position vector' 

    Returns
    -------
    ra : float
        right ascension (radians)    
    dec : float
        declination (radians)    
    pmr : float
        RA proper motion (radians/year)     
    pmd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, positive = receding)

    '''
    lib.iauPvstar.argtypes = [vector_double2,  c_double_p,  c_double_p,
                              c_double_p, c_double_p, c_double_p, c_double_p]
    
    lib.iauPvstar.restype  = c_int
 
    ra   = c_double()
    dec  = c_double()
    pmr  = c_double()
    pmd  = c_double()
    
    px   = c_double()
    rv   = c_double()
   
    j = lib.iauPvstar(pv,
                      byref(ra),  byref(dec),  byref(pmr),  byref(pmd),
                      byref(px),  byref(rv))

    if   j < 0:
        raise ValueError(pym_pvstar_msg[j])
#    elif j > 0:
#         ws.warn(iau_pvstar_msg[j], UserWarning, 2)
        
    return ra.value, dec.value, pmr.value, pmd.value, px.value, rv.value
 
 
def pymStarpm(ra1,  dec1, pmr1, pmd1,
              px1,  rv1,
              ep1a, ep1b, ep2a, ep2b):
    '''
    Star proper motion:  update star catalog data for space motion.

    Parameters
    ----------
    ra1 : float
        right ascension (radians), before    
    dec1 : float
        declination (radians), before    
    pmr1 : float
        RA proper motion (radians/year), before    
    pmd1 : float
        Dec proper motion (radians/year), before    
    px1 : float
        Dec proper motion (radians/year), before    
    rv1 : float
        radial velocity (km/s, +ve = receding), before    
    ep1a : float
        "before" epoch, part A    
    ep1b : float
        "before" epoch, part B    
    ep2a : float
        "after" epoch, part A    
    ep2b : float
        "after" epoch, part B

    Raises
    ------
    ValueError
        4: 'solution not converge',
        2: 'excessive velocity',
        1: 'distance overridden',
        0: 'no warnings or errors',
       -1: 'system error',

    Returns
    -------
    ra2 : float
        right ascension (radians), after    
    dec2 : float
        declination (radians), after    
    pmr2 : float
        RA proper motion (radians/year), after    
    pmd2 : float
        Dec proper motion (radians/year), after    
    px2 : float
        parallax (arcseconds), after    
    rv2 : flaot
        radial velocity (km/s, +ve = receding), after    

    '''
    lib.iauStarpm.argtypes = [c_double, c_double, c_double, c_double,
                              c_double, c_double,   
                              c_double, c_double, c_double, c_double, 
                              c_double_p, c_double_p, c_double_p, c_double_p,
                              c_double_p, c_double_p]
    
    lib.iauStarpm.restype  = c_int
 
    ra2   = c_double()
    dec2  = c_double()
    pmr2  = c_double()
    pmd2  = c_double()
    
    px2   = c_double()
    rv2   = c_double()
   
    j = lib.iauStarpm(ra1,  dec1, pmr1, pmd1,
                      px1,  rv1,
                      ep1a, ep1b,  ep2a, ep2b, 
                      byref(ra2),  byref(dec2),  byref(pmr2),  byref(pmd2),
                      byref(px2),  byref(rv2))
  

    if   j < 0:
         raise ValueError(pym_pmsafe_msg[j])
    elif j > 0:
         ws.warn(pym_pmsafe_msg[j], UserWarning, 2)
        
    return ra2.value, dec2.value, pmr2.value, pmd2.value, px2.value, rv2.value


pym_starpv_msg = {
                  0: 'no warnings',
                  1: 'distance overridden',
                  2: 'excessive speed',
                  4: 'solution did not converge '
                  }  

def pymStarpv(ra, dec, pmr, pmd, px, rv):
    '''
    Convert star catalog coordinates to position+velocity vector.

    Parameters
    ----------
    ra : float
        right ascension (radians)    
    dec : float
        declination (radians)    
    pmr : float
        RA proper motion (radians/year)     
    pmd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, positive = receding)

    Raises
    ------
    ValueError
        0: 'no warnings',
        1: 'distance overridden',
        2: 'excessive speed',
        4: 'solution did not converge '

    Returns
    -------
    pv : numpy.matrix(2,3)
        pv-vector (au, au/day)

    '''
    lib.iauStarpv.argtypes = [c_double, c_double,  
                              c_double, c_double, c_double, c_double,
                              vector_double2]
    
    lib.iauStarpv.restype  = c_int
 
    pv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    j = lib.iauStarpv(ra,  dec,  pmr,  pmd,   px,  rv, 
                      pv)

    if   j < 0:
         raise ValueError(pym_starpv_msg[j])
    elif j > 0:
         ws.warn(pym_starpv_msg[j], UserWarning, 2)
        
    return pv
    
