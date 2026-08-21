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

def pymBi00():
    '''
    Frame bias components of IAU 2000 precession-nutation models;  part
    of the Mathews-Herring-Buffett (MHB2000) nutation series, with
    additions.

    Returns
    -------
    dpsibi : float
        longitude and obliquity corrections    
    depsbi : float
        longitude and obliquity corrections    
    dra : flaot
        the ICRS RA of the J2000.0 mean equinox

    '''
    lib.iauBi00.argtypes = [c_double_p, c_double_p,  c_double_p]
    lib.iauBi00.restype  = None
    dpsibi   = c_double()
    depsbi   = c_double()
    dra      = c_double()
    
    lib.iauBi00(byref(dpsibi),  byref(depsbi),  byref(dra))
    
    return  dpsibi.value, depsbi.value, dra.value


#vector_double3  = nt.ndpointer(shape=(3,3), dtype=np.double, flags='C')
#using vector_double3 for rbpn[3][3]
    
def pymBpn2xy(rbpn):
    '''
    Extract from the bias-precession-nutation matrix the X,Y coordinates
    of the Celestial Intermediate Pole.

    Parameters
    ----------
    rbpn : numpy.matrix(3,3)
        celestial-to-true matrix

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole

    '''
    lib.iauBpn2xy.argtypes = [vector_double3, c_double_p,  c_double_p]
    lib.iauBpn2xy.restype  = None
    
    x   = c_double()
    y   = c_double()
    
    
    lib.iauBpn2xy(rbpn,  byref(x),  byref(y))
    
    return  x.value, y.value  

#version 2
#c_2d_double = nt.ndpointer(dtype=np.double, ndim=2, flags="C")
#using 2-dim array for rbpn[3][3]    
def pymBpn2xy_A(rbpn):
    
    lib.iauBpn2xy.argtypes = [c_2d_double, c_double_p,  c_double_p]
    lib.iauBpn2xy.restype  = None
    x   = c_double()
    y   = c_double()
    
    
    lib.iauBpn2xy(rbpn,  byref(x),  byref(y))
    
    return  x.value, y.value    

def pymC2i00a(date1, date2):
    '''
    Form the celestial-to-intermediate matrix for a given date using the
    IAU 2000A precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date   

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    lib.iauC2i00a.argtypes = [c_double,  c_double, vector_double3]
    lib.iauC2i00a.restype  = None
    
    rc2i = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauC2i00a(date1, date2, rc2i)
    
    return  rc2i 

#2023-05-06   

def pymC2i00b(date1, date2):
    '''
    Form the celestial-to-intermediate matrix for a given date using the
    IAU 2000B precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date   

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    lib.iauC2i00b.argtypes = [c_double,  c_double, vector_double3]
    lib.iauC2i00b.restype  = None
    
    rc2i = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauC2i00b(date1, date2, rc2i)
    
    return  rc2i 

def pymC2ibpn(date1, date2, rbpn):
    '''
    Form the celestial-to-intermediate matrix for a given date given
    the bias-precession-nutation matrix.  IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date    
    rbpn : numpy.matrix(3,3)
        celestial-to-true matrix

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    lib.iauC2ibpn.argtypes = [c_double,  c_double, vector_double3,
                              vector_double3]
    lib.iauC2ibpn.restype  = None
    
    rc2i = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauC2ibpn(date1, date2, rbpn, rc2i)
    
    return  rc2i 

def pymC2t00a(tta,  ttb,  uta,  utb,  xp,  yp):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1 and
    the polar motion, using the IAU 2000A precession-nutation model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole  

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    lib.iauC2t00a.argtypes = [c_double,  c_double, c_double,  c_double,
                              c_double,  c_double,
                              vector_double3]
    lib.iauC2t00a.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2t00a(tta,  ttb,  uta,  utb,  xp,  yp, rc2t)
    
    return  rc2t

def pymC2t00b(tta,  ttb,  uta,  utb,  xp,  yp):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1 and
    the polar motion, using the IAU 2000B precession-nutation model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole  

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix
        
    '''
    lib.iauC2t00b.argtypes = [c_double,  c_double, c_double,  c_double,
                              c_double,  c_double,
                              vector_double3]
    lib.iauC2t00b.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2t00b(tta,  ttb,  uta,  utb,  xp,  yp, rc2t)
    
    return  rc2t

def pymC2teqx(rbpn,  gst,  rpom):
    '''
    Assemble the celestial to terrestrial matrix from equinox-based
    components (the celestial-to-true matrix, the Greenwich Apparent
    Sidereal Time and the polar motion matrix).

    Parameters
    ----------
    rbpn : numpy.matrix(3,3)
        celestial-to-true matrix
    gst : float
        Greenwich (apparent) Sidereal Time (radians)
    rpom : numpy.matrix(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    lib.iauC2teqx.argtypes = [vector_double3,  c_double, vector_double3,
                              vector_double3]
    lib.iauC2teqx.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2teqx(rbpn,  gst,  rpom,  rc2t)
    
    return  rc2t

def pymC2tpe(tta,  ttb,  uta,  utb,  dpsi,  deps,  xp,  yp):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1,
    the nutation and the polar motion.  IAU 2000.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    dpsi : float
        nutation     
    deps : float
        nutation         
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole  

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    lib.iauC2tpe.argtypes = [c_double,  c_double, c_double,  c_double,
                             c_double,  c_double, c_double,  c_double,
                             vector_double3]
    lib.iauC2tpe.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2tpe(tta,  ttb,  uta,  utb, dpsi,  deps,  xp,  yp, rc2t)
    
    return  rc2t

#2023-05-09   
def pymC2txy(tta,  ttb,  uta,  utb,  x,  y,  xp,  yp):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1,
    the CIP coordinates and the polar motion.  IAU 2000.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole    

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    lib.iauC2txy.argtypes = [c_double,  c_double, c_double,  c_double,
                             c_double,  c_double, c_double,  c_double,
                             vector_double3]
    lib.iauC2txy.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2txy(tta,  ttb,  uta,  utb, x,  y,  xp,  yp, rc2t)
    
    return  rc2t


#def pymCal2jd(iyear, imon, iday):  has been done in time tools
    
#pym_CR: to be done in vector tools

#pym_DAT: has been done in time tools

def pymEe00(date1,  date2,   epsa,   dpsi):
    '''
    The equation of the equinoxes, compatible with IAU 2000 resolutions,
    given the nutation in longitude and the mean obliquity.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    epsa : float
        mean obliquity    
    dpsi : float
        nutation in longitude

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    lib.iauEe00.argtypes = [c_double,  c_double, c_double,  c_double]
    lib.iauEe00.restype  =  c_double
    
    return   lib.iauEe00(date1,  date2,   epsa,  dpsi)

def pymEe00a(date1,  date2):
    '''
    Equation of the equinoxes, compatible with IAU 2000 resolutions.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    lib.iauEe00a.argtypes = [c_double,  c_double]
    lib.iauEe00a.restype  =  c_double
    
    return   lib.iauEe00a(date1,  date2)

def pymEe00b(date1,  date2):
    '''
    Equation of the equinoxes, compatible with IAU 2000 resolutions but
    using the truncated nutation model IAU 2000B.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    lib.iauEe00b.argtypes = [c_double,  c_double]
    lib.iauEe00b.restype  =  c_double
    
    return   lib.iauEe00b(date1,  date2)

def pymEe06a(date1,  date2):
    '''
    Equation of the equinoxes, compatible with IAU 2000 resolutions and
    IAU 2006/2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    lib.iauEe06a.argtypes = [c_double,  c_double]
    lib.iauEe06a.restype  =  c_double
    
    return   lib.iauEe06a(date1,  date2)

def pymEect00(date1,  date2):
    '''
    Equation of the equinoxes complementary terms, consistent with
    IAU 2000 resolutions.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        complementary terms

    '''
    lib.iauEect00.argtypes = [c_double,  c_double]
    lib.iauEect00.restype  =  c_double
    
    return   lib.iauEect00(date1,  date2)

def pymEo06a(date1,  date2):
    '''
    Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the equation of the origins in radians

    '''
    lib.iauEo06a.argtypes = [c_double,  c_double]
    lib.iauEo06a.restype  =  c_double
    
    return   lib.iauEo06a(date1,  date2)

def pymEors(rnpb,  s):
    '''
    Equation of the origins, given the classical NPB matrix and the
    quantity s.

    Parameters
    ----------
    rnpb : numpy.matrix(3,3) 
        classical nutation x precession x bias matrix    
    s : float
        the quantity s (the CIO locator) in radians

    Returns
    -------
    function value : float
        the equation of the origins in radians

    '''
    lib.iauEors.argtypes = [vector_double3,  c_double]
    lib.iauEors.restype  =  c_double
    
    return   lib.iauEors(rnpb,  s)

def pymEqeq94(date1,  date2):
    '''
    Equation of the equinoxes, IAU 1994 model.

    Parameters
    ----------
    date1 : float
        TDB date    
    date2 : float
        TDB date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    lib.iauEqeq94.argtypes = [c_double,  c_double]
    lib.iauEqeq94.restype  =  c_double
    
    return   lib.iauEqeq94(date1,  date2)

def pymEra00(dj1,  dj2):
    '''
    Earth rotation angle (IAU 2000 model).

    Parameters
    ----------
    dj1 : float
        UT1 as a 2-part Julian Date    
    dj2 : float
        UT1 as a 2-part Julian Date

    Returns
    -------
    function value : float
        Earth rotation angle (radians), range 0-2pi

    '''
    lib.iauEra00.argtypes = [c_double,  c_double]
    lib.iauEra00.restype  =  c_double
    
    return   lib.iauEra00(dj1,  dj2)

def pymFad03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean elongation of the Moon from the Sun.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        D, radians

    '''
    lib.iauFad03.argtypes = [c_double]
    lib.iauFad03.restype  =  c_double
    
    return   lib.iauFad03(t)

#2023-05-10   
def pymFae03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Earth.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Earth, radians

    '''
    lib.iauFae03.argtypes = [c_double]
    lib.iauFae03.restype  =  c_double
    
    return   lib.iauFae03(t)

def pymFaf03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of the Moon minus mean longitude of the ascending
    node.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        F, radians

    '''
    lib.iauFaf03.argtypes = [c_double]
    lib.iauFaf03.restype  =  c_double
    
    return   lib.iauFaf03(t)

def pymFaju03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Jupiter.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Jupiter, radians

    '''
    lib.iauFaju03.argtypes = [c_double]
    lib.iauFaju03.restype  =  c_double
    
    return   lib.iauFaju03(t)


def pymFal03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean anomaly of the Moon.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        l, radians

    '''
    lib.iauFal03.argtypes = [c_double]
    lib.iauFal03.restype  =  c_double
    
    return   lib.iauFal03(t)

def pymFalp03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean anomaly of the Sun.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        l', radians

    '''
    lib.iauFalp03.argtypes = [c_double]
    lib.iauFalp03.restype  =  c_double
    
    return   lib.iauFalp03(t)

def pymFama03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Mars.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Mars, radians

    '''
    lib.iauFama03.argtypes = [c_double]
    lib.iauFama03.restype  =  c_double
    
    return   lib.iauFama03(t)

def pymFame03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Mercury.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Mercury, radians

    '''
    lib.iauFame03.argtypes = [c_double]
    lib.iauFame03.restype  =  c_double
    
    return   lib.iauFame03(t)

def pymFane03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Neptune.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Neptune, radians

    '''
    lib.iauFane03.argtypes = [c_double]
    lib.iauFane03.restype  =  c_double
    
    return   lib.iauFane03(t)

def pymFaom03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of the Moon's ascending node.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        Omega, radians 

    '''
    lib.iauFaom03.argtypes = [c_double]
    lib.iauFaom03.restype  =  c_double
    
    return   lib.iauFaom03(t)

def pymFapa03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    general accumulated precession in longitude.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        general precession in longitude, radians

    '''
    lib.iauFapa03.argtypes = [c_double]
    lib.iauFapa03.restype  =  c_double
    
    return   lib.iauFapa03(t)

def pymFasa03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Saturn.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Saturn, radians

    '''
    lib.iauFasa03.argtypes = [c_double]
    lib.iauFasa03.restype  =  c_double
    
    return   lib.iauFasa03(t)

def pymFaur03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Uranus.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Uranus, radians

    '''
    lib.iauFaur03.argtypes = [c_double]
    lib.iauFaur03.restype  =  c_double
    
    return   lib.iauFaur03(t)

def pymFave03(t):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Venus.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Venus, radians

    '''
    lib.iauFave03.argtypes = [c_double]
    lib.iauFave03.restype  =  c_double
    
    return   lib.iauFave03(t)

#2023-05-11
def pymFw2m(gamb,  phib,  psi,  eps):
    '''
    Form rotation matrix given the Fukushima-Williams angles.

    Parameters
    ----------
    gamb : float
        F-W angle gamma_bar (radians)    
    phib : float
        F-W angle phi_bar (radians)    
    psi : float
        F-W angle psi (radians)    
    eps : float
        F-W angle epsilon (radians)

    Returns
    -------
    r : numpy.matrix(3,3)
        rotation matrix

    '''
    lib.iauFw2m.argtypes = [c_double,  c_double, c_double, c_double,
                            vector_double3]
    lib.iauFw2m.restype  = None
    
    r  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauFw2m(gamb,  phib,  psi,  eps,  r)
    
    return  r 

#2023-05-20
def pymFw2xy(gamb,   phib,  psi,  eps):
    '''
    CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

    Parameters
    ----------
    gamb : float
        F-W angle gamma_bar (radians)    
    phib : float
        F-W angle phi_bar (radians)    
    psi : float
        F-W angle psi (radians)    
    eps : float
        F-W angle epsilon (radians)

    Returns
    -------
    x : float
        CIP unit vector X,Y     
    y : float
        CIP unit vector X,Y

    '''
    lib.iauFw2xy.argtypes =[c_double, c_double, c_double, c_double,
                            c_double_p, c_double_p]
    
    lib.iauFw2xy.restype  = None
    
    x    = c_double()
    y    = c_double()
    
    lib.iauFw2xy(gamb,   phib,  psi,  eps, 
                     byref(x),   byref(y))
    
    return  x.value, y.value  



#iauIR, vector matrix 
def pymNum00a(date1, date2):
    '''
    Form the matrix of nutation for a given date, IAU 2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    lib.iauNum00a.argtypes = [c_double,  c_double, vector_double3]
    lib.iauNum00a.restype  = None
    
    rmatn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauNum00a(date1, date2, rmatn)
    
    return  rmatn

def pymNum00b(date1, date2):
    '''
    Form the matrix of nutation for a given date, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    lib.iauNum00b.argtypes = [c_double,  c_double, vector_double3]
    lib.iauNum00b.restype  = None
    
    rmatn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauNum00b(date1, date2, rmatn)
    
    return  rmatn

def pymNum06a(date1, date2):
    '''
    Form the matrix of nutation for a given date, IAU 2006/2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    lib.iauNum06a.argtypes = [c_double,  c_double, vector_double3]
    lib.iauNum06a.restype  = None
    
    rmatn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauNum06a(date1, date2, rmatn)
    
    return  rmatn

def pymNumat(epsa,   dpsi,   deps):
    '''
    Form the matrix of nutation.

    Parameters
    ----------
    epsa : float 
        mean obliquity of date    
    dpsi : float
        nutation    
    deps : float
        nutation

    Returns
    -------
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    lib.iauNumat.argtypes = [c_double,  c_double,  c_double,  vector_double3]
    lib.iauNumat.restype  = None
    
    rmatn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauNumat(epsa,   dpsi,   deps,  rmatn)
    
    return  rmatn

def pymNut00a(date1,  date2):
    '''
    Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
    with free core nutation omitted).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary    
    deps : float
        nutation, luni-solar + planetary

    '''
    lib.iauNut00a.argtypes = [c_double,  c_double,  c_double_p,  c_double_p]
    lib.iauNut00a.restype  = None
    
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut00a(date1,  date2,  byref(dpsi),   byref(deps))
    
    return  dpsi.value, deps.value 

def pymNut00b(date1,  date2):
    '''
    Nutation, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary    
    deps : float
        nutation, luni-solar + planetary


    '''
    lib.iauNut00b.argtypes = [c_double,  c_double,  c_double_p,  c_double_p]
    lib.iauNut00b.restype  = None
    
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut00b(date1,  date2,  byref(dpsi),   byref(deps))
    
    return  dpsi.value, deps.value 

def pymNut06a(date1,  date2):
    '''
    IAU 2000A nutation with adjustments to match the IAU 2006
    precession.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary    
    deps : float
        nutation, luni-solar + planetary

    '''
    lib.iauNut06a.argtypes = [c_double,  c_double,  c_double_p,  c_double_p]
    lib.iauNut06a.restype  = None
    
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut06a(date1,  date2,  byref(dpsi),   byref(deps))
    
    return  dpsi.value, deps.value 

def pymNut80(date1,  date2):
    '''
    Nutation, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation in longitude (radians)    
    deps : float
        nutation in obliquity (radians)

    '''
    lib.iauNut80.argtypes = [c_double,  c_double,  c_double_p,  c_double_p]
    lib.iauNut80.restype  = None
    
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut80(date1,  date2,  byref(dpsi),   byref(deps))
    
    return  dpsi.value, deps.value 

def pymNutm80(date1, date2):
    '''
    Form the matrix of nutation for a given date, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TDB date      
    date2 : float
        TDB date

    Returns
    -------
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    lib.iauNutm80.argtypes = [c_double,  c_double, vector_double3]
    lib.iauNutm80.restype  = None
    
    rmatn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauNutm80(date1, date2, rmatn)
    
    return  rmatn

#2023-05-21
def pymObl06(date1, date2):
    '''
    Mean obliquity of the ecliptic, IAU 2006 precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        obliquity of the ecliptic (radians)

    '''
    lib.iauObl06.argtypes =[c_double, c_double]
    
    lib.iauObl06.restype  = c_double
    
    return  lib.iauObl06(date1, date2)

def pymObl80(date1, date2):
    '''
    Mean obliquity of the ecliptic, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        obliquity of the ecliptic (radians)

    '''
    lib.iauObl80.argtypes =[c_double, c_double]
    
    lib.iauObl80.restype  = c_double
    
    return  lib.iauObl80(date1, date2)

def pymP06e(date1, date2):
    '''
    Precession angles, IAU 2006, equinox based.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    eps0 : float
        epsilon_0    
    psia : float
        psi_A    
    oma : float
        omega_A    
    bpa : float
        P_A    
    bqa : float
        Q_A    
    pia : float
        pi_A    
    bpia : float
        Pi_A    
    epsa : float
        obliquity epsilon_A    
    chia : float
        chi_A     
    za : float
        z_A     
    zetaa : float
        zeta_A    
    thetaa : float 
        theta_A    
    pa : float
        p_A    
    gma : float
        F-W angle gamma_J2000    
    phi : float
        F-W angle phi_J2000    
    psi : float
        F-W angle psi_J2000    

    '''
    lib.iauP06e.argtypes =[c_double,   c_double,
                           c_double_p, c_double_p, c_double_p, c_double_p,
                           c_double_p, c_double_p, c_double_p, 
                           c_double_p, c_double_p, c_double_p, c_double_p, 
                           c_double_p, c_double_p,
                           c_double_p, c_double_p, c_double_p  
                           ]
    
    lib.iauP06e.restype  = None
    
    eps0  = c_double()
    psia  = c_double()
    oma   = c_double()
    bpa   = c_double()
  
    bqa   = c_double() 
    pia   = c_double() 
    bpia  = c_double()  
       
    epsa  = c_double()
    chia  = c_double()
    za    = c_double()
    zetaa = c_double()
    
    thetaa = c_double() 
    pa     = c_double()
    
    gam   = c_double() 
    phi   = c_double() 
    psi   = c_double()  
    
    lib.iauP06e(date1, date2, 
                byref(eps0), byref(psia), byref(oma), byref(bpa),
                byref(bqa),  byref(pia),  byref(bpia),
                byref(epsa), byref(chia), byref(za), byref(zetaa),
                byref(thetaa), byref(pa),
                byref(gam),  byref(phi),  byref(psi))
    
    return   eps0.value,  psia.value,   oma.value,  bpa.value,  \
             bqa.value,   pia.value,  bpia.value,               \
             epsa.value,  chia.value,    za.value, zetaa.value, \
             thetaa.value,   pa.value,                          \
             gam.value,   phi.value,  psi.value
                        
def pymPb06(date1,  date2):
    '''
    This function forms three Euler angles which implement general
    precession from epoch J2000.0, using the IAU 2006 model.  Frame
    bias (the offset between ICRS and mean J2000.0) is included.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    bzeta : float
        1st rotation: radians cw around z    
    ca : float
        3rd rotation: radians cw around z     
    btheta : float
        2nd rotation: radians ccw around y

    '''
    lib.iauPb06.argtypes = [c_double,   c_double,  
                            c_double_p, c_double_p,  c_double_p]
    lib.iauPb06.restype  = None
    
    bzeta  = c_double()
    bz     = c_double()
    btheta = c_double()
    
    lib.iauPb06(date1,  date2,  byref(bzeta),  byref(bz),  byref(btheta))
    
    return  bzeta.value,  bz.value,  btheta.value

#2023-05-26 
def pymPn06a(date1,  date2):
    '''
    Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based use
    indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation    
    deps : flaot
        nutation     
    epsa : flaot
        mean obliquity    
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix     
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    lib.iauPn06a.argtypes = [c_double,    c_double,  
                             c_double_p,  c_double_p,  c_double_p,  
                             vector_double3, vector_double3, vector_double3,
                             vector_double3, vector_double3 
                            ]
    lib.iauPn06a.restype  = None
    
    dpsi    = c_double()
    deps    = c_double()
    epsa    = c_double()
    
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rn   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPn06a(date1,  date2,  byref(dpsi),  byref(deps),  byref(epsa),  
                 rb,  rp,  rbp,  rn,  rbpn)
    
    return  dpsi.value,   deps.value,   epsa.value,  rb,  rp,  rbp,  rn,  rbpn     
           
def pymPnm00a(date1, date2):
    '''
    Form the matrix of precession-nutation for a given date (including
    frame bias), equinox based, IAU 2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbpn : numpy.matrix(3,3)
        bias-precession-nutation matrix

    '''
    lib.iauPnm00a.argtypes = [c_double,  c_double, vector_double3]
    lib.iauPnm00a.restype  = None
    
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPnm00a(date1, date2, rbpn)
    
    return  rbpn

def pymPnm00b(date1, date2):
    '''
    Form the matrix of precession-nutation for a given date (including
    frame bias), equinox-based, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbpn : numpy.matrix(3,3)
        bias-precession-nutation matrix

    '''
    lib.iauPnm00b.argtypes = [c_double,  c_double, vector_double3]
    lib.iauPnm00b.restype  = None
    
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPnm00b(date1, date2, rbpn)
    
    return  rbpn

def pymPnm06a(date1, date2):
    '''
    Form the matrix of precession-nutation for a given date (including
    frame bias), equinox based, IAU 2006 precession and IAU 2000A
    nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbpn : numpy.matrix(3,3)
        bias-precession-nutation matrix

    '''
    lib. iauPnm06a.argtypes = [c_double,  c_double, vector_double3]
    lib. iauPnm06a.restype  = None
    
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib. iauPnm06a(date1, date2, rbpn)
    
    return  rbpn

def pymPnm80(date1, date2):
    '''
    Form the matrix of precession/nutation for a given date, IAU 1976
    precession model, IAU 1980 nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : flaot
        TT as a 2-part Julian Date

    Returns
    -------
    rmatpn : numpy.matrix(3,3)
        combined precession/nutation matrix

    '''
    lib. iauPnm80.argtypes = [c_double,  c_double, vector_double3]
    lib. iauPnm80.restype  = None
    
    rmatpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib. iauPnm80(date1, date2, rmatpn)
    
    return  rmatpn

def pymPom00(xp,   yp,   sp):
    '''
    Form the matrix of polar motion for a given date, IAU 2000.

    Parameters
    ----------
    xp : float
        coordinates of the pole (radians)    
    yp : float
        coordinates of the pole (radians)    
    sp : float
        the TIO locator s' (radians)

    Returns
    -------
    rpom : numpy.matrix(3,3)
        polar-motion matrix

    '''
    lib. iauPom00.argtypes = [c_double,  c_double, c_double, vector_double3]
    lib. iauPom00.restype  = None
    
    rpom = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib. iauPom00(xp,   yp,   sp,  rpom)
    
    return  rpom

def pymPr00(date1, date2):
    '''
    Precession-rate part of the IAU 2000 precession-nutation models
    (part of MHB2000).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsipr : float 
        precession corrections    
    depspr : float
        precession corrections

    '''
    lib.iauPr00.argtypes = [c_double,  c_double,  c_double_p,  c_double_p]
    lib.iauPr00.restype  = None
    
    dpsipr   = c_double()
    depspr   = c_double()

#2023-08-06 update
    #lib.iauPr00(date1,  date2,  dpsipr,  depspr)
    lib.iauPr00(date1,  date2,  byref(dpsipr),  byref(depspr))
    
    return   dpsipr.value, depspr.value

def pymPrec76(date01, date02, date11, date12):
    '''
    IAU 1976 precession model.
    
    This function forms the three Euler angles which implement general
    precession between two dates, using the IAU 1976 model (as for the
    FK5 catalog).

    Parameters
    ----------
    date01 : float
        TDB starting date    
    date02 : float
        TDB starting date    
    date11 : float
        TDB ending date    
    date12 : flaot
        TDB ending date    

    Returns
    -------
    zeta : float
        1st rotation: radians cw around z    
    z : float
        3rd rotation: radians cw around z    
    theta : float
        2nd rotation: radians ccw around y

    '''
    lib.iauPrec76.argtypes = [c_double,  c_double, c_double,  c_double,  
                              c_double_p,  c_double_p, c_double_p]
    lib.iauPrec76.restype  = None
    
    zeta   = c_double()
    z      = c_double()
    theta  = c_double()
    
    
    lib.iauPrec76(date01,  date02,  date11,  date12, 
                  byref(zeta), byref(z), byref(theta))
    
    return   zeta.value, z.value, theta.value

#2023-05-27 
#'''
#pym_RX, RXP, RXR, RY, RZ, TR to be done in VM
#'''   

def pymS00(date1,   date2,   x,   y):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    coordinates.  Compatible with IAU 2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    x : float
        CIP coordinates    
    y : float
        CIP coordinates

    Returns
    -------
    Tfunction value : float
        the CIO locator s in radians

    '''
    lib.iauS00.argtypes = [c_double,  c_double, c_double,  c_double ]
    lib.iauS00.restype  = c_double
    
    return  lib.iauS00(date1,   date2,   x,   y)

def pymS00a(date1,   date2):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2000A
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    lib.iauS00a.argtypes = [c_double,  c_double]
    lib.iauS00a.restype  = c_double
    
    return  lib.iauS00a(date1,   date2)

def pymS00b(date1,   date2):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2000B
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    lib.iauS00b.argtypes = [c_double,  c_double]
    lib.iauS00b.restype  = c_double
    
    return  lib.iauS00b(date1,   date2)


def pymS06(date1,   date2,   x,   y):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    coordinates.  Compatible with IAU 2006/2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    x : float
        CIP coordinates    
    y : float
        CIP coordinates

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    lib.iauS06.argtypes = [c_double,  c_double, c_double,  c_double ]
    lib.iauS06.restype  = c_double
    
    return  lib.iauS06(date1,   date2,   x,   y)


def pymS06a(date1,   date2):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2006
    precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    lib.iauS06a.argtypes = [c_double,  c_double]
    lib.iauS06a.restype  =  c_double
    
    return  lib.iauS06a(date1,   date2)

def pymSp00(date1,   date2):
    '''
    The TIO locator s', positioning the Terrestrial Intermediate Origin
    on the equator of the Celestial Intermediate Pole.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the TIO locator s' in radians

    '''
    lib.iauSp00.argtypes = [c_double,  c_double]
    lib.iauSp00.restype  =  c_double
    
    return  lib.iauSp00(date1,   date2)

def pymXy06(date1,   date2):
    '''
    X,Y coordinates of celestial intermediate pole from series based
    on IAU 2006 precession and IAU 2000A nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date     

    Returns
    -------
    x : float
        CIP X,Y coordinates    
    y : float
        CIP X,Y coordinates

    '''
    lib.iauXy06.argtypes = [c_double,  c_double,  c_double_p,  c_double_p]
    lib.iauXy06.restype  =  None
    x      = c_double()
    y      = c_double()
    
    lib.iauXy06(date1,   date2,  byref(x),  byref(y))
    
    return  x.value,  y.value

def pymXys00a(date1,   date2):
    '''
    For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2000A
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s
        
    '''
    lib.iauXys00a.argtypes = [c_double,    c_double,  
                             c_double_p,  c_double_p,  c_double_p]
    lib.iauXys00a.restype  =  None
    
    x      = c_double()
    y      = c_double()
    s      = c_double()
    
    lib.iauXys00a(date1,   date2,  byref(x),  byref(y),  byref(s))
    
    return  x.value,  y.value,  s.value 

def pymXys00b(date1,   date2):
    '''
    For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2000B
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s

    '''
    lib.iauXys00b.argtypes = [c_double,    c_double,  
                             c_double_p,  c_double_p,  c_double_p]
    lib.iauXys00b.restype  =  None
    
    x      = c_double()
    y      = c_double()
    s      = c_double()
    
    lib.iauXys00b(date1,   date2,  byref(x),  byref(y),  byref(s))
    
    return  x.value,  y.value,  s.value 

def pymXys06a(date1,   date2):
    '''
    For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2006
    precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s

    '''
    lib.iauXys06a.argtypes = [c_double,    c_double,  
                             c_double_p,  c_double_p,  c_double_p]
    lib.iauXys06a.restype  =  None
    
    x      = c_double()
    y      = c_double()
    s      = c_double()
    
    lib.iauXys06a(date1,   date2,  byref(x),  byref(y),  byref(s))
    
    return  x.value,  y.value,  s.value 

def pymBp00(date1,  date2):
    '''
    Frame bias and precession, IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rb : numpy.matrix(3,3)
        frame bias matrix     
    rp : numpy.matrix(3,3)   
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix
    '''
    lib.iauBp00.argtypes = [c_double,    c_double,  
                            vector_double3, vector_double3, vector_double3]
    lib.iauBp00.restype  = None
 
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauBp00(date1,  date2,  rb,  rp,  rbp)
    
    return   rb,  rp,  rbp 


def pymBp06(date1,  date2):
    '''
    Frame bias and precession, IAU 2006.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix

    '''
    lib.iauBp06.argtypes = [c_double,    c_double,  
                            vector_double3, vector_double3, vector_double3]
    lib.iauBp06.restype  = None
 
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauBp06(date1,  date2,  rb,  rp,  rbp)
    
    return   rb,  rp,  rbp  


def pymC2i06a(date1, date2):
    '''
    Form the celestial-to-intermediate matrix for a given date using the
    IAU 2006 precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    lib.iauC2i06a.argtypes = [c_double,  c_double, vector_double3]
    lib.iauC2i06a.restype  = None
    
    rc2i = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauC2i06a(date1, date2, rc2i)
    
    return  rc2i 


def pymC2ixy(date1, date2, x, y):
    '''
    Form the celestial to intermediate-frame-of-date matrix for a given
    date when the CIP X,Y coordinates are known.  IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    lib.iauC2ixy.argtypes = [c_double,  c_double, c_double,  c_double, vector_double3]
    lib.iauC2ixy.restype  = None
    
    rc2i = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauC2ixy(date1, date2, x, y, rc2i)
    
    return  rc2i 


def pymC2ixys(x, y, s):
    '''
    Form the celestial to intermediate-frame-of-date matrix given the CIP
    X,Y and the CIO locator s.

    Parameters
    ----------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s 

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    lib.iauC2ixys.argtypes = [c_double,  c_double, c_double,  vector_double3]
    lib.iauC2ixys.restype  = None
    
    rc2i = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauC2ixys(x, y, s, rc2i)
    
    return  rc2i 


def pymC2t06a(tta,  ttb,  uta,  utb,  xp,  yp):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1 and
    the polar motion, using the IAU 2006/2000A precession-nutation
    model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date     
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date     
    xp : float
        coordinates of the pole (radians)    
    yp : float
        coordinates of the pole (radians)

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    lib.iauC2t06a.argtypes = [c_double,  c_double, c_double,  c_double,
                              c_double,  c_double,
                              vector_double3]
    lib.iauC2t06a.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2t06a(tta,  ttb,  uta,  utb,  xp,  yp, rc2t)
    
    return  rc2t


def pymC2tcio(rc2i, era, rpom):
    '''
    Assemble the celestial to terrestrial matrix from CIO-based
    components (the celestial-to-intermediate matrix, the Earth Rotation
    Angle and the polar motion matrix).

    Parameters
    ----------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix     
    era : float
        Earth rotation angle (radians)     
    rpom : numpy.matrix(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''    
    lib.iauC2tcio.argtypes = [vector_double3,  c_double, vector_double3,
                              vector_double3]
    
    lib.iauC2tcio.restype  = None
    
    rc2t = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauC2tcio(rc2i, era, rpom, rc2t)
    
    return  rc2t  

def pymLtp(epj):
    '''
    Long-term precession matrix.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)     

    Returns
    -------
    rp : numpy.matrix(3,3)
        precession matrix, J2000.0 to date

    '''
    lib.iauLtp.argtypes = [c_double, vector_double3]
    lib.iauLtp.restype  =  None
    
    rp = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauLtp(epj, rp)
    
    return rp 


def pymLtpb(epj):
    '''
    Long-term precession matrix, including ICRS frame bias.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    

    Returns
    -------
    rpb : numpy.matrix(3,3)
        precession-bias matrix, J2000.0 to date

    '''
    lib.iauLtpb.argtypes = [c_double, vector_double3]
    lib.iauLtpb.restype  =  None
    
    rpb = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauLtpb(epj, rpb)
    
    return rpb


def pymLtpecl(epj):
    '''
    Long-term precession of the ecliptic.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    
        
    Returns
    -------
    vec : numpy.matrix(1,3)
        ecliptic pole unit vector    

    '''
    lib.iauLtpecl.argtypes = [c_double, vector_double]
    lib.iauLtpecl.restype  =  None
    
    vec = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauLtpecl(epj, vec)
    
    return vec


def pymLtpequ(epj):
    '''
    Long-term precession of the equator.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)   

    Returns
    -------
    veq : numpy.matrix(1,3)
        equator pole unit vector

    '''
    lib.iauLtpequ.argtypes = [c_double, vector_double]
    lib.iauLtpequ.restype  =  None
    
    veq = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauLtpequ(epj, veq)
    
    return veq   

def pymPfw06(date1,   date2):
    '''
    Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    gamb : flaot
        F-W angle gamma_bar (radians)    
    phib : float
        F-W angle phi_bar (radians)    
    psib : float
        F-W angle psi_bar (radians)    
    epsa : float
        F-W angle epsilon_A (radians)

    '''
    lib.iauPfw06.argtypes = [c_double,    c_double,   c_double_p,   c_double_p,
                             c_double_p,  c_double_p]
    lib.iauPfw06.restype  =  None
    
    gamb     = c_double()
    phib     = c_double()
    psib     = c_double()
    epsa     = c_double()
 
    
    lib.iauPfw06(date1,   date2,   
                 byref(gamb),  byref(phib), byref(psib),  byref(epsa))
    
    return  gamb.value,  phib.value,  psib.value,  epsa.value 

def pymPmat00(date1, date2): 
    '''
    Precession matrix (including frame bias) from GCRS to a specified
    date, IAU 2000 model.

    Parameters
    ----------
    date1 : flaot
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbp : numpy.matrix(3,3)
        bias-precession matrix

    '''
    lib.iauPmat00.argtypes = [c_double,  c_double,  vector_double3]
    lib.iauPmat00.restype  =  None
 
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPmat00(date1, date2,  rbp)
    
    return  rbp  


def pymPmat06(date1, date2): 
    '''
    Precession matrix (including frame bias) from GCRS to a specified
    date, IAU 2006 model.

    Parameters
    ----------
    date1 : flaot
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbp : numpy.matrix(3,3)
        bias-precession matrix

    '''
    lib.iauPmat06.argtypes = [c_double,  c_double,  vector_double3]
    lib.iauPmat06.restype  =  None
 
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPmat06(date1, date2,  rbp)
    
    return  rbp    


def pymPmat76(date1, date2): 
    '''
    Precession matrix from J2000.0 to a specified date, IAU 1976 model.

    Parameters
    ----------
    date1 : float
        ending date, TT    
    date2 : float
        ending date, TT

    Returns
    -------
    rmatp : numpy.matrix(3,3)
        precession matrix, J2000.0 -> date1+date2

    '''    
    lib.iauPmat76.argtypes = [c_double,  c_double,  vector_double3]
    lib.iauPmat76.restype  =  None
 
    rmatp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPmat76(date1, date2,  rmatp)
    
    return  rmatp        


def pymPn00(date1,  date2,  dpsi,  deps):
    '''
    Precession-nutation, IAU 2000 model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based
    use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    dpsi : float
        nutation    
    deps : float
        nutation

    Returns
    -------
    epsa : flaot
        mean obliquity    
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix
        
    '''    
    lib.iauPn00.argtypes = [ c_double,    c_double,   c_double,   c_double, 
                             c_double_p,    
                             vector_double3, vector_double3, vector_double3,
                             vector_double3, vector_double3 
                            ]
    lib.iauPn00.restype  = None
    
    
    epsa    = c_double()
    
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rn   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPn00(date1,  date2,   dpsi,  deps,  byref(epsa),  
                 rb,  rp,  rbp,  rn,  rbpn)
    
    return   epsa.value,  rb,  rp,  rbp,  rn,  rbpn     


def pymPn00a(date1,  date2):
    '''
    Precession-nutation, IAU 2000A model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based
    use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    

    Returns
    -------
    dpsi : float
        nutation    
    deps : float
        nutation    
    epsa : flaot
        mean obliquity    
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    lib.iauPn00a.argtypes = [c_double,    c_double,  
                             c_double_p,  c_double_p,  c_double_p,  
                             vector_double3, vector_double3, vector_double3,
                             vector_double3, vector_double3 
                            ]
    lib.iauPn00a.restype  = None
    
    dpsi    = c_double()
    deps    = c_double()
    epsa    = c_double()
    
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rn   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPn00a(date1,  date2,  byref(dpsi),  byref(deps),  byref(epsa),  
                 rb,  rp,  rbp,  rn,  rbpn)
    
    return  dpsi.value,   deps.value,   epsa.value,  rb,  rp,  rbp,  rn,  rbpn    


def pymPn00b(date1,  date2):
    '''
    Precession-nutation, IAU 2000B model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based
    use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    

    Returns
    -------
    dpsi : float
        nutation    
    deps : float
        nutation    
    epsa : flaot
        mean obliquity    
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    lib.iauPn00b.argtypes = [c_double,    c_double,  
                             c_double_p,  c_double_p,  c_double_p,  
                             vector_double3, vector_double3, vector_double3,
                             vector_double3, vector_double3 
                            ]
    lib.iauPn00b.restype  = None
    
    dpsi    = c_double()
    deps    = c_double()
    epsa    = c_double()
    
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rn   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPn00b(date1,  date2,  byref(dpsi),  byref(deps),  byref(epsa),  
                 rb,  rp,  rbp,  rn,  rbpn)
    
    return  dpsi.value,   deps.value,   epsa.value,  rb,  rp,  rbp,  rn,  rbpn    


def pymPn06(date1,  date2,  dpsi,  deps):
    '''
    Precession-nutation, IAU 2006 model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based use
    indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date       
    dpsi : float
        nutation    
    deps : float
        nutation    

    Returns
    -------
    epsa : flaot
        mean obliquity    
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    lib.iauPn06.argtypes = [ c_double,    c_double,   c_double,   c_double, 
                             c_double_p,    
                             vector_double3, vector_double3, vector_double3,
                             vector_double3, vector_double3 
                            ]
    lib.iauPn06.restype  = None
    
    epsa    = c_double()
    
    rb   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rp   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbp  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rn   = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    rbpn = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauPn06(date1,  date2,   dpsi,  deps,  byref(epsa),  
                 rb,  rp,  rbp,  rn,  rbpn)
    
    return   epsa.value,  rb,  rp,  rbp,  rn,  rbpn   

def pymPlan94(date1, date2, npp): 
    '''
    Approximate heliocentric position and velocity of a nominated major
    planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
    Neptune (but not the Earth itself).

    Parameters
    ----------
    date1 : float
        TDB date part A    
    date2 : float
        TDB date part B    
    npp : int
        planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,
        5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)

    Returns
    -------
    pv : numpy.matrix(2,3)
        planet p,v (heliocentric, J2000.0, au,au/d)

    '''
    lib.iauPlan94.argtypes = [c_double, c_double, c_int, vector_double2]
    lib.iauPlan94.restype  =  c_int
 
    pv  = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauPlan94(date1, date2, npp, pv)
    
    return pv   

def pymMoon98(date1, date2): 
    '''
    Approximate geocentric position and velocity of the Moon.

    Parameters
    ----------
    date1 : float
        TT date part A    
    date2 : float
        TT date part A

    Returns
    -------
    pv : numpy.matrix(2,3)
        Moon p,v, GCRS (AU, AU/d)

    '''
    lib.iauMoon98.argtypes = [c_double, c_double, vector_double2]
    lib.iauMoon98.restype  =  None
 
    pv  = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauMoon98(date1, date2, pv)
    
    return  pv 