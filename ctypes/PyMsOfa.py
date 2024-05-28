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
#lib       = CDLL(libs_path)


print(pf.system().lower())
#'''
#for windows
if   pf.system().lower() == 'windows':
          libs_file  = 'libsofa_c.dll'

#for linux
elif pf.system().lower() == 'linux':
          dirname = os.getcwd()
          #lib = CDLL('./libsofa_c.so')

#for macOS
elif pf.system().lower() == 'darwin':
          pass


libs_path  = os.path.join(dirname, libs_file)
lib  = CDLL(libs_path)

#'''

c_int4 = c_int * 4

def pymA2af(ndp, angle):
    '''
    Decompose radians into degrees, arcminutes, arcseconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution      
    angle : float
        angle in radians

    Returns
    -------
    sign : bytes
        '+' or '-'    
    idmsf : tuple
        degrees, arcminutes, arcseconds, fraction 

         NDP         resolution
          :      ...0000 00 00
         -7         1000 00 00
         -6          100 00 00
         -5           10 00 00
         -4            1 00 00
         -3            0 10 00
         -2            0 01 00
         -1            0 00 10
          0            0 00 01
          1            0 00 00.1
          2            0 00 00.01
          3            0 00 00.001
          :            0 00 00.000...
    '''
    
     
    lib.iauA2af.argtypes = [c_int, c_double, POINTER(c_char), c_int4] 
    lib.iauA2af.restype  = None
    
    sign   = c_char()
    idmsf  = (c_int4)()
     
    lib.iauA2af(ndp, float(angle), byref(sign), idmsf)

    return sign.value, tuple([x for x in idmsf])


def pymPas(al, ap, bl, bp):
    '''
    Position-angle from spherical coordinates. 

    Parameters
    ----------
    al : float
        longitude of point A (e.g. RA) in radians     
    ap : float
        latitude of point A (e.g. Dec) in radians    
    bl : float
        longitude of point B     
    bp : float
        latitude of point B     

    Returns
    -------
    function value : float
        position angle of B with respect to A

    '''
    lib.iauPas.argtypes = [c_double, c_double, c_double, c_double]
    lib.iauPas.restype  = c_double
    
    return lib.iauPas(al, ap, bl, bp)

array_1d_double = nt.ndpointer(shape=(1,3), dtype=np.double, flags='C') 
vector_double   = nt.ndpointer(shape=(1,3), dtype=np.double, flags='C') 
array_double_2  = nt.ndpointer(shape=(1,2), dtype=np.double, flags='C') 

def pymS2c_A(theta, phi, c):
    
    lib.iauS2c.argtypes = [c_double, c_double, vector_double]
    lib.iauS2c.restype  = None
    lib.iauS2c(theta, phi, c)
    return 

def pymS2c(theta, phi):
    '''
    Convert spherical coordinates to Cartesian.

    Parameters
    ----------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)

    Returns
    -------
    c : numpy.matrix
        direction cosines

    '''

    lib.iauS2c.argtypes = [c_double, c_double, vector_double]
    lib.iauS2c.restype  = None
    c = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    lib.iauS2c(theta, phi, c)
    return c

#int iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
pym_cal2jd_msg = {
                 -1: 'bad year,  the year is simply valid from -4800 March 1',
                 -2: 'bad month, the month is not from 1 to 12',
                 -3: 'bad day,   the day is not related to the month'
                 } 

c_double_p = POINTER(c_double)   
def pymCal2jd_A(iyear, imon, iday, djm0, djm):
    
    lib.iauCal2jd.argtypes = [c_int, c_int, c_int, c_double_p,  c_double_p]
    lib.iauCal2jd.restype  = c_int
    
    j = lib.iauCal2jd(iyear, imon, iday, byref(djm0), byref(djm)) 
    if j != 0:
        raise ValueError(pym_cal2jd_msg[j])
    return    

def pymCal2jd(iyear, imon, iday):
    '''
    Gregorian Calendar to Julian Date.

    Parameters
    ----------
    iyear : int
        year in Gregorian calendar
    imon : int
        month in Gregorian calendar
    iday : int
        day in Gregorian calendar

    Raises
    ------
    ValueError
        -1: 'bad year,  the year is simply valid from -4800 March 1',
        -2: 'bad month, the month is not from 1 to 12',
        -3: 'bad day,   the day is not related to the month'

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5
    djm : float
        Modified Julian Date for 0 hrs

    '''
    lib.iauCal2jd.argtypes = [c_int, c_int, c_int, c_double_p,  c_double_p]
    lib.iauCal2jd.restype  = c_int    
    djm0 = c_double()
    djm  = c_double()    
    j    = lib.iauCal2jd(iyear, imon, iday, byref(djm0), byref(djm))
    if j != 0:
        raise ValueError(pym_cal2jd_msg[j]) 
    return  djm0.value, djm.value  
  
c_int_p = POINTER(c_int)
pym_jd2cal_msg ={
                -1: 'The valid date is -68569.5 (-4900 March 1) up to 1e9'
                }  

def pymJd2cal(dj1, dj2):
    '''
    Julian Date to Gregorian year, month, day, and fraction of a day.

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Raises
    ------
    ValueError
        -1: 'The valid date is -68569.5 (-4900 March 1) up to 1e9'

    Returns
    -------
    iy : int
        yaer    
    im : int
        month
    id : int
        day
    fd : float
        fraction of day

    '''
    lib.iauJd2cal.argtypes = [c_double, c_double, c_int_p, 
                             c_int_p, c_int_p, c_double_p ] 
    lib.iauJd2cal.restype  = c_int 
    iy   = c_int()
    im   = c_int()
    iday = c_int()
    fday = c_double()
    
    j = lib.iauJd2cal(dj1, dj2, byref(iy), byref(im), byref(iday), byref(fday))
    if j != 0:
        raise ValueError(pym_jd2cal_msg[j])
    return iy.value, im.value, iday.value, fday.value

#int iauJdcalf(int ndp, double dj1, double dj2, int iymdf[4])
def pymJdcalf(ndp, dj1, dj2):  
    '''
    Julian Date to Gregorian Calendar, expressed in a form convenient 
    for formatting messages:  rounded to a specified precision.

    Parameters
    ----------
    ndp : int
        number of decimal places of days in fraction    
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date    

    Returns
    -------
    iymdf : tuple
        year, month, day, fraction in Gregorian calendar

    '''
    lib.iauJdcalf.argtypes = [c_int, c_double, c_double]
    lib.iauJdcalf.restype  = c_int 
    idmsf  = (c_int4)()
    
    lib.iauJdcalf(ndp, dj1,  dj2, idmsf)
    return  tuple([x for x in idmsf])

#void iauD2tf(int ndp, double days, char *sign, int ihmsf[4])
def pymD2tf(ndp, days):
    '''
    Decompose days to hours, minutes, seconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution    
    days : float
        interval in days

    Returns
    -------
    sign : 'bytes'
        '+' or '-'    
    ihmsf : tuple    
        hours, minutes, seconds, fraction

         NDP         resolution
          :      ...0000 00 00
         -7         1000 00 00
         -6          100 00 00
         -5           10 00 00
         -4            1 00 00
         -3            0 10 00
         -2            0 01 00
         -1            0 00 10
          0            0 00 01
          1            0 00 00.1
          2            0 00 00.01
          3            0 00 00.001
          :            0 00 00.000...

    '''
    lib.iauD2tf.argtypes = [c_int, c_double, POINTER(c_char), c_int4] 
    lib.iauD2tf.restype  = None
    
    sign   = c_char()
    ihmsf  = (c_int4)()
    
    lib.iauD2tf(ndp, float(days), byref(sign), ihmsf)
    return sign.value, tuple([x for x in ihmsf])

#int iauD2dtf(const char *scale, int ndp, double d1, double d2,
#             int *iy, int *im, int *id, int ihmsf[4])
pym_d2dtf_msg = {
                 1: 'dubious year',
                -1: 'unacceptable date',
                }

def pymD2dtf(scale, ndp, d1, d2):
    '''
    Format for output a 2-part Julian Date (or in the case of UTC a 
    quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : ctypes.c_char_p
        time scale ID(Only the value "UTC" is significant)       
    ndp : int
        resolution    
    d1 : float
        time as a 2-part Julian Date    
    d2 : float
        time as a 2-part Julian Date

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    iy,im,id : int
        year, month, day in Gregorian calendar    
    ihmsf : tuple
        hours, minutes, seconds, fraction

         NDP         resolution
          :      ...0000 00 00
         -7         1000 00 00
         -6          100 00 00
         -5           10 00 00
         -4            1 00 00
         -3            0 10 00
         -2            0 01 00
         -1            0 00 10
          0            0 00 01
          1            0 00 00.1
          2            0 00 00.01
          3            0 00 00.001
          :            0 00 00.000...
    '''
    
    lib.iauD2dtf.argtypes = [c_wchar_p, c_int, c_double, c_double,
                             c_int_p, c_int_p, c_int_p, c_int4] 
    lib.iauD2dtf.restype  = c_int
    
    iy     = c_int()
    im     = c_int()
    iday   = c_int()
    ihmsf  =(c_int4)()
     
    j=lib.iauD2dtf(scale, ndp, d1, d2, byref(iy), byref(im), byref(iday), ihmsf)
    if j < 0:
         raise ValueError(pym_d2dtf_msg[j])
    elif j > 0:
         ws.warn(pym_d2dtf_msg[j], UserWarning, 2)
    #return iy.value, im.value, iday.value, tuple([x for x in ihmsf])
    return (iy.value, im.value, iday.value) + tuple([x for x in ihmsf])

#int iauDat(int iy, int im, int id, double fd, double *deltat)
#Note: current leap second is 2021, this should update and recompile iauDat for new date    
pym_dat_msg =   {
                  1: 'dubious year', 
                 -1: 'bad year,  the year is simply valid from -4800 March 1',
                 -2: 'bad month, the month is not from 1 to 12',
                 -3: 'bad day,   the day is not related to the month',
                 -4: 'bad fraction of day',
                 -5: 'internal error', 
                 } 

def pymDat(iy, im, iday, fday):
    '''
    For a given UTC date, calculate Delta(AT) = TAI-UTC.

    Parameters
    ----------
    iy : int
        year    
    im : int
        month    
    iday : int
        day    
    fday : float    
        fraction of day

    Raises
    ------
    ValueError
        1: 'dubious year', 
       -1: 'bad year,  the year is simply valid from -4800 March 1',
       -2: 'bad month, the month is not from 1 to 12',
       -3: 'bad day,   the day is not related to the month',
       -4: 'bad fraction of day',
       -5: 'internal error', 

    Returns
    -------
    deltat : float
        TAI minus UTC, seconds

    '''
    lib.iauDat.argtypes = [c_int, c_int, c_int, c_double]
    lib.iauDat.restype  = c_int
    deltat = c_double()
    j = lib.iauDat(iy, im, iday, fday, byref(deltat))
    if   j < 0:
          raise ValueError(pym_dat_msg[j])
    elif j > 0:
          ws.warn(pym_dat_msg[j], UserWarning, 2)    
    return deltat.value  

#double iauDtdb(double date1, double date2,
#               double ut, double elong, double u, double v)   

def pymDtdb(date1, date2, ut, elong, u, v):
    '''
    An approximation to TDB-TT, the difference between barycentric
    dynamical time and terrestrial time, for an observer on the Earth.

    Parameters
    ----------
    date1 : float
        date, TDB    
    date2 : float
        date, TDB    
    ut : float
        universal time (UT1, fraction of one day)    
    elong : float
        longitude (east positive, radians)     
    u : float
        distance from Earth spin axis (km)    
    v : float
        distance north of equatorial plane (km)

    Returns
    -------
    function value : float
        TDB-TT (seconds)

    '''
    lib.iauDtdb.argtypes = [c_double, c_double, c_double, 
                            c_double, c_double, c_double]
    lib.iauDtdb.restype  =  c_double
    
    return lib.iauDtdb(date1, date2, ut, elong, u, v)

#int iauDtf2d(const char *scale, int iy, int im, int id,           
#             int ihr, int imn, double sec, double *d1, double *d2)

pym_dtf2d_msg = {
                  3: 'both of next two',
                  2: 'time is after end of day',
                  1: 'dubious year', 
                 -1: 'bad year,  the year is simply valid from -4800 March 1',
                 -2: 'bad month, the month is not from 1 to 12',
                 -3: 'bad day,   the day is not related to the month',
                 -4: 'bad hour',
                 -5: 'bad minute',
                 -6: 'bad second', 
                } 

def pymDtf2d(scale,  iy,  im,  id,  ihr, imn,  sec):
    '''
    Encode date and time fields into 2-part Julian Date (or in the case
    of UTC a quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : ctypes.c_char_p
        time scale ID(Only the value "UTC" is significant)
    iy : int
        year in Gregorian calendar
    im : int
        month in Gregorian calendar
    id : int
        day in Gregorian calendar
    ihr : int
        hour
    imn : int
        minute
    sec : float
        seconds

    Raises
    ------
    ValueError
        3: 'both of next two',
        2: 'time is after end of day',
        1: 'dubious year', 
       -1: 'bad year,  the year is simply valid from -4800 March 1',
       -2: 'bad month, the month is not from 1 to 12',
       -3: 'bad day,   the day is not related to the month',
       -4: 'bad hour',
       -5: 'bad minute',
       -6: 'bad second', 

    Returns
    -------
    d1,d2 : float
        2-part Julian Date

    '''
    lib.iauDtf2d.argtypes = [c_wchar_p, c_int, c_int, c_int, 
                             c_int, c_int, c_double, c_double_p , 
                             c_double_p ]
    lib.iauDtf2d.restype  = c_int
    d1  = c_double()
    d2  = c_double()  
    j=lib.iauDtf2d(scale, iy, im, id, ihr, imn, sec, byref(d1), byref(d2))
    if j < 0:
         raise ValueError(pym_dtf2d_msg[j])
    elif j > 0:
         ws.warn(pym_dtf2d_msg[j], UserWarning, 2)
    return d1.value, d2.value

#double iauEpb(double dj1, double dj2)
def pymEpb(dj1, dj2):
    '''
    Julian Date to Besselian Epoch.

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Returns
    -------
    function value : float
        Besselian Epoch

    '''
    lib.iauEpb.argtypes = [c_double,  c_double]
    lib.iauEpb.restype  = c_double
   
    return lib.iauEpb(dj1, dj2)

#directly calculating for comparison
D1900 = 36524.68648E0 
DJ00  = 2451545.0E0 
DTY   = 365.242198781E0    
def pymEpb_A(dj1, dj2):
 
    return 1900.0E0 + ((dj1 - DJ00) + (dj2 + D1900)) / DTY

#void iauEpb2jd(double epb, double *djm0, double *djm)
def pymEpb2jd(epb):
    '''
    Besselian Epoch to Julian Date

    Parameters
    ----------
    epb : float
        Besselian Epoch

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5    
    djm : float
        Modified Julian Date

    '''
    lib.iauEpb2jd.argtypes = [c_double, c_double_p , c_double_p ]
    lib.iauEpb2jd.restype  = None
    djm0 = c_double()
    djm  = c_double()
    lib.iauEpb2jd(epb, byref(djm0), byref(djm))
    return djm0.value, djm.value

DJM0 = 2400000.5E0
def pymEpb2jd_A(epb):
    
    djm0 = DJM0
    djm  = 15019.81352E0 + (epb - 1900.0E0) * DTY 
    return djm0, djm

#double iauEpj(double dj1, double dj2)    
def pymEpj(dj1, dj2):
    '''
    Julian Date to Julian Epoch

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Returns
    -------
    function value : float
        Julian Epoch

    '''
    lib.iauEpj.argtypes = [c_double,  c_double]
    lib.iauEpj.restype  = c_double
   
    return lib.iauEpj(dj1, dj2)    

DJY = 365.25E0
def pymEpj_A(dj1, dj2):
    
    epj = 2000.0E0 + ((dj1 - DJ00) + dj2) / DJY 
    return epj   

def pymEpj2jd(epj):
    '''
    Julian Epoch to Julian Date

    Parameters
    ----------
    epj : float
        Julian Epoch (e.g. 1996.8)

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5     
    djm : float
        Modified Julian Date

    '''
    lib.iauEpj2jd.argtypes = [c_double, c_double_p , c_double_p ]
    lib.iauEpj2jd.restype  = None
    djm0 = c_double()
    djm  = c_double()
    lib.iauEpj2jd(epj, byref(djm0), byref(djm))
    
    return djm0.value, djm.value
 
DJM00 = 51544.5E0
def pymEpj2jd_A(epj):
    
    djm0 = DJM0
    djm  = DJM00 + (epj - 2000.0E0) * 365.25E0 
    return djm0, djm

#int iauTaitt(double tai1, double tai2, double *tt1, double *tt2)
def pymTaitt(tai1, tai2):
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Terrestrial Time, TT.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date    
    tai2 : float
        TAI as a 2-part Julian Date

    Returns
    -------
    tt1,tt2 : float
        TT as a 2-part Julian Date

    '''
    lib.iauTaitt.argtypes = [c_double, c_double, c_double_p, c_double_p]
    lib.iauTaitt.restype  = c_int
    tt1  = c_double()
    tt2  = c_double()
    j = lib.iauTaitt(tai1, tai2, byref(tt1), byref(tt2))
    return tt1.value, tt2.value

#int iauTaiut1(double tai1, double tai2, double dta,
#              double *ut11, double *ut12)
    
def pymTaiut1(tai1, tai2, dta):
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Universal Time, UT1.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date    
    tai2 : float
        TAI as a 2-part Julian Date    
    dta : float
        UT1-TAI in seconds

    Returns
    -------
    ut11,ut12 : float
        UT1 as a 2-part Julian Date

    '''
    lib.iauTaiut1.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauTaiut1.restype  = c_int
    ut11 = c_double()
    ut12 = c_double()
    j = lib.iauTaiut1(tai1, tai2, dta, byref(ut11), byref(ut12))
    return ut11.value, ut12.value

#int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2)
pym_taiutc_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }    
def pymTaiutc(tai1, tai2): 
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Coordinated Universal Time, UTC.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    '''
    lib.iauTaiutc.argtypes = [c_double, c_double,  c_double_p, c_double_p]
    lib.iauTaiutc.restype  = c_int
    utc1 = c_double()
    utc2 = c_double()
    j = lib.iauTaiutc(tai1, tai2, byref(utc1), byref(utc2))
    if   j < 0:
        raise ValueError(pym_taiutc_msg[j])
    elif j > 0:
        ws.warn(pym_taiutc_msg[j], UserWarning, 2)
    return utc1.value, utc2.value

#int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)
def pymTcbtdb(tcb1, tcb2):
    '''
    Time scale transformation:  Barycentric Coordinate Time, TCB, to
    Barycentric Dynamical Time, TDB.

    Parameters
    ----------
    tcb1 : float
        TCB as a 2-part Julian Date
    tcb2 : float
        TCB as a 2-part Julian Date

    Returns
    -------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 :float
        TDB as a 2-part Julian Date

    '''
    lib.iauTcbtdb.argtypes = [c_double, c_double,  c_double_p, c_double_p]
    lib.iauTcbtdb.restype  = c_int
    tdb1 = c_double()
    tdb2 = c_double()
    j    = lib.iauTcbtdb(tcb1, tcb2, byref(tdb1), byref(tdb2))
   
    return tdb1.value, tdb2.value   

#int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2) 
def pymTcgtt(tcg1, tcg2):  
    '''
    Time scale transformation:  Geocentric Coordinate Time, TCG, to
    Terrestrial Time, TT.

    Parameters
    ----------
    tcg1 : float
        TCG as a 2-part Julian Date
    tcg2 : float
        TCG as a 2-part Julian Date

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    '''    
    lib.iauTcgtt.argtypes = [c_double, c_double,  c_double_p, c_double_p]
    lib.iauTcgtt.restype  = c_int
    
    tt1 = c_double()
    tt2 = c_double()
    j   = lib.iauTcgtt(tcg1, tcg2, byref(tt1), byref(tt2))
   
    return tt1.value, tt2.value       


#int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)
def pymTdbtcb(tdb1, tdb2):  
    '''
    Time scale transformation:  Barycentric Dynamical Time, TDB, to
    Barycentric Coordinate Time, TCB.

    Parameters
    ----------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 : float
        TDB as a 2-part Julian Date

    Returns
    -------
    tcb1 : float
        TCB as a 2-part Julian Date
    tcb2 : float
        TCB as a 2-part Julian Date

    '''    
    lib.iauTdbtcb.argtypes = [c_double, c_double,  c_double_p, c_double_p]
    lib.iauTdbtcb.restype  = c_int
    
    tcb1 = c_double()
    tcb2 = c_double()
    j    = lib.iauTdbtcb(tdb1, tdb2, byref(tcb1), byref(tcb2))
   
    return tcb1.value, tcb2.value 


#int iauTdbtt(double tdb1, double tdb2, double dtr,
#             double *tt1, double *tt2 )      
def pymTdbtt(tdb1, tdb2, dtr):   
    '''
    Time scale transformation:  Barycentric Dynamical Time, TDB, to
    Terrestrial Time, TT.

    Parameters
    ----------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 : float
        TDB as a 2-part Julian Date
    dtr : float
        TDB-TT in seconds

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    '''
    lib.iauTdbtt.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauTdbtt.restype  = c_int
    tt1  = c_double()
    tt2  = c_double()
    j = lib.iauTdbtt(tdb1, tdb2, dtr, byref(tt1), byref(tt2))
    return tt1.value, tt2.value

#int iauTf2d(char s, int ihour, int imin, double sec, double *days)
pym_tf2d_msg = {
                1:'ihour outside range 0-23',
                2:'imin outside range 0-59',
                3:'sec outside range 0-59.999...'
                }    
def pymTf2d(s, ihour, imin, sec): 
    '''
    Convert hours, minutes, seconds to days

    Parameters
    ----------
    s : bytes
        sign:  '-' = negative, otherwise positive
    ihour : int
        hours
    imin : int
        minutes
    sec : float
        seconds
    
    Raises
    ------
    ValueError
        1:'ihour outside range 0-23',
        2:'imin outside range 0-59',
        3:'sec outside range 0-59.999...'

    Returns
    -------
    days : float
        interval in days

    '''
#Note: for single byte s, here using c_wchar that can work.  
    lib.iauTf2d.argtypes = [c_wchar, c_int, c_int, c_double, c_double_p]
    lib.iauTf2d.restype  = c_int
    days = c_double()
     
    j = lib.iauTf2d(s, ihour, imin, sec, byref(days))
    if j > 0:
        ws.warn(pym_tf2d_msg[j], UserWarning, 2)
    return days.value 

#int iauTttai(double tt1, double tt2, double *tai1, double *tai2)
def pymTttai(tt1, tt2):
    '''
    Time scale transformation:  Terrestrial Time, TT, to International
    Atomic Time, TAI.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    '''
    lib.iauTttai.argtypes = [c_double, c_double, c_double_p, c_double_p]
    lib.iauTttai.restype  = c_int
    tai1  = c_double()
    tai2  = c_double()
    j = lib.iauTttai(tt1, tt2, byref(tai1), byref(tai2))
    return tai1.value, tai2.value


#int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2)
def pymTttcg(tt1, tt2):  
    '''
    Time scale transformation:  Terrestrial Time, TT, to Geocentric
    Coordinate Time, TCG.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    tcg1 : float
        TCG as a 2-part Julian Date
    tcg2 : float
        TCG as a 2-part Julian Date

    '''
    
    lib.iauTttcg.argtypes = [c_double, c_double,  c_double_p, c_double_p]
    lib.iauTttcg.restype  = c_int
    
    tcg1 = c_double()
    tcg2 = c_double()
    j   = lib.iauTttcg(tt1, tt2, byref(tcg1), byref(tcg2))
   
    return tcg1.value, tcg2.value       

#int iauTttdb(double tt1, double tt2, double dtr,
#             double *tdb1, double *tdb2)
def pymTttdb(tt1, tt2, dtr):  
    '''
    Time scale transformation:  Terrestrial Time, TT, to Barycentric
    Dynamical Time, TDB.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date
    dtr : float
        TDB-TT in seconds

    Returns
    -------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 : float
        TDB as a 2-part Julian Date

    '''
    lib.iauTttdb.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauTttdb.restype  = c_int
    
    tdb1 = c_double()
    tdb2 = c_double()
    j   = lib.iauTttdb(tt1, tt2, dtr, byref(tdb1), byref(tdb2))
   
    return tdb1.value, tdb2.value  

#int iauTtut1(double tt1, double tt2, double dt,
#             double *ut11, double *ut12)         
def pymTtut1(tt1, tt2, dtr):  
    '''
    Time scale transformation:  Terrestrial Time, TT, to Universal Time,
    UT1.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date
    dtr : float
        TT-UT1 in seconds

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date

    '''
    lib.iauTtut1.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauTtut1.restype  = c_int
    
    ut11 = c_double()
    ut12 = c_double()
    j   = lib.iauTtut1(tt1, tt2, dtr, byref(ut11), byref(ut12))
   
    return ut11.value, ut12.value

#int iauUt1tai(double ut11, double ut12, double dta,
#              double *tai1, double *tai2)   
#The argument dta, i.e. UT1-TAI, is an observed quantity, and is
#     available from IERS tabulations.   
def pymUt1tai(ut11, ut12, dta):  
    '''
    Time scale transformation:  Universal Time, UT1, to International
    Atomic Time, TAI.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    dta : float
        UT1-TAI in seconds

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    '''
    lib.iauUt1tai.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauUt1tai.restype  = c_int
    
    tai1 = c_double()
    tai2 = c_double()
    j   = lib.iauUt1tai(ut11, ut12, dta, byref(tai1), byref(tai2))
   
    return tai1.value, tai2.value    

#int iauUt1tt(double ut11, double ut12, double dt,
#             double *tt1, double *tt2)
def pymUt1tt(ut11, ut12, dt):  
    '''
    Time scale transformation:  Universal Time, UT1, to Terrestrial
    Time, TT.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    dt : float
        TT-UT1 in seconds

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    '''
    lib.iauUt1tt.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauUt1tt.restype  = c_int
    
    tt1 = c_double()
    tt2 = c_double()
    j   = lib.iauUt1tt(ut11, ut12, dt, byref(tt1), byref(tt2))
   
    return tt1.value, tt2.value  

#int iauUt1utc(double ut11, double ut12, double dut1,
#              double *utc1, double *utc2)     
pym_ut1utc_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }       
def pymUt1utc(ut11, ut12, dut1):  
    '''
    Time scale transformation:  Universal Time, UT1, to Coordinated
    Universal Time, UTC.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    dut1 : float
        Delta UT1: UT1-UTC in seconds

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    '''
    lib.iauUt1utc.argtypes = [c_double, c_double, c_double, c_double_p, c_double_p]
    lib.iauUt1utc.restype  = c_int
    
    utc1 = c_double()
    utc2 = c_double()
    j    = lib.iauUt1utc(ut11, ut12, dut1, byref(utc1), byref(utc2))
    if   j < 0:
        raise ValueError(pym_ut1utc_msg[j])
    elif j > 0:
        ws.warn(pym_ut1utc_msg[j], UserWarning, 2)
    return utc1.value, utc2.value      

#int iauUtctai(double utc1, double utc2, double *tai1, double *tai2)
pym_utctai_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }    
def pymUtctai(utc1, utc2):  
    '''
    Time scale transformation:  Coordinated Universal Time, UTC, to
    International Atomic Time, TAI.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    '''
    lib.iauUtctai.argtypes = [c_double, c_double,  c_double_p, c_double_p]
    lib.iauUtctai.restype  = c_int
    tai1 = c_double()
    tai2 = c_double()
    j = lib.iauUtctai(utc1, utc2, byref(tai1), byref(tai2))
    if   j < 0:
        raise ValueError(pym_utctai_msg[j])
    elif j > 0:
        ws.warn(pym_utctai_msg[j], UserWarning, 2)
    return tai1.value, tai2.value  

#int iauUtcut1(double utc1, double utc2, double dut1,
#              double *ut11, double *ut12)
pym_utcut1_msg = pym_ut1utc_msg
def pymUtcut1(utc1, utc2, dut1):  
    '''
    Time scale transformation:  Coordinated Universal Time, UTC, to
    Universal Time, UT1.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date
    dut1 : float
        Delta UT1 = UT1-UTC in seconds

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date

    '''
    lib.iauUtcut1.argtypes = [c_double, c_double, c_double,  c_double_p, c_double_p]
    lib.iauUtcut1.restype  = c_int
    ut11 = c_double()
    ut12 = c_double()
    j = lib.iauUtcut1(utc1, utc2, dut1, byref(ut11), byref(ut12))
    if   j < 0:
        raise ValueError(pym_utcut1_msg[j])
    elif j > 0:
        ws.warn(pym_utcut1_msg[j], UserWarning, 2)
    return ut11.value, ut12.value    


#2023-03-08 SOFA Astrometry Tools 

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

        
    def printA(self):    
    #def printA(pymASTROM):
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
} pymLDBODY;                                                           
                                                                    
'''
#2023-03-12, pymLDBODY class from C structure 
class pymLDBODY(ctypes.Structure):
    _fields_ = [("bm:",    c_double),
                ("dl:",    c_double),
                ("pv:",    c_dbl32),
                ]
    
    def __init__(self, bm, dl, pv):
   
        self.bm  = bm
        self.dl  = dl
        self.pv  = pv 
        
    #@property     
    def printA(pymLDBODY):
         pv_c  = np.array([x for x in pymLDBODY.pv])
         print("bm:", pymLDBODY.bm)
         print("dl:", pymLDBODY.dl)
         print("pv:", pv_c)
    
  
#void iauAtciqn(double rc, double dc, double pr, double pd,
#               double px, double rv, pymASTROM *astrom,
#               int n, pymLDBODY b[], double *ri, double *di)

#using 2-dim array to replace pymLDBODY structure, it works
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
    b : pymLDBODY class
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
#it now work well for 2-dim of b rather than the structure pymLDBODY      
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
    b : pymLDBODY class
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
    b : pymLDBODY class
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

#2023-05-03   Sofa Tools for Earth Attitude 
###sofa_pn_f.pdf###

def pymAnp(a):
    '''
    Normalize angle into the range 0 <= a < 2pi.

    Parameters
    ----------
    a : float
        angle (radians)

    Returns
    -------
    function value : float
        angle in range 0-2pi

    '''
    lib.iauAnp.argtypes = [c_double]
    lib.iauAnp.restype  = c_double
    
    return lib.iauAnp(a)

#2023-05-04     
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

def pymGmst00(uta,  utb,  tta,  ttb):
    '''
    Greenwich mean sidereal time (model consistent with IAU 2000
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich mean sidereal time (radians)

    '''
    lib.iauGmst00.argtypes =[c_double, c_double, c_double, c_double]
    
    lib.iauGmst00.restype  = c_double
    
    return  lib.iauGmst00(uta,  utb,  tta,  ttb)

def pymGmst06(uta,  utb,  tta,  ttb):
    '''
    Greenwich mean sidereal time (consistent with IAU 2006 precession).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich mean sidereal time (radians)

    '''
    lib.iauGmst06.argtypes =[c_double, c_double, c_double, c_double]
    
    lib.iauGmst06.restype  = c_double
    
    return  lib.iauGmst06(uta,  utb,  tta,  ttb)

def pymGmst82(dj1,  dj2):
    '''
    Universal Time to Greenwich mean sidereal time (IAU 1982 model).

    Parameters
    ----------
    dj1 : float
        UT1 Julian Date    
    dj2 : float
        UT1 Julian Date

    Returns
    -------
    function value : float
        Greenwich mean sidereal time (radians)

    '''
    lib.iauGmst82.argtypes =[c_double, c_double]
    
    lib.iauGmst82.restype  = c_double
    
    return  lib.iauGmst82(dj1,  dj2)

def pymGst00a(uta,  utb,  tta,  ttb):
    '''
    Greenwich apparent sidereal time (consistent with IAU 2000
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    lib.iauGst00a.argtypes =[c_double, c_double, c_double, c_double]
    
    lib.iauGst00a.restype  = c_double
    
    return  lib.iauGst00a(uta,  utb,  tta,  ttb)

def pymGst00b(uta,  utb):
    '''
    Greenwich apparent sidereal time (consistent with IAU 2000
    resolutions but using the truncated nutation model IAU 2000B).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    lib.iauGst00b.argtypes =[c_double, c_double]
    
    lib.iauGst00b.restype  = c_double
    
    return  lib.iauGst00b(uta,  utb)

def pymGst06(uta,  utb, tta,  ttb,  rnpb):
    '''
    Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date     
    ttb : float
        TT as a 2-part Julian Date    
    rnpb : float
        nutation x precession x bias matrix

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    lib.iauGst06.argtypes = [c_double,  c_double, c_double,  c_double,
                             vector_double3]
    lib.iauGst06.restype  = c_double
    
    
    return  lib.iauGst06(uta,  utb, tta,  ttb,  rnpb)

def pymGst06a(uta,  utb, tta,  ttb):
    '''
    Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date     
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    lib.iauGst06a.argtypes = [c_double,  c_double, c_double,  c_double]
                              
    lib.iauGst06a.restype  = c_double
    
    
    return  lib.iauGst06a(uta,  utb, tta,  ttb)

def pymGst94(uta, utb):
    '''
    Greenwich apparent sidereal time (consistent with IAU 1982/94
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    lib.iauGst94.argtypes =[c_double, c_double]
    
    lib.iauGst94.restype  = c_double
    
    return  lib.iauGst94(uta, utb)

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

#sofa vector_matrix model 

def pymA2tf(ndp, angle):
    '''
    Decompose radians into hours, minutes, seconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution    
    angle : float
        angle in radians

    Returns
    -------
    sign : bytes  
        '+' or '-'    
    ihmsf : tuple
        
    NDP         resolution
     :      ...0000 00 00
    -7         1000 00 00
    -6          100 00 00
    -5           10 00 00
    -4            1 00 00
    -3            0 10 00
    -2            0 01 00
    -1            0 00 10
     0            0 00 01
     1            0 00 00.1
     2            0 00 00.01
     3            0 00 00.001
     :            0 00 00.000...
    '''
     
    lib.iauA2tf.argtypes = [c_int, c_double, POINTER(c_char), c_int4] 
    lib.iauA2tf.restype  = None
    
    sign   = c_char()
    idmsf  = (c_int4)()
     
    lib.iauA2tf(ndp, float(angle), byref(sign), idmsf)
    return sign.value, tuple([x for x in idmsf])    


pym_af2a_msg = {
                  1: 'ideg outside range 0-359',
                  2: 'iamin outside range 0-59',
                  3: 'asec outside range 0-59.999...',
                  }    

def pymAf2a(s,   ideg,   iamin,   asec):
    '''
    Convert degrees, arcminutes, arcseconds to radians.

    Parameters
    ----------
    s : bytes
        sign:  '-' = negative, otherwise positive    
    ideg : int
        degrees    
    iamin : int
        arcminutes    
    asec : float
        arcseconds    

    Raises
    ------
    ValueError
        1: 'ideg outside range 0-359',
        2: 'iamin outside range 0-59',
        3: 'asec outside range 0-59.999...',
        
    Returns
    -------
    rad : float
        angle in radians

    '''
    lib.iauAf2a.argtypes = [c_char,  c_int,  c_int,  c_double,  c_double_p] 
    lib.iauAf2a.restype  = c_int
    
    rad = c_double()
    
    j   = lib.iauAf2a(s,   ideg,   iamin,   asec,  byref(rad))
    if j > 0:
        ws.warn(pym_af2a_msg[j], UserWarning, 2)
        #raise ValueError(pym_af2a_msg[j])
    
    return rad.value  

def pymAnpm(a):
    '''
    Normalize angle into the range -pi <= a < +pi.

    Parameters
    ----------
    a : float
        angle (radians)

    Returns
    -------
    function value : float
        angle in range +/-pi

    '''
    lib.iauAnpm.argtypes = [c_double]
    lib.iauAnpm.restype  = c_double
    
    return lib.iauAnpm(a)

def pymC2s(p):
    '''
    P-vector to spherical coordinates.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)

    '''
    lib.iauC2s.argtypes = [vector_double, c_double_p, c_double_p]
    lib.iauC2s.restype  =  None

    theta = c_double()
    phi   = c_double()

    lib.iauC2s(p, byref(theta),  byref(phi))
    
    return theta.value, phi.value


def pymCp(p):
    '''
    Copy a p-vector.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector to be copied

    Returns
    -------
    c : numpy.matrix(1,3)
        copy

    '''
    lib.iauCp.argtypes = [vector_double, vector_double]
    lib.iauCp.restype  =  None

    c = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))

    lib.iauCp(p, c)
    
    return c

def pymCp_A(p):

    return p 
 
def pymCpv(pv):
    '''
    Copy a position/velocity vector.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        position/velocity vector to be copied

    Returns
    -------
    c : numpy.matrix(2,3)
        copy

    '''
    lib.iauCpv.argtypes = [vector_double2, vector_double2]
    lib.iauCpv.restype  =  None

    c = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))

    lib.iauCpv(pv, c)
    
    return c

def pymCpv_A(pv):

#    c = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
#    c = pv
    return pv

def pymCr(r):
    '''
    Copy an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix to be copied

    Returns
    -------
    c : numpy.matrix(3,3)
        r-matrix to be copied

    '''
    lib.iauCr.argtypes = [vector_double3, vector_double3]
    lib.iauCr.restype  =  None

    c = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))

    lib.iauCr(r, c)
    
    return c

def pymCr_A(r):
    return r

#2023-05-28 
'''
void iauIr(double r[3][3])
'''
def pymIr():
    '''
    Initialize an r-matrix to the identity matrix.

    Returns
    -------
    r : numpy.matrix(3,3)
        r-matrix

    '''
    lib.iauIr.argtypes = [vector_double3]
    lib.iauIr.restype  =  None

    r = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))

    lib.iauIr(r)
    
    return r

def pymP2pv(p):
    '''
    Extend a p-vector to a pv-vector by appending a zero velocity.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    pv : numpy.matrix(2,3)
        pv-vector

    '''
    lib.iauP2pv.argtypes = [vector_double, vector_double2]
    lib.iauP2pv.restype  =  None

    pv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))

    lib.iauP2pv(p, pv)
    
    return pv


def pymP2s(p):
    '''
    P-vector to spherical polar coordinates.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance

    '''
    lib.iauP2s.argtypes = [vector_double, c_double_p, c_double_p, c_double_p]
    lib.iauP2s.restype  =  None

    theta = c_double()
    phi   = c_double()
    r     = c_double()

    lib.iauP2s(p, byref(theta), byref(phi), byref(r))
    
    return theta.value, phi.value, r.value 

def pymPap(a, b):
    '''
    Position-angle from two p-vectors.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        direction of reference point    
    b : numpy.matrix(1,3)
        direction of point whose PA is required

    Returns
    -------
    function value : float
        position angle of b with respect to a (radians)

    '''
    lib.iauPap.argtypes = [vector_double, vector_double]
    lib.iauPap.restype  =  c_double
    
    return lib.iauPap(a, b)

#def pymPas(al, ap, bl, bp):     Line No. 106

def pymPdp(a, b):
    '''
    p-vector inner (=scalar=dot) product.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    function value : float
        a . b

    '''
    lib.iauPdp.argtypes = [vector_double, vector_double]
    lib.iauPdp.restype  =  c_double
    
    return lib.iauPdp(a, b)

#Two  vectors dot to have scalar product,  using np.dot
def pymPdp_A(a, b):
 
    return  np.dot(a, b)

def pymPm(p):
    '''
    Modulus of p-vector.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    function value : float
        modulus

    '''
    lib.iauPm.argtypes = [vector_double]
    lib.iauPm.restype  =  c_double
    
    return lib.iauPm(p)

def pymPm_A(p):
        
    return np.linalg.norm(p, ord=2)

def pymPmp(a, b):
    '''
    P-vector subtraction.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    amb : numpy.matrix(1,3)
        a - b

    '''
    lib.iauPmp.argtypes = [vector_double, vector_double, vector_double]
    lib.iauPmp.restype  =  None
    amb = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPmp(a, b, amb)
    
    return  amb

def pymPmp_A(a, b):

    return  a - b

def pymPn(p):
    '''
    Convert a p-vector into modulus and unit vector.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    r : float
        modulus    
    u : numpy.matrix(1,3)
        unit vector

    '''
    lib.iauPn.argtypes = [vector_double, c_double_p, vector_double]
    lib.iauPn.restype  =  None
    
    r = c_double()
    u = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPn(p, byref(r), u)
    
    return r.value, u

#Python 
def pymPn_A(p):
    
    r = np.linalg.norm(p, ord=2) 
    u = p/r 
    
    return r, u

def pymPpp(a, b):
    '''
    P-vector addition.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    apb : numpy.matrix(1,3)
        a + b

    '''
    lib.iauPpp.argtypes = [vector_double, vector_double, vector_double]
    lib.iauPpp.restype  =  None
    apb = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPpp(a, b, apb)
    
    return  apb

def pymPpp_A(a, b):

    return  a + b

def pymPpsp(a, s, b):
    '''
    P-vector plus scaled p-vector.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    s : float     
        scalar (multiplier for b)
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    apsb : numpy.matrix(1,3)
        a + s*b

    '''
    lib.iauPpsp.argtypes = [vector_double, c_double, vector_double, vector_double]
    lib.iauPpsp.restype  =  None
    
    apsb = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPpsp(a, s, b, apsb)
    
    return  apsb

def pymPpsp_A(a, s, b):
    
    return a + s*b


def pymPv2p(pv):
    '''
    Discard velocity component of a pv-vector.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    p : numpy.matrix(1,3)
        p-vector

    '''
    lib.iauPv2p.argtypes = [vector_double2, vector_double]
    lib.iauPv2p.restype  =  None

    p  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))

    lib.iauPv2p(pv, p)
    
    return p 

def pymPv2s(pv):
    '''
    Convert position/velocity from Cartesian to spherical coordinates.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance    
    td : float
        rate of change of theta    
    pd : float
        rate of change of phi    
    rd : float
        rate of change of r

    '''
    lib.iauPv2s.argtypes = [vector_double2, c_double_p,  c_double_p, c_double_p,
                            c_double_p,  c_double_p,  c_double_p ]
    lib.iauPv2s.restype  =  None

 
    theta = c_double()
    phi   = c_double()
    r     = c_double()
    
    td    = c_double()
    pd    = c_double()
    rd    = c_double()
    
    lib.iauPv2s(pv, byref(theta), byref(phi), byref(r),
                    byref(td),    byref(pd),  byref(rd))
    
    return theta.value,  phi.value,  r.value,  td.value, pd.value, rd.value

def pymPvdpv(a, b):
    '''
    Inner (=scalar=dot) product of two pv-vectors.

    Parameters
    ----------
    a : numpy.matrix(2,3)
        first pv-vector   
    b : numpy.matrix(2,3)
        second pv-vector   

    Returns
    -------
    adb : numpy.matrix(1,2)
        DESCRIPTION.

    '''
    lib.iauPvdpv.argtypes = [vector_double2, vector_double2, array_double_2]
    lib.iauPvdpv.restype  =  None

    adb  = np.asmatrix(np.zeros(shape=(1,2), dtype=float, order='C'))

    lib.iauPvdpv(a, b, adb)
    
    return adb

def pymPvm(pv):
    '''
    Modulus of pv-vector.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    r : float
        modulus of position component
    s : float
        modulus of velocity component

    '''
    lib.iauPvm.argtypes = [vector_double2, c_double_p, c_double_p]
    lib.iauPvm.restype  =  None
    
    r = c_double()
    s = c_double()
   
    lib.iauPvm(pv, byref(r), byref(s))
    
    return r.value, s.value

def pymPvmpv(a, b):
    '''
    Subtract one pv-vector from another.

    Parameters
    ----------
    a : numpy.matrix(2,3)
        first pv-vector    
    b : numpy.matrix(2,3)
        second pv-vector

    Returns
    -------
    amb : numpy.matrix(2,3)
        a - b

    '''
    lib.iauPvmpv.argtypes = [vector_double2, vector_double2, vector_double2]
    lib.iauPvmpv.restype  =  None
    
    amb = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauPvmpv(a, b, amb)
    
    return amb

def pymPvmpv_A(a, b):
    
    return a - b

def pymPvppv(a, b):
    '''
    Add one pv-vector to another.

    Parameters
    ----------
    a : numpy.matrix(2,3) 
        first pv-vector    
    b : numpy.matrix(2,3) 
        second pv-vector

    Returns
    -------
    apb : numpy.matrix(2,3) 
        a + b

    '''
    lib.iauPvppv.argtypes = [vector_double2, vector_double2, vector_double2]
    lib.iauPvppv.restype  =  None
    
    apb = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauPvppv(a, b, apb)
    
    return apb

def pymPvppv_A(a, b):
    
    return a + b

#2023-05-30

def pymPvu(dt, pv):
    '''
    Update a pv-vector.

    Parameters
    ----------
    dt : float
        time interval    
    pv : numpy.matrix(2,3) 
        pv-vector

    Returns
    -------
    upv : numpy.matrix(2,3) 
        p updated, v unchanged

    '''
    lib.iauPvu.argtypes = [c_double, vector_double2, vector_double2]
    lib.iauPvu.restype  =  None
    
    upv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauPvu(dt, pv, upv)
    
    return upv

def pymPvup(dt, pv):
    '''
    Update a pv-vector, discarding the velocity component.

    Parameters
    ----------
    dt : float
        time interval    
    pv : numpy.matrix(2,3) 
        pv-vector

    Returns
    -------
    p : numpy.matrix(1,3) 
        p-vector

    '''
    lib.iauPvup.argtypes = [c_double, vector_double2, vector_double]
    lib.iauPvup.restype  =  None
    
    p = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPvup(dt, pv, p)
    
    return p

def pymPvxpv(a, b):
    '''
    Outer (=vector=cross) product of two pv-vectors.

    Parameters
    ----------
    a : numpy.matrix(2,3) 
        first pv-vector    
    b : numpy.matrix(2,3) 
        second pv-vector

    Returns
    -------
    axb : numpy.matrix(2,3) 
        a x b

    '''
    lib.iauPvxpv.argtypes = [vector_double2, vector_double2, vector_double2]
    lib.iauPvxpv.restype  =  None
    
    axb = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauPvxpv(a, b, axb)
    
    return axb


def pymPxp(a, b):
    '''
    p-vector outer (=vector=cross) product.

    Parameters
    ----------
    a : numpy.matrix(1,3) 
        first p-vector    
    b : numpy.matrix(1,3) 
        second p-vector

    Returns
    -------
    axb : numpy.matrix(1,3) 
        a x b

    '''
    lib.iauPxp.argtypes = [vector_double, vector_double, vector_double]
    lib.iauPxp.restype  =  None
    
    axb = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauPxp(a, b, axb)
    
    return axb

def pymPxp_A(a, b):
    
    return np.cross(a, b)

#2023-05-31
def pymRm2v(r):
    '''
    Express an r-matrix as an r-vector.

    Parameters
    ----------
    r : numpy.matrix(3,3)  
        rotation matrix

    Returns
    -------
    w : numpy.matrix(1,3)  
        rotation vector

    '''
    lib.iauRm2v.argtypes = [vector_double3,  vector_double]
    lib.iauRm2v.restype  =  None
    
    w = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauRm2v(r, w)
    
    return w


def pymRv2m(w):
    '''
    Form the r-matrix corresponding to a given r-vector.

    Parameters
    ----------
    w : numpy.matrix(1,3)     
        rotation vector     

    Returns
    -------
    r : numpy.matrix(3,3)     
        rotation matrix    

    '''
    lib.iauRv2m.argtypes = [vector_double, vector_double3]
    lib.iauRv2m.restype  =  None
    
    r = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauRv2m(w, r)
    
    return r


def pymRx(phi, r):
    '''
    Rotate an r-matrix about the x-axis.

    Parameters
    ----------
    phi : float    
        angle (radians)    
    r : numpy.matrix(3,3)    
        r-matrix    

    Returns
    -------
    r : numpy.matrix(3,3)    
        r-matrix, rotated

    '''
    lib.iauRx.argtypes = [c_double, vector_double3]
    lib.iauRx.restype  =  None
    
#   r = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauRx(phi, r)
    
    return r

#Python codes for iauRx
def pymRx_A(phi, r):
    
    cos_phi = np.math.cos(phi)
    sin_phi = np.math.sin(phi)

# to be in accordance with SOFA definition for Rx
    mat_rx  = np.array(
                  [ [1,         0,          0],
                    [0,    cos_phi,   sin_phi],
                    [0,   -sin_phi,   cos_phi]
                  ]).reshape(3,3)
    
    return np.dot(mat_rx, r)


#2023-06-03

def pymRxp(r, p):
    '''
    Multiply a p-vector by an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    p : numpy.matrix(1,3)
        p-vector    

    Returns
    -------
    rp : numpy.matrix(3,3)
        r * p

    '''
    lib.iauRxp.argtypes = [vector_double3, vector_double, vector_double]
    lib.iauRxp.restype  =  None
    
    rp = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    lib.iauRxp(r, p, rp)
    
    return rp 

def pymRxp_A(r, p):
    
    return np.dot(r, p) 


def pymRxpv(r, pv):
    '''
    Multiply a pv-vector by an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    rpv : numpy.matrix(2,3)
        r * pv

    '''
    lib.iauRxpv.argtypes = [vector_double3, vector_double2, vector_double2]
    lib.iauRxpv.restype  =  None
    
    rpv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    lib.iauRxpv(r, pv, rpv)
    
    return rpv 

def pymRxpv_A(r, pv):

#transpose 3,2 => 2, 3  
    
    return np.dot(r, pv).T 



def pymRxr(a, b):
    '''
    Multiply two r-matrices.

    Parameters
    ----------
    a : numpy.matrix(3,3)
        first r-matrix    
    b : numpy.matrix(3,3)
        second r-matrix

    Returns
    -------
    atb : numpy.matrix(3,3)
        a * b

    '''
    lib.iauRxr.argtypes = [vector_double3, vector_double3, vector_double3]
    lib.iauRxr.restype  =  None
    
    atb = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    lib.iauRxr(a, b, atb)
    
    return atb

def pymRxr_A(a, b):
    
    return np.dot(a, b)



def pymRy(theta, r):
    '''
    Rotate an r-matrix about the y-axis.

    Parameters
    ----------
    theta : float
        angle (radians)    
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    r : numpy.matrix(3,3)
        r-matrix, rotated

    '''
    lib.iauRy.argtypes = [c_double, vector_double3]
    lib.iauRy.restype  =  None
    
#   r = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauRy(theta, r)
    
    return r

#Python codes for iauRy
def pymRy_A(theta, r):
    
    cos_theta = np.math.cos(theta)
    sin_theta = np.math.sin(theta)

# to be in accordance with SOFA definition for Ry
    mat_ry  = np.array(
                  [ [cos_theta,    0,      -sin_theta],
                    [0,            1,              0 ],
                    [sin_theta,    0,       cos_theta]
                  ]).reshape(3,3)
    
    return np.dot(mat_ry, r)



def pymRz(psi, r):
    '''
    Rotate an r-matrix about the z-axis.

    Parameters
    ----------
    psi : float
        angle (radians)    
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    r : numpy.matrix(3,3)
        r-matrix, rotated

    '''
    lib.iauRz.argtypes = [c_double, vector_double3]
    lib.iauRz.restype  =  None
    
#   r = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauRz(psi, r)
    
    return r

#Python codes for iauRz
def pymRz_A(psi, r):
    
    cos_psi = np.math.cos(psi)
    sin_psi = np.math.sin(psi)

# to be in accordance with SOFA definition for Ry
    mat_rz  = np.array(
                  [ [ cos_psi,     sin_psi,     0],
                    [-sin_psi,     cos_psi,     0],
                    [0,            0,           1]
                  ]).reshape(3,3)
    
    return np.dot(mat_rz, r)

#def pymS2c(theta, phi): see previously


def pymS2p(theta, phi, r):
    '''
    Convert spherical polar coordinates to p-vector.

    Parameters
    ----------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance

    Returns
    -------
    p : numpy.matrix(1,3)
        Cartesian coordinates

    '''
    lib.iauS2p.argtypes = [c_double, c_double, c_double, vector_double]
    lib.iauS2p.restype  = None
    
    p = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    lib.iauS2p(theta, phi, r, p)
    
    return p

#2023-06-04

def pymS2pv(theta, phi, r, td, pd, rd):
    '''
    Convert position/velocity from spherical to Cartesian coordinates.

    Parameters
    ----------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance    
    td : float
        rate of change of theta    
    pd : float
        rate of change of phi    
    rd : float
        rate of change of r

    Returns
    -------
    pv : numpy.matrix(2,3)
        pv-vector

    '''
    lib.iauS2pv.argtypes = [c_double, c_double, c_double, 
                            c_double, c_double, c_double, 
                            vector_double2]
    lib.iauS2pv.restype  = None
    
    pv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    lib.iauS2pv(theta, phi, r, td, pd, rd, pv)
    
    return pv


def pymS2xpv(s1, s2, pv):
    '''
    Multiply a pv-vector by two scalars.

    Parameters
    ----------
    s1 : float
        scalar to multiply position component by    
    s2 : float
        scalar to multiply velocity component by    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    spv : numpy.matrix(2,3)
        pv-vector: p scaled by s1, v scaled by s2

    '''
    lib.iauS2xpv.argtypes = [c_double, c_double, vector_double2, vector_double2]
    lib.iauS2xpv.restype  = None
    
    spv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C')) 
    lib.iauS2xpv(s1, s2, pv, spv)
    
    return spv

def pymS2xpv_A(s1, s2, pv):
   
    spv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C')) 
    spv[0] = s1 * pv[0]
    spv[1] = s2 * pv[1]
    
    return   spv


def pymSepp(a, b):
    '''
    Angular separation between two p-vectors.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector (not necessarily unit length)    
    b : numpy.matrix(1,3)
        second p-vector (not necessarily unit length)

    Returns
    -------
    function value : float
        angular separation (radians, always positive)

    '''
    lib.iauSepp.argtypes = [vector_double, vector_double]
    lib.iauSepp.restype  = c_double
    
    return lib.iauSepp(a, b)

#Python codes for iauSepp
def pymSepp_A(a, b):
     
    axb     = np.cross(a, b)
    sin_phi = np.linalg.norm(axb, ord=2)
    cos_phi = np.dot(a, b)
        
    return  np.math.atan2(sin_phi, cos_phi)


def pymSeps(al,   ap,   bl,   bp):
    '''
    Angular separation between two sets of spherical coordinates.

    Parameters
    ----------
    al : float
        first longitude (radians)    
    ap : float
        first latitude (radians)    
    bl : float
        second longitude (radians)    
    bp : float
        second latitude (radians)

    Returns
    -------
    function value : float
        angular separation (radians)

    '''
    lib.iauSeps.argtypes = [c_double, c_double, c_double, c_double]
    lib.iauSeps.restype  = c_double
    
    return lib.iauSeps(al,   ap,   bl,   bp)


def pymSxp(s, p):
    '''
    Multiply a p-vector by a scalar.

    Parameters
    ----------
    s : float
        scalar    
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    sp : numpy.matrix(1,3)
        s * p

    '''
    lib.iauSxp.argtypes = [c_double,  vector_double, vector_double]
    lib.iauSxp.restype  = None
    
    sp = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C')) 
    lib.iauSxp(s, p, sp)
    
    return sp 

 
def pymSxp_A(s, p):

    return   np.asmatrix(s*p)

  
def pymSxpv(s, pv):
    '''
    Multiply a pv-vector by a scalar.

    Parameters
    ----------
    s : float
        scalar    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    spv : numpy.matrix(2,3)
        s * pv

    '''
    lib.iauSxpv.argtypes = [c_double,  vector_double2, vector_double2]
    lib.iauSxpv.restype  = None
    
    spv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C')) 
    lib.iauSxpv(s, pv, spv)
    
    return spv 

pym_tf2a_msg = {
                1:'ihour outside range 0-23',
                2:'imin outside range 0-59',
                3:'sec outside range 0-59.999...'
                }   
 
def pymTf2a(s, ihour, imin, sec): 
    '''
    Convert hours, minutes, seconds to radians.

    Parameters
    ----------
    s : bytes  
        sign:  '-' = negative, otherwise positive      
    ihour : int    
        hours    
    imin : int    
        minutes    
    sec : float    
        seconds    
        
    Raises
    ------
    ValueError
        1:'ihour outside range 0-23',
        2:'imin outside range 0-59',
        3:'sec outside range 0-59.999...'

    Returns
    -------
    rad : float
        angle in radians

    '''
    #Note: for single byte s(b'+'), here using c_char that can work.  
    
    lib.iauTf2a.argtypes = [c_char, c_int, c_int, c_double, c_double_p]
    lib.iauTf2a.restype  = c_int
    
    rad = c_double()
     
    j = lib.iauTf2a(s, ihour, imin, sec, byref(rad))
    if j > 0:
        ws.warn(pym_tf2a_msg[j], UserWarning, 2)
        
    return rad.value 

def pymTf2a_A(s, ihour, imin, sec): 
    
#Note: for single byte s('+'), here using c_wchar that can work.  
    
    lib.iauTf2a.argtypes = [c_wchar, c_int, c_int, c_double, c_double_p]
    lib.iauTf2a.restype  = c_int
    
    rad = c_double()
     
    j = lib.iauTf2a(s, ihour, imin, sec, byref(rad))
    
    if j > 0:
        ws.warn(pym_tf2a_msg[j], UserWarning, 2)
        
    return rad.value 

#def pymTf2d(s, ihour, imin, sec): see previously

def pymTr(r):
    '''
    Transpose an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    rt : numpy.matrix(3,3)
        transpose

    '''
    lib.iauTr.argtypes = [vector_double3, vector_double3]
    lib.iauTr.restype  =  None
    
    rt = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauTr(r, rt)
    
    return rt

def pymTr_A(r):
    
    return r.T
 
def pymTrxp(r, p):
    '''
    Multiply a p-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    trp : numpy.matrix(1,3)
        r^T * p

    '''
    lib.iauTrxp.argtypes = [vector_double3, vector_double, vector_double]
    lib.iauTrxp.restype  =  None
    
    trp = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauTrxp(r, p, trp)
    
    return trp 

def pymTrxp_A(r, p):
   
#   trp  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    trp  = np.dot(r.T, p)
    
    return np.asmatrix(trp)

def pymTrxpv(r, pv):
    '''
    Multiply a pv-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    trpv : numpy.matrix(2,3)
        r^T * pv

    '''
    lib.iauTrxpv.argtypes = [vector_double3, vector_double2, vector_double2]
    lib.iauTrxpv.restype  =  None
    
    trpv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauTrxpv(r, pv, trpv)
    
    return trpv


def pymTrxpv_A(r, pv):
    
    return np.dot(r.T, pv.T).T 


def pymZp(p):
    '''
    Zero a p-vector

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    p : numpy.matrix(1,3)
        zero p-vector
    '''
    
    lib.iauZp.argtypes = [vector_double]
    lib.iauZp.restype  =  None
    
#   pr = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    
    lib.iauZp(p)
    
    return p

def pymZp_A(p):
   
    pt  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
     
    return pt

def pymZpv(pv):
    '''
    Zero a pv-vector

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector
        
    Returns
    -------
    pv : numpy.matrix(2,3)
        zero pv-vector

    '''
    lib.iauZpv.argtypes = [vector_double2]
    lib.iauZpv.restype  =  None
    
#   pr = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
    
    lib.iauZpv(pv)
    
    return pv

def pymZpv_A(pv):
   
    pvt  = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
     
    return pvt

def pymZr(r):
    '''
    Initialize an r-matrix to the null matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    r : numpy.matrix(3,3)
        null matrix

    '''
    lib.iauZr.argtypes = [vector_double3]
    lib.iauZr.restype  =  None
    
#   pr = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
    
    lib.iauZr(r)
    
    return r


def pymZr_A(r):
   
    rt  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
     
    return rt


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
    
