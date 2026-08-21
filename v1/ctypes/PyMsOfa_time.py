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














