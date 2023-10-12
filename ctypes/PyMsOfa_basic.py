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


#sofa vector_matrix model 


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
 
