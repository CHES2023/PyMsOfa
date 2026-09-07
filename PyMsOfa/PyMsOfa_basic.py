"""
Created on Sat Aug  9  04:15:10 2023
Done    on Mon Aug  26 08:25:12 2025
@author: Dr. Jianghui JI  (jijh@pmo.ac.cn)
Description: SOFA basic tools
"""

import numpy as np
try:                                  # installed as part of the PyMsOfa package
    from .sofa_const import *
except ImportError:                   # flat layout: all modules in one folder
    from sofa_const import *

#2025-08-25
 
def pymS2c(theta, phi):
    """
    Convert spherical coordinates to Cartesian unit vector(s).

    Parameters
    ----------
    theta : float or array_like
        Longitude angle(s) in radians.
    phi : float or array_like
        Latitude angle(s) in radians.

    Returns
    -------
    ndarray
        Cartesian unit vector(s), shape (..., 3)
    """
    theta = np.asarray(theta, dtype=float)
    phi   = np.asarray(phi,   dtype=float)

    cos_phi = np.cos(phi)

    return np.stack((
        cos_phi * np.cos(theta),
        cos_phi * np.sin(theta),
                  np.sin(phi)
    ), axis=-1)


#2025-08-25
def pymAnp(a):
    """
    Normalize angle(s) into the range 0 <= a < 2*pi.

    Parameters
    ----------
    a : float or array_like
        Angle(s) in radians.

    Returns
    -------
    float or ndarray
        Normalized angle(s) in [0, 2*pi).
        Scalar in -> scalar out; array in -> array out.
    """
    w = np.mod(a, D2PI)
    
    return w.item() if np.ndim(w) == 0 else w



def pymAnpm(a):
    """
    Normalize angle(s) into the range -pi <= a < +pi.

    Parameters
    ----------
    a : float or array_like
        Angle(s) in radians.

    Returns
    -------
    float or ndarray
        Normalized angle(s) in [-pi, +pi).
        Scalar in -> scalar out; array in -> array out.
    """
    w = np.mod(a + DPI, D2PI) - DPI
    
    return w.item() if np.ndim(w) == 0 else w


#2025-08-26  
def pymC2s(p):
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    p : array_like, shape (3,) or (..., 3)
        Cartesian vector(s).

    Returns
    -------
    theta : float or ndarray
        Longitude(s) in radians, range [-pi, +pi].
    phi : float or ndarray
        Latitude(s) in radians, range [-pi/2, +pi/2].
    """
    p = np.asarray(p, dtype=float)

    # reshape to (..., 3)
    if p.shape[-1] != 3:
        raise ValueError("Input must have last dimension of size 3")

    x = p[..., 0]
    y = p[..., 1]
    z = p[..., 2]

    rxy = np.hypot(x, y)

    # SOFA convention: atan2(0,0) -> 0
    theta = np.arctan2(y, x)
    theta = np.where(rxy == 0.0, 0.0, theta)

    # Latitude; for zero vector, return 0
    phi = np.arctan2(z, rxy)
    phi = np.where((rxy == 0.0) & (z == 0.0), 0.0, phi)

    # scalar in -> scalar out
    if np.ndim(theta) == 0:
        return theta.item(), phi.item()

    return theta, phi


def pymCp(p):
    """
    Copy a p-vector (3D vector).

    Parameters
    ----------
    p : array-like, shape (3,) or (..., 3)
        Input p-vector(s).

    Returns
    -------
    c : np.ndarray
        Copy of the p-vector(s), same shape as input.
    """
    return np.array(p, dtype=float, copy=True)


def pymCpv(pv):
    """
    Copy a position/velocity vector (2×3).

    Parameters
    ----------
    pv : array-like, shape (2,3) or (...,2,3)
        Position/velocity vector(s).

    Returns
    -------
    c : np.ndarray
        Copy of pv, same shape as input.
    """
    return np.array(pv, dtype=float, copy=True)



def pymCr(r):
    """
    Copy a 3x3 matrix  

    Parameters
    ----------
    r : array-like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    c : np.ndarray
        Deep copy of r, shape (3,3)
    """
    return np.array(r, dtype=float, copy=True)


def pymIr():
    """
    Initialize a 3x3 matrix to the identity matrix.

    Returns
    -------
    r : np.ndarray, shape (3,3)
        Identity matrix
    """
    return np.eye(3, dtype=float)


def pymP2pv(p):
    """
    Extend a 3D position vector to a position/velocity vector by appending zero velocity.

    Parameters
    ----------
    p : array-like, shape (3,)
        Input position vector.

    Returns
    -------
    pv : np.ndarray, shape (2,3)
        Position/velocity vector, velocity initialized to zero.
    """
    return np.array([p, [0.0, 0.0, 0.0]], dtype=float)


def pymP2s(p):
    """
    Convert a 3D Cartesian vector or multiple vectors to spherical polar coordinates.

    Parameters
    ----------
    p : array-like, shape (3,) or (N,3)
        Input Cartesian vector(s).

    Returns
    -------
    theta : float or np.ndarray
        Longitude angle(s) in radians.
    phi : float or np.ndarray
        Latitude angle(s) in radians.
    r : float or np.ndarray
        Radial distance(s).
    """
    p = np.asarray(p, dtype=float)
    if p.ndim == 1:
        p = p.reshape(1, 3)   
        scalar_input = True
    else:
        scalar_input = False

    x, y, z = p[:,0], p[:,1], p[:,2]

    # Radial distance
    r = np.linalg.norm(p, axis=1)

    # Longitude angle (theta), 0 if x=y=0
    theta = np.arctan2(y, x)
    theta[(x==0) & (y==0)] = 0.0

    # Latitude angle (phi), 0 if r=0
    rho = np.hypot(x, y)
    phi = np.arctan2(z, rho)
    phi[r==0] = 0.0

    if scalar_input:
        return float(theta[0]), float(phi[0]), float(r[0])
    
    return theta, phi, r



def pymPap(a, b):
    """
    Compute the position angle of vector b with respect to vector a.

    Parameters
    ----------
    a : array-like, shape (3,)
        Direction of reference point
    b : array-like, shape (3,)
        Direction of point whose PA is required

    Returns
    -------
    pa : float
        Position angle of b with respect to a (radians, -pi to +pi)
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    am = np.linalg.norm(a)
    bm = np.linalg.norm(b)

    # Handle null vectors
    if am == 0.0 or bm == 0.0:
        return 0.0

    # Unit vector of a
    au = a / am

    # North axis tangential from a
    xa, ya, za = a
    eta = np.array([-xa*za, -ya*za, xa**2 + ya**2])

    # East axis tangential from a
    xi = np.cross(eta, au)

    # Vector from a to b
    a2b = b - a

    # Resolve into components along north and east
    st = np.dot(a2b, xi)
    ct = np.dot(a2b, eta)

    # Handle degenerate case
    if st == 0.0 and ct == 0.0:
        ct = 1.0

    pa = np.arctan2(st, ct)
    
    return pa


def pymPdp(a, b):
    """
    Compute the dot (scalar) product of two 3D vectors.

    Parameters
    ----------
    a : array-like, shape (3,)
        First vector
    b : array-like, shape (3,)
        Second vector

    Returns
    -------
    w : float
        Dot product a . b
    """
    #a = np.asarray(a, dtype=float)
    #b = np.asarray(b, dtype=float)
    result = np.dot(a, b)
     
    return float(result) if result.ndim == 0 else result



def pymPm(p):
    """
    Calculate the modulus (Euclidean norm) of a 3D vector using NumPy.
    
    Parameters
    ----------
    p : array-like of 3 elements
        The 3D vector [x, y, z]
    
    Returns
    -------
    float
        The modulus (length) of the vector
    
    """
    # np.linalg.norm directly computes the Euclidean norm of the vector
    
    return np.linalg.norm(p)


def pymPmp(a, b):
    """
    Subtract two 3D vectors (p-vectors) using NumPy: amb = a - b

    Parameters
    ----------
    a, b : array-like of shape (3,)
        Input vectors

    Returns
    -------
    np.ndarray
        Difference vector a - b

    Example
    -------
    >>> pymPmp(np.array([1,2,3]), np.array([0.5,1.5,1]))
    array([0.5, 0.5, 2. ])
    """
    a = np.asarray(a)
    b = np.asarray(b)
    
    return a - b



def pymPn(p):
    """
    Convert a 3D vector into modulus and unit vector using NumPy.

    Parameters
    ----------
    p : array-like of 3 floats
        Input vector

    Returns
    -------
    r : float
        Modulus of the vector
    u : ndarray of 3 floats
        Unit vector

    Notes
    -----
    Return order follows the SOFA ``iauPn`` convention: ``(r, u)``
    (modulus first, unit vector second).
    """
    p = np.asarray(p, dtype=float)

    r = float(np.linalg.norm(p))

    if r == 0.0:
        return 0.0, np.zeros(3, dtype=float)

    u = p / r

    return r, u
    


def pymPpp(a, b):
    """
    P-vector addition: compute a + b for 3D vectors (or batches of vectors).

    Parameters
    ----------
    a : array-like, shape (3,) or (N,3)
        First vector(s).
    b : array-like, shape (3,) or (N,3)
        Second vector(s).

    Returns
    -------
    apb : np.ndarray
        Result of a + b, shape (3,) if input was (3,), otherwise (N,3).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    
    return a + b


def pymPpsp(a, s, b):
    """
    P-vector plus scaled p-vector: a + s*b  
    
    Parameters
    ----------
    a : array-like of shape (3,)
        First vector
    s : float
        Scalar multiplier
    b : array-like of shape (3,)
        Second vector
    
    Returns
    -------
    ndarray of shape (3,)
        Resulting vector
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    
    return a + s * b


def pymPv2p(pv):
    """
    Discard velocity component of a pv-vector (NumPy version).

    Parameters
    ----------
    pv : array-like, shape (2,3)
        pv[0] = position vector
        pv[1] = velocity vector

    Returns
    -------
    ndarray, shape (3,)
        Position vector only (independent copy of pv[0])
    """
    pv = np.asarray(pv, dtype=float)
    
    # Return an independent copy of pv[0]
    return np.array(pv[0], dtype=float, copy=True)


def pymPv2s(pv):
    """
    Convert position+velocity vector to spherical coordinates.

    Parameters
    ----------
    pv : ndarray, shape (2,3)
        Cartesian position and velocity vectors:
        pv[0] = [x, y, z]   (position)
        pv[1] = [xd, yd, zd] (velocity)

    Returns
    -------
    theta : float
        Longitude angle (radians), atan2(y, x)
    phi : float
        Latitude angle (radians), atan2(z, sqrt(x^2+y^2))
    r : float
        Radial distance (norm of position vector)
    td : float
        Time derivative of longitude (radians/unit time)
    pd : float
        Time derivative of latitude (radians/unit time)
    rd : float
        Radial velocity (same units as velocity components)
    """
    # Ensure input is a float array
    pv = np.array(pv, dtype=float, copy=True)
    pos, vel = pv
    x, y, z = pos
    xd, yd, zd = vel

    # Squared distance in XY-plane
    rxy2 = x*x + y*y
    # Squared distance in 3D
    r2 = rxy2 + z*z
    # True radial distance
    rtrue = np.sqrt(r2)

    # Handle zero position vector: use velocity as position
    if rtrue == 0.0:
        x, y, z = vel
        rxy2 = x*x + y*y
        r2 = rxy2 + z*z
        rtrue = np.sqrt(r2)

    # Distance in XY-plane
    rxy = np.sqrt(rxy2)
    # Dot product of position and velocity in XY-plane
    xyp = x*xd + y*yd

    if rxy2 != 0.0:
        # Longitude angle theta = atan2(y, x)
        theta = np.arctan2(y, x)
        # Latitude angle phi = atan2(z, sqrt(x^2+y^2))
        phi   = np.arctan2(z, rxy)
        # Longitude rate
        td    = (x*yd - y*xd) / rxy2
        # Latitude rate
        pd    = (zd*rxy2 - z*xyp) / (r2*rxy)
    else:
        # Near the pole, theta is set to zero
        theta = 0.0
        # phi only depends on z when rxy=0
        phi   = np.arctan2(z, rxy) if z != 0.0 else 0.0
        td    = 0.0
        pd    = 0.0

    # Radial distance
    r  = rtrue
    # Radial velocity (projection along position vector)
    rd = (xyp + z*zd) / rtrue if rtrue != 0.0 else 0.0

    return theta, phi, r, td, pd, rd


def pymPvdpv(a, b):
    """
    Inner (=scalar=dot) product of two pv-vectors.

    Parameters
    ----------
    a : ndarray, shape (2,3)
        First pv-vector (a[0]=position, a[1]=velocity)
    b : ndarray, shape (2,3)
        Second pv-vector (b[0]=position, b[1]=velocity)

    Returns
    -------
    adb : ndarray, shape (2,)
        Dot-product of pv-vectors:
        adb[0] = a_position . b_position
        adb[1] = a_position . b_velocity + a_velocity . b_position
    """
    a = np.array(a, dtype=float, copy=True)
    b = np.array(b, dtype=float, copy=True)

    # Constant (position) part
    adb0 = np.dot(a[0], b[0])

    # Velocity part
    adb1 = np.dot(a[0], b[1]) + np.dot(a[1], b[0])

    return np.array([adb0, adb1], dtype=float)


def pymPvm(pv):
    """
    Compute modulus of position and velocity components of a pv-vector
    
    Parameters
    ----------
    pv : array-like, shape (2,3)
        pv[0] = position vector
        pv[1] = velocity vector
    
    Returns
    -------
    r : float
        Modulus of position vector
    s : float
        Modulus of velocity vector
    """
    pv= np.asarray(pv, dtype=float)
    r = np.linalg.norm(pv[0])
    s = np.linalg.norm(pv[1])
    
    return r, s

def pymPvmpv(a, b):
    """
    Subtract two pv-vectors: amb = a - b (NumPy version)
    
    Parameters
    ----------
    a, b : array-like, shape (2,3)
        pv-vectors
    
    Returns
    -------
    amb : ndarray, shape (2,3)
        Resulting pv-vector
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    # amb = a - b
    return a - b

def pymPvppv(a, b):
    """
    Add two pv-vectors: apb = a + b  
    
    Parameters
    ----------
    a, b : array-like, shape (2,3)
        pv-vectors
    
    Returns
    -------
    apb : ndarray, shape (2,3)
        Resulting pv-vector
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    # apb = a + b
    return a + b

def pymPvu(dt, pv):
    """
    Update a pv-vector: upv = pv at time dt later (NumPy version)

    Parameters
    ----------
    dt : float
        Time interval
    pv : array-like, shape (2,3)
        Input pv-vector [position, velocity]

    Returns
    -------
    upv : ndarray, shape (2,3)
        Updated pv-vector
    """
    pv  = np.asarray(pv, dtype=float)
    pos = pv[0] + dt * pv[1]  
    vel = np.array(pv[1], dtype=float, copy=True)   
    
    return np.array([pos, vel])


def pymPvup(dt, pv):
    """
    Update a pv-vector, discarding the velocity component (NumPy version).
    
    Parameters
    ----------
    dt : float
        Time interval
    pv : array-like, shape (2,3)
        pv-vector [position, velocity]
    
    Returns
    -------
    p : ndarray, shape (3,)
        Updated position vector
    """
    pv = np.asarray(pv, dtype=float)
    # p = position + dt * velocity
    
    return pv[0] + dt * pv[1]


def pymPvxpv(a, b):
    """
    Outer (cross) product of two pv-vectors (NumPy version).

    Parameters
    ----------
    a, b : array-like, shape (2,3)
        Input pv-vectors

    Returns
    -------
    axb : ndarray, shape (2,3)
        Cross product pv-vector
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    
    # Position part
    pos = np.cross(a[0], b[0])
    
    # Velocity part: ap × bv + av × bp
    vel = np.cross(a[0], b[1]) + np.cross(a[1], b[0])
    
    return np.array([pos, vel])

 
def pymPxp(a, b):
    """
    Cross product of two 3D vectors  

    Parameters
    ----------
    a, b : array-like, shape (3,)
        Input vectors

    Returns
    -------
    ndarray, shape (3,)
        Cross product a x b
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    
    return np.cross(a, b)
 


#2025-08-27
def pymRm2v(r):
    """
    Express a rotation matrix as a rotation vector (axis * angle).

    Parameters
    ----------
    r : ndarray (3,3)
        Rotation matrix.

    Returns
    -------
    w : ndarray (3,)
        Rotation vector: direction = Euler axis, magnitude = rotation angle [rad].
    """
    r = np.asarray(r, dtype=float)

    # Extract antisymmetric part
    x, y, z = r[1,2] - r[2,1], r[2,0] - r[0,2], r[0,1] - r[1,0]
    s2 = np.linalg.norm([x, y, z])

    if s2 > 0.0:
        c2  = np.trace(r) - 1.0
        phi = np.arctan2(s2, c2)  # angle
        w = (phi / s2) * np.array([x, y, z])
    else:
        w = np.zeros(3)

    return w


def pymRv2m(w):
    """
    Convert a rotation vector to a rotation matrix.

    Parameters
    ----------
    w : array_like, shape (3,)
        Rotation vector. Its direction gives the Euler axis, 
        and its magnitude (radians) gives the rotation angle.

    Returns
    -------
    r : ndarray, shape (3, 3)
        Corresponding rotation matrix.

    Notes
    -----
    - If the input vector is [0, 0, 0], the identity matrix is returned.
    - Uses Rodrigues' rotation formula.
    """
    w = np.array(w, dtype=float, copy=True)
    phi = np.linalg.norm(w)

    # Identity if no rotation
    if phi == 0.0:
        return np.eye(3)

    # Normalize axis
    u = w / phi
    x, y, z = u

    c = np.cos(phi)
    s = np.sin(phi)
    f = 1.0 - c

    # Rodrigues' formula
    r = np.array([
        [c + x*x*f,   x*y*f + z*s, x*z*f - y*s],
        [y*x*f - z*s, c + y*y*f,   y*z*f + x*s],
        [z*x*f + y*s, z*y*f - x*s,   c + z*z*f]
    ])
    
    return r



def pymRxp(r, p):
    """
    Multiply a 3x3 rotation matrix with a 3D vector (rp = r * p).

    Parameters
    ----------
    r : array_like, shape (3, 3)
        Rotation matrix.
    p : array_like, shape (3,)
        Input vector.

    Returns
    -------
    rp : ndarray, shape (3,)
        Rotated vector (r * p).

    """
    r = np.array(r, dtype=float, copy=True).reshape(3, 3)
    p = np.array(p, dtype=float, copy=True).reshape(3)

    # Matrix-vector multiplication
    rp = r @ p
    
    return rp

 
def pymRxpv(r, pv):
    """
    Multiply a pv-vector by an r-matrix.
    
    Parameters
    ----------
    r : array_like, shape (3,3)
        Rotation matrix.
    pv : array_like, shape (2,3)
        Position-velocity vector (first row: position, second row: velocity).
    
    Returns
    -------
    rpv : ndarray, shape (2,3)
        Result of r * pv.
    
    """
    r   = np.array(r,  dtype=float, copy=True)
    pv  = np.array(pv, dtype=float, copy=True)
    
    rpv = np.dot(pv, r.T)  # Matrix-vector multiplication for both rows
    
    return rpv   


def pymRxr(a, b):
    """
    Multiply two 3x3 rotation matrices (r-matrices) using NumPy.
    
    Parameters
    ----------
    a : array_like, shape (3,3)
        First rotation matrix.
    b : array_like, shape (3,3)
        Second rotation matrix.
    
    Returns
    -------
    atb : ndarray, shape (3,3)
        The result of matrix multiplication a * b.
    """
    a = np.asarray(a, dtype=float).reshape(3, 3)
    b = np.asarray(b, dtype=float).reshape(3, 3)
    
    atb = np.dot(a, b)
    
    return atb


def pymRx(phi, r=None):
    """
    Rotate an r-matrix about the x-axis.

    Parameters
    ----------
    phi : float
        Rotation angle in radians.
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    r_out : ndarray, shape (3,3)
        Rotated matrix (original is not modified).

    Notes
    -----
    - Rotation is anticlockwise as seen looking towards the origin 
      from the positive x-axis.
    - Equivalent to left-multiplying r by the Rx(phi) matrix:

        [ 1      0      0 ]
        [ 0   cosφ   sinφ ]
        [ 0  -sinφ   cosφ ]
    """
    #r = np.array(r, dtype=float, copy=True)

    c = np.cos(phi)
    s = np.sin(phi)

    Rx = np.array([
        [1, 0, 0],
        [0, c, s],
        [0,-s, c]
    ])

    
    return Rx if r is None else Rx @ r


def pymRy(theta, r=None):
    """
    Rotate an r-matrix about the y-axis.

    Parameters
    ----------
    theta : float
        Rotation angle in radians.
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    r_out : ndarray, shape (3,3)
        Rotated matrix (original is not modified).

    Notes
    -----
    - Rotation is anticlockwise as seen looking towards the origin 
      from the positive y-axis.
    - Equivalent to left-multiplying r by the Ry(theta) matrix:

        [ cosθ   0  -sinθ ]
        [  0     1     0  ]
        [ sinθ   0   cosθ ]
    """
    #r = np.array(r, dtype=float, copy=True)

    c = np.cos(theta)
    s = np.sin(theta)

    Ry = np.array([
        [ c, 0, -s],
        [ 0, 1,  0],
        [ s, 0,  c]
    ])

    return Ry if r is None else Ry @ r  


def pymRz(psi, r=None):
    """
    Rotate an r-matrix about the z-axis.

    Parameters
    ----------
    psi : float
        Rotation angle in radians.
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    r_out : ndarray, shape (3,3)
        Rotated matrix (original is not modified).

    Notes
    -----
    - Rotation is anticlockwise as seen looking towards the origin 
      from the positive z-axis.
    - Equivalent to left-multiplying r by the Rz(psi) matrix:

        [ cosψ   sinψ   0 ]
        [-sinψ   cosψ   0 ]
        [  0      0     1 ]
    """
    #r = np.array(r, dtype=float, copy=True)

    c = np.cos(psi)
    s = np.sin(psi)

    Rz = np.array([
        [ c,  s, 0],
        [-s,  c, 0],
        [ 0,  0, 1]
    ])

    return Rz if r is None else Rz @ r  
 
'''
def pymRx(phi, r):
    """
    Rotate an r-matrix about the x-axis.

    Parameters
    ----------
    phi : float
        Rotation angle in radians.
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    r_out : ndarray, shape (3,3)
        Rotated matrix (original is not modified).

    Notes
    -----
    - Rotation is anticlockwise as seen looking towards the origin 
      from the positive x-axis.
    - Equivalent to left-multiplying r by the Rx(phi) matrix:

        [ 1      0      0 ]
        [ 0   cosφ   sinφ ]
        [ 0  -sinφ   cosφ ]
    """
    r = np.array(r, dtype=float, copy=True)

    c = np.cos(phi)
    s = np.sin(phi)

    Rx = np.array([
        [1, 0, 0],
        [0, c, s],
        [0,-s, c]
    ])

    return Rx @ r    

def pymRy(theta, r):
    """
    Rotate an r-matrix about the y-axis.

    Parameters
    ----------
    theta : float
        Rotation angle in radians.
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    r_out : ndarray, shape (3,3)
        Rotated matrix (original is not modified).

    Notes
    -----
    - Rotation is anticlockwise as seen looking towards the origin 
      from the positive y-axis.
    - Equivalent to left-multiplying r by the Ry(theta) matrix:

        [ cosθ   0  -sinθ ]
        [  0     1     0  ]
        [ sinθ   0   cosθ ]
    """
    r = np.array(r, dtype=float, copy=True)

    c = np.cos(theta)
    s = np.sin(theta)

    Ry = np.array([
        [ c, 0, -s],
        [ 0, 1,  0],
        [ s, 0,  c]
    ])

    return Ry @ r   


def pymRz(psi, r):
    """
    Rotate an r-matrix about the z-axis.

    Parameters
    ----------
    psi : float
        Rotation angle in radians.
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    r_out : ndarray, shape (3,3)
        Rotated matrix (original is not modified).

    Notes
    -----
    - Rotation is anticlockwise as seen looking towards the origin 
      from the positive z-axis.
    - Equivalent to left-multiplying r by the Rz(psi) matrix:

        [ cosψ   sinψ   0 ]
        [-sinψ   cosψ   0 ]
        [  0      0     1 ]
    """
    r = np.array(r, dtype=float, copy=True)

    c = np.cos(psi)
    s = np.sin(psi)

    Rz = np.array([
        [ c,  s, 0],
        [-s,  c, 0],
        [ 0,  0, 1]
    ])

    return Rz @ r   
'''

def pymS2p(theta, phi, r):
    """
    Convert spherical polar coordinates to a Cartesian p-vector.

    Parameters
    ----------
    theta : float
        Longitude angle in radians.
    phi : float
        Latitude angle in radians.
    r : float
        Radial distance.

    Returns
    -------
    p : ndarray, shape (3,)
        Cartesian coordinates corresponding to the spherical polar coordinates.

    Notes
    -----
    - p = r * [cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)]
    - Uses numpy for efficient computation.
    """
    cos_phi = np.cos(phi)
    
    p = r * np.array([
        cos_phi * np.cos(theta),
        cos_phi * np.sin(theta),
        np.sin(phi)
    ])
    
    return p


def pymS2pv(theta, phi, r, td, pd, rd):
    """
    Convert position and velocity from spherical to Cartesian coordinates.

    Parameters
    ----------
    theta : float
        Longitude angle in radians.
    phi : float
        Latitude angle in radians.
    r : float
        Radial distance.
    td : float
        Rate of change of theta (radians per unit time).
    pd : float
        Rate of change of phi (radians per unit time).
    rd : float
        Radial velocity (rate of change of r).

    Returns
    -------
    pv : ndarray, shape (2,3)
        Cartesian position and velocity vector.
        pv[0] = position [x, y, z]
        pv[1] = velocity [vx, vy, vz]

    Notes
    -----
    - Uses standard spherical to Cartesian conversion:
        x = r * cos(phi) * cos(theta)
        y = r * cos(phi) * sin(theta)
        z = r * sin(phi)
    - Velocity transformation accounts for theta_dot, phi_dot, and r_dot.
    """
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)
    rcp = r * cp

    # Position components
    x = rcp * ct
    y = rcp * st
    z = r * sp

    # Intermediate terms for velocity
    rpd = r * pd
    w = rpd*sp - cp*rd

    # Velocity components
    vx = -y*td - w*ct
    vy =  x*td - w*st
    vz = rpd*cp + sp*rd

    pv = np.array([[x, y, z],
                   [vx, vy, vz]], dtype=float)
    
    return pv

 

def pymS2xpv(s1, s2, pv):
    """
    Multiply a pv-vector by two scalars.

    Parameters
    ----------
    s1 : float
        Scalar to multiply the position component by.
    s2 : float
        Scalar to multiply the velocity component by.
    pv : array_like, shape (2,3)
        Input pv-vector [position; velocity].

    Returns
    -------
    spv : ndarray, shape (2,3)
        Scaled pv-vector: position scaled by s1, velocity scaled by s2.
    """
    pv = np.array(pv, dtype=float, copy=True)
    spv = np.empty_like(pv)
    spv[0] = s1 * pv[0]  # scale position
    spv[1] = s2 * pv[1]  # scale velocity
    
    return spv


def pymSepp(a, b):
    """
    Angular separation between two p-vectors.

    Parameters
    ----------
    a : array_like, shape (3,)
        First p-vector (not necessarily unit length).
    b : array_like, shape (3,)
        Second p-vector (not necessarily unit length).

    Returns
    -------
    s : float
        Angular separation in radians (always positive).

    Notes
    -----
    Uses cross and dot products for full accuracy for small and large angles.
    Returns 0 if either vector is null.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if np.all(a == 0) or np.all(b == 0):
        return 0.0

    # Sine of angle (magnitude of cross product)
    ss = np.linalg.norm(np.cross(a, b))
    # Cosine of angle (dot product)
    cs = np.dot(a, b)

    # Angle
    return np.arctan2(ss, cs)


def pymSeps(al, ap, bl, bp):
    """
    Angular separation between two points in spherical coordinates
    using pure NumPy.

    Parameters
    ----------
    al, ap : float
        Longitude and latitude of first point (radians)
    bl, bp : float
        Longitude and latitude of second point (radians)

    Returns
    -------
    s : float
        Angular separation in radians
    """
    # Unit vector
    a = np.array([np.cos(ap) * np.cos(al),
                  np.cos(ap) * np.sin(al),
                  np.sin(ap)])
    b = np.array([np.cos(bp) * np.cos(bl),
                  np.cos(bp) * np.sin(bl),
                  np.sin(bp)])
    
    # Magnitude of the cross product
    ss = np.linalg.norm(np.cross(a, b))
    
    # Dot product
    cs = np.dot(a, b)
    
    # arctan2 ensures numerical accuracy for both small and large angles
    return np.arctan2(ss, cs)



def pymSxp(s, p):
    """
    Multiply a p-vector by a scalar.

    Parameters
    ----------
    s : float
        Scalar multiplier
    p : array_like, shape (3,)
        Input p-vector

    Returns
    -------
    sp : ndarray, shape (3,)
        Resulting scaled p-vector
    """
    p  = np.array(p, dtype=float)
    sp = s * p
    
    return sp


def pymSxpv(s, pv):
    """
    Multiply a pv-vector by a scalar.

    Parameters
    ----------
    s : float
        Scalar multiplier
    pv : array_like, shape (2,3)
        Input pv-vector (position and velocity)

    Returns
    -------
    spv : ndarray, shape (2,3)
        Scaled pv-vector
    """
    pv  = np.array(pv, dtype=float)
    spv = s * pv
    
    return spv


def pymTr(r):
    """
    Transpose a 3x3 rotation matrix.

    Parameters
    ----------
    r : array_like, shape (3,3)
        Input rotation matrix.

    Returns
    -------
    rt : ndarray, shape (3,3)
        Transposed matrix.
    """
    r = np.array(r, dtype=float)
    
    return r.T


def pymTrxp(r, p):
    """
    Multiply a p-vector by the transpose of a 3x3 rotation matrix.

    Parameters
    ----------
    r : array_like, shape (3,3)
        Rotation matrix.
    p : array_like, shape (3,)
        Position vector.

    Returns
    -------
    trp : ndarray, shape (3,)
        Result of r^T * p.
    """
    r = np.array(r, dtype=float)
    p = np.array(p, dtype=float)
    
    trp = r.T @ p
    
    return trp



def pymTrxpv(r, pv):
    """
    Multiply a pv-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : ndarray, shape (3,3)
        Rotation matrix.
    pv : ndarray, shape (2,3)
        Position-velocity vector.

    Returns
    -------
    trpv : ndarray, shape (2,3)
        Result of r^T * pv.
    """

    r  = np.array(r,  dtype=float, copy=True)
    pv = np.array(pv, dtype=float, copy=True)

    # Transpose the matrix
    rt = r.T
    # Apply separately to position and velocity
    trpv = np.array([rt @ pv[0], rt @ pv[1]])
    
    return trpv


def pymZp():
    """
    Return a zero p-vector (3-element vector).

    Returns
    -------
    p : ndarray, shape (3,)
        Zero vector [0.0, 0.0, 0.0].
    """
    return np.zeros(3, dtype=float)


def pymZpv():
    """
    Return a zero pv-vector (2x3 array).

    Returns
    -------
    pv : ndarray, shape (2,3)
        Zero pv-vector.
    """
    return np.zeros((2, 3), dtype=float)


def pymZr():
    """
    Return a zero r-matrix (3x3 array).

    Returns
    -------
    r : ndarray, shape (3,3)
        Zero r-matrix.
    """
    return np.zeros((3, 3), dtype=float)
