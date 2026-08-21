# PyMsOfa 2.x Introduction

A **pure Python / NumPy** implementation of the IAU Standards of Fundamental
Astronomy (SOFA) service.  Version **2.0.0** provides all **247** SOFA routines
(naming `iauXxx` → `pymXxx`) and is based on SOFA release 2023-10-11.

> **What changed from 1.x** — one implementation instead of three (ctypes / cffi /
> pure python + a compiled C library); `import PyMsOfa as sf` exposes everything at
> the top level; invalid input raises `ValueError` (no trailing status codes, no
> `-1e9`/`None` sentinels); routines accept NumPy arrays as well as scalars;
> `pymASTROM` / `pymLDBODY` are ordinary Python classes.

```bash
pip install PyMsOfa
```

### Quick start

```python
import PyMsOfa as sf

sf.pymCal2jd(2024, 2, 27)                 # -> (2400000.5, 60367.0)
sf.pymD2tf(3, 0.5)                        # -> ('+', [12, 0, 0, 0])
sf.pymS2c([0.0, 0.7854], [0.0, 0.7854])   # shape (2, 3)  (array input)
sf.DAS2R                                   # 4.84813681109536e-06
```

### Citation

> Ji, Jiang-Hui, Tan, Dong-jie, Bao, Chun-hui, Huang, Xiu-min, Hu, Shoucun,
> Dong, Yao, Wang, Su. 2023, *PyMsOfa: A Python Package for the Standards of
> Fundamental Astronomy (SOFA) Service*, Research in Astronomy and Astrophysics,
> 23, 125015, doi:10.1088/1674-4527/ad0499

---

## 1. Basic tools

### 1.1 Copy a parameter

`pymCp(p)` : Copy a p-vector.(1*3)

`pymCpv(pv)` : Copy a position/velocity vector.(2*3)

`pymCr(r)` : Copy an r-matrix.(3*3)

```python
import PyMsOfa as sf

A = [3,4,5]
B = [[3,4,5],[4,5,6]]
C = [[3,4,5],[4,5,6],[5,6,7]]

print(sf.pymCp(A))
print(sf.pymCpv(B))
print(sf.pymCr(C))
```

```text
[3. 4. 5.]
[[3. 4. 5.]
 [4. 5. 6.]]
[[3. 4. 5.]
 [4. 5. 6.]
 [5. 6. 7.]]
```

### 1.2 Initialize / zero a parameter

`pymZp()` : Zero a p-vector.(1*3)

`pymZpv()` : Zero a pv-vector.(2*3)

`pymZr()` : Initialize an r-matrix to the null matrix.(3*3)

`pymIr()` : Initialize an r-matrix to the identity matrix.(3*3)

```python
import PyMsOfa as sf

A = [3,4,5]
B = [[3,4,5],[4,5,6]]
C = [[3,4,5],[4,5,6],[5,6,7]]

print(sf.pymZp())
print(sf.pymZpv())
print(sf.pymZr())
print(sf.pymIr())
```

```text
[0. 0. 0.]
[[0. 0. 0.]
 [0. 0. 0.]]
[[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]
[[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]
```

### 1.3 Normalize an angle

`pymAnp(a)` : Normalize angle into the range 0 <= a < 2pi.

`pymAnpm(a)` : Normalize angle into the range -pi <= a < +pi.

```python
import PyMsOfa as sf
import math as ma

A = 7/2*ma.pi

print(sf.pymAnp(A))
print(sf.pymAnpm(A))
```

```text
4.71238898038469
-1.5707963267948966
```

### 1.4 Parameter processing

`pymPm(p)` : Modulus of p-vector.

`pymPn(p)` : Convert a p-vector into modulus and unit vector.

`pymTr(r)` : Transpose an r-matrix.

`pymPv2p(pv)` : Discard velocity component of a pv-vector.

`pymPvm(pv)` : Modulus of pv-vector.

```python
import PyMsOfa as sf

A = [3,4,5]
B = [[3,4,5],[4,5,6]]
C = [[3,4,5],[3,4,5],[3,4,5]]

print(sf.pymPm(A))
print(sf.pymPn(A))
print(sf.pymTr(C))
print(sf.pymPv2p(B))
print(sf.pymPvm(B))
```

```text
7.0710678118654755
(7.0710678118654755, array([0.42426407, 0.56568542, 0.70710678]))
[[3. 3. 3.]
 [4. 4. 4.]
 [5. 5. 5.]]
[3. 4. 5.]
(7.0710678118654755, 8.774964387392123)
```

### 1.5 Parameter operations

`pymPpp(a, b)` : P-vector addition.

`pymPmp(a, b)` : P-vector subtraction.

`pymSxp(s, p)` : Multiply a p-vector by a scalar.

`pymPpsp(a, s, b)` : P-vector plus scaled p-vector.

`pymPdp(a, b)` : p-vector inner (=scalar=dot) product.

`pymPxp(a, b)` : p-vector outer (=vector=cross) product.

`pymRxp(r, p)` : Multiply a p-vector by an r-matrix.

`pymTrxp(r, p)` : Multiply a p-vector by the transpose of an r-matrix.

`pymRxpv(r, pv)` : Multiply a pv-vector by an r-matrix.

`pymTrxpv(r, pv)` : Multiply a pv-vector by the transpose of an r-matrix.

`pymRxr(a, b)` : Multiply two r-matrices.

`pymPvppv(a, b)` : Add one pv-vector to another.

`pymPvmpv(a, b)` : Subtract one pv-vector from another.

`pymPvdpv(a, b)` : Inner (=scalar=dot) product of two pv-vectors.

`pymPvxpv(a, b)` : Outer (=vector=cross) product of two pv-vectors.

`pymSxpv(s, pv)` : Multiply a pv-vector by a scalar.

`pymS2xpv(s1, s2, pv)` : Multiply a pv-vector by two scalars.

```python
import PyMsOfa as sf

A = [3,4,5]
B = [4,5,6]

print(sf.pymPpp(A,B))
print(sf.pymPmp(A,B))
print(sf.pymSxp(2,A))
print(sf.pymPpsp(A,2,B))
print(sf.pymPdp(A,B))
print(sf.pymPxp(A,B))
```

```text
[ 7.  9. 11.]
[-1 -1 -1]
[ 6.  8. 10.]
[11. 14. 17.]
62.0
[-1.  2. -1.]
```

```python
import PyMsOfa as sf

A = [3,4,5]
B = [[3,4,5],[4,5,6]]
C = [[3,4,5],[3,4,5],[3,4,5]]

print(sf.pymRxp(C,A))
print(sf.pymTrxp(C,A))
print(sf.pymRxpv(C,B))
print(sf.pymTrxpv(C,B))
print(sf.pymRxr(C,C))
```

```text
[50. 50. 50.]
[36. 48. 60.]
[[50. 50. 50.]
 [62. 62. 62.]]
[[36. 48. 60.]
 [45. 60. 75.]]
[[36. 48. 60.]
 [36. 48. 60.]
 [36. 48. 60.]]
```

```python
import PyMsOfa as sf

A = [[3,4,5],[3,4,5]]
B = [[3,4,5],[4,5,6]]

print(sf.pymPvppv(A,B))
print(sf.pymPvmpv(A,B))
print(sf.pymPvdpv(A,B))
print(sf.pymPvxpv(A,B))
print(sf.pymSxpv(2,A))
print(sf.pymS2xpv(2,3,A))
```

```text
[[ 6.  8. 10.]
 [ 7.  9. 11.]]
[[ 0.  0.  0.]
 [-1. -1. -1.]]
[ 50. 112.]
[[ 0.  0.  0.]
 [-1.  2. -1.]]
[[ 6.  8. 10.]
 [ 6.  8. 10.]]
[[ 6.  8. 10.]
 [ 9. 12. 15.]]
```

### 1.6 Matrix rotation about an axis

`pymRx(phi, r=None)` : Rotate an r-matrix about the x-axis.

`pymRy(theta, r=None)` : Rotate an r-matrix about the y-axis.

`pymRz(psi, r=None)` : Rotate an r-matrix about the z-axis.

```python
import PyMsOfa as sf
import math as ma

PHI = 1/4*ma.pi
A = [[1,1,1],[2,2,2],[2,2,2]]
B = [[2,2,2],[1,1,1],[2,2,2]]
C = [[2,2,2],[2,2,2],[1,1,1]]

print(sf.pymRx(PHI,A))
print(sf.pymRy(PHI,B))
print(sf.pymRz(PHI,C))
```

```text
[[1.         1.         1.        ]
 [2.82842712 2.82842712 2.82842712]
 [0.         0.         0.        ]]
[[0.         0.         0.        ]
 [1.         1.         1.        ]
 [2.82842712 2.82842712 2.82842712]]
[[2.82842712 2.82842712 2.82842712]
 [0.         0.         0.        ]
 [1.         1.         1.        ]]
```

### 1.7 Coordinate conversion

`pymC2s(p)` : P-vector to spherical coordinates.

`pymS2c(theta, phi)` : Convert spherical coordinates to Cartesian.

`pymP2s(p)` : P-vector to spherical polar coordinates (including radial distance).

`pymS2p(theta, phi, r)` : Convert spherical polar coordinates (including radial distance) to p-vector.

`pymPv2s(pv)` : Convert position/velocity from Cartesian to spherical coordinates (including latitude and longitude, radial distance, and the rate of change of all three).

`pymS2pv(theta, phi, r, td, pd, rd)` : Convert position/velocity from spherical to Cartesian coordinates.

`pymPvstar(pv)` : Convert star position+velocity vector to catalog coordinates (including RA, DEC, pmRA, pmDEC, plx and RV).

`pymStarpv(ra, dec, pmr, pmd, px, rv)` : Convert star catalog coordinates to position+velocity vector.

```python
import PyMsOfa as sf
import math as ma

A = [0.5, 0.5, ma.sqrt(2)/2]
THETA = ma.pi/4
PHI = ma.pi/4

print(sf.pymC2s(A))
print(sf.pymS2c(THETA, PHI))
```

```text
(0.7853981633974483, 0.7853981633974483)
[0.5        0.5        0.70710678]
```

```python
import PyMsOfa as sf
import math as ma

A = [0.5, 0.5, ma.sqrt(2)/2]
THETA = ma.pi/4
PHI = ma.pi/4
R = 1.0

print(sf.pymP2s(A))
print(sf.pymS2p(THETA, PHI, R))
```

```text
(0.7853981633974483, 0.7853981633974483, 1.0)
[0.5        0.5        0.70710678]
```

```python
import PyMsOfa as sf
import math as ma

PV = [[0.5, 0.5, ma.sqrt(2)/2],
      [1, 1, 1]]
THETA = ma.pi/4        #longitude angle (radians)
PHI = ma.pi/4        #latitude angle (radians)
R = 1.0            #radial distance
TD = 0.0        #rate of change of theta
PD = -0.293        #rate of change of phi
RD = 1.707        #rate of change of r

print(sf.pymPv2s(PV))
print(sf.pymS2pv(THETA, PHI, R, TD, PD, RD))
```

```text
(0.7853981633974483, 0.7853981633974483, 1.0, 0.0, -0.29289321881345254, 1.7071067811865475)
[[0.5        0.5        0.70710678]
 [1.         1.         0.99984899]]
```

```python
import PyMsOfa as sf
import math as ma

PV = [[0.5, 0.5, ma.sqrt(2)/2],
      [1, 1, ma.sqrt(2)]]
RA = ma.pi/4
DEC = ma.pi/4
PMRA = 0.0        #radians/year
PMDEC = -105.935    #radians/year
PLX = 206264.806    #arcsce
RV = 2941.778        #km/s,positive = receding

print(sf.pymPvstar(PV))
print(sf.pymStarpv(RA, DEC, PMRA, PMDEC ,PLX, RV))
```

```text
(0.7853981633974483, 0.7853981633974483, 0.0, -5.734762705342487e-14, 206264.80624709636, 3443.1425958340537)
[[0.5        0.5        0.70710678]
 [1.00000031 1.00000031 0.99999962]]
```

### 1.8 Rotation vector ↔ rotation matrix

`pymRv2m(w)` : Form the r-matrix corresponding to a given r-vector.

`pymRm2v(r)` : Express an r-matrix as an r-vector.

```python
import PyMsOfa as sf
import math as ma

W = [1,2,3]
R = [[0.00, -0.80, -0.60],
     [0.80, -0.36, 0.48],
     [0.60, 0.48, -0.64]]

print(sf.pymRv2m(W))
print(sf.pymRm2v(R))
```

```text
[[-0.69492056 -0.19200697  0.69297817]
 [ 0.71352099 -0.30378504  0.6313497 ]
 [ 0.08929286  0.93319235  0.34810748]]
[ 0.          1.41371669 -1.88495559]
```

### 1.9 Position relationships

`pymPas(al, ap, bl, bp)` : Position-angle from spherical coordinates.

`pymPap(a, b)` : Position-angle from two p-vectors.

`pymSeps(al, ap, bl, bp)` : Angular separation between two sets of spherical coordinates.

`pymSepp(a, b)` : Angular separation between two p-vectors.

```python
import PyMsOfa as sf
import math as ma

AL = 0.5 * ma.pi        #longitude of point A (e.g. RA) in radians 
AP = 0.8 * ma.pi        #latitude of point A (e.g. Dec) in radians  
BL = 0.4 * ma.pi        #longitude of point B  
BP = 0.6 * ma.pi        #latitude of point B 
A = [1, 2, 3]
B = [2, 3, 4]

print(sf.pymPas(AL, AP, BL, BP))
print(sf.pymPap(A, B))
print(sf.pymSeps(AL, AP, BL, BP))
print(sf.pymSepp(A, B))
```

```text
2.982899019958473
-2.389552564397458
0.6488468770597016
0.12186756768575521
```

### 1.10 Projection (tangent plane)

`pymTpors(xi, eta, a, b)` : In the tangent plane projection, given the rectangular coordinates of a star and its spherical coordinates, determine the spherical coordinates of the tangent point.

`pymTporv(xi, eta, v)` : In the tangent plane projection, given the rectangular coordinates of a star and its direction cosines, determine the direction cosines of the tangent point.

`pymTpsts(xi, eta, a0, b0)` : In the tangent plane projection, given the star's rectangular coordinates and the spherical coordinates of the tangent point, solve for the spherical coordinates of the star.

`pymTpstv(xi, eta, v0)` : In the tangent plane projection, given the star's rectangular coordinates and the direction cosines of the tangent point, solve for the direction cosines of the star.

`pymTpxes(a, b, a0, b0)` : In the tangent plane projection, given celestial spherical coordinates for a star and the tangent point, solve for the star's rectangular coordinates in the tangent plane.

`pymTpxev(v, v0)` : In the tangent plane projection, given celestial direction cosines for a star and the tangent point, solve for the star's rectangular coordinates in the tangent plane.

```python
import PyMsOfa as sf
import math as ma

XI = 1                    #rectangular coordinates of star image  
ETA = 2                    #rectangular coordinates of star image  
A = ma.pi/4                #star's spherical coordinates  
B = ma.pi/6                #star's spherical coordinates  
V = [0.2, 0.4, ma.sqrt(0.8)]        #star's direction cosines
A0 = ma.pi/4+ma.pi/180            #tangent point's spherical coordinates  
B0 = ma.pi/6-ma.pi/180            #tangent point's spherical coordinates  
V0 = [0.2, 0.401, ma.sqrt(0.7992)]    #tangent point's direction cosines

print(sf.pymTpors(XI, ETA, A, B))
print(sf.pymTporv(XI, ETA, V))
print(sf.pymTpsts(XI, ETA, A0, B0))
print(sf.pymTpstv(XI, ETA, V0))
print(sf.pymTpxes(A, B, A0, B0))
print(sf.pymTpxev(V, V0))
```

```text
(0.29451548510813685, -0.5275089774303864, 4.417873495276552, 1.4548041954319986)
(array([ 0.96490136, -0.04162585,  0.25929261]), array([ 0.49884201, -0.60859618,  0.61706348]))
(2.4683630662278397, 1.1482475839728823)
[-0.6094643 -0.3072788  0.7308446]
(-0.015118273995868148, 0.017521042566124613)
(-0.0004463207759768221, 0.001000448828478158)
```

### 1.11 Other functions

`pymP2pv(p)` : Extend a p-vector to a pv-vector by appending a zero velocity.

`pymPvu(dt, pv)` : Update a pv-vector.

`pymPvup(dt, pv)` : Update a pv-vector, discarding the velocity component.

```python
import PyMsOfa as sf

P = [1, 2, 3]
DT = 10                #time interval  
PV = [[1, 2, 3],
      [4, 5, 6]]

print(sf.pymP2pv(P))
print(sf.pymPvu(DT, PV))
print(sf.pymPvup(DT, PV))
```

```text
[[1. 2. 3.]
 [0. 0. 0.]]
[[41. 52. 63.]
 [ 4.  5.  6.]]
[41. 52. 63.]
```

## 2. Time

### 2.1 Sidereal time

`pymGmst82(dj1, dj2)` : Universal Time to Greenwich mean sidereal time (IAU 1982 model).

`pymGst94(uta, utb)` : Greenwich apparent sidereal time (consistent with IAU 1982/94 resolutions).

`pymGmst00(uta, utb, tta, ttb)` : Greenwich mean sidereal time (model consistent with IAU 2000 resolutions).

`pymGst00a(uta, utb, tta, ttb)` : Greenwich apparent sidereal time (consistent with IAU 2000 resolutions).

`pymGst00b(uta, utb)` : Greenwich apparent sidereal time (consistent with IAU 2000 resolutions but using the truncated nutation model IAU 2000B).

`pymGmst06(uta, utb, tta, ttb)` : Greenwich mean sidereal time (consistent with IAU 2006 precession).

`pymGst06(uta, utb, tta, ttb, rnpb)` : Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

`pymGst06a(uta, utb, tta, ttb)` : Greenwich apparent sidereal time (consistent with IAU 2000 and 2006 resolutions).

```python
import PyMsOfa as sf

DJ1 = 2400000.5
DJ2 = 52000.0
UTA = 2400000.5
UTB = 52200.0
TTA = 2400000.5
TTB = 52200.0
RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4], 
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
    #nutation x precession x bias matrix

print(sf.pymGmst82(DJ1, DJ2))
print(sf.pymGst94(UTA, UTB))
print(sf.pymGmst00(UTA, UTB, TTA, TTB))
print(sf.pymGst00a(UTA, UTB, TTA, TTB))
print(sf.pymGst00b(UTA, UTB))
print(sf.pymGmst06(UTA, UTB, TTA, TTB))
print(sf.pymGst06(UTA, UTB, TTA, TTB, RNPB))
print(sf.pymGst06a(UTA, UTB, TTA, TTB))
```

```text
3.3060549204240584
0.4633455111601113
0.4634280218387905
0.46334550139240926
0.4633455023818729
0.46342802175551745
0.4643595506102784
0.4633455012681201
```

### 2.2 Time-unit conversion

`pymD2tf(ndp, days)` : Decompose days to hours, minutes, seconds, fraction.

`pymTf2d(sign, ihour, imin, sec)` : Convert hours, minutes, seconds to days.

`pymA2af(ndp, angle)` : Decompose radians into degrees, arcminutes, arcseconds, fraction.
`pymA2tf(ndp, angle)` : Decompose radians into hours, minutes, seconds, fraction.


`pymAf2a(sign, ideg, iamin, asec)` : Convert degrees, arcminutes, arcseconds to radians.

`pymTf2a(sign, ihour, imin, sec)` : Convert hours, minutes, seconds to radians.

```python
import PyMsOfa as sf
import math as ma

NDP = 3        #resolution
DAYS = 0.5
S='+'
IHOUR = 5
IMIN = 24
SEC = 50.3
ANGLE = ma.pi/4
IDEG = 274
IAMIN = 35
ASEC = 42.38

print(sf.pymD2tf(NDP, DAYS))
print(sf.pymTf2d(S, IHOUR, IMIN, SEC))
print(sf.pymA2af(NDP, ANGLE))
print(sf.pymA2tf(NDP, ANGLE))
print(sf.pymAf2a(S, IDEG, IAMIN, ASEC))
print(sf.pymTf2a(S, IHOUR, IMIN, SEC))
```

```text
('+', [12, 0, 0, 0])
0.22558217592592592
('+', [45, 0, 0, 0])
('+', [3, 0, 0, 0])
4.792588701805817
1.4173746133393783
```

### 2.3 Date conversion

`pymCal2jd(iy, im, id)` : Gregorian Calendar to Julian Date.

`pymJd2cal(dj1, dj2)` : Julian Date to Gregorian year, month, day, and fraction of a day.

`pymJdcalf(ndp, dj1, dj2)` : Julian Date to Gregorian Calendar, expressed in a form convenient for formatting messages:  rounded to a specified precision.

`pymEpj(dj1, dj2)` : Julian Date to Julian Epoch

`pymEpj2jd(epj)` : Julian Epoch to Julian Date

`pymEpb(dj1, dj2)` : Julian Date to Besselian Epoch.

`pymEpb2jd(epb)` : Besselian Epoch to Julian Date

```python
import PyMsOfa as sf

IY = 2024
IM = 2
ID = 27
DJ1 = 2400000.5
DJ2 = 60367.0
NDP = 3                #resolution
EPJ = 2024.154688569473        #Julian Epoch
EPB = 2024.1564820038504    #Besselian Epoch

print(sf.pymCal2jd(IY, IM, ID))
print(sf.pymJd2cal(DJ1, DJ2))
print(sf.pymJdcalf(NDP, DJ1, DJ2))
print(sf.pymEpj(DJ1, DJ2))
print(sf.pymEpj2jd(EPJ))
print(sf.pymEpb(DJ1, DJ2))
print(sf.pymEpb2jd(EPB))
```

```text
(2400000.5, 60367.0)
(2024, 2, 27, 0.0)
[2024, 2, 27, 0]
2024.154688569473
(2400000.5, 60366.99999999998)
2024.1564820038504
(2400000.5, 60367.0)
```

### 2.4 Time-scale conversion

`pymDat(iy, im, id, fd)` : For a given UTC date, calculate Delta(AT) = TAI-UTC.

`pymDtf2d(scale, iy, im, id, ihr, imn, sec)` : Encode date and time fields into 2-part Julian Date (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).

`pymDtdb(date1, date2, ut, elong, u, v)` : An approximation to TDB-TT, the difference between barycentric dynamical time and terrestrial time, for an observer on the Earth.

`pymD2dtf(scale, ndp, d1, d2)` : Format for output a 2-part Julian Date (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).

`pymUtctai(utc1, utc2)` : Time scale transformation:  Coordinated Universal Time, UTC, to International Atomic Time, TAI.

`pymTaiutc(tai1, tai2)` : Time scale transformation:  International Atomic Time, TAI, to Coordinated Universal Time, UTC.

`pymTaiut1(tai1, tai2, dta)` : Time scale transformation:  International Atomic Time, TAI, to Universal Time, UT1.

`pymUt1tai(ut11, ut12, dta)` : Time scale transformation:  Universal Time, UT1, to International Atomic Time, TAI.

`pymUtcut1(utc1, utc2, dut1)` : Time scale transformation:  Coordinated Universal Time, UTC, to Universal Time, UT1.

`pymUt1utc(ut11, ut12, dut1)` : Time scale transformation:  Universal Time, UT1, to Coordinated Universal Time, UTC.

`pymTaitt(tai1, tai2)` : Time scale transformation:  International Atomic Time, TAI, to Terrestrial Time, TT.

`pymTttai(tt1, tt2)` : Time scale transformation:  Terrestrial Time, TT, to International Atomic Time, TAI.

`pymTttcg(tt1, tt2)` : Time scale transformation:  Terrestrial Time, TT, to Geocentric Coordinate Time, TCG.

`pymTcgtt(tcg1, tcg2)` : Time scale transformation:  Geocentric Coordinate Time, TCG, to Terrestrial Time, TT.

`pymTttdb(tt1, tt2, dtr)` : Time scale transformation:  Terrestrial Time, TT, to Barycentric Dynamical Time, TDB.

`pymTdbtt(tdb1, tdb2, dtr)` : Time scale transformation:  Barycentric Dynamical Time, TDB, to Terrestrial Time, TT.

`pymTtut1(tt1, tt2, dt)` : Time scale transformation:  Terrestrial Time, TT, to Universal Time, UT1.

`pymUt1tt(ut11, ut12, dt)` : Time scale transformation:  Universal Time, UT1, to Terrestrial Time, TT.

`pymTdbtcb(tdb1, tdb2)` : Time scale transformation:  Barycentric Dynamical Time, TDB, to Barycentric Coordinate Time, TCB.

`pymTcbtdb(tcb1, tcb2)` : Time scale transformation:  Barycentric Coordinate Time, TCB, to Barycentric Dynamical Time, TDB.

```python
import PyMsOfa as sf

NDP = 3
D1, D2 = 2400000.5, 60367.0
UTC1, UTC2 = 2400000.5, 60367.0
TAI1, TAI2 = 2400000.5, 60367.0
DTA = -32.6659
UT11, UT12 = 2400000.5, 60367.0
DUT1 = 0.3341
TT1, TT2 = 2400000.5, 60367.0
TCG1, TCG2 = 2400000.5, 60367.0
DTR = -0.000201
TDB1, TDB2 = 2400000.5, 60367.0
DT = 64.8499
TCB1, TCB2 = 2400000.5, 60367.0

print(sf.pymD2dtf("UTC", NDP, D1, D2))
print(sf.pymUtctai(UTC1, UTC2))
print(sf.pymTaiutc(TAI1, TAI2))
print(sf.pymTaiut1(TAI1, TAI2, DTA))
print(sf.pymUt1tai(UT11, UT12, DTA))
print(sf.pymUtcut1(UTC1, UTC2, DUT1))
print(sf.pymUt1utc(UT11, UT12, DUT1))
print(sf.pymTaitt(TAI1, TAI2))
print(sf.pymTttai(TT1, TT2))
print(sf.pymTttcg(TT1, TT2))
print(sf.pymTcgtt(TCG1, TCG2))
print(sf.pymTttdb(TT1, TT2, DTR))
print(sf.pymTdbtt(TDB1, TDB2, DTR))
print(sf.pymTtut1(TT1, TT2, DT))
print(sf.pymUt1tt(UT11, UT12, DT))
print(sf.pymTdbtcb(TDB1, TDB2))
print(sf.pymTcbtdb(TCB1, TCB2))
```

```text
(2024, 2, 27, [0, 0, 0, 0])
(2400000.5, 60367.00042824074)
(2400000.5, 60366.99957175926)
(2400000.5, 60366.999621922456)
(2400000.5, 60367.000378077544)
(2400000.5, 60367.000003866895)
(array(2400000.5), array(60366.99999613))
(2400000.5, 60367.0003725)
(2400000.5, 60366.9996275)
(2400000.5, 60367.000012003205)
(2400000.5, 60366.999987996795)
(2400000.5, 60366.99999999767)
(2400000.5, 60367.00000000233)
(2400000.5, 60366.99924942246)
(2400000.5, 60367.00075057754)
(2400000.5, 60367.00026704677)
(2400000.5, 60366.99973295323)
```

## 3. Precession, nutation and polar motion

### 3.1 Solar-system body parameters

`pymFapa03(t)` : Fundamental argument, IERS Conventions (2003) : general accumulated precession in longitude.

`pymFae03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Earth.

`pymFad03(t)` : Fundamental argument, IERS Conventions (2003) : mean elongation of the Moon from the Sun.

`pymFaom03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of the Moon's ascending node.

`pymFaf03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of the Moon minus mean longitude of the ascending node.

`pymFal03(t)` : Fundamental argument, IERS Conventions (2003) : mean anomaly of the Moon.

`pymFalp03(t)` : Fundamental argument, IERS Conventions (2003) : mean anomaly of the Sun.

`pymFame03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Mercury.

`pymFave03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Venus.

`pymFama03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Mars.

`pymFaju03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Jupiter.

`pymFasa03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Saturn.

`pymFaur03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Uranus.

`pymFane03(t)` : Fundamental argument, IERS Conventions (2003) : mean longitude of Neptune.

`pymPlan94(date1, date2, np_plt)` : Approximate heliocentric position and velocity of a nominated major planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or Neptune (but not the Earth itself).

`pymMoon98(date1, date2)` : Approximate geocentric position and velocity of the Moon.

```python
import PyMsOfa as sf

T = 0.8        #TDB, Julian centuries since J2000.0
DATE1 = 2400000.5
DATE2 = 60367.0
NP = 1        #planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars, 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)

print(sf.pymFapa03(T))
print(sf.pymFae03(T))
print(sf.pymFad03(T))
print(sf.pymFaom03(T))
print(sf.pymFaf03(T))
print(sf.pymFal03(T))
print(sf.pymFalp03(T))
print(sf.pymFame03(T))
print(sf.pymFave03(T))
print(sf.pymFama03(T))
print(sf.pymFaju03(T))
print(sf.pymFasa03(T))
print(sf.pymFaur03(T))
print(sf.pymFane03(T))
print(sf.pymPlan94(DATE1, DATE2, NP))
print(sf.pymMoon98(DATE1, DATE2))
```

```text
0.019508847622400002
1.7447137389131058
1.9467092053972077
-5.973618440951301
0.2597711366751447
5.132369751109105
6.226797973505531
5.417338184297677
3.4249004605338342
3.275506840277835
5.275711665202486
5.371574539440829
5.180636450180414
2.079343830860413
[[ 0.34673534 -0.14077236 -0.11113816]
 [ 0.00707931  0.02365909  0.01190497]]
[[-2.70343937e-03 -2.07161048e-04 -3.71702717e-05]
 [ 4.78785113e-05 -4.91093096e-04 -2.67211978e-04]]
```

### 3.2 Frame bias, precession, nutation and polar motion

`pymNumat(epsa, dpsi, deps)` : Form the matrix of nutation.

`pymLtp(epj)` : Long-term precession matrix.

`pymLtpb(epj)` : Long-term precession matrix, including ICRS frame bias.

`pymLtpequ(epj)` : Long-term precession of the equator.

`pymLtpecl(epj)` : Long-term precession of the ecliptic.

`pymPrec76(date01, date02, date11, date12)` : IAU 1976 precession model. This function forms the three Euler angles which implement general precession between two dates, using the IAU 1976 model (as for the FK5 catalog).

`pymPmat76(date1, date2)` : Precession matrix from J2000.0 to a specified date, IAU 1976 model.

`pymPnm80(date1, date2)` : Form the matrix of precession/nutation for a given date, IAU 1976 precession model, IAU 1980 nutation model.

`pymNut80(date1, date2)` : Nutation, IAU 1980 model.

`pymNutm80(date1, date2)` : Form the matrix of nutation for a given date, IAU 1980 model.

`pymObl80(date1, date2)` : Mean obliquity of the ecliptic, IAU 1980 model.

`pymEqeq94(date1, date2)` : Equation of the equinoxes, IAU 1994 model.

`pymBi00()` : Frame bias components of IAU 2000 precession-nutation models;  part of the Mathews-Herring-Buffett (MHB2000) nutation series, with additions.

`pymPr00(date1, date2)` : Precession-rate part of the IAU 2000 precession-nutation models (part of MHB2000).

`pymBp00(date1, date2)` : Frame bias and precession, IAU 2000.

`pymPmat00(date1, date2)` : Precession matrix (including frame bias) from GCRS to a specified date, IAU 2000 model.

`pymEe00(date1, date2, epsa, dpsi)` : The equation of the equinoxes, compatible with IAU 2000 resolutions, given the nutation in longitude and the mean obliquity.

`pymEe00a(date1, date2)` : Equation of the equinoxes, compatible with IAU 2000 resolutions.

`pymEe00b(date1, date2)` : Equation of the equinoxes, compatible with IAU 2000 resolutions but using the truncated nutation model IAU 2000B.

`pymEect00(date1, date2)` : Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions.

`pymNut00a(date1, date2)` : Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation with free core nutation omitted).

`pymNut00b(date1, date2)` : Nutation, IAU 2000B model.

`pymPn00(date1, date2, dpsi, deps)` : Precession-nutation, IAU 2000 model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

`pymPn00a(date1, date2)` : Precession-nutation, IAU 2000A model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

`pymNum00a(date1, date2)` : Form the matrix of nutation for a given date, IAU 2000A model.

`pymPn00b(date1, date2)` : Precession-nutation, IAU 2000B model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

`pymNum00b(date1, date2)` : Form the matrix of nutation for a given date, IAU 2000B model.

`pymEra00(dj1, dj2)` : Earth rotation angle (IAU 2000 model).

`pymSp00(date1, date2)` : The TIO locator s', positioning the Terrestrial Intermediate Origin on the equator of the Celestial Intermediate Pole.

`pymPom00(xp, yp, sp)` : Form the matrix of polar motion for a given date, IAU 2000.

`pymPnm00a(date1, date2)` : Form the matrix of precession-nutation for a given date (including frame bias), equinox based, IAU 2000A model.

`pymPnm00b(date1, date2)` : Form the matrix of precession-nutation for a given date (including frame bias), equinox-based, IAU 2000B model.

`pymS00(date1, date2, x, y)` : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, given the CIP's X,Y coordinates.  Compatible with IAU 2000A precession-nutation.

`pymS00a(date1, date2)` : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2000A precession-nutation model.

`pymS00b(date1, date2)` : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2000B precession-nutation model.

`pymXys00a(date1, date2)` : For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole and the CIO locator s, using the IAU 2000A precession-nutation model.

`pymXys00b(date1, date2)` : For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole and the CIO locator s, using the IAU 2000B precession-nutation model.

`pymBp06(date1, date2)` : Frame bias and precession, IAU 2006.

`pymNut06a(date1, date2)` : IAU 2000A nutation with adjustments to match the IAU 2006 precession.

`pymObl06(date1, date2)` : Mean obliquity of the ecliptic, IAU 2006 precession model.

`pymPfw06(date1, date2)` : Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

`pymEe06a(date1, date2)` : Equation of the equinoxes, compatible with IAU 2000 resolutions and IAU 2006/2000A precession-nutation.

`pymPn06(date1, date2, dpsi, deps)` : Precession-nutation, IAU 2006 model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

`pymPn06a(date1, date2)` : Precession-nutation, IAU 2006/2000A models:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

`pymNum06a(date1, date2)` : Form the matrix of nutation for a given date, IAU 2006/2000A model.

`pymPnm06a(date1, date2)` : Form the matrix of precession-nutation for a given date (including frame bias), equinox based, IAU 2006 precession and IAU 2000A nutation models.

`pymPb06(date1, date2)` : This function forms three Euler angles which implement general precession from epoch J2000.0, using the IAU 2006 model.  Frame bias (the offset between ICRS and mean J2000.0) is included.

`pymPmat06(date1, date2)` : Precession matrix (including frame bias) from GCRS to a specified date, IAU 2006 model.

`pymP06e(date1, date2)` : Precession angles, IAU 2006, equinox based.

`pymS06(date1, date2, x, y)` : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, given the CIP's X,Y coordinates.  Compatible with IAU 2006/2000A precession-nutation.

`pymS06a(date1, date2)` : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2006 precession and IAU 2000A nutation models.

`pymXy06(date1, date2)` : X,Y coordinates of celestial intermediate pole from series based on IAU 2006 precession and IAU 2000A nutation.

`pymXys06a(date1, date2)` : For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole and the CIO locator s, using the IAU 2006 precession and IAU 2000A nutation models.

`pymFw2m(gamb, phib, psi, eps)` : Form rotation matrix given the Fukushima-Williams angles.

`pymFw2xy(gamb, phib, psi, eps)` : CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

```python
import PyMsOfa as sf

EPSA =  0.4090789763356509900        #mean obliquity of date
DPSI = -0.9630909107115582393e-5    #nutation
DEPS =  0.4063239174001678826e-4    #nutation
EPJ = 1666.666                #Julian epoch (TT)

print(sf.pymNumat(EPSA, DPSI, DEPS))
print(sf.pymLtp(EPJ))
print(sf.pymLtpb(EPJ))
print(sf.pymLtpequ(EPJ))
print(sf.pymLtpecl(EPJ))
```

```text
[[ 1.00000000e+00  8.83623932e-06  3.83083345e-06]
 [-8.83608366e-06  9.99999999e-01 -4.06324087e-05]
 [-3.83119248e-06  4.06323748e-05  9.99999999e-01]]
[[ 0.99670441  0.07437802  0.03237624]
 [-0.07437803  0.99722939 -0.00120577]
 [-0.03237622 -0.00120629  0.99947502]]
[[ 0.99670442  0.07437795  0.03237633]
 [-0.07437796  0.99722939 -0.00120574]
 [-0.03237631 -0.00120632  0.99947502]]
[-0.03237622 -0.00120629  0.99947502]
[-5.73709164e-05 -3.98473355e-01  9.17179907e-01]
```

```python
import PyMsOfa as sf

DATE01 = 2400000.5
DATE02 = 33282.0
DATE11 = 2400000.5
DATE12 = 51544.0
DATE1 = 2400000.5
DATE2 = 50124.0

print(sf.pymPrec76(DATE01, DATE02, DATE11, DATE12))
print(sf.pymPmat76(DATE1, DATE2))
print(sf.pymPnm80(DATE1, DATE2))
print(sf.pymNut80(DATE1, DATE2))
print(sf.pymNutm80(DATE1, DATE2))
print(sf.pymObl80(DATE1, DATE2))
print(sf.pymEqeq94(DATE1, DATE2))
```

```text
(0.005588961642000161, 0.005589922365870681, 0.004858945471687296)
[[ 9.99999550e-01  8.69663160e-04  3.77915321e-04]
 [-8.69663160e-04  9.99999622e-01 -1.64328454e-07]
 [-3.77915321e-04 -1.64330651e-07  9.99999929e-01]]
[[ 9.99999583e-01  8.37365365e-04  3.63912174e-04]
 [-8.37380450e-04  9.99999649e-01  4.13020343e-05]
 [-3.63877462e-04 -4.16067501e-05  9.99999933e-01]]
(3.520279437486982e-05, -4.1454385835936036e-05)
[[ 9.99999999e-01 -3.22978088e-05 -1.40031525e-05]
 [ 3.22983892e-05  9.99999999e-01  4.14541597e-05]
 [ 1.40018136e-05 -4.14546119e-05  9.99999999e-01]]
0.40910163117239295
3.229357405996896e-05
```

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 53736
EPSA = 0.4090789763356509900        #mean obliquity
DPSI = -0.9630909107115582393e-5    #nutation in longitude

print(sf.pymBi00())
print(sf.pymPr00(DATE1, DATE2))
print(sf.pymBp00(DATE1, DATE2))
print(sf.pymPmat00(DATE1, DATE2))
print(sf.pymEe00(DATE1, DATE2, EPSA, DPSI))
print(sf.pymEe00a(DATE1, DATE2))
print(sf.pymEect00(DATE1, DATE2))
```

```text
(-2.0253091528350866e-07, -3.3060414542221477e-08, -7.078279744199226e-08)
(-8.716465172668348e-08, -7.342018386722812e-09)
(array([[ 1.00000000e+00, -7.07827974e-08,  8.05621715e-08],
       [ 7.07827948e-08,  1.00000000e+00,  3.30604145e-08],
       [-8.05621738e-08, -3.30604088e-08,  1.00000000e+00]]), array([[ 9.99998930e-01, -1.34164723e-03, -5.82988093e-04],
       [ 1.34164723e-03,  9.99999100e-01, -3.83744444e-07],
       [ 5.82988083e-04, -3.98420327e-07,  9.99999830e-01]]), array([[ 9.99998930e-01, -1.34171799e-03, -5.82907575e-04],
       [ 1.34171801e-03,  9.99999100e-01, -3.50575973e-07],
       [ 5.82907521e-04, -4.31521995e-07,  9.99999830e-01]]))
[[ 9.99998930e-01 -1.34171799e-03 -5.82907575e-04]
 [ 1.34171801e-03  9.99999100e-01 -3.50575973e-07]
 [ 5.82907521e-04 -4.31521995e-07  9.99999830e-01]]
-8.834193235367966e-06
-8.83419245922342e-06
2.0460850048851376e-09
```

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 53736.0
DPSI = -0.9632552291149335877e-5    #nutation
DEPS =  0.4063197106621141414e-4    #nutation

print(sf.pymNut00a(DATE1, DATE2))
print(sf.pymNut00b(DATE1, DATE2))
print(sf.pymPn00(DATE1, DATE2, DPSI, DEPS))
print(sf.pymPn00a(DATE1, DATE2))
print(sf.pymNum00a(DATE1, DATE2))
print(sf.pymPn00b(DATE1, DATE2))
print(sf.pymNum00b(DATE1, DATE2))
```

```text
(-9.630909107116423e-06, 4.063239174001665e-05)
(-9.632552291148323e-06, 4.063197106621164e-05)
(0.409079178940423, array([[ 1.00000000e+00, -7.07827974e-08,  8.05621715e-08],
       [ 7.07827948e-08,  1.00000000e+00,  3.30604145e-08],
       [-8.05621738e-08, -3.30604088e-08,  1.00000000e+00]]), array([[ 9.99998930e-01, -1.34164723e-03, -5.82988093e-04],
       [ 1.34164723e-03,  9.99999100e-01, -3.83744444e-07],
       [ 5.82988083e-04, -3.98420327e-07,  9.99999830e-01]]), array([[ 9.99998930e-01, -1.34171799e-03, -5.82907575e-04],
       [ 1.34171801e-03,  9.99999100e-01, -3.50575973e-07],
       [ 5.82907521e-04, -4.31521995e-07,  9.99999830e-01]]), array([[ 1.00000000e+00,  8.83774614e-06,  3.83148884e-06],
       [-8.83759046e-06,  9.99999999e-01, -4.06319880e-05],
       [-3.83184793e-06,  4.06319541e-05,  9.99999999e-01]]), array([[ 9.99998944e-01, -1.33288025e-03, -5.79076090e-04],
       [ 1.33285675e-03,  9.99999111e-01, -4.09774056e-05],
       [ 5.79130193e-04,  4.02055368e-05,  9.99999831e-01]]))
(-9.630909107116423e-06, 4.063239174001665e-05, 0.409079178940423, array([[ 1.00000000e+00, -7.07827974e-08,  8.05621715e-08],
       [ 7.07827948e-08,  1.00000000e+00,  3.30604145e-08],
       [-8.05621738e-08, -3.30604088e-08,  1.00000000e+00]]), array([[ 9.99998930e-01, -1.34164723e-03, -5.82988093e-04],
       [ 1.34164723e-03,  9.99999100e-01, -3.83744444e-07],
       [ 5.82988083e-04, -3.98420327e-07,  9.99999830e-01]]), array([[ 9.99998930e-01, -1.34171799e-03, -5.82907575e-04],
       [ 1.34171801e-03,  9.99999100e-01, -3.50575973e-07],
       [ 5.82907521e-04, -4.31521995e-07,  9.99999830e-01]]), array([[ 1.00000000e+00,  8.83623854e-06,  3.83083524e-06],
       [-8.83608288e-06,  9.99999999e-01, -4.06324087e-05],
       [-3.83119427e-06,  4.06323748e-05,  9.99999999e-01]]), array([[ 9.99998944e-01, -1.33288176e-03, -5.79076743e-04],
       [ 1.33285825e-03,  9.99999111e-01, -4.09778271e-05],
       [ 5.79130847e-04,  4.02059566e-05,  9.99999831e-01]]))
[[ 1.00000000e+00  8.83623854e-06  3.83083524e-06]
 [-8.83608288e-06  9.99999999e-01 -4.06324087e-05]
 [-3.83119427e-06  4.06323748e-05  9.99999999e-01]]
(-9.632552291148323e-06, 4.063197106621164e-05, 0.409079178940423, array([[ 1.00000000e+00, -7.07827974e-08,  8.05621715e-08],
       [ 7.07827948e-08,  1.00000000e+00,  3.30604145e-08],
       [-8.05621738e-08, -3.30604088e-08,  1.00000000e+00]]), array([[ 9.99998930e-01, -1.34164723e-03, -5.82988093e-04],
       [ 1.34164723e-03,  9.99999100e-01, -3.83744444e-07],
       [ 5.82988083e-04, -3.98420327e-07,  9.99999830e-01]]), array([[ 9.99998930e-01, -1.34171799e-03, -5.82907575e-04],
       [ 1.34171801e-03,  9.99999100e-01, -3.50575973e-07],
       [ 5.82907521e-04, -4.31521995e-07,  9.99999830e-01]]), array([[ 1.00000000e+00,  8.83774614e-06,  3.83148884e-06],
       [-8.83759046e-06,  9.99999999e-01, -4.06319880e-05],
       [-3.83184793e-06,  4.06319541e-05,  9.99999999e-01]]), array([[ 9.99998944e-01, -1.33288025e-03, -5.79076090e-04],
       [ 1.33285675e-03,  9.99999111e-01, -4.09774056e-05],
       [ 5.79130193e-04,  4.02055368e-05,  9.99999831e-01]]))
[[ 1.00000000e+00  8.83774614e-06  3.83148884e-06]
 [-8.83759046e-06  9.99999999e-01 -4.06319880e-05]
 [-3.83184793e-06  4.06319541e-05  9.99999999e-01]]
```

```python
import PyMsOfa as sf

DJ1 = 2400000.5
DJ2 = 54388.0
DATE1 = 2400000.5
DATE2 = 52541.0
XP = 2.55060238e-7            #coordinates of the pole (radians)
YP = 1.860359247e-6            #coordinates of the pole (radians)   
SP = -0.1367174580728891460e-10        #the TIO locator s' (radians)
X = 0.5791308486706011000e-3        #CIP coordinates 
Y = 0.4020579816732961219e-4        #CIP coordinates 

print(sf.pymEra00(DJ1, DJ2))
print(sf.pymSp00(DATE1, DATE2))
print(sf.pymPom00(XP, YP, SP))
print(sf.pymPnm00a(DATE1, DATE2))
print(sf.pymPnm00b(DATE1, DATE2))
print(sf.pymS00(DATE1, DATE2, X, Y))
print(sf.pymS00a(DATE1, DATE2))
print(sf.pymS00b(DATE1, DATE2))
print(sf.pymXys00a(DATE1, DATE2))
print(sf.pymXys00b(DATE1, DATE2))
```

```text
0.40228372400281387
-6.216698469981019e-12
[[ 1.00000000e+00 -1.36717458e-11  2.55060238e-07]
 [ 1.41462495e-11  1.00000000e+00 -1.86035925e-06]
 [-2.55060238e-07  1.86035925e-06  1.00000000e+00]]
[[ 9.99999829e-01 -5.35976002e-04 -2.32864778e-04]
 [ 5.35972044e-04  9.99999856e-01 -1.70602782e-05]
 [ 2.32873889e-04  1.69354663e-05  9.99999973e-01]]
[[ 9.99999829e-01 -5.35976719e-04 -2.32865089e-04]
 [ 5.35972761e-04  9.99999856e-01 -1.70612293e-05]
 [ 2.32874200e-04  1.69364170e-05  9.99999973e-01]]
-2.3077139553642807e-08
-1.3406844489207085e-08
-1.3406957829526553e-08
(0.00023287388885713626, 1.693546625002907e-05, -1.3406844489207085e-08)
(0.00023287420038471495, 1.6936416998432182e-05, -1.3406957829526553e-08)
```

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 50124.0
DPSI = -0.9632552291149335877e-5    #nutation 
DEPS =  0.4063197106621141414e-4    #nutation 

print(sf.pymBp06(DATE1, DATE2))
print(sf.pymNut06a(DATE1, DATE2))
print(sf.pymObl06(DATE1, DATE2))
print(sf.pymPfw06(DATE1, DATE2))
print(sf.pymEe06a(DATE1, DATE2))
print(sf.pymPn06(DATE1, DATE2, DPSI, DEPS))
print(sf.pymPn06a(DATE1, DATE2))
```

```text
(array([[ 1.00000000e+00, -7.07836896e-08,  8.05621398e-08],
       [ 7.07836869e-08,  1.00000000e+00,  3.30594373e-08],
       [-8.05621421e-08, -3.30594317e-08,  1.00000000e+00]]), array([[ 9.99999550e-01,  8.69611197e-04,  3.77892903e-04],
       [-8.69611195e-04,  9.99999622e-01, -1.69164593e-07],
       [-3.77892907e-04, -1.59455381e-07,  9.99999929e-01]]), array([[ 9.99999551e-01,  8.69540401e-04,  3.77973494e-04],
       [-8.69540411e-04,  9.99999622e-01, -1.36175226e-07],
       [-3.77973469e-04, -1.92488062e-07,  9.99999929e-01]]))
(3.520755263333626e-05, -4.145326347441321e-05)
0.409101431658115
(-2.2433875313873983e-06, 0.4091014602385099, -0.0009501953509247532, 0.409101431658115)
3.229806317150974e-05
(0.409101431658115, array([[ 1.00000000e+00, -7.07836896e-08,  8.05621398e-08],
       [ 7.07836869e-08,  1.00000000e+00,  3.30594373e-08],
       [-8.05621421e-08, -3.30594317e-08,  1.00000000e+00]]), array([[ 9.99999550e-01,  8.69611197e-04,  3.77892903e-04],
       [-8.69611195e-04,  9.99999622e-01, -1.69164593e-07],
       [-3.77892907e-04, -1.59455381e-07,  9.99999929e-01]]), array([[ 9.99999551e-01,  8.69540401e-04,  3.77973494e-04],
       [-8.69540411e-04,  9.99999622e-01, -1.36175226e-07],
       [-3.77973469e-04, -1.92488062e-07,  9.99999929e-01]]), array([[ 1.00000000e+00,  8.83766088e-06,  3.83168550e-06],
       [-8.83750519e-06,  9.99999999e-01, -4.06319880e-05],
       [-3.83204459e-06,  4.06319541e-05,  9.99999999e-01]]), array([[ 9.99999541e-01,  8.78378057e-04,  3.81805178e-04],
       [-8.78362554e-04,  9.99999613e-01, -4.07715007e-05],
       [-3.81840843e-04,  4.04361186e-05,  9.99999926e-01]]))
(3.520755263333626e-05, -4.145326347441321e-05, 0.409101431658115, array([[ 1.00000000e+00, -7.07836896e-08,  8.05621398e-08],
       [ 7.07836869e-08,  1.00000000e+00,  3.30594373e-08],
       [-8.05621421e-08, -3.30594317e-08,  1.00000000e+00]]), array([[ 9.99999550e-01,  8.69611197e-04,  3.77892903e-04],
       [-8.69611195e-04,  9.99999622e-01, -1.69164593e-07],
       [-3.77892907e-04, -1.59455381e-07,  9.99999929e-01]]), array([[ 9.99999551e-01,  8.69540401e-04,  3.77973494e-04],
       [-8.69540411e-04,  9.99999622e-01, -1.36175226e-07],
       [-3.77973469e-04, -1.92488062e-07,  9.99999929e-01]]), array([[ 9.99999999e-01, -3.23021772e-05, -1.40050388e-05],
       [ 3.23027577e-05,  9.99999999e-01,  4.14530373e-05],
       [ 1.40036998e-05, -4.14534896e-05,  9.99999999e-01]]), array([[ 9.99999583e-01,  8.37238238e-04,  3.63968460e-04],
       [-8.37253335e-04,  9.99999649e-01,  4.13290687e-05],
       [-3.63933730e-04, -4.16337852e-05,  9.99999933e-01]]))
```

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 50124.0
X = 0.5791308486706011000e-3        #CIP coordinates
Y = 0.4020579816732961219e-4        #CIP coordinates
GAMB = -0.2243387670997992368e-5    #F-W angle gamma_bar (radians)
PHIB =  0.4091014602391312982        #F-W angle phi_bar (radians)
PSI  = -0.9501954178013015092e-3    #F-W angle psi (radians)
EPS  =  0.4091014316587367472        #F-W angle epsilon (radians)

print(sf.pymNum06a(DATE1, DATE2))
print(sf.pymPnm06a(DATE1, DATE2))
print(sf.pymPb06(DATE1, DATE2))
print(sf.pymPmat06(DATE1, DATE2))
print(sf.pymP06e(DATE1, DATE2))
print(sf.pymS06(DATE1, DATE2, X, Y))
print(sf.pymXy06(DATE1, DATE2))
print(sf.pymXys06a(DATE1, DATE2))
print(sf.pymFw2m(GAMB, PHIB, PSI, EPS))
print(sf.pymFw2xy(GAMB, PHIB, PSI, EPS))
```

```text
[[ 9.99999999e-01 -3.23021772e-05 -1.40050388e-05]
 [ 3.23027577e-05  9.99999999e-01  4.14530373e-05]
 [ 1.40036998e-05 -4.14534896e-05  9.99999999e-01]]
[[ 9.99999583e-01  8.37238238e-04  3.63968460e-04]
 [-8.37253335e-04  9.99999649e-01  4.13290687e-05]
 [-3.63933730e-04 -4.16337852e-05  9.99999933e-01]]
(array(-0.00050926), array(-0.00036028), array(-0.00037797))
[[ 9.99999551e-01  8.69540401e-04  3.77973494e-04]
 [-8.69540411e-04  9.99999622e-01 -1.36175226e-07]
 [-3.77973469e-04 -1.92488062e-07  9.99999929e-01]]
(0.4090926006005829, -0.0009500121640916875, 0.4090926058345976, -7.903154196964373e-07, 8.826577573750554e-06, -8.861888526875804e-06, 3.0522926386722435, 0.409101431658115, -2.007869257151149e-06, -0.00044765239229031437, -0.00042195894411613553, -0.00037789294961742447, -0.0009481699832898824, -1.9867900599926174e-06, 0.40910142717887715, -0.0009499928243453799)
-7.773082059855285e-09
(-0.0003639337295849517, -4.163378483283794e-05)
(-0.00036393372975485804, -4.163378524520265e-05, -3.7068424268427255e-09)
[[ 9.99999551e-01  8.69540462e-04  3.77973520e-04]
 [-8.69540472e-04  9.99999622e-01 -1.36175250e-07]
 [-3.77973496e-04 -1.92488085e-07  9.99999929e-01]]
(-0.00037797349570340826, -1.9248808483033084e-07)
```

## 4. Coordinate-system transformations

`pymEform(n)` : Earth reference ellipsoids.(1 : WGS84 ; 2 : GRS80 ; 3 : WGS72)

`pymGd2gce(a, f, elong, phi, height)` : Transform geodetic coordinates to geocentric for a reference ellipsoid of specified form.

`pymGd2gc(n, elong, phi, height)` : Transform geodetic coordinates to geocentric using the specified reference ellipsoid.

`pymGc2gde(a, f, xyz)` : Transform geocentric coordinates to geodetic for a reference ellipsoid of specified form.

`pymGc2gd(n, xyz)` : Transform geocentric coordinates to geodetic using the specified reference ellipsoid.

`pymAe2hd(az, el, phi)` : Horizon to equatorial coordinates:  transform azimuth and altitude to hour angle and declination.

`pymHd2ae(ha, dec, phi)` : Equatorial to horizon coordinates:  transform hour angle and declination to azimuth and altitude.

`pymHd2pa(ha, dec, phi)` : Parallactic angle for a given hour angle and declination.

`pymLtecm(epj)` : ICRS equatorial to ecliptic rotation matrix, long-term.

`pymLteqec(epj, ra, dec)` : Transformation from ICRS equatorial coordinates to ecliptic coordinates (mean equinox and ecliptic of date) using a long-term precession model.

`pymLteceq(epj, dl, db)` : Transformation from ecliptic coordinates (mean equinox and ecliptic of date) to ICRS RA,Dec, using a long-term precession model.

`pymIcrs2g(dr, dd)` : Transformation from ICRS to Galactic Coordinates.

`pymG2icrs(dl, db)` : Transformation from Galactic Coordinates to ICRS.

`pymFk425(r1950, d1950, dr1950, dd1950, p1950, v1950)` : Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

`pymFk45z(r1950, d1950, bepoch)` : Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero proper motion in the FK5 system.

`pymFk524(r2000, d2000, dr2000, dd2000, p2000, v2000)` : Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

`pymFk54z(r2000, d2000, bepoch)` : Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero proper motion in FK5 and parallax.

`pymFk5hip()` : FK5 to Hipparcos rotation and spin.

`pymFk52h(r5, d5, dr5, dd5, px5, rv5)` : Transform FK5 (J2000.0) star data into the Hipparcos system.

`pymFk5hz(r5, d5, date1, date2)` : Transform an FK5 (J2000.0) star position into the system of the Hipparcos catalog, assuming zero Hipparcos proper motion.

`pymH2fk5(rh, dh, drh, ddh, pxh, rvh)` : Transform Hipparcos star data into the FK5 (J2000.0) system.

`pymHfk5z(rh, dh, date1, date2)` : Transform a Hipparcos star position into FK5 J2000.0, assuming = zero Hipparcos proper motion.

`pymEcm06(date1, date2)` : ICRS equatorial to ecliptic rotation matrix, IAU 2006.

`pymEceq06(date1, date2, dl, db)` : Transformation from ecliptic coordinates (mean equinox and ecliptic of date) to ICRS RA,Dec, using the IAU 2006 precession model.

`pymEqec06(date1, date2, dr, dd)` : Transformation from ICRS equatorial coordinates to ecliptic coordinates (mean equinox and ecliptic of date) using IAU 2006 precession model.

`pymStarpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b)` : Star proper motion:  update star catalog data for space motion.

`pymPmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b)` : Star proper motion:  update star catalog data for space motion, with special handling to handle the zero parallax case.

```python
import PyMsOfa as sf
import math as ma

N = 2                 # ellipsoid identifier
A = 6378136.0            #equatorial radius
F = 0.0033528            #flattening
ELONG = 3.1            #longitude (radians, east +ve)
PHI = -0.5            #latitude (geodetic, radians)
HEIGHT = 2500.0            #height above ellipsoid (geodetic)
XYZ = [2e6, 3e6, 5.244e6]    #geocentric vector
AZ = 5.5            #azimuth
EL = 1.1            #altitude (informally, elevation)
PHI = 0.7            #site latitude
HA = 1.1            #hour angle (local)
DEC = 1.2            #declination
PHI = 0.3            #site latitude
EPJ = 2024.154688569473        #Julian epoch (TT)
DR = ma.pi/3            #ICRS right ascension (radians) 
DD = ma.pi/4            #ICRS declination (radians)
DL = 1.5            #ecliptic longitude (radians)
DB = 0.6            #ecliptic latitude (radians)

print(sf.pymEform(N))
print(sf.pymGd2gce(A, F, ELONG, PHI, HEIGHT))
print(sf.pymGd2gc(N, ELONG, PHI, HEIGHT))
print(sf.pymGc2gde(A, F, XYZ))
print(sf.pymGc2gd(N, XYZ))
print(sf.pymAe2hd(AZ, EL, PHI))
print(sf.pymHd2ae(HA, DEC, PHI))
print(sf.pymHd2pa(HA, DEC, PHI))
print(sf.pymLtecm(EPJ))
print(sf.pymLteqec(EPJ, DR, DD))
print(sf.pymLteceq(EPJ, DL, DB))
```

```text
(6378137.0, 0.003352810681182319)
[-6092162.97213779   253535.44209086  1873536.4160412 ]
[-6092163.93259208   253535.48206175  1873536.67126909]
(0.982793723247329, 0.9716018377570411, 332.36862495798647)
(0.982793723247329, 0.9716018482060784, 331.41731754978827)
(0.40025365757863, 0.6070689395289157)
(5.916889243730066, 0.4472186304990486)
1.9062274280019955
[[ 9.99982658e-01 -5.40151427e-03 -2.34679824e-03]
 [ 5.88929492e-03  9.17488028e-01  3.97719794e-01]
 [ 4.87015078e-06 -3.97726717e-01  9.17503928e-01]]
(1.1796084210776348, 0.4172173738451268)
(1.4521514559004196, 1.0072288157204605)
```

```python
import PyMsOfa as sf
import math as ma

DR = ma.pi/3                #ICRS right ascension (radians)
DD = ma.pi/4                #ICRS declination (radians)
DL = ma.pi/6                #galactic longitude (radians)
DB = ma.pi/8                #galactic longitude (radians)
R1950 = 0.07626899753879587532        #B1950.0 RA (rad)
D1950 = -1.137405378399605780        #B1950.0 Dec (rad)
DR1950 = 0.1973749217849087460e-4    #B1950.0 proper motions (rad/trop.yr)
DD1950 = 0.5659714913272723189e-5    #B1950.0 proper motions (rad/trop.yr)
P1950 = 0.134                #parallax (arcsec)
V1950 = 8.7                #radial velocity (km/s, +ve = moving away)
R2000 = 0.01602284975382960982        #J2000.0 RA (rad)
D2000 = -0.1164347929099906024        #J2000.0 Dec (rad)
BEPOCH = 1954.677617625256806        #Besselian epoch (e.g. 1979.3)
DR2000 = 0.2019447755430472323e-4    #J2000.0 proper motions (rad/Jul.yr)
DD2000 = 0.3541563940505160433e-5    #J2000.0 proper motions (rad/Jul.yr)
P2000 = 0.1559                #parallax (arcsec)
V2000 = 86.87                #radial velocity (km/s, +ve = moving away)
R5  =  1.76779433            #RA (radians)     FK5, equinox J2000.0, epoch J2000.0
D5  = -0.2917517103            #Dec (radians)  FK5, equinox J2000.0, epoch J2000.0
DR5 = -1.91851572e-7            #proper motion in RA (dRA/dt, rad/Jyear)     FK5, equinox J2000.0, epoch J2000.0
DD5 = -5.8468475e-6            #proper motion in Dec (dDec/dt, rad/Jyear)    FK5, equinox J2000.0, epoch J2000.0
PX5 =  0.379210                #parallax (arcsec)    FK5, equinox J2000.0, epoch J2000.0
RV5 = -7.6                #radial velocity (km/s, positive = receding)    FK5, equinox J2000.0, epoch J2000.0

print(sf.pymIcrs2g(DR,DD))
print(sf.pymG2icrs(DL, DB))
print(sf.pymFk425(R1950, D1950, DR1950, DD1950, P1950, V1950))
print(sf.pymFk45z(R1950, D1950, BEPOCH))
print(sf.pymFk524(R2000, D2000, DR2000, DD2000, P2000, V2000))
print(sf.pymFk54z(R2000, D2000, BEPOCH))
print(sf.pymFk5hip())
print(sf.pymFk52h(R5, D5, DR5, DD5, PX5, RV5))
```

```text
(2.6938731316889215, -0.10583104427828438)
(4.562820525414904, 0.13283058755100807)
(0.08757989933556447, -1.132279113042092, 1.953670614474396e-05, 5.637686678659639e-06, 0.1339919950582768, 8.736999669183527)
(0.08660249114642907, -1.1325609304456188)
(0.0038357827667403516, -0.12147115357961684, 2.022271874289751e-05, 3.5660786624684333e-06, 0.15600796032308914, 86.73793467537503)
(0.004846756109282368, -0.12129386288333069, -1.1732616691106914e-08, 2.109122543091408e-08)
(array([[ 1.00000000e+00,  1.11022335e-07,  4.41180396e-08],
       [-1.11022331e-07,  1.00000000e+00, -9.64779250e-08],
       [-4.41180503e-08,  9.64779201e-08,  1.00000000e+00]]), array([-1.45444104e-09,  2.90888209e-09,  3.39369577e-09]))
(1.7677942262999478, -0.29175160705303915, -1.9618741256057196e-07, -5.845990517669392e-06, 0.37921000000000005, -7.600000094000026)
```

```python
import PyMsOfa as sf

R5 =  1.76779433        #FK5 RA (radians), equinox J2000.0, at date
D5 = -0.2917517103        #FK5 Dec (radians), equinox J2000.0, at date
DATE1 = 2400000.5
DATE2 = 54479.0
RH  =  1.767794352        #RA (radians)        Hipparcos, epoch J2000.0
DH  = -0.2917512594        #Dec (radians)        Hipparcos, epoch J2000.0
DRH = -2.76413026e-6        #proper motion in RA (dRA/dt, rad/Jyear)    Hipparcos, epoch J2000.0
DDH = -5.92994449e-6        #proper motion in Dec (dDec/dt, rad/Jyear)    Hipparcos, epoch J2000.0
PXH =  0.379210            #parallax (arcsec)      Hipparcos, epoch J2000.0
RVH = -7.6            #radial velocity (km/s, positive = receding)    Hipparcos, epoch J2000.0
DL = 5.1            #ecliptic longitude (radians)
DB = -0.9            #ecliptic latitude (radians)
DR = 1.234            #ICRS right ascension (radians)
DD = 0.987            #ICRS declination (radians)
RA1 =   0.01686756        #right ascension (radians), before
DEC1 = -1.093989828        #declination (radians), before
PMR1 = -1.78323516e-5        #RA proper motion (radians/year), before
PMD1 =  2.336024047e-6        #Dec proper motion (radians/year), before
PX1 =   0.74723            #parallax (arcseconds), before
RV1 = -21.6            #radial velocity (km/s, +ve = receding), before
EP1A = 2400000.5        #"before" epoch, part A
EP1B = 50083.0            #"before" epoch, part B
EP2A = 2400000.5        #"after" epoch, part A
EP2B = 53736.0            #"after" epoch, part B

print(sf.pymFk5hz(R5, D5, DATE1, DATE2))
print(sf.pymH2fk5(RH, DH, DRH, DDH, PXH, RVH))
print(sf.pymHfk5z(RH, DH, DATE1, DATE2))
print(sf.pymEcm06(DATE1, DATE2))
print(sf.pymEceq06(DATE1, DATE2, DL, DB))
print(sf.pymEqec06(DATE1, DATE2, DR, DD))
print(sf.pymStarpm(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))
print(sf.pymPmsafe(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))
```

```text
(1.767794191464424, -0.2917516001679885)
(1.7677944557000653, -0.2917513626469639, -2.759794502451121e-06, -5.930801409326283e-06, 0.37920999999999994, -7.600000130907113)
(1.767794490535581, -0.2917513695320114, 4.335890983539242e-09, -8.56964884123775e-10)
[[ 9.99998081e-01 -1.79659635e-03 -7.80558614e-04]
 [ 1.95883328e-03  9.17487623e-01  3.97759506e-01]
 [ 1.53958984e-06 -3.97760272e-01  9.17489382e-01]]
(5.535423121517396, -1.2462214232485105)
(1.3513428093792499, 0.592703434726813)
(0.01668919069414256, -1.093966454217128, -1.7836626821531766e-05, 2.33809291598399e-06, 0.7473533835317719, -21.599051704764175)
(0.01668919069414256, -1.093966454217128, -1.7836626821531766e-05, 2.33809291598399e-06, 0.7473533835317719, -21.599051704764175)
```

## 5. Astrometry parameters

### Astrometry parameters (general)

`pymEpv00(date1, date2)` : Earth position and velocity, heliocentric and barycentric, with respect to the Barycentric Celestial Reference System.

`pymAb(pnat, v, s, bm1)` : Apply aberration to transform natural direction into proper direction.

`pymPmpx(rc, dc, pr, pd, px, rv, pmt, pob)` : The position of the star after a certain time interval (PMT) is given based on the star's position, proper motion, parallax, and radial velocity and the observer's velocity relative to the Solar SystemBarycenter.

`pymRefco(phpa, tc, rh, wl)` : Determine the constants A and B in the atmospheric refraction model dZ = A tan Z + B tan^3 Z.

`pymPvtob(elong, phi, hm, xp, yp, sp, theta)` : Position and velocity of a terrestrial observing station.

`pymC2ixy(date1, date2, x, y)` : Form the celestial to intermediate-frame-of-date matrix for a given date when the CIP X,Y coordinates are known.  IAU 2000.

`pymC2ixys(x, y, s)` : Form the celestial to intermediate-frame-of-date matrix given the CIP X,Y and the CIO locator s.

`pymC2ibpn(date1, date2, rbpn)` : Form the celestial-to-intermediate matrix for a given date given the bias-precession-nutation matrix.  IAU 2000.

`pymC2i00a(date1, date2)` : Form the celestial-to-intermediate matrix for a given date using the IAU 2000A precession-nutation model.

`pymC2i00b(date1, date2)` : Form the celestial-to-intermediate matrix for a given date using the IAU 2000B precession-nutation model.

`pymC2i06a(date1, date2)` : Form the celestial-to-intermediate matrix for a given date using the IAU 2006 precession and IAU 2000A nutation models.

`pymEors(rnpb, s)` : Equation of the origins, given the classical NPB matrix and the quantity s.

`pymEo06a(date1, date2)` : Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

`pymBpn2xy(rbpn)` : Extract from the bias-precession-nutation matrix the X,Y coordinates of the Celestial Intermediate Pole.

`pymC2tcio(rc2i, era, rpom)` : Assemble the celestial to terrestrial matrix from CIO-based components (the celestial-to-intermediate matrix, the Earth Rotation Angle and the polar motion matrix).

`pymC2t00a(tta, ttb, uta, utb, xp, yp)` : Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2000A precession-nutation model.

`pymC2t00b(tta, ttb, uta, utb, xp, yp)` : Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2000B precession-nutation model.

`pymC2t06a(tta, ttb, uta, utb, xp, yp)` : Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2006/2000A precession-nutation model.

`pymC2teqx(rbpn, gst, rpom)` : Assemble the celestial to terrestrial matrix from equinox-based components (the celestial-to-true matrix, the Greenwich Apparent Sidereal Time and the polar motion matrix).

`pymC2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp)` : Form the celestial to terrestrial matrix given the date, the UT1, the nutation and the polar motion.  IAU 2000.

`pymC2txy(tta, ttb, uta, utb, x, y, xp, yp)` : Form the celestial to terrestrial matrix given the date, the UT1, the CIP coordinates and the polar motion.  IAU 2000.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 60367.0
PNAT = [-0.76321968546737951, -0.60869453983060384, -0.21676408580639883]    #natural direction to the source (unit vector)
V = [2.1044018893653786E-5, -8.9108923304429319E-5, -3.8633714797716569E-5]    #observer barycentric velocity in units of c
S = 0.99980921395708788            #distance between the Sun and the observer (au)
BM1 = 0.99999999506209258        #sqrt(1-|v|^2): reciprocal of Lorenz factor
RC = 1.234                #ICRS RA at catalog epoch (radians)
DC = 0.789                #ICRS Dec at catalog epoch (radians)
PR = 1e-5                #RA proper motion (radians/year)
PD = -2e-5                #Dec proper motion (radians/year)
PX = 1e-2                #parallax (arcsec)
RV = 10.0                #radial velocity (km/s, +ve if receding)
PMT = 8.75                #proper motion time interval (SSB, Julian years)
POB = [0.9, 0.4, 0.1]            #SSB to observer vector (au)
PHPA = 800.0                #pressure at the observer (hPa = millibar)
TC = 10.0                #ambient temperature at the observer (deg C)
RH = 0.9                #relative humidity at the observer (range 0-1)
WL = 0.4                #wavelength (micrometers)
ELONG = 2.0                #longitude (radians, east +ve)
PHI = 0.5                #latitude (geodetic, radians)
HM = 3000.0                #height above ref. ellipsoid (geodetic, m)
XP = 1e-6                #coordinates of the pole (radians)
YP = -0.5e-6                #coordinates of the pole (radians)
SP = 1e-8                #the TIO locator s' (radians)
THETA = 5.0                #Earth rotation angle (radians)

print(sf.pymEpv00(DATE1, DATE2))
print(sf.pymAb(PNAT, V, S, BM1))
print(sf.pymPmpx(RC, DC, PR, PD, PX, RV, PMT, POB))
print(sf.pymRefco(PHPA, TC, RH, WL))
print(sf.pymPvtob(ELONG, PHI, HM, XP, YP, SP, THETA))
```

```text
(array([[-0.91503673,  0.34697147,  0.15041286],
       [-0.00685123, -0.01464186, -0.00634635]]), array([[-0.92271488,  0.34386561,  0.14929229],
       [-0.00684596, -0.01464791, -0.00634904]]))
[-0.76316311 -0.60875531 -0.21679263]
[0.23281376 0.66510971 0.70952578]
(0.0002264949956241415, -2.598658261729344e-07)
[[ 4.22508137e+06  3.68194322e+06  3.04114940e+06]
 [-2.68491539e+02  3.08097798e+02  0.00000000e+00]]
```

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 60367.0
X = 0.5791308486706011000e-3        #Celestial Intermediate Pole
Y = 0.4020579816732961219e-4        #Celestial Intermediate Pole
S = -0.1220040848472271978e-7        #the CIO locator s
RBPN = [[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
        [2.516462370370876e-3, 9.999968329010883e-1, 4.006159587358310e-5],
        [1.093465510215479e-3, -4.281337229063151e-5, 9.999994012499173e-1]]
    #celestial-to-true matrix

print(sf.pymC2ixy(DATE1, DATE2, X, Y))
print(sf.pymC2ixys(X, Y, S))
print(sf.pymC2ibpn(DATE1, DATE2, RBPN))
print(sf.pymC2i00a(DATE1, DATE2))
print(sf.pymC2i00b(DATE1, DATE2))
print(sf.pymC2i06a(DATE1, DATE2))
```

```text
[[ 9.99999832e-01  4.05562419e-09 -5.79130849e-04]
 [-2.73400415e-08  9.99999999e-01 -4.02057891e-05]
 [ 5.79130849e-04  4.02057982e-05  9.99999831e-01]]
[[ 9.99999832e-01  5.58198490e-10 -5.79130849e-04]
 [-2.38426164e-08  9.99999999e-01 -4.02057911e-05]
 [ 5.79130849e-04  4.02057982e-05  9.99999831e-01]]
[[ 9.99999402e-01  4.05563220e-09 -1.09346551e-03]
 [ 4.27593162e-08  9.99999999e-01  4.28133511e-05]
 [ 1.09346551e-03 -4.28133723e-05  9.99999401e-01]]
[[ 9.99997265e-01  4.05556385e-09 -2.33899017e-03]
 [-9.36612932e-08  9.99999999e-01 -3.83094698e-05]
 [ 2.33899016e-03  3.83095841e-05  9.99997264e-01]]
[[ 9.99997265e-01  4.05556385e-09 -2.33899017e-03]
 [-9.36612835e-08  9.99999999e-01 -3.83094655e-05]
 [ 2.33899017e-03  3.83095798e-05  9.99997264e-01]]
[[ 9.99997265e-01  4.05597547e-09 -2.33898963e-03]
 [-9.36602025e-08  9.99999999e-01 -3.83088363e-05]
 [ 2.33898963e-03  3.83089506e-05  9.99997264e-01]]
```

```python
import PyMsOfa as sf

RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
    #classical nutation x precession x bias matrix
S = -0.1220040848472271978e-7    #the quantity s (the CIO locator) in radians
DATE1 = 2400000.5
DATE2 = 60367.0
RBPN = [[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
        [2.516462370370876e-3, 9.999968329010883e-1, 4.006159587358310e-5],
        [1.093465510215479e-3, -4.281337229063151e-5, 9.999994012499173e-1]]
    #celestial-to-true matrix

print(sf.pymEors(RNPB, S))
print(sf.pymEo06a(DATE1, DATE2))
print(sf.pymBpn2xy(RBPN))
```

```text
-0.0013328827151307448
-0.00538297699189191
(0.001093465510215479, -4.281337229063151e-05)
```

```python
import PyMsOfa as sf

RC2I = [[0.9999998323037164738, 0.5581526271714303683e-9, -0.5791308477073443903e-3],
        [-0.2384266227524722273e-7, 0.9999999991917404296, -0.4020594955030704125e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
    #celestial-to-intermediate matrix
ERA = 1.75283325530307        #Earth rotation angle (radians)
RPOM = [[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
        [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
        [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]
    #polar-motion matrix
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7    #coordinates of the pole
YP = 1.860359247e-6    #coordinates of the pole
RBPN = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
    #celestial-to-true matrix
GST = 1.754166138040730516        #Greenwich (apparent) Sidereal Time (radians)
DEPS = 0.4090789763356509900        #nutation
DPSI = -0.9630909107115582393e-5    #nutation
X = 0.5791308486706011000e-3        #coordinates of the pole (radians)
Y = 0.4020579816732961219e-4        #coordinates of the pole (radians)

print(sf.pymC2tcio(RC2I, ERA, RPOM))
print(sf.pymC2t00a(TTA, TTB, UTA, UTB, XP, YP))
print(sf.pymC2t00b(TTA, TTB, UTA, UTB, XP, YP))
print(sf.pymC2t06a(TTA, TTB, UTA, UTB, XP, YP))
print(sf.pymC2teqx(RBPN, GST, RPOM))
print(sf.pymC2tpe(TTA, TTB, UTA, UTB, DPSI, DEPS, XP, YP))
print(sf.pymC2txy(TTA, TTB, UTA, UTB, X, Y, XP, YP))
```

```text
[[-1.81033213e-01  9.83476981e-01  6.55553564e-05]
 [-9.83476813e-01 -1.81033220e-01  5.74980112e-04]
 [ 5.77347401e-04  3.96183239e-05  9.99999833e-01]]
[[-1.81033213e-01  9.83476981e-01  6.55553564e-05]
 [-9.83476813e-01 -1.81033220e-01  5.74980112e-04]
 [ 5.77347401e-04  3.96183239e-05  9.99999833e-01]]
[[-1.81033213e-01  9.83476981e-01  6.55556508e-05]
 [-9.83476813e-01 -1.81033220e-01  5.74979392e-04]
 [ 5.77346747e-04  3.96179041e-05  9.99999833e-01]]
[[-1.81033213e-01  9.83476981e-01  6.55555096e-05]
 [-9.83476813e-01 -1.81033220e-01  5.74980084e-04]
 [ 5.77347402e-04  3.96181683e-05  9.99999833e-01]]
[[-1.81033213e-01  9.83476981e-01  6.55553564e-05]
 [-9.83476813e-01 -1.81033220e-01  5.74980112e-04]
 [ 5.77347401e-04  3.96183239e-05  9.99999833e-01]]
[[-0.1813678   0.90234822 -0.39099029]
 [-0.98341476 -0.16598836  0.07309764]
 [ 0.00105969  0.39776319  0.91748751]]
[[-1.81033213e-01  9.83476981e-01  6.55555125e-05]
 [-9.83476813e-01 -1.81033220e-01  5.74980084e-04]
 [ 5.77347403e-04  3.96181655e-05  9.99999833e-01]]
```

### 5.1 Light deflection by solar-system bodies

`pymLd(bm, p, q, e, em, dlim)` : Apply light deflection by a solar-system body, as part of transforming coordinate direction into natural direction.

`pymLdn(n, bodies, ob, sc)` : For a star, apply light deflection by multiple solar-system bodies, as part of transforming coordinate direction into natural direction.

`pymLdsun(p, e, em)` : Deflection of starlight by the Sun.

```python
import PyMsOfa as sf

BM = 0.00028574                    #mass of the gravitating body (solar masses)
P = [-0.763276255, -0.608633767, -0.216735543]    #direction from observer to source (unit vector)
Q = [-0.763276255, -0.608633767, -0.216735543]    #direction from body to source (unit vector)
E = [0.76700421, 0.605629598, 0.211937094]    #direction from body to observer (unit vector)
EM = 8.91276983                    #distance from body to observer (au)
DLIM = 3e-10                    #deflection limiter
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]    #data for each of the n bodies
OB = [-0.974170437, -0.2115201, -0.0917583114]    #barycentric position of the observer (au)
SC = [-0.763276255, -0.608633767, -0.216735543]    #observer to star coord direction (unit vector)

print(sf.pymLd(BM, P, Q, E, EM, DLIM))
print(sf.pymLdn(N, [sf.pymLDBODY(b[0], b[1], [[b[2], b[3], b[4]], [b[5], b[6], b[7]]]) for b in B], OB, SC))
print(sf.pymLdsun(P, E, EM))
```

```text
[-0.76327625 -0.60863377 -0.21673554]
[-0.76327626 -0.60863376 -0.21673554]
[-0.76327589 -0.60863406 -0.21673601]
```

### 5.2 Transformations based on ASTROM

`pymApcg(date1, date2, ebpv, ehp)` : For a geocentric observer, prepare star-independent astrometry parameters for transformations between ICRS and GCRS coordinates. The Earth ephemeris is supplied by the caller.

`pymApcg13(date1, date2)` : For a geocentric observer, prepare star-independent astrometry parameters for transformations between ICRS and GCRS coordinates. The caller supplies the date, and SOFA models are used to predict the Earth ephemeris.

`pymApci(date1, date2, ebpv, ehp, x, y, s)` : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and geocentric CIRS coordinates.  The Earth ephemeris and CIP/CIO are supplied by the caller.

`pymApci13(date1, date2)` : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and geocentric CIRS coordinates.  The caller supplies the date, and SOFA models are used to predict the Earth ephemeris and CIP/CIO.

`pymApco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb)` : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and observed coordinates.  The caller supplies the Earth ephemeris, the Earth rotation information and the refraction constants as well as the site coordinates.

`pymApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)` : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and observed coordinates.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength, and SOFA models are used to obtain the Earth ephemeris, CIP/CIO and refraction constants.

`pymApcs(date1, date2, pv, ebpv, ehp)` : For an observer whose geocentric position and velocity are known, prepare star-independent astrometry parameters for transformations between ICRS and GCRS.  The Earth ephemeris is supplied by the caller.

`pymApcs13(date1, date2, pv)` : For an observer whose geocentric position and velocity are known, prepare star-independent astrometry parameters for transformations between ICRS and GCRS.  The Earth ephemeris is from SOFA models.

`pymAper(theta, astrom)` : In the star-independent astrometry parameters, update only the Earth rotation angle, supplied by the caller explicitly.

`pymAper13(ut11, ut12, astrom)` : In the star-independent astrometry parameters, update only the Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

`pymApio(sp, theta, elong, phi, hm, xp, yp, refa, refb)` : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between CIRS and observed coordinates.  The caller supplies the Earth orientation information and the refraction constants as well as the site coordinates.

`pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)` : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between CIRS and observed coordinates.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

```python
import PyMsOfa as sf

DATE1 = 2456165.5
DATE2 = 0.401182685
EBPV = [[0.901310875, -0.417402664, -0.180982288],
        [0.00742727954, 0.0140507459, 0.00609045792]]
EHP = [0.903358544, -0.415395237, -0.180084014]
X, Y, S = 0.0013122272, -2.92808623e-5, 3.05749468e-8
THETA = 3.14540971
ELONG, PHI, HM = -0.527800806, -1.2345856, 2738.0
XP, YP = 2.47230737e-7, 1.82640464e-6
SP = -3.01974337e-11
REFA, REFB = 0.000201418779, -2.36140831e-7
UTC1, UTC2 = 2456384.5, 0.969254051
DUT1 = 0.1550675
PHPA, TC, RH, WL = 731.0, 12.8, 0.59, 0.55
PV = [[-1836024.09, 1056607.72, -5998795.26],
      [-77.0361767, -133.310856, 0.0971855934]]
UT11, UT12 = 2456165.5, 0.401182685

a = sf.pymApcg(DATE1, DATE2, EBPV, EHP)
print(a.pmt, a.em)
a = sf.pymApcg13(DATE1, DATE2)
print(a.pmt, a.em)
a = sf.pymApci(DATE1, DATE2, EBPV, EHP, X, Y, S)
print(a.pmt)
a, eo = sf.pymApci13(DATE1, DATE2)
print(a.pmt, eo)
print(sf.pymApco(DATE1, DATE2, EBPV, EHP, X, Y, S, THETA, ELONG, PHI, HM, XP, YP, SP, REFA, REFB).along)
a13, eo13 = sf.pymApco13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL)
print(a13.along, eo13)
a = sf.pymApcs(DATE1, DATE2, PV, EBPV, EHP)
print(a.pmt, a.em)
a = sf.pymApcs13(DATE1, DATE2, PV)
print(a.pmt, a.em)
a = sf.pymAper(THETA, sf.pymASTROM())
a.along = 1.234
a = sf.pymAper(THETA, a)
print(a.eral)
a = sf.pymAper13(UT11, UT12, sf.pymASTROM())
print(a.eral)
a = sf.pymApio(SP, THETA, ELONG, PHI, HM, XP, YP, REFA, REFB)
print(a.along, a.eral)
a = sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL)
print(a.along, a.eral)
```

```text
12.651337940273786 1.0104652958110132
12.651337940273786 1.0104652959646598
12.651337940273786
12.651337940273786 -0.0029006187126573756
-0.5278008060296
-0.5278008060296 -0.003020548354802412
12.651337940273786 1.0104585674728983
12.651337940273786 1.0104585676265498
4.37940971
2.0822366617896932
-0.5278008060296 2.6176089039704
-0.5278008060296 2.6176089091896517
```

### 5.2.2 Applying ASTROM parameters

`pymAtccq(rc, dc, pr, pd, px, rv, astrom)` : Quick transformation of a star's ICRS catalog entry (epoch J2000.0) into ICRS astrometric place, given precomputed star-independent astrometry parameters.

`pymAtcc13(rc, dc, pr, pd, px, rv, date1, date2)` : Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS astrometric place.

`pymAtciq(rc, dc, pr, pd, px, rv, astrom)` : Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed star-independent astrometry parameters.

`pymAtci13(rc, dc, pr, pd, px, rv, date1, date2)` : Transform ICRS star data, epoch J2000.0, to CIRS.

`pymAtciqn(rc, dc, pr, pd, px, rv, astrom, n, bodies)` : Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed star-independent astrometry parameters plus a list of light-deflecting bodies.

`pymAtciqz(rc, dc, astrom)` : Quick ICRS to CIRS transformation, given precomputed star-independent astrometry parameters, and assuming zero parallax and proper motion.

`pymAticq(ri, di, astrom)` : Quick CIRS RA,Dec to ICRS astrometric place, given the star-independent astrometry parameters.

`pymAtic13(ri, di, date1, date2)` : Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

`pymAticqn(ri, di, astrom, n, bodies)` : Quick CIRS to ICRS astrometric place transformation, given the star-independent astrometry parameters plus a list of light-deflecting bodies.

`pymAtoc13(coord_type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)` : Observed place at a groundbased site to to ICRS astrometric RA,Dec. The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

`pymAtco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)` : ICRS RA,Dec to observed place.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

`pymAtoiq(coord_type, ob1, ob2, astrom)` : Quick observed place to CIRS, given the star-independent astrometry parameters.

`pymAtoi13(coord_type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)` : Observed place to CIRS.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

`pymAtioq(ri, di, astrom)` : Quick CIRS to observed place transformation.

`pymAtio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)` : CIRS RA,Dec to observed place.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

```python
import PyMsOfa as sf

RC, DC = 2.71, 0.174
PR, PD, PX, RV = 1e-5, 5e-6, 0.1, 55.0
DATE1, DATE2 = 2456165.5, 0.401182685
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98019819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692,
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303,
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
bodies = [sf.pymLDBODY(b[0], b[1], [[b[2], b[3], b[4]], [b[5], b[6], b[7]]])
          for b in B]
RI, DI = 2.710121572969038991, 0.1729371367218230438
UTC1, UTC2 = 2456384.5, 0.969254051
DUT1 = 0.1550675
ELONG, PHI, HM = -0.527800806, -1.2345856, 2738.0
XP, YP = 2.47230737e-7, 1.82640464e-6
PHPA, TC, RH, WL = 731.0, 12.8, 0.59, 0.55
OB1, OB2 = 2.710085107986886201, 0.1717653435758265198
TYPE = "R"

a, _ = sf.pymApci13(DATE1, DATE2)
print(sf.pymAtccq(RC, DC, PR, PD, PX, RV, a))
print(sf.pymAtcc13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))
print(sf.pymAtciq(RC, DC, PR, PD, PX, RV, a))
print(sf.pymAtci13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))
print(sf.pymAtciqn(RC, DC, PR, PD, PX, RV, a, len(bodies), bodies))
print(sf.pymAtciqz(RC, DC, a))
print(sf.pymAticq(RI, DI, a))
print(sf.pymAtic13(RI, DI, DATE1, DATE2))
print(sf.pymAticqn(RI, DI, a, len(bodies), bodies))
print(sf.pymAtoc13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
print(sf.pymAtco13(RC, DC, PR, PD, PX, RV, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
b = sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL)
print(sf.pymAtoiq(TYPE, OB1, OB2, b))
print(sf.pymAtoi13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
print(sf.pymAtioq(RI, DI, b))
print(sf.pymAtio13(RI, DI, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
```

```text
(2.710126504531372, 0.17406325376283505)
(2.710126504531372, 0.17406325376283505)
(2.7101215729686965, 0.17293713672195385)
(2.7101215729686965, 0.17293713672195385, -0.0029006187126573756)
(2.710122008104983, 0.17293719164927662)
(2.7099948992472567, 0.17287407209849318)
(2.7101265045317167, 0.1740632537627034)
(2.7101265045317167, 0.1740632537627034, -0.0029006187126573756)
(2.710126069960012, 0.17406319850087898)
(2.709956744659136, 0.17416965008984703)
(0.09251774485486736, 1.4076614052564997, -0.09265154431530948, 0.17166265600725286, 2.710260453504961, -0.003020548354802412)
(2.710121574447541, 0.1729371839116607)
(2.710121574447541, 0.1729371839116607)
(0.09233952224896334, 1.4077587045135502, -0.09247619879882912, 0.17176534357562356, 2.7100851079884807)
(0.09233952224896334, 1.4077587045135502, -0.09247619879882912, 0.17176534357562356, 2.7100851079884807)
```
