# PyMsOfa

This package is a Python package for the Standards of Fundamental Astronomy (SOFA) service of the International Astronomical Union (IAU). It implements the python package PyMsOfa for SOFA service in three ways:

(1) a python wrapper package based on a foreign function library for Python (ctypes),

(2) a python wrapper package with the foreign function interface for Python calling C code (cffi),

(3) a python package directly written in pure python codes from SOFA subroutines.

It implements all 247 functions in the SOFA service and is based on the latest version released on Oct 11, 2023.

This Python package can be suitable for the astrometric detection of habitable planets of the Closeby Habitable Exoplanet Survey ([CHES](https://doi.org/10.1088/1674-4527/ac77e4)) mission and for the frontier themes of black holes and dark matter related to astrometric calculations and other fields.

### Citation

To cite PyMsOfa in publications use:

> 1. Ji, Jiang-Hui, Tan, Dong-jie, Bao, Chun-hui, Huang, Xiu-min, Hu, Shoucun, Dong, Yao, Wang, Su. 2023, PyMsOfa: A Python Package for the Standards of Fundamental Astronomy (SOFA) Service, Research in Astronomy and Astrophysics, https://doi.org/10.1088/1674-4527/ad0499
> 2. Ji, Jiang-Hui, Li, Hai-Tao, Zhang, Jun-Bo, Fang, Liang, Li, Dong, Wang, Su, Cao, Yang, Deng, Lei, Li, Bao-Quan, Xian, Hao, Gao, Xiao-Dong, Zhang, Ang, Li, Fei, Liu, Jia-Cheng, Qi, Zhao-Xiang,  Jin, Sheng, Liu, Ya-Ning, Chen, Guo, Li, Ming-Tao, Dong, Yao, Zhu, Zi, and CHES Consortium. 2022, CHES: A Space-borne Astrometric Mission for the Detection of Habitable Planets of the Nearby Solar-type Stars, Research in Astronomy and Astrophysics, 22, 072003

A BibTeX entry for LaTeX users is

```bibtex
@ARTICLE{2023arXiv231008673J,
       author = {{Ji}, Jianghui and {Tan}, Dongjie and {Bao}, Chunhui and {Huang}, Xiumin and {Hu}, Shoucun and {Dong}, Yao and {Wang}, Su},
        title = "{PyMsOfa: A Python Package for the Standards of Fundamental Astronomy (SOFA) Service}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
         year = 2023,
        month = oct,
          eid = {arXiv:2310.08673},
        pages = {arXiv:2310.08673},
          doi = {10.48550/arXiv.2310.08673},
archivePrefix = {arXiv},
       eprint = {2310.08673},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv231008673J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@article{2022RAA....22g2003J,
       author = {{Ji}, Jiang-Hui and {Li}, Hai-Tao and {Zhang}, Jun-Bo and {Fang}, Liang and {Li}, Dong and {Wang}, Su and {Cao}, Yang and {Deng}, Lei and {Li}, Bao-Quan and {Xian}, Hao and {Gao}, Xiao-Dong and {Zhang}, Ang and {Li}, Fei and {Liu}, Jia-Cheng and {Qi}, Zhao-Xiang and {Jin}, Sheng and {Liu}, Ya-Ning and {Chen}, Guo and {Li}, Ming-Tao and {Dong}, Yao and {Zhu}, Zi and {CHES Consortium}},
        title = "{CHES: A Space-borne Astrometric Mission for the Detection of Habitable Planets of the Nearby Solar-type Stars}",
      journal = {Research in Astronomy and Astrophysics},
     keywords = {Astrometry and Celestial Mechanics, planets and satellites: detection, planets and satellites: terrestrial planets, stars: solar-type, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Astrophysics of Galaxies, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2022,
        month = jul,
       volume = {22},
       number = {7},
          eid = {072003},
        pages = {072003},
          doi = {10.1088/1674-4527/ac77e4},
archivePrefix = {arXiv},
       eprint = {2205.05645},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022RAA....22g2003J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

In this paper, the function and usage are briefly introduced, which can be used as a reference when using PyM-sOfa. For detailed features of SOFA, please refer to : [http://www.iausofa.org/cookbooks.html](http://www.iausofa.org/cookbooks.html)

---

### Contents

* Section 1 Basic
  * 1.1 Copy a paramater
  * 1.2 Initialize a paramater
  * 1.3 Normalize angle
  * 1.4 Paramater processing
  * 1.5 Parameter Operation
  * 1.6 Matrix rotation
  * 1.7 Coordinate transformation
  * 1.8 Conversion between rotation vector and rotation matrix
  * 1.9 Position relation
  * 1.10 Projection relationship
  * 1.11 Others
* Section 2 Time
  * 2.1 Sidereal time
  * 2.2 Time unit conversion
  * 2.3 Date conversion
  * 2.4 Time conversion
* Section 3 Precession, nutation, and polar shift
  * 3.1 Parameters of celestial bodies in the solar system
  * 3.2 Bias-precession-nutation matrix
* Section 4 Transformation of the coordinate system
* Section 5 Astrometric parameter
  * 5.1 Deflection of light by celestial bodies in the solar system
  * 5.2 Coordinate transformation based on astrometric parameter ASTROM
    * 5.2.1 Generate ASTROM parameters
  * 5.2.2 ASTROM parameters application

---

### Section 1   Basic

#### 1.1 Copy a paramater

pymCp(A) : Copy a p-vector.(1*3)

pymCpv(A) : Copy a position/velocity vector.(2*3)

pymCr(A) : Copy an r-matrix.(3*3)

```python
import PyMsOfa as sf

A = [3,4,5]
B = [[3,4,5],[4,5,6]]
C = [[3,4,5],[4,5,6],[5,6,7]]

print(sf.pymCp(A))
print(sf.pymCpv(B))
print(sf.pymCr(C))

>>> [3,4,5]
>>> [[3,4,5],
     [4,5,6]]
>>> [[3,4,5],
     [4,5,6],
     [5,6,7]]
```

#### 1.2 Initialize a paramater

pymZp(A) : Zero a p-vector.(1*3)

pymZpv(A) : Zero a pv-vector.(2*3)

pymZr(A) : Initialize an r-matrix to the null matrix.(3*3)

pymIr(A) : Initialize an r-matrix to the identity matrix.(3*3)

```python
import PyMsOfa as sf

A = [3,4,5]
B = [[3,4,5],[4,5,6]]
C = [[3,4,5],[4,5,6],[5,6,7]]

print(sf.pymZp(A))
print(sf.pymZpv(B))
print(sf.pymZr(C))
print(sf.pymIr(C))

>>> [0,0,0]
>>> [[0,0,0],
     [0,0,0]]
>>> [[0,0,0],
     [0,0,0],
     [0,0,0]]
>>> [[1,1,1],
     [1,1,1],
     [1,1,1]]
```

#### 1.3 Normalize angle

pymAnp(A) : Normalize angle into the range 0 <= a < 2pi.

pymAnpm(A) : Normalize angle into the range -pi <= a < +pi.

```python
import PyMsOfa as sf
import math as ma

A = 7/2*ma.pi

print(sf.pymAnp(A))
print(sf.pymAnpm(A))

>>> 4.712	#3/2 pi
>>> -1.571	#-1/2 pi
```

#### 1.4 Paramater processing

pymPm(A) : Modulus of p-vector.

pymPn(A) : Convert a p-vector into modulus and unit vector.

pymTr(A) : Transpose an r-matrix.

pymPv2p(A) : Discard velocity component of a pv-vector.

pymPvm(A) : Modulus of pv-vector.

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

>>> 7.071
>>> (7.071, [0.424, 0.566, 0.707])
>>> [[3, 3, 3], 
     [4, 4, 4], 
     [5, 5, 5]]
>>> [3, 4, 5]
>>> (7.071, 8.77)
```

#### 1.5 Parameter Operation

pymPpp(A, B) : P-vector addition.

pymPmp(A, B) : P-vector subtraction.

pymSxp(N, A) : Multiply a p-vector by a scalar.

pymPpsp(A, N, B) : P-vector plus scaled p-vector.

pymPdp(A, B) : p-vector inner (=scalar=dot) product.

pymPxp(A, B) : p-vector outer (=vector=cross) product.

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

>>> [7, 9, 11]
>>> [-1, -1, -1]
>>> [6, 8, 10]
>>> [11, 14, 17]
>>> 62.0
>>> [-1, 2, -1]
```

pymRxp(A, B) : Multiply a p-vector by an r-matrix.

pymTrxp(A, B) : Multiply a p-vector by the transpose of an r-matrix.

pymRxpv(A, B) : Multiply a pv-vector by an r-matrix.

pymTrxpv(A, B) : Multiply a pv-vector by the transpose of an r-matrix.

pymRxr(A, B) : Multiply two r-matrices.

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

>>> [50.0, 50.0, 50.0]
>>> [36.0, 48.0, 60.0]
>>> [[50.0, 50.0, 50.0], 
     [62.0, 62.0, 62.0]]
>>> [[36.0, 48.0, 60.0], 
     [45.0, 60.0, 75.0]]
>>> [[36.0, 48.0, 60.0], 
     [36.0, 48.0, 60.0], 
     [36.0, 48.0, 60.0]]
```

pymPvppv(A, B) : Add one pv-vector to another.

pymPvmpv(A, B) : Subtract one pv-vector from another.

pymPvdpv(A, B) : Inner (=scalar=dot) product of two pv-vectors.

pymPvxpv(A, B) : Outer (=vector=cross) product of two pv-vectors.

pymSxpv(N, A) : Multiply a pv-vector by a scalar.

pymS2xpv(N, M, A) : Multiply a pv-vector by two scalars.

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

>>> [[6, 8, 10], [7, 9, 11]]
>>> [[0, 0, 0], [-1, -1, -1]]
>>> [50.0, 112.0]
>>> [[0, 0, 0], [-1, 2, -1]]
>>> [[6, 8, 10], [6, 8, 10]]
>>> [[6, 8, 10], [9, 12, 15]]
```

#### 1.6 Matrix rotation

pymRx(PHI, A) : Rotate an r-matrix about the x-axis.

pymRy(PHI, A) : Rotate an r-matrix about the y-axis.

pymRz(PHI, A) : Rotate an r-matrix about the z-axis.

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

>>>[[1, 1, 1],
    [2.83, 2.83, 2.83], 
    [0.0, 0.0, 0.0]]
>>>[[0.0, 0.0, 0.0],
    [1, 1, 1],
    [2.83, 2.83, 2.83]]
>>>[[2.83, 2.83, 2.83],
    [0.0, 0.0, 0.0],
    [1, 1, 1]]
```

#### 1.7 Coordinate transformation

pymC2s(A) : P-vector to spherical coordinates.

pymS2c(THETA, PHI) : Convert spherical coordinates to Cartesian.

```python
import PyMsOfa as sf
import math as ma

A = [0.5, 0.5, ma.sqrt(2)/2]
THETA = ma.pi/4
PHI = ma.pi/4

print(sf.pymC2s(A))
print(sf.pymS2c(THETA, PHI))

>>>(0.785, 0.785)
>>>[0.5, 0.5, 0.707]
```

pymP2s(A) : P-vector to spherical polar coordinates (including radial distance).

pymS2p(THETA, PHI, R) : Convert spherical polar coordinates (including radial distance) to p-vector.

```python
import PyMsOfa as sf
import math as ma

A = [0.5, 0.5, ma.sqrt(2)/2]
THETA = ma.pi/4
PHI = ma.pi/4
R = 1.0

print(sf.pymP2s(A))
print(sf.pymS2p(THETA, PHI, R))

>>>(0.785, 0.785, 1.0)
>>>[0.5, 0.5, 0.707]
```

pymPv2s(PV) : Convert position/velocity from Cartesian to spherical coordinates (including latitude and longitude, radial distance, and the rate of change of all three).

pymS2pv(THETA, PHI, R, TD, PD, RD) : Convert position/velocity from spherical to Cartesian coordinates.

```python
import PyMsOfa as sf
import math as ma

PV = [[0.5, 0.5, ma.sqrt(2)/2],
      [1, 1, 1]]
THETA = ma.pi/4		#longitude angle (radians)
PHI = ma.pi/4		#latitude angle (radians)
R = 1.0			#radial distance
TD = 0.0		#rate of change of theta
PD = -0.293		#rate of change of phi
RD = 1.707		#rate of change of r

print(sf.pymPv2s(PV))
print(sf.pymS2pv(THETA, PHI, R, TD, PD, RD))

>>>(0.785, 0.785, 1.0, 0.0, -0.293, 1.707)
>>>[[0.5, 0.5, 0.707],
    [1.0, 1.0, 1.0]]
```

pymPvstar(PV) : Convert star position+velocity vector to catalog coordinates (including RA, DEC, pmRA, pmDEC, plx and RV).

pymStarpv(RA, DEC, PMRA, PMDEC, PLX, RV) : Convert star catalog coordinates to position+velocity vector.

```python
import PyMsOfa as sf
import math as ma

PV = [[0.5, 0.5, ma.sqrt(2)/2],
      [1, 1, ma.sqrt(2)]]
RA = ma.pi/4
DEC = ma.pi/4
PMRA = 0.0		#radians/year
PMDEC = -105.935	#radians/year
PLX = 206264.806	#arcsce
RV = 2941.778		#km/s,positive = receding

print(sf.pymPvstar(PV))
print(sf.pymStarpv(RA, DEC, PMRA, PMDEC ,PLX, RV))

>>>(0.785, 0.785, 0.0, -105.935, 206264.806, 2941.778, 0)#The last value is the check value, where 0 is no error.
>>>([[0.5, 0.5, 0.707], 
     [1.0, 1.0, 1.0]], 0)  				 #The last value is the check value, where 0 is no error.
```

#### 1.8 Conversion between rotation vector and rotation matrix

Rotation vector: The direction of the vector is the axis of rotation, and the magnitude of the vector is the Angle of rotation, which is rotated clockwise along the axis of rotation.

pymRv2m(W) : Form the r-matrix corresponding to a given r-vector.

pymRm2v(R) : Express an r-matrix as an r-vector.

```python
import PyMsOfa as sf
import math as ma

W = [1,2,3]
R = [[0.00, -0.80, -0.60],
     [0.80, -0.36, 0.48],
     [0.60, 0.48, -0.64]]

print(sf.pymRv2m(W))
print(sf.pymRm2v(R))

>>>[[-0.695, -0.192, 0.693], 
    [0.714, -0.304, 0.631], 
    [0.089, 0.933, 0.348]]
>>>[0.0, 1.414, -1.885]
```

#### 1.9 Position relation

pymPas(AL, AP, BL, BP) : Position-angle from spherical coordinates.

pymPap(A, B) : Position-angle from two p-vectors.

pymSeps(AL, AP, BL, BP) : Angular separation between two sets of spherical coordinates.

pymSepp(A,B) : Angular separation between two p-vectors.

```python
import PyMsOfa as sf
import math as ma

AL = 0.5 * ma.pi		#longitude of point A (e.g. RA) in radians 
AP = 0.8 * ma.pi		#latitude of point A (e.g. Dec) in radians  
BL = 0.4 * ma.pi		#longitude of point B  
BP = 0.6 * ma.pi		#latitude of point B 
A = [1, 2, 3]
B = [2, 3, 4]

print(sf.pymPas(AL, AP, BL, BP))
print(sf.pymPap(A, B))
print(sf.pymSeps(AL, AP, BL, BP))
print(sf.pymSepp(A, B))

>>>2.983
>>>-2.390
>>>0.649
>>>0.122
```

#### 1.10 Projection relationship

pymTpors(XI, ETA, A, B) : In the tangent plane projection, given the rectangular coordinates of a star and its spherical coordinates, determine the spherical coordinates of the tangent point.

pymTporv(XI, ETA, V) : In the tangent plane projection, given the rectangular coordinates of a star and its direction cosines, determine the direction cosines of the tangent point.

pymTpsts(XI, ETA, A0, B0) : In the tangent plane projection, given the star's rectangular coordinates and the spherical coordinates of the tangent point, solve for the spherical coordinates of the star.

pymTpstv(XI, ETA, V0) : In the tangent plane projection, given the star's rectangular coordinates and the direction cosines of the tangent point, solve for the direction cosines of the star.

pymTpxes(A, B, A0, B0) : In the tangent plane projection, given celestial spherical coordinates for a star and the tangent point, solve for the star's rectangular coordinates in the tangent plane.

pymTpxev(V, V0) : In the tangent plane projection, given celestial direction cosines for a star and the tangent point, solve for the star's rectangular coordinates in the tangent plane.

```python
import PyMsOfa as sf
import math as ma

XI = 1					#rectangular coordinates of star image  
ETA = 2					#rectangular coordinates of star image  
A = ma.pi/4				#star's spherical coordinates  
B = ma.pi/6				#star's spherical coordinates  
V = [0.2, 0.4, ma.sqrt(0.8)]		#star's direction cosines
A0 = ma.pi/4+ma.pi/180			#tangent point's spherical coordinates  
B0 = ma.pi/6-ma.pi/180			#tangent point's spherical coordinates  
V0 = [0.2, 0.401, ma.sqrt(0.7992)]	#tangent point's direction cosines

print(sf.pymTpors(XI, ETA, A, B))
print(sf.pymTporv(XI, ETA, V))
print(sf.pymTpsts(XI, ETA, A0, B0))
print(sf.pymTpstv(XI, ETA, V0))
print(sf.pymTpxes(A, B, A0, B0))
print(sf.pymTpxev(V, V0))

>>>(0.295, -0.528, 
    4.418, 1.455, 2)		#The last value is the number of solutions, where 2 = both solutions are useful.
>>>([0.965, -0.042, 0.260], 
    [0.499, -0.609, 0.617], 2)	#The last value is the number of solutions, where 2 = both solutions are useful.
>>>(2.468, 1.148)
>>>[-0.609, -0.307, 0.731]
>>>(-0.015, 0.018, 0)		#The last value is the check value, where 0 is OK.
>>>(-4.463e-4, 1.000e-3, 0)	#The last value is the check value, where 0 is OK.
```

#### 1.11 Others

pymP2pv(P) : Extend a p-vector to a pv-vector by appending a zero velocity.

pymPvu(DT, PV) : Update a pv-vector.

pymPvup(DT, PV) : Update a pv-vector, discarding the velocity component.

```python
import PyMsOfa as sf

P = [1, 2, 3]
DT = 10				#time interval  
PV = [[1, 2, 3],
      [4, 5, 6]]

print(sf.pymP2pv(P))
print(sf.pymPvu(DT, PV))
print(sf.pymPvup(DT, PV))

>>>[[1, 2, 3], 
    [0, 0, 0]]
>>>[[41, 52, 63], 
    [4, 5, 6]]
>>>[41, 52, 63]
```

### Section 2 Time

#### 2.1 Sidereal time

pymGmst82(DJ1, DJ2) : Universal Time to Greenwich mean sidereal time (IAU 1982 model).

pymGst94(UTA, UTB) : Greenwich apparent sidereal time (consistent with IAU 1982/94 resolutions).

pymGmst00(UTA, UTB, TTA, TTB) : Greenwich mean sidereal time (model consistent with IAU 2000 resolutions).

pymGst00a(UTA, UTB, TTA, TTB) : Greenwich apparent sidereal time (consistent with IAU 2000 resolutions).

pymGst00b(UTA, UTB) : Greenwich apparent sidereal time (consistent with IAU 2000 resolutions but using the truncated nutation model IAU 2000B).

pymGmst06(UTA, UTB, TTA, TTB) : Greenwich mean sidereal time (consistent with IAU 2006 precession).

pymGst06(UTA, UTB, TTA, TTB, RNPB) : Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

pymGst06a(UTA, UTB, TTA, TTB) : Greenwich apparent sidereal time (consistent with IAU 2000 and 2006 resolutions).

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

>>>3.3060549204240584
>>>0.4633455111601113
>>>0.4634280218387905
>>>0.46334550139240926
>>>0.4633455023818729
>>>0.46342802175551745
>>>0.4643595506102784
>>>0.4633455012681201
```

#### 2.2 Time unit conversion

pymD2tf(NDP, DAYS) : Decompose days to hours, minutes, seconds, fraction.

pymTf2d(S, IHOUR, IMIN, SEC) : Convert hours, minutes, seconds to days.

pymA2af(NDP, ANGLE) : Decompose radians into degrees, arcminutes, arcseconds, fraction.

pymAf2a(S, IDEG, IAMIN, ASEC) : Convert degrees, arcminutes, arcseconds to radians.

pymTf2a(S, IHOUR, IMIN, SEC) : Convert hours, minutes, seconds to radians.

```python
import PyMsOfa as sf
import math as ma

NDP = 3		#resolution
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
print(sf.pymAf2a(S, IDEG, IAMIN, ASEC))
print(sf.pymTf2a(S, IHOUR, IMIN, SEC))

>>>('+', [12, 0, 0, 0])
>>>(0.22558217592592592, 0)
>>>('+', [45, 0, 0, 0])
>>>(4.792588701805817, 0)
>>>(1.4173746133393783, 0)
```

#### 2.3 Date conversion

pymCal2jd(IY, IM, ID) : Gregorian Calendar to Julian Date.

pymJd2cal(DJ1, DJ2) : Julian Date to Gregorian year, month, day, and fraction of a day.

pymJdcalf(NDP, DJ1, DJ2) : Julian Date to Gregorian Calendar, expressed in a form convenient for formatting messages:  rounded to a specified precision.

pymEpj(DJ1, DJ2) : Julian Date to Julian Epoch

pymEpj2jd(EPJ) : Julian Epoch to Julian Date

pymEpb(DJ1, DJ2) : Julian Date to Besselian Epoch.

pymEpb2jd(EPB) : Besselian Epoch to Julian Date

```python
import PyMsOfa as sf

IY = 2024
IM = 2
ID = 27
DJ1 = 2400000.5
DJ2 = 60367.0
NDP = 3				#resolution
EPJ = 2024.154688569473		#Julian Epoch
EPB = 2024.1564820038504	#Besselian Epoch

print(sf.pymCal2jd(IY, IM, ID))
print(sf.pymJd2cal(DJ1, DJ2))
print(sf.pymJdcalf(NDP, DJ1, DJ2))
print(sf.pymEpj(DJ1, DJ2))
print(sf.pymEpj2jd(EPJ))
print(sf.pymEpb(DJ1, DJ2))
print(sf.pymEpb2jd(EPB))

>>>(2400000.5, 60367.0, 0)	#The last value is the check value, where 0 is OK.
>>>(2024, 2, 27, 0.0, 0)	#The last value is the check value, where 0 is OK.
>>>([2024, 2, 27, 0], 0)	#The last value is the check value, where 0 is OK.
>>>2024.154688569473
>>>(2400000.5, 60366.99999999998)
>>>2024.1564820038504
>>>(2400000.5, 60367.0)
```

#### 2.4 Time conversion

pymDat(IY, IM, ID, FD) : For a given UTC date, calculate Delta(AT) = TAI-UTC.

pymDtf2d(SCALE, IY, IM, ID, IHR, IMN, SEC) : Encode date and time fields into 2-part Julian Date (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).

pymDtdb(DATE1, DATE2, UT, ELONG, U, V) : An approximation to TDB-TT, the difference between barycentric dynamical time and terrestrial time, for an observer on the Earth.

```python
import PyMsOfa as sf

IY = 2024
IM = 2
ID = 27
FD = 0.6
SCALE = '+'
IHR = 14
IMN = 5
SEC = 30
DATE1 = 2400000.5
DATE2 = 60367.0
UT = 0.76543		#universal time (UT1, fraction of one day)
ELONG = 5.0123		#longitude (east positive, radians) 
U = 5525.242		#distance from Earth spin axis (km)
V = 3190.0		#distance north of equatorial plane (km)

print(sf.pymDat(IY, IM, ID, FD))
print(sf.pymDtf2d(SCALE,IY,IM,ID,IHR,IMN,SEC))
print(sf.pymDtdb(DATE1, DATE2, UT, ELONG, U, V))

>>>(37.0, 0)			#The last value is the check value, where 0 is OK.
>>>(2460367.5, 0.5871528, 0)	#The last value is the check value, where 0 is OK.
>>>0.0013094390627437918
```

pymD2dtf(SCALE, NDP, D1, D2) : Format for output a 2-part Julian Date (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).

pymUtctai(UTC1, UTC2) : Time scale transformation:  Coordinated Universal Time, UTC, to International Atomic Time, TAI.

pymTaiutc(TAI1, TAI2) : Time scale transformation:  International Atomic Time, TAI, to Coordinated Universal Time, UTC.

pymTaiut1(TAI1, TAI2, DTA) : Time scale transformation:  International Atomic Time, TAI, to Universal Time, UT1.

pymUt1tai(UT11, UT12, DTA) : Time scale transformation:  Universal Time, UT1, to International Atomic Time, TAI.

pymUtcut1(UTC1, UTC2, DUT1) : Time scale transformation:  Coordinated Universal Time, UTC, to Universal Time, UT1.

pymUt1utc(UT11, UT12, DUT1) : Time scale transformation:  Universal Time, UT1, to Coordinated Universal Time, UTC.

pymTaitt(TAI1, TAI2) : Time scale transformation:  International Atomic Time, TAI, to Terrestrial Time, TT.

pymTttai(TT1, TT2) : Time scale transformation:  Terrestrial Time, TT, to International Atomic Time, TAI.

pymTttcg(TT1, TT2) : Time scale transformation:  Terrestrial Time, TT, to Geocentric Coordinate Time, TCG.

pymTcgtt(TCG1, TCG2) : Time scale transformation:  Geocentric Coordinate Time, TCG, to Terrestrial Time, TT.

pymTttdb(TT1, TT2, DTR) : Time scale transformation:  Terrestrial Time, TT, to Barycentric Dynamical Time, TDB.

pymTdbtt(TDB1, TDB2, DTR) : Time scale transformation:  Barycentric Dynamical Time, TDB, to Terrestrial Time, TT.

pymTtut1(TT1, TT2, DT) : Time scale transformation:  Terrestrial Time, TT, to Universal Time, UT1.

pymUt1tt(UT11, UT12, DT) : Time scale transformation:  Universal Time, UT1, to Terrestrial Time, TT.

pymTdbtcb(TDB1, TDB2) : Time scale transformation:  Barycentric Dynamical Time, TDB, to Barycentric Coordinate Time, TCB.

pymTcbtdb(TCB1, TCB2) : Time scale transformation:  Barycentric Coordinate Time, TCB, to Barycentric Dynamical Time, TDB.

```python
		import PyMsOfa as sf

SCALE = '+'
NDP = 3			#resolution
D1 = 2400000.5
D2 = 60367.0
UTC1 = 2400000.5
UTC2 = 60367.0
TAI1 = 2400000.5
TAI2 = 60367.0
DTA = -32.6659		#UT1-TAI in seconds
UT11 = 2400000.5
UT12 = 60367.0
DUT1 = 0.3341		#Delta UT1: UT1-UTC in seconds
TT1 = 2400000.5
TT2 = 60367.0
TCG1 = 2400000.5
TCG2 = 60367.0
DTR = -0.000201		#TDB-TT in seconds
TDB1 = 2400000.5
TDB2 = 60367.0
DT = 64.8499		#TT-UT1 in seconds
TCB1 = 2400000.5
TCB2 = 60367.0

print(sf.pymD2dtf(SCALE,NDP,D1,D2))
print(sf.pymUtctai(UTC1, UTC2))
print(sf.pymTaiutc(TAI1,TAI2))
print(sf.pymTaiut1(TAI1,TAI2,DTA))
print(sf.pymUt1tai(UT11, UT12, DTA))
print(sf.pymUtcut1(UTC1, UTC2, DUT1))
print(sf.pymUt1utc(UT11, UT12, DUT1))
print(sf.pymTaitt(TAI1,TAI2))
print(sf.pymTttai(TT1, TT2))
print(sf.pymTttcg(TT1, TT2))
print(sf.pymTcgtt(TCG1, TCG2))
print(sf.pymTttdb(TT1, TT2, DTR))
print(sf.pymTdbtt(TDB1, TDB2, DTR))
print(sf.pymTtut1(TT1, TT2, DT))
print(sf.pymUt1tt(UT11, UT12, DT))
print(sf.pymTdbtcb(TDB1, TDB2))
print(sf.pymTcbtdb(TCB1, TCB2))

#The last value is the check value, where 0 is OK.
>>>(2024, 2, 27, [0, 0, 0, 0], 0)
>>>(2400000.5, 60367.00042824074, 0)
>>>(2400000.5, 60366.99957175926, 0)
>>>(2400000.5, 60366.999621922456, 0)
>>>(2400000.5, 60367.000378077544, 0)
>>>(2400000.5, 60367.000003866895, 0)
>>>(2400000.5, 60366.999996133105, 0)
>>>(2400000.5, 60367.0003725, 0)
>>>(2400000.5, 60366.9996275, 0)
>>>(2400000.5, 60367.000012003205, 0)
>>>(2400000.5, 60366.999987996795, 0)
>>>(2400000.5, 60366.99999999767, 0)
>>>(2400000.5, 60367.00000000233, 0)
>>>(2400000.5, 60366.99924942246, 0)
>>>(2400000.5, 60367.00075057754, 0)
>>>(2400000.5, 60367.00026704677, 0)
>>>(2400000.5, 60366.99973295323, 0)
```

### Section 3 Precession, nutation, and polar shift

#### 3.1 Parameters of celestial bodies in the solar system

pymFapa03(T) : Fundamental argument, IERS Conventions (2003) : general accumulated precession in longitude.

pymFae03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Earth.

pymFad03(T) : Fundamental argument, IERS Conventions (2003) : mean elongation of the Moon from the Sun.

pymFaom03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of the Moon's ascending node.

pymFaf03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of the Moon minus mean longitude of the ascending node.

pymFal03(T) : Fundamental argument, IERS Conventions (2003) : mean anomaly of the Moon.

pymFalp03(T) : Fundamental argument, IERS Conventions (2003) : mean anomaly of the Sun.

pymFame03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Mercury.

pymFave03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Venus.

pymFama03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Mars.

pymFaju03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Jupiter.

pymFasa03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Saturn.

pymFaur03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Uranus.

pymFane03(T) : Fundamental argument, IERS Conventions (2003) : mean longitude of Neptune.

pymPlan94(DATE1, DATE2, NP) : Approximate heliocentric position and velocity of a nominated major planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or Neptune (but not the Earth itself).

pymMoon98(DATE1, DATE2) : Approximate geocentric position and velocity of the Moon.

```python
import PyMsOfa as sf

T = 0.8		#TDB, Julian centuries since J2000.0
DATE1 = 2400000.5
DATE2 = 60367.0
NP = 1		#planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars, 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)

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

>>>0.019508847622400002
>>>1.7447137389131058
>>>1.9467092053972077
>>>0.30956686622828467
>>>0.2597711366751447
>>>5.132369751109105
>>>6.226797973505531
>>>5.417338184297677
>>>3.4249004605338342
>>>3.275506840277835
>>>5.275711665202486
>>>5.371574539440829
>>>5.180636450180414
>>>2.079343830860413
>>>([[0.34673534057893274, -0.1407723578655547, -0.1111381624431485], 
     [0.007079322757371702, 0.023659120706492093, 0.011904985268464919]], 0)	#The last value is the check value, where 0 is OK.
>>>[[-0.002703439370794801, -0.00020716104811251522, -3.717027171462571e-05], 
    [4.7878512048071815e-05, -0.0004910931056274818, -0.0002672119831026995]]
```

#### 3.2 Bias-precession-nutation matrix

pymNumat(EPSA, DPSI, DEPS) : Form the matrix of nutation.

pymLtp(EPJ) : Long-term precession matrix.

pymLtpb(EPJ) : Long-term precession matrix, including ICRS frame bias.

pymLtpequ(EPJ) : Long-term precession of the equator.

pymLtpecl(EPJ) : Long-term precession of the ecliptic.

```python
import PyMsOfa as sf

EPSA =  0.4090789763356509900		#mean obliquity of date
DPSI = -0.9630909107115582393e-5	#nutation
DEPS =  0.4063239174001678826e-4	#nutation
EPJ = 1666.666				#Julian epoch (TT)

print(sf.pymNumat(EPSA, DPSI, DEPS))
print(sf.pymLtp(EPJ))
print(sf.pymLtpb(EPJ))
print(sf.pymLtpequ(EPJ))
print(sf.pymLtpecl(EPJ))

>>>[[0.9999999999536228, 8.83623932023625e-06, 3.830833447458252e-06], 
    [-8.83608365701669e-06, 0.9999999991354653, -4.0632408653651364e-05], 
    [-3.8311924818333855e-06, 4.063237480217419e-05, 0.999999999167166]]
>>>[[0.9967044141159213, 0.07437801893193208, 0.03237624409345603], 
    [-0.07437802731819615, 0.9972293894454531, -0.0012057688427235401], 
    [-0.03237622482766575, -0.0012062860396976614, 0.9994750246704011]]
>>>[[0.996704416772327, 0.07437794731203337, 0.03237632684841626], 
    [-0.07437795663437174, 0.9972293947500012, -0.00120574186591119], 
    [-0.03237630543224665, -0.0012063167910765376, 0.9994750220222438]]
>>>[-0.03237622482766575, -0.0012062860396976614, 0.9994750246704011]
>>>[-5.737091641342703e-05, -0.39847335458761085, 0.917179907320409]
```

pymPrec76(DATE01, DATE02, DATE11, DATE12) : IAU 1976 precession model. This function forms the three Euler angles which implement general precession between two dates, using the IAU 1976 model (as for the FK5 catalog).

pymPmat76(DATE1, DATE2) : Precession matrix from J2000.0 to a specified date, IAU 1976 model.

pymPnm80(DATE1, DATE2) : Form the matrix of precession/nutation for a given date, IAU 1976 precession model, IAU 1980 nutation model.

pymNut80(DATE1, DATE2) : Nutation, IAU 1980 model.

pymNutm80(DATE1, DATE2) : Form the matrix of nutation for a given date, IAU 1980 model.

pymObl80(DATE1, DATE2) : Mean obliquity of the ecliptic, IAU 1980 model.

pymEqeq94(DATE1, DATE2) : Equation of the equinoxes, IAU 1994 model.

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

>>>(0.005588961642000161, 0.005589922365870681, 0.004858945471687296)
>>>[[0.9999995504328985, 0.000869663159726508, 0.0003779153208913893], 
    [-0.0008696631597269232, 0.9999996218429094, -1.643284544746013e-07], 
    [-0.00037791532089043395, -1.6433065147768536e-07, 0.9999999285899891]]
>>>[[0.9999995831935007, 0.0008373653649393289, 0.0003639121744485501], 
    [-0.0008373804499809821, 0.9999996485440003, 4.130203432601264e-05], 
    [-0.00036387746165638676, -4.16067500514884e-05, 0.9999999329310333]]
>>>(3.520279437486982e-05, -4.145438583593602e-05)
>>>[[0.9999999993803816, -3.22978087631124e-05, -1.4003152516090194e-05], 
    [3.229838922744821e-05, 0.9999999986191833, 4.145415968437316e-05], 
    [1.4001813618232485e-05, -4.14546119379966e-05, 0.9999999990427321]]
>>>0.40910163117239295
>>>3.229357405996896e-05
```

pymBi00() : Frame bias components of IAU 2000 precession-nutation models;  part of the Mathews-Herring-Buffett (MHB2000) nutation series, with additions.

pymPr00(DATE1, DATE2) : Precession-rate part of the IAU 2000 precession-nutation models (part of MHB2000).

pymBp00(DATE1, DATE2) : Frame bias and precession, IAU 2000.

pymPmat00(DATE1, DATE2) : Precession matrix (including frame bias) from GCRS to a specified date, IAU 2000 model.

pymEe00(DATE1, DATE2, EPSA, DPSI) : The equation of the equinoxes, compatible with IAU 2000 resolutions, given the nutation in longitude and the mean obliquity.

pymEe00a(DATE1, DATE2) : Equation of the equinoxes, compatible with IAU 2000 resolutions.

pymEe00b(DATE1, DATE2) : Equation of the equinoxes, compatible with IAU 2000 resolutions but using the truncated nutation model IAU 2000B.

pymEect00(DATE1, DATE2) : Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 53736
EPSA = 0.4090789763356509900		#mean obliquity
DPSI = -0.9630909107115582393e-5	#nutation in longitude

print(sf.pymBi00())
print(sf.pymPr00(DATE1, DATE2))
print(sf.pymBp00(DATE1, DATE2))
print(sf.pymPmat00(DATE1, DATE2))
print(sf.pymEe00(DATE1, DATE2, EPSA, DPSI))
print(sf.pymEe00a(DATE1, DATE2))
print(sf.pymEect00(DATE1, DATE2))

>>>(-2.0253091528350866e-07, -3.3060414542221477e-08, -7.078279744199226e-08)
>>>(-8.716465172668348e-08, -7.342018386722812e-09)
>>>([[0.9999999999999942, -7.078279744199198e-08, 8.056217146976134e-08], 
     [7.078279477857338e-08, 0.9999999999999969, 3.3060414542221364e-08], 
     [-8.056217380986972e-08, -3.306040883980552e-08, 0.9999999999999962]], 
    [[0.9999989300532289, -0.0013416472267918239, -0.0005829880927190295], 
     [0.0013416472310697584, 0.9999990999908749, -3.8374444417782927e-07], 
     [0.0005829880828740957, -3.9842032673043093e-07, 0.9999998300623537]], 
    [[0.9999989300052243, -0.0013417179902397033, -0.0005829075749891683], 
     [0.0013417180138317393, 0.9999990998959188, -3.505759733759999e-07], 
     [0.0005829075206857718, -4.3152199547940834e-07, 0.9999998301093036]])
>>>[[0.9999989300052243, -0.0013417179902397033, -0.0005829075749891683], 
    [0.0013417180138317393, 0.9999990998959188, -3.505759733759999e-07], 
    [0.0005829075206857718, -4.3152199547940834e-07, 0.9999998301093036]]
>>>-8.834193235367966e-06
>>>-8.83419245922342e-06
>>>2.0460850048851376e-09
```

pymNut00a(DATE1, DATE2) : Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation with free core nutation omitted).

pymNut00b(DATE1, DATE2) : Nutation, IAU 2000B model.

pymPn00(DATE1, DATE2, DPSI, DEPS) : Precession-nutation, IAU 2000 model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

pymPn00a(DATE1, DATE2) : Precession-nutation, IAU 2000A model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

pymNum00a(DATE1, DATE2) : Form the matrix of nutation for a given date, IAU 2000A model.

pymPn00b(DATE1, DATE2) : Precession-nutation, IAU 2000B model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

pymNum00b(DATE1, DATE2) : Form the matrix of nutation for a given date, IAU 2000B model.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 53736.0
DPSI = -0.9632552291149335877e-5	#nutation
DEPS =  0.4063197106621141414e-4	#nutation

print(sf.pymNut00a(DATE1, DATE2))
print(sf.pymNut00b(DATE1, DATE2))
print(sf.pymPn00(DATE1, DATE2, DPSI, DEPS))
print(sf.pymPn00a(DATE1, DATE2))
print(sf.pymNum00a(DATE1, DATE2))
print(sf.pymPn00b(DATE1, DATE2))
print(sf.pymNum00b(DATE1, DATE2))

>>>(-9.630909107116424e-06, 4.0632391740016646e-05)
>>>(-9.632552291148318e-06, 4.063197106621162e-05)
>>>(0.409079178940423, 
    [[0.9999999999999942, -7.078279744199198e-08, 8.056217146976134e-08], 
     [7.078279477857338e-08, 0.9999999999999969, 3.3060414542221364e-08], 
     [-8.056217380986972e-08, -3.306040883980552e-08, 0.9999999999999962]], 
    [[0.9999989300532289, -0.0013416472267918239, -0.0005829880927190295], 
     [0.0013416472310697584, 0.9999990999908749, -3.8374444417782927e-07], 
     [0.0005829880828740957, -3.9842032673043093e-07, 0.9999998300623537]],
    [[0.9999989300052243, -0.0013417179902397033, -0.0005829075749891683], 
     [0.0013417180138317393, 0.9999990998959188, -3.505759733759999e-07], 
     [0.0005829075206857718, -4.3152199547940834e-07, 0.9999998301093036]], 
    [[0.999999999953607, 8.83774614487214e-06, 3.83148883825259e-06], 
     [-8.837590456633198e-06, 0.9999999991354693, -4.063198798553991e-05],
     [-3.831847930135328e-06, 4.063195412257192e-05, 0.9999999991671806]], 
    [[0.9999989440499982, -0.001332880253640848, -0.0005790760898731086], 
     [0.001332856746979948, 0.9999991109064768, -4.0977405557194263e-05], 
     [0.0005791301929950204, 4.020553681376788e-05, 0.9999998314958529]])
>>>(-9.630909107116424e-06, 4.0632391740016646e-05, 0.409079178940423, 
    [[0.9999999999999942, -7.078279744199198e-08, 8.056217146976134e-08], 
     [7.078279477857338e-08, 0.9999999999999969, 3.3060414542221364e-08], 
     [-8.056217380986972e-08, -3.306040883980552e-08, 0.9999999999999962]], 
    [[0.9999989300532289, -0.0013416472267918239, -0.0005829880927190295], 
     [0.0013416472310697584, 0.9999990999908749, -3.8374444417782927e-07], 
     [0.0005829880828740957, -3.9842032673043093e-07, 0.9999998300623537]], 
    [[0.9999989300052243, -0.0013417179902397033, -0.0005829075749891683], 
     [0.0013417180138317393, 0.9999990998959188, -3.505759733759999e-07],
     [0.0005829075206857718, -4.3152199547940834e-07, 0.9999998301093036]], 
    [[0.9999999999536228, 8.836238544091704e-06, 3.830835237722761e-06], 
     [-8.836082880799401e-06, 0.9999999991354656, -4.063240865359585e-05], 
     [-3.831194272066356e-06, 4.0632374802118676e-05, 0.9999999991671661]],
    [[0.9999989440476101, -0.0013328817612400104, -0.0005790767434730081],
     [0.0013328582543089534, 0.9999991109044503, -4.097782710395611e-05], 
     [0.0005791308472168148, 4.02059566159112e-05, 0.9999998314954572]])
>>>[[0.9999999999536228, 8.836238544091704e-06, 3.830835237722761e-06], 
    [-8.836082880799401e-06, 0.9999999991354656, -4.063240865359585e-05], 
    [-3.831194272066356e-06, 4.0632374802118676e-05, 0.9999999991671661]]
>>>(-9.632552291148318e-06, 4.063197106621162e-05, 0.409079178940423, 
    [[0.9999999999999942, -7.078279744199198e-08, 8.056217146976134e-08], 
     [7.078279477857338e-08, 0.9999999999999969, 3.3060414542221364e-08], 
     [-8.056217380986972e-08, -3.306040883980552e-08, 0.9999999999999962]], 
    [[0.9999989300532289, -0.0013416472267918239, -0.0005829880927190295], 
     [0.0013416472310697584, 0.9999990999908749, -3.8374444417782927e-07], 
     [0.0005829880828740957, -3.9842032673043093e-07, 0.9999998300623537]], 
    [[0.9999989300052243, -0.0013417179902397033, -0.0005829075749891683],
     [0.0013417180138317393, 0.9999990998959188, -3.505759733759999e-07], 
     [0.0005829075206857718, -4.3152199547940834e-07, 0.9999998301093036]], 
    [[0.999999999953607, 8.837746144871207e-06, 3.8314888382521855e-06], 
     [-8.837590456632263e-06, 0.9999999991354693, -4.063198798553991e-05], 
     [-3.831847930134923e-06, 4.063195412257192e-05, 0.9999999991671806]], 
    [[0.9999989440499982, -0.0013328802536408488, -0.000579076089873109], 
     [0.0013328567469799491, 0.9999991109064768, -4.0977405557194263e-05], 
     [0.0005791301929950209, 4.020553681376788e-05, 0.9999998314958529]])
>>>[[0.999999999953607, 8.837746144871207e-06, 3.8314888382521855e-06],
    [-8.837590456632263e-06, 0.9999999991354693, -4.063198798553991e-05], 
    [-3.831847930134923e-06, 4.063195412257192e-05, 0.9999999991671806]]
```

pymEra00(DJ1, DJ2) : Earth rotation angle (IAU 2000 model).

pymSp00(DATE1, DATE2) : The TIO locator s', positioning the Terrestrial Intermediate Origin on the equator of the Celestial Intermediate Pole.

pymPom00(XP, YP, SP) : Form the matrix of polar motion for a given date, IAU 2000.

pymPnm00a(DATE1, DATE2) : Form the matrix of precession-nutation for a given date (including frame bias), equinox based, IAU 2000A model.

pymPnm00b(DATE1, DATE2) : Form the matrix of precession-nutation for a given date (including frame bias), equinox-based, IAU 2000B model.

pymS00(DATE1, DATE2, X, Y) : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, given the CIP's X,Y coordinates.  Compatible with IAU 2000A precession-nutation.

pymS00a(DATE1, DATE2) : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2000A precession-nutation model.

pymS00b(DATE1, DATE2) : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2000B precession-nutation model.

pymXys00a(DATE1, DATE2) : For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole and the CIO locator s, using the IAU 2000A precession-nutation model.

pymXys00b(DATE1, DATE2) : For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole and the CIO locator s, using the IAU 2000B precession-nutation model.

```python
import PyMsOfa as sf

DJ1 = 2400000.5
DJ2 = 54388.0
DATE1 = 2400000.5
DATE2 = 52541.0
XP = 2.55060238e-7			#coordinates of the pole (radians)
YP = 1.860359247e-6			#coordinates of the pole (radians)   
SP = -0.1367174580728891460e-10		#the TIO locator s' (radians)
X = 0.5791308486706011000e-3		#CIP coordinates 
Y = 0.4020579816732961219e-4		#CIP coordinates 

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

>>>0.40228372400281387
>>>-6.216698469981019e-12
>>>[[0.9999999999999675, -1.367174580728847e-11, 2.5506023799999727e-07], 
    [1.4146249479570297e-11, 0.9999999999982695, -1.8603592469988663e-06], 
    [-2.550602379741216e-07, 1.860359247002414e-06, 0.999999999998237]]
>>>[[0.9999998292518457, -0.000535976001895725, -0.00023286477843841014],
    [0.0005359720436122114, 0.9999998562214474, -1.70602781741906e-05],
    [0.00023287388885713628, 1.693546625001527e-05, 0.9999999727414706]]
>>>[[0.9999998292513887, -0.0005359767192496032, -0.0002328650894441306],
    [0.0005359727607393593, 0.9999998562210466, -1.706122925647015e-05],
    [0.00023287420038471497, 1.6936416998393357e-05, 0.999999972741382]]
>>>-2.3077139553642807e-08
>>>-1.3406844489205477e-08
>>>-1.3406957829522032e-08
>>>(0.00023287388885713628, 1.693546625001527e-05, -1.3406844489205477e-08)
>>>(0.00023287420038471497, 1.6936416998393357e-05, -1.3406957829522032e-08)
```

pymBp06(DATE1, DATE2) : Frame bias and precession, IAU 2006.

pymNut06a(DATE1, DATE2) : IAU 2000A nutation with adjustments to match the IAU 2006 precession.

pymObl06(DATE1, DATE2) : Mean obliquity of the ecliptic, IAU 2006 precession model.

pymPfw06(DATE1, DATE2) : Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

pymEe06a(DATE1, DATE2) : Equation of the equinoxes, compatible with IAU 2000 resolutions and IAU 2006/2000A precession-nutation.

pymPn06(DATE1, DATE2, DPSI, DEPS) : Precession-nutation, IAU 2006 model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

pymPn06a(DATE1, DATE2) : Precession-nutation, IAU 2006/2000A models:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 50124.0
DPSI = -0.9632552291149335877e-5	#nutation 
DEPS =  0.4063197106621141414e-4	#nutation 

print(sf.pymBp06(DATE1, DATE2))
print(sf.pymNut06a(DATE1, DATE2))
print(sf.pymObl06(DATE1, DATE2))
print(sf.pymPfw06(DATE1, DATE2))
print(sf.pymEe06a(DATE1, DATE2))
print(sf.pymPn06(DATE1, DATE2, DPSI, DEPS))
print(sf.pymPn06a(DATE1, DATE2))

>>>([[0.9999999999999941, -7.078368960971556e-08, 8.056213977613186e-08], 
     [7.078368694637676e-08, 0.9999999999999969, 3.3059437354321375e-08], 
     [-8.056214211620057e-08, -3.305943169218395e-08, 0.9999999999999962]], 
    [[0.9999995504865591, 0.0008696111966676076, 0.0003778929027311164], 
     [-0.0008696111948330861, 0.9999996218880989, -1.6916459330671076e-07], 
     [-0.00037789290695273693, -1.5945538146359447e-07, 0.9999999285984602]], 
    [[0.999999550517664, 0.0008695404005168889, 0.0003779734935835342],
     [-0.0008695404111592707, 0.9999996219496559, -1.3617522615083644e-07], 
     [-0.0003779734691003856, -1.9248806210070413e-07, 0.9999999285680072]])
>>>(3.5207552633336266e-05, -4.145326347441322e-05)
>>>0.409101431658115
>>>(-2.2433875313873983e-06, 0.4091014602385099, -0.0009501953509247532, 0.409101431658115)
>>>3.229806317150974e-05
>>>(0.409101431658115,
    [[0.9999999999999941, -7.078368960971556e-08, 8.056213977613186e-08], 
     [7.078368694637676e-08, 0.9999999999999969, 3.3059437354321375e-08], 
     [-8.056214211620057e-08, -3.305943169218395e-08, 0.9999999999999962]], 
    [[0.9999995504865591, 0.0008696111966676076, 0.0003778929027311164], 
     [-0.0008696111948330861, 0.9999996218880989, -1.6916459330671076e-07], 
     [-0.00037789290695273693, -1.5945538146359447e-07, 0.9999999285984602]], 
    [[0.999999550517664, 0.0008695404005168889, 0.0003779734935835342],
     [-0.0008695404111592707, 0.9999996219496559, -1.3617522615083644e-07],
     [-0.0003779734691003856, -1.9248806210070413e-07, 0.9999999285680072]], 
    [[0.9999999999536069, 8.83766088164433e-06, 3.831685501173984e-06], 
     [-8.837505185415067e-06, 0.9999999991354698, -4.0631987986276474e-05], 
     [-3.832044589592164e-06, 4.063195412177061e-05, 0.9999999991671797]], 
    [[0.9999995413382922, 0.0008783780572795583, 0.00038180517758999743], 
     [-0.0008783625538071836, 0.9999996134083795, -4.077150065256996e-05], 
     [-0.0003818408427788507, 4.0436118581343994e-05, 0.9999999262812428]])
>>>(3.5207552633336266e-05, -4.145326347441322e-05, 0.409101431658115, 
    [[0.9999999999999941, -7.078368960971556e-08, 8.056213977613186e-08], 
     [7.078368694637676e-08, 0.9999999999999969, 3.3059437354321375e-08], 
     [-8.056214211620057e-08, -3.305943169218395e-08, 0.9999999999999962]], 
    [[0.9999995504865591, 0.0008696111966676076, 0.0003778929027311164], 
     [-0.0008696111948330861, 0.9999996218880989, -1.6916459330671076e-07], 
     [-0.00037789290695273693, -1.5945538146359447e-07, 0.9999999285984602]], 
    [[0.999999550517664, 0.0008695404005168889, 0.0003779734935835342], 
     [-0.0008695404111592707, 0.9999996219496559, -1.3617522615083644e-07], 
     [-0.0003779734691003856, -1.9248806210070413e-07, 0.9999999285680072]], 
    [[0.9999999993802141, -3.230217715739815e-05, -1.4005038836398655e-05], 
     [3.230275768420892e-05, 0.9999999986190887, 4.1453037261891995e-05], 
     [1.4003699793705576e-05, -4.145348963758524e-05, 0.9999999990427523]], 
    [[0.9999995832794599, 0.0008372382377282141, 0.00036396845991203705], 
     [-0.0008372533349421212, 0.9999996486493189, 4.132906866105568e-05], 
     [-0.00036393372975485804, -4.163378524524308e-05, 0.9999999329094319]])
```

pymNum06a(DATE1, DATE2) : Form the matrix of nutation for a given date, IAU 2006/2000A model.

pymPnm06a(DATE1, DATE2) : Form the matrix of precession-nutation for a given date (including frame bias), equinox based, IAU 2006 precession and IAU 2000A nutation models.

pymPb06(DATE1, DATE2) : This function forms three Euler angles which implement general precession from epoch J2000.0, using the IAU 2006 model.  Frame bias (the offset between ICRS and mean J2000.0) is included.

pymPmat06(DATE1, DATE2) : Precession matrix (including frame bias) from GCRS to a specified date, IAU 2006 model.

pymP06e(DATE1, DATE2) : Precession angles, IAU 2006, equinox based.

pymS06(DATE1, DATE2, X, Y) : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, given the CIP's X,Y coordinates.  Compatible with IAU 2006/2000A precession-nutation.

pymS06a(DATE1, DATE2) : The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2006 precession and IAU 2000A nutation models.

pymXy06(DATE1, DATE2) : X,Y coordinates of celestial intermediate pole from series based on IAU 2006 precession and IAU 2000A nutation.

pymXys06a(DATE1, DATE2) : For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole and the CIO locator s, using the IAU 2006 precession and IAU 2000A nutation models.

pymFw2m(GAMB, PHIB, PSI, EPS) : Form rotation matrix given the Fukushima-Williams angles.

pymFw2xy(GAMB, PHIB, PSI, EPS) : CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 50124.0
X = 0.5791308486706011000e-3		#CIP coordinates
Y = 0.4020579816732961219e-4		#CIP coordinates
GAMB = -0.2243387670997992368e-5	#F-W angle gamma_bar (radians)
PHIB =  0.4091014602391312982		#F-W angle phi_bar (radians)
PSI  = -0.9501954178013015092e-3	#F-W angle psi (radians)
EPS  =  0.4091014316587367472		#F-W angle epsilon (radians)

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

>>>[[0.9999999993802141, -3.2302177157398035e-05, -1.400503883639873e-05], 
    [3.2302757684209126e-05, 0.9999999986190887, 4.145303726182936e-05],
    [1.4003699793705669e-05, -4.1453489637577334e-05, 0.9999999990427522]]
>>>[[0.9999995832794599, 0.0008372382377282141, 0.00036396845991203705], 
    [-0.0008372533349421212, 0.9999996486493189, 4.132906866105568e-05], 
    [-0.00036393372975485804, -4.163378524524308e-05, 0.9999999329094319]]
>>>(-0.0005092633773821592, -0.00036027716908913136, -0.00037797352711374776)
>>>[[0.999999550517664, 0.0008695404005168889, 0.0003779734935835342], 
    [-0.0008695404111592707, 0.9999996219496559, -1.3617522615083644e-07], 
    [-0.0003779734691003856, -1.9248806210070413e-07, 0.9999999285680072]]
>>>(0.4090926006005829, -0.0009500121640916875, 0.4090926058345976, 
    -7.903154196964373e-07, 8.826577573750554e-06, -8.861888526875804e-06, 
    3.0522926386722435, 0.409101431658115, -2.007869257151149e-06, 
    -0.00044765239229031437, -0.00042195894411613553, -0.00037789294961742447, 
    -0.0009481699832898824, -1.9867900599926174e-06, 0.40910142717887715, -0.0009499928243453799)
>>>-7.773082059855285e-09
>>>(-0.00036393372958495174, -4.163378483283794e-05)
>>>(-0.00036393372975485804, -4.163378524524308e-05, -3.706842426850084e-09)
>>>[[0.9999995505176007, 0.0008695404617348192, 0.00037797352018655825], 
    [-0.0008695404723772016, 0.9999996219496027, -1.3617524957654226e-07], 
    [-0.00037797349570340826, -1.9248808486027613e-07, 0.9999999285679972]]
>>>(-0.00037797349570340826, -1.9248808486027613e-07)
```

### Section 4 Transformation of the coordinate system

pymEform(N) : Earth reference ellipsoids.(1 : WGS84 ; 2 : GRS80 ; 3 : WGS72)

pymGd2gce(A, F, ELONG, PHI, HEIGHT) : Transform geodetic coordinates to geocentric for a reference ellipsoid of specified form.

pymGd2gc(N, ELONG, PHI, HEIGHT) : Transform geodetic coordinates to geocentric using the specified reference ellipsoid.

pymGc2gde(A, F, XYZ) : Transform geocentric coordinates to geodetic for a reference ellipsoid of specified form.

pymGc2gd(N, XYZ) : Transform geocentric coordinates to geodetic using the specified reference ellipsoid.

pymAe2hd(AZ, EL, PHI) : Horizon to equatorial coordinates:  transform azimuth and altitude to hour angle and declination.

pymHd2ae(HA, DEC, PHI) : Equatorial to horizon coordinates:  transform hour angle and declination to azimuth and altitude.

pymHd2pa(HA, DEC, PHI) : Parallactic angle for a given hour angle and declination.

pymLtecm(EPJ) : ICRS equatorial to ecliptic rotation matrix, long-term.

pymLteqec(EPJ, DR, DD) : Transformation from ICRS equatorial coordinates to ecliptic coordinates (mean equinox and ecliptic of date) using a long-term precession model.

pymLteceq(EPJ, DL, DB) : Transformation from ecliptic coordinates (mean equinox and ecliptic of date) to ICRS RA,Dec, using a long-term precession model.

```python
import PyMsOfa as sf
import math as ma

N = 2 				# ellipsoid identifier
A = 6378136.0			#equatorial radius
F = 0.0033528			#flattening
ELONG = 3.1			#longitude (radians, east +ve)
PHI = -0.5			#latitude (geodetic, radians)
HEIGHT = 2500.0			#height above ellipsoid (geodetic)
XYZ = [2e6, 3e6, 5.244e6]	#geocentric vector
AZ = 5.5			#azimuth
EL = 1.1			#altitude (informally, elevation)
PHI = 0.7			#site latitude
HA = 1.1			#hour angle (local)
DEC = 1.2			#declination
PHI = 0.3			#site latitude
EPJ = 2024.154688569473		#Julian epoch (TT)
DR = ma.pi/3			#ICRS right ascension (radians) 
DD = ma.pi/4			#ICRS declination (radians)
DL = 1.5			#ecliptic longitude (radians)
DB = 0.6			#ecliptic latitude (radians)

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

>>>(6378137.0, 0.003352810681182319, 0)					#The last value is the check value, where 0 is OK.
>>>([-6092162.972137791, 253535.44209085943, 1873536.416041204], 0)	#The last value is the check value, where 0 is OK.
>>>([-6092163.932592076, 253535.4820617537, 1873536.6712690855], 0)	#The last value is the check value, where 0 is OK.
>>>(0.982793723247329, 0.9716018377570411, 332.36862495798647, 0)	#The last value is the check value, where 0 is OK.
>>>(0.982793723247329, 0.9716018482060784, 331.41731754978827, 0)	#The last value is the check value, where 0 is OK.
>>>(0.40025365757863, 0.6070689395289157)
>>>(5.916889243730066, 0.4472186304990486)
>>>1.9062274280019955
>>>[[0.9999826579404173, -0.005401514272291263, -0.0023467982435674247], 
    [0.0058892949248147515, 0.917488028189624, 0.39771979374178257], 
    [4.870150781661126e-06, -0.39772671744841026, 0.9175039281681606]]
>>>(1.1796084210776348, 0.4172173738451268)
>>>(1.4521514559004196, 1.0072288157204605)
```

pymIcrs2g(DR, DD) : Transformation from ICRS to Galactic Coordinates.

pymG2icrs(DL, DB) : Transformation from Galactic Coordinates to ICRS.

pymFk425(R1950, D1950, DR1950, DD1950, P1950, V1950) : Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

pymFk45z(R1950, D1950, BEPOCH) : Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero proper motion in the FK5 system.

pymFk524(R2000, D2000, DR2000, DD2000, P2000, V2000) : Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

pymFk54z(R2000, D2000, BEPOCH) : Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero proper motion in FK5 and parallax.

pymFk5hip() : FK5 to Hipparcos rotation and spin.

pymFk52h(R5, D5, DR5, DD5, PX5, RV5) : Transform FK5 (J2000.0) star data into the Hipparcos system.

```python
import PyMsOfa as sf
import math as ma

DR = ma.pi/3				#ICRS right ascension (radians)
DD = ma.pi/4				#ICRS declination (radians)
DL = ma.pi/6				#galactic longitude (radians)
DB = ma.pi/8				#galactic longitude (radians)
R1950 = 0.07626899753879587532		#B1950.0 RA (rad)
D1950 = -1.137405378399605780		#B1950.0 Dec (rad)
DR1950 = 0.1973749217849087460e-4	#B1950.0 proper motions (rad/trop.yr)
DD1950 = 0.5659714913272723189e-5	#B1950.0 proper motions (rad/trop.yr)
P1950 = 0.134				#parallax (arcsec)
V1950 = 8.7				#radial velocity (km/s, +ve = moving away)
R2000 = 0.01602284975382960982		#J2000.0 RA (rad)
D2000 = -0.1164347929099906024		#J2000.0 Dec (rad)
BEPOCH = 1954.677617625256806		#Besselian epoch (e.g. 1979.3)
DR2000 = 0.2019447755430472323e-4	#J2000.0 proper motions (rad/Jul.yr)
DD2000 = 0.3541563940505160433e-5	#J2000.0 proper motions (rad/Jul.yr)
P2000 = 0.1559				#parallax (arcsec)
V2000 = 86.87				#radial velocity (km/s, +ve = moving away)
R5  =  1.76779433			#RA (radians) 	FK5, equinox J2000.0, epoch J2000.0
D5  = -0.2917517103			#Dec (radians)  FK5, equinox J2000.0, epoch J2000.0
DR5 = -1.91851572e-7			#proper motion in RA (dRA/dt, rad/Jyear) 	FK5, equinox J2000.0, epoch J2000.0
DD5 = -5.8468475e-6			#proper motion in Dec (dDec/dt, rad/Jyear)	FK5, equinox J2000.0, epoch J2000.0
PX5 =  0.379210				#parallax (arcsec)	FK5, equinox J2000.0, epoch J2000.0
RV5 = -7.6				#radial velocity (km/s, positive = receding)	FK5, equinox J2000.0, epoch J2000.0

print(sf.pymIcrs2g(DR,DD))
print(sf.pymG2icrs(DL, DB))
print(sf.pymFk425(R1950, D1950, DR1950, DD1950, P1950, V1950))
print(sf.pymFk45z(R1950, D1950, BEPOCH))
print(sf.pymFk524(R2000, D2000, DR2000, DD2000, P2000, V2000))
print(sf.pymFk54z(R2000, D2000, BEPOCH))
print(sf.pymFk5hip())
print(sf.pymFk52h(R5, D5, DR5, DD5, PX5, RV5))

>>>(2.6938731316889215, -0.10583104427828448)
>>>(4.562820525414904, 0.13283058755100796)
>>>(0.08757989933556447, -1.132279113042092, 1.953670614474396e-05, 
    5.63768667865964e-06, 0.1339919950582768, 8.73699966918353)
>>>(0.08660249025427784, -1.1325609302769601)
>>>(0.0038357827667403507, -0.12147115357961685, 2.0222718742897512e-05, 
    3.566078662468435e-06, 0.15600796032308914, 86.73793467537503)
>>>(0.004846756109282368, -0.12129386288333069, -1.1732616691106914e-08, 2.109122543091408e-08)
>>>([[0.9999999999999929, 1.1102233509835495e-07, 4.411803963527301e-08], 
     [-1.1102233084981164e-07, 0.9999999999999892, -9.647792498531089e-08], 
     [-4.4118050326662195e-08, 9.647792009628368e-08, 0.9999999999999943]], 
    [-1.4544410433286077e-09, 2.9088820866572155e-09, 3.3936957677667517e-09])
>>>(1.7677942262999478, -0.2917516070530392, -1.9618741256057225e-07,
   -5.845990517669392e-06, 0.37921000000000005, -7.600000094000024)
```

pymFk5hz(R5, D5, DATE1, DATE2) : Transform an FK5 (J2000.0) star position into the system of the Hipparcos catalog, assuming zero Hipparcos proper motion.

pymH2fk5(RH, DH, DRH, DDH, PXH, RVH) : Transform Hipparcos star data into the FK5 (J2000.0) system.

pymHfk5z(RH, DH, DATE1, DATE2) : Transform a Hipparcos star position into FK5 J2000.0, assuming = zero Hipparcos proper motion.

pymEcm06(DATE1, DATE2) : ICRS equatorial to ecliptic rotation matrix, IAU 2006.

pymEceq06(DATE1, DATE2, DL, DB) : Transformation from ecliptic coordinates (mean equinox and ecliptic of date) to ICRS RA,Dec, using the IAU 2006 precession model.

pymEqec06(DATE1, DATE2, DR, DD) : Transformation from ICRS equatorial coordinates to ecliptic coordinates (mean equinox and ecliptic of date) using IAU 2006 precession model.

pymStarpm(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B) : Star proper motion:  update star catalog data for space motion.

pymPmsafe(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B) : Star proper motion:  update star catalog data for space motion, with special handling to handle the zero parallax case.

```python
import PyMsOfa as sf

R5 =  1.76779433		#FK5 RA (radians), equinox J2000.0, at date
D5 = -0.2917517103		#FK5 Dec (radians), equinox J2000.0, at date
DATE1 = 2400000.5
DATE2 = 54479.0
RH  =  1.767794352		#RA (radians)		Hipparcos, epoch J2000.0
DH  = -0.2917512594		#Dec (radians)		Hipparcos, epoch J2000.0
DRH = -2.76413026e-6		#proper motion in RA (dRA/dt, rad/Jyear)	Hipparcos, epoch J2000.0
DDH = -5.92994449e-6		#proper motion in Dec (dDec/dt, rad/Jyear)	Hipparcos, epoch J2000.0
PXH =  0.379210			#parallax (arcsec)  	Hipparcos, epoch J2000.0
RVH = -7.6			#radial velocity (km/s, positive = receding)	Hipparcos, epoch J2000.0
DL = 5.1			#ecliptic longitude (radians)
DB = -0.9			#ecliptic latitude (radians)
DR = 1.234			#ICRS right ascension (radians)
DD = 0.987			#ICRS declination (radians)
RA1 =   0.01686756		#right ascension (radians), before
DEC1 = -1.093989828		#declination (radians), before
PMR1 = -1.78323516e-5		#RA proper motion (radians/year), before
PMD1 =  2.336024047e-6		#Dec proper motion (radians/year), before
PX1 =   0.74723			#parallax (arcseconds), before
RV1 = -21.6			#radial velocity (km/s, +ve = receding), before
EP1A = 2400000.5		#"before" epoch, part A
EP1B = 50083.0			#"before" epoch, part B
EP2A = 2400000.5		#"after" epoch, part A
EP2B = 53736.0			#"after" epoch, part B

print(sf.pymFk5hz(R5, D5, DATE1, DATE2))
print(sf.pymH2fk5(RH, DH, DRH, DDH, PXH, RVH))
print(sf.pymHfk5z(RH, DH, DATE1, DATE2))
print(sf.pymEcm06(DATE1, DATE2))
print(sf.pymEceq06(DATE1, DATE2, DL, DB))
print(sf.pymEqec06(DATE1, DATE2, DR, DD))
print(sf.pymStarpm(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))
print(sf.pymPmsafe(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))

>>>(1.7677942259515924, -0.2917516069841887)
>>>(1.7677944557000653, -0.2917513626469639, -2.759794502451122e-06,
   -5.9308014093262826e-06, 0.37920999999999994, -7.600000130907111)
>>>(1.767794490535581, -0.2917513695320114, 4.335890983539242e-09, -8.56964884123775e-10)
>>>[[0.9999980814830622, -0.0017965963501775856, -0.0007805586135427632], 
    [0.0019588332813974653, 0.9174876229215057, 0.3977595061818394],
    [1.539589835524444e-06, -0.39776027205768, 0.9174893819386849]]
>>>(5.535423121517396, -1.2462214232485105)
>>>(1.3513428093792499, 0.5927034347268125)
>>>(0.01668919069414256, -1.093966454217128, -1.7836626821531766e-05, 
    2.3380929159839913e-06, 0.7473533835317719, -21.599051704764175, 0)		#The last value is the check value, where 0 is OK.
>>>(0.01668919069414256, -1.093966454217128, -1.7836626821531766e-05, 
    2.3380929159839913e-06, 0.7473533835317719, -21.599051704764175, 0)		#The last value is the check value, where 0 is OK.
```

### Section 5 Astrometric parameter

pymEpv00(DATE1, DATE2) : Earth position and velocity, heliocentric and barycentric, with respect to the Barycentric Celestial Reference System.

pymAb(PNAT, V, S, BM1) : Apply aberration to transform natural direction into proper direction.

pymPmpx(RC, DC, PR, PD, PX, RV, PMT, POB):The position of the star after a certain time interval (PMT) is given based on the star's position, proper motion, parallax, and radial velocity and the observer's velocity relative to the Solar SystemBarycenter.

pymRefco(PHPA, TC, RH, WL) : Determine the constants A and B in the atmospheric refraction model dZ = A tan Z + B tan^3 Z.

pymPvtob(ELONG, PHI, HM, XP, YP, SP, THETA) : Position and velocity of a terrestrial observing station.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 60367.0
PNAT = [-0.76321968546737951, -0.60869453983060384, -0.21676408580639883]	#natural direction to the source (unit vector)
V = [2.1044018893653786E-5, -8.9108923304429319E-5, -3.8633714797716569E-5]	#observer barycentric velocity in units of c
S = 0.99980921395708788			#distance between the Sun and the observer (au)
BM1 = 0.99999999506209258		#sqrt(1-|v|^2): reciprocal of Lorenz factor
RC = 1.234				#ICRS RA at catalog epoch (radians)
DC = 0.789				#ICRS Dec at catalog epoch (radians)
PR = 1e-5				#RA proper motion (radians/year)
PD = -2e-5				#Dec proper motion (radians/year)
PX = 1e-2				#parallax (arcsec)
RV = 10.0				#radial velocity (km/s, +ve if receding)
PMT = 8.75				#proper motion time interval (SSB, Julian years)
POB = [0.9, 0.4, 0.1]			#SSB to observer vector (au)
PHPA = 800.0				#pressure at the observer (hPa = millibar)
TC = 10.0				#ambient temperature at the observer (deg C)
RH = 0.9				#relative humidity at the observer (range 0-1)
WL = 0.4				#wavelength (micrometers)
ELONG = 2.0				#longitude (radians, east +ve)
PHI = 0.5				#latitude (geodetic, radians)
HM = 3000.0				#height above ref. ellipsoid (geodetic, m)
XP = 1e-6				#coordinates of the pole (radians)
YP = -0.5e-6				#coordinates of the pole (radians)
SP = 1e-8				#the TIO locator s' (radians)
THETA = 5.0				#Earth rotation angle (radians)

print(sf.pymEpv00(DATE1, DATE2))
print(sf.pymAb(PNAT, V, S, BM1))
print(sf.pymPmpx(RC, DC, PR, PD, PX, RV, PMT, POB))
print(sf.pymRefco(PHPA, TC, RH, WL))
print(sf.pymPvtob(ELONG, PHI, HM, XP, YP, SP, THETA))

>>>([[-0.9150367263184106, 0.34697147073735446, 0.15041286078980728], 
     [-0.006851234370024847, -0.014641858672797554, -0.006346351691059997]], 
    [[-0.9227148846354889, 0.3438656140831509, 0.1492922882346728], 
     [-0.006845957992232367, -0.014647909564465206, -0.006349037736469475]], 0)	#The last value is the check value, where 0 is OK.
>>>[-0.7631631094219556, -0.6087553082505591, -0.21679262693684712]
>>>[0.23281376239603083, 0.6651097085397856, 0.7095257765896361]
>>>(0.0002264949956241415, -2.598658261729344e-07)
>>>[[4225081.36707116, 3681943.2158561987, 3041149.3992412607], 
    [-268.49153893659997, 308.0977983288904, 0.0]]
```

pymC2ixy(DATE1, DATE2, X, Y) : Form the celestial to intermediate-frame-of-date matrix for a given date when the CIP X,Y coordinates are known.  IAU 2000.

pymC2ixys(X, Y, S) : Form the celestial to intermediate-frame-of-date matrix given the CIP X,Y and the CIO locator s.

pymC2ibpn(DATE1, DATE2, RBPN) : Form the celestial-to-intermediate matrix for a given date given the bias-precession-nutation matrix.  IAU 2000.

pymC2i00a(DATE1, DATE2) : Form the celestial-to-intermediate matrix for a given date using the IAU 2000A precession-nutation model.

pymC2i00b(DATE1, DATE2) : Form the celestial-to-intermediate matrix for a given date using the IAU 2000B precession-nutation model.

pymC2i06a(DATE1, DATE2) : Form the celestial-to-intermediate matrix for a given date using the IAU 2006 precession and IAU 2000A nutation models.

```python
import PyMsOfa as sf

DATE1 = 2400000.5
DATE2 = 60367.0
X = 0.5791308486706011000e-3		#Celestial Intermediate Pole
Y = 0.4020579816732961219e-4		#Celestial Intermediate Pole
S = -0.1220040848472271978e-7		#the CIO locator s
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

>>>[[0.9999998323037157, 4.055624169829741e-09, -0.0005791308493017451],
    [-2.734004153626124e-08, 0.9999999991917469, -4.020578907622958e-05], 
    [0.0005791308486706012, 4.020579816732961e-05, 0.9999998314954628]]
>>>[[0.9999998323037157, 5.581984874325485e-10, -0.0005791308491611283], 
    [-2.3842616436731134e-08, 0.9999999991917469, -4.020579110169669e-05], 
    [0.0005791308486706012, 4.020579816732961e-05, 0.9999998314954628]]
>>>[[0.9999994021664094, 4.055632212007776e-09, -0.0010934655110439967], 
    [4.275931623359819e-08, 0.9999999990835077, 4.2813351130053356e-05],
    [0.001093465510215479, -4.28133722906315e-05, 0.9999994012499173]]
>>>[[0.9999972645587599, 4.055563843086141e-09, -0.0023389901661903423], 
    [-9.366129322466166e-08, 0.9999999992661879, -3.8309469819015054e-05], 
    [0.0023389901643185965, 3.830958409855523e-05, 0.9999972638249501]]
>>>[[0.9999972645587417, 4.055563843086141e-09, -0.002338990173985668], 
    [-9.366128352061853e-08, 0.9999999992661881, -3.8309465543798595e-05],
    [0.0023389901721139225, 3.830957982332781e-05, 0.9999972638249321]]
>>>[[0.9999972645600166, 4.055975468680861e-09, -0.0023389896289374846],
    [-9.366020251033724e-08, 0.9999999992662121, -3.830883631870127e-05], 
    [0.0023389896270657826, 3.8308950597421e-05, 0.9999972638262311]]
```

pymEors(RNPB, S) : Equation of the origins, given the classical NPB matrix and the quantity s.

pymEo06a(DATE1, DATE2) : Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

pymBpn2xy(RBPN) : Extract from the bias-precession-nutation matrix the X,Y coordinates of the Celestial Intermediate Pole.

```python
import PyMsOfa as sf

RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
	#classical nutation x precession x bias matrix
S = -0.1220040848472271978e-7	#the quantity s (the CIO locator) in radians
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

pymC2tcio(RC2I, ERA, RPOM) : Assemble the celestial to terrestrial matrix from CIO-based components (the celestial-to-intermediate matrix, the Earth Rotation Angle and the polar motion matrix).

pymC2t00a(TTA, TTB, UTA, UTB, XP, YP) : Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2000A precession-nutation model.

pymC2t00b(TTA, TTB, UTA, UTB, XP, YP) : Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2000B precession-nutation model.

pymC2t06a(TTA, TTB, UTA, UTB, XP, YP) : Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2006/2000A precession-nutation model.

pymC2teqx(RBPN, GST, RPOM) : Assemble the celestial to terrestrial matrix from equinox-based components (the celestial-to-true matrix, the Greenwich Apparent Sidereal Time and the polar motion matrix).

pymC2tpe(TTA, TTB, UTA, UTB, DPSI, DEPS, XP, YP) : Form the celestial to terrestrial matrix given the date, the UT1, the nutation and the polar motion.  IAU 2000.

pymC2txy(TTA, TTB, UTA, UTB, X, Y, XP, YP) : Form the celestial to terrestrial matrix given the date, the UT1, the CIP coordinates and the polar motion.  IAU 2000.

```python
import PyMsOfa as sf

RC2I = [[0.9999998323037164738, 0.5581526271714303683e-9, -0.5791308477073443903e-3],
        [-0.2384266227524722273e-7, 0.9999999991917404296, -0.4020594955030704125e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
	#celestial-to-intermediate matrix
ERA = 1.75283325530307		#Earth rotation angle (radians)
RPOM = [[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
        [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
        [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]
	#polar-motion matrix
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7	#coordinates of the pole
YP = 1.860359247e-6	#coordinates of the pole
RBPN = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
	#celestial-to-true matrix
GST = 1.754166138040730516		#Greenwich (apparent) Sidereal Time (radians)
DEPS = 0.4090789763356509900		#nutation
DPSI = -0.9630909107115582393e-5	#nutation
X = 0.5791308486706011000e-3		#coordinates of the pole (radians)
Y = 0.4020579816732961219e-4		#coordinates of the pole (radians)

print(sf.pymC2tcio(RC2I, ERA, RPOM))
print(sf.pymC2t00a(TTA, TTB, UTA, UTB, XP, YP))
print(sf.pymC2t00b(TTA, TTB, UTA, UTB, XP, YP))
print(sf.pymC2t06a(TTA, TTB, UTA, UTB, XP, YP))
print(sf.pymC2teqx(RBPN, GST, RPOM))
print(sf.pymC2tpe(TTA, TTB, UTA, UTB, DPSI, DEPS, XP, YP))
print(sf.pymC2txy(TTA, TTB, UTA, UTB, X, Y, XP, YP))

>>>[[-0.18103321283071114, 0.9834769806938469, 6.555535638685471e-05], 
    [-0.9834768134135996, -0.18103322036494493, 0.0005749801116141106], 
    [0.0005773474014081407, 3.9618323917726584e-05, 0.9999998325501692]]
>>>[[-0.18103321283071463, 0.9834769806938464, 6.555535638688495e-05], 
    [-0.983476813413599, -0.18103322036494846, 0.0005749801116141047], 
    [0.0005773474014081402, 3.961832391769783e-05, 0.9999998325501692]]
>>>[[-0.1810332128439643, 0.9834769806913878, 6.555565082455188e-05],
    [-0.9834768134115444, -0.1810332203783966, 0.0005749793922030075], 
    [0.0005773467471863534, 3.961790411553012e-05, 0.9999998325505636]]
>>>[[-0.1810332128305861, 0.9834769806938598, 6.555550962990436e-05], 
    [-0.983476813413622, -0.1810332203649095, 0.0005749800844905738], 
    [0.0005773474024748544, 3.961816829640605e-05, 0.9999998325501748]]
>>>[[-0.18103321285286855, 0.9834769806897684, 6.555535639982635e-05],
    [-0.9834768134095211, -0.1810332203871024, 0.0005749801116126439], 
    [0.000577347401408154, 3.961832391768641e-05, 0.9999998325501692]]
>>>[[-0.1813677995762994, 0.9023482206891685, -0.3909902938641089], 
    [-0.9834147641476809, -0.16598836354349625, 0.07309763898042682],
    [0.0010596854306732149, 0.39776318556050766, 0.9174875068792734]]
>>>[[-0.1810332128306243, 0.9834769806938525, 6.555551248057454e-05],
    [-0.9834768134136147, -0.1810332203649493, 0.0005749800843594143], 
    [0.0005773474028619265, 3.9618165469116246e-05, 0.9999998325501747]]
```

#### 5.1 Deflection of light by celestial bodies in the solar system

pymLd(BM, P, Q, E, EM, DLIM) : Apply light deflection by a solar-system body, as part of transforming coordinate direction into natural direction.

pymLdn(N, B, OB, SC) : For a star, apply light deflection by multiple solar-system bodies, as part of transforming coordinate direction into natural direction.

pymLdsun(P, E, EM) : Deflection of starlight by the Sun.

```python
import PyMsOfa as sf

BM = 0.00028574					#mass of the gravitating body (solar masses)
P = [-0.763276255, -0.608633767, -0.216735543]	#direction from observer to source (unit vector)
Q = [-0.763276255, -0.608633767, -0.216735543]	#direction from body to source (unit vector)
E = [0.76700421, 0.605629598, 0.211937094]	#direction from body to observer (unit vector)
EM = 8.91276983					#distance from body to observer (au)
DLIM = 3e-10					#deflection limiter
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]	#data for each of the n bodies
OB = [-0.974170437, -0.2115201, -0.0917583114]	#barycentric position of the observer (au)
SC = [-0.763276255, -0.608633767, -0.216735543]	#observer to star coord direction (unit vector)

print(sf.pymLd(BM, P, Q, E, EM, DLIM))
print(sf.pymLdn(N, B, OB, SC))
print(sf.pymLdsun(P, E, EM))

>>>[-0.7632762548968159, -0.6086337670823763, -0.2167355431320547]
>>>[-0.7632762579693334, -0.6086337636093003, -0.21673554206463283]
>>>[-0.7632758938883696, -0.6086340552909992, -0.21673600514983796]
```

#### 5.2 Coordinate transformation based on astrometric parameter ASTROM

```
astrom : list(30)
        star-independent astrometry parameters   
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)
```

##### 5.2.1 Generate ASTROM parameters

The parameters produced by these functions are required in the parallax, light deflection and aberration parts of the astrometric transformation chain.

pymApcg(DATE1, DATE2, EBPV, EHP, ASTROM) : For a geocentric observer, prepare star-independent astrometry parameters for transformations between ICRS and GCRS coordinates. The Earth ephemeris is supplied by the caller.

pymApcg13(DATE1, DATE2, ASTROM) : For a geocentric observer, prepare star-independent astrometry parameters for transformations between ICRS and GCRS coordinates. The caller supplies the date, and SOFA models are used to predict the Earth ephemeris.

pymApci(DATE1, DATE2, EBPV, EHP, X, Y, S, ASTROM) : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and geocentric CIRS coordinates.  The Earth ephemeris and CIP/CIO are supplied by the caller.

pymApci13(DATE1, DATE2, ASTROM) : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and geocentric CIRS coordinates.  The caller supplies the date, and SOFA models are used to predict the Earth ephemeris and CIP/CIO.

pymApco(DATE1, DATE2, EBPV, EHP, X, Y, S, THETA, ELONG, PHI, HM, XP, YP, SP, REFA, REFB) : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and observed coordinates.  The caller supplies the Earth ephemeris, the Earth rotation information and the refraction constants as well as the site coordinates.

pymApco13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL) : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between ICRS and observed coordinates.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength, and SOFA models are used to obtain the Earth ephemeris, CIP/CIO and refraction constants.

pymApcs(DATE1, DATE2, PV, EBPV, EHP, ASTROM) : For an observer whose geocentric position and velocity are known, prepare star-independent astrometry parameters for transformations between ICRS and GCRS.  The Earth ephemeris is supplied by the caller.

pymApcs13(DATE1, DATE2, PV, ASTROM) : For an observer whose geocentric position and velocity are known, prepare star-independent astrometry parameters for transformations between ICRS and GCRS.  The Earth ephemeris is from SOFA models.

pymAper(THETA, ASTROM) : In the star-independent astrometry parameters, update only the Earth rotation angle, supplied by the caller explicitly.

pymAper13(UT11, UT12, ASTROM) : In the star-independent astrometry parameters, update only the Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

pymApio(SP, THETA, ELONG, PHI, HM, XP, YP, REFA, REFB, ASTROM) : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between CIRS and observed coordinates.  The caller supplies the Earth orientation information and the refraction constants as well as the site coordinates.

pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM) : For a terrestrial observer, prepare star-independent astrometry parameters for transformations between CIRS and observed coordinates.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

```python
import PyMsOfa as sf

DATE1 = 2456165.5
DATE2 = 0.401182685
EBPV = [[0.901310875, -0.417402664, -0.180982288],
        [0.00742727954, 0.0140507459, 0.00609045792]]	#Earth barycentric pos/vel (au, au/day)
EHP = [0.903358544, -0.415395237, -0.180084014]		#Earth heliocentric position (au)
X =  0.0013122272			#CIP X,Y (components of unit vector)
Y = -2.92808623e-5			#CIP X,Y (components of unit vector)
S =  3.05749468e-8			#the CIO locator s (radians)
THETA = 3.14540971			#Earth rotation angle (radians)
ELONG = -0.527800806			#longitude (radians, east +ve)
PHI = -1.2345856			#latitude (geodetic, radians)
HM = 2738.0				#height above ellipsoid (m, geodetic)
XP = 2.47230737e-7			#polar motion coordinates (radians)
YP = 1.82640464e-6			#polar motion coordinates (radians)
SP = -3.01974337e-11			#the TIO locator s' (radians)
REFA = 0.000201418779			#refraction constant A (radians)
REFB = -2.36140831e-7			#refraction constant B (radians)
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675			#UT1-UTC (seconds)
PHPA = 731.0				#pressure at the observer (hPa = mB)
TC = 12.8				#ambient temperature at the observer (deg C)
RH = 0.59				#relative humidity at the observer (range 0-1)
WL = 0.55				#wavelength (micrometers)
PV = [[-1836024.09, 1056607.72, -5998795.26],
      [-77.0361767, -133.310856, 0.0971855934]]
UT11 = 2456165.5
UT12 = 0.401182685

ASTROM=[0 for i in range(30)]
print(sf.pymApcg(DATE1, DATE2, EBPV, EHP, ASTROM))
ASTROM=[0 for i in range(30)]
print(sf.pymApcg13(DATE1, DATE2, ASTROM))
ASTROM=[0 for i in range(30)]
print(sf.pymApci(DATE1, DATE2, EBPV, EHP, X, Y, S, ASTROM))
ASTROM=[0 for i in range(30)]
print(sf.pymApci13(DATE1, DATE2, ASTROM))
print(sf.pymApco(DATE1, DATE2, EBPV, EHP, X, Y, S, THETA, ELONG, PHI, HM, XP, YP, SP, REFA, REFB))
print(sf.pymApco13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
ASTROM=[0 for i in range(30)]
print(sf.pymApcs(DATE1, DATE2, PV, EBPV, EHP, ASTROM))
ASTROM=[0 for i in range(30)]
print(sf.pymApcs13(DATE1, DATE2, PV, ASTROM))
ASTROM=[0 for i in range(30)]
ASTROM[21] = 1.234
print(sf.pymAper(THETA, ASTROM))
ASTROM=[0 for i in range(30)]
ASTROM[21] = 1.234
print(sf.pymAper13(UT11, UT12, ASTROM))
ASTROM=[0 for i in range(30)]
print(sf.pymApio(SP, THETA, ELONG, PHI, HM, XP, YP, REFA, REFB, ASTROM))
ASTROM=[0 for i in range(30)]
print(sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM))

>>>[12.651337940273786, 0.901310875, -0.417402664, -0.180982288, 0.8940025429324143, 
    -0.41109302686798177, -0.178218900487287, 1.0104652958110132, 4.289638913597693e-05, 8.11503405158132e-05, 
    3.517555136380563e-05, 0.9999999951686013, 1, 0, 0,
    0, 1, 0, 0, 0, 
    1, 0, 0, 0, 0, 
    0, 0, 0, 0, 0]
>>>[12.651337940273786, 0.9013108747340665, -0.41740266404059817, -0.18098228778677578, 0.8940025429255563, 
    -0.4110930268331781, -0.17821890060197, 1.0104652959646592, 4.2896389129412893e-05, 8.115034032405105e-05,
    3.517555135536498e-05, 0.9999999951686013, 1, 0, 0,
    0, 1, 0, 0, 0, 
    1, 0, 0, 0, 0, 
    0, 0, 0, 0, 0]
>>>[12.651337940273786, 0.901310875, -0.417402664, -0.180982288, 0.8940025429324143, 
    -0.41109302686798177, -0.178218900487287, 1.0104652958110132, 4.289638913597693e-05, 8.11503405158132e-05,
    3.517555136380563e-05, 0.9999999951686013, 0.9999991390295159, -1.1363366539396402e-08, -0.0013122272008952603, 
    4.978650072762214e-08, 0.9999999995713156, 2.9280822178723155e-05, 0.0013122272000000003, -2.9280862300000004e-05, 
    0.9999991386008323, 0, 0, 0, 0,
    0, 0, 0, 0, 0]
>>>([12.651337940273786, 0.9013108747340665, -0.41740266404059817, -0.18098228778677578, 0.8940025429255563, 
     -0.4110930268331781, -0.17821890060197, 1.0104652959646592, 4.2896389129412893e-05, 8.115034032405105e-05, 
     3.517555135536498e-05, 0.9999999951686013, 0.999999206037676, -1.2822919871757765e-08, -0.0012601285716613744, 
     4.124244859862736e-08, 0.9999999997456834, 2.255285422955367e-05, 0.0012601285710517097, -2.2552888294224964e-05, 
     0.9999992057833604, 0, 0, 0, 0,
     0, 0, 0, 0, 0],
     -0.002900618712657376)
>>>[12.651337940273786, 0.9012986019369676, -0.41739560101353446, -0.1810223874694953, 0.8939963497922233,
    -0.41108877432989205, -0.1782597716203172, 1.0104585674728062, 4.263942410653276e-05, 8.070566336357042e-05, 
    3.51758755400555e-05, 0.9999999952155666, 0.9999991390295159, -1.1363366539396402e-08, -0.0013122272008952603, 
    4.978650072762214e-08, 0.9999999995713156, 2.9280822178723155e-05, 0.0013122272000000003, -2.9280862300000004e-05, 
    0.9999991386008323, 0.5201666932102418, 2.542004084904101e-07, 1.8254476433615963e-06, -0.9440115679003211, 
    0.32991235149714754, 0.0, -2.6176089039693444, 0.000201418779, -2.36140831e-07]
>>>([13.252484686224758, -0.9741827107320855, -0.21151301904898623, -0.09179840189497385, -0.9736425572586899, 
     -0.2092452121603479, -0.09075578153886282, 0.9998233240913919, 2.0787049945205506e-05, -8.955360133238837e-05, 
     -3.863338993055873e-05, 0.999999995027756, 0.9999991390295148, -1.1363366528988061e-08, -0.0013122272017455534, 
     4.9786500755377716e-08, 0.9999999995713155, 2.9280822188467148e-05, 0.0013122272008502932, -2.928086230974403e-05, 
     0.9999991386008312, 0.5201666827717384, 2.5420041801788083e-07, 1.8254476420348602e-06, -0.9440115679003211, 
     0.32991235149714754, 0.0, -2.617608909188596, 0.0002014187785940397, -2.3614083149436963e-07], 
     -0.0030205483548024123, 0)						#The last value is the check value, where 0 is OK.
>>>[12.651337940273786, 0.9012986019370117, -0.4173956010136416, -0.1810223874695441, 0.8939963497921855, 
    -0.41108877432996066, -0.17825977162034926, 1.0104585674728983, 4.263942411019233e-05, 8.07056633651959e-05, 
    3.517587554005092e-05, 0.9999999952155666, 1, 0, 0,
    0, 1, 0, 0, 0, 
    1, 0, 0, 0, 0, 
    0, 0, 0, 0, 0]
>>>[12.651337940273786, 0.9012986016710782, -0.4173956010542398, -0.1810223872563199, 0.8939963497853238, 
    -0.41108877429515533, -0.17825977173502586, 1.0104585676265494, 4.263942410362829e-05, 8.070566317343375e-05, 
    3.517587553161026e-05, 0.9999999952155666, 1, 0, 0, 
    0, 1, 0, 0, 0, 
    1, 0, 0, 0, 0, 
    0, 0, 0, 0, 0]
>>>[0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 1.234, 0, 0, 0,
    0, 0, 4.37940971, 0, 0]
>>>[0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 1.234, 0, 0, 0,
    0, 0, 3.3162366617896932, 0, 0]
>>>[0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0.5201666932102418, 2.542004084904101e-07, 1.8254476433615963e-06, -0.9440115679003211,
    0.32991235149714754, 5.135843661699914e-07, -2.6176089039693444, 0.000201418779, -2.36140831e-07]
>>>([0, 0, 0, 0, 0,
     0, 0, 0, 0, 0,
     0, 0, 0, 0, 0,
     0, 0, 0, 0, 0,
     0, 0.5201666827717384, 2.5420041801788083e-07, 1.8254476420348602e-06, -0.9440115679003211,
     0.32991235149714754, 5.135843661699914e-07, -2.617608909188596, 0.0002014187785940397, -2.3614083149436963e-07],
     0)								#The last value is the check value, where 0 is OK.
```

##### 5.2.2 ASTROM parameters application

Use of these functions are appropriate when efficiency is important and where many star positions are to be transformed for one date.  The star-independent parameters can be obtained by calling one of the functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].

pymAtccq(RC, DC, PR, PD, PX, RV, ASTROM) : Quick transformation of a star's ICRS catalog entry (epoch J2000.0) into ICRS astrometric place, given precomputed star-independent astrometry parameters.

pymAtcc13(RC, DC, PR, PD, PX, RV, DATE1, DATE2) : Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS astrometric place.

pymAtciq(RC, DC, PR, PD, PX, RV, ASTROM) : Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed star-independent astrometry parameters.

pymAtci13(RC, DC, PR, PD, PX, RV, DATE1, DATE2) : Transform ICRS star data, epoch J2000.0, to CIRS.

pymAtciqn(RC, DC, PR, PD, PX, RV, ASTROM, N, B) : Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed star-independent astrometry parameters plus a list of light-deflecting bodies.

pymAtciqz(RC, DC, ASTROM) : Quick ICRS to CIRS transformation, given precomputed star-independent astrometry parameters, and assuming zero parallax and proper motion.

pymAticq(RI, DI, ASTROM) : Quick CIRS RA,Dec to ICRS astrometric place, given the star-independent astrometry parameters.

pymAtic13(RI, DI, DATE1, DATE2) : Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

pymAticqn(RI, DI, ASTROM, N, B) : Quick CIRS to ICRS astrometric place transformation, given the star-independent astrometry parameters plus a list of light-deflecting bodies.

pymAtoc13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL) : Observed place at a groundbased site to to ICRS astrometric RA,Dec. The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

pymAtco13(RC, DC, PR, PD, PX, RV, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL) : ICRS RA,Dec to observed place.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

pymAtoiq(TYPE, OB1, OB2, ASTROM) : Quick observed place to CIRS, given the star-independent astrometry parameters.

pymAtoi13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL) : Observed place to CIRS.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

pymAtioq(RI, DI, ASTROM) : Quick CIRS to observed place transformation.

pymAtio13(RI, DI, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL) : CIRS RA,Dec to observed place.  The caller supplies UTC, site coordinates, ambient air conditions and observing wavelength.

```python
import PyMsOfa as sf

RC = 2.71		#ICRS right ascension at J2000.0 (radians)
DC = 0.174		#ICRS declination at J2000.0 (radians)
PR = 1e-5		#RA proper motion (radians/year)
PD = 5e-6		#Dec proper motion (radians/year)
PX = 0.1		#parallax (arcsec)
RV = 55.0		#radial velocity (km/s, +ve if receding)
DATE1 = 2456165.5
DATE2 = 0.401182685
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
	#data for each of the n bodies
RI = 2.710121572969038991	#CIRS geocentric RA,Dec (radians)
DI = 0.1729371367218230438	#CIRS geocentric RA,Dec (radians)
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675		#UT1-UTC (seconds)
ELONG = -0.527800806		#longitude (radians, east +ve)
PHI = -1.2345856		#geodetic latitude (radians)
HM = 2738.0			#height above ellipsoid (m, geodetic)
XP = 2.47230737e-7		#polar motion coordinates (radians)
YP = 1.82640464e-6		#polar motion coordinates (radians)
PHPA = 731.0			#pressure at the observer (hPa = mB)
TC = 12.8			#ambient temperature at the observer (deg C)
RH = 0.59			#relative humidity at the observer (range 0-1)
WL = 0.55			#wavelength (micrometers)
OB1 = 2.710085107986886201	#observed Az, HA or RA (radians; Az is N=0,E=90)
OB2 = 0.1717653435758265198	#observed ZD or Dec (radians)
TYPE = 'R'			#type of coordinates - "R", "H" or "A" 

ASTROM=[0 for i in range(30)]
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtccq(RC, DC, PR, PD, PX, RV, ASTROM))
print(sf.pymAtcc13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))
ASTROM=[0 for i in range(30)]
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtciq(RC, DC, PR, PD, PX, RV, ASTROM))
print(sf.pymAtci13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))
ASTROM=[0 for i in range(30)]
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtciqn(RC, DC, PR, PD, PX, RV, ASTROM, N, B))
ASTROM=[0 for i in range(30)]
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtciqz(RC, DC, ASTROM))
ASTROM=[0 for i in range(30)]
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAticq(RI, DI, ASTROM))
print(sf.pymAtic13(RI, DI, DATE1, DATE2))
ASTROM=[0 for i in range(30)]
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAticqn(RI, DI, ASTROM, N, B))
print(sf.pymAtoc13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
print(sf.pymAtco13(RC, DC, PR, PD, PX, RV, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
ASTROM=[0 for i in range(30)]
ASTROM, J = sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM)
print(sf.pymAtoiq(TYPE, OB1, OB2, ASTROM))
print(sf.pymAtoi13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
ASTROM=[0 for i in range(30)]
ASTROM, J = sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM)
print(sf.pymAtioq(RI, DI, ASTROM))
print(sf.pymAtio13(RI, DI, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

>>>(2.7101265045313725, 0.17406325376283502)
>>>(2.7101265045313725, 0.17406325376283502)
>>>(2.7101215729686965, 0.17293713672195396)
>>>(2.7101215729686965, 0.17293713672195396, -0.002900618712657376)
>>>(2.710122008104983, 0.1729371916492768)
>>>(2.7099948992472567, 0.17287407209849318)
>>>(2.7101265045317167, 0.17406325376270346)
>>>(2.7101265045317167, 0.17406325376270346, -0.002900618712657376)
>>>(2.710126069960012, 0.174063198500879)
>>>(2.7093338203209707, 0.17505968077815165, 0)			#The last value is the check value, where 0 is OK.
>>>(5.348237013786584, 1.543277300890511, 0.9546737967220833, 
    0.17070036612674977, 2.7109026012689066, -0.0030205483548024123, 0)	#The last value is the check value, where 0 is OK.
>>>(2.7094994558424434, 0.17382752250461816)
>>>(2.7094994558424434, 0.17382752250461816, 0)			#The last value is the check value, where 0 is OK.
>>>(5.348098241468909, 1.5434235119717519, 0.9548493265835329, 
    0.17080380611316373, 2.710727071407457)
>>>(5.348098241468909, 1.5434235119717519, 0.9548493265835329, 
    0.17080380611316373, 2.710727071407457, 0)			#The last value is the check value, where 0 is OK.
```
