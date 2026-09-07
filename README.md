# PyMsOfa

[![arXiv](https://img.shields.io/badge/arxiv-2310.08673-b31b1b.svg)](https://arxiv.org/abs/2310.08673) | [Paper](https://doi.org/10.1088/1674-4527/ad0499) | ![Python](https://img.shields.io/badge/Python-3.8%2B-green.svg) | ![License](https://img.shields.io/badge/license-MIT-blue.svg)

`PyMsOfa` is a Python package for the **Standards of Fundamental Astronomy
(SOFA)** service of the International Astronomical Union (IAU).  It implements
all **247** SOFA routines and is based on the SOFA release of **2023 October 11**.

**Version 2.x is a complete rewrite in pure Python / NumPy.**
There is no C source, no shared library, no compiler and no `ctypes` / `cffi`
layer any more: `pip install PyMsOfa` now works identically on Windows, Linux
and macOS, and the routines accept NumPy arrays as well as scalars.

## Installation

```bash
pip install PyMsOfa
```

The only dependency is NumPy.

## Quick start

Every routine and every SOFA constant is available directly from the top level
of the package:

```python
import PyMsOfa as sf

# --- time -------------------------------------------------------------
djm0, djm = sf.pymCal2jd(2003, 6, 1)          # Gregorian calendar -> MJD
print(djm0, djm)                              # 2400000.5 52791.0

sign, ihmsf = sf.pymD2tf(3, 0.5)              # days -> hours/min/sec
print(sign, ihmsf)                            # + [12, 0, 0, 0]

tai1, tai2 = sf.pymUtctai(2453750.5, 0.892100694)
tt1,  tt2  = sf.pymTaitt(tai1, tai2)

# --- vectors / spherical ---------------------------------------------
sf.pymS2c(0.3, 0.4)                           # -> array of shape (3,)
sf.pymS2c([0.3, 0.5], [0.4, 0.6])             # -> array of shape (2, 3)   <- new in 2.x

# --- astrometry -------------------------------------------------------
ri, di, eo = sf.pymAtci13(2.71, 0.174,
                          -354.45e-3, 595.35e-3, 164.99e-3, 0.0,
                          2456165.5, 0.401182685)

# --- constants (from the SOFA header sofam.h) -------------------------
print(sf.DAS2R, sf.DJ00, sf.DAU, sf.CMPS)
```

Thematic sub-modules are still available if you prefer to keep the namespaces
apart:

```python
from PyMsOfa import PyMsOfa_time             as t
from PyMsOfa import PyMsOfa_earth_attitude   as e
from PyMsOfa import PyMsOfa_astrometry       as a
from PyMsOfa import PyMsOfa_basic            as b
from PyMsOfa import sofa_const               as c
```

## Package layout

```
PyMsOfa/
├── __init__.py                    # flat API: all 247 routines + all constants
├── sofa_const.py                  # SOFA constants (sofam.h) + dint/dnint/dsign/gmax/gmin
├── PyMsOfa_basic.py               #  51 routines - vector/matrix & angle utilities
├── PyMsOfa_time.py                #  33 routines - calendars & time scales
├── PyMsOfa_earth_attitude.py      #  96 routines - precession, nutation, Earth rotation
└── PyMsOfa_astrometry.py          #  72 routines - astrometry, ephemerides, geodesy
```

## What changed from 1.x to 2.x

`2.0.0` is a **breaking** release.  The routine names are unchanged (`pymXxx`,
one per SOFA `iauXxx`), but the calling convention was modernised.

| | 1.1.6 | 2.0.0 |
|---|---|---|
| Implementation | 3 parallel back ends: `PyMsOfa.ctypes`, `PyMsOfa.cffi`, `PyMsOfa.python` | one pure Python / NumPy implementation |
| C library | `sofa_a.c` compiled into `libsofa_c`; `ctypes`/`cffi` back ends unusable on Windows from PyPI | none |
| Import | `from PyMsOfa import python as sf` | `import PyMsOfa as sf` |
| Argument names | upper case (`DATE1`, `IY`, `THETA`) | lower case, matching the SOFA C prototypes (`date1`, `iy`, `theta`) |
| Error reporting | trailing integer status `J` in the return tuple | `ValueError` / `Warning`-style exceptions; the status is no longer returned |
| Vectorisation | scalars only | scalars **and** NumPy arrays (broadcasting) for the vector/angle routines |
| `iauGc2gd` | `pymGC2GD` in the pure-python back end | `pymGc2gd` everywhere |
| Helper `_A` variants (`pymRz_A`, `pymS2c_A`, …) | present in the `ctypes`/`cffi` back ends | removed — the single routine handles both cases |
| `ASTROM` / `LDBODY` | `ctypes.Structure` subclasses | plain Python classes with the same field names, holding NumPy arrays |

### Migration examples

```python
# ---- 1.x -----------------------------------------------------------
from PyMsOfa import python as sf
djm0, djm, j = sf.pymCal2jd(2003, 6, 1)
if j < 0:
    raise ValueError("bad date")

# ---- 2.x -----------------------------------------------------------
import PyMsOfa as sf
djm0, djm = sf.pymCal2jd(2003, 6, 1)      # raises ValueError on a bad date
```

```python
# ---- 1.x: the ASTROM context had to be pre-allocated ---------------
astrom = sf.pymASTROM()
astrom, eo = sf.pymApci13(date1, date2, astrom)

# ---- 2.x -----------------------------------------------------------
astrom, eo = sf.pymApci13(date1, date2)
```

Routines whose return tuple lost its trailing status value: `pymAf2a`,
`pymApco13`, `pymApio13`, `pymAtco13`, `pymAtio13`, `pymAtoc13`, `pymAtoi13`,
`pymCal2jd`, `pymD2dtf`, `pymDat`, `pymDtf2d`, `pymEform`, `pymEpv00`,
`pymGc2gde`, `pymGd2gc`, `pymGd2gce`, `pymJd2cal`, `pymJdcalf`, `pymPlan94`,
`pymPmsafe`, `pymPvstar`, `pymStarpm`, `pymStarpv`, `pymTaitt`, `pymTaiut1`,
`pymTaiutc`, `pymTcbtdb`, `pymTcgtt`, `pymTdbtcb`, `pymTdbtt`, `pymTf2a`,
`pymTf2d`, `pymTpors`, `pymTporv`, `pymTpxes`, `pymTpxev`, `pymTttai`,
`pymTttcg`, `pymTttdb`, `pymTtut1`, `pymUt1tai`, `pymUt1tt`, `pymUt1utc`,
`pymUtctai`, `pymUtcut1`.

Routines that now raise instead of signalling: `pymAb`, `pymAf2a`, `pymC2s`,
`pymCal2jd`, `pymD2dtf`, `pymDat`, `pymDtf2d`, `pymEform`, `pymEpv00`,
`pymGc2gd`, `pymGc2gde`, `pymGd2gc`, `pymGd2gce`, `pymJd2cal`, `pymJdcalf`,
`pymPlan94`, `pymPmsafe`, `pymPvstar`, `pymStarpm`, `pymStarpv`, `pymTaiutc`,
`pymTf2a`, `pymTf2d`, `pymTpxes`, `pymTpxev`, `pymUtcut1`.

In version 2.0.0 **every** routine reports an invalid input by raising
`ValueError` — there are no trailing status codes, no `-1e9` sentinel values and
no `None` placeholders anywhere in the public API.

If you need the old behaviour, pin the previous release:

```bash
pip install "PyMsOfa==1.1.6"
```

## Applications

This package is suitable for the astrometric detection of habitable planets of
the Closeby Habitable Exoplanet Survey
([CHES](https://doi.org/10.1088/1674-4527/ac77e4)) mission, and for frontier
themes of black holes and dark matter related to astrometric calculations and
other fields.

## Licence and acknowledgement

`PyMsOfa` is **not** a part of the SOFA routines; it is an independent Python
implementation of the algorithms published by the IAU SOFA Board.  It is
neither distributed, supported nor endorsed by the International Astronomical
Union.

In addition to `PyMsOfa`'s MIT licence, any use of this module should comply
with [SOFA's licence and terms of use](http://www.iausofa.org/tandc.html).
Especially, but not exclusively, any published work or commercial product
including results achieved by using `PyMsOfa` shall acknowledge that the SOFA
algorithms were used to obtain those results.

## Citation

To cite `PyMsOfa` in publications use:

1. Ji, Jiang-Hui, Tan, Dong-jie, Bao, Chun-hui, Huang, Xiu-min, Hu, Shoucun,
   Dong, Yao, Wang, Su. 2023, *PyMsOfa: A Python Package for the Standards of
   Fundamental Astronomy (SOFA) Service*, Research in Astronomy and
   Astrophysics, 23, 125015,
   doi:[10.1088/1674-4527/ad0499](https://doi.org/10.1088/1674-4527/ad0499)

2. Ji, Jiang-Hui, Li, Hai-Tao, Zhang, Jun-Bo, Fang, Liang, Li, Dong, Wang, Su,
   Cao, Yang, Deng, Lei, Li, Bao-Quan, Xian, Hao, Gao, Xiao-Dong, Zhang, Ang,
   Li, Fei, Liu, Jia-Cheng, Qi, Zhao-Xiang, Jin, Sheng, Liu, Ya-Ning, Chen,
   Guo, Li, Ming-Tao, Dong, Yao, Zhu, Zi, and CHES Consortium. 2022, *CHES: A
   Space-borne Astrometric Mission for the Detection of Habitable Planets of
   the Nearby Solar-type Stars*, Research in Astronomy and Astrophysics, 22,
   072003,
   doi:[10.1088/1674-4527/ac77e4](https://doi.org/10.1088/1674-4527/ac77e4)
