# -*- coding: utf-8 -*-
"""
PyMsOfa_extension
=================

Pure Python / NumPy subroutines for the updated **IAU 2006J2/2000AR26**
precession-nutation model, written as an *extension* of the PyMsOfa package.

The routines follow the style of the *PyMsOfa v2* package (function names
``pym*``, two-part Julian Date arguments ``(date1, date2)``, NumPy
vectorisation).  Everything that already exists in the base package --
SOFA constants, the 14 IERS 2003 fundamental arguments, the rotation
matrices and matrix helpers, and the IAU 2000A nutation ``pymNut00a`` -- is
imported from :mod:`sofa_const` and :mod:`PyMsOfa_earth_attitude` rather than
being redefined here.

References
----------
* Liu, J.-C. & Huang, C.-L. 2025, A&A, 703, L21
  ("The IAU 2006 precession quantities with an improved Earth's J2
  long-term variation") -- the IAU 2006J2 precession.
* Liu, J.-C. et al. 2026, "Precession-nutation quantities compatible with
  the IAU 2006J2 precession model" -- the IAU 2000AR26 nutation and the
  complete set of operational parameters.
* Ferrándiz, J. M., Navarro, J. F., Martínez-Belda, M. C., Escapa, A.,
  Getino, J. 2018, A&A, 618, A69 -- the complete planetary Oppolzer terms.

Data
----
The series coefficients (X, Y, s, s + XY/2, EO) are stored in the companion
data module :mod:`iau2006j2_data`, generated from ``X_coeff.txt``,
``Y_coeff.txt``, ``s_coeff.txt``, ``sPlusXY2_coeff.txt`` and ``EO_coeff.txt``.

Notes
-----
* Fundamental arguments are the 14 IERS Conventions (2003) arguments
  (l, l', F, D, Om and the eight planetary longitudes plus the general
  precession pA), imported from :mod:`PyMsOfa_earth_attitude`.
* CIP/CIO series coefficients are in **microarcsecond** (converted with
  ``1e-6 * DAS2R``).
* The CIO locator ``s`` is evaluated through the compact ``s + XY/2``
  series (mirroring ``iauS06``) and corrected by ``-X*Y/2``.

Routines
--------
pymXy06J2(date1, date2)        -> x, y        : CIP coordinates (radians)
pymS06J2(date1, date2, x, y)   -> s           : CIO locator (radians)
pymS06J2direct(date1, date2)   -> s           : CIO locator, direct series
pymXys06J2a(date1, date2)      -> x, y, s     : CIP + CIO locator
pymEo06J2a(date1, date2)       -> eo          : equation of the origins
pymNut00aR26(date1, date2)     -> dpsi, deps  : IAU 2000AR26 nutation
pymOppolzer(date1, date2)      -> dpsi, deps  : planetary Oppolzer terms
pymP06J2(date1, date2)         -> psia, oma, pa, epsa, chia : precession
pymPfw06J2(date1, date2)       -> gamb, phib, psib, epsa   : F-W angles
pymObl06J2(date1, date2)       -> epsa        : mean obliquity (IAU 2006J2)
pymPnm06J2a(date1, date2)      -> rnpb        : bias-precession-nutation matrix
pymC2i06J2a(date1, date2)      -> rc2i        : celestial-to-intermediate matrix
"""

import numpy as np

try:                                     # installed inside the package
    from .iau2006j2_data import *
    from .sofa_const import DAS2R, DJ00, DJC
    from .PyMsOfa_earth_attitude import (
        pymFal03, pymFalp03, pymFaf03, pymFad03, pymFaom03,
        pymFame03, pymFave03, pymFae03, pymFama03, pymFaju03,
        pymFasa03, pymFaur03, pymFane03, pymFapa03,
        pymNut00a, pymFw2m, pymBpn2xy, pymC2ixys,
    )
except ImportError:                      # flat layout: all modules in one folder
    from iau2006j2_data import *
    from sofa_const import DAS2R, DJ00, DJC
    from PyMsOfa_earth_attitude import (
        pymFal03, pymFalp03, pymFaf03, pymFad03, pymFaom03,
        pymFame03, pymFave03, pymFae03, pymFama03, pymFaju03,
        pymFasa03, pymFaur03, pymFane03, pymFapa03,
        pymNut00a, pymFw2m, pymBpn2xy, pymC2ixys,
    )

__all__ = [
    "pymXy06J2", "pymS06J2", "pymS06J2direct", "pymXys06J2a", "pymEo06J2a",
    "pymNut00aR26", "pymOppolzer",
    "pymP06J2", "pymPfw06J2", "pymObl06J2",
    "pymPnm06J2a", "pymC2i06J2a",
]

UAS2R = 1e-6 * DAS2R                       # microarcseconds -> radians

# ---------------------------------------------------------------------------
# Planetary Oppolzer terms for the Earth's figure axis (complete table).
#
# Ferrándiz et al. 2018, A&A 618, A69, Table 1.  Coefficients (in uas) of the
# largest Oppolzer terms of planetary origin.  Each row is
#   [nVe, nE, nMa, nJ, nPa,  dpsi_sin, dpsi_cos, deps_cos, deps_sin]
# and the argument is  nVe*LVe + nE*LE + nMa*LMa + nJ*LJ + nPa*pA
# (the remaining fundamental arguments have zero coefficients).
# dpsi = dpsi_sin*sin(arg) + dpsi_cos*cos(arg)
# deps = deps_cos*cos(arg) + deps_sin*sin(arg)
# ---------------------------------------------------------------------------
OPPOLZER_TERMS = np.array([
    # nVe  nE  nMa  nJ  nPa   dpsi_sin dpsi_cos deps_cos deps_sin
    [  0,   1,   0,  -1,   0,     1,      0,       0,       0],  # LE-LJ          (ind. Moon)
    [  2,  -4,   0,   0,   2,    -4,      0,       2,       0],  # 2LVe-4LE+2pA   (ind. Sun)
    [  0,   1,   0,  -2,   0,     3,    -35,       1,      14],  # LE-2LJ         (ind. Sun)
    [  0,   3,  -4,   0,   0,    -8,      1,      -3,       0],  # 3LE-4LMa       (ind. Sun)
    [  3,  -4,   0,   0,   0,    -7,     21,      -3,      -8],  # 3LVe-4LE       (ind. Sun)
    [  0,   1,   0,  -1,   0,   -47,     -8,     -17,       3],  # LE-LJ          (ind. Sun)
    [  0,   2,  -2,   0,   0,   -11,     -1,      -4,       0],  # 2LE-2LMa       (ind. Sun)
    [  2,  -4,   0,   0,   2,   -14,     -1,       6,      -1],  # 2LVe-4LE+2pA   (dir. Venus)
    [  0,   1,   0,  -1,   0,     4,      1,       1,       0],  # LE-LJ          (dir. Jupiter)
])

# Series tables as NumPy arrays: each row is [14 multipliers, S, C].
_X_SERIES = {j: np.asarray(r, dtype=float) for j, r in X_SERIES.items()}
_Y_SERIES = {j: np.asarray(r, dtype=float) for j, r in Y_SERIES.items()}
_S_SERIES = {j: np.asarray(r, dtype=float) for j, r in S_SERIES.items()}
_SP_SERIES = {j: np.asarray(r, dtype=float) for j, r in SPLUSXY2_SERIES.items()}
_EO_SERIES = {j: np.asarray(r, dtype=float) for j, r in EO_SERIES.items()}


def pymXy06J2(date1, date2):
    """
    X, Y coordinates of the celestial intermediate pole, IAU 2006J2/2000AR26
    model (series-based).

    Parameters
    ----------
    date1, date2 : float
        TT as a 2-part Julian Date.

    Returns
    -------
    x, y : float
        CIP X, Y coordinates (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    fa = np.array([pymFal03(t), pymFalp03(t), pymFaf03(t), pymFad03(t),
                   pymFaom03(t), pymFame03(t), pymFave03(t), pymFae03(t),
                   pymFama03(t), pymFaju03(t), pymFasa03(t), pymFaur03(t),
                   pymFane03(t), pymFapa03(t)])

    # X = polynomial + Poisson series (microarcsecond), then to radians
    x = 0.0
    for c in reversed(X_POLY):
        x = x * t + c
    for j, r in _X_SERIES.items():
        a = r[:, :14] @ fa
        x += np.sum(r[:, 14] * np.sin(a) + r[:, 15] * np.cos(a)) * (t ** j)

    # Y = polynomial + Poisson series (microarcsecond), then to radians
    y = 0.0
    for c in reversed(Y_POLY):
        y = y * t + c
    for j, r in _Y_SERIES.items():
        a = r[:, :14] @ fa
        y += np.sum(r[:, 14] * np.sin(a) + r[:, 15] * np.cos(a)) * (t ** j)

    return x * UAS2R, y * UAS2R


def pymS06J2(date1, date2, x, y):
    """
    The CIO locator s, given the CIP X, Y coordinates, IAU 2006J2/2000AR26
    model.  The series is for ``s + X*Y/2`` (mirroring ``iauS06``); the result
    is corrected by ``-X*Y/2``.

    Returns
    -------
    s : float
        CIO locator s (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    fa = np.array([pymFal03(t), pymFalp03(t), pymFaf03(t), pymFad03(t),
                   pymFaom03(t), pymFame03(t), pymFave03(t), pymFae03(t),
                   pymFama03(t), pymFaju03(t), pymFasa03(t), pymFaur03(t),
                   pymFane03(t), pymFapa03(t)])

    s = 0.0
    for c in reversed(SPLUSXY2_POLY):
        s = s * t + c
    for j, r in _SP_SERIES.items():
        a = r[:, :14] @ fa
        s += np.sum(r[:, 14] * np.sin(a) + r[:, 15] * np.cos(a)) * (t ** j)

    return s * UAS2R - x * y / 2.0


def pymS06J2direct(date1, date2):
    """
    The CIO locator s from its direct series (``s_coeff.txt``).

    Returns
    -------
    s : float
        CIO locator s (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    fa = np.array([pymFal03(t), pymFalp03(t), pymFaf03(t), pymFad03(t),
                   pymFaom03(t), pymFame03(t), pymFave03(t), pymFae03(t),
                   pymFama03(t), pymFaju03(t), pymFasa03(t), pymFaur03(t),
                   pymFane03(t), pymFapa03(t)])

    s = 0.0
    for c in reversed(S_POLY):
        s = s * t + c
    for j, r in _S_SERIES.items():
        a = r[:, :14] @ fa
        s += np.sum(r[:, 14] * np.sin(a) + r[:, 15] * np.cos(a)) * (t ** j)

    return s * UAS2R


def pymXys06J2a(date1, date2):
    """
    X, Y coordinates of the CIP and the CIO locator s, IAU 2006J2/2000AR26 model.

    Returns
    -------
    x, y, s : float
        CIP coordinates and CIO locator (radians).
    """
    rnpb = pymPnm06J2a(date1, date2)
    x, y = pymBpn2xy(rnpb)
    s = pymS06J2(date1, date2, x, y)
    return x, y, s


def pymEo06J2a(date1, date2):
    """
    Equation of the origins, IAU 2006J2/2000AR26 model.  Assembled as
    ``EO = [EO series] - dpsi(R26) * cos(epsA)``.

    Returns
    -------
    eo : float
        Equation of the origins (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    fa = np.array([pymFal03(t), pymFalp03(t), pymFaf03(t), pymFad03(t),
                   pymFaom03(t), pymFame03(t), pymFave03(t), pymFae03(t),
                   pymFama03(t), pymFaju03(t), pymFasa03(t), pymFaur03(t),
                   pymFane03(t), pymFapa03(t)])

    eo = 0.0
    for c in reversed(EO_POLY):
        eo = eo * t + c
    for j, r in _EO_SERIES.items():
        a = r[:, :14] @ fa
        eo += np.sum(r[:, 14] * np.sin(a) + r[:, 15] * np.cos(a)) * (t ** j)

    dpsi, _ = pymNut00aR26(date1, date2)
    return eo * UAS2R - dpsi * np.cos(pymObl06J2(date1, date2))


def pymOppolzer(date1, date2):
    """
    Planetary Oppolzer terms for the Earth's figure axis (complete table of
    Ferrándiz et al. 2018, A&A 618, A69).

    Returns
    -------
    dpsi, deps : float
        Oppolzer corrections to the nutation in longitude and obliquity
        (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    ve = pymFave03(t)
    ea = pymFae03(t)
    ma = pymFama03(t)
    ju = pymFaju03(t)
    pa = pymFapa03(t)
    arg = (OPPOLZER_TERMS[:, 0] * ve + OPPOLZER_TERMS[:, 1] * ea +
           OPPOLZER_TERMS[:, 2] * ma + OPPOLZER_TERMS[:, 3] * ju +
           OPPOLZER_TERMS[:, 4] * pa)

    dpsi = (OPPOLZER_TERMS[:, 5] @ np.sin(arg) +
            OPPOLZER_TERMS[:, 6] @ np.cos(arg)) * UAS2R
    deps = (OPPOLZER_TERMS[:, 7] @ np.cos(arg) +
            OPPOLZER_TERMS[:, 8] @ np.sin(arg)) * UAS2R
    return dpsi, deps


def pymNut00aR26(date1, date2):
    """
    Nutation, IAU 2000AR26 model (IAU 2000A nutation adjusted for the
    IAU 2006J2 precession)::

        dpsi = dpsi(IAU2000A) + d_eps0(dpsi) + d_J2(dpsi) + d_Opz(dpsi)
        deps = deps(IAU2000A)           + d_J2(deps) + d_Opz(deps)

    Returns
    -------
    dpsi, deps : float
        Nutation in longitude and obliquity (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    dpsi00, deps00 = pymNut00a(date1, date2)

    om = pymFaom03(t)
    f = pymFaf03(t)
    d = pymFad03(t)
    lp = pymFalp03(t)

    # Eq. (25): obliquity-epsilon0 adjustment (uas)
    d_eps0_psi = -8.1 * np.sin(om) - 0.6 * np.sin(2 * f - 2 * d + 2 * om)

    # Eq. (32): parabolic J2 adjustment (uas)
    d_j2_psi = (1.2 * np.sin(om) * t +
                (-189.2 * np.sin(om) - 14.5 * np.sin(2 * f - 2 * d + 2 * om)
                 - 2.5 * np.sin(2 * f + 2 * om) + 2.3 * np.sin(2 * om)
                 + 1.6 * np.sin(lp)) * t * t)
    d_j2_eps = (-0.7 * np.cos(om) * t +
                (101.3 * np.cos(om) + 6.3 * np.cos(2 * f - 2 * d + 2 * om)
                 + 1.1 * np.cos(2 * f + 2 * om) - 1.0 * np.cos(2 * om)) * t * t)

    # Eq. (33): complete planetary Oppolzer terms
    d_opz_psi, d_opz_eps = pymOppolzer(date1, date2)

    dpsi = dpsi00 + (d_eps0_psi + d_j2_psi) * UAS2R + d_opz_psi
    deps = deps00 + d_j2_eps * UAS2R + d_opz_eps
    return dpsi, deps


def pymObl06J2(date1, date2):
    """
    Mean obliquity of the ecliptic, IAU 2006J2 model (Liu et al. 2026, Eq. 21).

    Returns
    -------
    epsa : float
        Mean obliquity of date (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC
    epsa = (84381.406 +
           (-46.836734 +
           (-0.0001936 +
           ( 0.00200004 +
           (-0.000000602 +
           ( 0.000000011) * t) * t) * t) * t) * t) * DAS2R
    return epsa


def pymP06J2(date1, date2):
    """
    IAU 2006J2 precession quantities (Liu & Huang 2025, Eqs. 20-21).

    Returns
    -------
    psia, oma : float
        Precession in longitude and obliquity (radians).
    pa, epsa, chia : float
        General precession, mean obliquity and ecliptic precession (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC

    psia = (5038.482041 +
           (-1.07182 +
           (0.01754827 +
           (0.000126577 +
           (-0.000000103) * t) * t) * t) * t) * t * DAS2R

    oma = (84381.406 +
          (-0.025754 +
          (0.0512625 +
          (-0.0077249 +
          (-0.000000245 +
          (0.000000260) * t) * t) * t) * t) * t) * DAS2R

    pa = (5028.796900 +
         (1.1125525 +
         (0.0187702 +
         (-0.000019662 +
         (-0.000000017) * t) * t) * t) * t) * t * DAS2R

    epsa = pymObl06J2(date1, date2)

    chia = (10.556240 +
           (-2.3813876 +
           (-0.00121400 +
           (0.000159277 +
           (-0.000000087) * t) * t) * t) * t) * t * DAS2R

    return psia, oma, pa, epsa, chia


def pymPfw06J2(date1, date2):
    """
    Fukushima-Williams angles for frame bias and precession, IAU 2006J2
    model (Liu et al. 2026, Eq. 44).

    Returns
    -------
    gamb, phib, psib, epsa : float
        The four F-W angles (radians).
    """
    t = ((date1 - DJ00) + date2) / DJC

    gamb = (-0.052928 +
           (10.556239 +
           (0.493244 +
           (-0.0003096 +
           (-0.0000033116 +
           (0.0000000013) * t) * t) * t) * t) * t) * DAS2R

    phib = (84381.412819 +
           (-46.810980 +
           (0.0511146 +
           (0.0005299 +
           (-0.0000003175 +
           (0.0000000185) * t) * t) * t) * t) * t) * DAS2R

    psib = (-0.041775 +
           (5038.482019 +
           (1.565603 +
           (0.0185079 +
           (-0.0000227596 +
           (-0.0000000164) * t) * t) * t) * t) * t) * DAS2R

    epsa = pymObl06J2(date1, date2)

    return gamb, phib, psib, epsa


def pymPnm06J2a(date1, date2):
    """
    Form the classical bias-precession-nutation matrix, IAU 2006J2/2000AR26 model.

    Returns
    -------
    rnpb : ndarray, shape (3, 3)
        Bias-precession-nutation matrix.
    """
    gamb, phib, psib, epsa = pymPfw06J2(date1, date2)
    dpsi, deps = pymNut00aR26(date1, date2)
    return pymFw2m(gamb, phib, psib + dpsi, epsa + deps)


def pymC2i06J2a(date1, date2):
    """
    Form the celestial-to-intermediate matrix, IAU 2006J2/2000AR26 model.

    Returns
    -------
    rc2i : ndarray, shape (3, 3)
        Celestial-to-intermediate matrix.
    """
    rbpn = pymPnm06J2a(date1, date2)
    x, y = pymBpn2xy(rbpn)
    s = pymS06J2(date1, date2, x, y)
    return pymC2ixys(x, y, s)
