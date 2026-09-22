# -*- coding: utf-8 -*-
"""
PyMsOfa_extension
=================

Pure Python / NumPy subroutines for the updated **IAU 2006J2 / IAU 2000AR26**
precession-nutation model.

The routines follow the style of the *PyMsOfa v2* package: function names
``pym*``, two-part Julian Date arguments ``(date1, date2)``, NumPy
vectorisation and SOFA-style docstrings.

References
----------
* Liu, J.-C. & Huang, C.-L. 2025, A&A, 703, L21
  ("The IAU 2006 precession quantities with an improved Earth's J2
  long-term variation") -- the IAU 2006J2 precession.
* Liu, J.-C. et al. 2026 (in prep.), "Precession-nutation quantities
  compatible with the IAU 2006J2 precession model" -- the IAU 2000AR26
  nutation and the complete set of operational parameters.
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
  precession pA), as used by the SOFA routines.
* CIP/CIO series coefficients are in **microarcsecond** (converted with
  ``1e-6 * DAS2R``); the nutation series uses 0.1 microarcsecond.
* The CIO locator ``s`` is evaluated through the compact ``s + XY/2``
  series (mirroring ``iauS06``) and corrected by ``-X*Y/2``.
* The equation of the origins is assembled as
  ``EO = [EO series] - dpsi(R26) * cos(epsA)``.

Routines
--------
pymXy06J2(date1, date2)        -> x, y        : CIP coordinates (radians)
pymS06J2(date1, date2, x, y)   -> s           : CIO locator (radians)
pymS06J2direct(date1, date2)   -> s           : CIO locator, direct series
pymXys06J2a(date1, date2)      -> x, y, s     : CIP + CIO locator
pymEo06J2a(date1, date2)       -> eo          : equation of the origins
pymNut00a(date1, date2)        -> dpsi, deps  : IAU 2000A nutation
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
except ImportError:                      # flat layout: both files in one folder
    from iau2006j2_data import *

__all__ = [
    "pymXy06J2", "pymS06J2", "pymS06J2direct", "pymXys06J2a", "pymEo06J2a",
    "pymNut00aR26", "pymOppolzer",
    "pymP06J2", "pymPfw06J2", "pymObl06J2",
    "pymPnm06J2a", "pymC2i06J2a",
]

# ---------------------------------------------------------------------------
# SOFA constants (sofa.h / sofam.h)
# ---------------------------------------------------------------------------
D2PI = 6.283185307179586476925287          # 2 pi
DAS2R = 4.848136811095359935899141e-6      # arcseconds -> radians
DJ00 = 2451545.0                           # Julian Date of J2000.0
DJC = 36525.0                              # days per Julian century
TURNAS = 1296000.0                         # arcseconds in a full circle

UAS2R = 1e-6 * DAS2R                       # microarcseconds -> radians

# ---------------------------------------------------------------------------
# Series tables: each block j is pre-split into (multipliers, S, C), where
# multipliers is (M, 14) and S, C are the sine/cosine coefficients (M,).
# ---------------------------------------------------------------------------
def _split(series_dict):
    return {j: (np.asarray(rows)[:, :14],
                np.asarray(rows)[:, 14],
                np.asarray(rows)[:, 15])
            for j, rows in series_dict.items()}


_X_SERIES = _split(X_SERIES)
_Y_SERIES = _split(Y_SERIES)
_S_SERIES = _split(S_SERIES)
_SP_SERIES = _split(SPLUSXY2_SERIES)
_EO_SERIES = _split(EO_SERIES)

_NUT_LS = np.asarray(NUT_LS, dtype=float)
_NUT_PL = np.asarray(NUT_PL, dtype=float)


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


# ---------------------------------------------------------------------------
# Fundamental arguments, IERS Conventions (2003)  (SOFA iauFal03 ... iauFapa03)
# ---------------------------------------------------------------------------
def pymFal03(t):
    """Mean anomaly of the Moon (radians)."""
    return np.fmod(485868.249036 +
                   t * (1717915923.2178 +
                   t * (31.8792 +
                   t * (0.051635 +
                   t * (-0.00024470)))), TURNAS) * DAS2R


def pymFalp03(t):
    """Mean anomaly of the Sun (radians)."""
    return np.fmod(1287104.793048 +
                   t * (129596581.0481 +
                   t * (-0.5532 +
                   t * (0.000136 +
                   t * (-0.00001149)))), TURNAS) * DAS2R


def pymFaf03(t):
    """Mean longitude of the Moon minus that of the ascending node (radians)."""
    return np.fmod(335779.526232 +
                   t * (1739527262.8478 +
                   t * (-12.7512 +
                   t * (-0.001037 +
                   t * (0.00000417)))), TURNAS) * DAS2R


def pymFad03(t):
    """Mean elongation of the Moon from the Sun (radians)."""
    return np.fmod(1072260.703692 +
                   t * (1602961601.2090 +
                   t * (-6.3706 +
                   t * (0.006593 +
                   t * (-0.00003169)))), TURNAS) * DAS2R


def pymFaom03(t):
    """Mean longitude of the ascending node of the Moon (radians)."""
    return np.fmod(450160.398036 +
                   t * (-6962890.5431 +
                   t * (7.4722 +
                   t * (0.007702 +
                   t * (-0.00005939)))), TURNAS) * DAS2R


def pymFame03(t):
    """Mean longitude of Mercury (radians)."""
    return np.fmod(4.402608842 + 2608.7903141574 * t, D2PI)


def pymFave03(t):
    """Mean longitude of Venus (radians)."""
    return np.fmod(3.176146697 + 1021.3285546211 * t, D2PI)


def pymFae03(t):
    """Mean longitude of the Earth (radians)."""
    return np.fmod(1.753470314 + 628.3075849991 * t, D2PI)


def pymFama03(t):
    """Mean longitude of Mars (radians)."""
    return np.fmod(6.203480913 + 334.0612426700 * t, D2PI)


def pymFaju03(t):
    """Mean longitude of Jupiter (radians)."""
    return np.fmod(0.599546497 + 52.9690962641 * t, D2PI)


def pymFasa03(t):
    """Mean longitude of Saturn (radians)."""
    return np.fmod(0.874016757 + 21.3299104960 * t, D2PI)


def pymFaur03(t):
    """Mean longitude of Uranus (radians)."""
    return np.fmod(5.481293872 + 7.4781598567 * t, D2PI)


def pymFane03(t):
    """Mean longitude of Neptune (radians)."""
    return np.fmod(5.311886287 + 3.8133035638 * t, D2PI)


def pymFapa03(t):
    """General accumulated precession in longitude (radians)."""
    return (0.024381750 + 0.00000538691 * t) * t


def pymFundArgs(t):
    """
    The 14 fundamental arguments used by the CIP/CIO series, in the order of
    the coefficient tables::

        [ l, l', F, D, Om, LMe, LVe, LE, LMa, LJ, LSa, LU, LNe, pA ]

    Parameters
    ----------
    t : float or ndarray
        Julian centuries since J2000.0 (TT).

    Returns
    -------
    ndarray, shape (14,) or (14, N)
        Fundamental arguments in radians.
    """
    return np.array([
        pymFal03(t), pymFalp03(t), pymFaf03(t), pymFad03(t), pymFaom03(t),
        pymFame03(t), pymFave03(t), pymFae03(t), pymFama03(t), pymFaju03(t),
        pymFasa03(t), pymFaur03(t), pymFane03(t), pymFapa03(t),
    ])


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _t(date1, date2):
    """Julian centuries since J2000.0 (TT) from a 2-part Julian Date."""
    return ((date1 - DJ00) + date2) / DJC


def _polyval(p, t):
    """Horner evaluation of the polynomial p = [c0, c1, ..., cn]."""
    val = 0.0
    for c in reversed(p):
        val = val * t + c
    return val


def _pn_series(t, poly, series, fa):
    """
    Evaluate a precession-nutation series, returning the value in
    microarcseconds.

    ``series`` maps the power j of t to a tuple ``(mult, S, C)`` where
    ``mult`` is an (M, 14) array of fundamental-argument multipliers and
    ``S``, ``C`` the sine/cosine coefficients.  Works for scalar ``t`` or a
    NumPy array of epochs.
    """
    val = _polyval(poly, t)
    for j in sorted(series):
        mult, S, C = series[j]
        arg = mult @ fa
        val = val + (S @ np.sin(arg) + C @ np.cos(arg)) * (t ** j)
    return val


# ---------------------------------------------------------------------------
# Rotation matrices and matrix helpers (SOFA iauIr / iauRx / iauRy / iauRz)
# ---------------------------------------------------------------------------
def pymIr():
    """Identity rotation matrix, shape (3, 3)."""
    return np.eye(3)


def pymRx(phi):
    """Rotation about the x-axis by angle phi (radians)."""
    s, c = np.sin(phi), np.cos(phi)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, c, s],
                     [0.0, -s, c]])


def pymRy(theta):
    """Rotation about the y-axis by angle theta (radians)."""
    s, c = np.sin(theta), np.cos(theta)
    return np.array([[c, 0.0, -s],
                     [0.0, 1.0, 0.0],
                     [s, 0.0, c]])


def pymRz(psi):
    """Rotation about the z-axis by angle psi (radians)."""
    s, c = np.sin(psi), np.cos(psi)
    return np.array([[c, s, 0.0],
                     [-s, c, 0.0],
                     [0.0, 0.0, 1.0]])


def pymFw2m(gamb, phib, psi, eps):
    """Rotation matrix from Fukushima-Williams angles:  R1(-eps).R3(-psi).R1(phib).R3(gamb)."""
    return pymRx(-eps) @ pymRz(-psi) @ pymRx(phib) @ pymRz(gamb)


def pymFw2xy(gamb, phib, psi, eps):
    """CIP X, Y from Fukushima-Williams bias-precession-nutation angles."""
    return pymBpn2xy(pymFw2m(gamb, phib, psi, eps))


def pymBpn2xy(rbpn):
    """Extract the CIP X, Y coordinates from a bias-precession-nutation matrix."""
    rbpn = np.asarray(rbpn)
    return rbpn[2, 0], rbpn[2, 1]


def pymC2ixys(x, y, s):
    """Form the celestial-to-intermediate matrix given CIP X, Y and CIO locator s."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)

    r2 = x * x + y * y
    e = np.where(r2 != 0.0, np.arctan2(y, x), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(r2 <= 1.0, np.arctan(np.sqrt(r2 / (1.0 - r2))),
                     np.pi / 2.0)

    # Rotation about z and y, broadcasting over any trailing batch dims.
    def _r3(a):
        c, sn = np.cos(a), np.sin(a)
        r = np.zeros(np.shape(a) + (3, 3))
        r[..., 0, 0] = c
        r[..., 0, 1] = sn
        r[..., 1, 0] = -sn
        r[..., 1, 1] = c
        r[..., 2, 2] = 1.0
        return r

    def _r2(a):
        c, sn = np.cos(a), np.sin(a)
        r = np.zeros(np.shape(a) + (3, 3))
        r[..., 0, 0] = c
        r[..., 0, 2] = -sn
        r[..., 1, 1] = 1.0
        r[..., 2, 0] = sn
        r[..., 2, 2] = c
        return r

    return _r3(-(e + s)) @ _r2(d) @ _r3(e)


def pymC2ixy(date1, date2, x, y):
    """Form the celestial-to-intermediate matrix given CIP X, Y (IAU 2006J2)."""
    return pymC2ixys(x, y, pymS06J2(date1, date2, x, y))


def pymEors(rnpb, s):
    """Equation of the origins, given the classical NPB matrix and CIO locator s."""
    rnpb = np.asarray(rnpb, dtype=float)
    x = rnpb[..., 2, 0]
    y = rnpb[..., 2, 1]
    z = rnpb[..., 2, 2]
    ax = x / (1.0 + z)
    xs = 1.0 - ax * x
    ys = -ax * y
    zs = -x
    p = rnpb[..., 0, 0] * xs + rnpb[..., 0, 1] * ys + rnpb[..., 0, 2] * zs
    q = rnpb[..., 1, 0] * xs + rnpb[..., 1, 1] * ys + rnpb[..., 1, 2] * zs
    return np.where((p != 0.0) | (q != 0.0), s - np.arctan2(q, p), s)


# ---------------------------------------------------------------------------
# IAU 2000A nutation  (SOFA iauNut00a)
# ---------------------------------------------------------------------------
def pymNut00a(date1, date2):
    """
    Nutation, IAU 2000A model (MHB2000).

    Parameters
    ----------
    date1, date2 : float
        TT as a 2-part Julian Date.

    Returns
    -------
    dpsi, deps : float
        Nutation in longitude and obliquity (radians).
    """
    t = _t(date1, date2)
    U2R = DAS2R / 1e7                     # 0.1 microarcsecond -> radian

    # Luni-solar fundamental arguments (IERS 2003).
    el = pymFal03(t)
    elp = pymFalp03(t)
    f = pymFaf03(t)
    d = pymFad03(t)
    om = pymFaom03(t)

    # Planetary fundamental arguments (MHB2000 for l, F, D, Om, Neptune).
    al = np.fmod(2.35555598 + 8328.6914269554 * t, D2PI)
    af = np.fmod(1.627905234 + 8433.466158131 * t, D2PI)
    ad = np.fmod(5.198466741 + 7771.3771468121 * t, D2PI)
    aom = np.fmod(2.18243920 - 33.757045 * t, D2PI)
    alme = pymFame03(t)
    alve = pymFave03(t)
    alea = pymFae03(t)
    alma = pymFama03(t)
    alju = pymFaju03(t)
    alsa = pymFasa03(t)
    alur = pymFaur03(t)
    alne = np.fmod(5.321159000 + 3.8127774000 * t, D2PI)
    apa = pymFapa03(t)

    # Luni-solar nutation.
    arg = np.fmod(_NUT_LS[:, 0] * el + _NUT_LS[:, 1] * elp +
                  _NUT_LS[:, 2] * f + _NUT_LS[:, 3] * d +
                  _NUT_LS[:, 4] * om, D2PI)
    sarg = np.sin(arg)
    carg = np.cos(arg)
    dpsils = ((_NUT_LS[:, 5] + _NUT_LS[:, 6] * t) @ sarg +
              _NUT_LS[:, 7] @ carg) * U2R
    depsls = ((_NUT_LS[:, 8] + _NUT_LS[:, 9] * t) @ carg +
              _NUT_LS[:, 10] @ sarg) * U2R

    # Planetary nutation.
    arg = np.fmod(_NUT_PL[:, 0] * al + _NUT_PL[:, 1] * af +
                  _NUT_PL[:, 2] * ad + _NUT_PL[:, 3] * aom +
                  _NUT_PL[:, 4] * alme + _NUT_PL[:, 5] * alve +
                  _NUT_PL[:, 6] * alea + _NUT_PL[:, 7] * alma +
                  _NUT_PL[:, 8] * alju + _NUT_PL[:, 9] * alsa +
                  _NUT_PL[:, 10] * alur + _NUT_PL[:, 11] * alne +
                  _NUT_PL[:, 12] * apa, D2PI)
    sarg = np.sin(arg)
    carg = np.cos(arg)
    dpsipl = (_NUT_PL[:, 13] @ sarg + _NUT_PL[:, 14] @ carg) * U2R
    depspl = (_NUT_PL[:, 15] @ sarg + _NUT_PL[:, 16] @ carg) * U2R

    return dpsils + dpsipl, depsls + depspl


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
    t = _t(date1, date2)
    ve = pymFave03(t)
    ea = pymFae03(t)
    ma = pymFama03(t)
    ju = pymFaju03(t)
    pa = pymFapa03(t)

    arg = (OPPOLZER_TERMS[:, 0] * ve + OPPOLZER_TERMS[:, 1] * ea +
           OPPOLZER_TERMS[:, 2] * ma + OPPOLZER_TERMS[:, 3] * ju +
           OPPOLZER_TERMS[:, 4] * pa)
    sarg = np.sin(arg)
    carg = np.cos(arg)

    dpsi = (OPPOLZER_TERMS[:, 5] @ sarg + OPPOLZER_TERMS[:, 6] @ carg) * UAS2R
    deps = (OPPOLZER_TERMS[:, 7] @ carg + OPPOLZER_TERMS[:, 8] @ sarg) * UAS2R
    return dpsi, deps


def pymNut00aR26(date1, date2):
    """
    Nutation, IAU 2000AR26 model (IAU 2000A nutation adjusted for the
    IAU 2006J2 precession).

    Following Liu et al. (2026, Eqs. 25, 32-34), the adjustments are::

        dpsi = dpsi(IAU2000A) + d_eps0(dpsi) + d_J2(dpsi) + d_Opz(dpsi)
        deps = deps(IAU2000A)           + d_J2(deps) + d_Opz(deps)

    where ``d_eps0`` is the obliquity-epsilon0 correction, ``d_J2`` the
    parabolic J2 (Poisson) correction and ``d_Opz`` the complete planetary
    Oppolzer terms.

    Returns
    -------
    dpsi, deps : float
        Nutation in longitude and obliquity (radians).
    """
    t = _t(date1, date2)
    dpsi00, deps00 = pymNut00a(date1, date2)

    om = pymFaom03(t)
    f = pymFaf03(t)
    d = pymFad03(t)
    lp = pymFalp03(t)

    # Eq. (25): obliquity-epsilon0 adjustment (uas).
    d_eps0_psi = -8.1 * np.sin(om) - 0.6 * np.sin(2 * f - 2 * d + 2 * om)

    # Eq. (32): parabolic J2 adjustment (uas).
    d_j2_psi = (1.2 * np.sin(om) * t +
                (-189.2 * np.sin(om) - 14.5 * np.sin(2 * f - 2 * d + 2 * om)
                 - 2.5 * np.sin(2 * f + 2 * om) + 2.3 * np.sin(2 * om)
                 + 1.6 * np.sin(lp)) * t * t)
    d_j2_eps = (-0.7 * np.cos(om) * t +
                (101.3 * np.cos(om) + 6.3 * np.cos(2 * f - 2 * d + 2 * om)
                 + 1.1 * np.cos(2 * f + 2 * om) - 1.0 * np.cos(2 * om)) * t * t)

    # Eq. (33): complete planetary Oppolzer terms.
    d_opz_psi, d_opz_eps = pymOppolzer(date1, date2)

    dpsi = dpsi00 + (d_eps0_psi + d_j2_psi) * UAS2R + d_opz_psi
    deps = deps00 + d_j2_eps * UAS2R + d_opz_eps
    return dpsi, deps


# ---------------------------------------------------------------------------
# IAU 2006J2 precession quantities
# ---------------------------------------------------------------------------
def pymObl06J2(date1, date2):
    """
    Mean obliquity of the ecliptic, IAU 2006J2 model (Liu et al. 2026, Eq. 21).

    Returns
    -------
    epsa : float
        Mean obliquity of date (radians).
    """
    t = _t(date1, date2)
    return _polyval([84381.406, -46.836734, -0.0001936, 0.00200004,
                     -0.000000602, 0.000000011], t) * DAS2R


def pymP06J2(date1, date2):
    """
    IAU 2006J2 precession quantities (Liu & Huang 2025, Eqs. 20-21).

    Returns
    -------
    psia, oma : float
        Precession in longitude and obliquity (radians), relative to the
        ecliptic of epoch.
    pa, epsa, chia : float
        General precession in longitude, mean obliquity of date and the
        precession of the ecliptic along the mean equator (radians).
    """
    t = _t(date1, date2)

    psia = _polyval([0.0, 5038.482041, -1.07182, 0.01754827,
                     0.000126577, -0.000000103], t) * DAS2R
    oma = _polyval([84381.406, -0.025754, 0.0512625, -0.0077249,
                    -0.000000245, 0.000000260], t) * DAS2R
    pa = _polyval([0.0, 5028.796900, 1.1125525, 0.0187702,
                   -0.000019662, -0.000000017], t) * DAS2R
    epsa = _polyval([84381.406, -46.836734, -0.0001936, 0.00200004,
                     -0.000000602, 0.000000011], t) * DAS2R
    chia = _polyval([0.0, 10.556240, -2.3813876, -0.00121400,
                     0.000159277, -0.000000087], t) * DAS2R
    return psia, oma, pa, epsa, chia


def pymPfw06J2(date1, date2):
    """
    Fukushima-Williams angles for frame bias and precession, IAU 2006J2
    model (Liu et al. 2026, Eq. 44).

    Returns
    -------
    gamb, phib, psib, epsa : float
        The four F-W angles (radians).  ``epsa`` is the mean obliquity.
    """
    t = _t(date1, date2)

    gamb = _polyval([-0.052928, 10.556239, 0.493244, -0.0003096,
                     -0.0000033116, 0.0000000013], t) * DAS2R
    phib = _polyval([84381.412819, -46.810980, 0.0511146, 0.0005299,
                     -0.0000003175, 0.0000000185], t) * DAS2R
    psib = _polyval([-0.041775, 5038.482019, 1.565603, 0.0185079,
                     -0.0000227596, -0.0000000164], t) * DAS2R
    epsa = pymObl06J2(date1, date2)

    return gamb, phib, psib, epsa


# ---------------------------------------------------------------------------
# IAU 2006J2 / IAU 2000AR26 CIP and CIO quantities
# ---------------------------------------------------------------------------
def _xy(t, fa):
    """CIP X, Y (radians) from the precomputed century t and arguments fa."""
    x = _pn_series(t, X_POLY, _X_SERIES, fa) * UAS2R
    y = _pn_series(t, Y_POLY, _Y_SERIES, fa) * UAS2R
    return x, y


def pymXy06J2(date1, date2):
    """
    X, Y coordinates of the celestial intermediate pole, IAU 2006J2 /
    IAU 2000AR26 model (series-based).

    Parameters
    ----------
    date1, date2 : float
        TT as a 2-part Julian Date.

    Returns
    -------
    x, y : float
        CIP X, Y coordinates (radians).
    """
    t = _t(date1, date2)
    return _xy(t, pymFundArgs(t))


def _s(t, fa, x, y):
    """CIO locator s (radians) via the s+XY/2 series, given X, Y."""
    spxy = _pn_series(t, SPLUSXY2_POLY, _SP_SERIES, fa) * UAS2R
    return spxy - x * y / 2.0


def pymS06J2(date1, date2, x, y):
    """
    The CIO locator s, given the CIP X, Y coordinates, IAU 2006J2 /
    IAU 2000AR26 model.

    The series is actually for ``s + X*Y/2`` (more compact than a direct
    series for s); the result is corrected by ``-X*Y/2``, exactly as in the
    SOFA routine ``iauS06``.

    Returns
    -------
    s : float
        CIO locator s (radians).
    """
    t = _t(date1, date2)
    return _s(t, pymFundArgs(t), x, y)


def pymS06J2direct(date1, date2):
    """
    The CIO locator s from its direct series (``s_coeff.txt``).

    This is an alternative to :func:`pymS06J2`; the two agree to within the
    truncation accuracy of the series (~0.1 uas).

    Returns
    -------
    s : float
        CIO locator s (radians).
    """
    t = _t(date1, date2)
    return _pn_series(t, S_POLY, _S_SERIES, pymFundArgs(t)) * UAS2R


def pymXys06J2a(date1, date2):
    """
    X, Y coordinates of the CIP and the CIO locator s, IAU 2006J2 /
    IAU 2000AR26 model.

    Returns
    -------
    x, y, s : float
        CIP coordinates and CIO locator (radians).
    """
    t = _t(date1, date2)
    fa = pymFundArgs(t)
    x, y = _xy(t, fa)
    return x, y, _s(t, fa, x, y)


def pymEo06J2a(date1, date2):
    """
    Equation of the origins, IAU 2006J2 / IAU 2000AR26 model.

    Following Liu et al. (2026, Eqs. 63-66), the EO is assembled as::

        EO = [EO series: polynomial + complementary terms] - dpsi * cos(epsA)

    where ``dpsi`` is the IAU 2000AR26 nutation in longitude and ``epsA``
    the mean obliquity of date.

    Returns
    -------
    eo : float
        Equation of the origins (radians).
    """
    t = _t(date1, date2)
    eo_series = _pn_series(t, EO_POLY, _EO_SERIES, pymFundArgs(t)) * UAS2R
    dpsi, _ = pymNut00aR26(date1, date2)
    return eo_series - dpsi * np.cos(pymObl06J2(date1, date2))


# ---------------------------------------------------------------------------
# Matrices
# ---------------------------------------------------------------------------
def pymPnm06J2a(date1, date2):
    """
    Form the classical bias-precession-nutation matrix, IAU 2006J2 /
    IAU 2000AR26 model.

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
    Form the celestial-to-intermediate matrix, IAU 2006J2 / IAU 2000AR26 model.

    Returns
    -------
    rc2i : ndarray, shape (3, 3)
        Celestial-to-intermediate matrix.
    """
    x, y, s = pymXys06J2a(date1, date2)
    return pymC2ixys(x, y, s)
