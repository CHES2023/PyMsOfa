# -*- coding: utf-8 -*-
"""
PyMsOfa -- a Python package for the Standards of Fundamental Astronomy (SOFA)
service of the International Astronomical Union (IAU).

Version 2.x is a *pure Python / NumPy* implementation: no C library, no
compilation, no ``ctypes``/``cffi`` shared object.  All 247 SOFA routines are
available directly from the top level of the package::

    import PyMsOfa as sf

    sf.pymD2tf(3, 0.5)                     # scalar call
    sf.pymS2c([0.1, 0.2], [0.3, 0.4])      # vectorised call (NumPy arrays)
    print(sf.DAS2R, sf.DJ00)               # SOFA constants from sofam.h

Thematic sub-modules remain importable if you prefer them::

    from PyMsOfa import PyMsOfa_time as t

To cite PyMsOfa in publications use:

  Ji, Jiang-Hui, Tan, Dong-jie, Bao, Chun-hui, Huang, Xiu-min, Hu, Shoucun,
  Dong, Yao, Wang, Su. 2023, PyMsOfa: A Python Package for the Standards of
  Fundamental Astronomy (SOFA) Service, Research in Astronomy and Astrophysics,
  23, 125015, doi:10.1088/1674-4527/ad0499

NOTE FOR MAINTAINERS
--------------------
This file is hand-written.  When bumping the version, change BOTH
``__version__`` below and ``version`` in ``pyproject.toml``.
"""

__version__ = "2.0.0"
__author__ = "Ji, Jianghui"
__license__ = "MIT"

# SOFA release the algorithms were transcribed from
SOFA_RELEASE = "2023-10-11"

# --- thematic modules ------------------------------------------------------
from . import sofa_const                    # noqa: F401
from . import PyMsOfa_basic                 # noqa: F401
from . import PyMsOfa_time                  # noqa: F401
from . import PyMsOfa_earth_attitude        # noqa: F401
from . import PyMsOfa_astrometry            # noqa: F401

# --- flat API: every constant and every pym* routine at the top level ------
from .sofa_const import *                   # noqa: F401,F403
from .PyMsOfa_basic import *                # noqa: F401,F403
from .PyMsOfa_time import *                 # noqa: F401,F403
from .PyMsOfa_earth_attitude import *       # noqa: F401,F403
from .PyMsOfa_astrometry import *           # noqa: F401,F403


def _build_all():
    """Public API = SOFA routines + SOFA constants, without leaked helpers."""
    import types

    skip = {"np", "numpy", "types"}
    out = []
    for name, obj in list(globals().items()):
        if name.startswith("_") or name in skip:
            continue
        if isinstance(obj, types.ModuleType):
            continue
        out.append(name)
    return sorted(out)


__all__ = _build_all() + [
    "__version__",
    "SOFA_RELEASE",
    "sofa_const",
    "PyMsOfa_basic",
    "PyMsOfa_time",
    "PyMsOfa_earth_attitude",
    "PyMsOfa_astrometry",
]

del _build_all
