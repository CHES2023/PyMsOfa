# -*- coding: utf-8 -*-
"""Convenience alias module: ``PyMsOfa_astrometry`` re-exports ``PyMsOfa_astrometry_n``."""
try:
    from .PyMsOfa_astrometry_n import *          # noqa: F401,F403
except ImportError:
    from PyMsOfa_astrometry_n import *           # noqa: F401,F403
