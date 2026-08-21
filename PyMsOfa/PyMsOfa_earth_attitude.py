# -*- coding: utf-8 -*-
"""Convenience alias module: ``PyMsOfa_earth_attitude`` re-exports ``PyMsOfa_earth_attitude_n``."""
try:
    from .PyMsOfa_earth_attitude_n import *          # noqa: F401,F403
except ImportError:
    from PyMsOfa_earth_attitude_n import *           # noqa: F401,F403
