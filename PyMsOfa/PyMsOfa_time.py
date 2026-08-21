# -*- coding: utf-8 -*-
"""Convenience alias module: ``PyMsOfa_time`` re-exports ``PyMsOfa_time_n``."""
try:
    from .PyMsOfa_time_n import *          # noqa: F401,F403
except ImportError:
    from PyMsOfa_time_n import *           # noqa: F401,F403
