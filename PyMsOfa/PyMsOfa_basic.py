# -*- coding: utf-8 -*-
"""Convenience alias module: ``PyMsOfa_basic`` re-exports ``PyMsOfa_basic_n``."""
try:
    from .PyMsOfa_basic_n import *          # noqa: F401,F403
except ImportError:
    from PyMsOfa_basic_n import *           # noqa: F401,F403
