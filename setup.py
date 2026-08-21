# -*- coding: utf-8 -*-
"""
Compatibility shim.

All packaging metadata lives in ``pyproject.toml`` (PEP 621).  This file only
exists so that the legacy commands still work, e.g.::

    python setup.py sdist bdist_wheel

The recommended, modern invocation is::

    python -m build
"""
from setuptools import setup

setup()
