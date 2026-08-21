# Contributing to PyMsOfa

Thanks for contributing!  This is a short guide for building, testing and
releasing PyMsOfa 2.x.

## Repository layout

```
PyMsOfa/                 # the package (pure Python + NumPy, no C library)
docs/                    # English introduction (md) + Chinese introduction (pdf/tex)
examples/                # runnable usage example
pyproject.toml           # packaging metadata (PEP 621)
setup.py                 # thin setuptools shim
MANIFEST.in
LICENSE
README.md
```

## Building the package

```bash
python -m pip install build
python -m build
```

This produces `dist/PyMsOfa-2.0.0-py3-none-any.whl` (a pure-Python wheel, valid
on Windows / Linux / macOS) and `dist/PyMsOfa-2.0.0.tar.gz` (the source
distribution).  No C compiler is required.

Validate the artifacts:

```bash
python -m twine check dist/*
```

## Running the tests / smoke checks

```bash
python -m pip install -e .
python examples/example_usage.py
```

## Releasing to PyPI

```bash
# test first (versions on PyPI can never be reused)
python -m twine upload --repository testpypi dist/*

# production
python -m twine upload dist/*
```

Use a PyPI API token (`pypi-...`).  Before uploading, bump the version in **two
places** and keep them in sync:

1. `pyproject.toml`  →  `version = "2.0.0"`
2. `PyMsOfa/__init__.py`  →  `__version__ = "2.0.0"`

## API conventions (2.x)

* Every SOFA routine is `iauXxx` → `pymXxx`, and is exposed at the top level:
  `import PyMsOfa as sf; sf.pymCal2jd(...)`.
* Invalid input raises `ValueError` — there are no trailing status codes, no
  `-1e9` / `None` sentinels anywhere in the public API.
* Routines accept NumPy arrays as well as scalars.
* `pymASTROM` / `pymLDBODY` are ordinary Python classes.
