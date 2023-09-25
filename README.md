# PyMsOfa
[![arXiv](https://img.shields.io/badge/arxiv-2200.0000-b31b1b.svg)](https://arxiv.org/abs/2200.0000) | [Paper](https://www.raa-journal.org/) |![Python](https://img.shields.io/badge/Python-3.0-green.svg)

This package is a Python wrapper for the Standards of Fundamental Astronomy (SOFA) service of the International Astronomical Union (IAU). It implements all 247 functions in the SOFA service and is based on the latest version released on May 12, 2021.

This Python package can be suitable for the astrometric detection of habitable planets of the Closeby Habitable Exoplanet Survey (CHES) mission and for the frontier themes of black holes and dark matter related to astrometric calculations and other fields.

## Structure
```bash
├─ C
    ├─makefile
    ├─sofa.h
    ├─sofa_a.c
    ├─sofam.h
                 # All configurations live there
├─ Example
    ├─makefile               # Notebooks to get started with multimodal datasets and models
├─ eval.py                 # Eval script
├─ insall.sh               # Installation script for DeepViewAgg
├─ scripts                 # Some scripts to help manage the project
├─ torch_points3d
    ├─ core                # Core components
    ├─ datasets            # All code related to datasets
    ├─ metrics             # All metrics and trackers
    ├─ models              # All models
    ├─ modules             # Basic modules that can be used in a modular way
    ├─ utils               # Various utils
    └─ visualization       # Visualization
└─ train.py                # Main script to launch a training
```

## Installation

Put the files in the C folder in the same directory.

### For windows

Use the following command in CMD.
```
mingw32-make
```
If this fails, please install [TDM-GCC](https://jmeubank.github.io/tdm-gcc/download/).

Then put the newly generated libsofa_c.so file in the same directory as PyMsOfa.py in the Python folder.

### For macOS

Use the following command in Terminal.
```
make
```
Then put the newly generated libsofa_c.so file in the same directory as PyMsOfa.py in the Python folder.

### For Linux
Use the following command in Terminal.
```
make
```
And then
```
ln -s libsofa_c.so /usr/lib/libsofa_c.so
```

## Documentation

For package documentation please refer to the underlying SOFA documentation at:
[(http://www.iausofa.org/cookbooks.html)](http://www.iausofa.org/cookbooks.html)

## Description

PyMsOfa is a Python module for accessing the [International Astronomical Union](https://www.iau.org/)’s [SOFA library](http://www.iausofa.org/) from Python. SOFA (Standards of Fundamental Astronomy) is a set of algorithms and procedures that implement standard models used in fundamental astronomy.

PyMsOfa is not a port of SOFA routines but a wrapper around the SOFA C library. Thus, no calculations are made into the PyMsOfa software, they are all delegated to the underlying SOFA C library.

PyMsOfa is neither distributed, supported nor endorsed by the International Astronomical Union. In addition to PyMsOfa’s license, any use of this module should comply with [SOFA’s license and terms of use](http://www.iausofa.org/tandc.html). Especially, but not exclusively, any published work or commercial products including results achieved by using PyMsOfa shall acknowledge that the SOFA software was used to obtain those results.


## Reference
In case you use all or part of the present code, please include a citation to the following paper:

```
@article{ji2022ches,
	title={CHES: a space-borne astrometric mission for the detection of habitable planets of the nearby solar-type stars},
	author={Ji, Jiang-Hui and Li, Hai-Tao and Zhang, Jun-Bo and Fang, Liang and Li, Dong and Wang, Su and Cao, Yang and Deng, Lei and Li, Bao-Quan and Xian, Hao and others},
	journal={Research in Astronomy and Astrophysics},
	volume={22},
	number={7},
	pages={072003},
	year={2022},
	publisher={IOP Publishing}
}
```
