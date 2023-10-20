# PyMsOfa
[![arXiv](https://img.shields.io/badge/arxiv-2310.08673-b31b1b.svg)](https://arxiv.org/abs/2310.08673) | [Paper](https://doi.org/10.1088/1674-4527/ad0499) |![Python](https://img.shields.io/badge/Python-3.0-green.svg)

This package is a Python package for the Standards of Fundamental Astronomy (SOFA) service of the International Astronomical Union (IAU). It implements the python package PyMsOfa for SOFA service in three ways: 

(1) a python wrapper package based on a foreign function library for Python (ctypes), 

(2) a python wrapper package with the foreign function interface for Python calling C code (cffi), 

(3) a python package directly written in pure python codes from SOFA subroutines. 

It implements all 247 functions in the SOFA service and is based on the latest version released on Oct 11, 2023.

This Python package can be suitable for the astrometric detection of habitable planets of the Closeby Habitable Exoplanet Survey (CHES) mission and for the frontier themes of black holes and dark matter related to astrometric calculations and other fields.

## Structure
```
├─ C
    ├─makefile			#makefile
    ├─sofa.h			#header file
    ├─sofa_a.c			#integrated all 247 routines in SOFA service 
    ├─sofam.h			#header file
├─ Example
    ├─Application in CHES	               
    	├─Fig4.png		#Figure 4 in paper
    	├─Fig4.py		#Programs to generate Figure 4
    	├─input.xls		#Input parameters obtained from Gaia DR3
    ├─Precession-nutation.py	#Programs on precession-nutation
    ├─coordinate.py		#Programs on coordinate
    ├─time.py			#Programs on time
├─ cffi                         #Cffi version, run with the libsofa_c.so generated in the C folder
    ├─PyMsOfa.py		#Contains all 247 routines
    ├─PyMsOfa_astrometry.py	#Astrometry module of PyMsOfa.py
    ├─PyMsOfa_basic.py		#Basic module of PyMsOfa.py
    ├─PyMsOfa_earth_attitude.py	#Earth attitude module of PyMsOfa.py
    ├─PyMsOfa_t.py		#Test file
    ├─PyMsOfa_time.py           #Time module of PyMsOfa.py
├─ ctypes                       #Ctypes version, run with the libsofa_c.so generated in the C folder
    ├─PyMsOfa.py		#Contains all 247 routines
    ├─PyMsOfa_astrometry.py	#Astrometry module of PyMsOfa.py
    ├─PyMsOfa_basic.py		#Basic module of PyMsOfa.py
    ├─PyMsOfa_earth_attitude.py	#Earth attitude module of PyMsOfa.py
    ├─PyMsOfa_t.py		#Test file
    ├─PyMsOfa_time.py           #Time module of PyMsOfa.py
├─ Python			#Python version, pure python codes
    ├─PyMsOfa.py		#Contains all 247 routines
    ├─PyMsOfa_astrometry.py	#Astrometry module of PyMsOfa.py
    ├─PyMsOfa_basic.py		#Basic module of PyMsOfa.py
    ├─PyMsOfa_earth_attitude.py	#Earth attitude module of PyMsOfa.py
    ├─PyMsOfa_t.py		#Test file
    ├─PyMsOfa_time.py           #Time module of PyMsOfa.py
├─ LICENSE.md                 
├─ README.md
```

## Installation

Put the files in the C folder in the same directory.

### For windows

Use the following command in CMD.
```
mingw32-make
```
If this doesn't work, please install [TDM-GCC](https://jmeubank.github.io/tdm-gcc/download/).

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

PyMsOfa is not a part of SOFA routines but a Python package for the SOFA C library. Thus, no calculations are made into the PyMsOfa package based on ctypes and cffi interface, which are all delegated to the underlying SOFA C library.

PyMsOfa is neither distributed, supported nor endorsed by the International Astronomical Union. In addition to PyMsOfa’s license, any use of this module should comply with [SOFA’s license and terms of use](http://www.iausofa.org/tandc.html). Especially, but not exclusively, any published work or commercial products including results achieved by using PyMsOfa shall acknowledge that the SOFA software was used to obtain those results.


## Citation

To cite PyMsOfa in publications use:
> 1.    Ji, Jiang-Hui, Tan, Dong-jie, Bao, Chun-hui, Huang, Xiu-min, Hu, Shoucun, Dong, Yao, Wang, Su. 2023, PyMsOfa: A Python Package for the Standards of Fundamental Astronomy (SOFA) Service, Research in Astronomy and Astrophysics, https://doi.org/10.1088/1674-4527/ad0499
> 2.	Ji, Jiang-Hui, Li, Hai-Tao, Zhang, Jun-Bo, Fang, Liang, Li, Dong, Wang, Su, Cao, Yang, Deng, Lei, Li, Bao-Quan, Xian, Hao, Gao, Xiao-Dong, Zhang, Ang, Li, Fei, Liu, Jia-Cheng, Qi, Zhao-Xiang,  Jin, Sheng, Liu, Ya-Ning, Chen, Guo, Li, Ming-Tao, Dong, Yao, Zhu, Zi, and CHES Consortium. 2022, CHES: A Space-borne Astrometric Mission for the Detection of Habitable Planets of the Nearby Solar-type Stars, Research in Astronomy and Astrophysics, 22, 072003


A BibTeX entry for LaTeX users is
```bibtex
@ARTICLE{2023arXiv231008673J,
       author = {{Ji}, Jianghui and {Tan}, Dongjie and {Bao}, Chunhui and {Huang}, Xiumin and {Hu}, Shoucun and {Dong}, Yao and {Wang}, Su},
        title = "{PyMsOfa: A Python Package for the Standards of Fundamental Astronomy (SOFA) Service}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
         year = 2023,
        month = oct,
          eid = {arXiv:2310.08673},
        pages = {arXiv:2310.08673},
          doi = {10.48550/arXiv.2310.08673},
archivePrefix = {arXiv},
       eprint = {2310.08673},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv231008673J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@article{2022RAA....22g2003J,
       author = {{Ji}, Jiang-Hui and {Li}, Hai-Tao and {Zhang}, Jun-Bo and {Fang}, Liang and {Li}, Dong and {Wang}, Su and {Cao}, Yang and {Deng}, Lei and {Li}, Bao-Quan and {Xian}, Hao and {Gao}, Xiao-Dong and {Zhang}, Ang and {Li}, Fei and {Liu}, Jia-Cheng and {Qi}, Zhao-Xiang and {Jin}, Sheng and {Liu}, Ya-Ning and {Chen}, Guo and {Li}, Ming-Tao and {Dong}, Yao and {Zhu}, Zi and {CHES Consortium}},
        title = "{CHES: A Space-borne Astrometric Mission for the Detection of Habitable Planets of the Nearby Solar-type Stars}",
      journal = {Research in Astronomy and Astrophysics},
     keywords = {Astrometry and Celestial Mechanics, planets and satellites: detection, planets and satellites: terrestrial planets, stars: solar-type, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Astrophysics of Galaxies, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2022,
        month = jul,
       volume = {22},
       number = {7},
          eid = {072003},
        pages = {072003},
          doi = {10.1088/1674-4527/ac77e4},
archivePrefix = {arXiv},
       eprint = {2205.05645},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022RAA....22g2003J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
