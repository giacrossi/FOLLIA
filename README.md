<a name="top"></a>

# FOLLIA

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-alpha-orange.svg)]()
[![Build Status](https://travis-ci.org/giacombum/FOLLIA.svg?branch=master)](https://travis-ci.org/giacombum/FOLLIA)
[![Coverage Status](https://img.shields.io/codecov/c/github/giacombum/FOLLIA.svg)](http://codecov.io/github/giacombum/FOLLIA?branch=master)
[![GitHub issues](https://img.shields.io/github/issues/giacombum/FOLLIA.svg)]()

### FOLLIA, Fortran Library for Lagrange InterpolAtion

- FOLLIA is a pure Fortran (KISS) library for computing Lagrange interpolations;
- FOLLIA is Fortran 2003+ standard compliant;
- FOLLIA is OOP designed;
- FOLLIA is a Free, Open Source Project.

#### Table of Contents

+ [What is FOLLIA?](#what-is-wenoof?)
	+ [What is Lagrange Interpolation?](#what-is-lagrange?)
+ [Main features](#main-features)
+ [Status](#status)
+ [Copyrights](#copyrights)
+ [Documentation](#documentation)

## What is FOLLIA?

FOLLIA is a Modern Fortran library for the evaluation of Lagrange coefficients and interpolation on uniform and non-uniform one-dimensional grids.

### What is Lagrange Interpolation?

Lagrange Interpolation is a data fitting technique based on [Lagrange polynonials](https://en.wikipedia.org/wiki/Lagrange_polynomial). For a given set of distinct points and numbers, the Lagrange polynomial is the polynomial of lowest degree that assumes at each point the corresponding value (i.e. the functions coincide at each point). The interpolating polynomial of the least degree is unique, however, and since it can be arrived at through multiple methods, referring to "the Lagrange polynomial" is perhaps not as correct as referring to "the Lagrange form" of that unique polynomial.

FOLLIA is designed to provide a KISS, Object Oriented Fortran API for computing Lagrange interpolation

Go to [Top](#top)

## Main features

FOLLIA is aimed to be a KISS-pure-Fortran library for computing Lagrange interpolation, it being:

+ [x] Pure Fortran implementation;
+ [x] KISS and user-friendly:
  + [x] simple API;
  + [ ] easy building and porting on heterogeneous architectures;
+ [ ] efficient:
  + [ ] high scalability on parallel architectures:
    + [ ] support for shared memory multi/many cores architecture;
    + [ ] support for distributed memory cluster;
    + [ ] support for GPGPU/accelerators device;
+ [ ] well documented:
  + [x] clear documentation of schemes implementations;
  + [x] complete API reference;
  + [ ] comprehensive wiki:
+ [ ] collaborative developed;
+ [x] FOSS licensed;

Any feature request is welcome.

Go to [Top](#top)

## Copyrights

FOLLIA is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to FOLLIA is welcome, feel free to select the license that best matches your soul!

### Externals libraries

FOLLIA uses some external libraries (placed into the *external* subdirectory of the root project) for the testing suite. These library maybe distributed under different licensing system with respect the FOLLIA one, please refer to their own licenses.

Go to [Top](#top)

Go to [Top](#top)
