[![GitHub release; latest by date](https://img.shields.io/github/v/release/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/releases)
[![Test Status](https://img.shields.io/github/actions/workflow/status/SETI/rms-starcat/run-tests.yml?branch=main)](https://github.com/SETI/rms-starcat/actions)
[![Documentation Status](https://readthedocs.org/projects/rms-starcat/badge/?version=latest)](https://rms-starcat.readthedocs.io/en/latest/?badge=latest)
[![Code coverage](https://img.shields.io/codecov/c/github/SETI/rms-starcat/main?logo=codecov)](https://codecov.io/gh/SETI/rms-starcat)
<br />
[![PyPI - Version](https://img.shields.io/pypi/v/rms-starcat)](https://pypi.org/project/rms-starcat)
[![PyPI - Format](https://img.shields.io/pypi/format/rms-starcat)](https://pypi.org/project/rms-starcat)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/rms-starcat)](https://pypi.org/project/rms-starcat)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/rms-starcat)](https://pypi.org/project/rms-starcat)
<br />
[![GitHub commits since latest release](https://img.shields.io/github/commits-since/SETI/rms-starcat/latest)](https://github.com/SETI/rms-starcat/commits/main/)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/m/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/commits/main/)
[![GitHub last commit](https://img.shields.io/github/last-commit/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/commits/main/)
<br />
[![Number of GitHub open issues](https://img.shields.io/github/issues-raw/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/issues)
[![Number of GitHub closed issues](https://img.shields.io/github/issues-closed-raw/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/issues)
[![Number of GitHub open pull requests](https://img.shields.io/github/issues-pr-raw/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/pulls)
[![Number of GitHub closed pull requests](https://img.shields.io/github/issues-pr-closed-raw/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/pulls)
<br />
![GitHub License](https://img.shields.io/github/license/SETI/rms-starcat)
[![Number of GitHub stars](https://img.shields.io/github/stars/SETI/rms-starcat)](https://github.com/SETI/rms-starcat/stargazers)
![GitHub forks](https://img.shields.io/github/forks/SETI/rms-starcat)

# Introduction

`starcat` is a set of routines for converting between NumPy floating point and complex
scalars/arrays and starcat-format single- and double-precision floats.

`starcat` is a product of the [PDS Ring-Moon Systems Node](https://pds-rings.seti.org).

# Installation

The `starcat` module is available via the `rms-starcat` package on PyPI and can be installed with:

```sh
pip install rms-starcat
```

# Getting Started

The `starcat` module provides two functions for converting *from* starcat-format floats:

- [`from_starcat32`](https://rms-starcat.readthedocs.io/en/latest/module.html#starcat.from_starcat32):
  Interpret a series of bytes or NumPy array as one or more starcat single-precision floats
  and convert them to a NumPy float or complex scalar or array.
- [`from_starcat64`](https://rms-starcat.readthedocs.io/en/latest/module.html#starcat.from_starcat64):
  Interpret a series of bytes NumPy array as one or more starcat double-precision floats and
  convert them to a NumPy float or complex scalar or array.

and two functions for converting *to* starcat-format floats::

- [`to_starcat32`](https://rms-starcat.readthedocs.io/en/latest/module.html#starcat.to_starcat32):
  Convert a NumPy float or complex scalar or array to a NumPy array containing the
  binary representation of starcat single-precision floats. Such an array can not be
  used for arithmetic operations since it is not in IEEE 754 format.
- [`to_starcat32_bytes`](https://rms-starcat.readthedocs.io/en/latest/module.html#starcat.to_starcat32_bytes):
  Convert a NumPy float or complex scalar or array to a Python `bytes` object containing
  the binary representation of starcat single-precision floats.

Note that there are no functions to convert a NumPy array to starcat double-precision format.

Details of each function are available in the [module documentation](https://rms-starcat.readthedocs.io/en/latest/module.html).

Basic operation is as follows:

```python
import starcat
b = starcat.to_starcat32([1., 2., 3.])
print(f'b = {b!r}')
ba = starcat.to_starcat32_bytes([1., 2., 3.])
print(f'ba = {ba!r}')
v = starcat.from_starcat32(b)
print(f'v = {v!r}')
va = starcat.from_starcat32(ba)
print(f'va = {va!r}')
```

yields:

```python
b = array([2.3138e-41, 2.3318e-41, 2.3407e-41], dtype=float32)
ba = b'\x80@\x00\x00\x00A\x00\x00@A\x00\x00'
v = array([1., 2., 3.], dtype=float32)
va = array([1., 2., 3.], dtype=float32)
```

As NASA data products stored as starcat-format floats are often provided in JPL's VICAR file
format, you may also be interested in the `rms-vicar` package
([documentation](https://rms-vicar.readthedocs.io/en/latest)).

# Contributing

Information on contributing to this package can be found in the
[Contributing Guide](https://github.com/SETI/rms-starcat/blob/main/CONTRIBUTING.md).

# Links

- [Documentation](https://rms-starcat.readthedocs.io)
- [Repository](https://github.com/SETI/rms-starcat)
- [Issue tracker](https://github.com/SETI/rms-starcat/issues)
- [PyPi](https://pypi.org/project/rms-starcat)

# Licensing

This code is licensed under the [Apache License v2.0](https://github.com/SETI/rms-starcat/blob/main/LICENSE).
