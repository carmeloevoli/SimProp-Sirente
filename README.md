# SimProp - beta version

[![Build Status](https://github.com/carmeloevoli/SimProp-Sirente/actions/workflows/ci.yml/badge.svg)](https://github.com/carmeloevoli/SimProp-Sirente/actions)
![GitHub](https://img.shields.io/github/license/carmeloevoli/SimProp-Sirente)
[![Coverage Status](https://coveralls.io/repos/github/carmeloevoli/SimProp-beta/badge.svg)](https://coveralls.io/github/carmeloevoli/SimProp-beta)

`SimProp` is a Monte Carlo simulation code of ultra-high energy cosmic ray propagation. `SimProp-Sirente` is based on `SimProp-v2r4` but still in its initial development.

The stable release (and older versions) of `SimProp` can be found at 
https://augeraq.sites.lngs.infn.it

Enquires about `SimProp` must be addressed to the SimProp development group at
SimProp-dev@aquila.infn.it

## Obtaining SIMPROP source

The preferred way is to clone the source from the git repository:

```sh
git clone https://github.com/carmeloevoli/SimProp-Sirente.git
```

The other method is to download a ZIP file from [the GitHub page](https://github.com/carmeloevoli/SimProp-Sirente) by clicking "Clone and download" and then "Download ZIP".

## Install

`SimProp` can be installed on GNU/Linux and macOS (OS X), while other operating systems have not been tested and are generally not supported.

### Dependencies and requirements

Required:
- **CMake** (`cmake`) is necessary to configure and build the source code (ver. 3.1+)
- **GCC** (`g++`) or **Clang** (`clang++`) are supported compilers; as modern C++ is employed (C++14), at least [GCC 6.1](https://gcc.gnu.org/projects/cxx-status.html#cxx14) / [Clang 3.4](https://clang.llvm.org/cxx_status.html) should be considered
- **GNU Scientific Library (GSL)** (`gsl`) is mandatory for numerical integration and for special functions

Optional:
- **Git** is needed if one wants to clone and keep in sync the source code from the git repository (**recommended**)
- **LCOV** is used with `gcov` (GCC) to generate the code coverage reports

Provided with the source:
- **Google Test** is employed as a framework for unit tests [v1.10.0](https://github.com/google/googletest)
- **PLOG** is a C++ logging library [v1.1.5](https://github.com/SergiusTheBest/plog)
- **NamedType** is a C++ library for the easy implementation of StrongTypes [master](https://github.com/joboccara/NamedType)

### Install on GNU/Linux

Required packages to build on RHEL/CentOS/Fedora systems:
```sh
dnf install git cmake g++ gsl-devel
```

To build it with cmake:
```sh
cd SimProp-Sirente
mkdir build
cd build
cmake -DENABLE_TESTING=On ..
make -j
make test
```

### Install on macOS (OS X)

To obtain the required and optional packages, one can use [Homebrew](https://brew.sh):
```sh
brew install cmake gsl 
brew install gcc # only if GCC is desired
brew link --overwrite gcc # same as above
```

To build it with cmake:
```sh
cd SimProp-Sirente
mkdir build
cd build
cmake -DENABLE_TESTING=On ..
make -j
make test
```

Or for GCC:
```sh
export GCC_BREW_PATH=$(brew --cellar gcc)/$(brew info --json gcc | jq -r '.[0].installed[0].version');
export CC=$GCC_BREW_PATH/bin/gcc-9
export CXX=$GCC_BREW_PATH/bin/g++-9
cmake -DENABLE_TESTING=On ..
make -j
make test
```
