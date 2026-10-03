# Ariadne Algebra

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/algebra/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/algebra/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/win.yml)
[![Coverage Status](https://github.com/ariadne-cps/algebra/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/coverage.yml)

Ariadne Algebra is the standalone C++20 algebra layer used by Ariadne. It provides linear-algebra containers, power-series and expansion structures, differential-algebra types, graded algebras and related operations over the validated numeric types supplied by the lower Ariadne stack.

## Dependencies

Algebra depends directly on [ariadne-cps/interval](https://github.com/ariadne-cps/interval), included as a Git submodule. Interval supplies Numeric, Foundation and Utility transitively.

The intended repository chain is:

```text
utility -> foundation -> numeric -> interval -> algebra -> function
```

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/algebra.git
cd algebra
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel
ctest --output-on-failure
```

A C++20 compiler, CMake, GMP and MPFR are required.

## Python bindings

When a compatible Python development environment is available, the build produces a standalone `pyariadne` module containing the lower Interval/Numeric/Foundation bindings plus Algebra linear-algebra and differentiation bindings. The public `pyariadne-algebra` interface consumes only the direct lower `pyariadne-interval` interface.

## Tutorials

- C++: `tutorials/tutorial_algebra/tutorial_algebra.cpp`
- Python: `python/tutorials/tutorial_algebra.py`

Both cover linear algebra and differential algebra, the two main user-facing parts of this layer.

## License

Algebra is released under the GNU General Public License v3.0.
