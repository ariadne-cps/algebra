# Ariadne Algebra

[![Format](https://github.com/ariadne-cps/algebra/actions/workflows/format.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/format.yml)
[![macOS Status](https://github.com/ariadne-cps/algebra/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/algebra/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/algebra/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/algebra/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/algebra/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/algebra/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/algebra/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/algebra)

Ariadne Algebra is the standalone C++20 algebra layer used by Ariadne. It provides linear-algebra containers, sparse expansions, power series, power-basis and Chebyshev polynomial representations, differential-algebra types, graded algebras and related operations over the validated numeric types supplied by the lower Ariadne stack.
Although independently maintained and usable, this repository is also part of the [Ariadne repository structure](https://github.com/ariadne-cps/ariadne/blob/main/doc/structure.md).

## Features

- **Linear algebra:** matrices, vectors, decompositions and linear-system operations.
- **Expansions and series:** sparse expansions and power series.
- **Polynomials:** power-basis and Chebyshev representations.
- **Differential and graded algebra:** differential-algebra types and graded structures.
- **Spectral analysis:** `Spectrum` for real/complex eigenvalues, `SchurDecomposition` for real Schur forms and invariant subspaces, and a unified `EigenDecomposition` for real/complex eigenvectors.
- **Validated computation:** interval eigenpair certification and invariant-subspace certificates for coordinate spaces and invariant graphs, including local Krawczyk verification.

## Dependencies

Algebra depends directly on [ariadne-cps/interval](https://github.com/ariadne-cps/interval), included as a Git submodule. Interval supplies Numeric, Foundation and Utility transitively.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/algebra.git
cd algebra
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel
ctest --output-on-failure
```

A C++20 compiler and CMake are required.

## Eigenvalue benchmark

An optional benchmark measures complete `EigenDecomposition<P,PR>`
calls in double and multiple precision, including the outcome of
validated certification.

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-release --parallel --target benchmarks
./build-release/benchmarks/benchmark_eigen_decomposition --sizes 3,8,16 --bits 53,128 --iterations 3 --samples 7 --warmups 1
```

Use `--csv` for machine-readable results and `--profile` for algorithm
and verification timings. CI builds the benchmark but does not execute it.
See `benchmarks/benchmark_eigen_decomposition.cpp` for workloads and options.

## Python bindings

When a compatible Python development environment is available, the build produces a standalone `pyariadne` module containing the lower Interval/Numeric/Foundation bindings plus Algebra bindings for linear algebra, sparse expansions and series, polynomial representations, and differential algebra. The public `pyariadne-algebra` interface consumes only the direct lower `pyariadne-interval` interface.
