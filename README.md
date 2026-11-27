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

## Eigenvalue and eigenvector benchmark

The `eigen` branch provides an optional benchmark of complete
`EigenDecomposition<P,PR>` calls, **including eigenvectors**, and in the
validated paradigm the certification attempt. Shared timing, warm-ups and
reporting use the auxiliary `ariadne-cps/benchmarking` gitlink at `ss/bmk`.
Benchmarking does not become a production dependency.

```bash
git submodule update --init --recursive
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-release --parallel --target benchmarks
./build-release/benchmarks/benchmark_eigen_decomposition --sizes 3,8,16 --bits 53,128 --iterations 3 --samples 7 --warmups 1 --summary median,mad,p95
./build-release/benchmarks/benchmark_eigen_decomposition --sizes 3,8 --bits 128 --iterations 3 --samples 7 --warmups 1 --csv > eigen-benchmark.csv
./build-release/benchmarks/benchmark_eigen_decomposition --sizes 3,8 --bits 128 --iterations 3 --samples 7 --warmups 1 --profile --csv > eigen-profile.csv
```

The Ubuntu, macOS and Windows Release CI jobs **only compile** the
`benchmarks` target; they never execute benchmark binaries. Run them
manually in a controlled Release environment for performance measurements
and A/B comparisons.

The deterministic workload catalogue includes dense 2x2, upper-triangular,
dense symmetric with separated disks, and dense symmetric with overlapping
original Gershgorin disks (Jacobi and validated preconditioning). Two additional
4x4 non-normal workloads use an exact rational matrix with a known simple
real spectrum {1,3,5,7}. The point-matrix workload is now available for
both ApproximateTag (shifted Householder QR) and ValidatedTag (rigorous
eigenpair certificate); the small interval-perturbation case is
Validated-only. Their fixed order four is independent of `--sizes`. Options:
`--sizes` (orders 3–64), `--bits` (MP precision), `--iterations`,
`--samples` (odd, >=3), `--warmups`, `--summary` and `--csv`.
DP53 and MP at the selected bit widths are measured with both paradigms.

An optional `--profile` switch instruments the eigendecomposition itself
and adds the **actual algorithm strategy**, candidate-search median time and
interval-verification median time (nanoseconds per call). These stage values
are grouped into disjoint regions; dispatch and minor bookkeeping are not
accounted for, so they need not sum to end-to-end time. The instrumentation
adds overhead; compare profiled builds against other profiled builds, not
against the normal uninstrumented mode. Warm-ups are excluded from all
statistics. The default CSV columns and normal timing mode are unchanged.

Display and CSV reports round all timing values to whole nanoseconds (and
variance to whole ns^2); internal samples and statistics retain their original
floating-point precision. Sub-nanosecond fractions are not displayed.

Results report nanoseconds per complete call and the semantic outcome.
A fast `INCONCLUSIVE` outcome is not equivalent to a successful validated
decomposition. Compare only identical workload, size, precision and outcome.
Matrix construction is excluded from the timer. This does not measure
eigenvalues alone.

## Python bindings

When a compatible Python development environment is available, the build produces a standalone `pyariadne` module containing the lower Interval/Numeric/Foundation bindings plus Algebra bindings for linear algebra, sparse expansions and series, polynomial representations, and differential algebra. The public `pyariadne-algebra` interface consumes only the direct lower `pyariadne-interval` interface.
