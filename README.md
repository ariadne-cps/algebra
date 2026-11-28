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

## Spectral analysis and invariant subspaces

Ariadne Algebra provides spectral algorithms for real matrices in double and
multiple precision, with both approximate computation and rigorous interval
validation. The APIs separate three mathematical questions: finding the
eigenvalue multiset, computing a Schur representation, and constructing a
complete eigenvector basis. Invariant-subspace operations additionally support
matrices that are not diagonalizable.

| API | Purpose |
| --- | --- |
| `Spectrum<P, PR>` | Computes the complete eigenvalue multiset, including real and complex eigenvalues and algebraic multiplicities, without requiring eigenvectors. |
| `SchurDecomposition<P, PR>` | Computes a real orthogonal Schur factorization, `A Q ≈ Q T`, with real 1×1 and conjugate-pair 2×2 diagonal blocks. |
| `EigenDecomposition<P, PR>` | Computes a complete real or complex eigenvector basis when diagonalization can be established, automatically choosing the appropriate representation. |
| `SchurDecompositionData<P, PR>` | Extracts real invariant subspaces, reorders Schur blocks, and selects complete spectral groups by block indices or real-axis-centered disks. |
| `certified_coordinate_invariant_subspace` | Certifies coordinate invariant subspaces for entire interval-matrix families using exact block-zero structure. |
| `certified_invariant_graph` | Certifies a fixed, possibly non-coordinate invariant graph over an interval family. |
| `certified_invariant_graph_near` | Uses interval Krawczyk inclusion and contraction to certify a locally unique invariant graph for each matrix in an interval family. |

### Eigenvalues and diagonalization

`Spectrum` returns the eigenvalue multiset independently of
diagonalizability. It exploits the matrix's structural decomposition into
strongly connected components, handles scalar and 2×2 blocks directly, and
uses the real Schur algorithm for dense approximate spectra. Its validated
path provides outward complex enclosures for supported structural blocks.

`EigenDecomposition` presents **one public interface** for real input
matrices with real or nonreal eigenvalues. It preserves specialized real
paths—including symmetric Jacobi, triangular structure, and real
Gershgorin/Krawczyk certification—and invokes complex eigenvector
construction only when a nonreal spectrum requires it. On success,
`EigenDecompositionResult<P, PR>` holds a tagged real or complex result;
`real_data()` and `complex_data()` access the corresponding factors.
The real path retains real matrices and arithmetic without automatically
promoting them to complex matrices.

These operations report `SUCCESS`, `INCONCLUSIVE`, or
`UNSUPPORTED_DIMENSION` through their respective status types.
`INCONCLUSIVE` means that the requested numerical result or rigorous
certificate was not established; it is not a claim of mathematical
nonexistence.

### Real Schur form and spectral subspaces

The real Schur implementation combines Householder Hessenberg reduction
with shifted QR iterations, including Francis implicit double shifts for
conjugate pairs. Its orthogonal factor and real quasi-upper-triangular
matrix provide a basis for working with invariant subspaces independently
of whether an eigenvector basis exists.

`SchurDecompositionData` supplies
`leading_invariant_subspace`, `reorder_blocks`,
`selected_invariant_subspace`, and
`invariant_subspace_in_real_disk`. These operations retain complete
1×1 or 2×2 Schur blocks and produce a real basis `U` with restriction
`B` satisfying `A U ≈ U B`. The real-axis-centered disk selection
keeps conjugate eigenvalue pairs together. Orthogonal block reordering
uses small real Sylvester equations and verifies the resulting factors.

### Validated invariant subspaces

For `FloatBounds<PR>` interval-matrix families, the invariant-subspace
facilities provide mathematical certificates rather than relying on small
floating-point residuals:

- **Coordinate subspaces:** prove `A[Sᶜ,S] = 0` uniformly over the
  interval family and return the corresponding exact coordinate basis
  and interval restriction.
- **Fixed invariant graphs:** verify the matrix Riccati identity for
  a supplied singleton graph `X`, with basis `[I;X]`, uniformly over
  the input family.
- **Locally unique invariant graphs:** verify strict interval
  Krawczyk inclusion and contraction for a supplied graph center,
  candidate box, and preconditioner. The resulting enclosure contains
  a locally unique graph for each matrix represented by the input
  family.

The verified graph operations reuse Ariadne Algebra's own
`validated_root_contractors.hpp` primitives and preserve the module
boundary with equation-solving. Their certificates explicitly distinguish
a fixed common subspace from a family of matrix-dependent invariant
subspaces.

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

An optional benchmark measures complete
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
