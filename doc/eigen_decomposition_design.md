# Eigenvalue and eigenvector computation — design journal

## Scope and workflow

Development branch: `eigen` in `ariadne-cps/algebra`, branched from `main` at `1b85f7e864044b697be2408c3575181a75f36d4b`.
Do not merge into `main` until implementation, numerical and semantic tests, and review are complete.

**Journal rule:** Every development commit must include an update to this document describing its design decisions, code changes, tests run, test results and outstanding issues. The first commit creates this journal. No successful build or test is to be claimed without running it.

## Purpose

Implement reusable approximate and validated spectral methods in Algebra. Evolution's first-order PDE solver currently takes user-supplied characteristic eigenvalues/eigenvectors. Algebra must remain independent of PDE assumptions; integration with Evolution is a later project.

## Public API direction (first restricted implementation)

- `EigenDecomposition<P, PR>`: algorithm/solver object, with `P` in `ApproximateTag` or `ValidatedTag` and `PR` denoting numeric precision.
- `EigenDecompositionResult<P, PR>`: corresponding result type. `ApproximateTag` returns approximations without certification. `ValidatedTag` provides mathematically guaranteed enclosures **only when certified**, and an explicit inconclusive status otherwise.
- Support real and complex eigenvalues eventually; do not suggest every real matrix has a real eigenbasis.
- A failed certificate must never be reported as proof that no decomposition exists.
- Prefer an optional validated result over a success-looking empty or invalid enclosure.
- Input matrices should use the existing Ariadne number/matrix conventions. The exact call signatures, storage of complex outputs, convergence controls and precision traits remain to be decided after checking available primitives.

## Mathematical contracts

For a real square matrix `M`, a diagonalization `M T = T D` requires `T` invertible. Some matrices are defective; some have nonreal eigenvalues. A general eigensolver should ideally expose a Schur-based result even if diagonalization fails.

For the first increment, target **simple real eigenvalues and their normalized eigenvectors** for orders 1 and 2 only. This is a deliberate limited-scope milestone, not a general eigenvalue decomposition.

Approximate stage:
1. Current 2x2 milestone uses the closed-form quadratic discriminant, not a general QR algorithm. A stable general backend for larger matrices remains to be designed.
2. Measure residuals and account for convergence and ill-conditioning, without treating small numerical residuals as proof.

Validated stage:
1. For the present 2x2 case, apply directly validated interval arithmetic to the exact quadratic formula; a general validated eigenpair solver seeded by approximate candidates is future work.
2. Fix eigenvector scaling with an explicit normalization constraint.
3. Apply an interval inclusion technique (e.g. interval Newton / Krawczyk) to the eigenpair equations `Mv - lambda v = 0` and the normalization equation.
4. Certify the unique eigenpair in a specific interval box when the inclusion succeeds. If it does not, return inconclusive.
5. For a full real diagonalization, certify all eigenpairs, distinctness / independence, and invertibility of the eigenvector matrix. Repeated eigenvalues, clusters and invariant subspaces are later milestones.

For interval input matrices, distinguish certifying an eigenpair for one fixed matrix from certifying a consistent family of eigenpairs for **every** matrix represented by the input interval matrix. Document the quantifiers in the final API.

## Repository reconnaissance

- `include/algebra/matrix.hpp`: matrix operations, inverses and linear solves.
- `include/algebra/matrix.tpl.hpp`: orthogonal decomposition and other template implementations.
- `include/algebra/symmetric_matrix.hpp`: positive definiteness, validated Cholesky and SPD solves.
- `include/foundation/paradigm.hpp`: `ApproximateTag`, `ValidatedTag` and paradigm relationships.
- `include/numeric/float.decl.hpp`: `FloatApproximation<PR>` and `FloatBounds<PR>` for DP/MP.
- No general non-symmetric eigensolver was identified during initial code search.

## Proposed milestones

1. **Design and baseline:** this journal, API/contract specification, test matrix catalogue; preserve the baseline.
2. **Approximate real eigenpairs:** implement a clearly delimited first solver with unit tests and convergence/failure handling.
3. **Validated real simple eigenpairs:** implement rigorous inclusion and tests that prove enclosure, plus ambiguous/degenerate rejection or inconclusive status.
4. **Full diagonalization and complex spectra:** extend representation, repeated/clustered cases, numerical stability and certificates; decide Schur support.
5. **Integration readiness:** public API, documentation, examples, DP/MP coverage, CI confirmation and eventual Evolution consumption. No premature main merge.

## Test plan

- Scalar and diagonal matrices, nontrivial symmetric and nonsymmetric real matrices.
- Matrices with rational entries but irrational eigenvalues; non-normal and poorly conditioned examples.
- Positive/negative/zero eigenvalues, eigenvectors up to sign and scaling, permutations.
- Complex conjugate eigenvalues, defective Jordan blocks, repeated and clustered eigenvalues: fail or report unsupported correctly until implemented.
- Validate `M v = lambda v` and invertibility; use independent reference eigenvalues or exact constructions, never only self-consistency.
- For validated outputs, check mathematical enclosure semantics rather than merely a numerical tolerance.
- Test both double and multiple precision when supported, including failure paths.

## Change log

### 2026-11-27 — initial design commit

- Created the `eigen` branch from the current `main` baseline.
- Added this design and implementation journal.
- No eigensolver implementation or runtime tests yet.
- Next: settle result types and the first self-contained algorithm; add tests in the same implementation commit.

### 2026-11-27 — first restricted implementation (this commit)

- Added public header-only `EigenDecomposition<P,PR>` and a status-bearing `EigenDecompositionResult<P,PR>` with an optional successful data payload.
- The initial algorithm handles 1x1 and 2x2 real matrices with simple real spectrum. It computes both approximate and outward-rounded validated eigenvalues and normalized eigenvectors using the quadratic discriminant, and checks basis invertibility.
- Validated success has a **family-wise enclosure** contract: for every real matrix represented by the input intervals, its two ascending eigenvalues and normalized eigenvectors belong to the respective output intervals. `INCONCLUSIVE` is not a proof of absence of diagonalization.
- Added tests for real symmetric and nonsymmetric, triangular, diagonal and 1x1 cases, irrational eigenvalues, interval families, complex and repeated spectrum rejection, and DP/MP variants. Registered the new test and public header.
- All new C++ files use repository GPL headers with copyright 2026 Luca Geretti.
- **Verification status:** Code and tests committed for CI; a complete local build cannot run in the available environment (GitHub clone/dependencies unavailable). Passing tests are NOT claimed. Revisit compiler/runtime diagnostics before broadening scope.
- **Out of scope:** general dimension, complex spectra, clusters, Schur, QR, Krawczyk, and end-to-end PDE integration.

### 2026-11-27 — fix clang -Wshadow in the eigendecomposition (this commit)

- CI macOS/arm64 with Clang and `-Werror -Wshadow` rejected the local `const X two` in `EigenDecomposition<P,PR>::operator()`; `Ariadne::two` is already declared in `numeric/twoexp.hpp`.
- Renamed the local constant to `two_value` and updated both eigenvalue formulas. This is a naming-only correction; no change to numerical behavior or test coverage.
- **Verification:** inspected the reported diagnostic and the relevant source; checked the replacement locations. Full compilation and regression results for this new commit remain pending; the earlier CI build was stopped before the test binary could be produced.

### 2026-11-27 — stable 2x2 eigenvalue computation (this commit)

- **Baseline CI verified:** commit `16fdf00` passed Ubuntu, macOS, Windows, sanitizers and coverage jobs. `Format` failed; the external Codecov project/patch checks also report failure, independently of the passing coverage jobs.
- **Numerical change:** for a positive or negative provable trace, compute the dominant real eigenvalue using the cancellation-free branch of the quadratic formula, then recover the other from the characteristic determinant. When the trace sign is not known, retain the direct quadratic formula. A quotient is used only when the denominator is provably nonzero.
- **Regression coverage added:** real 2x2 matrices with widely separated eigenvalues, testing both positive and negative dominant eigenvalues; check the small eigenvalue explicitly.
- **Mathematical contract:** interval arithmetic still encloses the exact characteristic determinant and both eigenvalues for every represented matrix. This improves conditioning of the computed candidate but does not guarantee that every validated input will pass the final basis-invertibility certificate.
- **Verification at commit creation:** source edits inspected and new cases registered, but new commit has not yet been compiled or executed. Do not claim tests pass for this change.
- **Outstanding:** `check-format` uses clang-format 23 from the shared Build submodule. The current environment cannot execute that binary; resolve the format failure separately without changing the numerical algorithm. General-order eigendecomposition and fully validated eigenvectors remain future work.

### 2026-11-27 — general matrix determinant API (this commit)

- Introduced `determinant(const Matrix<X>&)` in `algebra/determinant.hpp` with a separate template implementation and `Copyright 2026 Luca Geretti` headers. Registered it in the umbrella header.
- The **approximate** implementation performs partial-row-pivot Gaussian elimination in O(n^3), tracking row permutation parity. This is a numerical approximation; no certified error bound is claimed.
- The **validated** implementation applies the Faddeev-LeVerrier characteristic-polynomial recurrence in O(n^4), avoiding uncertain pivots and division by input-dependent intervals. All operations are outward-rounded, thus the result encloses the determinant of **every** represented matrix, including singular matrices. Overestimation is expected for wide interval inputs.
- The empty matrix convention is `det([])=1`. Nonsquare input violates the precondition.
- Added the standalone `test_determinant` target, testing 0x0, 1x1, 3x3, 4x4 row permutation, singularity, a 2x2 interval matrix family and 6x6 diagonal matrix, in DP and MP.
- **Verification:** Source and test definitions reviewed, but this commit has not yet run CI or a local build. No pass is claimed. The preceding `e4fa739` commit passed Ubuntu, macOS, Sanitizers and Coverage CI; Windows was still running at inspection.
- Per user request, `check-format` is explicitly deferred to pre-merge work; no formatting intervention is required during implementation.
- **Open:** stronger numerical pivots and scaling for ill-conditioned approximate determinants, tight validated determinant enclosures, and general-order eigenvector computation.

### 2026-11-27 — triangular matrices of arbitrary dimension (this commit)

- Extended `EigenDecomposition<P,PR>` to square real **upper/lower triangular** matrices of order >= 3 when the diagonal eigenvalue intervals admit a strict common increasing order (or numeric strict order for ApproximateTag). The current 1x1/2x2 dense algorithm is unchanged.
- Computed normalized eigenvectors using backward substitution for upper triangular inputs and forward substitution for lower triangular inputs. Each column uses the corresponding diagonal position as an exact unit pivot; entries beyond the triangular support are exact zero.
- **Validated contract:** all involved denominators exclude zero for every represented matrix, each recurrence uses outward-rounded operations and hence encloses the normalized eigenvector family. With columns ordered by disjoint diagonal eigenvalue intervals, the eigenvector matrix is a column permutation of a unit triangular matrix; invertibility is structural and does not depend on an overestimated interval determinant.
- Diagonal intervals that cannot be ordered consistently return `INCONCLUSIVE`. Dense higher-order matrices continue to return `UNSUPPORTED_DIMENSION`; complex or repeated eigenvalues remain outside the milestone.
- Added DP/MP tests for upper and lower triangular 4x4 matrices, matrix-eigenpair residual inclusion, strictly separated interval diagonal families, overlapping intervals and unsupported dense 3x3 inputs.
- **Verification:** code and tests published for CI, with runtime results still pending at commit time. `check-format` is intentionally deferred until pre-merge, as requested.

### 2026-11-27 — guard zero-order eigendecomposition (this commit)

- During review of the triangular extension, identified that the prior early dimension guard excluded order zero, while the new branch only rejected dense order >2. Added an explicit order-zero `UNSUPPORTED_DIMENSION` return before any matrix indexing and a regression test.
- No change to the supported nonempty matrices. `determinant` intentionally treats the empty square matrix differently, returning one by convention.
- **Verification:** source inspection and regression added; runtime CI remains to be confirmed. `check-format` intentionally deferred until pre-merge.

### 2026-11-27 — dense symmetric approximate eigendecomposition (this commit)

- Extended the dimension-general solver to **real symmetric dense matrices** (order >=3) for `ApproximateTag`, using a maximum-off-diagonal-pivot Jacobi rotation algorithm.
- Require exactly matching stored symmetric matrix entries; arbitrary dense nonsymmetric matrices are still unsupported. Do not silently symmetrize user inputs.
- Compute plane rotations with a stable tangent formula; diagonalize iteratively, sort real eigenvalues, reject unresolved repeated/very close eigenvalues, and return `INCONCLUSIVE` for failed iteration or residual checks.
- Jacobi convergence, pivot selection and residual acceptance are currently assessed with double-precision estimates even if the underlying numeric type is MP. Thus **MP is supported as an input/output type, but this implementation makes no claim of exploiting full MP accuracy**. A precision-native convergence policy is necessary for a mature generic eigensolver.
- ValidatedTag deliberately remains restricted to existing 1x1, 2x2 and triangular workflows: an approximate Jacobi eigenpair must **not** be returned as a validated result.
- Added DP and MP regression tests for a symmetric tridiagonal 3x3 with exact known spectrum `2,2±sqrt(2)`, a 4x4 dense symmetric construction with known spectrum `1,2,3,4`, and rejection of nonsymmetric/uncertified cases.
- **CI status at preparation:** earlier commit `0bd47a6` has passed Ubuntu/macOS/Sanitizers/Coverage; newer branch CI still ongoing. Tests for this new implementation have not yet been executed. `check-format` deferred per user instruction.
- **Next:** address CI findings, test near-degenerate and ill-conditioned symmetric systems, then introduce a robust nonsymmetric approximate backend and verified simple-eigenpair certification.

### 2026-11-27 — native-precision Jacobi convergence and residuals (this commit)

- Confirmed `ac3598e` passes Ubuntu, macOS, Sanitizers and Coverage workflows; Windows remained in progress when checked. Per user instruction, ignore `check-format` until the merge.
- Removed binary64 `get_d()` conversions from the Jacobi algorithm: finite checks, matrix scale, pivot choice, convergence threshold, rotation sign, eigenvalue sorting, gap rejection and residual checks now use the native raw DP/MP representation.
- Derived thresholds from `RawFloatType<PR>::eps(precision)`, so an MP run can resolve spectral couplings far below binary64 roundoff. These remain *heuristic* numerical acceptance thresholds; no validated guarantee is implied.
- Updated the shared `ApproximateTag` sign/zero tests to use raw arithmetic rather than `get_d()`; extremely small or large finite MP values can no longer be mistaken for zero or infinity solely during a double conversion.
- Added an MP160 regression with a dense symmetric matrix containing an off-diagonal coupling `2^-75`, demonstrating detection of a nontrivial eigenvector and checking residuals against an MP-dependent budget. Existing DP and MP tests remain.
- **Verification:** inspected edits and test coverage, but no CI execution for this new commit is claimed. Further stress tests (tight clusters, large exponents, finite overflow) are warranted before using the algorithm in production.

### 2026-11-27 — validated dense eigenpairs through separated Gershgorin disks (this commit)

- Added `EigenDecomposition<ValidatedTag,PR>` support for **dense real interval matrices of arbitrary dimension** under the sufficient condition that their outward-rounded Gershgorin disk real projections are pairwise strictly separated. This includes nonsymmetric matrices; the criterion itself proves every eigenvalue is simple and real for each represented matrix.
- For every disk i, enclose its unique real eigenvalue by the interval projection `A_ii + [-r_i,r_i]`, with `r_i = sum_{j != i} |A_ij|`. Gershgorin's second theorem proves the eigenvalue count, not a small floating-point residual.
- Normalize the associated eigenvector to have component i equal to one and solve the remaining `(lambda I - A_rest)v_rest = A_rest,i` equations with outward-rounded interval Gaussian elimination. Abort inconclusively if any interval pivot contains zero; pointwise elimination is otherwise inclusion-preserving for all members and eigenvalues in the disk.
- Verify the interval eigenvector columns as a matrix `V=I+E` **before** sorting: strict `||E||_1<1` proves the eigenvector basis invertible for every matrix; sorting only permutes the columns.
- Return `INCONCLUSIVE` when the disks are not separated, an interval linear solve cannot be established, or basis invertibility is not certified; do not claim a negative mathematical result.
- Registered DP/MP regression cases for a symmetric 3x3 interval family, a dense nonsymmetric 4x4 system, and overlapping Gershgorin disks. Updated formerly unsupported dense-validated tests to the more accurate `INCONCLUSIVE` status.
- **Status:** Source inspection completed; compilation and test runs for this new commit remain pending. The preceding native-MP Jacobi commit `760f2d9` was still running the CI at time of design. `check-format` is deferred until pre-merge as requested.
- **Limitations:** The criterion is deliberately conservative: a symmetric dense matrix with clustered spectra may have an excellent simple real eigendecomposition but fail the disk-separation test. A general approximate-seeded interval-Newton/Krawczyk certification remains future work.

### 2026-11-27 — verified similarity preconditioning for dense validated eigenpairs (this commit)

- **Baseline verified:** `306ed2e` passed Ubuntu, macOS, Windows, Coverage and Sanitizers. `check-format` remains deferred until pre-merge; the Codecov project/patch status checks are separate.
- Introduced `_preconditioned` for dense `ValidatedTag` matrices when direct interval Gershgorin disk isolation fails.
- Construct a symmetric **approximate** seed matrix from input interval midpoints and run the ApproximateTag Jacobi algorithm to obtain a candidate basis S. It is only a choice of *fixed coordinates*, not a certificate of the input matrix (which may be nonsymmetric).
- Regard S as a point matrix with exact dyadic coefficients. Apply verified interval Gauss-Jordan to enclose S^{-1}; return `INCONCLUSIVE` when any required pivot includes zero. In particular, do not replace inverse by an unchecked approximate transpose of the Jacobi output.
- Form the interval enclosure of the exact similarity B=S^{-1}AS. Apply the previously certified Gershgorin / eigenvector isolation to B, requiring real isolated eigenvalues, verified eigenvector interval solves and invertibility of the transformed eigenvector basis W.
- Return the original eigenvalues and transformed eigenvectors V=SW. The basis V is invertible because both S and W are; unlike the direct Gershgorin path, its columns need not have a component equal to one.
- This is a conservative but rigorous sufficient method, not a general verified QR algorithm. It may succeed even for some nonsymmetric dense matrices with overlapping original Gershgorin disks.
- Added DP/MP tests for the symmetric tridiagonal matrix with known spectrum 2, 2±sqrt(2), a nonsymmetric perturbation and an uncertain interval family that must remain inconclusive.
- **Verification:** repository source and tests reviewed and committed; this new commit has not yet passed a CI build. Do not interpret any unverified test as already passing.
- Remaining: evaluate interval overestimation from Gauss-Jordan and transformed matrix products; design a general interval-Newton/Krawczyk refinement for difficult isolated eigenpairs.

### 2026-11-27 — eigenvalue/eigenvector performance benchmark (this commit)

- Added optional auxiliary `ariadne-cps/benchmarking` submodule `ss/bmk` at `8dca1e7eb605a378852e574632c9e9b82641f582`, mirroring Numeric. Verified shared Build `5e6444f` and Utility `341736b` pins against Interval/Numeric/Foundation and Benchmarking/Logging.
- Registered a separate `benchmarks` target via Build's existing `ariadne_setup_benchmarking` root orchestration, and linked only benchmark executable to `ariadne-benchmarking` plus any additional Logging objects (following Numeric's standalone linking pattern). Production Algebra remains independent of Benchmarking.
- Introduced `benchmarks/benchmark_eigen_decomposition.cpp` with deterministic test matrices and four distinct paths: dense 2x2, upper triangular, well-separated dense, and overlapping-disk rotated dense. Both Approximate and Validated paradigms; DP53 and configurable MP bits.
- Uses shared `Runner`, `Statistics`, `Report`, CLI parser and optimization barrier. Inputs and pilot status are computed before timing; each measured operation is the entire eigenvalue/eigenvector decomposition. Check statuses after every round, and mark inconclusive results explicitly.
- Stable CSV schema: `backend,workload,n,bits,status,median_ns,mad_ns,min_ns,max_ns,p95_ns`. Compare only rows with the same semantic outcome; no false speedups by treating an inconclusive attempt as successful. Configurable sizes, bit precisions, repetitions, samples, warmups and summary columns.
- Added Release benchmark build/run examples to README. No performance measurements are claimed for this commit: benchmark build and execution remain to be verified in CI or on a dedicated machine.
- `check-format` remains deferred to pre-merge as requested. Next measure actual timings, identify per-path bottlenecks and add A/B implementations without changing the fixed workloads.
