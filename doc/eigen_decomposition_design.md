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
