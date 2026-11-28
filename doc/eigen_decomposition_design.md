# Eigenvalue and eigenvector computation — design journal

## Scope and workflow

Development branch: `eigen` in `ariadne-cps/algebra`, branched from `main` at `1b85f7e864044b697be2408c3575181a75f36d4b`.
Do not merge into `main` until implementation, numerical and semantic tests, and review are complete.

**Journal rule:** Every development commit must include an update to this document describing its design decisions, code changes, tests run, test results and outstanding issues. The first commit creates this journal. No successful build or test is to be claimed without running it.

## Purpose

Implement reusable approximate and validated spectral methods in Algebra. Evolution's first-order PDE solver currently takes user-supplied characteristic eigenvalues/eigenvectors. Algebra must remain independent of PDE assumptions; integration with Evolution is a later project.

## Current target architecture — approved 2026-11-27

**This section supersedes the historical "Public API direction", "Mathematical
contracts", and "Proposed milestones" below.** There are no downstream users
of the experimental spectral API. We are free to change, replace, or delete the
current classes and result types without compatibility adapters or deprecation.

The intended Algebra responsibilities are independent and layered:

1. **Spectrum<P,PR>** computes all eigenvalues of a real square matrix, with
   algebraic multiplicities, in a **complex-valued representation**. The result
   never implies existence of an eigenvector basis. Order of entries is not
   part of the public contract: Schur block order or structural coordinate
   order is permitted, and interval eigenvalue boxes may overlap. If the
   input is an interval family, SUCCESS means that for *every* represented
   real matrix, its eigenvalue multiset has a labeling enclosed by the
   returned complex rectangles, with multiplicities. ApproximateTag is
   not a certificate; ValidatedTag must use interval inclusion proofs.
2. **SchurDecomposition<P,PR>** (future) represents A=Q*T*Q^T for real
   matrices, with orthogonal Q and upper real quasi-triangular T (1x1
   real blocks and 2x2 conjugate-pair blocks). It is the general numerical
   backend for Spectrum and works for matrices that are defective or
   have nonreal eigenvalues. Initially implement an approximate real
   Hessenberg reduction and implicit double-shift Francis QR/deflation.
   A validated Schur contract needs separate rigorous analysis; do not
   pretend a floating Q and T prove an interval identity.
3. **EigenDecomposition<P,PR>** is requested only when a full eigenvector
   basis (over an explicitly specified real or complex scalar field) is
   desired. SUCCESS must guarantee the required basis/invertibility in
   the requested paradigm. A failed diagonalization is not a failed
   spectrum. Existing real eigensolvers/interval certificates are valuable
   implementation primitives but the current public shapes are temporary.
   Do not preserve them merely for backwards compatibility.

**One public eigendecomposition API and the real fast path (decision 2026-11-28).**
There must be exactly one user-facing `EigenDecomposition<P,PR>` for a
real input matrix: the caller must not have to know in advance whether
the eigenvalues are real or complex. The interim
`ComplexEigenDecomposition` implementation is not a second permanent
public solver. A complete real basis is a successful special case of
the unified operation. The intended success payload is a tagged
real/complex representation (e.g. `std::variant<RealEigenData,
ComplexEigenData>`), **not** eager promotion of every real eigenbasis to
complex coefficients.

- Preserve specialized real paths (symmetric Jacobi; triangular/block
  structure; existing validated Gershgorin/Krawczyk/similarity pathways)
  where they are applicable. The generic real Schur backend runs only
  real-arithmetic Householder/Hessenberg/Francis/Wilkinson iterations.
- When Schur yields only real scalar blocks, construct/refine/check
  eigenvectors and their inverse in the **real** number system. Avoid
  allocating dense complex matrices or performing complex arithmetic
  on this overwhelmingly common path.
- Only use complex back-substitution, complex vectors and complex
  inversion when nonreal conjugate 2x2 Schur blocks require them.
  Defective/unresolved clusters remain INCONCLUSIVE; an approximate
  factorization must never be passed off as a validated certificate.
- Measure the performance impact before declaring completion. In
  DP and MP, compare wall time, allocations, and iteration counts
  against the last proven real implementation for symmetric,
  triangular, non-normal but real-spectrum, and mixed-spectrum
  examples. Do not assume the new Schur path is cost-free: the
  change in real QR algorithms can itself affect speed. Execute
  benchmarks locally only; CI may compile but must not run them.
- Remove/replace transitional public names and storage types during
  final unification: there are no downstream compatibility clients.

**Status semantics:** SUCCESS gives the complete requested result (with
the requested proof strength), INCONCLUSIVE means not established and never
asserts nonexistence, and UNSUPPORTED_DIMENSION handles order zero. Future
handling of known non-diagonalizability must be distinguished from failure to
certify it. For validated interval families never infer correlation between
independently represented entries. Repeated enclosures and crossing
eigenvalues require no artificial fixed sorting of eigenvalue boxes.

**Implementation roadmap, in order:**
- Stage A: stabilize CI baseline and implement independent Spectrum result
  for 1x1/2x2 real input with complex pairs and repeated defective spectra,
  plus triangular matrices of arbitrary order. Exercise DP/MP approximate
  and validated families. Keep dense order >=3 explicitly inconclusive.
- Stage B: introduce general approximate real Schur (Householder Hessenberg,
  implicit double-shift QR, deflation and nonconvergence reporting), with
  tests for nonreal pairs, Jordan cases, clusters and non-normal matrices.
  Use Schur to supply Spectrum without requesting a diagonalization.
- Stage C: extract complex eigenvectors when possible; define exact result
  types for complete diagonalization and invariant subspaces, retaining
  the already proved real validated pathways as internal building blocks.
- Stage D: verified complex isolated eigenvalues and spectral clusters,
  including multiplicities/invariant subspace enclosures where possible.
  Failure to certify remains INCONCLUSIVE, never a false proof.
- Stage E: hard numerical/semantic regressions, DP/MP precision studies,
  benchmarks performed locally only (CI compile-only), coverage and format
  review, and only then consider merging into main.

**CI baseline investigation:** commit 919e3da passed functional macOS,
Windows, Coverage and Sanitizers. Its sole Ubuntu x86 GCC 13 Release
failure occurred during "Set up Linux dependencies" before CMake, any
build or any test. This is not evidence of an RTL/C++ regression; the
exact package/network error is not exposed in the available job metadata.
No code or workflow change is justified solely by that result.
The check-format and external Codecov policy failures remain separate,
with format deferred until pre-merge as previously agreed.

## Historical public API direction (superseded by current target architecture)

- `EigenDecomposition<P, PR>`: algorithm/solver object, with `P` in `ApproximateTag` or `ValidatedTag` and `PR` denoting numeric precision.
- `EigenDecompositionResult<P, PR>`: corresponding result type. `ApproximateTag` returns approximations without certification. `ValidatedTag` provides mathematically guaranteed enclosures **only when certified**, and an explicit inconclusive status otherwise.
- Support real and complex eigenvalues eventually; do not suggest every real matrix has a real eigenbasis.
- A failed certificate must never be reported as proof that no decomposition exists.
- Prefer an optional validated result over a success-looking empty or invalid enclosure.
- Input matrices should use the existing Ariadne number/matrix conventions. The exact call signatures, storage of complex outputs, convergence controls and precision traits remain to be decided after checking available primitives.

## Historical mathematical contracts (superseded by current target architecture)

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

## Historical proposed milestones (superseded)

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

### 2026-11-27 — restore eigendecomposition regression fidelity after similarity preconditioning (this commit)

- Diagnosed Debug CI failure from `07ea706` and `b8a197d`: the pre-existing 3x3 example `[[0,1,0],[1,0,0],[0,0,0]]` was still asserted to return `INCONCLUSIVE`. Once the validated similarity path was introduced, it can certify this real simple spectrum. This is an **outdated test expectation**, not an indication that certification must be weakened.
- Preserved the original example, now asserting `SUCCESS`, enclosures for its independently known real eigenvalues `-1, 0, +1`, and interval overlap of `A*V` with `V*diag(lambda)` for every component. A failure in the new certificate will now be reported as an algorithm regression, not hidden by accepting either status.
- Added a distinct 3x3 matrix with repeated spectrum `-1, 1, 1` to keep explicit `INCONCLUSIVE` coverage for an unsupported simple-spectrum contract.
- The CI log's `last check [1.230]: residual.raw() <= tolerance.raw()` is test context; the failing assertion is `no_dimension.status == Status::INCONCLUSIVE`. Compared the MP residual test between `306ed2e` (where CI passed) and current branch: its body is unchanged. Do **not** relax the MP residual tolerance or alter its scenario.
- Reviewed `ariadne-cps/equation-solving` Krawczyk and interval-Newton implementations. Their public solvers depend on Algebra, Function and function-patch types. Keep `algebra` free of dependency cycles. A future shared numerical Krawczyk kernel may be placed in Algebra as a matrix/vector interval primitive, with the higher-level equation-solving solver wrapping it; do not import its public API into Algebra.
- **Verification:** exact source edit and test assertions reviewed; new CI on this commit remains pending. The previous Release jobs compiled; Debug tests failed on the stale expectation. `check-format` remains intentionally deferred until pre-merge.

### 2026-11-27 — minimal shared Newton/Krawczyk numerical kernel (this commit)

- Reviewed the existing `equation-solving` implementations of `IntervalNewtonSolver::step` and `KrawczykSolver::step`. Their shared numerical operations depend only on a center, residual vector, interval Jacobian and (for Krawczyk) a caller-chosen preconditioner. Their higher-level function evaluation, Jacobian generation, iteration, patch handling and root certification remain outside this layer.
- Introduced exactly **two stateless functions** in `algebra/validated_root_contractors.hpp`: `interval_newton_image(center, residual_at_center, jacobian)` and `krawczyk_image(center, residual_at_center, jacobian, box, preconditioner)`. They use existing `Matrix`, `Vector`, `inverse` and arithmetic operations; there is no new dependency, solver hierarchy, result type, or iteration policy.
- Both functions compute **interval images**, not proofs of existence/uniqueness. The calling solver must check the applicable hypotheses, enclosure, contraction and strict inclusion; an invalid inverse may fail. The Krawczyk function intentionally takes the preconditioner as an input rather than selecting an approximate inverse internally.
- Added standalone DP/MP tests of a coupled 2x2 linear system (known exact root), a nonlinear scalar function with independently known root `sqrt(2)`, and a zero preconditioner to demonstrate that image computation alone cannot certify a root.
- Registered the header in the Algebra umbrella and the test in the existing test target list. No changes to `equation-solving` and no dependency cycle.
- **Verification status:** source and interfaces inspected and changes committed for CI; local C++ compilation and test execution not performed. As of inspection, the preceding regression fix `ace1e92` had passed macOS, Ubuntu, Sanitizers and Coverage; Windows remained in progress. `check-format` remains deferred until merge.
- **Next step:** once the DP/MP tests pass, have `equation-solving` call these two helpers from its existing solver steps (in that repository and its own tests) before considering any richer shared abstraction. Implement eigenpair-specific validated Newton/Krawczyk refinement only once the shared operation is exercised in both consumers.

### 2026-11-27 — route eigendecomposition benchmark output through Logging (this commit)

- Reviewed Numeric benchmarks and shared Benchmarking `Report`: Numeric uses `LOGGING_PRINTLN` for informational and error output, and `Report::print()` to present tables. The Algebra eigendecomposition benchmark previously wrote directly to `std::cout`/`std::cerr`, bypassing the Logger.
- Replaced the table banner and configuration summary with `LOGGING_PRINTLN`, then invoked `Bench::Report::print()` for the tabular report. Replaced error output in `main` with `LOGGING_PRINTLN`. Removed the no-longer-needed `<iostream>` include.
- Kept the established CSV column names and order, constructing each CSV record with `std::ostringstream` at 17-digit precision and printing each complete record through `LOGGING_PRINTLN`. This keeps CSV data together per line and avoids introducing extra diagnostics in CSV mode; CSV output now follows the Logger's configured destination and verbosity.
- No benchmark timing, workload construction, numeric semantics, status checks, or benchmarking linkage changed. All presentation remains outside measured callbacks.
- **Verification:** source references and exact replacements inspected; compilation/CLI output for this commit remain to be checked. `check-format` remains deferred until pre-merge.

### 2026-11-27 — compile and exercise eigendecomposition benchmarks in Ubuntu Release CI (this commit)

- Audited the Algebra and Numeric Ubuntu workflows. The earlier Algebra release jobs built the library and bindings but **never compiled the new benchmark target**; in contrast, Numeric explicitly built `benchmarks`. Therefore a green Algebra CI was not evidence that `benchmark_eigen_decomposition` compiled or produced output.
- Mirrored Numeric's Release-only `cmake --build ... --target benchmarks` step in Algebra's Ubuntu workflow (GCC/Clang, x86/arm).
- Added a short Release execution exercising both formatted reporting and CSV with `--sizes 3 --bits 80 --iterations 1 --samples 3 --warmups 0`. This exercises DP and MP, ApproximateTag and ValidatedTag, and each matrix-workload family, while keeping timings outside correctness tests.
- Documented the smoke check in README. CI timing results are **not performance benchmarks**: accurate A/B measurements require controlled Release runs on the same stable machine.
- **Verification:** inspected workflow syntax and CMake/build target contracts; GitHub Actions results for this commit still pending. `check-format` remains deferred until the merge; no formatting or workflow changes related to Format.

### 2026-11-27 — enforce compile-only benchmark CI (this commit)

- User requirement: **CI must compile benchmarks, never execute benchmark binaries**. The preceding `a176b81` introduced a Release smoke execution step; this was explicitly rejected as contrary to the benchmarking policy.
- Removed the `Smoke-test Release eigen benchmark` step in `.github/workflows/ubuntu.yml`. Kept the `Build Release eigen benchmarks` step and its `cmake --build ... --target benchmarks` command in the Release matrix, following Numeric.
- Updated README to state that benchmark runs are manual and outside CI. No changes to benchmark algorithms, compiler configuration, functional tests or formatting workflows.
- **Verification:** inspected the resulting workflow and README source; new CI results pending. The preceding commit may already have executed its smoke step before this corrective commit, but future runs of the updated workflow contain no benchmark invocation. `check-format` stays deferred.

### 2026-11-27 — reusable Krawczyk root-certification predicate (this commit)

- Baseline `09e7057`: Ubuntu Release successfully compiled the optional benchmark target **without running benchmarks**; Ubuntu, macOS, sanitizers and coverage succeeded, while Windows Debug and Format remained pending/deferred at inspection.
- Kept the shared numerical API minimal and header-only. Added `krawczyk_certifies_unique_root(center,residual_at_center,jacobian,box,preconditioner)`, alongside the existing Newton and Krawczyk image operations; no solver classes, function abstractions, new repositories or additional dependencies.
- Mathematical contract: for a fixed continuously differentiable real `f`, the caller must provide a **point** center, an enclosure of `f(center)`, and an interval Jacobian enclosing `Df` for every point of the candidate box (and every uncertain parameter). The method checks strict containment of the Krawczyk image inside the box and an outward-rounded bound `||I-C*J||_infinity < 1`. These jointly give a self-map and a contraction, hence a unique zero in the box for each fixed function, and rule out singular preconditioners. Mere image inclusion without contraction is not accepted as a uniqueness proof.
- The certificate returns false (inconclusive) for nonpoint center, nonfinite or noninterior data, failed contraction, or failed strict containment. An image alone is still not a certificate; the preconditions about user-supplied residual/Jacobian are essential and are not inferred.
- Extended the DP and MP tests: coupled 2x2 linear system and `x^2-2` must pass the Krawczyk unique-root check; a zero preconditioner, a wide Jacobian containing a critical point, and a nonpoint center must fail. Existing image-enclosure tests remain unchanged.
- This algebraic predicate is immediately reusable by an eigenpair-specific Krawczyk contractor; defer updating `equation-solving` itself until a compatible Algebra revision is available through its pinned Function/Algebra dependency chain, avoiding broken builds or a dependency cycle.
- **Verification:** source and contracts reviewed and committed, new CI to be checked. `check-format` remains deferred until pre-merge. Benchmarks are compiled but **never executed** in CI.

### 2026-11-27 — compile benchmarks on every Release CI platform (this commit)

- Corrected the previous Ubuntu-only scope after checking the actual `ariadne-cps/numeric` workflows: Numeric compiles `benchmarks` in **Ubuntu, macOS and Windows Release** (including MSVC x64/ARM64), always as a build-only target; no CI benchmark execution occurs.
- Added to Algebra `.github/workflows/macos.yml` and `.github/workflows/win.yml` the same Release-only `Build Release benchmarks` step that Numeric uses, immediately after the main library/bindings compilation. macOS uses `cmake --build "${{ runner.workspace }}/build" --parallel "$NPROC" --target benchmarks`; Windows uses `cmake --build build --parallel --target benchmarks`.
- Left the already verified Ubuntu Release build step unchanged, and did not add execution, smoke tests, benchmarking measurements or changes to any Debug, Format, Sanitizers or Coverage workflow.
- Updated README to describe the three-platform, compile-only policy. This is a CI coverage change only: functional code and workload definitions are unchanged.
- **Verification:** checked workflow content against Numeric and confirmed required steps and Release conditions in the committed source. Full macOS/Windows compile results for this commit remain pending. Ignore `check-format` until pre-merge as requested.

### 2026-11-27 — Krawczyk fallback for real eigenpairs of non-normal dense matrices (this commit)

- Added a third validated dense fallback `_krawczyk` after direct Gershgorin and verified-similarity Gershgorin fail. It uses `krawczyk_certifies_unique_root` from the small, existing Algebra root kernel; no dependency on `equation-solving`, no new solver hierarchy.
- For every Jacobi seed of the **symmetrized midpoint**, normalize the candidate eigenvector at its largest-magnitude component. Solve the square real system `F_i(lambda,v)=(Av)_i-lambda*v_i=0` with that component fixed to one; it has exactly n unknowns and n equations.
- Approximate Newton refinement operates on the actual, possibly nonsymmetric midpoint matrix to construct *starting guesses only*. Its results do not certify anything. A point preconditioner is computed from the approximate analytic Jacobian and reinterpreted as a matrix of exact floating-point constants in interval arithmetic.
- Build a verified interval residual at the exact center and an interval Jacobian enclosing derivatives for every input matrix and every unknown in a candidate box: `dF_i/dlambda=-v_i`, `dF_i/dv_j=A_ij-lambda*delta_ij`. Choose tentative radii from preconditioned residuals and precision-dependent slack, then attempt the strict Krawczyk contraction/inclusion test. Failure is `INCONCLUSIVE`, never evidence of a missing root.
- Each successful Krawczyk check proves an actual real eigenpair for **every** represented matrix. Require all n eigenvalue intervals pairwise strictly separated; this ensures the roots are simple, consistently sortable, exhaust the spectrum, and yield an invertible real eigenvector basis for each member. Sorting occurs **after** certification because the symmetric seed order need not agree with the nonsymmetric actual spectrum.
- Added DP and MP regression tests of a 4x4 non-normal matrix obtained by a known **integer unimodular similarity** of `diag(1,3,5,7)`. Original and symmetrically preconditioned Gershgorin disks overlap. Tests require a successful certificate, the independently known eigenvalue enclosures, separated ordering, interval eigenpair residual overlap, and a small nonzero interval input perturbation.
- **Verification:** implementation and proof obligations reviewed; CI tests on this commit remain pending. Prior `d606a77` macOS/Windows Release benchmark builds were still in progress at the time of this edit. No benchmark binary is executed in CI. `check-format` remains deferred.
- **Limitations:** candidate discovery still depends on Jacobi for the symmetric seed; the numerical Newton iteration may not reach the correct non-normal eigenpair or may converge multiple seeds to the same root, and rigorous Krawczyk may fail due dependency/rounding. Such cases return `INCONCLUSIVE`. Wider intervals and clustered or complex spectra remain open.

### 2026-11-27 — benchmark non-normal validated eigendecomposition (this commit)

- Verified `295ba45` Krawczyk tests have passed Ubuntu, macOS, Coverage and Sanitizers. At review, Windows ARM64 jobs were still running. `check-format` is intentionally deferred.
- Existing benchmark catalogue had no nonsymmetric dense matrices, so it could not measure the new Krawczyk eigenpair fallback. Added two **Validated-only** fixed-order 4x4 workloads: `dense_nonnormal_4x4` (point coefficients) and `dense_nonnormal_interval_4x4` (a small interval perturbation).
- Reused the exact input coefficients from the Krawczyk regression (known eigenvalues 1,3,5,7) and perturbation of `A[0][0]` by `[0,1]/1048576`. Preserved existing workload names and generated inputs for baseline comparisons. Both new workloads require a successful validated certificate at pilot measurement; inconclusive attempts abort instead of masquerading as improved performance.
- Do not add the non-normal cases to `ApproximateTag`: the current approximate dense solver supports symmetric inputs only, so `UNSUPPORTED_DIMENSION` would make the benchmark incomparable. Added two configurations per validated backend and updated Runner's expected configuration count accordingly.
- These cases measure **complete calls**, including direct Gershgorin, preconditioned Gershgorin and fallback attempts. Their execution does not by itself prove which internal strategy succeeded; future path-specific profiling should use explicit instrumentation rather than inference from matrix structure.
- Documented fixed 4x4 workloads and result interpretation in README. No changed CI policy: benchmark binaries are compiled on Ubuntu/macOS/Windows Release but **never executed in CI**.
- **Verification:** source changes and dimension/count consistency inspected. CI compilation of the new benchmark commit is pending; no timing numbers claimed. Next run controlled external Release A/B measurements and consider exposing algorithm strategy metrics without contaminating timed sections.

### 2026-11-27 — optional per-stage eigendecomposition profiling (this commit)

- Baseline `802edc2` has passed Ubuntu, macOS, Windows, Sanitizers and Coverage. Formatting is intentionally deferred until pre-merge, and benchmark executables remain compile-only in every Release CI matrix.
- Added an **optional** `EigenDecompositionDiagnostics` and strategy enum to the Algebra eigendecomposition API, with a two-argument overload `solver(A,&diagnostics)`. The existing `solver(A)` overload delegates with a null pointer; no instrumentation is active and no clocks are read in normal mode.
- Diagnostics identify the *actual* successful strategy (or last attempted one for inconclusive results): analytic, triangular, Jacobi, direct Gershgorin, verified similarity, or eigenpair Krawczyk. Added scoped `steady_clock` counters for disjoint approximate candidate-generation/Newton and verified interval-certification stages. The counters exclude some dispatch/bookkeeping work and must not be interpreted as additive exact partitions of the complete runtime.
- Added `--profile` to the existing performance executable, preserving existing default benchmark/CSV output. In profiling mode, `Runner` collects per-round per-call candidate and verification nanoseconds, excluding warm-ups with the same policy used for total timings; it reports the median by strategy in both the table and optional CSV columns. Timer overhead is present **only** in profile mode, so compare like with like.
- Added DP/MP functional regression checks requiring identical uninstrumented/instrumented statuses and eigenvalue representations for direct validated Gershgorin and approximate Jacobi. Profiling does not alter eigenpair calculations or interval acceptance criteria.
- README documents profile mode and an optional manual Release CSV command. The CI continues to **compile** benchmark binaries on Ubuntu/macOS/Windows Release and never runs them.
- **Verification:** source inspections and coherent API/callback pairing completed. New commit CI is pending; no elapsed timings or success for new instrumentation claimed. Future use should first compare per-stage medians across controlled external builds, then target verified inverse construction and interval dependency inflation.

### 2026-11-27 — first dense nonsymmetric ApproximateTag eigensolver (this commit)

- Established baseline: previous `c176683` passed Ubuntu, macOS, Windows, Sanitizers and Coverage. Format remains deferred per the user's instruction.
- Added an approximate eigenpair backend for **dense nonsymmetric real square matrices of order >=3** with a fully resolved simple real spectrum. This is a basic functional milestone, not a general complex eigensolver and not a validated result.
- Implemented shifted real Householder QR iterations with full trailing-row deflation, accumulating approximate orthogonal similarities. When a trailing real 2x2 block has non-real eigenvalues, or QR does not converge within its finite iteration budget, return `INCONCLUSIVE`. Following deflation, each QR similarity updates the active-to-deflated cross block as well as the active square block, preserving consistency for the Schur eigenvectors.
- Back-substitute approximate upper-Schur eigenvectors, map them to original coordinates, then run Newton refinement on the actual matrix using the square normalized system `(A-lambda I)v=0` with a selected component of v fixed to one.
- Require finite inputs and intermediate values, numerically separated sorted eigenvalues, small original-input residuals, and numerical invertibility of the eigenvector basis. These are numerical acceptance tests only; no interval enclosure is claimed. Unsupported/unconverged/complex/repeated cases return `INCONCLUSIVE`.
- Symmetric ApproximateTag matrices continue to use Jacobi. Triangular matrices continue to use the dedicated recurrence. The new backend is selected only when both structural fast paths do not apply and Jacobi rejects asymmetry.
- Registered the strategy as `SHIFTED_QR` for the existing opt-in diagnostics and human-readable benchmark reporting, without further profiling work.
- Added DP/MP tests for two 3x3/4x4 nonsymmetric real matrices with independently known spectra, residuals and strategy dispatch; tested the 3x3 previously unsupported matrix and separate complex-spectrum/repeated-eigenvalue rejection scenarios. Updated the benchmark to measure the same non-normal 4x4 point workload in both ApproximateTag and ValidatedTag while keeping the interval workload validated-only.
- **Verification:** Householder QR and real-shift deflation behavior were checked in a separate numeric prototype on the chosen reference matrices; source/test changes are committed for C++ CI. No local C++ build result or new CI pass is claimed yet. Benchmark executables remain **compiled only**, never run by CI.
- **Limitations:** shifted dense QR without Hessenberg reduction and full real 2x2 Schur handling is computationally expensive and may not converge for strongly nonnormal or clustered cases. This deliberately restricted simple-real backend will later need more robust QR/Schur and optional complex eigenspectrum support; optimize only after baseline features stabilize.

### 2026-11-27 — repeated real eigenvalues with complete eigenvector bases (this commit)

- Confirmed preceding `dbee294` passes Ubuntu, macOS, Windows, Sanitizers and Coverage; Format remains deferred until pre-merge. All three Release CI platforms continue to **compile benchmarks without running them**.
- Corrected an overrestrictive semantic assumption: **repeated eigenvalues do not imply a defective matrix**. For exactly diagonal matrices of any positive order (including 2x2), introduced a coordinate-eigenvector fast path that accepts repeated real eigenvalues and returns a sorted basis permutation. This bypasses the 2x2 discriminant test, which previously rejected scalar matrices.
- In ValidatedTag, interval diagonal coefficients are independent: accepting the sorted result requires strictly ordered disjoint intervals, except for *identical singleton* diagonal values (provably equal for all represented matrices). Equal non-singleton interval descriptions such as [2,3] and [2,3] do not prove the values coincide or have fixed ordering; correctly report INCONCLUSIVE in that case.
- Extended ApproximateTag Jacobi to symmetric matrices with repeated eigenvalues. Removed an unnecessary spectral-gap rejection: the Jacobi rotation basis remains orthogonal even for higher-multiplicity eigenspaces. Added a native DP/MP orthogonality residual check as an explicit numerical invertibility safeguard alongside the existing eigenpair residual.
- Preserved the original symmetric 3x3 repeated-spectrum test previously asserted inconclusive: it now requires SUCCESS with eigenvalues {0,2,2}. Added an independent defective upper-triangular Jordan block requiring INCONCLUSIVE, so real-spectrum defectiveness is still covered rather than hidden.
- Added DP and MP regression scenarios for repeated 2x2 diagonal scalar matrices, a 4x4 validated diagonal family with repeated *point* eigenvalues plus one uncertain interval, an unorderable diagonal interval family, and a non-diagonal symmetric matrix with spectrum {1,3,3}; check eigenvalue ordering, eigenvector residuals and Jacobi dispatch.
- Non-diagonal triangular repeated spectra and general dense nonsymmetric repeated spectra remain INCONCLUSIVE unless a separate validated/full-rank eigenspace method is added. This deliberately avoids claiming a basis without independent eigenvectors.
- **Verification:** source edit and regression scenarios inspected; compilation and tests for this new commit remain pending. No profiling or benchmark runtime work, no changes to Format.

### 2026-11-27 — repeated eigenvalues in non-diagonal triangular matrices (this commit)

- Baseline: `49a8c75` passes Ubuntu, macOS, Windows, Sanitizers and Coverage; `check-format` is deliberately deferred. The user prioritizes essential eigendecomposition capabilities over profiling/optimization work.
- Generalized the upper- and lower-triangular eigenvector back/forward substitution to allow **exactly repeated diagonal eigenvalues** when each row equation associated with a zero denominator is identically satisfied. Assign zero to that free component and retain a fixed unit pivot in the corresponding eigenvector. A nonzero or unprovably-zero numerator returns `INCONCLUSIVE`; no division by zero is attempted.
- Sorting allows identical real diagonal scalars in ApproximateTag and provably equal **singleton** diagonal coefficients in ValidatedTag. Distinct interval entries must be strictly separated. Repeated nonsingleton interval descriptions are independent uncertain coefficients and cannot be assumed equal or consistently ordered across the matrix family.
- The resulting eigenvector matrix is a column permutation of a **unit triangular** matrix with exact structural zeros, therefore it is invertible for each member of the interval family, even when the eigenvalues repeat. A defective Jordan block still returns INCONCLUSIVE.
- Added DP/MP tests of a non-diagonal upper triangular 3x3 with diagonal (1,2,1), a compatible repeated eigenpair, independently known spectrum (1,1,2), residuals in both paradigms, transpose (lower triangular), ambiguous repeated interval descriptions, and an incompatible Jordan block. Every successful validated result still encloses a consistently ordered real eigenbasis for every matrix member.
- No changes to benchmark execution policy (Release CI **compiles** benchmark binaries only), no new tooling or dependencies, and no new profiling work.
- **Verification:** source edit and proof obligations reviewed; CI for this new commit remains pending. Any errors found must be classified as RTL/implementation issues versus test expectation rather than suppressed to make regression pass.

### 2026-11-27 — consolidate eigendecomposition template definitions (this commit)

- Following a request to simplify the file structure, consolidated the entire Householder QR/approximate Newton implementation from `include/algebra/eigen_nonsymmetric.tpl.hpp` into the existing `include/algebra/eigen_decomposition.tpl.hpp`. Deleted the additional template header and its internal include.
- Preserved the QR function definition and its Doxygen contract **verbatim** (apart from surrounding file namespace/guard wrappers). The primary template already includes `<exception>` and `<vector>`, so no new includes are needed.
- The public `eigen_decomposition.hpp`, algorithm dispatch, test source, build targets, and benchmarks remain unchanged. This is a pure source-organization change, not a new spectral algorithm or an optimization.
- Verified before commit that the independent QR definition appears exactly once in the merged source and the secondary include reference is gone. CI compilation/test outcomes for the resulting commit must still be checked.
- Prior commit `2d8b694` had passed Ubuntu, macOS, Windows, Sanitizers and Coverage; `check-format` stays intentionally deferred until pre-merge. Benchmark executables are compiled in Release CI on all three OSes but **never run in CI**.

### 2026-11-27 — display nanoseconds as integers in eigen benchmarks (this commit)

- User requirement: the eigendecomposition benchmark should not display fractional nanoseconds. Keep measurement precision and timing statistics unchanged; round **only the displayed output** to the nearest integer nanosecond (and variance to integer ns^2).
- `Benchmarking::Report` uses `Table::number(value)`, which defaults to three decimal places. Avoid changing the shared benchmarking library or the arithmetic of the stored samples. Use `Benchmarking::Table` directly for this benchmark, retaining its same configurable summary-field headings, grouping, column order and formatting through `LOGGING_PRINTLN`.
- The normal table and optional profiling columns use `Table::number(value,0u)`; CSV writes `std::fixed << std::setprecision(0)` for all time statistics and optional stage medians. The existing CSV headers and schema are unchanged; count remains an integer and variance retains the `ns^2` unit.
- Documented that display rounding does not imply nanosecond-level measurement accuracy; only the report representation is changed. No workload definitions, timers, acceptance criteria or algorithms altered.
- **Verification:** inspected source and the shared `Report`/`Table` formatting behavior, preserved requested summary-field support and logging. The new commit CI is pending. `check-format` remains deferred; Release CI compiles benchmark binaries without executing them.

### 2026-11-27 — exact structural block-diagonal eigendecomposition (this commit)

- The existing 1x1, 2x2, triangular and dense solvers did not exploit an **interleaved block-diagonal structure**; in ValidatedTag they could fail to certify repeated simple block eigenvalues even when a complete basis exists.
- Added a fast structural graph decomposition for matrices of order >2: vertices are coordinate indices, and an edge joins i and j whenever either A_ij or A_ji is **not provably exactly zero**. Connected components identify independent principal blocks up to a simultaneous coordinate permutation. Nonsingleton zero-containing interval entries are **not** treated as structural zeros.
- Each block is solved recursively by the **existing** eigendecomposition API (no new numerical solver). Assemble its eigenvectors into their original coordinate support with exact zero entries outside the block, then globally sort eigenpairs. For ValidatedTag, sorting requires provably disjoint eigenvalue enclosures or identical singleton values; overlapping uncertain eigenvalues remain INCONCLUSIVE.
- Mathematical invariant: each successful block supplies an invertible eigenvector basis for every matrix in its input family; the global basis is block diagonal up to row/column permutations and is therefore invertible. Repeated eigenvalues from separate blocks are allowed, while a defective block fails as before.
- Added a BLOCK_DIAGONAL diagnostic strategy identifier and the corresponding string for optional manual benchmark reporting. Block-local candidate/verification timings are accumulated if requested; no new profiling functionality or benchmark workload added.
- Added DP/MP Approximate/Validated regressions with two **interleaved** coupled 2x2 blocks (spectrum {2,4,4,6}), interval-uncertain separated blocks, overlapping interval-order uncertainty, and a defective Jordan subblock. Check globally sorted eigenvalues, eigenvector residuals and the actual dispatch strategy.
- **Verification:** source checked and submitted for CI. Do not claim tests passed until workflows for the new commit complete. Current baseline c06e170 passed Ubuntu/macOS/Sanitizers, while Windows and Coverage were still running when reviewed. `check-format` remains deferred; benchmark binaries are compiled but **never executed in CI**.

- Preserved independent Jacobi regression coverage: replaced the earlier block-decomposable symmetric repeated-spectrum example with fully connected `4I-J` (eigenvalues `1,4,4`). Otherwise the new structural dispatch would bypass Jacobi, unintentionally reducing algorithm-specific test coverage.

### 2026-11-27 — fix stale repeated-eigenvalue expectation after block decomposition (this commit)

- CI regression introduced by `bb896e2`: Release builds succeed, but Debug/Coverage/Sanitizer tests fail in `eigen_decomposition`. The failing assertion is `repeated_spectrum.status == INCONCLUSIVE` in `test_real_examples`; the preceding `validated(v).status == INCONCLUSIVE` in the failure report is the *last successful check context*, not evidence of a second failed assertion.
- The tested matrix `[[0,1,0],[1,0,0],[0,0,1]]` has eigenvalues `{-1,1,1}` and is fully diagonalizable. The independent `2x2` and `1x1` blocks yield a complete eigenvector basis, including the duplicated value `1`. `SUCCESS` is therefore the mathematically correct outcome after adding blockwise decomposition; demanding `INCONCLUSIVE` was an obsolete test expectation.
- Updated this **same** test matrix to require `SUCCESS`, presence of data, rigorous enclosures for its independently known sorted eigenvalues and interval overlap of the full eigenpair residuals `A*V` and `V*D`. Did not weaken any numerical tolerance, discard the scenario, or modify the block decomposition algorithm.
- Added a separate defective block-family counterexample `[[0,1,0],[0,0,0],[0,0,1]]` with Jordan eigenvalue zero; this still requires `INCONCLUSIVE` and no data. The existing overlapping independent-interval ordering test `validated(v)` remains untouched.
- **Verification:** inspected source and mathematical spectra, committed corrected assertions and journal. Compilation and test results for the corrected commit remain pending. `check-format` stays deferred; CI Release compiles, but never runs, benchmarks.

### 2026-11-27 — expose the diagonal eigenvalue matrix in eigendecomposition results (this commit)

- Baseline `8d86fcb` has passed Ubuntu, macOS, Windows, Sanitizers and Coverage. `check-format` remains deferred until pre-merge; external Codecov statuses do not change functional test results.
- `EigenDecompositionData<P,PR>` now exposes `diagonal_matrix()`, returning the existing `DiagonalMatrix<NumericType>` constructed from the already-computed, sorted `eigenvalues`. No new representation, solver or dependency is introduced, and no numerical algorithms are changed.
- This supplies the missing concrete diagonal factor `D` alongside eigenvector matrix `V` for the basic `A V = V D` eigendecomposition relation. The user can use existing `inverse(V)` independently if appropriate, rather than always computing a potentially over-enclosed interval inverse as part of an otherwise successful validated result.
- In ValidatedTag, `D` contains rigorous eigenvalue **enclosures**; `AV` and `VD` evaluated independently in interval arithmetic may differ as intervals because each entry can be over-enclosed. Do not claim equality of interval matrices; the pointwise eigenpair relation holds for coherent selections for every input member.
- Added DP/MP test coverage for the new API in approximate 2x2, validated 2x2, and interleaved block-diagonal 4x4 with repeated eigenvalues; check exact off-diagonal zeros and the appropriate approximate residual or rigorous interval overlap of `AV` and `VD`.
- **Verification:** source-level interface and test changes committed; CI for this new commit remains pending. Benchmark binaries are still **compiled only** in Release CI across all platforms and never executed there.

### 2026-11-27 — optionally recover inverse eigenvector basis for full diagonal factorization (this commit)

- Follow-up to the `diagonal_matrix()` API: add `EigenDecompositionData<P,PR>::inverse_eigenvector_matrix()` returning `std::optional<Matrix<NumericType>>`. It delegates to the existing matrix `inverse()` implementation and catches only `SingularMatrixException` to represent failure to establish a nonzero pivot. No second factorization type, new numerical algorithm, or cached inverse is introduced.
- For ApproximateTag, the returned basis inverse is numerical and not a certificate. For ValidatedTag, successful interval elimination encloses the inverses of all independently varying interval matrices described by `eigenvectors`, so callers can form an enclosing `V*D*inverse(V)`. This is **stronger** than the pointwise, correlated eigenbasis-existence guarantee of `SUCCESS`; an inconclusive interval inverse does not invalidate the eigenpair certificate.
- Tests in DP and MP require a usable inverse for the exact 2x2 approximate and validated symmetric examples and the interleaved 4x4 repeated-spectrum validated example. Approximate reconstructions are checked against numerical tolerances; validated reconstructions are required to contain every coefficient of the original interval matrix. The pre-existing `A*V=V*D` tests are retained.
- No formatting changes or benchmarking changes. Benchmark executables are still **only compiled** in Ubuntu/macOS/Windows Release CI, never executed. `check-format` remains postponed until the merge.
- **Verification:** source and mathematical contract reviewed and committed; actual new-commit CI compilation and tests remain to be checked.

### 2026-11-27 — reuse actual-matrix QR eigenpairs as validated Krawczyk seeds (this commit)

- Baseline `e96bd52` has passed all functional workflows (Ubuntu, macOS, Windows, Sanitizers and Coverage). Format remains intentionally deferred, with external Codecov checks separate from tests.
- Previously, validated dense Krawczyk certification constructed its **only** approximate seed by applying symmetric Jacobi to `(mid(A)+mid(A)^T)/2`. For a non-normal matrix this solves a different eigenproblem, and Newton can converge poorly or collapse multiple eigenpair seeds.
- The Krawczyk fallback now tries `EigenDecomposition<ApproximateTag,PR>` on the **actual midpoint matrix first**: symmetric inputs automatically use Jacobi, and nonsymmetric inputs use the already available real shifted-QR/Newton backend. The resulting approximate eigenpairs are only candidate centers, never certification evidence.
- Refactored the existing outward-rounded normalized eigenpair Krawczyk checks into a local attempt function (same equations, interval Jacobian, norm/contraction predicate, radius expansion and pairwise eigenvalue separation). If actual-midpoint candidates fail or cannot be certified, retry the **original** symmetrized-Jacobi candidates; skip the duplicate attempt when the midpoint is exactly symmetric. This preserves previously supported cases instead of changing one heuristic for another.
- The mathematical `SUCCESS` conditions are unchanged: every accepted root is certified for every matrix represented by the original interval input, with a consistent sorted simple real spectrum and pointwise invertible eigenbasis. Failures remain `INCONCLUSIVE`.
- Added DP/MP validated regression tests with a fully coupled nonsymmetric 3x3 matrix with independently known spectrum `{1,3,5}`, and its interval family formed by a small nonnegative perturbation. Check eigenvalue enclosures and interval `A*V`/`V*D` overlap. Existing non-normal 4x4 tests and all prior pathways remain active.
- No new dependencies, separate template files, tuning/profiling work or benchmark executions in CI. Benchmark executables remain **compile-only** on Ubuntu/macOS/Windows Release.
- **Verification:** source semantics and fallback preservation reviewed; CI results for this commit remain pending. `check-format` remains deferred until pre-merge.

### 2026-11-27 — eliminate duplicate local half in the QR/Krawczyk fallback (this commit)

- The `79c2148` eigendecomposition refactor introduced two declarations of `const AX half` in the same `_krawczyk` function scope: one at the start of the method and another just before the symmetrized-Jacobi fallback. Clang fails to compile the header with `redefinition of 'half'` (reported by the macOS Release benchmark build, and other functional CI builds).
- Removed the **unused outer declaration**. The sole retained `half` is declared immediately before constructing the fallback symmetric seed, its only use. The actual-midpoint QR candidate path, the fallback and all outward-rounded proof obligations are unchanged.
- Verified the new source contains exactly one local `const AX half` in `_krawczyk`, with the fallback use intact; this is a **C++ compilation fix**, not a test-expectation or algorithm change. New CI outcomes must be checked before claiming compilation/test success.
- `check-format` remains deferred until pre-merge. Release CI compiles benchmark binaries but never executes them.

### 2026-11-27 — detect triangular matrices up to a coordinate permutation (this commit)

- Baseline `b2a225a`: Ubuntu/macOS/Sanitizers/Coverage jobs succeeded, Windows x64 succeeded and ARM64 was still running at inspection; `check-format` deliberately deferred. The `half` redeclaration regression is resolved in compiled workflows.
- Essential structural eigendecomposition extension: for a nontriangular connected matrix, build a **directed graph** of structurally nonzero off-diagonal coefficients (row i -> column j). Perform a deterministic topological sort. If the graph is acyclic, the simultaneous permutation of rows and columns gives an **upper triangular** matrix. Pass it to the existing triangular solver; map eigenvector rows back to original coordinates, retaining its sorted eigenvalue order.
- Structural absence requires a *provably exact zero*. A floating interval containing both zero and nonzero values still contributes a directed edge. The graph test never assumes that distinct uncertain coefficients are correlated.
- The triangular solver already certifies a unit-pivot eigenvector matrix (pointwise for every ValidatedTag input member) and handles compatible repeated point eigenvalues. A cyclic pattern, or an inconclusive permuted triangular solve, **falls through** to the existing dense algorithms without claiming a missing eigenbasis.
- Reused the existing TRIANGULAR diagnostics status and propagated per-block timing counters when requested. No new solver classes, algorithms, dependencies or profiling features were introduced.
- Added DP/MP Approximate/Validated regression for a connected, interleaved 4x4 acyclic matrix that is neither upper nor lower triangular in stored order, with known spectrum `{1,3,5,7}`. The validated test also uses an uncertain off-diagonal coefficient `[0,1]`, checking the same spectral ordering and interval overlap of `A*V` and `V*D`.
- **Verification:** source and graph-permutation argument inspected; new CI results pending. Benchmark executables are compiled only in Release CI on Ubuntu/macOS/Windows and **never run**. `check-format` remains deferred until pre-merge.

### 2026-11-27 — restore genuine dense QR coverage after structural triangularization (this commit)

- CI regression for `180f998`: Debug/Coverage tests fail at `diagnostics.strategy == SHIFTED_QR`. The old nonnormal 4x4 example `[[1,2,0,0],[0,3,0,0],[-8,10,5,0],[-6,2,0,7]]` has an **acyclic** directed off-diagonal graph. The new `TRIANGULAR` path is mathematically correct; requiring QR was an obsolete test expectation.
- Preserved this exact original 4x4 matrix, its known spectrum `{1,3,5,7}` and residual checks in the approximate test, now explicitly expecting `TRIANGULAR` instead of removing or weakening its coverage.
- Introduced a separate truly cyclic 4x4 nonnormal matrix `[[3,0,-2,-2],[-1,5,0,-2],[-2,2,4,-1],[0,-2,-1,4]]` from `S*(H*diag(1,3,5,7)*H^T/4)*S^{-1}`, with integer order-four Hadamard `H` and unimodular shear `S=I+E_01`. The spectrum is independently known, and its directed graph contains cycles; it must use actual `SHIFTED_QR`, preserving high-order QR coverage.
- The original validated `test_krawczyk_nonnormal` and labeled nonnormal benchmark also used that acyclic fixture, silently bypassing dense validation. Both now use the cyclic matrix and existing interval perturbation. Validated regression requires `SUCCESS`, real ordered eigenvalue enclosures, interval residual compatibility, and a nonstructural strategy, but does not require a particular valid dense certificate (Gershgorin, verified similarity or Krawczyk).
- Benchmarks continue to be **compiled but never executed** in CI. No numerical solver behavior, precision tolerance, dependencies or formatting workflows have changed. `check-format` remains deferred.
- **Verification:** graph and characteristic spectrum checked by exact similarity; new CI compilation and tests pending, including possible newly revealed failures from the previously bypassed dense path. Do not suppress such failures to make tests pass.

### 2026-11-27 — certify nonnormal dense matrices using an actual-midpoint QR similarity (this commit)

- Regression on `b6089e7`: the genuinely cyclic nonnormal 4x4 matrix (real simple spectrum `{1,3,5,7}`) returned `INCONCLUSIVE` from the validated eigendecomposition, not a wrong eigenvalue. CI Debug/Coverage/Sanitizer failures occurred at `result.status == SUCCESS` in the restored dense-regression test; merely changing that assertion would conceal an essential unsupported case.
- Found the existing `_preconditioned` certificate was **only** seeded by eigenvectors of the symmetric part of the midpoint; for a nonnormal input those eigenvectors need not even approximate the actual eigenspaces. The other dense fallback `_krawczyk` already tries QR eigenvectors of the actual midpoint, but a valid QR eigenbasis can also support a simpler, rigorous **similarity plus Gershgorin** certificate.
- Extended `_preconditioned` to propose a fixed basis from the **actual midpoint eigendecomposition first** (nonsymmetric shifted QR / symmetric Jacobi). Reused its existing validated interval Gauss–Jordan check of `S^{-1}`, the verified enclosure of `B=S^{-1}AS`, and the strictly separated Gershgorin/eigenvector basis checks, unchanged, as a reusable local attempt. This yields a pointwise real diagonalization certificate for each represented interval matrix when successful.
- Preserved the original symmetrized-Jacobi candidate as a fallback if QR fails to converge or its similarity certificate is inconclusive; symmetric midpoints skip the redundant repeat. The existing `_krawczyk` fallback remains in place after both similarity strategies fail. No approximation is mistaken for certification; all failed attempts return `INCONCLUSIVE`.
- Renamed the misleading test label `test_krawczyk_nonnormal` to `test_dense_nonnormal_validated`, because the required contract is successful **dense validated spectral certification**, not necessarily a particular internal certificate. Preserved the cyclic point and uncertain-family matrices, exact eigenvalue enclosures, strict separation, true eigenpair residual overlap, and nonstructural strategy guard. No successful test claimed until CI completes.
- **Verification:** independent numerical calculation shows the cyclic 4x4 has simple eigenvalues and a well-conditioned point eigenvector basis; the strict interval criterion and CI still decide success. No benchmark algorithm changes, no benchmark execution in Release CI, no format changes; `check-format` remains deferred.


### 2026-11-27 — structural strongly-connected block-triangular eigendecomposition (this commit)

- Continued on `eigen` after `0b58156` (actual-midpoint QR similarity). At inspection the baseline functional macOS/Coverage/Ubuntu jobs passed; remaining Windows/Sanitizer jobs were still running. Format remains intentionally deferred until pre-merge.
- Added an exact structural decomposition between the independent-block and permutation-triangular paths. Compute the reachability relation of the directed graph of provably nonzero coefficients, partition into strongly connected components (SCCs), and topologically order the condensation DAG. An interval entry merely containing zero still creates an edge. This new path applies only to genuine cyclic SCCs coupled to other components; the existing all-singleton permutation-triangular and completely connected dense paths remain intact.
- Solve each SCC's principal diagonal block recursively with the existing approximate/validated eigendecomposition. For an eigenpair from block k, lift its eigenvector in reverse topological order by solving `(B_i-lambda I)v_i = -sum_{j>i} A_ij v_j` in every preceding block. Use the existing `inverse` implementation; in ValidatedTag, the outward-rounded interval inverses enclose each inverse for every represented coefficient/eigenvalue selection when pivots are provably nonzero. Any inconclusive block, pivot, or globally ambiguous eigenvalue ordering falls through to the original dense solver.
- The lifted eigenvector matrix is upper block triangular after permutation, with invertible diagonal eigenvector bases from the recursive successes. Thus its columns form an invertible basis for every input family member. Globally sort the eigenpairs only when validated intervals are strictly separated or identical exact singletons, preserving the existing fixed-order contract. Introduced the diagnostic `BLOCK_TRIANGULAR` and its benchmark string; existing Release CI continues to compile benchmarks without running them.
- Added DP and MP tests for an interleaved, directed-cyclic 2x2 block coupled to a 1x1 block with known spectrum {1,3,5}. The tests require actual structural strategy selection, approximate residuals, rigorous eigenvalue enclosures and interval eigenpair residual compatibility for both a point matrix and a nontrivial uncertain family. A coupled defective repeated-eigenvalue example must remain inconclusive.
- **Verification:** algorithms and changes inspected at source level. Local C++ compilation was unavailable because the environment cannot resolve github.com to retrieve the repository and dependencies. New CI results must be checked; no passing new-commit tests are claimed. `check-format` remains deferred.


### 2026-11-27 — certify compatible repeated eigenvalues across coupled cyclic blocks (this commit)

- Baseline `7ed2f18`: Ubuntu Release/Debug, macOS Release/Debug/Coverage, Ubuntu Coverage and ThreadSanitizer passed at inspection; Windows and AddressSanitizer/UBSan were still in progress. No new-commit CI results claimed. `check-format` continues to be intentionally deferred.
- The original block-triangular lifting used `inverse(B_i-lambda I)`, which necessarily fails at a repeated eigenvalue shared with an earlier diagonal SCC, even when the coupling is exactly compatible and the full matrix remains diagonalizable. Such an inverse failure is not a proof of defectiveness.
- After that inverse throws `SingularMatrixException`, reuse the already established eigenbasis `V_i` of the earlier diagonal block and attempt an independent interval-certified `V_i^{-1}`. Express the right-hand side in spectral coordinates `z=V_i^{-1} rhs`. For each eigenvalue `mu_j`, solve `(mu_j-lambda)w_j=-z_j` if the difference is provably nonzero; if the difference is exactly zero and the corresponding `z_j` is exactly zero, set the free spectral component to zero. Otherwise the certificate is inconclusive and existing dense fallback is preserved.
- For every interval member, the spectral basis `V_i` is invertible and the certified components enclose a coherent solution of the shifted block equations. All earlier block vectors remain lifted without changing their own diagonal block coordinates, retaining the invertibility proof for the assembled block-upper-triangular eigenbasis. This method does not infer algebraic correlation between independently represented interval entries.
- Added DP/MP approximate and validated tests for the coupled interleaved 3x3 matrix `[[2,1,1],[0,3,0],[1,-1,2]]` with independently known real eigenvalues `{1,3,3}`, eigenpair residuals and full-basis inverses. Changing the coupling to `[1,1]` creates a defective eigenvalue 3 and must yield `INCONCLUSIVE`; a nontrivial independent coefficient interval family containing defective members must also remain inconclusive. Retained previous distinct-spectrum SCC, QR and Gershgorin scenarios.
- **Verification:** mathematical projection equations and source edits reviewed; CI compilation and tests for this commit remain pending. Benchmarks are compiled only, not executed in Release CI; no formatting changes.


### 2026-11-27 — nonfinite input guard and multi-SCC block-triangular regression (this commit)

- Baseline `e561368`: all Ubuntu, macOS, Coverage, AddressSanitizer/UBSan and ThreadSanitizer jobs passed at inspection. Windows x86 Debug/Release and ARM64 Release passed; ARM64 Debug remained in progress. Codecov patch/project thresholds are coverage reports, not functional regression. `check-format` is still deliberately deferred until pre-merge.
- Added a paradigm-native finite-coefficient check for every input matrix before any structural shortcut. Previously order-one analytic matrices and off-diagonal entries traversing block-triangular paths could permit NaN/infinity with an invalid `SUCCESS`. Positive order matrices with nonfinite approximate coefficients or interval endpoints now return `INCONCLUSIVE` without data. Preserve zero-order `UNSUPPORTED_DIMENSION` and the nonsquare precondition.
- Added DP/MP tests constructing raw NaN and infinity for both approximate and validated 1x1 inputs, and off-diagonal infinity in a 3x3 triangular input. They exercise early rejection and guard against apparently successful invalid eigenvector bases.
- Added DP/MP tests for **two distinct cyclic SCCs** with one-way coupling, rather than only a 2x2 cycle plus singleton: interleaved 4x4 `[[2,1,1,2],[0,6,0,1],[1,3,2,4],[0,1,0,6]]`, whose exact characteristic polynomial was independently checked as `(lambda-1)(lambda-3)(lambda-5)(lambda-7)`. Assert real ordered eigenpairs, actual BLOCK_TRIANGULAR dispatch, approximate residuals and validated interval eigenpair residual overlap. Perturb both the first SCC diagonal and one cross-coupling entry independently by `[0,2^-20]` to exercise family-wide certification.
- No changes to benchmark execution (build-only in Release CI), numerical heuristics, dependencies or formatting. **Verification:** exact spectrum checked independently and source edits reviewed; new commit CI compilation and tests remain pending. No unrun test is claimed as passed.


### 2026-11-27 — reject nonfinite eigenpairs from finite-input arithmetic overflow (this commit)

- Baseline `0b74c45` passed every functional CI workflow: Ubuntu, macOS, Windows, Coverage, and Sanitizers. Codecov coverage thresholds and the intentionally deferred `check-format` remain separate from functional results.
- Identified a distinct postcondition gap after the finite-input precheck: the upper/lower triangular and structurally lifted paths can overflow while calculating normalized eigenvectors from otherwise finite real matrix entries, and could then return `SUCCESS` containing infinity. Finite inputs do not imply finite intermediate arithmetic.
- Moved the existing strategy dispatch into a private `_attempt` while preserving the two public `operator()` overloads and recursive calls through the public interface. The public diagnostics overload now checks every successful eigenvalue and eigenvector element for finite raw values (approximate) or finite interval endpoints (validated), as well as result dimensions. If any output is invalid it returns `INCONCLUSIVE` with no data. It makes no new certification assertion, and does not change the numerical solvers, ordering, or selection of approximate candidates.
- Added a binary64 approximate/validated 3x3 upper-triangular regression with exact diagonal `{1,2,3}` and off-diagonal couplings `2^600` at (0,1) and (1,2). Input coefficients are provably finite while unit-pivot back-substitution for eigenvalue 3 requires a component of magnitude `2^1200/2`, exceeding binary64 range. Both paradigms must now report `INCONCLUSIVE`, retain `TRIANGULAR` diagnostics, and return no data. Reducing both couplings to one must preserve `SUCCESS`. The positive-case DP/MP structural and dense regressions remain intact.
- **Verification:** source edit and independent exponent reasoning checked. Full new-commit CI compilation/test results have not yet been observed; no claim of passing this new regression is made. Release CI compiles benchmark binaries without running them; `check-format` remains deferred.


### 2026-11-27 — separate independent complex Spectrum from diagonalization (this commit)

- Reviewed the 919e3da CI jobs: all functional workflows other than Ubuntu GCC 13 Release passed, and that job stopped at Linux dependency installation before configuration/compilation. No unverified claim of a solver regression, no weakening of numerical tests, and no changes to Format.
- Superseded the original real-diagonalization-first API plan with Spectrum, SchurDecomposition and EigenDecomposition as distinct requested mathematical results. Explicitly removed all backwards-compatibility commitments: there are no API users yet. Kept the earlier journal sections as labelled historical records, not current mandates.
- Implemented initial independent header-only Spectrum<P,PR> and SpectrumResult<P,PR>. Output is a vector of Ariadne Complex<FloatType<P,PR>> values, including algebraic repetitions; no eigenvector basis is needed. 1x1, 2x2 and upper/lower triangular real inputs of arbitrary positive order are supported. For 2x2, use a stable real quadratic branch or negative-discriminant conjugate pair, and allow exactly repeated roots of defective matrices. Reject uncertain discriminant sign or nonfinite results with INCONCLUSIVE.
- For ValidatedTag, every reported complex rectangle encloses the corresponding root for every real matrix of the interval family, using outward arithmetic. Triangular eigenvalues preserve diagonal order (including independently overlapping diagonal intervals), rather than imposing impossible family-wide sorting constraints. Other dense orders >=3 remain INCONCLUSIVE until Schur is implemented: not a false claim of generality.
- Registered standalone DP/MP tests for a rotation with conjugate complex eigenvalues (including an interval family), a defective Jordan 2x2 with a repeated real eigenvalue, higher-order defective triangular spectra and overlapping uncertain diagonal intervals, plus honest inconclusive dense fallback. Exposed Spectrum in the Algebra umbrella header. Existing EigenDecomposition remains available for now only as a working internal starting point, **not** for compatibility preservation.
- **Verification:** source and mathematical contracts inspected. New CI compile/test results are pending; do not claim success yet. Release benchmark targets remain compile-only in CI, check-format deferred.

### 2026-11-27 — fix Spectrum local name shadowing (this commit)

- The first Spectrum CI run for `38ec509` reported a Clang 16 compilation failure in `tests/test_spectrum.cpp`: the local `const X two` in `include/algebra/spectrum.hpp` shadows the pre-existing namespace-level `Ariadne::two` from `numeric/twoexp.hpp`. With `-Wshadow -Werror`, this is a hard compilation failure before tests execute.
- Renamed the one local constant to `two_value` and updated all seven division uses in the 2x2 spectrum formulas. The arithmetic, interval branches, public contracts and test expectations are unchanged. No weakening of warnings or test coverage.
- **Verification:** inspected the reported compiler error and checked that the modified source contains no remaining standalone `two` identifiers in `spectrum.hpp`, with eight occurrences of `two_value` (one definition, seven uses). New CI has not yet been observed; this is a source-level correction, not a claimed passing build. Formatting remains deferred until pre-merge.

### 2026-11-27 — fix GCC installed umbrella-header compilation (this commit)

- Investigated the attached GitHub Actions log archive from the 2364be5 pull-request merge build. The two failing jobs were Ubuntu x86 and ARM GCC 13 Release. Debug test_spectrum and other tests passed, and Release library/benchmark builds completed, but the installed tutorial failed compiling tutorial_algebra.cpp through ariadne-algebra.hpp.
- Root diagnostic: after many algebra headers, the new late inclusion of algebra/spectrum.hpp causes numeric/complex.hpp to re-declare Ariadne::TranscendentalType already seen in numeric/operators.hpp. GCC diagnoses the dependent decltype(sin(declval<X>())) aliases as conflicting due to different unqualified name lookup contexts. The many subsequent std::deque and std::mt19937 failures are cascades, not separate STL defects. This is independent of the fixed local two shadow.
- Move algebra/spectrum.hpp to the very first inclusion in ariadne-algebra.hpp, matching the minimal standalone test_spectrum.cpp inclusion context that already compiled. Do not edit the numeric dependency, silence diagnostics, hide Spectrum or loosen any test. This is only an include-order repair.
- **Verification:** GCC error location and previous successful Debug tests inspected; the source now includes Spectrum exactly once, before algebra/algebra.hpp. New CI installed-tutorial results are pending; do not claim the upstream numeric alias duplication itself is eliminated. check-format remains deferred.


### 2026-11-28 — structural 1x1/2x2 spectra beyond triangular matrices (this commit)

- Verified baseline `58fbbc4`: all functional GitHub workflows (Ubuntu, macOS, Windows, Coverage and Sanitizers) passed. In particular the installed tutorial compiled on Ubuntu GCC 13 Release, confirming the umbrella-header ordering fix. Format remains intentionally deferred; external Codecov threshold failures do not signify a failed C++ test.
- Extended `Spectrum<P,PR>` from stored triangular matrices to **permutationally upper block-triangular matrices whose strongly connected components have order at most two**. The directed structural graph has edge i->j for every input coefficient not *provably exactly zero* (interval coefficients containing zero are not treated as absent). Compute transitive reachability, partition into SCCs, build their condensation DAG and use a deterministic topological ordering. The matrix is permutation-similar to upper block-triangular form, so its characteristic polynomial equals the product of diagonal-block polynomials; its spectrum is their multiset union, with algebraic multiplicity, independent of inter-block couplings.
- For singleton SCCs, take the corresponding diagonal real value. For doubleton SCCs, recursively apply the existing 2x2 real/complex spectral formulas, including rigorous outward interval enclosures and repeated roots. Do not construct or require any eigenvector basis. Any SCC of order three or more or an inconclusive 2x2 certificate remains `INCONCLUSIVE`; **no claim of general dense Schur or QR**.
- Added DP/MP approximate and validated regressions on a non-triangular interleaved 4x4 with rotation SCC {0,2}, scalar nodes {1,3}, one-way cross-block couplings, and independent exact spectrum {-i,+i,3,5}. Validate both point and interval-family enclosures with uncertainty in a cross edge and one diagonal coefficient. Adding a closing directed edge to make the full graph strongly connected must return `INCONCLUSIVE` instead of falsely reusing the structural eigenvalues.
- Existing dense-cycle refusal, Jordan repetitions, triangular spectra and diagonalization tests are preserved. Implementation only changes Spectrum's previously inconclusive structural branch; SchurDecomposition is the next algorithmic phase, which must implement real Hessenberg/Francis QR for general dense inputs. Benchmark execution remains disabled in CI.
- **Verification:** source and graph/spectral arguments reviewed, previous CI confirmed; new-commit build and tests pending. No new passing test claimed.

### 2026-11-28 — initial approximate real Schur and dense complex Spectrum (this commit)

- Baseline `0a217dd` passed all functional workflows (Ubuntu, macOS, Windows, Coverage, Sanitizers) on GCC/Clang/MSVC, including the DP/MP SCC tests. Formatting remains deliberately deferred, and the Codecov coverage threshold failures are separate from functional tests.
- Introduced the independent `SchurDecomposition<P,PR>` and `SchurStatus` API. On approximate success the factors `Q`, `T` and a 1x1/2x2 block partition satisfy a checked small residual `A*Q-Q*T` and approximate orthogonality of `Q`. Schur does not demand an eigenvector basis. A 2x2 block must have a strictly negative discriminant (conjugate complex pair). A successful result is an **approximation**, never a certificate.
- Initial numerical backend: Householder similarity to Hessenberg form, then explicitly shifted real Householder QR; a single Wilkinson shift is used for real trailing roots and a real quadratic polynomial `T*T-tr(B)*T+det(B)*I` for a conjugate trailing pair. QR transformations are applied to the entire matrix (including coupling into deflated blocks), and accumulated in `Q`. Small fill-in is removed only within precision-dependent tolerances; all outputs must meet finite, quasi-triangular, orthogonality and original-matrix residual checks. This is explicitly **not yet** an optimized implicit Francis QR implementation. Bounded iterations and numerical problems return `INCONCLUSIVE`.
- `Spectrum<ApproximateTag,PR>` now tries Schur for strongly connected components larger than two instead of treating every such dense matrix as automatically inconclusive. The accepted Schur block spectrum is assembled from scalar diagonal elements and the existing 2x2 complex formula. `ValidatedTag` for general dense complex spectra and Schur decomposition remains `INCONCLUSIVE`; the approximate Q does not certify interval spectra.
- Added a standalone `test_schur_decomposition` target with DP/MP regression: fully strongly connected integer 3x3 `[[4,5,0],[3,4,-4],[-4,0,1]]`, with independently known characteristic polynomial `(lambda-9)(lambda^2+9)`, requiring one complex pair, one real value, orthogonality and low residual. Also tests fully coupled symmetric 3x3 with real Schur blocks, defective triangular Jordan block, and honest validated inconclusiveness. Registered target in tests/CMakeLists.txt. Existing spectral and diagonalization regressions preserved.
- **Verification:** source-level inspection and numerical double-shift prototype only. New C++ compilation and CI tests for this commit remain pending; no passed build is claimed. The next task is to assess CI, then replace explicit polynomial QR by robust implicit Francis bulge chasing and increase conditioning/cluster coverage.


### 2026-11-28 — Doxygen method derivations and historical provenance (this commit)

- Responding to the requirement that API-level documentation explain **which methods are used, where they come from, and how the code implements them**, expanded Doxygen in schur_decomposition.hpp, spectrum.hpp, eigen_decomposition.hpp, and the _nonsymmetric implementation. Not only the journal: generated Doxygen now exposes the mathematical descriptions where developers discover the API and private algorithm methods.
- Documented **Householder reflections** U=I-2vv^T/(v^Tv), their sign choice and left/right orthogonal similarities, and the **upper Hessenberg** structure (zeros below the first subdiagonal). The explanation states why the reduction preserves eigenvalues and why Householder transformations are preferred to nonorthogonal triangularization for numerical stability.
- Derived real shifted QR (T-mu*I=QR; T'=RQ+mu*I), the Wilkinson-type trailing-block shift, and the **explicit** real quadratic double shift p(T)=T^2-tr(B)T+det(B)I for conjugate complex roots. Explicitly distinguished the existing polynomial QR from the still-unimplemented implicit Francis bulge chase; the documentation must not misattribute that algorithm to the present source.
- Explained trailing 1x1/2x2 deflation, the numerical (non-certified) thresholds, finite checks, original-input residual A*Q-Q*T and approximate orthogonality. Distinguished the older EigenDecomposition QR (no preliminary Hessenberg reduction, Newton-refined real eigenbasis) from Schur (complete real/complex spectral form without an eigenbasis).
- Added method-level provenance for Jacobi orthogonal rotations, Gershgorin localization, and Krawczyk interval root inclusion. Bibliographic references checked against the original Householder (JACM 1958), Francis (Computer Journal 1961/1962), Gershgorin (1931), and Krawczyk (Computing 1969) publications; Wilkinson's 1965 treatise identifies the Hessenberg and QR background.
- No algorithms, public signatures, numerical thresholds, tests or build settings changed. Latest baseline ecd31f3 has had failing Debug C++ test steps and Coverage tests in GitHub Actions; these remain an outstanding numerical/test failure to investigate separately, not something that documentation changes can fix. Format and Codecov threshold checks are still separate/deferred.
- **Verification:** reviewed source against the equations and bibliographic records. The actual new commit CI results are pending; no claim of passing those tests.


### 2026-11-28 — audit Schur block boundaries after Debug out-of-range failures (this commit)

- The `ecd31f3` Schur changes introduced a runtime precondition failure in DP `test_schur_decomposition` and `test_spectrum`: an attempted element `(3,2)` on a 3x3 approximate matrix. The input matrix printed by the exception resembles a converged 1x1 plus conjugate 2x2 real Schur form. The preceding Doxygen-only commit `d572e35` left this failure unchanged.
- Explicitly audit the final Schur block partition using `remaining=n-covered` before every subdiagonal probe; handle a trailing singleton separately and check the partition covers exactly n coordinates. Keep the strict real Schur 2x2 negative-discriminant condition, finite checks and residual acceptance criteria.
- Independently validate the Schur factor dimensions and each block width and offset in `Spectrum::_schur_spectrum` **before** indexing any 1x1 or 2x2 block. Require the returned eigenvalue count to equal the input dimension. Extend the existing regression with a block-width boundary invariant, without relaxing any SUCCESS requirements, numerical tolerances or certified contracts.
- **Verification:** the provided error and affected source have been inspected, but a precise throwing statement has not yet been isolated by a stack trace. Consequently this is a defensive boundary fix and regression-strengthening pass, **not** a verified root-cause solution. New CI remains to be evaluated; if failures persist, instrument the execution path to locate the exact offending access. Format remains deferred.

### 2026-11-28 — stop unsigned Schur deflation underflow and update dense-spectrum regressions (this commit)

- The follow-up macOS log after `c8ebe01` still showed an access to matrix element `(0,18446744073709551615)` for a 3x3 approximate Schur matrix. This second coordinate is `SizeType(-1)`, not a meaningful matrix index. The prior Schur-block parsing guards did not solve it. The implicated accesses are to the subdiagonal of an active QR block, such as `T[active-1][active-2]` or `T[active-2][active-3]`. The exact throw-site is not instrumented; no false claim that the stack frame was observed.
- Strengthened the QR active-block termination before any `active-2` expression. For a final real 2x2 block, added a direct orthogonal Schur rotation built from a real eigenvector of the trailing matrix. The two algebraically valid eigenvector candidates are `(b,lambda-a)` and `(lambda-d,c)`; normalize the larger candidate and apply the full similarity to `T` and `Q`, including couplings to deflated blocks. This avoids unnecessary QR iterations on a terminal real 2x2 and also handles exactly repeated defective real roots without requiring a full eigenbasis. The original finite, orthogonality, original-input residual and block-structure acceptance checks remain in force.
- Two `test_spectrum` expectations were stale since `Spectrum<ApproximateTag>` acquired the Schur fallback. A strongly connected 3x3 and a closed-cycle 4x4 must now be tested against independently computed characteristic polynomials and their real/complex eigenvalues, not required to be INCONCLUSIVE just because they are no longer structurally triangular. Replaced each stale assertion with unordered complete-spectrum accuracy checks against independently calculated numerical roots. All `ValidatedTag` dense complex-family tests still require INCONCLUSIVE. These tests preserve intended coverage and increase the required functionality; they do not allow either outcome just to pass.
- Verified the 3x3 characteristic polynomial `x^3-9*x^2+24*x-21` and the 4x4 polynomial `x^4-8*x^3+16*x^2-9*x+16` independently using exact integer symbolic matrices. Test roots are compared with fixed `1e-8` numerical tolerance, not with loosened existing conditions.
- **Verification:** source logic and exact reference polynomials checked, but the compiled C++ regression and CI results for the new commit are pending. If the reported out-of-bounds access remains, locate the precise stack frame before further changes. No format run or benchmark execution, consistent with the documented pre-merge phase.


### 2026-11-28 — implicit Francis bulge chasing for complex shifts (this commit)

- Baseline `758cc5c` passed all 18 functional jobs on Ubuntu, macOS and Windows, including Debug, Release, Coverage and Sanitizers. This confirms correction of the preceding Schur/Spectrum runtime failures. Separate Codecov thresholds remain below target (project 92.84%, patch 80.18% versus 98.83%); check-format remains intentionally deferred until the pre-merge phase.
- Replaced the **explicit** double-shift polynomial QR step for a nonreal trailing 2x2 block with the genuine **implicit Francis real double-shift** algorithm. For the active upper Hessenberg H, compute only the first three entries of `p(H)e_1`, where `p(H)=H^2-tr(B)H+det(B)I` for the trailing 2x2 B. A three-entry Householder reflection begins a bulge, and subsequent 2/3-entry Householder reflections chase it along the subdiagonal. Each reflection applies a full orthogonal similarity to T and is accumulated in Q, preserving deflated trailing-block couplings without constructing H^2 or a full shift polynomial.
- Retained the Householder QR implementation for the single real Wilkinson-type shift. Left the finite-value checks, active-block deflation, final-real-2x2 reduction, 1x1/2x2 Schur structural requirements and original-input residual/orthogonality acceptance tests unchanged. The validated paradigm still refuses to misrepresent approximate factors as interval certificates. This is an initial implicit double shift, *not* a complete aggressive-early-deflation or multishift LAPACK implementation.
- Updated Doxygen on the Schur class and new helper with the implicit Q theorem, exact first-column formula, bulge formation/chasing argument, and Francis (1961-1962) and Wilkinson (1965) references. SchurDecomposition is now included explicitly in the public Algebra umbrella header.
- Added a DP/MP regression for a fully coupled 4x4 integer matrix orthogonally similar to real rotation-generator blocks of frequencies 1 and 3. Its independently known characteristic polynomial is `(lambda^2+1)(lambda^2+9)`, with spectrum `{-3i,-i,+i,+3i}`. Require two 2x2 conjugate Schur blocks, complete approximate Spectrum, and separate 1e-8 checks on A*Q-Q*T and Q^T Q-I. Existing 3x3 complex, symmetric, defective Jordan, dense and validated tests remain unchanged.
- Prototyped the bulge-chasing step in NumPy (not Ariadne C++), with the 3x3 and 4x4 exact-reference cases and 400 random matrices of orders 3, 4, 5 and 8: the prototype converged with small residuals. **C++ build and test results for this commit are pending; passing the prototype does not establish C++ correctness.** Check CI before proceeding.


### 2026-11-28 — fix GCC -Wcomment on Francis Doxygen equations (this commit)

- The new implicit Francis commit `7ff5c50` fails GCC 13 Coverage/Debug compilation with `-Werror=comment` in `include/algebra/schur_decomposition.hpp`: the Doxygen `split` LaTeX expression contains three lines ending in `\\`. C/C++ line splicing treats each trailing backslash as a multi-line comment continuation, even when intended as a TeX row separator.
- Replaced that aligned equation block with four separate Doxygen display-math expressions for x, y, z, and the trailing 2x2 trace/determinant. This keeps the complete mathematical derivation and references in the generated API documentation without any physical C++ source line ending in a backslash. The new mathematical formulas have exactly the same meanings.
- No C++ execution paths, API types, numerical tolerances or tests changed; no warning is disabled. Source inspection confirms no remaining backslash-terminated physical line in the header. The actual compiler and test results for the new commit are pending. Formatting remains deferred until pre-merge.


### 2026-11-28 — preserve SCC reducibility when applying Schur (this commit)

- Baseline `73100e5` passed **all 18 functional CI jobs** across Ubuntu/macOS/Windows, including Debug, Release, DP/MP, Coverage, ASan/UBSan and TSan. The new Francis implicit double shift and Doxygen GCC -Wcomment correction are functionally verified. Separate thresholds remain below target (Codecov project 92.69%, diff 80.06% versus 98.83%); check-format remains deferred to pre-merge.
- Found a structural inefficiency in `Spectrum<P,PR>`: discovering any strongly connected component (SCC) larger than 2 previously returned `_schur_spectrum(A)` on the **entire** input matrix, discarding the condensation graph and demanding an unnecessarily large numerical Schur solve. This could make one irreducible dense SCC determine success/failure for unrelated scalar and 2x2 SCCs, contradicting the intended blockwise API architecture.
- Retained all SCCs through the Warshall-style Boolean transitive closure, extracted each diagonal principal submatrix from the topological condensation order, and solved 1x1 structurally, 2x2 analytically, and **each order >=3 SCC independently** through approximate Schur. The output is their multiset union. A failure of any required block remains `INCONCLUSIVE`; validated dense blocks still cannot be certified by floating QR and stay `INCONCLUSIVE`. Added a top-level check that the output contains exactly n eigenvalues, counted with algebraic multiplicity.
- Expanded API-level Doxygen to derive the block upper triangular determinant factorization, explain exactly zero/interval graph edges and the independent SCC backend, and cite Warshall (1962), `A Theorem on Boolean Matrices`, JACM 9(1), 11-12, doi:10.1145/321105.321107; and Kressner (2005), `Numerical Methods for General and Structured Eigenvalue Problems`, Section 3.3.1, doi:10.1007/3-540-28502-4. The cited publications and relevant sections were checked.
- Added a mixed 6x6 DP/MP test with **interleaved** irreducible blocks of sizes 3, 2 and 1 and one-way cross-block couplings. Their independently known characteristic polynomial is `(lambda-9)(lambda^2+9)(lambda^2+1)(lambda-7)`; the test requires the complete multiset `{9, 7, -3i, 3i, -i, i}` without relying on the output order. Interval-family inputs are checked to remain conservatively inconclusive for the unverified large component. Existing dense 3x3/4x4 and validated structural test expectations are unchanged.
- **Verification:** mathematical block-factorization and source transformations inspected; the new GitHub Actions compilation/tests are pending and will be distinguished from the prior proven CI baseline. No benchmark run, formatting modification, test weakening, or hidden backward compatibility adapter was introduced.


### 2026-11-28 — repair macro-safe validated Spectrum assertions (this commit)

- The new mixed-SCC test in `bcad2f2` failed across Debug/Coverage/Sanitizers/Windows test builds, **before executing tests**. The supplied AppleClang diagnostic identifies two preprocessor errors at lines 379 and 386 of `test_spectrum.cpp`: `ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(context, expression)` misinterprets the comma in `Spectrum<ValidatedTag, PR>(precision)(v)` as an additional macro argument.
- Introduced one typed `const Spectrum<ValidatedTag, PR> validated(precision)` outside both assertions and call `validated(v)` inside the existing macro. This is a preprocessor syntax repair: neither expected `INCONCLUSIVE` status nor inputs, numeric thresholds, spectral algorithms or test coverage were changed.
- The Release builds in the previous run passed because the Debug test target was not compiled there; the failures do **not** establish an SCC numerical bug. The new commit's actual compiled/tested results are pending. Keep the `eigen` branch, Doxygen algorithm provenance, deferred check-format and compile-only CI benchmarks.

### 2026-11-28 — route real nonsymmetric diagonalization through Schur (this commit)

- Resolved first the `bcad2f2` Debug compiler regression in `2d9c9c3`: the legacy two-argument assertion macro parsed the comma inside `Spectrum<ValidatedTag,PR>` as a third argument. That repair only binds a local `validated` solver and leaves both expected validated INCONCLUSIVE outcomes unchanged. CI for the repair is still pending as this refactor is authored.
- Eliminated the duplicated shifted real Householder QR iteration in `EigenDecomposition::_nonsymmetric`. That implementation formerly managed its own active blocks, Wilkinson shifts, QR orthogonal factors and deflation before extracting its real Schur-like triangular form. It now calls the shared `SchurDecomposition<ApproximateTag,PR>` backend for Householder Hessenberg reduction and Francis/Wilkinson QR.
- The mathematical contract differs: `SchurDecomposition` may succeed with nonreal 2x2 conjugate blocks, while the current **real** `EigenDecomposition` succeeds only when the Schur block partition contains exactly n scalar 1x1 blocks. A 2x2 complex block therefore leaves the diagonalization INCONCLUSIVE (not a claim the spectrum is absent). The diagonal/back-substitution seed, transformation by Q, ORIGINAL-matrix Newton refinement, eigenvalue separation checks, original residual checks and eigenvector-basis inversion/identity check remain unchanged. The validated Gershgorin, verified similarity and Krawczyk paths are unaffected; approximate Schur is not a rigorous interval certificate.
- Updated EigenDecomposition's API-level Doxygen and inline method documentation to identify Schur as the sole dense QR backend and link to the Householder/Hessenberg/Francis/Wilkinson mathematical references. No backwards-compatibility wrappers and no benchmark execution; CI benchmarks remain compile-only.
- **Remaining campaign after this refactor:** (1) secure green functional CI for SCC change plus refactor and expand mixed real/complex/defective/clustered cases and DP/MP convergence; (2) complex eigenvector extraction / complete complex diagonalization, invariant-subspace representation for defective matrices, and explicit real vs complex contracts; (3) rigorous validated complex isolated-eigenvalue and cluster/invariant-subspace certification, distinguish undecidable cases; (4) numerical stability improvements for near-defective/non-normal and scaled matrices, iteration and conditioning diagnostics, performance profiling; (5) coverage threshold investigation, check-format/pre-merge formatting, installation/tutorial API verification, documentation and full final CI. This is a substantial outstanding campaign, not only a formatting exercise. Avoid declaring completion before these features and tests are implemented.
- **Verification:** code structure and mathematical API separation reviewed; C++ build and runtime tests of this commit pending. Test expectations, convergence tolerances and certification requirements were not relaxed.

### 2026-11-28 — first approximate complex diagonalization on real Schur factors (this commit)

- Reused the now centralized approximate `SchurDecomposition` backend: no second QR implementation was introduced. Its real quasi-triangular factors are traversed in 1x1 scalar and 2x2 nonreal-conjugate blocks. For a 2x2 block with eigenvalue lambda, an initial eigenvector is `(b,lambda-a)`; for each earlier diagonal Schur block, solve the triangular 1x1 or complex 2x2 system for the coupling to later blocks. Lift these Schur eigenvectors to the original real matrix using the orthogonal Q factor.
- Introduced an explicitly named `ComplexEigenDecomposition<P,PR>`, distinct from the existing real-only `EigenDecomposition<P,PR>`. Its successful data contain `vector<Complex<X>>` eigenvalues and paired `Matrix<Complex<X>>` eigenvectors / complex inverse. ApproximateTag accepts only numerically separated complex/real eigenvalues, finite intermediates, small residuals `A*V-V*D` against the ORIGINAL real matrix and a numerically invertible full complex basis. A custom pivoted complex Gauss-Jordan inversion checks BOTH `V*V^{-1}` and `V^{-1}*V`. Every inconclusive operation returns empty data. Repeated or defective cases remain inconclusive (this initial API is not a general repeated-semisimple eigensolver).
- ValidatedTag intentionally remains INCONCLUSIVE, since a small floating QR residual and complex inverse are not an interval proof. No attempt is made to assert that all matrices in a validated input interval family are uniformly diagonalizable.
- The public Doxygen explains the real Schur block back-substitution, the 2x2 complex eigenvector construction, lift by Q, complex pivoting and invertibility checks, and references Wilkinson (1965), Francis (1961-62) and SchurDecomposition's Householder/Hessenberg derivation. Public header is included in ariadne-algebra.hpp; the new test target is registered in tests/CMakeLists.txt.
- DP/MP regression cases include scalar and 2x2 rotation, fully strongly connected 3x3 with exact characteristic polynomial `(lambda-9)(lambda^2+9)`, the 4x4 two-pair spectrum `(lambda^2+1)(lambda^2+9)`, a dense simple-real case, the defective 3x3 Jordan refusal, validated refusal and zero-order handling. Independent test assertions verify original-matrix complex eigenpair residuals AND both inverse products, with 1e-8 absolute tolerances. Existing real decomposition tests remain unchanged.
- **Verification:** implementation and tests inspected, source blobs staged for GitHub commit, but C++ compilation and CI correctness are NOT yet established. Continue by examining all CI jobs and correcting compiler/algorithm defects without loosening numerical expectations. Existing format deferral and compile-only CI benchmark policy remain unchanged.

**Remaining campaign:** After the complex diagonalization first increment, (1) stabilize all functional jobs and check DP/MP precision and conditioning; (2) represent defective matrices and spectral invariant subspaces independently of diagonalizability; (3) extend ValidatedTag to isolated complex eigenvalues, clusters and certified invariant subspaces, including multiplicities; (4) improve convergence and scaling, add near-defective and non-normal regression tests, document numeric limitations; (5) final Doxygen/tutorials/coverage threshold analysis/format pre-merge/CI review, and only then merge. A complete campaign therefore remains materially larger than this first complex API.


### 2026-11-28 — repair complex matrix storage; converge on one eigendecomposition API (this commit)

- The first complex eigendecomposition build `1f531e6` failed when instantiating `Matrix<Complex<FloatApproximation<PR>>>`: `UniformArray` requires scalar `characteristics`, which Ariadne `Complex<X>` does not provide, and checked-assignment matrix element proxies are unsuitable as operands of several of `Complex<X>`'s arithmetic overloads. The supplied Clang logs show nested `Complex<Complex<X>>`, ambiguous multiplication and assignment errors under both DP and MP. Release builds did not instantiate the new Debug-only test executable.
- Introduced a minimal explicitly initialized `ComplexMatrix<X>` storage type with actual `Complex<X>&` indexed elements, replacing only the unsupported `Matrix<Complex<X>>` values in complex eigendecomposition results, the Schur eigenvector array, and pivoted Gauss-Jordan inversion. It preserves row/column-size inspection and `[][]` access required by the existing regression tests. No changes to numeric's generic `Matrix` or `Complex` are made, and no numerical thresholds, tests, or Schur solver algorithms are altered. This is a narrowly scoped compile-integration repair; the new CI result remains to be checked.
- **Final API decision:** the separate public `ComplexEigenDecomposition` class is a *temporary implementation boundary*, NOT a permanent choice that a user must make by guessing the matrix spectrum. The target public `EigenDecomposition` must automatically handle both real and complex eigenvalues for a real input matrix, with complex-valued eigenvalues/eigenvectors and a complete basis when diagonalizable. Users ask separately for `Spectrum` (all eigenvalues regardless of diagonalizability) or `SchurDecomposition` (real orthogonal quasi-triangular factorization). Existing real-only approximate/validated diagonalization code should be reused *internally*, not retained as a competing public class. There are no compatibility clients; no backwards-compatibility wrapper is required. A validated result must not claim that approximate complex Schur proves a complete interval eigenbasis.
- **Verification:** compile errors and source type requirements inspected; no C++ build of this change confirmed yet. Prior failing checks are not categorized as numerical failures. Once CI is green, unify the public API and add tests proving that the caller does not need to know the spectral field in advance. Keep check-format deferred to the premerge stage; benchmarks compile-only in CI.

### 2026-11-28 — repair GCC Root name shadow and capture the real fast-path contract (this commit)

- The latest `3f4ae85` job results show successful Clang, AppleClang and MSVC Debug tests, while GCC 13 Debug/Coverage fail in `test_complex_eigen_decomposition.cpp` compilation with `-Werror=shadow`. The local type `struct Root` in `ComplexEigenDecomposition::operator()` shadows `Ariadne::Root` from `numeric/operators.hpp:506`.
- Renamed that local record type to `SpectralRootEntry` and changed its vector declaration and range-for type. Verified the three original uses are renamed and no independent uses of the old type remain. Neither the eigenpair algorithm, convergence checks nor test expectations have changed, and no compiler warning is disabled.
- Promoted the **single public EigenDecomposition API and real fast-path** choice into the authoritative design section. A final unified result must allow real and complex representations without mandatory promotion or allocation of complex eigenvector matrices. Reuse all existing real specializations and verified real interval paths. Complex eigenvectors/inverse are computed only when needed. Before closing the campaign, benchmark DP/MP latency, allocations and iterations for the existing real-spectrum workload; Schur's different real QR can itself cause a regression even without complex arithmetic. Benchmarks remain compile-only in CI.
- **Verification:** GCC error and source type references inspected; new commit build/test and numerical performance comparisons remain pending. The fix is not described as CI-passing before the new checks complete. The temporary ComplexEigenDecomposition class remains for now solely while the common engine is being stabilized. Pre-merge formatting remains deferred.


### 2026-11-28 — unify public EigenDecomposition with zero-cost real-success dispatch (this commit)

- The preceding `916a837` build completed the Ubuntu, macOS, Coverage, and Sanitizers workflows successfully, including GCC Debug/Coverage and DP/MP complex eigenvector tests. Three Windows jobs were still running at the initial check; no claim is made about their later conclusion without another inspection. Check-format and Codecov thresholds remain separate issues.
- **Final public API architecture:** there is now one `EigenDecomposition<P,PR>` class whose successful payload is `std::variant<RealData,ComplexData>` in an optional result. The caller does not select a real versus complex eigensolver. `real_data()` and `complex_data()` inspect the actual representation. A numerical real success immediately returns the original real data with no allocation of complex matrices. This preserves the common real-case fast path, specialized symmetric Jacobi, structural and block triangular algorithms, and the validated real Gershgorin/Krawczyk/similarity certification routes.
- The existing real implementation, including its numerical/validated methods, was moved into `Ariadne::detail::RealEigenDecomposition` with its previous internal result/data contracts. The existing complex implementation was moved into `Ariadne::detail::ComplexEigenDecomposition`. Both are **implementation details**, not competing public solver choices; no backward compatibility wrappers are maintained. Existing semantic/regression tests are explicitly directed to the internal real engine to retain fine-grained coverage of its historical contracts and strategies, and the complex engine's existing direct tests remain internal.
- The public `EigenDecomposition` facade calls the real engine first. Only on an inconclusive approximate real result does it inspect a real Schur form for nonreal 2x2 conjugate blocks. If no such block is found, it returns inconclusive without complex allocation. Otherwise it uses the internal complex Schur eigenvector routine, which checks original-matrix `AV-VD`, complex basis invertibility, and two inverse identities. The new diagnostics strategy `COMPLEX_SCHUR` records the complex fallback. A failure never invents a basis or changes the meaning of validated real certification; ValidatedTag does not attempt unproved complex certification. The current fallback performs a second Schur reduction in the complex engine, which is an optimization opportunity for later.
- Added DP/MP **public-interface tests** for both real-only successful data (no complex representation, no COMPLEX_SCHUR strategy), conjugate complex rotation (same public class, complex data with original residual checks), defective Jordan rejection, validated diagonal-family real certificates, and validated complex-family refusal. All original real and complex algorithmic tests remain unchanged apart from addressing internal implementation names. The benchmark harness now recognizes COMPLEX_SCHUR and still exercises the public `EigenDecomposition` interface. Benchmarks remain compile-only in CI.
- Moved the public mathematical-contract explanation into the unified EigenDecomposition Doxygen, including links to Schur's Householder/Hessenberg/Francis/Wilkinson references. Removed the redundant umbrella inclusion of the now-internal complex implementation, which is instead an implementation dependency of the public unified header.
- **Verification:** source-level API transformations, type references, and test design inspected. This commit has not been C++ compiled or tested yet; monitor GCC, Clang, AppleClang, MSVC, sanitizers, and coverage before considering the unification stable. Do not weaken test expectations to make a build green, and keep formatting deferred to the pre-merge phase.

**Remaining campaign after public API unification:** stabilize all new CI jobs and repair genuine implementation issues; remove the duplicated real Schur work on the uncommon complex fallback; extend complex eigendecomposition to well-conditioned repeated semisimple cases, and distinguish certified non-diagonalizability from numerical inconclusiveness; represent invariant subspaces for defective and clustered matrices; implement validated complex isolated eigenvalue/cluster certificates; benchmark DP/MP throughput and allocations against the last known real path, especially symmetric/triangular/non-normal real workloads; add precision, conditioning and scaled/near-defective regressions; finish API tutorials, documentation, coverage and format review before a final premerge CI.


### 2026-11-28 — avoid duplicate Schur in the complex fallback (this commit)

- The unified public EigenDecomposition in `7ae432a` has passed Ubuntu, macOS, Coverage and Sanitizers, including DP/MP tests of both real and complex paths; Windows CI was still running at the time of the first check. Format (deferred) and external Codecov coverage thresholds remain separate.
- Removed one redundant real Schur QR computation on the **uncommon complex fallback**. Previously the public dispatcher computed Schur to detect a nonreal 2x2 diagonal block, then discarded its factors; the internal complex eigenvector engine recomputed Schur for back-substitution. The engine now accepts `const SchurDecompositionData<P,PR>&` as an overload, using the supplied Q, T and block partition directly while retaining the exact same finite checks, eigenvalue separation, original-input complex residual verification, complex basis pivoting and two inverse-identity tests. The convenience internal overload that computes Schur independently remains for focused algorithmic testing only.
- The public facade passes the just-computed Schur factors directly to the internal complex routine. On a **successful real path** the dispatcher returns before Schur detection or complex matrix allocation: no additional overhead or changed real semantics. An inconclusive real nonsymmetric path may already have performed a separate Schur attempt; avoiding that earlier duplication is future work and is not claimed here.
- Expanded Doxygen of the internal engine, its overload and public facade to explain factor reuse in terms of the Schur/Householder/Hessenberg/Francis QR chain, preserving attribution to the original methods. In the DP/MP complex regression, compare the original internal entry point against the explicit precomputed-Schur overload for every existing successful test matrix (scalar, rotation, mixed 3x3, two-pair 4x4 and real spectrum). Both must return SUCCESS with the same ordered eigenvalues; the original complex A*V=V*D and inverse-basis checks remain unchanged. No reduced tolerances and no test weakening.
- **Verification:** source transformations and Doxygen reviewed; the new commit's C++ CI compilation and runtime tests are pending. No changes to CI benchmark execution, format policy, validated certificates or the public API.

### 2026-11-28 — forward nonsymmetric real Schur factors to complex fallback (this commit)

- Baseline `de8856d`: all **18 functional CI jobs passed** (Ubuntu x86/ARM, macOS ARM, Windows x86/ARM, Debug/Release, Coverage and Sanitizers). Codecov coverage thresholds remain unmet, and check-format stays deferred to pre-merge. This verifies the prior reuse from nonreal detection into complex extraction.
- Removed the remaining redundant QR for the standard **dense nonsymmetric matrix with a nonreal pair**. The internal real method `_nonsymmetric` already computes a valid real Schur Q,T, then returns INCONCLUSIVE on its 2x2 conjugate block. Instead of dropping those factors, it now moves them to an optional output-only cache carried through `_attempt` and the internal real solver's three-argument overload. The public unified EigenDecomposition reuses that candidate for complex Schur-block back substitution and inverse checks. It does not recompute Hessenberg/Francis/Wilkinson QR.
- On **successful real paths** (symmetric Jacobi, structural or dense real), this cache is not populated: there is no complex matrix allocation, no Schur-factor copy and no complex arithmetic. A cache supplied by a caller is reset on entry. For other inconclusive paths that have not run nonsymmetric Schur (e.g. a 2x2 analytic real failure), the facade still computes Schur once to detect conjugate blocks and hands the factors to the complex solver. ValidatedTag remains on established certified real paths and never treats numerical Schur as a proof of a complex interval eigenspace.
- Added DP/MP regression requiring the internal real method to return INCONCLUSIVE **with Schur factors** on a strongly connected real 3x3 having exact spectrum `{9,-3i,3i}`. Checked its 1x1/2x2 block partition, successful complex extraction from the saved factors, and public unified success with strategy COMPLEX_SCHUR. A following successful symmetric real call must clear any previous cached factors and keep the public data real. Existing numerical tolerances, residual and inverse checks, and test scenarios are unchanged.
- Expanded API Doxygen with the move-only-on-failure ownership handoff and its relation to the mathematical Householder/Hessenberg/Francis Schur solver. No additional public solver class, backward-compatibility wrapper, format run or CI benchmark execution. **Verification:** source reviewed and GitHub CI of this commit pending; do not claim test success until the new jobs complete.


### 2026-11-28 — initial leading invariant subspaces from real Schur (this commit)

- Baseline `7ccbef4`: at the time of this increment, Ubuntu, macOS and Coverage functional workflows have passed; Windows x86/ARM and one sanitizer job were still running and are **not** claimed green. Codecov coverage and check-format remain separate outstanding checks.
- Introduced `SchurInvariantSubspaceData<P,PR>` with two real matrices: an n-by-k orthonormal basis `U` and a k-by-k restriction `B` for the approximate invariant-subspace relation `A*U ≈ U*B`. The `SchurDecompositionData::leading_invariant_subspace(block_count, precision)` method extracts only the first `block_count` **complete** diagonal Schur blocks from the existing real Schur factors: `U=Q(:,0:k)`, `B=T(0:k,0:k)`. It never reruns QR or allocates complex matrices. This operation does not require any individual eigenvectors, distinct eigenvalues, or diagonalizability.
- Mathematical basis: if `A*Q≈Q*T` and T is block upper triangular with a leading block B and exact zero lower-left coupling, the first k columns of Q satisfy `A*Q1≈Q1*B`. A real 2x2 Schur block associated with a nonreal conjugate pair is **indivisible**, representing a two-dimensional REAL invariant subspace. A defective Jordan 3x3 also supports nested leading invariant subspaces even though it has no full eigenbasis. Doxygen states the real Schur theorem, formula, and links to Householder/Hessenberg/Francis/Wilkinson origins (Wilkinson 1965, Chapters 6 and 8).
- The extractor validates n-by-n Q/T dimensions, a nonempty and complete 1x1/2x2 block partition, requested prefix count and zero lower-left T coupling at the chosen boundary. Malformed metadata and attempts to split a coupled 2x2 pair yield no result. It does **not** claim to select arbitrary eigenvalue clusters; Schur block reordering is a separate future operation. ValidatedTag remains unsupported, because approximate factors cannot establish rigorous interval invariant subspaces.
- Added DP/MP regressions on a defective 3x3 Jordan block (each of three prefix dimensions), a 2x2 pure rotation (only the full real invariant plane), and a fully coupled 3x3 with spectrum `{9,-3i,+3i}`. Tests check `A*U-U*B`, `U^T U-I`, and `U^T A U-B` against original inputs with 1e-8 residual limits, reject zero/excessive block_count, malformed block sizes, and an attempted split of a nonreal Schur pair. Existing spectrum, eigenbasis and validated tests are untouched.
- **Verification:** source-level mathematics and bounds inspected. New C++ compilation and runtime CI of this increment are pending; do not claim a pass until the checks finish. Benchmarks are compile-only in CI, and check-format remains deferred to premerge. Once CI is green, add Schur block reordering for arbitrary isolated clusters, investigate defective/near-defective stability and design rigorous interval invariant-subspace certificates.
