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
