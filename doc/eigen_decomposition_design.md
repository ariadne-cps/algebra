# Eigenvalue and eigenvector computation — design journal

## Scope and workflow

Development branch: `eigen` in `ariadne-cps/algebra`, branched from `main` at `1b85f7e864044b697be2408c3575181a75f36d4b`.
Do not merge into `main` until implementation, numerical and semantic tests, and review are complete.

**Journal rule:** Every development commit must include an update to this document describing its design decisions, code changes, tests run, test results and outstanding issues. The first commit creates this journal. No successful build or test is to be claimed without running it.

## Purpose

Implement reusable approximate and validated spectral methods in Algebra. Evolution's first-order PDE solver currently takes user-supplied characteristic eigenvalues/eigenvectors. Algebra must remain independent of PDE assumptions; integration with Evolution is a later project.

## Public API direction (proposed; not yet implemented)

- `EigenDecomposition<P, PR>`: algorithm/solver object, with `P` in `ApproximateTag` or `ValidatedTag` and `PR` denoting numeric precision.
- `EigenDecompositionResult<P, PR>`: corresponding result type. `ApproximateTag` returns approximations without certification. `ValidatedTag` provides mathematically guaranteed enclosures **only when certified**, and an explicit inconclusive status otherwise.
- Support real and complex eigenvalues eventually; do not suggest every real matrix has a real eigenbasis.
- A failed certificate must never be reported as proof that no decomposition exists.
- Prefer an optional validated result over a success-looking empty or invalid enclosure.
- Input matrices should use the existing Ariadne number/matrix conventions. The exact call signatures, storage of complex outputs, convergence controls and precision traits remain to be decided after checking available primitives.

## Mathematical contracts

For a real square matrix `M`, a diagonalization `M T = T D` requires `T` invertible. Some matrices are defective; some have nonreal eigenvalues. A general eigensolver should ideally expose a Schur-based result even if diagonalization fails.

For the first increment, target **simple real eigenvalues and their normalized eigenvectors**. This is a deliberate limited-scope milestone, not a general eigenvalue decomposition.

Approximate stage:
1. Produce candidate eigenpairs using a numerically stable standard method (evaluate existing Algebra capabilities, then decide whether to use an external backend or implement a specific algorithm).
2. Measure residuals and account for convergence and ill-conditioning, without treating small numerical residuals as proof.

Validated stage:
1. Start from an approximate candidate.
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
