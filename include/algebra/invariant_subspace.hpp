/***************************************************************************
 *            algebra/invariant_subspace.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Algebra.
 *
 *  Ariadne Algebra is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Algebra is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/invariant_subspace.hpp
 *  @brief Exact structural certification of coordinate invariant subspaces.
 */
#ifndef ARIADNE_INVARIANT_SUBSPACE_HPP
#define ARIADNE_INVARIANT_SUBSPACE_HPP

#include <exception>
#include <optional>
#include <utility>
#include <vector>

#include "algebra/matrix.hpp"
#include "algebra/schur_decomposition.hpp"
#include "algebra/validated_root_contractors.hpp"

namespace Ariadne {

//! \brief A certified coordinate subspace for a family of real matrices.
//! \details Every member A of the input interval matrix family preserves
//! the REAL subspace spanned by the indicated coordinate basis columns.
//! The n-by-k matrix U has EXACT coordinate columns (zeros and ones),
//! so U^T U=I and A*U=U*B hold for every member A, where B is its
//! corresponding k-by-k restriction. The interval matrix stored in
//! restriction encloses all such B; the same quantified A occurs on
//! both sides of the equality. This is a structural, not numerical,
//! certificate and does not assert anything about individual eigenpairs.
template<class PR>
struct CertifiedCoordinateInvariantSubspaceData {
    using X = FloatBounds<PR>;
    Matrix<X> basis;             //!< Exact n-by-k coordinate basis.
    Matrix<X> restriction;       //!< Interval family of k-by-k restrictions.
    std::vector<SizeType> coordinates; //!< Basis-column coordinate indices.
};

//! \brief Certify a coordinate invariant subspace for an interval family.
//! \param A Square interval matrix describing independent entry bounds.
//! \param coordinates Nonempty distinct coordinate indices, in the
//! requested order. Ordering changes the coordinate basis, not its span.
//! \param precision Precision for EXACT zero and one basis entries.
//! \return Certified basis and interval restriction if the whole family
//! preserves the subspace; no value otherwise.
//!
//! \details For a set of coordinates S, let U consist of the columns
//! e_j of the identity for j in S, in the requested order. Then
//! \f[ A U = U B \quad\Longleftrightarrow\quad
//!     A_{S^c,S}=0,\qquad B=A_{S,S}. \f]
//! The equality is algebraic and requires no eigenvalue separation,
//! no diagonalizability, no Schur computation and NO floating tolerance.
//! It holds uniformly for EVERY point matrix in the input interval
//! family exactly when all cross-block coefficients A(i,j), with
//! i outside S and j inside S, are the singleton interval [0,0].
//! If a cross-block interval merely CONTAINS zero, invariance for
//! the whole family is NOT proven and the function returns no result.
//!
//! This works for nonconsecutive coordinate sets, repeated eigenvalues,
//! defective Jordan matrices and blocks with nonreal conjugate spectra.
//! A full coordinate set always spans the full real space, provided
//! the supplied matrix is finite. An empty set is deliberately not
//! returned. Finite bounds are required throughout the input.
//!
//! This is an elementary block-invariance theorem (see the reducible
//! block triangular forms in Wilkinson, The Algebraic Eigenvalue
//! Problem, 1965, and Stewart & Sun, Matrix Perturbation Theory,
//! 1990, for the broader invariant-subspace context). It is distinct
//! from the approximate Schur invariant subspaces, which cannot be
//! promoted to validated certificates by checking small residuals.
//! Arbitrary NON-coordinate invariant subspace certificates and
//! rigorous spectral cluster projectors remain future work.
template<class PR>
std::optional<CertifiedCoordinateInvariantSubspaceData<PR>>
certified_coordinate_invariant_subspace(
        const Matrix<FloatBounds<PR>>& A,
        const std::vector<SizeType>& coordinates, PR precision) {
    using X = FloatBounds<PR>;
    ARIADNE_PRECONDITION(A.row_size() == A.column_size());
    const SizeType n = A.row_size();
    const SizeType k = coordinates.size();
    if (n == 0u || k == 0u || k > n) {
        return std::nullopt;
    }

    const X exact_zero(0u, precision);
    const X exact_one(1u, precision);
    std::vector<bool> included(n, false);
    for (const SizeType j: coordinates) {
        if (j >= n || included[j]) { return std::nullopt; }
        included[j] = true;
    }

    // A nonsingleton interval containing zero is NOT an exact zero
    // coefficient. All cross-block coefficients must be identically
    // zero for every matrix in the interval family.
    for (SizeType row = 0u; row != n; ++row) {
        for (SizeType col = 0u; col != n; ++col) {
            const X& coefficient = A[row][col];
            if (!is_finite(coefficient.lower_raw()) ||
                !is_finite(coefficient.upper_raw()) ||
                coefficient.lower_raw() > coefficient.upper_raw()) {
                return std::nullopt;
            }
            if (!included[row] && included[col] &&
                (coefficient.lower_raw() != exact_zero.lower_raw() ||
                 coefficient.upper_raw() != exact_zero.upper_raw())) {
                return std::nullopt;
            }
        }
    }

    Matrix<X> basis(n, k, exact_zero);
    Matrix<X> restriction(k, k, exact_zero);
    for (SizeType col = 0u; col != k; ++col) {
        basis[coordinates[col]][col] = exact_one;
        for (SizeType row = 0u; row != k; ++row) {
            restriction[row][col] =
                A[coordinates[row]][coordinates[col]];
        }
    }
    return CertifiedCoordinateInvariantSubspaceData<PR>{
        std::move(basis), std::move(restriction), coordinates};
}

//! \brief Rigorous invariant graph subspace, with exact fixed basis.
//! \details U=[I;X] is an n-by-k FULL-COLUMN-RANK point basis, because
//! its leading k-by-k block is the identity. X consists of singleton
//! intervals whose values are the same for every input matrix.
//! For each A in the interval input family, A*U=U*B(A) holds exactly,
//! where B(A) is the restriction in these fixed graph coordinates.
//! The stored interval restriction encloses all these B(A); its
//! elements must not be treated as independently chosen witnesses.
template<class PR>
struct CertifiedInvariantGraphData {
    using X = FloatBounds<PR>;
    Matrix<X> basis;             //!< Exact fixed graph columns [I;X].
    Matrix<X> restriction;       //!< Enclosure for B(A) for each A.
    Matrix<X> graph;             //!< Fixed singleton (n-k)-by-k graph.
};

//! \brief Certify that a prescribed non-coordinate graph is invariant.
//! \param A Finite n-by-n real interval matrix representing a family.
//! \param graph Finite singleton (n-k)-by-k matrix, with 0<k<n.
//! \param precision Precision used to build the exact identity basis.
//! \return A rigorous, full-column-rank graph subspace and interval
//! restriction, or no value if the interval proof is inconclusive.
//!
//! \details Split the coordinates into the first k and the last n-k:
//! \f[ A=\begin{pmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{pmatrix},
//! \qquad U=\begin{pmatrix}I_k\\X\end{pmatrix}. \f]
//! For a FIXED real point matrix X, define
//! \f[ B(A)=A_{11}+A_{12}X. \f]
//! The graph is invariant for each A precisely when
//! \f[ A_{21}+A_{22}X-X(A_{11}+A_{12}X)=0. \f]
//! Evaluate this expression with directed-rounding interval arithmetic
//! on the ORIGINAL input family, using no midpoint or floating Schur
//! residual. SUCCESS requires the interval enclosure of EVERY entry
//! to be the exact singleton [0,0]. That condition implies the
//! equation for every represented A, including all dependent uses
//! of the same coefficient within the expression. The converse need
//! not be detected: natural interval evaluation loses correlations,
//! so valid graphs may be reported INCONCLUSIVE.
//!
//! The fixed graph is supplied as singleton intervals, not an interval
//! family of candidate bases. If graph entries are uncertain, returning
//! an invariant subspace for ALL matrices would require a quantified
//! existence proof and a coherent X(A); this routine never makes
//! that unsupported assertion. U has full column rank independently
//! of the spectrum because it contains I_k. B(A) may vary with A.
//!
//! In contrast to coordinate-subspace certification, this method can
//! establish non-axis-aligned invariant lines and higher-dimensional
//! graphs without requiring real eigenvalue simplicity, complex
//! eigenvectors or diagonalizability. It is nevertheless a SUFFICIENT,
//! conservative proof for a caller-provided FIXED graph, not a
//! validated algorithm for discovering unknown invariant graphs or
//! spectral-cluster enclosures.
//!
//! \par Mathematical provenance
//! Graph subspaces and the matrix Riccati invariance equation are
//! standard in invariant-subspace perturbation theory; see
//! G. W. Stewart and J.-G. Sun, Matrix Perturbation Theory,
//! Academic Press (1990), chapters on invariant subspaces; see
//! also J. H. Wilkinson, The Algebraic Eigenvalue Problem (1965).
//! The rigorous certificate here is the exact graph-block identity
//! combined with outward interval evaluation of its residual.
template<class PR>
std::optional<CertifiedInvariantGraphData<PR>>
certified_invariant_graph(const Matrix<FloatBounds<PR>>& A,
                          const Matrix<FloatBounds<PR>>& graph,
                          PR precision) {
    using X = FloatBounds<PR>;
    ARIADNE_PRECONDITION(A.row_size() == A.column_size());
    const SizeType n = A.row_size();
    const SizeType k = graph.column_size();
    if (n < 2u || k == 0u || k >= n ||
        graph.row_size() != n - k) {
        return std::nullopt;
    }
    const SizeType m = n - k;
    const X exact_zero(0u, precision);
    const X exact_one(1u, precision);

    const auto valid_bounds = [](const X& value) {
        return is_finite(value.lower_raw()) &&
               is_finite(value.upper_raw()) &&
               value.lower_raw() <= value.upper_raw();
    };
    for (SizeType row = 0u; row != n; ++row) {
        for (SizeType col = 0u; col != n; ++col) {
            if (!valid_bounds(A[row][col])) {
                return std::nullopt;
            }
        }
    }
    for (SizeType row = 0u; row != m; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            const X& value = graph[row][col];
            if (!valid_bounds(value) ||
                value.lower_raw() != value.upper_raw()) {
                return std::nullopt;
            }
        }
    }

    // B(A)=A11+A12*X. The interval matrix encloses restrictions
    // coherently induced by each point matrix A in the input family.
    Matrix<X> restriction(k, k, exact_zero);
    for (SizeType row = 0u; row != k; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            X entry(A[row][col]);
            for (SizeType j = 0u; j != m; ++j) {
                entry += A[row][k + j] * graph[j][col];
            }
            if (!valid_bounds(entry)) { return std::nullopt; }
            restriction[row][col] = entry;
        }
    }

    // An EXACT singleton-zero interval extension is a rigorous
    // sufficient condition for every represented A. An interval
    // merely containing zero cannot be used as a proof.
    for (SizeType row = 0u; row != m; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            X residual(A[k + row][col]);
            for (SizeType j = 0u; j != m; ++j) {
                residual += A[k + row][k + j] * graph[j][col];
            }
            for (SizeType j = 0u; j != k; ++j) {
                residual -= graph[row][j] * restriction[j][col];
            }
            if (!valid_bounds(residual) ||
                residual.lower_raw() != exact_zero.lower_raw() ||
                residual.upper_raw() != exact_zero.upper_raw()) {
                return std::nullopt;
            }
        }
    }

    Matrix<X> basis(n, k, exact_zero);
    for (SizeType col = 0u; col != k; ++col) {
        basis[col][col] = exact_one;
    }
    for (SizeType row = 0u; row != m; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            basis[k + row][col] = graph[row][col];
        }
    }
    return CertifiedInvariantGraphData<PR>{
        std::move(basis), std::move(restriction), graph};
}

//! \brief A uniformly certified family of invariant graph subspaces.
//! \details For EVERY point matrix A from the input interval family,
//! exactly one graph X(A) exists in graph_enclosure within the
//! specified Krawczyk box, and A*[I;X(A)]=[I;X(A)]*B(A).
//! Unlike CertifiedInvariantGraphData, the basis here may depend on A.
//! The displayed interval matrices are correlated enclosures, NOT
//! independent choices of U, B and A that satisfy literal interval
//! matrix equality.
template<class PR>
struct CertifiedInvariantGraphEnclosureData {
    using X = FloatBounds<PR>;
    Matrix<X> graph_enclosure;       //!< Unique local graph X(A), n-k by k.
    Matrix<X> basis_enclosure;       //!< [I;X(A)] in each family member.
    Matrix<X> restriction_enclosure; //!< A11+A12*X(A), enclosed for all A.
};

//! \brief Prove local existence of an invariant graph for each interval A.
//! \param A Original square interval family, with independent entries.
//! \param graph_center Singleton (n-k)-by-k approximate graph matrix.
//! \param graph_box Finite interval search box, containing every center
//! entry in its STRICT interior; each graph entry is a separate unknown.
//! \param preconditioner k(n-k)-by-k(n-k) interval enclosure of a
//! fixed numerical approximate inverse of the graph Riccati Jacobian.
//! \param precision Precision for verified interval arithmetic.
//! \return Uniformly verified graph enclosures, or inconclusive.
//!
//! \details The unknown X is a real (n-k)-by-k matrix, vectorized in
//! row-major order, and
//! \f[ F_A(X)=A_{21}+A_{22}X-X(A_{11}+A_{12}X). \f]
//! For a variable X_ab, the (i,j) Jacobian component is
//! \f[ \frac{\partial F_{ij}}{\partial X_{ab}}=
//!   \delta_{jb}(A_{22,ia}-(XA_{12})_{ia})
//!   -\delta_{ia}(A_{11}+A_{12}X)_{bj}. \f]
//! These formulas generate interval enclosures of F_A(X0), and of
//! D_X F_A throughout the whole graph_box, for EVERY input A.
//! The existing krawczyk_certifies_unique_root operation then checks
//! (a) strict enclosure of the Krawczyk image inside the graph box
//! and (b) strict infinity-norm contraction ||I-C*J|| < 1.
//! By the parameter-uniform interval fixed-point argument, every
//! represented A admits ONE real graph X(A) in the specified box,
//! and that graph is unique IN THIS BOX. No spectral simplicity,
//! eigenvector basis or non-defectiveness is assumed.
//!
//! The graph is NOT asserted to be the same for different A, and
//! neither this proof nor the point-graph structural certificate
//! establishes spectral-cluster separation or eigenvalue multiplicity.
//! All failures are inconclusive (never a nonexistence proof).
//! Input graph boxes are caller-provided: this routine CERTIFIES
//! candidates but does not yet discover or optimize a midpoint graph.
//! This preserves the conceptual separation of candidate generation
//! and rigorous interval proof, as in the existing validated eigenpair
//! Krawczyk pipeline.
//!
//! \par Mathematical provenance
//! The matrix Riccati graph equation is standard in invariant
//! subspace perturbation theory (G. W. Stewart and J.-G. Sun,
//! Matrix Perturbation Theory, Academic Press, 1990).
//! The validated fixed-point step builds on R. Krawczyk,
//! "Newton-Algorithmen zur Bestimmung von Nullstellen mit
//! Fehlerschranken", Computing 4, 187-201 (1969).
//! The shared contraction-and-inclusion theorem is documented
//! on krawczyk_certifies_unique_root.
template<class PR>
std::optional<CertifiedInvariantGraphEnclosureData<PR>>
certified_invariant_graph_near(
        const Matrix<FloatBounds<PR>>& A,
        const Matrix<FloatBounds<PR>>& graph_center,
        const Matrix<FloatBounds<PR>>& graph_box,
        const Matrix<FloatBounds<PR>>& preconditioner,
        PR precision) {
    using X = FloatBounds<PR>;
    ARIADNE_PRECONDITION(A.row_size() == A.column_size());
    const SizeType n = A.row_size();
    const SizeType k = graph_center.column_size();
    if (n < 2u || k == 0u || k >= n) {
        return std::nullopt;
    }
    const SizeType m = n - k;
    const SizeType count = m * k;
    if (graph_center.row_size() != m ||
        graph_box.row_size() != m ||
        graph_box.column_size() != k ||
        preconditioner.row_size() != count ||
        preconditioner.column_size() != count) {
        return std::nullopt;
    }
    const X zero(0u, precision);
    const X one(1u, precision);
    const auto finite = [](const X& x) {
        return is_finite(x.lower_raw()) &&
            is_finite(x.upper_raw()) &&
            x.lower_raw() <= x.upper_raw();
    };
    for (SizeType i = 0u; i != n; ++i) {
        for (SizeType j = 0u; j != n; ++j) {
            if (!finite(A[i][j])) { return std::nullopt; }
        }
    }
    for (SizeType i = 0u; i != count; ++i) {
        for (SizeType j = 0u; j != count; ++j) {
            if (!finite(preconditioner[i][j])) { return std::nullopt; }
        }
    }

    Vector<X> center(count, zero);
    Vector<X> box(count, zero);
    for (SizeType i = 0u; i != m; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            const SizeType index = i * k + j;
            const X& at = graph_center[i][j];
            const X& domain = graph_box[i][j];
            if (!finite(at) || !finite(domain) ||
                at.lower_raw() != at.upper_raw() ||
                !(domain.lower_raw() < at.lower_raw() &&
                  at.upper_raw() < domain.upper_raw())) {
                return std::nullopt;
            }
            center[index] = at;
            box[index] = domain;
        }
    }

    // B0=A11+A12*X0 and Bbox=A11+A12*Xbox. Both are computed
    // from the ORIGINAL interval family (no midpoint substitution).
    Matrix<X> base_center(k, k, zero);
    Matrix<X> base_box(k, k, zero);
    for (SizeType i = 0u; i != k; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            X c(A[i][j]);
            X b(A[i][j]);
            for (SizeType a = 0u; a != m; ++a) {
                c += A[i][k + a] * graph_center[a][j];
                b += A[i][k + a] * graph_box[a][j];
            }
            if (!finite(c) || !finite(b)) { return std::nullopt; }
            base_center[i][j] = c;
            base_box[i][j] = b;
        }
    }

    Vector<X> residual(count, zero);
    Matrix<X> jacobian(count, count, zero);
    for (SizeType i = 0u; i != m; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            const SizeType row = i * k + j;
            X value(A[k + i][j]);
            for (SizeType a = 0u; a != m; ++a) {
                value += A[k + i][k + a] * graph_center[a][j];
            }
            for (SizeType b = 0u; b != k; ++b) {
                value -= graph_center[i][b] * base_center[b][j];
            }
            if (!finite(value)) { return std::nullopt; }
            residual[row] = value;

            // dF_ij / dX_ab =
            // delta_jb*(A22_ia - (X*A12)_ia) - delta_ia*B_bj.
            for (SizeType a = 0u; a != m; ++a) {
                for (SizeType b = 0u; b != k; ++b) {
                    X derivative(zero);
                    if (j == b) {
                        derivative += A[k + i][k + a];
                        for (SizeType c = 0u; c != k; ++c) {
                            derivative -= graph_box[i][c] * A[c][k + a];
                        }
                    }
                    if (i == a) {
                        derivative -= base_box[b][j];
                    }
                    if (!finite(derivative)) { return std::nullopt; }
                    jacobian[row][a * k + b] = derivative;
                }
            }
        }
    }

    // THIS is the proof gate. A small midpoint residual alone
    // is never sufficient to create a validated result.
    if (!krawczyk_certifies_unique_root(
            center, residual, jacobian, box, preconditioner)) {
        return std::nullopt;
    }

    Matrix<X> basis(n, k, zero);
    for (SizeType i = 0u; i != k; ++i) {
        basis[i][i] = one;
    }
    for (SizeType i = 0u; i != m; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            basis[k + i][j] = graph_box[i][j];
        }
    }
    return CertifiedInvariantGraphEnclosureData<PR>{
        graph_box, std::move(basis), std::move(base_box)};
}

//! \brief Approximate a Schur-selected invariant subspace as a graph.
//! \param A Real point matrix used for the real Schur decomposition.
//! \param selected_blocks Complete diagonal Schur blocks, indexed in
//! their original real Schur ordering.
//! \param precision Working double or multiple precision.
//! \return Approximate (n-k)-by-k graph X, or no value if selection
//! is unstable, the leading chart is singular, or the residual fails.
//!
//! \details The real Schur solver produces an orthonormal basis
//! Q1 (n-by-k) of the selected invariant subspace after orthogonal
//! block reordering. Partition Q1=[R;S] into its top k and bottom
//! n-k rows. If R is nonsingular then the SAME subspace is the graph
//! of X=S*R^{-1}, because Q1*R^{-1}=[I;X]. The numerical inverse
//! and original-matrix Riccati residual are both checked before
//! returning the candidate; a nearly vertical chart is INCONCLUSIVE.
//!
//! This is candidate generation ONLY, in real approximate arithmetic,
//! even if the selected Schur block represents a conjugate pair.
//! It is not a validated invariant-subspace enclosure. A validated
//! proof can be attempted with certified_invariant_graph_near,
//! using a caller-supplied box and preconditioner. Schur block swaps
//! follow Bai-Demmel (LAWN 54 / LAPACK DTREXC); see also Stewart and
//! Sun, Matrix Perturbation Theory (1990), on invariant graph charts.
template<class PR>
std::optional<Matrix<FloatApproximation<PR>>>
approximate_invariant_graph_from_schur(
        const Matrix<FloatApproximation<PR>>& A,
        const std::vector<SizeType>& selected_blocks,
        PR precision) {
    using AX = FloatApproximation<PR>;
    using Raw = typename AX::RawType;
    ARIADNE_PRECONDITION(A.row_size() == A.column_size());
    const SizeType n = A.row_size();
    if (n < 2u || selected_blocks.empty()) { return std::nullopt; }
    const AX zero(0u, precision);
    const AX one(1u, precision);
    const auto finite = [](const AX& x) { return is_finite(x.raw()); };
    AX scale(one);
    for (SizeType i = 0u; i != n; ++i) {
        for (SizeType j = 0u; j != n; ++j) {
            if (!finite(A[i][j])) { return std::nullopt; }
            if (abs(A[i][j]).raw() > scale.raw()) {
                scale = abs(A[i][j]);
            }
        }
    }
    const auto schur = SchurDecomposition<ApproximateTag, PR>(precision)(A);
    if (schur.status != SchurStatus::SUCCESS || !schur.data) {
        return std::nullopt;
    }
    const auto subspace = schur.data->selected_invariant_subspace(
        selected_blocks, precision);
    if (!subspace) { return std::nullopt; }
    const SizeType k = subspace->basis.column_size();
    if (k == 0u || k >= n || subspace->basis.row_size() != n) {
        return std::nullopt;
    }
    const SizeType m = n - k;
    Matrix<AX> leading(k, k, zero);
    Matrix<AX> trailing(m, k, zero);
    for (SizeType i = 0u; i != k; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            leading[i][j] = subspace->basis[i][j];
        }
    }
    for (SizeType i = 0u; i != m; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            trailing[i][j] = subspace->basis[k + i][j];
        }
    }

    Matrix<AX> inverse_chart(k, k, zero);
    try {
        inverse_chart = inverse(leading);
    } catch (const std::exception&) {
        return std::nullopt;
    }
    const Matrix<AX> identity_check = leading * inverse_chart;
    const AX eps(Raw::eps(precision));
    const AX chart_limit =
        AX(4096u, precision) * eps *
        AX(static_cast<unsigned>(n * n), precision);
    for (SizeType i = 0u; i != k; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            if (!finite(inverse_chart[i][j]) ||
                !finite(identity_check[i][j]) ||
                abs(identity_check[i][j] -
                    (i == j ? one : zero)).raw() > chart_limit.raw()) {
                return std::nullopt;
            }
        }
    }
    Matrix<AX> graph = trailing * inverse_chart;
    AX graph_scale(one);
    for (SizeType i = 0u; i != m; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            if (!finite(graph[i][j])) { return std::nullopt; }
            if (abs(graph[i][j]).raw() > graph_scale.raw()) {
                graph_scale = abs(graph[i][j]);
            }
        }
    }
    // Reject a numerically vertical chart even if the matrix
    // inverse happened not to raise an exception.
    const AX chart_scale_limit = one / (AX(128u, precision) * eps);
    if (!finite(chart_scale_limit) ||
        graph_scale.raw() > chart_scale_limit.raw()) {
        return std::nullopt;
    }

    Matrix<AX> restriction(k, k, zero);
    for (SizeType i = 0u; i != k; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            AX value(A[i][j]);
            for (SizeType a = 0u; a != m; ++a) {
                value += A[i][k + a] * graph[a][j];
            }
            if (!finite(value)) { return std::nullopt; }
            restriction[i][j] = value;
        }
    }
    const AX residual_limit =
        AX(8192u, precision) * eps * scale * graph_scale *
        AX(static_cast<unsigned>(n * n), precision);
    if (!finite(residual_limit)) { return std::nullopt; }
    for (SizeType i = 0u; i != m; ++i) {
        for (SizeType j = 0u; j != k; ++j) {
            AX residual(A[k + i][j]);
            for (SizeType a = 0u; a != m; ++a) {
                residual += A[k + i][k + a] * graph[a][j];
            }
            for (SizeType b = 0u; b != k; ++b) {
                residual -= graph[i][b] * restriction[b][j];
            }
            if (!finite(residual) ||
                abs(residual).raw() > residual_limit.raw()) {
                return std::nullopt;
            }
        }
    }
    return graph;
}

//! \brief Generate a real Schur graph candidate and certify its enclosure.
//! \param A ORIGINAL interval matrix family.
//! \param selected_blocks Whole blocks of the midpoint real Schur form.
//! \param graph_box Supplied interval box for the selected invariant
//! graph; it must strictly contain the generated singleton center.
//! \param preconditioner Supplied interval preconditioner for the
//! row-major graph Riccati Jacobian.
//! \param precision Precision for approximate and verified arithmetic.
//!
//! \details Convert the interval matrix midpoint to a real approximate
//! matrix, compute and reorder real Schur, form the graph chart X0,
//! then call certified_invariant_graph_near on the ORIGINAL interval
//! matrix A. A small approximate residual is NEVER a certificate:
//! SUCCESS requires the same strict Krawczyk inclusion and contraction
//! for every represented A as the underlying validated verifier.
//! The box and preconditioner remain caller-supplied; this function
//! does not run an iterative equation solver or depend on the
//! equation-solving repository. Failure is inconclusive.
//! Mathematical sources: real Schur reordering (Bai-Demmel, 1993),
//! invariant graph/Riccati equations (Stewart-Sun, 1990), and
//! Krawczyk's validated inclusion theorem (1969).
template<class PR>
std::optional<CertifiedInvariantGraphEnclosureData<PR>>
certified_invariant_graph_from_schur(
        const Matrix<FloatBounds<PR>>& A,
        const std::vector<SizeType>& selected_blocks,
        const Matrix<FloatBounds<PR>>& graph_box,
        const Matrix<FloatBounds<PR>>& preconditioner,
        PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    ARIADNE_PRECONDITION(A.row_size() == A.column_size());
    const SizeType n = A.row_size();
    if (n < 2u || selected_blocks.empty()) { return std::nullopt; }
    const AX zero(0u, precision);
    Matrix<AX> midpoint(n, n, zero);
    for (SizeType i = 0u; i != n; ++i) {
        for (SizeType j = 0u; j != n; ++j) {
            if (!is_finite(A[i][j].lower_raw()) ||
                !is_finite(A[i][j].upper_raw()) ||
                A[i][j].lower_raw() > A[i][j].upper_raw()) {
                return std::nullopt;
            }
            midpoint[i][j] = AX(A[i][j]);
            if (!is_finite(midpoint[i][j].raw())) {
                return std::nullopt;
            }
        }
    }
    const auto candidate = approximate_invariant_graph_from_schur(
        midpoint, selected_blocks, precision);
    if (!candidate || graph_box.row_size() != candidate->row_size() ||
        graph_box.column_size() != candidate->column_size()) {
        return std::nullopt;
    }
    const VX verified_zero(0u, precision);
    Matrix<VX> center(candidate->row_size(), candidate->column_size(),
                      verified_zero);
    for (SizeType i = 0u; i != center.row_size(); ++i) {
        for (SizeType j = 0u; j != center.column_size(); ++j) {
            center[i][j] = VX((*candidate)[i][j].raw());
        }
    }
    return certified_invariant_graph_near(
        A, center, graph_box, preconditioner, precision);
}

} // namespace Ariadne

#endif // ARIADNE_INVARIANT_SUBSPACE_HPP
