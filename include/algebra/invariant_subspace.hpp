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

#include <optional>
#include <utility>
#include <vector>

#include "algebra/matrix.hpp"

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

} // namespace Ariadne

#endif // ARIADNE_INVARIANT_SUBSPACE_HPP
