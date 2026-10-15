/***************************************************************************
 *            algebra/positive_definite_solve.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file algebra/positive_definite_solve.hpp
 *  \brief Validated solution enclosures for symmetric positive-definite systems.
 */

#ifndef ARIADNE_POSITIVE_DEFINITE_SOLVE_HPP
#define ARIADNE_POSITIVE_DEFINITE_SOLVE_HPP

#include <optional>

#include "algebra/cholesky_decomposition.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

//! \ingroup LinearAlgebraSubModule
//! \brief Result of a validated symmetric positive-definite linear solve.
//! \tparam T Either Matrix<FloatDPBounds> or Vector<FloatDPBounds>.
//! \details A solution is present exactly when positive definiteness was
//! certified. It contains every exact solution for every admissible real
//! symmetric coefficient matrix and every admissible right-hand side.
//! An absent solution is not a claim that the linear system is unsolvable:
//! the positive-definiteness precondition could not be established.
template<class T> struct ValidatedPositiveDefiniteSolveResult {
    PositiveDefinitenessStatus status; //!< Classification of the coefficient matrix.
    std::optional<T> solution;         //!< Validated solution enclosure if certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Enclose the solutions of AX=B using validated Cholesky factorisation.
//! \param A Symmetric interval coefficient matrix.
//! \param B Interval right-hand sides, with exactly A.size() rows.
//! \return POSITIVE_DEFINITE together with an enclosure of every solution
//! if Cholesky can be certified; otherwise the failure or indeterminate
//! status of the factorisation, without a solution.
//! \pre B.row_size()==A.size().
//! \details For every real symmetric a in A and every real b in B,
//! the returned X includes the unique solution x of ax=b componentwise.
//! The algorithm first encloses a=L*L^T via validated Cholesky, then
//! applies interval forward and backward substitution. All arithmetic
//! is outward rounded. Correlations between entries are not retained;
//! the returned interval matrix may be wider than the optimal enclosure.
//! No inverse is constructed, and no numerical rank threshold is used.
//! Zero-dimensional matrices and a zero-column right-hand side are valid.
//! A non-SPD status does not imply that AX=B is inconsistent.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, Cholesky factorization
//! and triangular substitution.
ValidatedPositiveDefiniteSolveResult<Matrix<FloatDPBounds>>
validated_positive_definite_solve(
    const SymmetricMatrix<FloatDPBounds>& A,
    const Matrix<FloatDPBounds>& B);

//! \ingroup LinearAlgebraSubModule
//! \brief Enclose the solutions of Ax=b via validated Cholesky.
//! \param A Symmetric interval coefficient matrix.
//! \param b Interval right-hand side with exactly A.size() elements.
//! \return A certified interval solution iff the input is certified SPD;
//! otherwise the factorisation status and no solution.
//! \pre b.size()==A.size().
//! \details Equivalent to the matrix right-hand-side overload with
//! one column, including its universal enclosure guarantee.
ValidatedPositiveDefiniteSolveResult<Vector<FloatDPBounds>>
validated_positive_definite_solve(
    const SymmetricMatrix<FloatDPBounds>& A,
    const Vector<FloatDPBounds>& b);

} // namespace Ariadne

#endif
