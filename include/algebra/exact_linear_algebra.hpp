/***************************************************************************
 *            algebra/exact_linear_algebra.hpp
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

/*! \file algebra/exact_linear_algebra.hpp
 *  \brief Exact matrix rank, null space and symmetric definiteness.
 */

#ifndef ARIADNE_EXACT_LINEAR_ALGEBRA_HPP
#define ARIADNE_EXACT_LINEAR_ALGEBRA_HPP

#include "numeric/rational.hpp"
#include "algebra/matrix.hpp"

namespace Ariadne {

//! \ingroup LinearAlgebraSubModule
//! \brief Exact right null-space and rank of a rational matrix.
//! \details For an m-by-n matrix A, basis is n-by-(n-rank). Its columns
//! are linearly independent and satisfy A*basis=0 exactly.
struct ExactNullSpaceResult {
    Matrix<Rational> basis; //!< A basis of the right null space.
    SizeType rank;          //!< Exact row and column rank.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Compute a basis of the right null space of a rational matrix.
//! \details Applies Gauss-Jordan elimination over the exact rational field;
//! a pivot is accepted if and only if it is mathematically nonzero.
//! The algorithm supports singular, rectangular and empty matrices;
//! zero rows and zero columns are retained in the output dimensions.
//! There are no numerical rank thresholds.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, Gaussian elimination.
ExactNullSpaceResult exact_null_space(const Matrix<Rational>& A);

//! \ingroup LinearAlgebraSubModule
//! \brief Definiteness of the quadratic form associated with a symmetric matrix.
//! \details INDEFINITE means *not positive semidefinite*: this includes
//! negative-definite matrices. The empty matrix is POSITIVE_DEFINITE
//! according to the vacuous condition on nonzero vectors.
enum class ExactSymmetricDefiniteness {
    POSITIVE_DEFINITE,      //!< Every nonzero x has x^T A x > 0.
    POSITIVE_SEMIDEFINITE,  //!< Every x has x^T A x >= 0, but A is singular.
    INDEFINITE              //!< Some x has x^T A x < 0.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Classify the definiteness of a symmetric rational matrix exactly.
//! \details Requires a square and symmetric matrix (enforced by preconditions).
//! Successive exact Schur complements eliminate positive pivots; a zero
//! pivot requires an entirely zero trailing row for positive semidefiniteness.
//! The computation never uses floating-point eigenvalues or tolerances.
//! \par Reference
//! R. A. Horn and C. R. Johnson, \em Matrix \em Analysis, 2nd ed.,
//! Cambridge University Press, 2012, Schur complements and definiteness.
ExactSymmetricDefiniteness exact_symmetric_definiteness(const Matrix<Rational>& A);

} // namespace Ariadne

#endif
