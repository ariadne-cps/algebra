/***************************************************************************
 *            algebra/cholesky_decomposition.hpp
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

/*! \file algebra/cholesky_decomposition.hpp
 *  \brief Validated interval Cholesky decomposition.
 */

#ifndef ARIADNE_CHOLESKY_DECOMPOSITION_HPP
#define ARIADNE_CHOLESKY_DECOMPOSITION_HPP

#include <optional>

#include "algebra/positive_definiteness.hpp"

namespace Ariadne {

//! \ingroup LinearAlgebraSubModule
//! \brief The result of an attempted validated interval Cholesky decomposition.
//! \details The factor is present if and only if the status is
//! POSITIVE_DEFINITE. Otherwise no complete validated factor is available.
//! A successful lower matrix encloses the lower-triangular Cholesky factor
//! of every real symmetric matrix represented by the interval input.
//! In particular, the result is not a pointwise equality between interval
//! matrices: L*transpose(L) is an interval enclosure of the input family.
struct ValidatedCholeskyResult {
    PositiveDefinitenessStatus status; //!< Rigorous universal classification.
    std::optional<Matrix<FloatDPBounds>> lower; //!< Validated lower factor, if certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Compute a lower-triangular Cholesky enclosure with outward rounding.
//! \param A Symmetric matrix with coefficients represented as FloatDPBounds.
//! \return A certified Cholesky enclosure, a proof of failure of positive
//! definiteness for all represented matrices, or INDETERMINATE.
//! \details For each represented real symmetric matrix a, when the status
//! is POSITIVE_DEFINITE, the unique Cholesky factor l (a=l*l^T,
//! with positive diagonal) belongs componentwise to the returned
//! interval matrix L. The upper triangular portion of L is exactly zero.
//!
//! The algorithm computes successive Schur pivots
//! d_i=a_ii-sum_{k<i} l_ik^2, and l_ij is formed by the corresponding
//! triangular recurrence. Each interval operation encloses its exact
//! real counterpart. A pivot known strictly positive permits sqrt(d_i);
//! a pivot known nonpositive, with all preceding pivots certified positive,
//! proves that no represented matrix is positive definite; otherwise the
//! method returns INDETERMINATE without dividing by an uncertain pivot.
//!
//! No floating-point tolerance is used. Dependency overestimation may
//! yield INDETERMINATE even for a family of positive-definite matrices.
//! An empty matrix is positive definite by convention.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, Cholesky factorization.
ValidatedCholeskyResult validated_cholesky_decomposition(
    const SymmetricMatrix<FloatDPBounds>& A);

} // namespace Ariadne

#endif
