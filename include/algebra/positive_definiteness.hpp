/***************************************************************************
 *            algebra/positive_definiteness.hpp
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

/*! \file algebra/positive_definiteness.hpp
 *  \brief Validated certification of positive definiteness.
 */

#ifndef ARIADNE_POSITIVE_DEFINITENESS_HPP
#define ARIADNE_POSITIVE_DEFINITENESS_HPP

#include "numeric/numeric.hpp"
#include "algebra/matrix.hpp"
#include "algebra/symmetric_matrix.hpp"

namespace Ariadne {

//! \ingroup LinearAlgebraSubModule
//! \brief Rigorous three-way classification of positive definiteness.
//! \details For an interval symmetric matrix, the classification refers
//! to *every* real symmetric matrix represented by its entries.
//! Indeterminate is not a claim that some represented matrix is singular.
enum class PositiveDefinitenessStatus {
    POSITIVE_DEFINITE,      //!< Every represented real matrix is positive definite.
    NOT_POSITIVE_DEFINITE,  //!< No represented real matrix is positive definite.
    INDETERMINATE           //!< Neither assertion could be certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Certify positive definiteness using an interval LDL^T elimination.
//! \param A An interval symmetric matrix. Off-diagonal entries use the same
//! interval in both symmetric positions; each interval encloses all
//! possible exact coefficient values.
//! \return POSITIVE_DEFINITE if all successive interval Schur pivots
//! are strictly positive; NOT_POSITIVE_DEFINITE if a pivot is provably
//! nonpositive after all preceding pivots were proved positive;
//! otherwise INDETERMINATE.
//! \details Arithmetic on FloatDPBounds uses directed outward rounding.
//! The conclusion is universally quantified over all represented
//! symmetric real matrices. An undecidable pivot is *never* rounded
//! or treated as zero to select an arbitrary numerical rank.
//! The zero-dimensional matrix is positive definite by convention.
//! This is a sufficient interval test: dependency overestimation can
//! lead to INDETERMINATE even when every represented matrix is positive
//! definite. The method does not certify positive semidefiniteness
//! with zero pivots, and it does not compute eigenvalues.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, LDL^T factorization
//! and symmetric Schur complements.
PositiveDefinitenessStatus validated_positive_definiteness(
    const SymmetricMatrix<FloatDPBounds>& A);

} // namespace Ariadne

#endif
