/***************************************************************************
 *            algebra/linear_independence.hpp
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

/*! \file algebra/linear_independence.hpp
 *  \brief Validated certificates of linear independence and dependence.
 */

#ifndef ARIADNE_LINEAR_INDEPENDENCE_HPP
#define ARIADNE_LINEAR_INDEPENDENCE_HPP

#include "algebra/positive_definiteness.hpp"

namespace Ariadne {

//! \ingroup LinearAlgebraSubModule
//! \brief Rigorous classification of independence of interval matrix rows or columns.
//! \details INDEPENDENT means *every* real matrix represented by the input
//! has independent selected vectors; DEPENDENT means *every* represented
//! matrix has dependent selected vectors. INDETERMINATE means neither
//! universal claim was established. A zero-length family of vectors is
//! independent by convention. The classification is not a numerical
//! estimate of rank and uses no singular-value or pivot threshold.
enum class LinearIndependenceStatus {
    INDEPENDENT,    //!< All admissible real matrices have independent selected vectors.
    DEPENDENT,      //!< All admissible real matrices have dependent selected vectors.
    INDETERMINATE   //!< Neither statement could be certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Certify linear independence of the rows of an interval matrix.
//! \param A An m-by-n matrix with validated interval coefficients.
//! \return INDEPENDENT if every admissible real matrix has row rank m;
//! DEPENDENT if every admissible real matrix has row rank less than m;
//! otherwise INDETERMINATE.
//! \details Computes an interval enclosure of the Gram matrix A*A^T.
//! Each real Gram matrix is positive semidefinite and is positive definite
//! exactly when the original rows are linearly independent. The Gram
//! matrix is passed to validated_positive_definiteness(). If m>n,
//! dependence is certified immediately from dimension. For m=0,
//! independence holds vacuously even if n>0. Squared diagonal terms use
//! interval sqr() to avoid the dependency error of multiplying the same
//! interval by itself. Interval overestimation can make the answer
//! INDETERMINATE even when the real rows are all independent.
//! No numerical rank threshold is introduced.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, normal equations,
//! Gram matrices and full-rank conditions.
LinearIndependenceStatus validated_independent_rows(
    const Matrix<FloatDPBounds>& A);

//! \ingroup LinearAlgebraSubModule
//! \brief Certify linear independence of the columns of an interval matrix.
//! \param A An m-by-n matrix with validated interval coefficients.
//! \return INDEPENDENT if every admissible real matrix has column rank n;
//! DEPENDENT if every admissible real matrix has column rank less than n;
//! otherwise INDETERMINATE.
//! \details Computes an interval enclosure of A^T*A and certifies its
//! positive definiteness. If n>m, dependence follows from dimension.
//! Zero columns are independent vacuously. This is a sound, possibly
//! inconclusive test; it is not a general-purpose null-space computation.
//! \see validated_independent_rows
LinearIndependenceStatus validated_independent_columns(
    const Matrix<FloatDPBounds>& A);

} // namespace Ariadne

#endif
