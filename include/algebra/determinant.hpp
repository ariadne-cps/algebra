/***************************************************************************
 *            algebra/determinant.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/determinant.hpp
 *  @brief Approximate and validated determinants of square matrices.
 */

#ifndef ARIADNE_DETERMINANT_HPP
#define ARIADNE_DETERMINANT_HPP

#include "algebra/matrix.hpp"
#include "foundation/paradigm.hpp"

namespace Ariadne {

//! \brief Compute the determinant of a square approximate or validated matrix.
//! \details For ApproximateTag, uses pivoted Gaussian elimination (O(n^3))
//! without any enclosure guarantee. For ValidatedTag, uses the division-free
//! Faddeev-LeVerrier recurrence (O(n^4)) in outward-rounded interval
//! arithmetic. The latter encloses the determinant of every matrix in the
//! represented interval family, even when singular or containing zero pivots.
//! Both methods return one for the empty square matrix (0x0).
//! \pre A is square.
//! \warning The validated method may overestimate substantially for
//! high dimensions or interval matrices with strong dependencies.
template<class X>
X determinant(const Matrix<X>& A);

} // namespace Ariadne

#include "algebra/determinant.tpl.hpp"

#endif // ARIADNE_DETERMINANT_HPP
