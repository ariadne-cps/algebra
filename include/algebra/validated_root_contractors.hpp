/***************************************************************************
 *            algebra/validated_root_contractors.hpp
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

/** @file algebra/validated_root_contractors.hpp
 *  @brief Minimal interval Newton and Krawczyk image operations for vectors and matrices.
 */

#ifndef ARIADNE_VALIDATED_ROOT_CONTRACTORS_HPP
#define ARIADNE_VALIDATED_ROOT_CONTRACTORS_HPP

#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

//! \brief Compute one interval Newton image from numerical input.
//! \details Returns x0 - J(X)^(-1) f(x0). All inputs must be interval
//! enclosures; J is assumed to enclose the Jacobian throughout the
//! candidate box X (not passed because it is not needed by this formula).
//! A valid interval inverse is required. This helper performs no iteration,
//! intersection, existence/uniqueness test or user-function evaluation.
//! \pre Dimensions are consistent and J is square.
//! \warning An image alone is not a proof that a root exists.
//! \warning An inverse that cannot be established may fail; callers are
//! responsible for treating such failures as inconclusive.
template<class X>
Vector<X> interval_newton_image(const Vector<X>& center,
                                const Vector<X>& residual_at_center,
                                const Matrix<X>& jacobian) {
    const SizeType n = center.size();
    ARIADNE_PRECONDITION(residual_at_center.size() == n);
    ARIADNE_PRECONDITION(jacobian.row_size() == n);
    ARIADNE_PRECONDITION(jacobian.column_size() == n);

    const Matrix<X> inverse_jacobian = inverse(jacobian);
    const Vector<X> correction = inverse_jacobian * residual_at_center;
    return center - correction;
}

//! \brief Compute one Krawczyk image using a caller-supplied preconditioner.
//! \details K(X) = x0 - C*f(x0) + (I-C*J(X))*(X-x0).
//! The matrix C is chosen by the caller, typically as an approximate inverse
//! to the midpoint Jacobian; it need not itself be a validated inverse for
//! the image computation. In a root-existence/uniqueness proof, additional
//! hypotheses on C and strict inclusion of K(X) in the interior of X must
//! be verified by the higher-level solver.
//! The operation has no state, iteration, function evaluation or dependency
//! on equation-solving.
//! \pre All vector dimensions agree and J and C are n-by-n matrices.
//! \warning An image alone is not a proof that a root exists.
template<class X>
Vector<X> krawczyk_image(const Vector<X>& center,
                         const Vector<X>& residual_at_center,
                         const Matrix<X>& jacobian,
                         const Vector<X>& box,
                         const Matrix<X>& preconditioner) {
    const SizeType n = center.size();
    ARIADNE_PRECONDITION(residual_at_center.size() == n);
    ARIADNE_PRECONDITION(box.size() == n);
    ARIADNE_PRECONDITION(jacobian.row_size() == n);
    ARIADNE_PRECONDITION(jacobian.column_size() == n);
    ARIADNE_PRECONDITION(preconditioner.row_size() == n);
    ARIADNE_PRECONDITION(preconditioner.column_size() == n);

    const Matrix<X> identity = Matrix<X>::identity(n, center.zero_element());
    const Vector<X> correction = preconditioner * residual_at_center;
    const Matrix<X> defect = identity - preconditioner * jacobian;
    const Vector<X> radius = box - center;
    const Vector<X> remainder = defect * radius;
    return center - correction + remainder;
}

} // namespace Ariadne

#endif // ARIADNE_VALIDATED_ROOT_CONTRACTORS_HPP
