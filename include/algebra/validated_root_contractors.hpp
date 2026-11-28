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
#include <type_traits>

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

//! \brief Verify an inclusion-and-contraction certificate for one root.
//! \details For each fixed continuously differentiable real function f,
//! assuming residual_at_center encloses f(center) and jacobian encloses
//! Df(x) for every x in box, the conditions K(box) strictly inside box
//! and ||I-C*jacobian||_infinity < 1 ensure a unique root in box.
//! Both conditions are checked with outward-rounded interval arithmetic.
//! This is valid also for families of functions satisfying the assumptions.
//! \pre Dimensions are consistent and box is a nonempty convex interval box.
//! \warning The caller must establish the stated residual and Jacobian
//! enclosure assumptions; this operation cannot infer them from matrices.
//! \warning Failure returns false (inconclusive), not proof of no root.
template<class X>
bool krawczyk_certifies_unique_root(const Vector<X>& center,
                                   const Vector<X>& residual_at_center,
                                   const Matrix<X>& jacobian,
                                   const Vector<X>& box,
                                   const Matrix<X>& preconditioner) {
    static_assert(std::is_same_v<typename X::Paradigm, ValidatedTag>,
                  "Krawczyk root certification requires validated interval arithmetic");
    const SizeType n = center.size();
    ARIADNE_PRECONDITION(residual_at_center.size() == n);
    ARIADNE_PRECONDITION(box.size() == n);
    ARIADNE_PRECONDITION(jacobian.row_size() == n);
    ARIADNE_PRECONDITION(jacobian.column_size() == n);
    ARIADNE_PRECONDITION(preconditioner.row_size() == n);
    ARIADNE_PRECONDITION(preconditioner.column_size() == n);
    if (n == 0u) {
        return false;
    }

    const X zero = center.zero_element();
    const X one = zero + X(1u, center.element_characteristics());
    for (SizeType i = 0u; i != n; ++i) {
        // The expansion point is fixed, not an unknown interval variable.
        if (center[i].lower_raw() != center[i].upper_raw() ||
            !is_finite(center[i].lower_raw()) ||
            !is_finite(box[i].lower_raw()) ||
            !is_finite(box[i].upper_raw()) ||
            !is_finite(residual_at_center[i].lower_raw()) ||
            !is_finite(residual_at_center[i].upper_raw()) ||
            !(box[i].lower_raw() < center[i].lower_raw() &&
              center[i].upper_raw() < box[i].upper_raw())) {
            return false;
        }
    }

    const Matrix<X> defect =
        Matrix<X>::identity(n, center.element_characteristics()) -
        preconditioner * jacobian;
    for (SizeType row = 0u; row != n; ++row) {
        X row_sum(zero);
        for (SizeType col = 0u; col != n; ++col) {
            if (!is_finite(defect[row][col].lower_raw()) ||
                !is_finite(defect[row][col].upper_raw())) {
                return false;
            }
            row_sum += abs(defect[row][col]);
        }
        if (!(row_sum.upper_raw() < one.lower_raw())) {
            return false;
        }
    }

    const Vector<X> image =
        krawczyk_image(center, residual_at_center, jacobian, box, preconditioner);
    for (SizeType i = 0u; i != n; ++i) {
        if (!is_finite(image[i].lower_raw()) ||
            !is_finite(image[i].upper_raw()) ||
            !(box[i].lower_raw() < image[i].lower_raw() &&
              image[i].upper_raw() < box[i].upper_raw())) {
            return false;
        }
    }
    return true;
}

} // namespace Ariadne

#endif // ARIADNE_VALIDATED_ROOT_CONTRACTORS_HPP
