/***************************************************************************
 *            algebra/determinant.tpl.hpp
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

/** @file algebra/determinant.tpl.hpp
 *  @brief Template implementation of approximate and validated matrix determinants.
 */

#ifndef ARIADNE_DETERMINANT_TPL_HPP
#define ARIADNE_DETERMINANT_TPL_HPP

#include <type_traits>
#include <utility>

namespace Ariadne {

template<class X>
X determinant(const Matrix<X>& A) {
    static_assert(std::is_same_v<typename X::Paradigm, ApproximateTag> ||
                  std::is_same_v<typename X::Paradigm, ValidatedTag>,
                  "determinant currently accepts approximate and validated floating-point matrices");
    ARIADNE_PRECONDITION(A.row_size() == A.column_size());

    const SizeType n = A.row_size();
    const auto precision = A.element_characteristics();
    const X zero(0u, precision);
    const X one(1u, precision);
    if (n == 0u) {
        return one;
    }

    if constexpr (std::is_same_v<typename X::Paradigm, ApproximateTag>) {
        // Partial row pivoting avoids dividing by zero whenever the
        // remaining column has a nonzero approximate entry.
        Matrix<X> U(A);
        X result(one);
        bool odd_permutation = false;
        for (SizeType k = 0u; k != n; ++k) {
            SizeType pivot = k;
            for (SizeType i = k + 1u; i != n; ++i) {
                if (decide(abs(U[i][k]) > abs(U[pivot][k]))) {
                    pivot = i;
                }
            }
            if (U[pivot][k].raw() == zero.raw()) {
                return zero;
            }
            if (pivot != k) {
                for (SizeType j = 0u; j != n; ++j) {
                    const X original = U[k][j];
                    U[k][j] = U[pivot][j];
                    U[pivot][j] = original;
                }
                odd_permutation = !odd_permutation;
            }
            const X diagonal = U[k][k];
            result *= diagonal;
            for (SizeType i = k + 1u; i != n; ++i) {
                const X factor = U[i][k] / diagonal;
                U[i][k] = zero;
                for (SizeType j = k + 1u; j != n; ++j) {
                    U[i][j] -= factor * U[k][j];
                }
            }
        }
        return odd_permutation ? -result : result;
    } else {
        // Coefficient c_n of det(tI-A) is (-1)^n det(A).
        // B_0=I; B_k=A B_(k-1)+c_k I and
        // c_k= -trace(A B_(k-1))/k. All operations enclose the
        // corresponding exact matrices for every interval input member.
        Matrix<X> B = Matrix<X>::identity(n, precision);
        X coefficient(zero);
        for (SizeType k = 1u; k <= n; ++k) {
            Matrix<X> product = A * B;
            X trace(zero);
            for (SizeType i = 0u; i != n; ++i) {
                trace += product[i][i];
            }
            coefficient = -trace / X(static_cast<unsigned>(k), precision);
            if (k != n) {
                for (SizeType i = 0u; i != n; ++i) {
                    product[i][i] += coefficient;
                }
                B = std::move(product);
            }
        }
        return n % 2u == 0u ? coefficient : -coefficient;
    }
}

} // namespace Ariadne

#endif // ARIADNE_DETERMINANT_TPL_HPP
