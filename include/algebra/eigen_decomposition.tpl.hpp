/***************************************************************************
 *            algebra/eigen_decomposition.tpl.hpp
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

/** @file algebra/eigen_decomposition.tpl.hpp
 *  @brief Template implementation of real eigendecomposition.
 */

#ifndef ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP
#define ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP

#include <vector>
#include <type_traits>
#include <utility>

namespace Ariadne {

template<class P, class PR>
bool EigenDecomposition<P, PR>::_positive(const NumericType& x) const {
    if constexpr (std::is_same_v<P, ValidatedTag>) {
        return x.lower_raw() > NumericType(0u, _precision).lower_raw();
    } else {
        return is_finite(x.raw()) && x.raw() > NumericType(0u, _precision).raw();
    }
}

template<class P, class PR>
bool EigenDecomposition<P, PR>::_negative(const NumericType& x) const {
    if constexpr (std::is_same_v<P, ValidatedTag>) {
        return x.upper_raw() < NumericType(0u, _precision).upper_raw();
    } else {
        return is_finite(x.raw()) && x.raw() < NumericType(0u, _precision).raw();
    }
}

template<class P, class PR>
bool EigenDecomposition<P, PR>::_zero(const NumericType& x) const {
    if constexpr (std::is_same_v<P, ValidatedTag>) {
        const auto zero = NumericType(0u, _precision).lower_raw();
        return x.lower_raw() == zero && x.upper_raw() == zero;
    } else {
        return is_finite(x.raw()) && x.raw() == NumericType(0u, _precision).raw();
    }
}

template<class P, class PR>
bool EigenDecomposition<P, PR>::_nonzero(const NumericType& x) const {
    return this->_positive(x) || this->_negative(x);
}

template<class P, class PR>
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::operator()(const Matrix<NumericType>& M) const {
    using X = NumericType;
    using Status = EigenDecompositionStatus;
    ARIADNE_PRECONDITION(M.row_size() == M.column_size());

    const SizeType n = M.row_size();
    if (n == 0u) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    }
    const X zero(0u, _precision);
    const X one(1u, _precision);
    Vector<X> values(n, zero);
    Matrix<X> vectors(n, n, zero);

    if (n > 2u) {
        bool upper_triangular = true;
        bool lower_triangular = true;
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                if (row > col && !this->_zero(M[row][col])) {
                    upper_triangular = false;
                }
                if (row < col && !this->_zero(M[row][col])) {
                    lower_triangular = false;
                }
            }
        }
        if (!upper_triangular && !lower_triangular) {
            if constexpr (std::is_same_v<P, ApproximateTag>) {
                return this->_symmetric(M);
            } else {
                return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
            }
        }

        // An interval-independent ordering is necessary to identify the
        // eigenpair belonging to each output column across the entire
        // represented family of real matrices.
        std::vector<SizeType> order;
        order.reserve(n);
        for (SizeType row = 0u; row != n; ++row) {
            SizeType position = 0u;
            for (; position != order.size(); ++position) {
                const X difference = M[row][row] - M[order[position]][order[position]];
                if (this->_negative(difference)) {
                    break;
                }
                if (!this->_positive(difference)) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
            order.insert(order.begin() + static_cast<std::ptrdiff_t>(position), row);
        }

        for (SizeType col = 0u; col != n; ++col) {
            const SizeType pivot = order[col];
            const X eigenvalue = M[pivot][pivot];
            values[col] = eigenvalue;
            vectors[pivot][col] = one;

            if (upper_triangular) {
                // The eigenvector has zero components above the pivot index;
                // solve the remaining equations by backward substitution.
                for (SizeType row = pivot; row > 0u;) {
                    --row;
                    X numerator(zero);
                    for (SizeType k = row + 1u; k <= pivot; ++k) {
                        numerator += M[row][k] * vectors[k][col];
                    }
                    const X divisor = M[row][row] - eigenvalue;
                    if (!this->_nonzero(divisor)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    vectors[row][col] = -numerator / divisor;
                }
            } else {
                // Lower triangular analogue: forward substitution.
                for (SizeType row = pivot + 1u; row != n; ++row) {
                    X numerator(zero);
                    for (SizeType k = pivot; k < row; ++k) {
                        numerator += M[row][k] * vectors[k][col];
                    }
                    const X divisor = M[row][row] - eigenvalue;
                    if (!this->_nonzero(divisor)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    vectors[row][col] = -numerator / divisor;
                }
            }
        }

        // The eigenvector matrix is triangular after a column permutation,
        // with exact unit pivots and exact structural zero entries. This
        // proves invertibility without a pessimistic interval determinant.
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(values), std::move(vectors)}};
    }

    if (n == 1u) {
        values[0u] = M[0u][0u];
        vectors[0u][0u] = one;
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(values), std::move(vectors)}};
    }

    const X& a = M[0u][0u];
    const X& b = M[0u][1u];
    const X& c = M[1u][0u];
    const X& d = M[1u][1u];
    const X two_value(2u, _precision);
    const X four(4u, _precision);

    // sqr rather than multiplication is necessary for interval inputs:
    // (a-d)*(a-d) may include negative values due to dependency.
    const X discriminant = sqr(a - d) + four * b * c;
    if (!this->_positive(discriminant)) {
        return {Status::INCONCLUSIVE, std::nullopt};
    }

    const X root = sqrt(discriminant);
    const X trace = a + d;
    const X characteristic_determinant = a * d - b * c;

    // Compute the dominant eigenvalue without cancellation and derive the
    // other from their product. With interval inputs, the quotient still
    // encloses the corresponding exact eigenvalue for every input matrix.
    // A trace interval crossing zero uses the direct formulas instead.
    if (this->_positive(trace)) {
        values[1u] = (trace + root) / two_value;
        if (!this->_nonzero(values[1u])) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        values[0u] = characteristic_determinant / values[1u];
    } else if (this->_negative(trace)) {
        values[0u] = (trace - root) / two_value;
        if (!this->_nonzero(values[0u])) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        values[1u] = characteristic_determinant / values[0u];
    } else {
        values[0u] = (trace - root) / two_value;
        values[1u] = (trace + root) / two_value;
    }

    // A diagonal matrix has a canonical exact eigenbasis even if both
    // off-diagonal coefficients are identically zero.
    if (this->_zero(b) && this->_zero(c)) {
        if (this->_negative(a - d)) {
            vectors[0u][0u] = one;
            vectors[1u][1u] = one;
        } else if (this->_positive(a - d)) {
            vectors[1u][0u] = one;
            vectors[0u][1u] = one;
        } else {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
    } else if (this->_nonzero(b)) {
        // v = (1, (lambda-a)/b). Each column encloses the normalized
        // exact eigenvector for every represented input matrix.
        for (SizeType j = 0u; j != 2u; ++j) {
            vectors[0u][j] = one;
            vectors[1u][j] = (values[j] - a) / b;
        }
    } else if (this->_nonzero(c)) {
        // v = ((lambda-d)/c, 1), when b is not provably nonzero.
        for (SizeType j = 0u; j != 2u; ++j) {
            vectors[0u][j] = (values[j] - d) / c;
            vectors[1u][j] = one;
        }
    } else {
        // A wide interval input may contain both zero and nonzero
        // off-diagonals. Do not choose an unverified normalization.
        return {Status::INCONCLUSIVE, std::nullopt};
    }

    // For validated arithmetic, certification that the interval determinant
    // excludes zero implies invertibility of every enclosed eigenvector
    // matrix. A failed test is inconclusive, not proof of singularity.
    const X eigenvector_determinant = vectors[0u][0u] * vectors[1u][1u]
                                     - vectors[0u][1u] * vectors[1u][0u];
    if (!this->_nonzero(eigenvector_determinant)) {
        return {Status::INCONCLUSIVE, std::nullopt};
    }

    return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(values), std::move(vectors)}};
}


/** @brief Diagonalize an approximately symmetric real matrix using Jacobi rotations. */
template<class P, class PR>
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::_symmetric(const Matrix<NumericType>& input) const {
    using Status = EigenDecompositionStatus;
    if constexpr (!std::is_same_v<P, ApproximateTag>) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    } else {
        using X = NumericType;
        const SizeType n = input.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
        const X two_value(2u, _precision);

        // Require equal stored symmetric entries; silently symmetrizing
        // a genuinely nonsymmetric matrix would give incorrect eigenpairs.
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = row + 1u; col != n; ++col) {
                if (input[row][col].raw() != input[col][row].raw()) {
                    return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
                }
            }
        }

        // Work at the input precision, not via get_d(): an MP value may
        // be unrepresentable or indistinguishable from zero in binary64.
        using Raw = typename X::RawType;
        const X epsilon(Raw::eps(_precision));
        X scale(zero);
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                if (!is_finite(input[row][col].raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                const X magnitude = abs(input[row][col]);
                if (magnitude.raw() > scale.raw()) {
                    scale = magnitude;
                }
            }
        }
        if (scale.raw() == zero.raw()) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }

        Matrix<X> rotated(input);
        Matrix<X> basis = Matrix<X>::identity(n, _precision);
        const X tolerance = X(64u, _precision) * epsilon * scale;
        const SizeType iteration_limit = 64u * n * n;
        bool converged = false;

        for (SizeType iteration = 0u; iteration != iteration_limit; ++iteration) {
            SizeType pivot_row = 0u;
            SizeType pivot_col = 1u;
            X off_diagonal_maximum(zero);

            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = row + 1u; col != n; ++col) {
                    if (!is_finite(rotated[row][col].raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    const X magnitude = abs(rotated[row][col]);
                    if (magnitude.raw() > off_diagonal_maximum.raw()) {
                        off_diagonal_maximum = magnitude;
                        pivot_row = row;
                        pivot_col = col;
                    }
                }
            }
            if (off_diagonal_maximum.raw() <= tolerance.raw()) {
                converged = true;
                break;
            }

            const X app = rotated[pivot_row][pivot_row];
            const X aqq = rotated[pivot_col][pivot_col];
            const X apq = rotated[pivot_row][pivot_col];
            const X tau = (aqq - app) / (two_value * apq);
            if (!is_finite(tau.raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            // The stable Jacobi tangent avoids cancellation in the
            // quadratic formula for the plane rotation.
            const X tangent = (tau.raw() < zero.raw() ? -one : one) /
                              (abs(tau) + sqrt(one + sqr(tau)));
            const X cosine = one / sqrt(one + sqr(tangent));
            const X sine = tangent * cosine;

            for (SizeType row = 0u; row != n; ++row) {
                if (row != pivot_row && row != pivot_col) {
                    const X rp = rotated[row][pivot_row];
                    const X rq = rotated[row][pivot_col];
                    const X updated_p = cosine * rp - sine * rq;
                    const X updated_q = sine * rp + cosine * rq;
                    rotated[row][pivot_row] = updated_p;
                    rotated[pivot_row][row] = updated_p;
                    rotated[row][pivot_col] = updated_q;
                    rotated[pivot_col][row] = updated_q;
                }
                const X bp = basis[row][pivot_row];
                const X bq = basis[row][pivot_col];
                basis[row][pivot_row] = cosine * bp - sine * bq;
                basis[row][pivot_col] = sine * bp + cosine * bq;
            }
            rotated[pivot_row][pivot_row] = app - tangent * apq;
            rotated[pivot_col][pivot_col] = aqq + tangent * apq;
            rotated[pivot_row][pivot_col] = zero;
            rotated[pivot_col][pivot_row] = zero;
        }
        if (!converged) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }

        // Sort the columns, and enforce simple eigenvalues for this
        // implementation milestone.
        std::vector<SizeType> ordering;
        ordering.reserve(n);
        for (SizeType row = 0u; row != n; ++row) {
            const X eigenvalue = rotated[row][row];
            if (!is_finite(eigenvalue.raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            SizeType pos = 0u;
            for (; pos != ordering.size(); ++pos) {
                if (eigenvalue.raw() < rotated[ordering[pos]][ordering[pos]].raw()) {
                    break;
                }
            }
            ordering.insert(ordering.begin() + static_cast<std::ptrdiff_t>(pos), row);
        }

        Vector<X> eigenvalues(n, zero);
        Matrix<X> eigenvectors(n, n, zero);
        for (SizeType col = 0u; col != n; ++col) {
            const SizeType source = ordering[col];
            eigenvalues[col] = rotated[source][source];
            if (col > 0u) {
                const X gap = eigenvalues[col] - eigenvalues[col - 1u];
                if (!is_finite(gap.raw()) ||
                    gap.raw() <= (X(128u, _precision) * epsilon * scale).raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
            for (SizeType row = 0u; row != n; ++row) {
                eigenvectors[row][col] = basis[row][source];
            }
        }

        // Small residual and near-orthogonality are numerical quality
        // checks only, not a validated enclosure or mathematical proof.
        const X budget = X(4096u, _precision) * epsilon * scale *
                         X(static_cast<unsigned>(n), _precision);
        for (SizeType col = 0u; col != n; ++col) {
            for (SizeType row = 0u; row != n; ++row) {
                X action(zero);
                for (SizeType k = 0u; k != n; ++k) {
                    action += input[row][k] * eigenvectors[k][col];
                }
                const X residual = abs(action - eigenvalues[col] * eigenvectors[row][col]);
                if (!is_finite(residual.raw()) || residual.raw() > budget.raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(eigenvalues), std::move(eigenvectors)}};
    }
}

} // namespace Ariadne

#endif // ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP
