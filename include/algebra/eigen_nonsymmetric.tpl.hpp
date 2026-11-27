/***************************************************************************
 *            algebra/eigen_nonsymmetric.tpl.hpp
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

/** @file algebra/eigen_nonsymmetric.tpl.hpp
 *  @brief Approximate simple real eigenpairs of nonsymmetric dense matrices.
 */

#ifndef ARIADNE_EIGEN_NONSYMMETRIC_TPL_HPP
#define ARIADNE_EIGEN_NONSYMMETRIC_TPL_HPP

#include <exception>
#include <vector>

namespace Ariadne {

/** @brief Compute a real nonsymmetric eigendecomposition when numerically resolvable.
 *  @details Real-shift Householder QR with trailing-row deflation gives initial
 *  eigenvalues and approximate Schur vectors. Eigenpairs are subsequently
 *  refined with Newton's method on the normalized equations
 *  (A-lambda I)v=0, v[pivot]=1. All decisions use the input DP/MP precision.
 *  This is an approximate method, NOT a proof of diagonalizability.
 */
template<class P, class PR>
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::_nonsymmetric(const Matrix<NumericType>& input) const {
    using Status = EigenDecompositionStatus;
    if constexpr (!std::is_same_v<P, ApproximateTag>) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    } else {
        using X = NumericType;
        using Raw = typename X::RawType;
        const SizeType n = input.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
        const X two_value(2u, _precision);
        const X epsilon(Raw::eps(_precision));

        X scale(zero);
        for (SizeType i = 0u; i != n; ++i) {
            for (SizeType j = 0u; j != n; ++j) {
                if (!is_finite(input[i][j].raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                const X magnitude = abs(input[i][j]);
                if (magnitude.raw() > scale.raw()) {
                    scale = magnitude;
                }
            }
        }
        if (scale.raw() == zero.raw()) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }

        const X tolerance = X(128u, _precision) * epsilon * scale;
        const X separation = X(256u, _precision) * epsilon * scale;
        const X residual_budget = X(8192u, _precision) * epsilon * scale *
                                  X(static_cast<unsigned>(n), _precision);
        Matrix<X> triangular(input);
        Matrix<X> orthogonal = Matrix<X>::identity(n, _precision);
        SizeType active = n;
        const SizeType iteration_limit = 256u * n * n;

        // The orthogonal similarity A=O*T*O^T is approximate. The
        // later Newton solve is performed on the ORIGINAL input.
        for (SizeType iteration = 0u; active > 1u &&
                iteration != iteration_limit; ++iteration) {
            bool last_row_small = true;
            for (SizeType j = 0u; j + 1u < active; ++j) {
                const X magnitude = abs(triangular[active - 1u][j]);
                if (!is_finite(magnitude.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                if (magnitude.raw() > tolerance.raw()) {
                    last_row_small = false;
                }
            }
            if (last_row_small) {
                for (SizeType j = 0u; j + 1u < active; ++j) {
                    triangular[active - 1u][j] = zero;
                }
                --active;
                continue;
            }

            // A real Wilkinson shift near the bottom-right diagonal.
            // Complex or repeated final real 2x2 blocks cannot yield
            // the requested complete SIMPLE REAL decomposition.
            const X& a = triangular[active - 2u][active - 2u];
            const X& b = triangular[active - 2u][active - 1u];
            const X& c = triangular[active - 1u][active - 2u];
            const X& d = triangular[active - 1u][active - 1u];
            const X discriminant = sqr(a - d) + X(4u, _precision) * b * c;
            if (!is_finite(discriminant.raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            X shift(d);
            if (discriminant.raw() > zero.raw()) {
                const X root = sqrt(discriminant);
                const X low = (a + d - root) / two_value;
                const X high = (a + d + root) / two_value;
                shift = abs(low - d).raw() < abs(high - d).raw() ? low : high;
            } else if (active == 2u) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            // Numerically stable Householder QR of T_active - shift*I.
            // The QR factorization exists even when the shift is exactly
            // an eigenvalue, unlike normalizing Gram-Schmidt columns.
            Matrix<X> R(active, active, zero);
            Matrix<X> Q = Matrix<X>::identity(active, _precision);
            for (SizeType i = 0u; i != active; ++i) {
                for (SizeType j = 0u; j != active; ++j) {
                    R[i][j] = triangular[i][j] - (i == j ? shift : zero);
                }
            }
            for (SizeType k = 0u; k + 1u < active; ++k) {
                X norm_squared(zero);
                for (SizeType i = k; i != active; ++i) {
                    norm_squared += sqr(R[i][k]);
                }
                if (!is_finite(norm_squared.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                if (norm_squared.raw() == zero.raw()) {
                    continue;
                }
                const X magnitude = sqrt(norm_squared);
                const X alpha = R[k][k].raw() < zero.raw() ? magnitude : -magnitude;
                std::vector<X> reflector(active - k, zero);
                reflector[0u] = R[k][k] - alpha;
                X reflector_squared(sqr(reflector[0u]));
                for (SizeType i = k + 1u; i != active; ++i) {
                    reflector[i - k] = R[i][k];
                    reflector_squared += sqr(reflector[i - k]);
                }
                if (!is_finite(reflector_squared.raw()) ||
                    !(reflector_squared.raw() > zero.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                const X beta = two_value / reflector_squared;
                if (!is_finite(beta.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }

                // R <- H_k R, Q <- Q H_k, where H_k is symmetric.
                for (SizeType col = k; col != active; ++col) {
                    X dot(zero);
                    for (SizeType row = k; row != active; ++row) {
                        dot += reflector[row - k] * R[row][col];
                    }
                    const X multiple = beta * dot;
                    for (SizeType row = k; row != active; ++row) {
                        R[row][col] -= multiple * reflector[row - k];
                    }
                }
                for (SizeType row = 0u; row != active; ++row) {
                    X dot(zero);
                    for (SizeType col = k; col != active; ++col) {
                        dot += Q[row][col] * reflector[col - k];
                    }
                    const X multiple = beta * dot;
                    for (SizeType col = k; col != active; ++col) {
                        Q[row][col] -= multiple * reflector[col - k];
                    }
                }
            }

            Matrix<X> next = R * Q;
            for (SizeType i = 0u; i != active; ++i) {
                next[i][i] += shift;
            }
            Matrix<X> updated_columns(n, active, zero);
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != active; ++j) {
                    for (SizeType k = 0u; k != active; ++k) {
                        updated_columns[i][j] += orthogonal[i][k] * Q[k][j];
                    }
                }
            }
            // A QR step acts on the entire leading block by similarity.
            // Transform the active-to-deflated upper-right block too:
            // T(0:active,active:n) <- Q^T*T(0:active,active:n).
            // Omitting this update makes Schur eigenvectors inconsistent
            // after the first deflation, even when diagonal values converge.
            for (SizeType col = active; col != n; ++col) {
                std::vector<X> transformed(active, zero);
                for (SizeType row = 0u; row != active; ++row) {
                    for (SizeType k = 0u; k != active; ++k) {
                        transformed[row] += Q[k][row] * triangular[k][col];
                    }
                }
                for (SizeType row = 0u; row != active; ++row) {
                    triangular[row][col] = transformed[row];
                }
            }
            for (SizeType i = 0u; i != active; ++i) {
                for (SizeType j = 0u; j != active; ++j) {
                    triangular[i][j] = next[i][j];
                }
            }
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != active; ++j) {
                    orthogonal[i][j] = updated_columns[i][j];
                }
            }
        }
        if (active > 1u) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }

        Vector<X> eigenvalues(n, zero);
        Matrix<X> eigenvectors(n, n, zero);
        for (SizeType col = 0u; col != n; ++col) {
            const X initial_value = triangular[col][col];
            if (!is_finite(initial_value.raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            // Back-substitute in the approximate upper Schur form.
            Vector<X> schur_vector(n, zero);
            schur_vector[col] = one;
            for (SizeType row = col; row > 0u;) {
                --row;
                X numerator(zero);
                for (SizeType k = row + 1u; k <= col; ++k) {
                    numerator += triangular[row][k] * schur_vector[k];
                }
                const X divisor = initial_value - triangular[row][row];
                if (!is_finite(divisor.raw()) ||
                    abs(divisor).raw() <= separation.raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                schur_vector[row] = numerator / divisor;
            }
            Vector<X> initial_vector(n, zero);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType k = 0u; k != n; ++k) {
                    initial_vector[row] += orthogonal[row][k] * schur_vector[k];
                }
            }
            SizeType pivot = 0u;
            for (SizeType row = 1u; row != n; ++row) {
                if (abs(initial_vector[row]).raw() >
                    abs(initial_vector[pivot]).raw()) {
                    pivot = row;
                }
            }
            if (!is_finite(initial_vector[pivot].raw()) ||
                initial_vector[pivot].raw() == zero.raw()) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            std::vector<SizeType> free_rows;
            free_rows.reserve(n - 1u);
            for (SizeType row = 0u; row != n; ++row) {
                if (row != pivot) {
                    free_rows.push_back(row);
                }
            }

            // Newton's square system in n real unknowns is
            // F_i(lambda,v)=(Av)_i-lambda*v_i, with v_pivot=1.
            Vector<X> z(n, zero);
            z[0u] = initial_value;
            for (SizeType k = 0u; k != n - 1u; ++k) {
                z[k + 1u] = initial_vector[free_rows[k]] /
                             initial_vector[pivot];
            }
            for (SizeType step = 0u; step != 12u; ++step) {
                Vector<X> v(n, zero);
                v[pivot] = one;
                for (SizeType k = 0u; k != n - 1u; ++k) {
                    v[free_rows[k]] = z[k + 1u];
                }
                Vector<X> residual(n, zero);
                Matrix<X> jacobian(n, n, zero);
                for (SizeType row = 0u; row != n; ++row) {
                    for (SizeType k = 0u; k != n; ++k) {
                        residual[row] += input[row][k] * v[k];
                    }
                    residual[row] -= z[0u] * v[row];
                    jacobian[row][0u] = -v[row];
                    for (SizeType k = 0u; k != n - 1u; ++k) {
                        const SizeType j = free_rows[k];
                        jacobian[row][k + 1u] = input[row][j] -
                            (row == j ? z[0u] : zero);
                    }
                }
                Matrix<X> inverted(n, n, zero);
                try {
                    inverted = inverse(jacobian);
                } catch (const std::exception&) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                const Vector<X> correction = inverted * residual;
                for (SizeType k = 0u; k != n; ++k) {
                    z[k] -= correction[k];
                    if (!is_finite(z[k].raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            }
            eigenvalues[col] = z[0u];
            eigenvectors[pivot][col] = one;
            for (SizeType k = 0u; k != n - 1u; ++k) {
                eigenvectors[free_rows[k]][col] = z[k + 1u];
            }
        }

        // Sort the refined approximate eigenpairs and require separation.
        std::vector<SizeType> order;
        order.reserve(n);
        for (SizeType col = 0u; col != n; ++col) {
            if (!is_finite(eigenvalues[col].raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            SizeType position = 0u;
            for (; position != order.size(); ++position) {
                if (eigenvalues[col].raw() < eigenvalues[order[position]].raw()) {
                    break;
                }
            }
            order.insert(order.begin() + static_cast<std::ptrdiff_t>(position), col);
        }
        Vector<X> sorted_values(n, zero);
        Matrix<X> sorted_vectors(n, n, zero);
        for (SizeType col = 0u; col != n; ++col) {
            const SizeType source = order[col];
            sorted_values[col] = eigenvalues[source];
            if (col != 0u) {
                const X gap = sorted_values[col] - sorted_values[col - 1u];
                if (!is_finite(gap.raw()) ||
                    gap.raw() <= separation.raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
            for (SizeType row = 0u; row != n; ++row) {
                sorted_vectors[row][col] = eigenvectors[row][source];
            }
        }

        // Check each eigenpair against the ORIGINAL input, not the Schur
        // surrogate, then reject ill-conditioned or singular vector bases.
        for (SizeType col = 0u; col != n; ++col) {
            X vector_scale(one);
            for (SizeType row = 0u; row != n; ++row) {
                const X magnitude = abs(sorted_vectors[row][col]);
                if (!is_finite(magnitude.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                if (magnitude.raw() > vector_scale.raw()) {
                    vector_scale = magnitude;
                }
            }
            for (SizeType row = 0u; row != n; ++row) {
                X action(zero);
                for (SizeType k = 0u; k != n; ++k) {
                    action += input[row][k] * sorted_vectors[k][col];
                }
                const X residual = abs(action -
                    sorted_values[col] * sorted_vectors[row][col]);
                if (!is_finite(residual.raw()) ||
                    residual.raw() > (residual_budget * vector_scale).raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }

        Matrix<X> inverse_vectors(n, n, zero);
        try {
            inverse_vectors = inverse(sorted_vectors);
        } catch (const std::exception&) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        const Matrix<X> product = sorted_vectors * inverse_vectors;
        const X basis_budget = X(8192u, _precision) * epsilon *
                               X(static_cast<unsigned>(n), _precision);
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                const X difference = abs(product[row][col] -
                    (row == col ? one : zero));
                if (!is_finite(difference.raw()) ||
                    difference.raw() > basis_budget.raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{
            std::move(sorted_values), std::move(sorted_vectors)}};
    }
}

} // namespace Ariadne

#endif // ARIADNE_EIGEN_NONSYMMETRIC_TPL_HPP
