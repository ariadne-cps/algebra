/***************************************************************************
 *            algebra/schur_decomposition.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/schur_decomposition.hpp
 *  @brief Approximate real Schur form by Hessenberg and shifted QR.
 */
#ifndef ARIADNE_SCHUR_DECOMPOSITION_HPP
#define ARIADNE_SCHUR_DECOMPOSITION_HPP

#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "algebra/matrix.hpp"

namespace Ariadne {

//! \brief Status of a Schur computation; INCONCLUSIVE is not nonexistence.
enum class SchurStatus { SUCCESS, INCONCLUSIVE, UNSUPPORTED_DIMENSION };

//! \brief Approximate A*Q=Q*T, with real quasi-triangular T.
//! \details Blocks have sizes one or two. A two-block is reserved for
//! a strictly nonreal conjugate pair. This is NOT a validated certificate.
template<class P, class PR>
struct SchurDecompositionData {
    using X = FloatType<P, PR>;
    Matrix<X> orthogonal;
    Matrix<X> quasi_triangular;
    std::vector<SizeType> block_sizes;
};

//! \brief Status and optional factors on SUCCESS.
template<class P, class PR>
struct SchurDecompositionResult {
    SchurStatus status;
    std::optional<SchurDecompositionData<P, PR>> data;
};

//! \brief Conservative first real Schur backend.
//! \details ApproximateTag: Householder Hessenberg reduction followed by
//! real Wilkinson single shifts or explicit quadratic double shifts.
//! All successful outputs pass finite, structure, orthogonality and
//! original-matrix residual checks. An implicit Francis bulge-chasing
//! implementation will replace the explicit polynomial QR in future.
//! ValidatedTag is not yet supported and always returns INCONCLUSIVE.
template<class P, class PR>
class SchurDecomposition {
  public:
    using X = FloatType<P, PR>;
    using Result = SchurDecompositionResult<P, PR>;
    explicit SchurDecomposition(PR pr) : _precision(pr) {}

    Result operator()(const Matrix<X>& A) const {
        ARIADNE_PRECONDITION(A.row_size() == A.column_size());
        const SizeType n = A.row_size();
        if (n == 0u) {
            return {SchurStatus::UNSUPPORTED_DIMENSION, std::nullopt};
        }
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return {SchurStatus::INCONCLUSIVE, std::nullopt};
        } else {
            using Raw = typename X::RawType;
            const X zero(0u, _precision);
            const X one(1u, _precision);
            const X two_value(2u, _precision);
            const X eps(Raw::eps(_precision));
            X scale(one);
            if (!this->_finite_matrix(A)) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != n; ++j) {
                    const X mag = abs(A[i][j]);
                    if (mag.raw() > scale.raw()) { scale = mag; }
                }
            }
            const X tolerance = X(128u, _precision) * eps * scale;
            const X residual_limit = X(32768u, _precision) * eps * scale *
                X(static_cast<unsigned>(n * n), _precision);
            const X orthogonal_limit = X(32768u, _precision) * eps *
                X(static_cast<unsigned>(n * n), _precision);
            if (!this->_finite(tolerance) || !this->_finite(residual_limit) ||
                !this->_finite(orthogonal_limit)) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            Matrix<X> T(A);
            Matrix<X> Q = Matrix<X>::identity(n, _precision);

            // Initial orthogonal Householder reduction to Hessenberg form.
            for (SizeType k = 0u; k + 2u < n; ++k) {
                X norm(zero);
                for (SizeType i = k + 1u; i != n; ++i) {
                    norm += sqr(T[i][k]);
                }
                if (!this->_finite(norm)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                if (norm.raw() == zero.raw()) { continue; }
                const X mag = sqrt(norm);
                const X alpha = T[k + 1u][k].raw() < zero.raw() ? mag : -mag;
                std::vector<X> v(n - k - 1u, zero);
                v[0u] = T[k + 1u][k] - alpha;
                X vnorm(sqr(v[0u]));
                for (SizeType i = k + 2u; i != n; ++i) {
                    v[i - k - 1u] = T[i][k];
                    vnorm += sqr(v[i - k - 1u]);
                }
                if (!this->_finite(vnorm) || !(vnorm.raw() > zero.raw())) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                const X beta = two_value / vnorm;
                for (SizeType j = k; j != n; ++j) {
                    X dot(zero);
                    for (SizeType i = k + 1u; i != n; ++i) {
                        dot += v[i - k - 1u] * T[i][j];
                    }
                    for (SizeType i = k + 1u; i != n; ++i) {
                        T[i][j] -= beta * dot * v[i - k - 1u];
                    }
                }
                for (SizeType i = 0u; i != n; ++i) {
                    X dot(zero);
                    for (SizeType j = k + 1u; j != n; ++j) {
                        dot += T[i][j] * v[j - k - 1u];
                    }
                    for (SizeType j = k + 1u; j != n; ++j) {
                        T[i][j] -= beta * dot * v[j - k - 1u];
                    }
                    X qdot(zero);
                    for (SizeType j = k + 1u; j != n; ++j) {
                        qdot += Q[i][j] * v[j - k - 1u];
                    }
                    for (SizeType j = k + 1u; j != n; ++j) {
                        Q[i][j] -= beta * qdot * v[j - k - 1u];
                    }
                }
            }
            if (!this->_finite_matrix(T) || !this->_finite_matrix(Q)) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }

            SizeType active = n;
            const SizeType limit = 512u * n * n;
            for (SizeType iter = 0u; active != 0u && iter != limit; ++iter) {
                if (active == 1u) { --active; continue; }
                if (abs(T[active - 1u][active - 2u]).raw() <= tolerance.raw()) {
                    T[active - 1u][active - 2u] = zero;
                    --active;
                    continue;
                }
                const X& a = T[active - 2u][active - 2u];
                const X& b = T[active - 2u][active - 1u];
                const X& c = T[active - 1u][active - 2u];
                const X& d = T[active - 1u][active - 1u];
                const X disc = sqr(a - d) + X(4u, _precision) * b * c;
                if (!this->_finite(disc)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                const bool pair = disc.raw() < zero.raw();
                if (pair && active == 2u) { active = 0u; break; }
                if (pair && active > 2u &&
                    abs(T[active - 2u][active - 3u]).raw() <= tolerance.raw()) {
                    T[active - 2u][active - 3u] = zero;
                    active -= 2u;
                    continue;
                }

                Matrix<X> polynomial(active, active, zero);
                if (pair) {
                    // Explicit double shift p(H)=H^2-tr(B)H+det(B)I;
                    // use purely real arithmetic for a conjugate pair.
                    const X trace = a + d;
                    const X det = a * d - b * c;
                    for (SizeType i = 0u; i != active; ++i) {
                        for (SizeType j = 0u; j != active; ++j) {
                            X product(zero);
                            for (SizeType k = 0u; k != active; ++k) {
                                product += T[i][k] * T[k][j];
                            }
                            polynomial[i][j] = product - trace * T[i][j] +
                                (i == j ? det : zero);
                        }
                    }
                } else {
                    const X root = sqrt(disc);
                    const X low = (a + d - root) / two_value;
                    const X high = (a + d + root) / two_value;
                    const X shift = abs(low - d).raw() < abs(high - d).raw()
                        ? low : high;
                    for (SizeType i = 0u; i != active; ++i) {
                        for (SizeType j = 0u; j != active; ++j) {
                            polynomial[i][j] = T[i][j] -
                                (i == j ? shift : zero);
                        }
                    }
                }
                const auto small_Q = this->_qr(polynomial);
                if (!small_Q) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                Matrix<X> step = Matrix<X>::identity(n, _precision);
                for (SizeType i = 0u; i != active; ++i) {
                    for (SizeType j = 0u; j != active; ++j) {
                        step[i][j] = (*small_Q)[i][j];
                    }
                }
                // Full similarity includes couplings to deflated blocks.
                T = transpose(step) * (T * step);
                Q = Q * step;
                if (!this->_finite_matrix(T) || !this->_finite_matrix(Q)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                for (SizeType col = 0u; col + 2u < n; ++col) {
                    for (SizeType row = col + 2u; row != n; ++row) {
                        if (abs(T[row][col]).raw() > tolerance.raw()) {
                            return {SchurStatus::INCONCLUSIVE, std::nullopt};
                        }
                        T[row][col] = zero;
                    }
                }
            }
            if (active != 0u) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            std::vector<SizeType> blocks;
            for (SizeType i = 0u; i < n;) {
                if (i + 1u < n &&
                    abs(T[i + 1u][i]).raw() > tolerance.raw()) {
                    if (i + 2u < n &&
                        abs(T[i + 2u][i + 1u]).raw() > tolerance.raw()) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    const X disc = sqr(T[i][i] - T[i + 1u][i + 1u]) +
                        X(4u, _precision) * T[i][i + 1u] * T[i + 1u][i];
                    if (!this->_finite(disc) || !(disc.raw() < zero.raw())) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    blocks.push_back(2u);
                    i += 2u;
                } else {
                    if (i + 1u < n) { T[i + 1u][i] = zero; }
                    blocks.push_back(1u);
                    ++i;
                }
            }

            const Matrix<X> A_Q = A * Q;
            const Matrix<X> Q_T = Q * T;
            const Matrix<X> Q_Q = transpose(Q) * Q;
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != n; ++j) {
                    const X error = abs(A_Q[i][j] - Q_T[i][j]);
                    const X orth_error = abs(Q_Q[i][j] -
                        (i == j ? one : zero));
                    if (!this->_finite(error) || !this->_finite(orth_error) ||
                        error.raw() > residual_limit.raw() ||
                        orth_error.raw() > orthogonal_limit.raw()) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                }
            }
            return {SchurStatus::SUCCESS, SchurDecompositionData<P, PR>{
                std::move(Q), std::move(T), std::move(blocks)}};
        }
    }

  private:
    bool _finite(const X& value) const {
        if constexpr (std::is_same_v<P, ApproximateTag>) {
            return is_finite(value.raw());
        } else {
            return is_finite(value.lower_raw()) && is_finite(value.upper_raw());
        }
    }

    bool _finite_matrix(const Matrix<X>& A) const {
        for (SizeType i = 0u; i != A.row_size(); ++i) {
            for (SizeType j = 0u; j != A.column_size(); ++j) {
                if (!this->_finite(A[i][j])) { return false; }
            }
        }
        return true;
    }

    // Full Householder QR, returning its orthogonal factor.
    // Zero columns in rank-deficient shift polynomials are allowed.
    std::optional<Matrix<X>> _qr(const Matrix<X>& input) const {
        const SizeType n = input.row_size();
        const X zero(0u, _precision);
        const X two_value(2u, _precision);
        Matrix<X> R(input);
        Matrix<X> Q = Matrix<X>::identity(n, _precision);
        for (SizeType k = 0u; k + 1u < n; ++k) {
            X norm(zero);
            for (SizeType i = k; i != n; ++i) {
                norm += sqr(R[i][k]);
            }
            if (!this->_finite(norm)) { return std::nullopt; }
            if (norm.raw() == zero.raw()) { continue; }
            const X mag = sqrt(norm);
            const X alpha = R[k][k].raw() < zero.raw() ? mag : -mag;
            std::vector<X> v(n - k, zero);
            v[0u] = R[k][k] - alpha;
            X vnorm(sqr(v[0u]));
            for (SizeType i = k + 1u; i != n; ++i) {
                v[i - k] = R[i][k];
                vnorm += sqr(v[i - k]);
            }
            if (!this->_finite(vnorm) || !(vnorm.raw() > zero.raw())) {
                return std::nullopt;
            }
            const X beta = two_value / vnorm;
            for (SizeType j = k; j != n; ++j) {
                X dot(zero);
                for (SizeType i = k; i != n; ++i) { dot += v[i - k] * R[i][j]; }
                for (SizeType i = k; i != n; ++i) {
                    R[i][j] -= beta * dot * v[i - k];
                }
            }
            for (SizeType i = 0u; i != n; ++i) {
                X dot(zero);
                for (SizeType j = k; j != n; ++j) { dot += Q[i][j] * v[j - k]; }
                for (SizeType j = k; j != n; ++j) {
                    Q[i][j] -= beta * dot * v[j - k];
                }
            }
        }
        if (!this->_finite_matrix(Q)) { return std::nullopt; }
        return Q;
    }

    PR _precision;
};

} // namespace Ariadne

#endif // ARIADNE_SCHUR_DECOMPOSITION_HPP
