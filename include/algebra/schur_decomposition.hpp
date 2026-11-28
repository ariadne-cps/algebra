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
 *  @brief Real Schur reduction using Householder reflections and shifted QR.
 *
 *  @details Mathematical background and historical references are documented
 *  on SchurDecomposition. The present double shift is explicit: it must not
 *  be confused with Francis's implicit double-shift algorithm.
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

//! \brief Real Schur factors approximating A*Q = Q*T.
//! \details Q is approximately orthogonal; T is real upper quasi-triangular.
//! Its diagonal 1x1 blocks hold real eigenvalues, while 2x2 blocks hold
//! complex-conjugate eigenvalue pairs (negative 2x2 discriminant).
//! A real Schur form exists even if A is defective and has no eigenbasis.
//! The factors are floating-point approximations, NOT interval certificates.
//! \see SchurDecomposition for the underlying reductions and references.
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

//! \brief Approximate a real Schur decomposition of a real square matrix.
//!
//! \details The real Schur theorem gives an orthogonal Q and a real
//! quasi-upper-triangular T such that
//! \f[ A = Q T Q^{\mathsf T}, \qquad Q^{\mathsf T}Q=I. \f]
//! Scalar diagonal blocks encode real eigenvalues; irreducible 2x2
//! diagonal blocks encode conjugate nonreal eigenvalues. Schur reduction
//! does not assume diagonalizability and does not produce eigenvectors.
//!
//! \par Step 1: Householder reflections and Hessenberg reduction
//! A Householder reflector is an orthogonal, symmetric matrix
//! \f[ U = I - 2vv^{\mathsf T}/(v^{\mathsf T}v). \f]
//! Choosing v from the portion of column k below the diagonal makes
//! the components below the first subdiagonal vanish after applying
//! the similarity T <- U*T*U, while accumulating Q <- Q*U. Repeating
//! this for k=0,...,n-3 yields an upper Hessenberg matrix H, i.e.
//! H(i,j)=0 when i>j+1. Orthogonal similarity preserves eigenvalues
//! and is preferable to a nonorthogonal triangularization because
//! it does not amplify errors through an ill-conditioned basis.
//! The reflection construction derives from Householder's 1958
//! orthogonal triangularization; its use for Hessenberg eigenvalue
//! reduction is described by Wilkinson (1965, Chapter 6).
//!
//! \par Step 2: shifted QR transformations
//! For a real shift mu, form H-mu*I=U*R by Householder QR, then set
//! \f[ H_{\mathrm{new}} = R U+\mu I = U^{\mathsf T} H U. \f]
//! The shift is chosen as the real eigenvalue of the trailing 2x2
//! block closest to its bottom-right diagonal entry (a Wilkinson-type
//! shift). When that 2x2 has a negative discriminant, its eigenvalues
//! are a conjugate pair and one instead forms the REAL quadratic
//! polynomial
//! \f[ p(H)=H^2-\operatorname{tr}(B)H+\det(B)I, \f]
//! where B is the trailing 2x2 block. Factoring p(H)=U*R and applying
//! H <- U^T*H*U is an explicit real double-shift QR step. This is
//! based on the QR iteration of Francis (1961-1962), NOT its implicit
//! bulge-chasing implementation. Forming H^2 explicitly is slower
//! and less robust than implicit Francis QR, which remains future work.
//!
//! \par Step 3: deflation and acceptance criteria
//! A sufficiently small subdiagonal entry is set to zero using a
//! precision/scale-dependent tolerance, splitting off a real 1x1
//! block or a nonreal 2x2 block. Iterations are bounded. Before
//! SUCCESS, finite factors, the 1x1/2x2 block structure, approximate
//! orthogonality and the ORIGINAL-input residual A*Q-Q*T are checked.
//! The tolerances are numerical heuristics, not proof of a backward
//! error bound in the ValidatedTag sense. ApproximateTag can return
//! INCONCLUSIVE on clusters, overflows and lack of convergence;
//! ValidatedTag currently always returns INCONCLUSIVE.
//!
//! \par Historical sources
//! - A. S. Householder, "Unitary Triangularization of a Nonsymmetric
//!   Matrix", Journal of the ACM 5(4), 339-342 (1958).
//!   https://doi.org/10.1145/320941.320947
//! - J. G. F. Francis, "The QR Transformation", Parts 1 and 2,
//!   The Computer Journal 4(3), 265-271 (1961) and 4(4), 332-345 (1962).
//!   https://doi.org/10.1093/comjnl/4.3.265
//!   https://doi.org/10.1093/comjnl/4.4.332
//! - J. H. Wilkinson, The Algebraic Eigenvalue Problem, Clarendon
//!   Press, Oxford (1965), Chapters 6 and 8 (Hessenberg and QR).
//!
//! \see Spectrum for extracting the eigenvalue multiset without an
//! eigenvector basis; EigenDecomposition for the stronger requirement
//! of a complete real diagonalization.
template<class P, class PR>
class SchurDecomposition {
  public:
    using X = FloatType<P, PR>;
    using Result = SchurDecompositionResult<P, PR>;
    //! \brief Select the floating-point precision used by the QR iteration.
    explicit SchurDecomposition(PR pr) : _precision(pr) {}

    //! \brief Attempt approximate real Schur reduction of A.
    //! \param A Real square matrix (non-square input violates precondition).
    //! \return Schur factors only on SUCCESS. Zero order is unsupported;
    //! the validated paradigm is deliberately inconclusive.
    //! \details First uses orthogonal Householder similarities to obtain
    //! Hessenberg form, then performs shifted QR with deflation and
    //! checks A*Q-Q*T and Q^T*Q-I. See the class documentation for
    //! derivation, numerical limitations and historical references.
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

            /** \brief Householder reduction to upper Hessenberg form.
             *  \details Each reflector annihilates the tail of a column
             *  by choosing its sign to avoid cancellation. Apply the
             *  reflector on both sides (orthogonal similarity), and
             *  multiply it into Q. See Householder (1958), Wilkinson (1965).
             */
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
                    /** \brief Explicit REAL quadratic QR double shift.
                     *  \details p(H)=H^2-tr(B)*H+det(B)*I for the final
                     *  2x2 block B. Its two shifts are complex conjugates,
                     *  so p has real coefficients. This construction
                     *  follows the Francis QR framework, but it is NOT
                     *  the implicit Francis bulge-chasing algorithm.
                     */
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
                    // Wilkinson-type single shift: choose the nearer
                    // real root of the trailing 2x2 characteristic
                    // polynomial. See Wilkinson (1965), Chapter 8.
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
            // A Schur partition must exhaust the matrix with disjoint
            // 1x1/2x2 blocks, without probing beyond its final row.
            if (T.row_size() != n || T.column_size() != n) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            std::vector<SizeType> blocks;
            SizeType covered = 0u;
            while (covered < n) {
                const SizeType remaining = n - covered;
                if (remaining == 1u) {
                    blocks.push_back(1u);
                    ++covered;
                    break;
                }
                const SizeType i = covered;
                const X first_subdiagonal = abs(T[i + 1u][i]);
                if (!this->_finite(first_subdiagonal)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                if (first_subdiagonal.raw() > tolerance.raw()) {
                    if (remaining > 2u) {
                        const X second_subdiagonal = abs(T[i + 2u][i + 1u]);
                        if (!this->_finite(second_subdiagonal) ||
                            second_subdiagonal.raw() > tolerance.raw()) {
                            return {SchurStatus::INCONCLUSIVE, std::nullopt};
                        }
                    }
                    const X disc = sqr(T[i][i] - T[i + 1u][i + 1u]) +
                        X(4u, _precision) * T[i][i + 1u] * T[i + 1u][i];
                    if (!this->_finite(disc) || !(disc.raw() < zero.raw())) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    blocks.push_back(2u);
                    covered += 2u;
                } else {
                    T[i + 1u][i] = zero;
                    blocks.push_back(1u);
                    ++covered;
                }
            }
            if (covered != n) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
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

    //! \brief Orthogonal factor of a full Householder QR factorization.
    //! \details Repeated reflections of the form
    //! \f$U=I-2vv^{\mathsf T}/(v^{\mathsf T}v)\f$ reduce an arbitrary
    //! shifted matrix (or real quadratic shift polynomial) to
    //! triangular form R. The accumulated Q satisfies input=Q*R
    //! approximately. A rank-deficient column with zero norm is
    //! skipped, not normalized. The reflector sign reduces cancellation.
    //! This is Householder QR, not classical Gram-Schmidt; see
    //! Householder (1958) and Wilkinson (1965), Chapter 4.
    //! \return Orthogonal Q, or no value if an intermediate value
    //! is not finite or a reflector cannot be formed.
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
