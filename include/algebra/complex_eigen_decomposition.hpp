/***************************************************************************
 *            algebra/complex_eigen_decomposition.hpp
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

/** @file algebra/complex_eigen_decomposition.hpp
 *  @brief Approximate complex eigenvector bases of real matrices.
 */
#ifndef ARIADNE_COMPLEX_EIGEN_DECOMPOSITION_HPP
#define ARIADNE_COMPLEX_EIGEN_DECOMPOSITION_HPP

#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "algebra/matrix.hpp"
#include "algebra/spectrum.hpp"
#include "numeric/complex.hpp"

namespace Ariadne {

//! \brief Dense, explicitly initialized storage of complex entries.
//! \details The generic Matrix<Complex<X>> currently requires scalar
//! characteristics metadata not provided by numeric/complex.hpp.
//! Its mutable element proxies also introduce ambiguous complex
//! operators. Here each indexed element is an actual Complex<X>&,
//! so arithmetic and assignment retain their intended scalar type.
//! This storage type is not a second QR or eigendecomposition backend.
//! It may be replaced once Matrix supports complex coefficients.
template<class X>
class ComplexMatrix {
  public:
    using ValueType = Complex<X>;
    using Row = std::vector<ValueType>;

    //! \brief Initialize every entry in an m-by-n matrix.
    ComplexMatrix(SizeType m, SizeType n, const ValueType& value)
        : _rows(m, Row(n, value)), _columns(n) {}

    //! \brief Number of rows.
    SizeType row_size() const { return _rows.size(); }
    //! \brief Number of columns.
    SizeType column_size() const { return _columns; }

    //! \brief Mutable row access with ordinary Complex<X> references.
    Row& operator[](SizeType i) { return _rows.at(i); }
    //! \brief Read-only row access.
    const Row& operator[](SizeType i) const { return _rows.at(i); }

  private:
    std::vector<Row> _rows;
    SizeType _columns;
};

//! \brief Outcome of complex diagonalization; failure is not nonexistence.
enum class ComplexEigenDecompositionStatus {
    SUCCESS,              //!< Complete complex eigenbasis, numerically checked.
    INCONCLUSIVE,         //!< A complete well-conditioned basis was not established.
    UNSUPPORTED_DIMENSION //!< Order-zero matrices are not supported.
};

//! \brief Approximate A*V = V*D over the complex numbers.
//! \details The eigenvalues (including conjugate pairs) correspond to
//! columns of V. The inverse of V is provided and checked numerically.
//! The factors are approximate: no interval/eigenvalue certification
//! follows from a small residual alone.
template<class P, class PR>
struct ComplexEigenDecompositionData {
    using X = FloatType<P, PR>;
    using CX = Complex<X>;
    std::vector<CX> eigenvalues;
    ComplexMatrix<X> eigenvectors;
    ComplexMatrix<X> inverse_eigenvectors;
};

//! \brief Status and optional diagonalization factors.
template<class P, class PR>
struct ComplexEigenDecompositionResult {
    ComplexEigenDecompositionStatus status;
    std::optional<ComplexEigenDecompositionData<P, PR>> data;
};

//! \brief Compute a complete complex eigenvector basis of a REAL matrix.
//!
//! \details This initial ApproximateTag implementation is for spectra
//! whose complex eigenvalues are numerically simple and well separated.
//! It is a distinct scalar-field result from the older real-only
//! EigenDecomposition API; there is NO second QR implementation here.
//!
//! \par Method derived from the real Schur form
//! First compute A*Q = Q*T with SchurDecomposition (Householder
//! Hessenberg reduction and Francis/Wilkinson shifted QR). A real
//! Schur diagonal block is either 1x1 with a real eigenvalue or 2x2
//! with a pair of conjugate complex eigenvalues. The eigenvalues
//! are extracted with Spectrum's 2x2 quadratic formula.
//!
//! For a scalar diagonal block T_kk=lambda, initialize w_k=1.
//! For a 2x2 block B=[a b;c d] and an eigenvalue lambda,
//! initialize the corresponding two entries using
//! \f$w_k=b,\quad w_{k+1}=\lambda-a\f$. Earlier Schur blocks
//! are solved by block back-substitution in
//! \f[(T-\lambda I)w=0.\f]
//! For an earlier 2x2 block, this requires an ordinary complex
//! 2x2 linear solve; two distinct Schur blocks with nearly
//! coincident eigenvalues are deliberately rejected rather than
//! dividing by an unreliable determinant. Lift to the original
//! coordinates by v=Q*w and scale each vector.
//!
//! \par Numerical acceptance and limitations
//! Before SUCCESS, verify all eigenvalues are separated, all complex
//! entries are finite, A*V-V*D is small for the ORIGINAL real matrix,
//! the complex eigenvector matrix V is invertible by pivoted complex
//! Gauss-Jordan elimination, and both V*V^{-1} and V^{-1}*V are close
//! to identity. These are approximate numerical checks, not proofs.
//! Repeated/defective and poorly conditioned spectra may be
//! INCONCLUSIVE even when a basis mathematically exists. ValidatedTag
//! currently returns INCONCLUSIVE: no floating Schur transform or
//! complex Gauss-Jordan inverse is a rigorous interval certificate.
//!
//! \par Mathematical sources
//! - J. H. Wilkinson, The Algebraic Eigenvalue Problem (1965),
//!   Chapters 6-8: Hessenberg, Schur/QR and eigenvectors.
//! - J. G. F. Francis, "The QR Transformation", Parts 1 and 2,
//!   The Computer Journal 4 (1961-1962).
//! - See SchurDecomposition for the historical references and Doxygen
//!   derivation of Householder reflections and implicit QR.
template<class P, class PR>
class ComplexEigenDecomposition {
  public:
    using X = FloatType<P, PR>;
    using CX = Complex<X>;
    using Result = ComplexEigenDecompositionResult<P, PR>;

    //! \brief Select the requested real coefficient precision.
    explicit ComplexEigenDecomposition(PR precision) : _precision(precision) {}

    //! \brief Attempt complex diagonalization of the real square matrix A.
    Result operator()(const Matrix<X>& A) const {
        using Status = ComplexEigenDecompositionStatus;
        ARIADNE_PRECONDITION(A.row_size() == A.column_size());
        const SizeType n = A.row_size();
        if (n == 0u) {
            return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
        }
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return {Status::INCONCLUSIVE, std::nullopt};
        } else {
            using Raw = typename X::RawType;
            const X zero(0u, _precision);
            const X one(1u, _precision);
            const CX complex_zero(zero, zero);
            const CX complex_one(one, zero);
            const X eps(Raw::eps(_precision));
            X scale(one);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    if (!this->_finite_real(A[row][col])) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    const X size = abs(A[row][col]);
                    if (size.raw() > scale.raw()) { scale = size; }
                }
            }
            const X separation = X(256u, _precision) * eps * scale;
            const X residual_limit = X(32768u, _precision) * eps * scale *
                X(static_cast<unsigned>(n * n), _precision);
            const X basis_limit = X(32768u, _precision) * eps *
                X(static_cast<unsigned>(n * n), _precision);
            const X pivot_limit = X(256u, _precision) * eps *
                X(static_cast<unsigned>(n), _precision);
            if (!this->_finite_real(residual_limit) ||
                !this->_finite_real(separation) ||
                !this->_finite_real(basis_limit) ||
                !this->_finite_real(pivot_limit)) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            const auto schur = SchurDecomposition<P, PR>(_precision)(A);
            if (schur.status != SchurStatus::SUCCESS || !schur.data) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            const Matrix<X>& T = schur.data->quasi_triangular;
            const Matrix<X>& Q = schur.data->orthogonal;
            if (T.row_size() != n || T.column_size() != n ||
                Q.row_size() != n || Q.column_size() != n) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            struct Block {
                SizeType start;
                SizeType width;
            };
            struct Root {
                SizeType block;
                CX eigenvalue;
            };
            std::vector<Block> blocks;
            std::vector<Root> roots;
            SizeType offset = 0u;
            for (const SizeType width: schur.data->block_sizes) {
                if ((width != 1u && width != 2u) ||
                    offset >= n || width > n - offset) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                const SizeType id = blocks.size();
                blocks.push_back({offset, width});
                if (width == 1u) {
                    roots.push_back({id, CX(T[offset][offset], zero)});
                } else {
                    Matrix<X> subblock(2u, 2u, zero);
                    for (SizeType row = 0u; row != 2u; ++row) {
                        for (SizeType col = 0u; col != 2u; ++col) {
                            subblock[row][col] = T[offset + row][offset + col];
                        }
                    }
                    const auto pair = Spectrum<P, PR>(_precision)(subblock);
                    if (pair.status != SpectralStatus::SUCCESS || !pair.data ||
                        pair.data->eigenvalues.size() != 2u) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    for (const CX& root: pair.data->eigenvalues) {
                        roots.push_back({id, root});
                    }
                }
                offset += width;
            }
            if (offset != n || roots.size() != n) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            // A simple-spectrum contract makes every earlier Schur
            // block nonsingular at each target eigenvalue. Do not
            // confuse an unresolved cluster with a valid eigenbasis.
            const X min_gap_squared = sqr(separation);
            for (SizeType i = 0u; i != n; ++i) {
                if (!this->_finite(roots[i].eigenvalue)) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                for (SizeType j = 0u; j != i; ++j) {
                    const X distance_squared = this->_norm2(
                        roots[i].eigenvalue - roots[j].eigenvalue);
                    if (!this->_finite_real(distance_squared) ||
                        distance_squared.raw() <= min_gap_squared.raw()) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            }

            ComplexMatrix<X> V(n, n, complex_zero);
            for (SizeType col = 0u; col != n; ++col) {
                const CX& lambda = roots[col].eigenvalue;
                const SizeType id = roots[col].block;
                const Block& target = blocks[id];
                std::vector<CX> w(n, complex_zero);
                if (target.width == 1u) {
                    w[target.start] = complex_one;
                } else {
                    const SizeType k = target.start;
                    // The vector (b,lambda-a) solves the first
                    // row of (B-lambda I)w=0 identically. Since
                    // the discriminant is negative, b is nonzero.
                    w[k] = CX(T[k][k + 1u], zero);
                    w[k + 1u] = lambda - CX(T[k][k], zero);
                }

                for (SizeType before = id; before > 0u;) {
                    --before;
                    const Block& block = blocks[before];
                    const SizeType k = block.start;
                    const SizeType end = k + block.width;
                    CX rhs0(complex_zero);
                    CX rhs1(complex_zero);
                    for (SizeType j = end; j != n; ++j) {
                        rhs0 += CX(T[k][j], zero) * w[j];
                        if (block.width == 2u) {
                            rhs1 += CX(T[k + 1u][j], zero) * w[j];
                        }
                    }
                    if (block.width == 1u) {
                        const CX divisor = lambda - CX(T[k][k], zero);
                        const X size = this->_norm2(divisor);
                        if (!this->_finite_real(size) ||
                            size.raw() <= min_gap_squared.raw()) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                        w[k] = rhs0 / divisor;
                    } else {
                        const CX lambda_a = lambda - CX(T[k][k], zero);
                        const CX lambda_d = lambda - CX(T[k + 1u][k + 1u], zero);
                        const CX b(T[k][k + 1u], zero);
                        const CX c(T[k + 1u][k], zero);
                        const CX determinant = lambda_a * lambda_d - b * c;
                        const X size = this->_norm2(determinant);
                        // A 2x2 determinant has eigenvalue^2 units.
                        const X threshold = sqr(min_gap_squared);
                        if (!this->_finite_real(size) ||
                            size.raw() <= threshold.raw()) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                        w[k] = (lambda_d * rhs0 + b * rhs1) / determinant;
                        w[k + 1u] = (c * rhs0 + lambda_a * rhs1) / determinant;
                    }
                }

                // Lift the Schur eigenvector through the REAL
                // orthogonal similarity into the original coordinates.
                X largest(zero);
                std::vector<CX> lifted(n, complex_zero);
                for (SizeType row = 0u; row != n; ++row) {
                    for (SizeType k = 0u; k != n; ++k) {
                        lifted[row] += CX(Q[row][k], zero) * w[k];
                    }
                    const X size = this->_norm2(lifted[row]);
                    if (!this->_finite_real(size)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    if (size.raw() > largest.raw()) { largest = size; }
                }
                if (!(largest.raw() > zero.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                const X magnitude = sqrt(largest);
                for (SizeType row = 0u; row != n; ++row) {
                    V[row][col] = lifted[row] / CX(magnitude, zero);
                    if (!this->_finite(V[row][col])) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            }

            const X residual_squared = sqr(residual_limit);
            for (SizeType col = 0u; col != n; ++col) {
                for (SizeType row = 0u; row != n; ++row) {
                    CX action(complex_zero);
                    for (SizeType k = 0u; k != n; ++k) {
                        action += CX(A[row][k], zero) * V[k][col];
                    }
                    const CX residual = action - roots[col].eigenvalue * V[row][col];
                    const X size = this->_norm2(residual);
                    if (!this->_finite_real(size) ||
                        size.raw() > residual_squared.raw()) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            }

            auto inverse = this->_inverse(V, pivot_limit);
            if (!inverse) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            const X basis_squared = sqr(basis_limit);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    CX first(complex_zero);
                    CX second(complex_zero);
                    for (SizeType k = 0u; k != n; ++k) {
                        first += V[row][k] * (*inverse)[k][col];
                        second += (*inverse)[row][k] * V[k][col];
                    }
                    const CX identity = row == col ? complex_one : complex_zero;
                    const X first_error = this->_norm2(first - identity);
                    const X second_error = this->_norm2(second - identity);
                    if (!this->_finite_real(first_error) ||
                        !this->_finite_real(second_error) ||
                        first_error.raw() > basis_squared.raw() ||
                        second_error.raw() > basis_squared.raw()) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            }

            std::vector<CX> eigenvalues;
            eigenvalues.reserve(n);
            for (const Root& root: roots) {
                eigenvalues.push_back(root.eigenvalue);
            }
            return {Status::SUCCESS, ComplexEigenDecompositionData<P, PR>{
                std::move(eigenvalues), std::move(V), std::move(*inverse)}};
        }
    }

  private:
    bool _finite_real(const X& value) const {
        if constexpr (std::is_same_v<P, ApproximateTag>) {
            return is_finite(value.raw());
        } else {
            return is_finite(value.lower_raw()) && is_finite(value.upper_raw());
        }
    }

    bool _finite(const CX& value) const {
        return this->_finite_real(value.real_part()) &&
               this->_finite_real(value.imaginary_part());
    }

    X _norm2(const CX& value) const {
        return sqr(value.real_part()) + sqr(value.imaginary_part());
    }

    //! \brief Compute and check an approximate complex matrix inverse.
    //! \details Gauss-Jordan elimination with maximal-modulus-squared
    //! row pivoting implements invertibility testing without relying
    //! on real ordered comparisons for complex numbers. A pivot
    //! below the precision-relative threshold is inconclusive.
    //! The caller also verifies both inverse residuals.
    std::optional<ComplexMatrix<X>> _inverse(const ComplexMatrix<X>& V,
                                       const X& pivot_limit) const {
        const SizeType n = V.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
        const CX complex_zero(zero, zero);
        const CX complex_one(one, zero);
        ComplexMatrix<X> work(V);
        ComplexMatrix<X> inverse(n, n, complex_zero);
        for (SizeType i = 0u; i != n; ++i) {
            inverse[i][i] = complex_one;
        }
        const X threshold = sqr(pivot_limit);
        for (SizeType col = 0u; col != n; ++col) {
            SizeType pivot = col;
            X largest(zero);
            for (SizeType row = col; row != n; ++row) {
                const X size = this->_norm2(work[row][col]);
                if (!this->_finite_real(size)) { return std::nullopt; }
                if (size.raw() > largest.raw()) {
                    largest = size;
                    pivot = row;
                }
            }
            if (largest.raw() <= threshold.raw()) { return std::nullopt; }
            if (pivot != col) {
                for (SizeType j = 0u; j != n; ++j) {
                    std::swap(work[pivot][j], work[col][j]);
                    std::swap(inverse[pivot][j], inverse[col][j]);
                }
            }
            const CX diagonal = work[col][col];
            for (SizeType j = 0u; j != n; ++j) {
                work[col][j] /= diagonal;
                inverse[col][j] /= diagonal;
            }
            for (SizeType row = 0u; row != n; ++row) {
                if (row == col) { continue; }
                const CX factor = work[row][col];
                for (SizeType j = 0u; j != n; ++j) {
                    work[row][j] -= factor * work[col][j];
                    inverse[row][j] -= factor * inverse[col][j];
                }
            }
        }
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                if (!this->_finite(work[row][col]) ||
                    !this->_finite(inverse[row][col])) {
                    return std::nullopt;
                }
            }
        }
        return inverse;
    }

    PR _precision;
};

} // namespace Ariadne

#endif // ARIADNE_COMPLEX_EIGEN_DECOMPOSITION_HPP
