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
 *  @brief Template implementations of real and complex eigendecomposition.
 */

#ifndef ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP
#define ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP

#include <algorithm>
#include <vector>
#include <type_traits>
#include <utility>
#include <exception>

#include "algebra/validated_root_contractors.hpp"

namespace Ariadne {
namespace detail {

// Complex diagonalization is an implementation detail of the unified solver.
// Keep its algorithm in this template implementation header.

template<class P, class PR>
auto ComplexEigenDecomposition<P, PR>::operator()(const Matrix<X>& A) const -> Result {
        using Status = ComplexEigenDecompositionStatus;
        ARIADNE_PRECONDITION(A.row_size() == A.column_size());
        if (A.row_size() == 0u) {
            return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
        }
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return {Status::INCONCLUSIVE, std::nullopt};
        } else {
            const auto schur = SchurDecomposition<P, PR>(_precision)(A);
            if (schur.status != SchurStatus::SUCCESS || !schur.data) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            return (*this)(A, *schur.data);
        }
    }

template<class P, class PR>
auto ComplexEigenDecomposition<P, PR>::operator()(const Matrix<X>& A,
                      const SchurDecompositionData<P, PR>& schur_data) const -> Result {
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

            const Matrix<X>& T = schur_data.quasi_triangular;
            const Matrix<X>& Q = schur_data.orthogonal;
            if (T.row_size() != n || T.column_size() != n ||
                Q.row_size() != n || Q.column_size() != n) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            struct Block {
                SizeType start;
                SizeType width;
            };
            struct SpectralRootEntry {
                SizeType block;
                CX eigenvalue;
            };
            std::vector<Block> blocks;
            std::vector<SpectralRootEntry> roots;
            SizeType offset = 0u;
            for (const SizeType width: schur_data.block_sizes) {
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
            for (const SpectralRootEntry& root: roots) {
                eigenvalues.push_back(root.eigenvalue);
            }
            return {Status::SUCCESS, ComplexEigenDecompositionData<P, PR>{
                std::move(eigenvalues), std::move(V), std::move(*inverse)}};
        }
    }

template<class P, class PR>
bool ComplexEigenDecomposition<P, PR>::_finite_real(const X& value) const {
        if constexpr (std::is_same_v<P, ApproximateTag>) {
            return is_finite(value.raw());
        } else {
            return is_finite(value.lower_raw()) && is_finite(value.upper_raw());
        }
    }

template<class P, class PR>
bool ComplexEigenDecomposition<P, PR>::_finite(const CX& value) const {
        return this->_finite_real(value.real_part()) &&
               this->_finite_real(value.imaginary_part());
    }

template<class P, class PR>
auto ComplexEigenDecomposition<P, PR>::_norm2(const CX& value) const -> X {
        return sqr(value.real_part()) + sqr(value.imaginary_part());
    }

template<class P, class PR>
auto ComplexEigenDecomposition<P, PR>::_inverse(const ComplexMatrix<X>& V,
                                       const X& pivot_limit) const -> std::optional<ComplexMatrix<X>> {
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


//! \brief Account for one optional, non-overlapping algorithm stage.
class EigenDecompositionSectionTimer {
  public:
    explicit EigenDecompositionSectionTimer(std::chrono::nanoseconds* counter)
        : _counter(counter), _start(counter ? Clock::now() : Clock::time_point{}) {}
    ~EigenDecompositionSectionTimer() { this->stop(); }
    void stop() {
        if (_counter) {
            *_counter += std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - _start);
            _counter = nullptr;
        }
    }
    EigenDecompositionSectionTimer(const EigenDecompositionSectionTimer&) = delete;
    EigenDecompositionSectionTimer& operator=(const EigenDecompositionSectionTimer&) = delete;

  private:
    using Clock = std::chrono::steady_clock;
    std::chrono::nanoseconds* _counter;
    Clock::time_point _start;
};

template<class P, class PR>
bool RealEigenDecomposition<P, PR>::_positive(const NumericType& x) const {
    if constexpr (std::is_same_v<P, ValidatedTag>) {
        return x.lower_raw() > NumericType(0u, _precision).lower_raw();
    } else {
        return is_finite(x.raw()) && x.raw() > NumericType(0u, _precision).raw();
    }
}

template<class P, class PR>
bool RealEigenDecomposition<P, PR>::_negative(const NumericType& x) const {
    if constexpr (std::is_same_v<P, ValidatedTag>) {
        return x.upper_raw() < NumericType(0u, _precision).upper_raw();
    } else {
        return is_finite(x.raw()) && x.raw() < NumericType(0u, _precision).raw();
    }
}

template<class P, class PR>
bool RealEigenDecomposition<P, PR>::_zero(const NumericType& x) const {
    if constexpr (std::is_same_v<P, ValidatedTag>) {
        const auto zero = NumericType(0u, _precision).lower_raw();
        return x.lower_raw() == zero && x.upper_raw() == zero;
    } else {
        return is_finite(x.raw()) && x.raw() == NumericType(0u, _precision).raw();
    }
}

template<class P, class PR>
bool RealEigenDecomposition<P, PR>::_nonzero(const NumericType& x) const {
    return this->_positive(x) || this->_negative(x);
}

template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::operator()(const Matrix<NumericType>& M) const {
    return this->operator()(M, nullptr);
}

template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::operator()(
        const Matrix<NumericType>& M,
        EigenDecompositionDiagnostics* diagnostics) const {
    return this->operator()(M, diagnostics, nullptr);
}

template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::operator()(
        const Matrix<NumericType>& M,
        EigenDecompositionDiagnostics* diagnostics,
        std::optional<SchurDecompositionData<P, PR>>*
            complex_candidate) const {
    using Status = EigenDecompositionStatus;
    if (complex_candidate) {
        complex_candidate->reset();
    }
    auto result = this->_attempt(M, diagnostics, complex_candidate);
    if (result.status != Status::SUCCESS) {
        return result;
    }
    // Every algorithm (including recursive structural paths) passes
    // through this gate. Finite inputs alone cannot prevent arithmetic
    // overflow while lifting or normalizing an eigenvector.
    if (!result.data ||
        result.data->eigenvalues.size() != M.row_size() ||
        result.data->eigenvectors.row_size() != M.row_size() ||
        result.data->eigenvectors.column_size() != M.column_size()) {
        return {Status::INCONCLUSIVE, std::nullopt};
    }
    const auto finite = [](const NumericType& value) {
        if constexpr (std::is_same_v<P, ValidatedTag>) {
            return is_finite(value.lower_raw()) && is_finite(value.upper_raw());
        } else {
            return is_finite(value.raw());
        }
    };
    const SizeType n = M.row_size();
    for (SizeType col = 0u; col != n; ++col) {
        if (!finite(result.data->eigenvalues[col])) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        for (SizeType row = 0u; row != n; ++row) {
            if (!finite(result.data->eigenvectors[row][col])) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
        }
    }
    return result;
}

template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::_attempt(
        const Matrix<NumericType>& M,
        EigenDecompositionDiagnostics* diagnostics,
        std::optional<SchurDecompositionData<P, PR>>*
            complex_candidate) const {
    if (diagnostics) {
        *diagnostics = {};
    }
    using X = NumericType;
    using Status = EigenDecompositionStatus;
    ARIADNE_PRECONDITION(M.row_size() == M.column_size());

    const SizeType n = M.row_size();
    if (n == 0u) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    }
    // Check all input entries before any structural shortcut, including
    // the 1x1 and block-triangular paths. Only finite real eigenpairs
    // belong to the result contract.
    for (SizeType row = 0u; row != n; ++row) {
        for (SizeType col = 0u; col != n; ++col) {
            if constexpr (std::is_same_v<P, ValidatedTag>) {
                if (!is_finite(M[row][col].lower_raw()) ||
                    !is_finite(M[row][col].upper_raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            } else if (!is_finite(M[row][col].raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
        }
    }
    const X zero(0u, _precision);
    const X one(1u, _precision);
    Vector<X> values(n, zero);
    Matrix<X> vectors(n, n, zero);

    // Every real diagonal matrix has the canonical coordinate eigenbasis,
    // including matrices with repeated eigenvalues. In ValidatedTag the
    // diagonal entries represent *independent* interval coordinates:
    // an ascending assignment is sound only for strictly separated
    // intervals, or equal singleton values known to coincide exactly.
    if (n > 1u) {
        bool diagonal = true;
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                if (row != col && !this->_zero(M[row][col])) {
                    diagonal = false;
                }
            }
        }
        if (diagonal) {
            std::vector<SizeType> order;
            order.reserve(n);
            for (SizeType row = 0u; row != n; ++row) {
                const X& current = M[row][row];
                if constexpr (std::is_same_v<P, ValidatedTag>) {
                    if (!is_finite(current.lower_raw()) ||
                        !is_finite(current.upper_raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                } else {
                    if (!is_finite(current.raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
                SizeType position = 0u;
                for (; position != order.size(); ++position) {
                    const X& previous = M[order[position]][order[position]];
                    if constexpr (std::is_same_v<P, ValidatedTag>) {
                        if (current.upper_raw() < previous.lower_raw()) {
                            break;
                        }
                        const bool separated_after =
                            current.lower_raw() > previous.upper_raw();
                        const bool repeated_point =
                            current.lower_raw() == current.upper_raw() &&
                            previous.lower_raw() == previous.upper_raw() &&
                            current.lower_raw() == previous.lower_raw();
                        if (!separated_after && !repeated_point) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                    } else {
                        if (current.raw() < previous.raw()) {
                            break;
                        }
                    }
                }
                order.insert(order.begin() + static_cast<std::ptrdiff_t>(position), row);
            }
            for (SizeType col = 0u; col != n; ++col) {
                const SizeType pivot = order[col];
                values[col] = M[pivot][pivot];
                vectors[pivot][col] = one;
            }
            if (diagnostics) {
                diagnostics->strategy = n == 2u
                    ? EigenDecompositionStrategy::ANALYTIC
                    : EigenDecompositionStrategy::TRIANGULAR;
            }
            return {Status::SUCCESS, EigenDecompositionData<P, PR>{
                std::move(values), std::move(vectors)}};
        }
    }

    // Detect a block-diagonal structure up to a simultaneous row and
    // column permutation. An off-block entry must be *exactly zero*;
    // intervals containing zero but also nonzero values connect blocks.
    // The undirected nonzero pattern gives the finest independent blocks.
    if (n > 2u) {
        std::vector<bool> reached(n, false);
        std::vector<std::vector<SizeType>> components;
        for (SizeType start = 0u; start != n; ++start) {
            if (reached[start]) {
                continue;
            }
            std::vector<SizeType> component{start};
            reached[start] = true;
            for (SizeType index = 0u; index != component.size(); ++index) {
                const SizeType row = component[index];
                for (SizeType col = 0u; col != n; ++col) {
                    if (!reached[col] &&
                        (!this->_zero(M[row][col]) || !this->_zero(M[col][row]))) {
                        reached[col] = true;
                        component.push_back(col);
                    }
                }
            }
            std::sort(component.begin(), component.end());
            components.push_back(std::move(component));
        }
        if (components.size() > 1u) {
            if (diagnostics) {
                diagnostics->strategy = EigenDecompositionStrategy::BLOCK_DIAGONAL;
            }
            SizeType offset = 0u;
            for (const auto& component: components) {
                const SizeType m = component.size();
                Matrix<X> block(m, m, zero);
                for (SizeType row = 0u; row != m; ++row) {
                    for (SizeType col = 0u; col != m; ++col) {
                        block[row][col] = M[component[row]][component[col]];
                    }
                }
                EigenDecompositionDiagnostics subdiagnostics;
                const auto subresult = this->operator()(
                    block, diagnostics ? &subdiagnostics : nullptr);
                if (diagnostics) {
                    diagnostics->candidate_time += subdiagnostics.candidate_time;
                    diagnostics->verification_time += subdiagnostics.verification_time;
                }
                if (subresult.status != Status::SUCCESS || !subresult.data) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                for (SizeType col = 0u; col != m; ++col) {
                    values[offset + col] = subresult.data->eigenvalues[col];
                    for (SizeType row = 0u; row != m; ++row) {
                        vectors[component[row]][offset + col] =
                            subresult.data->eigenvectors[row][col];
                    }
                }
                offset += m;
            }

            // Sort all block eigenpairs without silently imposing a fixed
            // ordering on independently varying interval coefficients.
            std::vector<SizeType> order;
            order.reserve(n);
            for (SizeType col = 0u; col != n; ++col) {
                const X& current_value = values[col];
                if constexpr (std::is_same_v<P, ValidatedTag>) {
                    if (!is_finite(current_value.lower_raw()) ||
                        !is_finite(current_value.upper_raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                } else if (!is_finite(current_value.raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                SizeType position = 0u;
                for (; position != order.size(); ++position) {
                    const X& prior = values[order[position]];
                    if constexpr (std::is_same_v<P, ValidatedTag>) {
                        if (current_value.upper_raw() < prior.lower_raw()) {
                            break;
                        }
                        const bool strictly_after =
                            current_value.lower_raw() > prior.upper_raw();
                        const bool repeated_point =
                            current_value.lower_raw() == current_value.upper_raw() &&
                            prior.lower_raw() == prior.upper_raw() &&
                            current_value.lower_raw() == prior.lower_raw();
                        if (!strictly_after && !repeated_point) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                    } else if (current_value.raw() < prior.raw()) {
                        break;
                    }
                }
                order.insert(order.begin() +
                             static_cast<std::ptrdiff_t>(position), col);
            }
            Vector<X> sorted_values(n, zero);
            Matrix<X> sorted_vectors(n, n, zero);
            for (SizeType col = 0u; col != n; ++col) {
                const SizeType source = order[col];
                sorted_values[col] = values[source];
                for (SizeType row = 0u; row != n; ++row) {
                    sorted_vectors[row][col] = vectors[row][source];
                }
            }

            // Every block contributes an invertible eigenvector basis.
            // The global basis is block diagonal up to permutations and
            // thus invertible, even for eigenvalues repeated across blocks.
            return {Status::SUCCESS, EigenDecompositionData<P, PR>{
                std::move(sorted_values), std::move(sorted_vectors)}};
        }
    }

    // A directed cycle can be confined to a diagonal block even when
    // the entire matrix is not permutationally triangular. The strongly
    // connected components of the structural nonzero graph give an upper
    // block-triangular form without treating uncertain zeros as absent.
    if (n > 2u) {
        std::vector<std::vector<bool>> reachable(
            n, std::vector<bool>(n, false));
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                reachable[row][col] =
                    row == col || !this->_zero(M[row][col]);
            }
        }
        for (SizeType via = 0u; via != n; ++via) {
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    if (reachable[row][via] && reachable[via][col]) {
                        reachable[row][col] = true;
                    }
                }
            }
        }

        std::vector<SizeType> component_of(n, n);
        std::vector<std::vector<SizeType>> components;
        for (SizeType start = 0u; start != n; ++start) {
            if (component_of[start] != n) {
                continue;
            }
            const SizeType id = components.size();
            std::vector<SizeType> component;
            for (SizeType vertex = 0u; vertex != n; ++vertex) {
                if (component_of[vertex] == n &&
                    reachable[start][vertex] && reachable[vertex][start]) {
                    component_of[vertex] = id;
                    component.push_back(vertex);
                }
            }
            components.push_back(std::move(component));
        }

        // Acyclic singleton components are already handled by the
        // permutation-triangular solver below. Restrict this path to
        // genuinely cyclic diagonal blocks, with at least two blocks.
        bool cyclic_block = false;
        for (const auto& component: components) {
            cyclic_block = cyclic_block || component.size() > 1u;
        }
        if (components.size() > 1u && cyclic_block) {
            const SizeType count = components.size();
            std::vector<std::vector<bool>> edges(
                count, std::vector<bool>(count, false));
            std::vector<SizeType> indegree(count, 0u);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    const SizeType from = component_of[row];
                    const SizeType to = component_of[col];
                    if (from != to && !this->_zero(M[row][col]) &&
                        !edges[from][to]) {
                        edges[from][to] = true;
                        ++indegree[to];
                    }
                }
            }

            std::vector<SizeType> order;
            std::vector<bool> placed(count, false);
            for (SizeType k = 0u; k != count; ++k) {
                SizeType next = count;
                for (SizeType id = 0u; id != count; ++id) {
                    if (!placed[id] && indegree[id] == 0u) {
                        next = id;
                        break;
                    }
                }
                if (next == count) {
                    break; // The condensation graph must be acyclic.
                }
                placed[next] = true;
                order.push_back(next);
                for (SizeType id = 0u; id != count; ++id) {
                    if (edges[next][id]) {
                        --indegree[id];
                    }
                }
            }

            if (order.size() == count) {
                if (diagnostics) {
                    diagnostics->strategy = EigenDecompositionStrategy::BLOCK_TRIANGULAR;
                }
                bool viable = true;
                std::vector<ResultType> solved_blocks;
                for (SizeType k = 0u; k != count && viable; ++k) {
                    const auto& indices = components[order[k]];
                    const SizeType m = indices.size();
                    Matrix<X> diagonal_block(m, m, zero);
                    for (SizeType row = 0u; row != m; ++row) {
                        for (SizeType col = 0u; col != m; ++col) {
                            diagonal_block[row][col] =
                                M[indices[row]][indices[col]];
                        }
                    }
                    EigenDecompositionDiagnostics nested_diagnostics;
                    auto answer = this->operator()(
                        diagonal_block, diagnostics ? &nested_diagnostics : nullptr);
                    if (diagnostics) {
                        diagnostics->candidate_time += nested_diagnostics.candidate_time;
                        diagnostics->verification_time += nested_diagnostics.verification_time;
                    }
                    if (answer.status != Status::SUCCESS || !answer.data) {
                        viable = false;
                    } else {
                        solved_blocks.push_back(std::move(answer));
                    }
                }

                if (viable) {
                    EigenDecompositionSectionTimer lifting_time(
                        diagnostics ? (std::is_same_v<P, ValidatedTag>
                            ? &diagnostics->verification_time
                            : &diagnostics->candidate_time) : nullptr);
                    SizeType output_col = 0u;
                    for (SizeType k = 0u; k != count && viable; ++k) {
                        const auto& indices = components[order[k]];
                        const auto& data = *solved_blocks[k].data;
                        for (SizeType eigen_col = 0u;
                             eigen_col != indices.size() && viable; ++eigen_col) {
                            const X lambda = data.eigenvalues[eigen_col];
                            Vector<X> lifted(n, zero);
                            for (SizeType row = 0u; row != indices.size(); ++row) {
                                lifted[indices[row]] =
                                    data.eigenvectors[row][eigen_col];
                            }

                            // For an upper block-triangular matrix, an
                            // eigenvector of block k is lifted into blocks
                            // preceding k via (B_i-lambda*I)v_i = -sum A_ij*v_j.
                            for (SizeType previous = k; previous > 0u && viable;) {
                                --previous;
                                const auto& earlier = components[order[previous]];
                                const SizeType m = earlier.size();
                                Matrix<X> shifted(m, m, zero);
                                Vector<X> rhs(m, zero);
                                for (SizeType row = 0u; row != m; ++row) {
                                    const SizeType original = earlier[row];
                                    for (SizeType col = 0u; col != m; ++col) {
                                        shifted[row][col] =
                                            M[original][earlier[col]];
                                    }
                                    shifted[row][row] -= lambda;
                                    for (SizeType col = 0u; col != n; ++col) {
                                        // All components in the current block
                                        // are still zero; only later blocks
                                        // contribute to the right-hand side.
                                        rhs[row] += M[original][col] * lifted[col];
                                    }
                                }

                                try {
                                    const Matrix<X> shifted_inverse = inverse(shifted);
                                    for (SizeType row = 0u; row != m; ++row) {
                                        X entry(zero);
                                        for (SizeType col = 0u; col != m; ++col) {
                                            entry -= shifted_inverse[row][col] * rhs[col];
                                        }
                                        lifted[earlier[row]] = entry;
                                    }
                                } catch (const SingularMatrixException&) {
                                    // A singular shifted block may still
                                    // admit a solution when the coupling
                                    // has zero projection onto every
                                    // eigenspace with eigenvalue lambda.
                                    // Use its already certified eigenbasis,
                                    // not a pseudoinverse of an interval
                                    // matrix. Unproved compatibility is
                                    // inconclusive and falls through.
                                    const auto& earlier_data =
                                        *solved_blocks[previous].data;
                                    const auto inverse_basis =
                                        earlier_data.inverse_eigenvector_matrix();
                                    if (!inverse_basis) {
                                        viable = false;
                                    } else {
                                        Vector<X> coefficients(m, zero);
                                        for (SizeType eig = 0u; eig != m; ++eig) {
                                            X projection(zero);
                                            for (SizeType row = 0u; row != m; ++row) {
                                                projection +=
                                                    (*inverse_basis)[eig][row] * rhs[row];
                                            }
                                            const X difference =
                                                earlier_data.eigenvalues[eig] - lambda;
                                            if (this->_nonzero(difference)) {
                                                coefficients[eig] = -projection / difference;
                                            } else if (this->_zero(difference) &&
                                                       this->_zero(projection)) {
                                                // Compatible repeated exact
                                                // eigenvalue: set this free
                                                // eigenvector coordinate to 0.
                                                coefficients[eig] = zero;
                                            } else {
                                                viable = false;
                                                break;
                                            }
                                        }
                                        if (viable) {
                                            for (SizeType row = 0u; row != m; ++row) {
                                                X entry(zero);
                                                for (SizeType eig = 0u; eig != m; ++eig) {
                                                    entry += earlier_data.eigenvectors[row][eig] *
                                                        coefficients[eig];
                                                }
                                                lifted[earlier[row]] = entry;
                                            }
                                        }
                                    }
                                }
                            }
                            if (viable) {
                                values[output_col] = lambda;
                                for (SizeType row = 0u; row != n; ++row) {
                                    vectors[row][output_col] = lifted[row];
                                }
                                ++output_col;
                            }
                        }
                    }
                }

                if (viable) {
                    // Preserve the documented ascending order for *every*
                    // member of an interval family. Overlapping independent
                    // eigenvalue ranges cannot be assigned fixed columns.
                    std::vector<SizeType> eigen_order;
                    for (SizeType col = 0u; col != n && viable; ++col) {
                        const X& current = values[col];
                        if constexpr (std::is_same_v<P, ValidatedTag>) {
                            if (!is_finite(current.lower_raw()) ||
                                !is_finite(current.upper_raw())) {
                                viable = false;
                                break;
                            }
                        } else if (!is_finite(current.raw())) {
                            viable = false;
                            break;
                        }
                        SizeType position = 0u;
                        for (; position != eigen_order.size(); ++position) {
                            const X& previous = values[eigen_order[position]];
                            if constexpr (std::is_same_v<P, ValidatedTag>) {
                                if (current.upper_raw() < previous.lower_raw()) {
                                    break;
                                }
                                const bool strictly_after =
                                    current.lower_raw() > previous.upper_raw();
                                const bool same_point =
                                    current.lower_raw() == current.upper_raw() &&
                                    previous.lower_raw() == previous.upper_raw() &&
                                    current.lower_raw() == previous.lower_raw();
                                if (!strictly_after && !same_point) {
                                    viable = false;
                                    break;
                                }
                            } else if (current.raw() < previous.raw()) {
                                break;
                            }
                        }
                        if (viable) {
                            eigen_order.insert(eigen_order.begin() +
                                static_cast<std::ptrdiff_t>(position), col);
                        }
                    }

                    if (viable) {
                        Vector<X> sorted_values(n, zero);
                        Matrix<X> sorted_vectors(n, n, zero);
                        for (SizeType col = 0u; col != n; ++col) {
                            const SizeType source = eigen_order[col];
                            sorted_values[col] = values[source];
                            for (SizeType row = 0u; row != n; ++row) {
                                sorted_vectors[row][col] = vectors[row][source];
                            }
                        }
                        // Before lifting, each diagonal block supplies an
                        // invertible basis. Lifting only adds components in
                        // earlier blocks; hence the full basis is invertible.
                        return {Status::SUCCESS, EigenDecompositionData<P, PR>{
                            std::move(sorted_values), std::move(sorted_vectors)}};
                    }
                }
                // An inconclusive block or shifted inverse must not block
                // the existing generic dense eigendecomposition fallback.
            }
        }
    }

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
            // A connected directed acyclic off-diagonal pattern is upper
            // triangular after a simultaneous coordinate permutation.
            // Only exactly zero entries are absent from this graph:
            // an uncertain interval containing zero is still an edge.
            std::vector<SizeType> indegree(n, 0u);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    if (row != col && !this->_zero(M[row][col])) {
                        ++indegree[col];
                    }
                }
            }
            std::vector<bool> placed(n, false);
            std::vector<SizeType> permutation;
            permutation.reserve(n);
            for (SizeType k = 0u; k != n; ++k) {
                SizeType next = n;
                for (SizeType row = 0u; row != n; ++row) {
                    if (!placed[row] && indegree[row] == 0u) {
                        next = row;
                        break;
                    }
                }
                if (next == n) {
                    // A directed cycle prevents triangularization by a
                    // permutation. Preserve the generic dense fallback.
                    break;
                }
                placed[next] = true;
                permutation.push_back(next);
                for (SizeType col = 0u; col != n; ++col) {
                    if (!placed[col] && !this->_zero(M[next][col])) {
                        --indegree[col];
                    }
                }
            }
            if (permutation.size() == n) {
                Matrix<X> permuted(n, n, zero);
                for (SizeType row = 0u; row != n; ++row) {
                    for (SizeType col = 0u; col != n; ++col) {
                        permuted[row][col] =
                            M[permutation[row]][permutation[col]];
                    }
                }
                EigenDecompositionDiagnostics nested_diagnostics;
                const auto triangular = this->operator()(
                    permuted, diagnostics ? &nested_diagnostics : nullptr);
                if (diagnostics) {
                    diagnostics->candidate_time += nested_diagnostics.candidate_time;
                    diagnostics->verification_time += nested_diagnostics.verification_time;
                }
                if (triangular.status == Status::SUCCESS && triangular.data) {
                    // Undo the row permutation on the eigenvectors. The
                    // eigenvalues are already sorted by the triangular
                    // solver; columns must retain that ordering.
                    for (SizeType col = 0u; col != n; ++col) {
                        values[col] = triangular.data->eigenvalues[col];
                        for (SizeType row = 0u; row != n; ++row) {
                            vectors[permutation[row]][col] =
                                triangular.data->eigenvectors[row][col];
                        }
                    }
                    if (diagnostics) {
                        diagnostics->strategy = EigenDecompositionStrategy::TRIANGULAR;
                    }
                    return {Status::SUCCESS, EigenDecompositionData<P, PR>{
                        std::move(values), std::move(vectors)}};
                }
                // Inconclusive triangular substitution is not a proof that
                // another solver cannot handle the matrix; fall through.
            }
            if constexpr (std::is_same_v<P, ApproximateTag>) {
                if (diagnostics) {
                    diagnostics->strategy = EigenDecompositionStrategy::JACOBI;
                }
                EigenDecompositionSectionTimer candidate(
                    diagnostics ? &diagnostics->candidate_time : nullptr);
                const auto symmetric = this->_symmetric(M);
                if (symmetric.status != Status::UNSUPPORTED_DIMENSION) {
                    return symmetric;
                }
                if (diagnostics) {
                    diagnostics->strategy = EigenDecompositionStrategy::SHIFTED_QR;
                }
                return this->_nonsymmetric(M, complex_candidate);
            } else {
                if (diagnostics) {
                    diagnostics->strategy = EigenDecompositionStrategy::GERSHGORIN;
                }
                ResultType direct;
                {
                    EigenDecompositionSectionTimer verified(
                        diagnostics ? &diagnostics->verification_time : nullptr);
                    direct = this->_gershgorin(M);
                }
                if (direct.status == Status::SUCCESS) {
                    return direct;
                }
                if (diagnostics) {
                    diagnostics->strategy = EigenDecompositionStrategy::SIMILARITY;
                }
                const auto preconditioned = this->_preconditioned(M, diagnostics);
                if (preconditioned.status == Status::SUCCESS) {
                    return preconditioned;
                }
                if (diagnostics) {
                    diagnostics->strategy = EigenDecompositionStrategy::KRAWCZYK;
                }
                return this->_krawczyk(M, diagnostics);
            }
        }

        // An interval-independent ordering is necessary to identify the
        // eigenpair belonging to each output column across the entire
        // represented family of real matrices.
        std::vector<SizeType> order;
        order.reserve(n);
        for (SizeType row = 0u; row != n; ++row) {
            const X& current_value = M[row][row];
            if constexpr (std::is_same_v<P, ValidatedTag>) {
                if (!is_finite(current_value.lower_raw()) ||
                    !is_finite(current_value.upper_raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            } else if (!is_finite(current_value.raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            SizeType position = 0u;
            for (; position != order.size(); ++position) {
                const X& prior_value = M[order[position]][order[position]];
                if constexpr (std::is_same_v<P, ValidatedTag>) {
                    if (current_value.upper_raw() < prior_value.lower_raw()) {
                        break;
                    }
                    const bool ordered_after =
                        current_value.lower_raw() > prior_value.upper_raw();
                    const bool exactly_repeated =
                        current_value.lower_raw() == current_value.upper_raw() &&
                        prior_value.lower_raw() == prior_value.upper_raw() &&
                        current_value.lower_raw() == prior_value.lower_raw();
                    // A repeated interval description does not identify a
                    // common parameter: independent diagonal entries may
                    // represent different eigenvalues for different members.
                    if (!ordered_after && !exactly_repeated) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                } else {
                    if (current_value.raw() < prior_value.raw()) {
                        break;
                    }
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
                    if (this->_nonzero(divisor)) {
                        vectors[row][col] = -numerator / divisor;
                    } else if (this->_zero(divisor) && this->_zero(numerator)) {
                        // Equal diagonal eigenvalues are permitted when this
                        // dependent row equation is identically satisfied.
                        // Assign zero to the free component, retaining the
                        // exact unit-pivot triangular eigenvector basis.
                        vectors[row][col] = zero;
                    } else {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            } else {
                // Lower triangular analogue: forward substitution.
                for (SizeType row = pivot + 1u; row != n; ++row) {
                    X numerator(zero);
                    for (SizeType k = pivot; k < row; ++k) {
                        numerator += M[row][k] * vectors[k][col];
                    }
                    const X divisor = M[row][row] - eigenvalue;
                    if (this->_nonzero(divisor)) {
                        vectors[row][col] = -numerator / divisor;
                    } else if (this->_zero(divisor) && this->_zero(numerator)) {
                        // Equal diagonal eigenvalues are permitted when this
                        // dependent row equation is identically satisfied.
                        // Assign zero to the free component, retaining the
                        // exact unit-pivot triangular eigenvector basis.
                        vectors[row][col] = zero;
                    } else {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
            }
        }

        // Even with repeated diagonal eigenvalues, the computed eigenvector
        // matrix is a column permutation of a unit triangular matrix.
        // Its exact pivots establish invertibility for every interval member.
        // Repeated-eigenvalue equations that cannot be shown identically
        // compatible have already returned INCONCLUSIVE.
        if (diagnostics) {
            diagnostics->strategy = EigenDecompositionStrategy::TRIANGULAR;
        }
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(values), std::move(vectors)}};
    }

    if (diagnostics) {
        diagnostics->strategy = EigenDecompositionStrategy::ANALYTIC;
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

    // Exactly diagonal matrices (also with repeated eigenvalues) have
    // already been handled by the dimension-independent coordinate basis.
    if (this->_nonzero(b)) {
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
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::_symmetric(const Matrix<NumericType>& input) const {
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

        // Sort the eigenpairs. Jacobi maintains an orthogonal eigenbasis
        // even when eigenvalues repeat; simplicity is not required.
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
        // Orthogonality also certifies numerical linear independence when
        // a symmetric eigenspace has multiplicity >1. It is checked in the
        // input arithmetic, without collapsing MP values to binary64.
        const X orthogonality_budget =
            X(4096u, _precision) * epsilon *
            X(static_cast<unsigned>(n), _precision);
        for (SizeType col = 0u; col != n; ++col) {
            for (SizeType other = 0u; other != n; ++other) {
                X dot(zero);
                for (SizeType row = 0u; row != n; ++row) {
                    dot += eigenvectors[row][col] * eigenvectors[row][other];
                }
                const X difference = abs(dot - (col == other ? one : zero));
                if (!is_finite(difference.raw()) ||
                    difference.raw() > orthogonality_budget.raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(eigenvalues), std::move(eigenvectors)}};
    }
}


/** @brief Certify dense interval eigenpairs using separated Gershgorin disks. */
template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::_gershgorin(const Matrix<NumericType>& input) const {
    using Status = EigenDecompositionStatus;
    if constexpr (!std::is_same_v<P, ValidatedTag>) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    } else {
        using X = NumericType;
        const SizeType n = input.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
        Vector<X> disks(n, zero);

        // Every entry of the input interval matrix varies independently.
        // The real projection of Gershgorin disk i is contained in
        // A_ii + [-sum_{j != i} |A_ij|, +sum_{j != i} |A_ij|].
        for (SizeType row = 0u; row != n; ++row) {
            X radius(zero);
            for (SizeType col = 0u; col != n; ++col) {
                if (!is_finite(input[row][col].lower_raw()) ||
                    !is_finite(input[row][col].upper_raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                if (col != row) {
                    radius += abs(input[row][col]);
                }
            }
            if (!is_finite(radius.upper_raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            disks[row] = input[row][row] +
                         X(-radius.upper_raw(), radius.upper_raw());
        }

        // Strictly separated *real projections* imply separated complex
        // Gershgorin disks for every matrix in the input family. Each
        // disk contains one eigenvalue with algebraic multiplicity one.
        // It must be real: nonreal conjugate eigenvalues would both
        // belong to the same isolated disk of a real matrix.
        std::vector<SizeType> order;
        order.reserve(n);
        for (SizeType row = 0u; row != n; ++row) {
            SizeType position = 0u;
            for (; position != order.size(); ++position) {
                const SizeType previous = order[position];
                if (disks[row].upper_raw() < disks[previous].lower_raw()) {
                    break;
                }
                if (!(disks[row].lower_raw() > disks[previous].upper_raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
            order.insert(order.begin() + static_cast<std::ptrdiff_t>(position), row);
        }

        // For the eigenvalue in disk i, the principal submatrix obtained by
        // deleting row and column i has no eigenvalue in that disk. Thus the
        // ith eigenvector component is nonzero and may be fixed to exactly
        // one. Solve the remaining equations with interval Gaussian
        // elimination; if an interval pivot contains zero, do not certify.
        Matrix<X> unsorted(n, n, zero);
        const SizeType m = n - 1u;
        for (SizeType pivot = 0u; pivot != n; ++pivot) {
            unsorted[pivot][pivot] = one;
            const X eigenvalue = disks[pivot];
            Matrix<X> system(m, m, zero);
            Vector<X> right(m, zero);

            for (SizeType r = 0u; r != m; ++r) {
                const SizeType original_row = r < pivot ? r : r + 1u;
                right[r] = input[original_row][pivot];
                for (SizeType c = 0u; c != m; ++c) {
                    const SizeType original_col = c < pivot ? c : c + 1u;
                    system[r][c] = (r == c ? eigenvalue : zero) -
                                   input[original_row][original_col];
                }
            }

            for (SizeType k = 0u; k != m; ++k) {
                const X diagonal = system[k][k];
                if (!this->_nonzero(diagonal)) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                for (SizeType r = k + 1u; r != m; ++r) {
                    const X multiplier = system[r][k] / diagonal;
                    system[r][k] = zero;
                    for (SizeType c = k + 1u; c != m; ++c) {
                        system[r][c] -= multiplier * system[k][c];
                    }
                    right[r] -= multiplier * right[k];
                }
            }

            Vector<X> solution(m, zero);
            for (SizeType r = m; r > 0u;) {
                --r;
                X numerator = right[r];
                for (SizeType c = r + 1u; c != m; ++c) {
                    numerator -= system[r][c] * solution[c];
                }
                const X diagonal = system[r][r];
                if (!this->_nonzero(diagonal)) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                solution[r] = numerator / diagonal;
            }
            for (SizeType r = 0u; r != m; ++r) {
                const SizeType original_row = r < pivot ? r : r + 1u;
                unsorted[original_row][pivot] = solution[r];
            }
        }

        // V=I+E before column ordering. Strict ||E||_1 < 1 certifies
        // that every possible V is invertible (Neumann series).
        for (SizeType col = 0u; col != n; ++col) {
            X off_diagonal_sum(zero);
            for (SizeType row = 0u; row != n; ++row) {
                if (row != col) {
                    off_diagonal_sum += abs(unsorted[row][col]);
                }
            }
            if (!(off_diagonal_sum.upper_raw() < one.lower_raw())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
        }

        Vector<X> values(n, zero);
        Matrix<X> vectors(n, n, zero);
        for (SizeType col = 0u; col != n; ++col) {
            const SizeType source = order[col];
            values[col] = disks[source];
            for (SizeType row = 0u; row != n; ++row) {
                vectors[row][col] = unsorted[row][source];
            }
        }
        return {Status::SUCCESS, EigenDecompositionData<P, PR>{std::move(values), std::move(vectors)}};
    }
}


/** @brief Certify a dense interval matrix by a validated similarity transform. */
template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::_preconditioned(
        const Matrix<NumericType>& A,
        EigenDecompositionDiagnostics* diagnostics) const {
    using Status = EigenDecompositionStatus;
    if constexpr (!std::is_same_v<P, ValidatedTag>) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    } else {
        using X = NumericType;
        using ApproximateX = FloatApproximation<PR>;
        EigenDecompositionSectionTimer candidate_timer(
            diagnostics ? &diagnostics->candidate_time : nullptr);
        const SizeType n = A.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
        // Try actual-midpoint eigenvectors as a fixed *candidate* basis.
        // Non-symmetric midpoints use shifted QR; symmetric ones use Jacobi.
        // Nothing is accepted until S^{-1} and the Gershgorin eigenpairs
        // have been verified using the ORIGINAL interval family.
        const ApproximateX azero(0u, _precision);
        Matrix<ApproximateX> midpoint(n, n, azero);
        bool symmetric_midpoint = true;
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                midpoint[row][col] = ApproximateX(A[row][col]);
                if (!is_finite(midpoint[row][col].raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = row + 1u; col != n; ++col) {
                if (midpoint[row][col].raw() != midpoint[col][row].raw()) {
                    symmetric_midpoint = false;
                }
            }
        }
        RealEigenDecomposition<ApproximateTag, PR> approximate(_precision);
        const auto midpoint_candidate = approximate(midpoint);
        candidate_timer.stop();

        // All certificate operations are independent of the approximate
        // algorithm that proposed S; failures are inconclusive and may
        // be retried with another seed.
        const auto certify_similarity = [&](const auto& candidate) -> ResultType {
            if (candidate.status != Status::SUCCESS || !candidate.data) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            EigenDecompositionSectionTimer verified_timer(
                diagnostics ? &diagnostics->verification_time : nullptr);
            // S is a point matrix: its floating-point coefficients are treated
            // as exact dyadic constants when constructing the interval matrix.
            // Nothing relies on the computed Jacobi columns being *exactly*
            // orthogonal. We explicitly certify a valid inverse of S.
            Matrix<X> similarity(n, n, zero);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    const auto& entry = candidate.data->eigenvectors[row][col];
                    if (!is_finite(entry.raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    similarity[row][col] = X(entry.raw());
                }
            }

            // Interval Gauss-Jordan reduction encloses S^{-1} since every
            // pivot is verified nonzero and S is a fixed point matrix.
            Matrix<X> working(similarity);
            Matrix<X> inverse = Matrix<X>::identity(n, _precision);
            for (SizeType pivot = 0u; pivot != n; ++pivot) {
                SizeType chosen = pivot;
                while (chosen != n && !this->_nonzero(working[chosen][pivot])) {
                    ++chosen;
                }
                if (chosen == n) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                if (chosen != pivot) {
                    for (SizeType col = 0u; col != n; ++col) {
                        const X original = working[pivot][col];
                        working[pivot][col] = working[chosen][col];
                        working[chosen][col] = original;
                        const X inverse_original = inverse[pivot][col];
                        inverse[pivot][col] = inverse[chosen][col];
                        inverse[chosen][col] = inverse_original;
                    }
                }

                const X divisor = working[pivot][pivot];
                if (!this->_nonzero(divisor)) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                for (SizeType col = pivot + 1u; col != n; ++col) {
                    working[pivot][col] /= divisor;
                }
                working[pivot][pivot] = one;
                for (SizeType col = 0u; col != n; ++col) {
                    inverse[pivot][col] /= divisor;
                }

                for (SizeType row = 0u; row != n; ++row) {
                    if (row == pivot) {
                        continue;
                    }
                    const X multiplier = working[row][pivot];
                    for (SizeType col = pivot + 1u; col != n; ++col) {
                        working[row][col] -= multiplier * working[pivot][col];
                    }
                    working[row][pivot] = zero;
                    for (SizeType col = 0u; col != n; ++col) {
                        inverse[row][col] -= multiplier * inverse[pivot][col];
                    }
                }
            }

            // The exact similarity B = S^{-1} A S lies inside the interval
            // evaluation below for every input matrix A. Gershgorin certifies
            // simple real eigenvalues and a basis W of B when it succeeds.
            const Matrix<X> transformed = inverse * A * similarity;
            auto transformed_result = this->_gershgorin(transformed);
            if (transformed_result.status != Status::SUCCESS ||
                !transformed_result.data) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }

            // V=S W is the corresponding eigenvector basis of A.
            // Both S and W are certified invertible; the columns do not
            // necessarily retain the unit-pivot normalization of W.
            Matrix<X> vectors = similarity * transformed_result.data->eigenvectors;
            return {Status::SUCCESS,
                    EigenDecompositionData<P, PR>{
                        std::move(transformed_result.data->eigenvalues),
                        std::move(vectors)}};
        };

        if (midpoint_candidate.status == Status::SUCCESS && midpoint_candidate.data) {
            const auto direct_result = certify_similarity(midpoint_candidate);
            if (direct_result.status == Status::SUCCESS) {
                return direct_result;
            }
        }
        if (symmetric_midpoint) {
            // Jacobi has already supplied the symmetric midpoint basis.
            return {Status::INCONCLUSIVE, std::nullopt};
        }

        // Retain the old symmetrized-Jacobi route if QR fails or its
        // fixed-basis Gershgorin certificate is too conservative.
        EigenDecompositionSectionTimer fallback_timer(
            diagnostics ? &diagnostics->candidate_time : nullptr);
        const ApproximateX half = ApproximateX(1u, _precision) /
                                  ApproximateX(2u, _precision);
        Matrix<ApproximateX> symmetric_seed(n, n, azero);
        for (SizeType row = 0u; row != n; ++row) {
            symmetric_seed[row][row] = midpoint[row][row];
            for (SizeType col = row + 1u; col != n; ++col) {
                const ApproximateX center =
                    (midpoint[row][col] + midpoint[col][row]) * half;
                symmetric_seed[row][col] = center;
                symmetric_seed[col][row] = center;
            }
        }
        const auto symmetric_candidate = approximate(symmetric_seed);
        fallback_timer.stop();
        return certify_similarity(symmetric_candidate);
    }
}


/** @brief Verify normalized real eigenpairs by interval Krawczyk inclusion. */
template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::_krawczyk(
        const Matrix<NumericType>& A,
        EigenDecompositionDiagnostics* diagnostics) const {
    using Status = EigenDecompositionStatus;
    if constexpr (!std::is_same_v<P, ValidatedTag>) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    } else {
        using X = NumericType;
        using AX = FloatApproximation<PR>;
        using Raw = typename AX::RawType;
        EigenDecompositionSectionTimer candidate_timer(
            diagnostics ? &diagnostics->candidate_time : nullptr);
        const SizeType n = A.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
        const AX approximate_zero(0u, _precision);
        const AX approximate_one(1u, _precision);

        // QR on the actual midpoint provides real nonsymmetric eigenpairs
        // when available. This is only candidate generation: every proof
        // still uses outward-rounded operations on the ORIGINAL intervals.
        Matrix<AX> midpoint(n, n, approximate_zero);
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                midpoint[row][col] = AX(A[row][col]);
                if (!is_finite(midpoint[row][col].raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        bool symmetric_midpoint = true;
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = row + 1u; col != n; ++col) {
                if (midpoint[row][col].raw() != midpoint[col][row].raw()) {
                    symmetric_midpoint = false;
                }
            }
        }
        const RealEigenDecomposition<ApproximateTag, PR> approximate(_precision);
        const auto midpoint_candidates = approximate(midpoint);
        candidate_timer.stop();

        // One self-contained attempt. Approximate candidates and Newton
        // iterates determine only the trial boxes: acceptance requires the
        // validated Krawczyk predicate for every real eigenpair.
        const auto certify_candidates = [&](const auto& candidates) -> ResultType {
            if (candidates.status != Status::SUCCESS || !candidates.data) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            Vector<X> unordered_values(n, zero);
            Matrix<X> unordered_vectors(n, n, zero);
            const X eps(Raw::eps(_precision));

            for (SizeType col = 0u; col != n; ++col) {
                EigenDecompositionSectionTimer newton_timer(
                    diagnostics ? &diagnostics->candidate_time : nullptr);
                SizeType pivot = 0u;
                for (SizeType row = 1u; row != n; ++row) {
                    if (abs(candidates.data->eigenvectors[row][col]).raw() >
                        abs(candidates.data->eigenvectors[pivot][col]).raw()) {
                        pivot = row;
                    }
                }
                const AX pivot_value = candidates.data->eigenvectors[pivot][col];
                if (!is_finite(pivot_value.raw()) ||
                    pivot_value.raw() == approximate_zero.raw()) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }

                std::vector<SizeType> other_rows;
                other_rows.reserve(n - 1u);
                for (SizeType row = 0u; row != n; ++row) {
                    if (row != pivot) {
                        other_rows.push_back(row);
                    }
                }

                // Unknown z = (lambda, v[other_rows]); v[pivot]=1. Hence there
                // are exactly n equations in n unknowns:
                // F_i(z)=(A*v)_i-lambda*v_i.
                Vector<AX> z(n, approximate_zero);
                z[0u] = candidates.data->eigenvalues[col];
                for (SizeType k = 0u; k != n - 1u; ++k) {
                    z[k + 1u] = candidates.data->eigenvectors[other_rows[k]][col] /
                                 pivot_value;
                }

                const auto approximate_system =
                    [&](const Vector<AX>& at, Vector<AX>& residual,
                        Matrix<AX>& jacobian) {
                    Vector<AX> v(n, approximate_zero);
                    v[pivot] = approximate_one;
                    for (SizeType k = 0u; k != n - 1u; ++k) {
                        v[other_rows[k]] = at[k + 1u];
                    }
                    for (SizeType row = 0u; row != n; ++row) {
                        AX result(approximate_zero);
                        for (SizeType j = 0u; j != n; ++j) {
                            result += midpoint[row][j] * v[j];
                        }
                        residual[row] = result - at[0u] * v[row];
                        jacobian[row][0u] = -v[row];
                        for (SizeType k = 0u; k != n - 1u; ++k) {
                            const SizeType original_col = other_rows[k];
                            jacobian[row][k + 1u] = midpoint[row][original_col] -
                                (row == original_col ? at[0u] : approximate_zero);
                        }
                    }
                };

                // Newton steps improve the candidate but carry no certificate.
                // The interval Krawczyk check below decides whether to accept.
                for (SizeType step = 0u; step != 12u; ++step) {
                    Vector<AX> residual(n, approximate_zero);
                    Matrix<AX> jacobian(n, n, approximate_zero);
                    approximate_system(z, residual, jacobian);
                    Matrix<AX> point_inverse(n, n, approximate_zero);
                    try {
                        point_inverse = inverse(jacobian);
                    } catch (const std::exception&) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    const Vector<AX> correction = point_inverse * residual;
                    for (SizeType k = 0u; k != n; ++k) {
                        z[k] -= correction[k];
                        if (!is_finite(z[k].raw())) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                    }
                }

                Vector<AX> approximate_residual(n, approximate_zero);
                Matrix<AX> approximate_jacobian(n, n, approximate_zero);
                approximate_system(z, approximate_residual, approximate_jacobian);
                Matrix<AX> point_inverse(n, n, approximate_zero);
                try {
                    point_inverse = inverse(approximate_jacobian);
                } catch (const std::exception&) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }

                newton_timer.stop();
                EigenDecompositionSectionTimer verified_timer(
                    diagnostics ? &diagnostics->verification_time : nullptr);
                Vector<X> center(n, zero);
                Matrix<X> preconditioner(n, n, zero);
                for (SizeType k = 0u; k != n; ++k) {
                    if (!is_finite(z[k].raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    center[k] = X(z[k].raw());
                    for (SizeType j = 0u; j != n; ++j) {
                        if (!is_finite(point_inverse[k][j].raw())) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                        preconditioner[k][j] = X(point_inverse[k][j].raw());
                    }
                }

                Vector<X> exact_vector(n, zero);
                exact_vector[pivot] = one;
                for (SizeType k = 0u; k != n - 1u; ++k) {
                    exact_vector[other_rows[k]] = center[k + 1u];
                }
                Vector<X> residual(n, zero);
                for (SizeType row = 0u; row != n; ++row) {
                    X result(zero);
                    for (SizeType j = 0u; j != n; ++j) {
                        result += A[row][j] * exact_vector[j];
                    }
                    residual[row] = result - center[0u] * exact_vector[row];
                }

                // An enclosure of C*F(center), plus precision-dependent slack,
                // seeds an interval box wide enough to contain a root if one
                // can be verified. The chosen radius is NOT a proof.
                const Vector<X> predicted_correction = preconditioner * residual;
                Vector<X> base_radii(n, zero);
                for (SizeType k = 0u; k != n; ++k) {
                    const X magnitude = abs(center[k]);
                    const X coordinate_scale =
                        magnitude.upper_raw() > one.upper_raw() ? magnitude : one;
                    base_radii[k] = X(4u, _precision) * abs(predicted_correction[k]) +
                                    X(4096u, _precision) * eps * coordinate_scale;
                    if (!is_finite(base_radii[k].upper_raw()) ||
                        !(base_radii[k].upper_raw() > zero.upper_raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }

                bool certified = false;
                Vector<X> certified_box(n, zero);
                for (SizeType attempt = 0u; attempt != 4u; ++attempt) {
                    Vector<X> box(n, zero);
                    for (SizeType k = 0u; k != n; ++k) {
                        const X radius = X(1u << attempt, _precision) * base_radii[k];
                        if (!is_finite(radius.upper_raw())) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                        box[k] = center[k] +
                                 X(-radius.upper_raw(), radius.upper_raw());
                    }

                    // Analytic interval Jacobian valid for every A in the
                    // input family and every candidate (lambda,v) in box.
                    Vector<X> interval_vector(n, zero);
                    interval_vector[pivot] = one;
                    for (SizeType k = 0u; k != n - 1u; ++k) {
                        interval_vector[other_rows[k]] = box[k + 1u];
                    }
                    Matrix<X> jacobian(n, n, zero);
                    for (SizeType row = 0u; row != n; ++row) {
                        jacobian[row][0u] = -interval_vector[row];
                        for (SizeType k = 0u; k != n - 1u; ++k) {
                            const SizeType original_col = other_rows[k];
                            jacobian[row][k + 1u] = A[row][original_col] -
                                (row == original_col ? box[0u] : zero);
                        }
                    }

                    if (krawczyk_certifies_unique_root(
                            center, residual, jacobian, box, preconditioner)) {
                        certified_box = std::move(box);
                        certified = true;
                        break;
                    }
                }
                if (!certified) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
                unordered_values[col] = certified_box[0u];
                unordered_vectors[pivot][col] = one;
                for (SizeType k = 0u; k != n - 1u; ++k) {
                    unordered_vectors[other_rows[k]][col] = certified_box[k + 1u];
                }
            }

            EigenDecompositionSectionTimer ordering_timer(
                diagnostics ? &diagnostics->verification_time : nullptr);
            // Each verified box contains a real normalized eigenpair for
            // every input matrix. Pairwise disjoint real eigenvalue boxes
            // imply exactly n simple real eigenvalues and an invertible
            // eigenvector basis, pointwise for every represented matrix.
            std::vector<SizeType> order;
            order.reserve(n);
            for (SizeType col = 0u; col != n; ++col) {
                SizeType position = 0u;
                for (; position != order.size(); ++position) {
                    if (unordered_values[col].upper_raw() <
                        unordered_values[order[position]].lower_raw()) {
                        break;
                    }
                    if (!(unordered_values[col].lower_raw() >
                          unordered_values[order[position]].upper_raw())) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                }
                order.insert(order.begin() + static_cast<std::ptrdiff_t>(position), col);
            }

            Vector<X> values(n, zero);
            Matrix<X> vectors(n, n, zero);
            for (SizeType col = 0u; col != n; ++col) {
                const SizeType original = order[col];
                values[col] = unordered_values[original];
                for (SizeType row = 0u; row != n; ++row) {
                    vectors[row][col] = unordered_vectors[row][original];
                }
            }
            return {Status::SUCCESS,
                    EigenDecompositionData<P, PR>{std::move(values), std::move(vectors)}};
        };

        if (midpoint_candidates.status == Status::SUCCESS && midpoint_candidates.data) {
            const auto direct_result = certify_candidates(midpoint_candidates);
            if (direct_result.status == Status::SUCCESS) {
                return direct_result;
            }
        }
        // When the midpoint is symmetric, Jacobi has already been tried:
        // repeating the same candidates cannot improve certification.
        if (symmetric_midpoint) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }

        // Preserve the original symmetrized-Jacobi seed as a genuine
        // fallback: QR may fail or its refined eigenpair boxes may be
        // inconclusive even though another candidate basis can certify.
        EigenDecompositionSectionTimer fallback_candidate_timer(
            diagnostics ? &diagnostics->candidate_time : nullptr);
        const AX half = approximate_one / AX(2u, _precision);
        Matrix<AX> symmetric_seed(n, n, approximate_zero);
        for (SizeType row = 0u; row != n; ++row) {
            symmetric_seed[row][row] = midpoint[row][row];
            for (SizeType col = row + 1u; col != n; ++col) {
                const AX value = (midpoint[row][col] + midpoint[col][row]) * half;
                symmetric_seed[row][col] = value;
                symmetric_seed[col][row] = value;
            }
        }
        const auto symmetric_candidates = approximate(symmetric_seed);
        fallback_candidate_timer.stop();
        return certify_candidates(symmetric_candidates);
    }
}

/** @brief Compute a real nonsymmetric eigendecomposition when numerically resolvable.
 *  @details First delegate Hessenberg reduction and the Francis/Wilkinson
 *  QR iteration to SchurDecomposition; reject nonreal Schur 2x2 blocks
 *  because this API requires a complete REAL eigenvector basis.
 *  For an all-scalar real Schur partition, back-substitute triangular
 *  eigenvectors, lift them by the orthogonal Schur factor and apply
 *  Newton refinement to (A-lambda I)v=0 with v[pivot]=1 on the ORIGINAL
 *  matrix. Require separated eigenvalues, small original residuals,
 *  and an invertible eigenvector matrix before SUCCESS. No duplicated
 *  QR iteration or implicit validated certificate exists here.
 *  See SchurDecomposition Doxygen for Householder (1958), Francis
 *  (1961-1962) and Wilkinson (1965), with mathematical derivations.
 */
template<class P, class PR>
typename RealEigenDecomposition<P, PR>::ResultType
RealEigenDecomposition<P, PR>::_nonsymmetric(
        const Matrix<NumericType>& input,
        std::optional<SchurDecompositionData<P, PR>>*
            complex_candidate) const {
    using Status = EigenDecompositionStatus;
    if constexpr (!std::is_same_v<P, ApproximateTag>) {
        return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
    } else {
        using X = NumericType;
        using Raw = typename X::RawType;
        const SizeType n = input.row_size();
        const X zero(0u, _precision);
        const X one(1u, _precision);
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

        const X separation = X(256u, _precision) * epsilon * scale;
        const X residual_budget = X(8192u, _precision) * epsilon * scale *
                                  X(static_cast<unsigned>(n), _precision);

        // The unique QR eigensolver now lives in SchurDecomposition.
        // Do not maintain a second full-matrix real QR iteration here.
        auto schur = SchurDecomposition<P, PR>(_precision)(input);
        if (schur.status != SchurStatus::SUCCESS || !schur.data) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        const auto& blocks = schur.data->block_sizes;
        SizeType block_order = 0u;
        for (const SizeType width: blocks) {
            if (width != 1u) {
                // Reuse the nonreal Schur factors when the public
                // facade must fall back to complex eigenvectors.
                // Move only on failure; real successes are unaffected.
                if (width == 2u && complex_candidate) {
                    *complex_candidate = std::move(*schur.data);
                }
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            ++block_order;
        }
        if (block_order != n) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        const Matrix<X>& triangular = schur.data->quasi_triangular;
        const Matrix<X>& orthogonal = schur.data->orthogonal;

        // Schur provides a real triangular seed. Back-substitution,
        // original-matrix Newton correction, real eigenpair residual
        // checks and verified numerical basis invertibility are retained
        // below, with unchanged tolerances and failure semantics.

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

} // namespace detail
} // namespace Ariadne

#endif // ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP
