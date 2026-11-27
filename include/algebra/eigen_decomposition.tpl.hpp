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
#include <exception>

#include "algebra/validated_root_contractors.hpp"

namespace Ariadne {

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
    return this->operator()(M, nullptr);
}

template<class P, class PR>
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::operator()(
        const Matrix<NumericType>& M,
        EigenDecompositionDiagnostics* diagnostics) const {
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
                return this->_nonsymmetric(M);
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


/** @brief Certify dense interval eigenpairs using separated Gershgorin disks. */
template<class P, class PR>
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::_gershgorin(const Matrix<NumericType>& input) const {
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
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::_preconditioned(
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
        const ApproximateX half = ApproximateX(1u, _precision) /
                                  ApproximateX(2u, _precision);

        // Candidate basis only: use a symmetric matrix built from interval
        // midpoints. Even if the original family contains nonsymmetric
        // matrices, the later interval certificate applies to every member.
        Matrix<ApproximateX> seed(n, n, ApproximateX(0u, _precision));
        for (SizeType row = 0u; row != n; ++row) {
            seed[row][row] = ApproximateX(A[row][row]);
            for (SizeType col = row + 1u; col != n; ++col) {
                const ApproximateX center = (ApproximateX(A[row][col]) +
                                             ApproximateX(A[col][row])) * half;
                seed[row][col] = center;
                seed[col][row] = center;
            }
        }

        EigenDecomposition<ApproximateTag, PR> approximate(_precision);
        const auto candidate = approximate(seed);
        if (candidate.status != Status::SUCCESS || !candidate.data) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        candidate_timer.stop();
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
    }
}


/** @brief Verify normalized real eigenpairs by interval Krawczyk inclusion. */
template<class P, class PR>
typename EigenDecomposition<P, PR>::ResultType
EigenDecomposition<P, PR>::_krawczyk(
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
        const AX half = approximate_one / AX(2u, _precision);

        // Jacobi supplies only starting guesses. All subsequent proof
        // obligations use outward-rounded interval arithmetic.
        Matrix<AX> midpoint(n, n, approximate_zero);
        Matrix<AX> seed(n, n, approximate_zero);
        for (SizeType row = 0u; row != n; ++row) {
            for (SizeType col = 0u; col != n; ++col) {
                midpoint[row][col] = AX(A[row][col]);
                if (!is_finite(midpoint[row][col].raw())) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        for (SizeType row = 0u; row != n; ++row) {
            seed[row][row] = midpoint[row][row];
            for (SizeType col = row + 1u; col != n; ++col) {
                const AX symmetric_entry =
                    (midpoint[row][col] + midpoint[col][row]) * half;
                seed[row][col] = symmetric_entry;
                seed[col][row] = symmetric_entry;
            }
        }

        const EigenDecomposition<ApproximateTag, PR> approximate(_precision);
        const auto candidates = approximate(seed);
        if (candidates.status != Status::SUCCESS || !candidates.data) {
            return {Status::INCONCLUSIVE, std::nullopt};
        }
        candidate_timer.stop();

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
    }
}

} // namespace Ariadne

#include "algebra/eigen_nonsymmetric.tpl.hpp"

#endif // ARIADNE_EIGEN_DECOMPOSITION_TPL_HPP
