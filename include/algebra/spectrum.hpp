/***************************************************************************
 *            algebra/spectrum.hpp
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

/** @file algebra/spectrum.hpp
 *  @brief Complete real-matrix spectra in complex coordinates where established.
 */
#ifndef ARIADNE_SPECTRUM_HPP
#define ARIADNE_SPECTRUM_HPP

#include <optional>
#include <type_traits>
#include <vector>
#include <utility>

#include "algebra/matrix.hpp"
#include "numeric/complex.hpp"

namespace Ariadne {

//! \brief Whether a complete spectrum was obtained.
enum class SpectralStatus {
    SUCCESS, //!< All eigenvalues, counted with algebraic multiplicity.
    INCONCLUSIVE, //!< No complete spectrum established, not a nonexistence proof.
    UNSUPPORTED_DIMENSION //!< Zero-order matrix.
};

//! \brief A complete complex spectrum of a real matrix, including repetitions.
//! \details No eigenvector basis is implied. Ordering is implementation-
//! dependent; triangular inputs retain coordinate order. In ValidatedTag,
//! for each matrix in the input interval family, a labeling of its entire
//! spectrum exists within the corresponding complex rectangular enclosures.
template<class P, class PR>
struct SpectrumData {
    using RealType = FloatType<P, PR>;
    std::vector<Complex<RealType>> eigenvalues;
};

//! \brief Status plus optional spectral data (only on success).
template<class P, class PR>
struct SpectrumResult {
    SpectralStatus status;
    std::optional<SpectrumData<P, PR>> data;
};

//! \brief Calculate eigenvalues without requiring diagonalizability.
//! \details Handles scalar, 2x2 and structurally block-upper-triangular
//! matrices whose strongly connected diagonal blocks have order at most two.
//! General dense higher-order strongly connected blocks await Schur.
//! ValidatedTag requires a discriminant with verified sign or exact zero.
template<class P, class PR>
class Spectrum {
  public:
    using X = FloatType<P, PR>;
    using ResultType = SpectrumResult<P, PR>;
    explicit Spectrum(PR precision) : _precision(precision) {}

    ResultType operator()(const Matrix<X>& A) const {
        using Status = SpectralStatus;
        ARIADNE_PRECONDITION(A.row_size() == A.column_size());
        const SizeType n = A.row_size();
        if (n == 0u) {
            return {Status::UNSUPPORTED_DIMENSION, std::nullopt};
        }
        const X zero(0u, _precision);
        const X two_value(2u, _precision);
        for (SizeType i = 0u; i != n; ++i) {
            for (SizeType j = 0u; j != n; ++j) {
                if (!this->_finite(A[i][j])) {
                    return {Status::INCONCLUSIVE, std::nullopt};
                }
            }
        }
        SpectrumData<P, PR> result;
        if (n == 1u) {
            result.eigenvalues.emplace_back(A[0u][0u], zero);
        } else if (n == 2u) {
            const X& a = A[0u][0u];
            const X& b = A[0u][1u];
            const X& c = A[1u][0u];
            const X& d = A[1u][1u];
            // sqr() avoids spurious negative interval products.
            const X discriminant = sqr(a - d) + X(4u, _precision) * b * c;
            const X trace = a + d;
            if (!this->_finite(discriminant) || !this->_finite(trace)) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            if (this->_negative(discriminant)) {
                const X real_part = trace / two_value;
                const X imaginary_part = sqrt(-discriminant) / two_value;
                result.eigenvalues.emplace_back(real_part, -imaginary_part);
                result.eigenvalues.emplace_back(real_part, imaginary_part);
            } else if (this->_positive(discriminant)) {
                const X root = sqrt(discriminant);
                const X determinant = a * d - b * c;
                X lower(zero), upper(zero);
                if (this->_positive(trace)) {
                    upper = (trace + root) / two_value;
                    if (!this->_nonzero(upper)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    lower = determinant / upper;
                } else if (this->_negative(trace)) {
                    lower = (trace - root) / two_value;
                    if (!this->_nonzero(lower)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    upper = determinant / lower;
                } else {
                    lower = (trace - root) / two_value;
                    upper = (trace + root) / two_value;
                }
                result.eigenvalues.emplace_back(lower, zero);
                result.eigenvalues.emplace_back(upper, zero);
            } else if (this->_zero(discriminant)) {
                const X repeated = trace / two_value;
                result.eigenvalues.emplace_back(repeated, zero);
                result.eigenvalues.emplace_back(repeated, zero);
            } else {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
        } else {
            bool upper = true, lower = true;
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != n; ++j) {
                    if (i > j && !this->_zero(A[i][j])) { upper = false; }
                    if (i < j && !this->_zero(A[i][j])) { lower = false; }
                }
            }
            if (upper || lower) {
                for (SizeType i = 0u; i != n; ++i) {
                    result.eigenvalues.emplace_back(A[i][i], zero);
                }
            } else {
                // A_ij is a structural edge i -> j unless it is
                // provably exactly zero. A row/column permutation
                // topologically orders the strongly connected components,
                // producing an upper block-triangular matrix. The
                // characteristic polynomial is the product of the
                // characteristic polynomials of its diagonal blocks.
                std::vector<std::vector<bool>> reachable(
                    n, std::vector<bool>(n, false));
                for (SizeType i = 0u; i != n; ++i) {
                    for (SizeType j = 0u; j != n; ++j) {
                        reachable[i][j] = i == j || !this->_zero(A[i][j]);
                    }
                }
                for (SizeType through = 0u; through != n; ++through) {
                    for (SizeType i = 0u; i != n; ++i) {
                        for (SizeType j = 0u; j != n; ++j) {
                            if (reachable[i][through] && reachable[through][j]) {
                                reachable[i][j] = true;
                            }
                        }
                    }
                }

                std::vector<SizeType> component_of(n, n);
                std::vector<std::vector<SizeType>> components;
                for (SizeType i = 0u; i != n; ++i) {
                    if (component_of[i] != n) {
                        continue;
                    }
                    const SizeType id = components.size();
                    std::vector<SizeType> component;
                    for (SizeType j = 0u; j != n; ++j) {
                        if (component_of[j] == n &&
                            reachable[i][j] && reachable[j][i]) {
                            component_of[j] = id;
                            component.push_back(j);
                        }
                    }
                    if (component.size() > 2u) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    components.push_back(std::move(component));
                }

                const SizeType count = components.size();
                std::vector<std::vector<bool>> edges(
                    count, std::vector<bool>(count, false));
                std::vector<SizeType> indegree(count, 0u);
                for (SizeType i = 0u; i != n; ++i) {
                    for (SizeType j = 0u; j != n; ++j) {
                        const SizeType from = component_of[i];
                        const SizeType to = component_of[j];
                        if (from != to && !this->_zero(A[i][j]) &&
                            !edges[from][to]) {
                            edges[from][to] = true;
                            ++indegree[to];
                        }
                    }
                }
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
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    placed[next] = true;
                    for (SizeType id = 0u; id != count; ++id) {
                        if (edges[next][id]) {
                            --indegree[id];
                        }
                    }
                    const auto& indices = components[next];
                    if (indices.size() == 1u) {
                        result.eigenvalues.emplace_back(
                            A[indices[0u]][indices[0u]], zero);
                    } else {
                        Matrix<X> block(2u, 2u, zero);
                        for (SizeType i = 0u; i != 2u; ++i) {
                            for (SizeType j = 0u; j != 2u; ++j) {
                                block[i][j] = A[indices[i]][indices[j]];
                            }
                        }
                        // Reuse exactly the small-order approximate or
                        // outward-rounded validated spectral contract.
                        const auto solved = (*this)(block);
                        if (solved.status != Status::SUCCESS || !solved.data) {
                            return {Status::INCONCLUSIVE, std::nullopt};
                        }
                        for (const auto& lambda: solved.data->eigenvalues) {
                            result.eigenvalues.push_back(lambda);
                        }
                    }
                }
            }
        }
        for (const auto& lambda: result.eigenvalues) {
            if (!this->_finite(lambda.real_part()) ||
                !this->_finite(lambda.imaginary_part())) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
        }
        return {Status::SUCCESS, std::move(result)};
    }

  private:
    bool _finite(const X& x) const {
        if constexpr (std::is_same_v<P, ValidatedTag>) {
            return is_finite(x.lower_raw()) && is_finite(x.upper_raw());
        } else {
            return is_finite(x.raw());
        }
    }
    bool _positive(const X& x) const {
        const X z(0u, _precision);
        if constexpr (std::is_same_v<P, ValidatedTag>) {
            return x.lower_raw() > z.lower_raw();
        } else {
            return x.raw() > z.raw();
        }
    }
    bool _negative(const X& x) const {
        const X z(0u, _precision);
        if constexpr (std::is_same_v<P, ValidatedTag>) {
            return x.upper_raw() < z.upper_raw();
        } else {
            return x.raw() < z.raw();
        }
    }
    bool _zero(const X& x) const {
        const X z(0u, _precision);
        if constexpr (std::is_same_v<P, ValidatedTag>) {
            return x.lower_raw() == z.lower_raw() &&
                   x.upper_raw() == z.upper_raw();
        } else {
            return x.raw() == z.raw();
        }
    }
    bool _nonzero(const X& x) const {
        return this->_positive(x) || this->_negative(x);
    }
    PR _precision;
};

} // namespace Ariadne

#endif // ARIADNE_SPECTRUM_HPP
