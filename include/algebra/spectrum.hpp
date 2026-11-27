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
//! \details This first stage handles 1x1, 2x2 and structurally upper/lower
//! triangular matrices. General dense higher-order matrices await Schur.
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
        const X two(2u, _precision);
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
                const X real_part = trace / two;
                const X imaginary_part = sqrt(-discriminant) / two;
                result.eigenvalues.emplace_back(real_part, -imaginary_part);
                result.eigenvalues.emplace_back(real_part, imaginary_part);
            } else if (this->_positive(discriminant)) {
                const X root = sqrt(discriminant);
                const X determinant = a * d - b * c;
                X lower(zero), upper(zero);
                if (this->_positive(trace)) {
                    upper = (trace + root) / two;
                    if (!this->_nonzero(upper)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    lower = determinant / upper;
                } else if (this->_negative(trace)) {
                    lower = (trace - root) / two;
                    if (!this->_nonzero(lower)) {
                        return {Status::INCONCLUSIVE, std::nullopt};
                    }
                    upper = determinant / lower;
                } else {
                    lower = (trace - root) / two;
                    upper = (trace + root) / two;
                }
                result.eigenvalues.emplace_back(lower, zero);
                result.eigenvalues.emplace_back(upper, zero);
            } else if (this->_zero(discriminant)) {
                const X repeated = trace / two;
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
            if (!upper && !lower) {
                return {Status::INCONCLUSIVE, std::nullopt};
            }
            for (SizeType i = 0u; i != n; ++i) {
                result.eigenvalues.emplace_back(A[i][i], zero);
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
