/***************************************************************************
 *            tests/test_schur_decomposition.cpp
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

/** @file tests/test_schur_decomposition.cpp
 *  @brief Real quasi-triangular factors and independent complex spectra.
 */
#include <cmath>

#include "algebra/schur_decomposition.hpp"
#include "algebra/spectrum.hpp"
#include "algebra_testing_context.hpp"

using namespace Ariadne;

namespace {

template<class PR>
Void test_schur(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    const AX zero_value(0u, precision);
    SchurDecomposition<ApproximateTag, PR> solver(precision);
    SchurDecomposition<ValidatedTag, PR> no_certificate(precision);
    Spectrum<ApproximateTag, PR> spectrum(precision);
    Spectrum<ValidatedTag, PR> rigorous(precision);

    // Similarity of a rotation block and a scalar, with exactly
    // known characteristic polynomial (lambda-9)(lambda^2+9).
    // All three vertices belong to one strongly connected component.
    const Int entries[3][3] = {{4, 5, 0}, {3, 4, -4}, {-4, 0, 1}};
    Matrix<AX> A(3u, 3u, zero_value);
    Matrix<VX> V(3u, 3u, VX(0u, precision));
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            A[i][j] = AX(entries[i][j], precision);
            V[i][j] = VX(entries[i][j], precision);
        }
    }
    const auto result = solver(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        result.status == SchurStatus::SUCCESS && result.data.has_value());
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        no_certificate(V).status == SchurStatus::INCONCLUSIVE);
    if (result.data) {
        const auto& data = *result.data;
        SizeType dimensions = 0u, complex_blocks = 0u;
        for (const auto width: data.block_sizes) {
            // Every Schur block must fit in the original matrix.
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                width > 0u && width <= 3u - dimensions);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                width == 1u || width == 2u);
            dimensions += width;
            if (width == 2u) { ++complex_blocks; }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            dimensions == 3u && complex_blocks == 1u);
        const Matrix<AX> aq = A * data.orthogonal;
        const Matrix<AX> qt = data.orthogonal * data.quasi_triangular;
        const Matrix<AX> orth = transpose(data.orthogonal) * data.orthogonal;
        for (SizeType i = 0u; i != 3u; ++i) {
            for (SizeType j = 0u; j != 3u; ++j) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs((aq[i][j] - qt[i][j]).get_d()) < 1e-8 &&
                    std::abs(orth[i][j].get_d() -
                        (i == j ? 1.0 : 0.0)) < 1e-8);
            }
        }
    }
    const auto spec = spectrum(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        spec.status == SpectralStatus::SUCCESS && spec.data.has_value());
    if (spec.data) {
        SizeType real_count = 0u, positive_count = 0u, negative_count = 0u;
        for (const auto& lambda: spec.data->eigenvalues) {
            const double re = lambda.real_part().get_d();
            const double im = lambda.imaginary_part().get_d();
            if (std::abs(re - 9.0) < 1e-8 && std::abs(im) < 1e-8) {
                ++real_count;
            } else if (std::abs(re) < 1e-8 &&
                       std::abs(im - 3.0) < 1e-8) {
                ++positive_count;
            } else if (std::abs(re) < 1e-8 &&
                       std::abs(im + 3.0) < 1e-8) {
                ++negative_count;
            }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            real_count == 1u && positive_count == 1u && negative_count == 1u);
    }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        rigorous(V).status == SpectralStatus::INCONCLUSIVE);

    // Distinct real eigenvalues from a fully coupled symmetric matrix.
    const Int symmetric[3][3] = {{2, 1, 1}, {1, 3, 1}, {1, 1, 4}};
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            A[i][j] = AX(symmetric[i][j], precision);
        }
    }
    const auto symmetric_result = solver(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        symmetric_result.status == SchurStatus::SUCCESS &&
        symmetric_result.data.has_value());
    if (symmetric_result.data) {
        for (const auto width: symmetric_result.data->block_sizes) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                width == 1u);
        }
    }

    // A defective Jordan matrix still possesses a real Schur form.
    Matrix<AX> jordan(3u, 3u, zero_value);
    for (SizeType i = 0u; i != 3u; ++i) {
        jordan[i][i] = AX(2u, precision);
    }
    jordan[0u][1u] = AX(1u, precision);
    jordan[1u][2u] = AX(1u, precision);
    const auto defective = solver(jordan);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        defective.status == SchurStatus::SUCCESS &&
        defective.data.has_value());
    const auto roots = spectrum(jordan);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        roots.status == SpectralStatus::SUCCESS && roots.data.has_value());
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("schur_decomposition");
    suite.add("approximate real Schur and dense complex spectra",
        [](Ariadne::Testing::Context& ctx) {
            AlgebraTesting::CaseScope capture(ctx);
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_schur(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_schur(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
