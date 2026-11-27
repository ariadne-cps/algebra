/***************************************************************************
 *            tests/test_determinant.cpp
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

/** @file tests/test_determinant.cpp
 *  @brief Regression tests for approximate and validated matrix determinants.
 */

#include <cmath>

#include "algebra/determinant.hpp"
#include "algebra_testing_context.hpp"
#include "numeric/numeric.hpp"

using namespace Ariadne;

namespace {

//! \brief Check that a validated determinant encloses an independently known integer.
template<class PR>
bool contains_integer(const FloatBounds<PR>& result, Int value, PR precision) {
    const FloatBounds<PR> expected(value, precision);
    return result.lower_raw() <= expected.lower_raw() &&
           result.upper_raw() >= expected.upper_raw();
}

//! \brief Exercise general dimensions, singularity and odd permutations.
template<class PR>
Void test_determinants(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;

    Matrix<AX> empty(0u, 0u, AX(0u, precision));
    Matrix<VX> validated_empty(0u, 0u, VX(0u, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), determinant(empty).get_d() == 1.0);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     contains_integer(determinant(validated_empty), 1, precision));

    Matrix<AX> scalar(1u, 1u, AX(0u, precision));
    scalar[0u][0u] = AX(-3, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), determinant(scalar).get_d() == -3.0);

    Matrix<AX> nonsymmetric(3u, 3u, AX(0u, precision));
    Matrix<VX> interval(3u, 3u, VX(0u, precision));
    const Int entries[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}};
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            nonsymmetric[i][j] = AX(entries[i][j], precision);
            interval[i][j] = VX(entries[i][j], precision);
        }
    }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     std::abs(determinant(nonsymmetric).get_d() + 3.0) < 1e-10);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     contains_integer(determinant(interval), -3, precision));

    // A row transposition changes the determinant's sign.
    Matrix<AX> permutation(4u, 4u, AX(0u, precision));
    for (SizeType i = 0u; i != 4u; ++i) {
        permutation[i][i] = AX(1u, precision);
    }
    permutation[0u][0u] = AX(0u, precision);
    permutation[1u][1u] = AX(0u, precision);
    permutation[0u][1u] = AX(1u, precision);
    permutation[1u][0u] = AX(1u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     determinant(permutation).get_d() == -1.0);

    // Singular inputs must yield an enclosure containing zero.
    for (SizeType j = 0u; j != 3u; ++j) {
        nonsymmetric[2u][j] = nonsymmetric[1u][j];
        interval[2u][j] = interval[1u][j];
    }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), determinant(nonsymmetric).get_d() == 0.0);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     contains_integer(determinant(interval), 0, precision));

    // Interval input denotes a whole family, not a single midpoint matrix.
    Matrix<VX> family(2u, 2u, VX(0u, precision));
    family[0u][0u] = VX(1, 2, precision);
    family[1u][1u] = VX(3, 4, precision);
    const VX family_det = determinant(family);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), contains_integer(family_det, 3, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), contains_integer(family_det, 8, precision));

    // Nontrivial dimension and fully separable known value (6!).
    Matrix<AX> diagonal(6u, 6u, AX(0u, precision));
    Matrix<VX> validated_diagonal(6u, 6u, VX(0u, precision));
    for (SizeType i = 0u; i != 6u; ++i) {
        diagonal[i][i] = AX(static_cast<unsigned>(i + 1u), precision);
        validated_diagonal[i][i] = VX(static_cast<unsigned>(i + 1u), precision);
    }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), determinant(diagonal).get_d() == 720.0);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     contains_integer(determinant(validated_diagonal), 720, precision));
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("determinant");
    suite.add("approximate and validated", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_determinants(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_determinants(MultiplePrecision(80u)));
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
