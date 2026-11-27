/***************************************************************************
 *            tests/test_eigen_decomposition.cpp
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

/** @file tests/test_eigen_decomposition.cpp
 *  @brief Test real approximate and validated eigendecomposition contracts.
 */

#include <cmath>

#include "algebra/eigen_decomposition.hpp"
#include "algebra_testing_context.hpp"
#include "numeric/numeric.hpp"

using namespace Ariadne;

namespace {

using Status = EigenDecompositionStatus;

// Check inclusion against exact integers using stored directed endpoints.
// A small floating point residual alone is not a validated certificate.
template<class PR>
bool encloses(const FloatBounds<PR>& interval, Int integer, PR pr) {
    const FloatBounds<PR> exact(integer, pr);
    return interval.lower_raw() <= exact.lower_raw()
        && interval.upper_raw() >= exact.upper_raw();
}

template<class PR>
Void test_real_examples(PR pr) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;

    Matrix<AX> a(2u, 2u, AX(0u, pr));
    a[0u][0u] = AX(2u, pr);
    a[0u][1u] = AX(1u, pr);
    a[1u][0u] = AX(1u, pr);
    a[1u][1u] = AX(2u, pr);

    EigenDecomposition<ApproximateTag, PR> approximate(pr);
    const auto approximation = approximate(a);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), approximation.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), approximation.data.has_value());
    if (approximation.data) {
        const auto& values = approximation.data->eigenvalues;
        const auto& vectors = approximation.data->eigenvectors;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(values[0u].get_d() - 1.0) < 1e-12);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(values[1u].get_d() - 3.0) < 1e-12);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(vectors[1u][0u].get_d() + 1.0) < 1e-12);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(vectors[1u][1u].get_d() - 1.0) < 1e-12);
    }

    Matrix<VX> v(2u, 2u, VX(0u, pr));
    v[0u][0u] = VX(2u, pr);
    v[0u][1u] = VX(1u, pr);
    v[1u][0u] = VX(1u, pr);
    v[1u][1u] = VX(2u, pr);

    EigenDecomposition<ValidatedTag, PR> validated(pr);
    const auto certificate = validated(v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), certificate.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), certificate.data.has_value());
    if (certificate.data) {
        const auto& eigen = *certificate.data;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(eigen.eigenvalues[0u], 1, pr));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(eigen.eigenvalues[1u], 3, pr));
        const auto product = v * eigen.eigenvectors;
        for (SizeType row = 0u; row != 2u; ++row) {
            for (SizeType col = 0u; col != 2u; ++col) {
                const auto expected = eigen.eigenvalues[col] * eigen.eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    product[row][col].lower_raw() <= expected.upper_raw()
                    && product[row][col].upper_raw() >= expected.lower_raw());
            }
        }
    }

    // Nonsymmetric, real distinct spectrum; b!=0 branch.
    v[0u][0u] = VX(2u, pr);
    v[0u][1u] = VX(1u, pr);
    v[1u][0u] = VX(0u, pr);
    v[1u][1u] = VX(3u, pr);
    const auto triangular = validated(v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), triangular.status == Status::SUCCESS);
    if (triangular.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(triangular.data->eigenvalues[0u], 2, pr));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(triangular.data->eigenvalues[1u], 3, pr));
    }

    // b=0, c!=0 branch and diagonal branch (including a zero eigenvalue).
    v[0u][0u] = VX(1u, pr);
    v[0u][1u] = VX(0u, pr);
    v[1u][0u] = VX(1u, pr);
    v[1u][1u] = VX(2u, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), validated(v).status == Status::SUCCESS);
    v[1u][0u] = VX(0u, pr);
    v[0u][0u] = VX(0u, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), validated(v).status == Status::SUCCESS);

    // A 1x1 matrix is decomposable even with an interval entry.
    Matrix<VX> scalar(1u, 1u, VX(0u, pr));
    scalar[0u][0u] = VX(1, 2, pr);
    const auto one_by_one = validated(scalar);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), one_by_one.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), one_by_one.data.has_value());
    if (one_by_one.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            one_by_one.data->eigenvalues[0u].lower_raw() <= VX(1u, pr).lower_raw());
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            one_by_one.data->eigenvalues[0u].upper_raw() >= VX(2u, pr).upper_raw());
    }

    Matrix<VX> unsupported(3u, 3u, VX(0u, pr));
    const auto no_dimension = validated(unsupported);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), no_dimension.status == Status::UNSUPPORTED_DIMENSION);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !no_dimension.data.has_value());
}

Void test_validated_uncertainty() {
    using X = FloatDPBounds;
    EigenDecomposition<ValidatedTag, DoublePrecision> solver(dp);
    Matrix<X> M(2u, 2u, X(0u, dp));

    // Irrational eigenvalues +/- sqrt(2), checked through exact squares.
    M[0u][1u] = X(1u, dp);
    M[1u][0u] = X(2u, dp);
    auto irrational = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), irrational.status == Status::SUCCESS);
    if (irrational.data) {
        for (SizeType j = 0u; j != 2u; ++j) {
            const X square = sqr(irrational.data->eigenvalues[j]);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(square, 2, dp));
        }
    }

    // A nontrivial interval matrix family: every member has two real
    // separated eigenvalues with a certified eigenvector basis.
    M[0u][0u] = X(2, 3, dp);
    M[0u][1u] = X(1u, dp);
    M[1u][0u] = X(1u, dp);
    M[1u][1u] = X(2u, dp);
    const auto family = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), family.status == Status::SUCCESS);

    // A rotation has nonreal eigenvalues; a Jordan block has a double
    // eigenvalue. Neither may claim a real simple eigendecomposition.
    M[0u][0u] = X(0u, dp);
    M[0u][1u] = X(-1, dp);
    M[1u][0u] = X(1u, dp);
    M[1u][1u] = X(0u, dp);
    const auto rotation = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), rotation.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !rotation.data.has_value());

    M[0u][0u] = X(1u, dp);
    M[0u][1u] = X(1u, dp);
    M[1u][0u] = X(0u, dp);
    M[1u][1u] = X(1u, dp);
    const auto jordan = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), jordan.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !jordan.data.has_value());
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("eigen_decomposition");
    suite.add("real 1x1 and 2x2 cases", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_real_examples(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_real_examples(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_validated_uncertainty());
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
