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

    // Avoid losing the small eigenvalue to subtraction cancellation.
    // Both examples have one eigenvalue very close to +1 or -1, while
    // the other has magnitude 1e16.
    Matrix<AX> ill_conditioned(2u, 2u, AX(0u, pr));
    ill_conditioned[0u][0u] = AX(1e16, pr);
    ill_conditioned[0u][1u] = AX(1u, pr);
    ill_conditioned[1u][0u] = AX(1u, pr);
    ill_conditioned[1u][1u] = AX(1u, pr);
    const auto stable_positive = approximate(ill_conditioned);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), stable_positive.status == Status::SUCCESS);
    if (stable_positive.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         std::abs(stable_positive.data->eigenvalues[0u].get_d() - 1.0) < 1e-6);
    }

    ill_conditioned[0u][0u] = AX(-1e16, pr);
    ill_conditioned[1u][1u] = AX(-1, pr);
    const auto stable_negative = approximate(ill_conditioned);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), stable_negative.status == Status::SUCCESS);
    if (stable_negative.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         std::abs(stable_negative.data->eigenvalues[1u].get_d() + 1.0) < 1e-6);
    }

    // Dense symmetric 3x3 Jacobi eigenpairs: exactly 2, 2 +/- sqrt(2).
    // Test the known spectrum independently and check M*V=V*D.
    Matrix<AX> symmetric(3u, 3u, AX(0u, pr));
    symmetric[0u][0u] = AX(2u, pr);
    symmetric[1u][1u] = AX(2u, pr);
    symmetric[2u][2u] = AX(2u, pr);
    symmetric[0u][1u] = AX(1u, pr);
    symmetric[1u][0u] = AX(1u, pr);
    symmetric[1u][2u] = AX(1u, pr);
    symmetric[2u][1u] = AX(1u, pr);
    const auto three = approximate(symmetric);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), three.status == Status::SUCCESS);
    if (three.data) {
        const double expected[3] = {2.0 - std::sqrt(2.0), 2.0, 2.0 + std::sqrt(2.0)};
        for (SizeType col = 0u; col != 3u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             std::abs(three.data->eigenvalues[col].get_d() - expected[col]) < 1e-10);
            for (SizeType row = 0u; row != 3u; ++row) {
                AX action(0u, pr);
                for (SizeType k = 0u; k != 3u; ++k) {
                    action += symmetric[row][k] * three.data->eigenvectors[k][col];
                }
                const AX residual = action - three.data->eigenvalues[col] * three.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                                 std::abs(residual.get_d()) < 1e-10);
            }
        }
    }

    // A dense symmetric 4x4 matrix with an independently constructed
    // orthogonal eigenbasis and known spectrum {1,2,3,4}.
    Matrix<AX> four_by_four(4u, 4u, AX(0u, pr));
    const double entries_4[4][4] = {
        {2.5, -0.5, -1.0, 0.0},
        {-0.5, 2.5, 0.0, -1.0},
        {-1.0, 0.0, 2.5, -0.5},
        {0.0, -1.0, -0.5, 2.5}
    };
    for (SizeType row = 0u; row != 4u; ++row) {
        for (SizeType col = 0u; col != 4u; ++col) {
            four_by_four[row][col] = AX(entries_4[row][col], pr);
        }
    }
    const auto four = approximate(four_by_four);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), four.status == Status::SUCCESS);
    if (four.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             std::abs(four.data->eigenvalues[col].get_d() - static_cast<double>(col + 1u)) < 1e-10);
        }
    }

    // Approximate nonsymmetric dense matrices remain unsupported. Validated
    // dense symmetric matrices must never be passed off as certified.
    symmetric[0u][1u] = AX(3u, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     approximate(symmetric).status == Status::UNSUPPORTED_DIMENSION);
    Matrix<VX> symmetric_bounds(3u, 3u, VX(0u, pr));
    symmetric_bounds[0u][0u] = VX(2u, pr);
    symmetric_bounds[1u][1u] = VX(2u, pr);
    symmetric_bounds[2u][2u] = VX(2u, pr);
    symmetric_bounds[0u][1u] = VX(1u, pr);
    symmetric_bounds[1u][0u] = VX(1u, pr);
    symmetric_bounds[1u][2u] = VX(1u, pr);
    symmetric_bounds[2u][1u] = VX(1u, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     validated(symmetric_bounds).status == Status::UNSUPPORTED_DIMENSION);

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

    // Real upper triangular matrix of order four. Its eigenvalues are
    // the diagonal entries, independent of the off-diagonal coefficients.
    Matrix<VX> upper(4u, 4u, VX(0u, pr));
    upper[0u][0u] = VX(4u, pr);
    upper[1u][1u] = VX(1u, pr);
    upper[2u][2u] = VX(3u, pr);
    upper[3u][3u] = VX(2u, pr);
    upper[0u][1u] = VX(2u, pr);
    upper[0u][2u] = VX(1u, pr);
    upper[1u][2u] = VX(-1, pr);
    upper[2u][3u] = VX(1u, pr);
    const auto upper_result = validated(upper);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), upper_result.status == Status::SUCCESS);
    if (upper_result.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             encloses(upper_result.data->eigenvalues[col], static_cast<Int>(col + 1u), pr));
        }
        const auto lhs = upper * upper_result.data->eigenvectors;
        for (SizeType row = 0u; row != 4u; ++row) {
            for (SizeType col = 0u; col != 4u; ++col) {
                const VX rhs = upper_result.data->eigenvalues[col] * upper_result.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                                 lhs[row][col].lower_raw() <= rhs.upper_raw() &&
                                                 lhs[row][col].upper_raw() >= rhs.lower_raw());
            }
        }
    }

    // Lower triangular case with the same pairwise distinct spectrum.
    Matrix<VX> lower(4u, 4u, VX(0u, pr));
    for (SizeType row = 0u; row != 4u; ++row) {
        for (SizeType col = 0u; col != 4u; ++col) {
            lower[row][col] = upper[col][row];
        }
    }
    const auto lower_result = validated(lower);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), lower_result.status == Status::SUCCESS);
    if (lower_result.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             encloses(lower_result.data->eigenvalues[col], static_cast<Int>(col + 1u), pr));
        }
    }

    // Wide diagonal intervals prevent a consistent eigenvalue ordering.
    Matrix<VX> overlapping(3u, 3u, VX(0u, pr));
    overlapping[0u][0u] = VX(0, 2, pr);
    overlapping[1u][1u] = VX(1, 3, pr);
    overlapping[2u][2u] = VX(4, 5, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     validated(overlapping).status == Status::INCONCLUSIVE);

    // An interval family with strictly separated diagonal values.
    overlapping[0u][0u] = VX(0, 1, pr);
    overlapping[1u][1u] = VX(2, 3, pr);
    const auto separated = validated(overlapping);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), separated.status == Status::SUCCESS);

    // A dense 3x3 matrix still requires a general spectral backend.
    // Dimension zero has no individual eigenpairs and is deliberately
    // excluded from this diagonalization API.
    Matrix<VX> empty(0u, 0u, VX(0u, pr));
    const auto empty_result = validated(empty);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), empty_result.status == Status::UNSUPPORTED_DIMENSION);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !empty_result.data.has_value());

    Matrix<VX> unsupported(3u, 3u, VX(0u, pr));
    unsupported[0u][1u] = VX(1u, pr);
    unsupported[1u][0u] = VX(1u, pr);

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
