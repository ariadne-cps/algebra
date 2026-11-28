/***************************************************************************
 *            tests/test_complex_eigen_decomposition.cpp
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

/** @file tests/test_complex_eigen_decomposition.cpp
 *  @brief Complex Schur eigenvectors, inverse bases and defective refusal.
 */
#include <cmath>

#include "algebra/complex_eigen_decomposition.hpp"
#include "algebra_testing_context.hpp"

using namespace Ariadne;

namespace {

template<class PR>
Void check_complex_factorization(const Matrix<FloatApproximation<PR>>& A,
                                  PR precision) {
    using X = FloatApproximation<PR>;
    using CX = Complex<X>;
    const SizeType n = A.row_size();
    const auto result = ComplexEigenDecomposition<ApproximateTag, PR>(precision)(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        result.status == ComplexEigenDecompositionStatus::SUCCESS &&
        result.data.has_value());
    if (!result.data) {
        return;
    }
    const auto& data = *result.data;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        data.eigenvalues.size() == n &&
        data.eigenvectors.row_size() == n &&
        data.eigenvectors.column_size() == n &&
        data.inverse_eigenvectors.row_size() == n &&
        data.inverse_eigenvectors.column_size() == n);
    for (SizeType i = 0u; i != n; ++i) {
        for (SizeType j = 0u; j != n; ++j) {
            CX action(X(0u, precision), X(0u, precision));
            CX inverse_left(X(0u, precision), X(0u, precision));
            CX inverse_right(X(0u, precision), X(0u, precision));
            for (SizeType k = 0u; k != n; ++k) {
                action += CX(A[i][k]) * data.eigenvectors[k][j];
                inverse_left += data.eigenvectors[i][k] *
                    data.inverse_eigenvectors[k][j];
                inverse_right += data.inverse_eigenvectors[i][k] *
                    data.eigenvectors[k][j];
            }
            const CX residual = action -
                data.eigenvalues[j] * data.eigenvectors[i][j];
            const double re = residual.real_part().get_d();
            const double im = residual.imaginary_part().get_d();
            const double left_real = inverse_left.real_part().get_d();
            const double left_imaginary = inverse_left.imaginary_part().get_d();
            const double right_real = inverse_right.real_part().get_d();
            const double right_imaginary = inverse_right.imaginary_part().get_d();
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(re) < 1e-8 && std::abs(im) < 1e-8 &&
                std::abs(left_real - (i == j ? 1.0 : 0.0)) < 1e-8 &&
                std::abs(left_imaginary) < 1e-8 &&
                std::abs(right_real - (i == j ? 1.0 : 0.0)) < 1e-8 &&
                std::abs(right_imaginary) < 1e-8);
        }
    }
}

template<class PR>
Void test_complex_diagonalization(PR precision) {
    using X = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    using Status = ComplexEigenDecompositionStatus;
    const X zero(0u, precision);
    ComplexEigenDecomposition<ApproximateTag, PR> solver(precision);

    Matrix<X> scalar(1u, 1u, zero);
    scalar[0u][0u] = X(7u, precision);
    check_complex_factorization(scalar, precision);

    // Real rotation: eigenvalues -i,+i and no complete REAL basis.
    Matrix<X> rotation(2u, 2u, zero);
    rotation[0u][1u] = X(-1, precision);
    rotation[1u][0u] = X(1u, precision);
    check_complex_factorization(rotation, precision);
    const auto rotation_result = solver(rotation);
    if (rotation_result.data) {
        const auto& eigenvalues = rotation_result.data->eigenvalues;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            eigenvalues.size() == 2u);
        if (eigenvalues.size() == 2u) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(eigenvalues[0u].real_part().get_d()) < 1e-8 &&
                std::abs(eigenvalues[1u].real_part().get_d()) < 1e-8 &&
                std::abs(eigenvalues[0u].imaginary_part().get_d() + 1.0) < 1e-8 &&
                std::abs(eigenvalues[1u].imaginary_part().get_d() - 1.0) < 1e-8);
        }
    }

    // Fully connected, nonsymmetric real matrix with one real root
    // and a complex conjugate pair. Exact characteristic polynomial
    // is (lambda-9)(lambda^2+9).
    const Int three[3][3] = {{4, 5, 0}, {3, 4, -4}, {-4, 0, 1}};
    Matrix<X> a3(3u, 3u, zero);
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            a3[i][j] = X(three[i][j], precision);
        }
    }
    check_complex_factorization(a3, precision);
    const auto a3_result = solver(a3);
    if (a3_result.data) {
        SizeType real_hits = 0u, negative_hits = 0u, positive_hits = 0u;
        for (const auto& lambda: a3_result.data->eigenvalues) {
            const double re = lambda.real_part().get_d();
            const double im = lambda.imaginary_part().get_d();
            if (std::abs(re - 9.0) < 1e-8 && std::abs(im) < 1e-8) {
                ++real_hits;
            } else if (std::abs(re) < 1e-8 && std::abs(im + 3.0) < 1e-8) {
                ++negative_hits;
            } else if (std::abs(re) < 1e-8 && std::abs(im - 3.0) < 1e-8) {
                ++positive_hits;
            }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            real_hits == 1u && negative_hits == 1u && positive_hits == 1u);
    }

    // Two complex pairs in one fully connected real 4x4 matrix.
    // Exact characteristic polynomial (lambda^2+1)(lambda^2+9).
    const Int four[4][4] = {
        {0, 0, -1, 2}, {0, 0, -2, 1},
        {1, 2, 0, 0}, {-2, -1, 0, 0}
    };
    Matrix<X> a4(4u, 4u, zero);
    for (SizeType i = 0u; i != 4u; ++i) {
        for (SizeType j = 0u; j != 4u; ++j) {
            a4[i][j] = X(four[i][j], precision);
        }
    }
    check_complex_factorization(a4, precision);

    // A dense nonsymmetric system with a simple REAL spectrum
    // also belongs to the complex-diagonalization domain.
    const Int real_entries[3][3] = {{6, 1, 0}, {1, 3, 1}, {0, 1, 0}};
    Matrix<X> a_real(3u, 3u, zero);
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            a_real[i][j] = X(real_entries[i][j], precision);
        }
    }
    check_complex_factorization(a_real, precision);

    // A defective Jordan block has a repeated eigenvalue but no
    // full complex eigenbasis: never report numerical success.
    Matrix<X> jordan(3u, 3u, zero);
    for (SizeType i = 0u; i != 3u; ++i) {
        jordan[i][i] = X(2u, precision);
    }
    jordan[0u][1u] = X(1u, precision);
    jordan[1u][2u] = X(1u, precision);
    const auto defective = solver(jordan);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        defective.status == Status::INCONCLUSIVE &&
        !defective.data.has_value());

    // ValidatedTag is not allowed to reuse a floating inverse as proof.
    Matrix<VX> validated_rotation(2u, 2u, VX(0u, precision));
    validated_rotation[0u][1u] = VX(-1, precision);
    validated_rotation[1u][0u] = VX(1u, precision);
    const ComplexEigenDecomposition<ValidatedTag, PR> validated(precision);
    const auto uncertified = validated(validated_rotation);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        uncertified.status == Status::INCONCLUSIVE &&
        !uncertified.data.has_value());

    Matrix<X> empty(0u, 0u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        solver(empty).status == Status::UNSUPPORTED_DIMENSION);
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("complex_eigen_decomposition");
    suite.add("complex Schur eigenvectors and inverse bases",
        [](Ariadne::Testing::Context& ctx) {
            AlgebraTesting::CaseScope capture(ctx);
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_complex_diagonalization(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_complex_diagonalization(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
