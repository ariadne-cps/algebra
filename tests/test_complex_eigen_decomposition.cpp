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

#include "algebra/eigen_decomposition.hpp"
#include "algebra_testing_context.hpp"

using namespace Ariadne;

namespace {

template<class PR>
Void check_complex_factorization(const Matrix<FloatApproximation<PR>>& A,
                                  PR precision) {
    using X = FloatApproximation<PR>;
    using CX = Complex<X>;
    const SizeType n = A.row_size();
    const detail::ComplexEigenDecomposition<ApproximateTag, PR> solver(precision);
    // Compare the direct method with the shared-Schur overload: both
    // must construct a complete complex eigenbasis, and neither may
    // bypass the original-input residual and inverse-basis checks.
    const auto schur = SchurDecomposition<ApproximateTag, PR>(precision)(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        schur.status == SchurStatus::SUCCESS && schur.data.has_value());
    const auto result = solver(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        result.status == detail::ComplexEigenDecompositionStatus::SUCCESS &&
        result.data.has_value());
    if (schur.data) {
        const auto reused = solver(A, *schur.data);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            reused.status == detail::ComplexEigenDecompositionStatus::SUCCESS &&
            reused.data.has_value());
        if (reused.data && result.data) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                reused.data->eigenvalues.size() == result.data->eigenvalues.size());
            for (SizeType k = 0u; k != result.data->eigenvalues.size(); ++k) {
                const auto& first = reused.data->eigenvalues[k];
                const auto& second = result.data->eigenvalues[k];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs((first.real_part() - second.real_part()).get_d()) < 1e-8 &&
                    std::abs((first.imaginary_part() - second.imaginary_part()).get_d()) < 1e-8);
            }
        }
    }
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
    using Status = detail::ComplexEigenDecompositionStatus;
    const X zero(0u, precision);
    detail::ComplexEigenDecomposition<ApproximateTag, PR> solver(precision);

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
    const detail::ComplexEigenDecomposition<ValidatedTag, PR> validated(precision);
    const auto uncertified = validated(validated_rotation);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        uncertified.status == Status::INCONCLUSIVE &&
        !uncertified.data.has_value());

    Matrix<X> empty(0u, 0u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        solver(empty).status == Status::UNSUPPORTED_DIMENSION);
}

//! \brief Reject malformed Schur metadata without manufacturing an eigenbasis.
template<class PR>
Void test_invalid_complex_schur_factors(PR precision) {
    using X = FloatApproximation<PR>;
    using Status = detail::ComplexEigenDecompositionStatus;
    Matrix<X> rotation(2u, 2u, X(0u, precision));
    rotation[0u][1u] = X(-1, precision);
    rotation[1u][0u] = X(1u, precision);
    const auto original = SchurDecomposition<ApproximateTag, PR>(precision)(rotation);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        original.status == SchurStatus::SUCCESS && original.data.has_value());
    if (!original.data) { return; }
    const detail::ComplexEigenDecomposition<ApproximateTag, PR> solver(precision);
    auto reject = [&](const auto& factors) {
        const auto result = solver(rotation, factors);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            result.status == Status::INCONCLUSIVE && !result.data.has_value());
    };
    auto invalid_T = *original.data;
    invalid_T.quasi_triangular = Matrix<X>(1u, 1u, X(0u, precision));
    reject(invalid_T);
    auto invalid_Q = *original.data;
    invalid_Q.orthogonal = Matrix<X>(1u, 1u, X(0u, precision));
    reject(invalid_Q);
    auto empty = *original.data;
    empty.block_sizes.clear();
    reject(empty);
    auto width = *original.data;
    width.block_sizes = {3u};
    reject(width);
    auto incomplete = *original.data;
    incomplete.block_sizes = {1u};
    reject(incomplete);
    auto overfull = *original.data;
    overfull.block_sizes = {1u, 1u, 1u};
    reject(overfull);
}

//! \brief Verify automatic real/complex dispatch through one public API.
template<class PR>
Void test_unified_eigen_decomposition(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    using Status = EigenDecompositionStatus;
    const AX zero(0u, precision);
    const EigenDecomposition<ApproximateTag, PR> solver(precision);

    // A symmetric matrix with only real eigenvalues must use
    // the real representation and must never allocate a
    // ComplexMatrix on the successful fast path.
    Matrix<AX> real_symmetric(2u, 2u, zero);
    real_symmetric[0u][0u] = AX(2u, precision);
    real_symmetric[1u][1u] = AX(3u, precision);
    real_symmetric[0u][1u] = AX(1u, precision);
    real_symmetric[1u][0u] = AX(1u, precision);
    EigenDecompositionDiagnostics diagnostics;
    const auto real_result = solver(real_symmetric, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        real_result.status == Status::SUCCESS &&
        real_result.real_data() != nullptr &&
        real_result.complex_data() == nullptr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy != EigenDecompositionStrategy::COMPLEX_SCHUR);

    // The same class is also used for matrices with nonreal
    // eigenvalues. These require complex vectors and inverse.
    Matrix<AX> rotation(2u, 2u, zero);
    rotation[0u][1u] = AX(-1, precision);
    rotation[1u][0u] = AX(1u, precision);
    const auto complex_result = solver(rotation, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        complex_result.status == Status::SUCCESS &&
        complex_result.real_data() == nullptr &&
        complex_result.complex_data() != nullptr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy == EigenDecompositionStrategy::COMPLEX_SCHUR);
    if (complex_result.complex_data()) {
        const auto& data = *complex_result.complex_data();
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            data.eigenvalues.size() == 2u);
        for (SizeType i = 0u; i != 2u; ++i) {
            for (SizeType j = 0u; j != 2u; ++j) {
                Complex<AX> action(AX(0u, precision), AX(0u, precision));
                for (SizeType k = 0u; k != 2u; ++k) {
                    action += Complex<AX>(AX(rotation[i][k])) * data.eigenvectors[k][j];
                }
                const auto residual = action - data.eigenvalues[j] *
                    data.eigenvectors[i][j];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs(residual.real_part().get_d()) < 1e-8 &&
                    std::abs(residual.imaginary_part().get_d()) < 1e-8);
            }
        }
    }

    // A defective Jordan block does not have a complex basis:
    // neither result representation may claim success.
    Matrix<AX> jordan(2u, 2u, zero);
    jordan[0u][0u] = AX(2u, precision);
    jordan[0u][1u] = AX(1u, precision);
    jordan[1u][1u] = AX(2u, precision);
    const auto defective = solver(jordan);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        defective.status == Status::INCONCLUSIVE &&
        !defective.data.has_value());

    // Preserve validated real-basis certificates as real data.
    Matrix<VX> certified(2u, 2u, VX(0u, precision));
    certified[0u][0u] = VX(1u, precision);
    certified[1u][1u] = VX(3u, precision);
    const EigenDecomposition<ValidatedTag, PR> validated(precision);
    const auto valid_result = validated(certified);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        valid_result.status == Status::SUCCESS &&
        valid_result.real_data() != nullptr &&
        valid_result.complex_data() == nullptr);

    Matrix<VX> validated_rotation(2u, 2u, VX(0u, precision));
    validated_rotation[0u][1u] = VX(-1, precision);
    validated_rotation[1u][0u] = VX(1u, precision);
    const auto unverified = validated(validated_rotation);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        unverified.status == Status::INCONCLUSIVE &&
        !unverified.data.has_value());
}


//! \brief Reuse real nonsymmetric Schur factors on complex fallback.
template<class PR>
Void test_cached_complex_schur_fallback(PR precision) {
    using AX = FloatApproximation<PR>;
    using Status = EigenDecompositionStatus;
    // A strongly connected nonsymmetric 3x3 with spectrum
    // {9,-3i,+3i}; the real QR path reaches Schur's 2x2 block.
    const Int entries[3][3] = {{4, 5, 0}, {3, 4, -4}, {-4, 0, 1}};
    Matrix<AX> A(3u, 3u, AX(0u, precision));
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            A[i][j] = AX(entries[i][j], precision);
        }
    }
    std::optional<SchurDecompositionData<ApproximateTag, PR>> cached;
    const detail::RealEigenDecomposition<ApproximateTag, PR> real(precision);
    const auto real_result = real(A, nullptr, &cached);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        real_result.status == Status::INCONCLUSIVE &&
        !real_result.data.has_value() && cached.has_value());
    if (cached) {
        SizeType covered = 0u, pairs = 0u;
        for (const auto width: cached->block_sizes) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                width > 0u && width <= 3u - covered);
            covered += width;
            if (width == 2u) { ++pairs; }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            covered == 3u && pairs == 1u);
        const auto extracted =
            detail::ComplexEigenDecomposition<ApproximateTag, PR>(precision)(
                A, *cached);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            extracted.status ==
                detail::ComplexEigenDecompositionStatus::SUCCESS &&
            extracted.data.has_value());
    }

    EigenDecompositionDiagnostics diagnostics;
    const auto unified = EigenDecomposition<ApproximateTag, PR>(precision)(
        A, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        unified.status == Status::SUCCESS &&
        unified.real_data() == nullptr &&
        unified.complex_data() != nullptr &&
        diagnostics.strategy == EigenDecompositionStrategy::COMPLEX_SCHUR);

    // The cache is cleared by the next invocation and is NOT
    // populated by a real success (here a symmetric example).
    const Int real_entries[3][3] = {
        {6, 1, 0}, {1, 3, 1}, {0, 1, 0}
    };
    Matrix<AX> real_matrix(3u, 3u, AX(0u, precision));
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = 0u; j != 3u; ++j) {
            real_matrix[i][j] = AX(real_entries[i][j], precision);
        }
    }
    const auto real_success = real(real_matrix, nullptr, &cached);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        real_success.status == Status::SUCCESS &&
        real_success.data.has_value() && !cached.has_value());
    const auto unified_real =
        EigenDecomposition<ApproximateTag, PR>(precision)(real_matrix);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        unified_real.status == Status::SUCCESS &&
        unified_real.real_data() != nullptr &&
        unified_real.complex_data() == nullptr);
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
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_unified_eigen_decomposition(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_unified_eigen_decomposition(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_invalid_complex_schur_factors(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_invalid_complex_schur_factors(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_cached_complex_schur_fallback(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_cached_complex_schur_fallback(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
