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

//! \brief Two coupled complex conjugate pairs for implicit Francis QR.
template<class PR>
Void test_francis_two_pairs(PR precision) {
    using AX = FloatApproximation<PR>;
    const AX zero_value(0u, precision);
    // U*diag([[0,-1],[1,0]], [[0,-3],[3,0]])*U with
    // U=I-(1/2)*ones(4,4), an exactly orthogonal matrix.
    // The strongly connected integer input has characteristic
    // polynomial (lambda^2+1)(lambda^2+9).
    const Int coefficients[4][4] = {
        {0, 0, -1, 2},
        {0, 0, -2, 1},
        {1, 2, 0, 0},
        {-2, -1, 0, 0}
    };
    Matrix<AX> A(4u, 4u, zero_value);
    for (SizeType row = 0u; row != 4u; ++row) {
        for (SizeType col = 0u; col != 4u; ++col) {
            A[row][col] = AX(coefficients[row][col], precision);
        }
    }

    const auto schur = SchurDecomposition<ApproximateTag, PR>(precision)(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        schur.status == SchurStatus::SUCCESS && schur.data.has_value());
    if (schur.data) {
        const auto& data = *schur.data;
        SizeType covered = 0u, conjugate_pairs = 0u;
        for (const auto width: data.block_sizes) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                width > 0u && width <= 4u - covered);
            covered += width;
            if (width == 2u) { ++conjugate_pairs; }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            covered == 4u && conjugate_pairs == 2u);
        const Matrix<AX> aq = A * data.orthogonal;
        const Matrix<AX> qt = data.orthogonal * data.quasi_triangular;
        const Matrix<AX> orth = transpose(data.orthogonal) * data.orthogonal;
        for (SizeType row = 0u; row != 4u; ++row) {
            for (SizeType col = 0u; col != 4u; ++col) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs((aq[row][col] - qt[row][col]).get_d()) < 1e-8 &&
                    std::abs(orth[row][col].get_d() -
                        (row == col ? 1.0 : 0.0)) < 1e-8);
            }
        }
    }

    const auto roots = Spectrum<ApproximateTag, PR>(precision)(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        roots.status == SpectralStatus::SUCCESS && roots.data.has_value());
    if (roots.data) {
        const double expected[4] = {-3.0, -1.0, 1.0, 3.0};
        SizeType hits[4] = {0u, 0u, 0u, 0u};
        for (const auto& lambda: roots.data->eigenvalues) {
            const double re = lambda.real_part().get_d();
            const double im = lambda.imaginary_part().get_d();
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(re) < 1e-8);
            for (SizeType i = 0u; i != 4u; ++i) {
                if (std::abs(im - expected[i]) < 1e-8) {
                    ++hits[i];
                }
            }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            roots.data->eigenvalues.size() == 4u);
        for (SizeType i = 0u; i != 4u; ++i) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                hits[i] == 1u);
        }
    }
}

//! \brief Verify the invariant-subspace identity without diagonalizability.
//! \details A U = U B and U^T U = I follow from a complete leading
//! block partition in the real Schur theorem (Wilkinson, 1965).
template<class PR>
Void check_schur_invariant_prefix(
        const Matrix<FloatApproximation<PR>>& A,
        const SchurDecompositionData<ApproximateTag, PR>& factors,
        SizeType block_count, PR precision) {
    const auto subspace = factors.leading_invariant_subspace(
        block_count, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        subspace.has_value());
    if (!subspace) { return; }
    const SizeType n = A.row_size();
    const SizeType k = subspace->basis.column_size();
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        k > 0u && k <= n &&
        subspace->basis.row_size() == n &&
        subspace->restriction.row_size() == k &&
        subspace->restriction.column_size() == k);
    const Matrix<FloatApproximation<PR>> action = A * subspace->basis;
    const Matrix<FloatApproximation<PR>> represented =
        subspace->basis * subspace->restriction;
    const Matrix<FloatApproximation<PR>> gram =
        transpose(subspace->basis) * subspace->basis;
    const Matrix<FloatApproximation<PR>> compressed =
        transpose(subspace->basis) * action;
    for (SizeType row = 0u; row != n; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs((action[row][col] - represented[row][col]).get_d()) <
                    1e-8);
        }
    }
    for (SizeType row = 0u; row != k; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(gram[row][col].get_d() -
                    (row == col ? 1.0 : 0.0)) < 1e-8 &&
                std::abs((compressed[row][col] -
                    subspace->restriction[row][col]).get_d()) < 1e-8);
        }
    }
}

//! \brief Defective, nonreal and malformed Schur prefix cases.
template<class PR>
Void test_leading_invariant_subspaces(PR precision) {
    using AX = FloatApproximation<PR>;
    const AX zero_value(0u, precision);
    SchurDecomposition<ApproximateTag, PR> solver(precision);

    // A defective 3x3 Jordan block has just ONE independent
    // eigenvector, yet three nested invariant Schur prefixes.
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
    if (defective.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            defective.data->block_sizes.size() == 3u);
        for (SizeType blocks = 1u; blocks <= 3u; ++blocks) {
            check_schur_invariant_prefix(jordan, *defective.data,
                                         blocks, precision);
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            !defective.data->leading_invariant_subspace(0u, precision) &&
            !defective.data->leading_invariant_subspace(4u, precision));

        // The partition is checked independently, not trusted blindly.
        auto malformed = *defective.data;
        malformed.block_sizes = {1u};
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            !malformed.leading_invariant_subspace(1u, precision));
        malformed.block_sizes = {1u, 0u, 2u};
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            !malformed.leading_invariant_subspace(1u, precision));
    }

    // A 2x2 rotation has NO invariant real line. It is one
    // indivisible real Schur block with conjugate eigenvalues.
    Matrix<AX> rotation(2u, 2u, zero_value);
    rotation[0u][1u] = AX(-1, precision);
    rotation[1u][0u] = AX(1u, precision);
    const auto nonreal = solver(rotation);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        nonreal.status == SchurStatus::SUCCESS &&
        nonreal.data.has_value());
    if (nonreal.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            nonreal.data->block_sizes.size() == 1u &&
            nonreal.data->block_sizes[0u] == 2u);
        check_schur_invariant_prefix(rotation, *nonreal.data, 1u, precision);
        auto invalid = *nonreal.data;
        invalid.block_sizes = {1u, 1u};
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            !invalid.leading_invariant_subspace(1u, precision));
    }

    // Fully coupled 3x3, orthogonally reduced by Schur. A first
    // block may have width 1 or 2; either is invariant, but an
    // incomplete pair is NOT. Also verify the whole space.
    const Int coefficients[3][3] = {
        {4, 5, 0}, {3, 4, -4}, {-4, 0, 1}
    };
    Matrix<AX> mixed(3u, 3u, zero_value);
    for (SizeType row = 0u; row != 3u; ++row) {
        for (SizeType col = 0u; col != 3u; ++col) {
            mixed[row][col] = AX(coefficients[row][col], precision);
        }
    }
    const auto mixed_schur = solver(mixed);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        mixed_schur.status == SchurStatus::SUCCESS &&
        mixed_schur.data.has_value());
    if (mixed_schur.data) {
        const auto& factors = *mixed_schur.data;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            factors.block_sizes.size() == 2u);
        check_schur_invariant_prefix(mixed, factors, 1u, precision);
        check_schur_invariant_prefix(mixed, factors, 2u, precision);

        SizeType pair_start = 0u;
        for (const auto width: factors.block_sizes) {
            if (width == 2u) { break; }
            pair_start += width;
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            pair_start < 3u);
        auto split = factors;
        split.block_sizes = {1u, 1u, 1u};
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            !split.leading_invariant_subspace(pair_start + 1u, precision));
    }
}

//! \brief Reorder adjacent Schur blocks and select nonleading subspaces.
//! \details Small Sylvester exchanges follow orthogonal real Schur
//! reordering (Bai-Demmel, LAPACK Working Note 54 / DTREXC).
template<class PR>
Void test_reordered_schur_subspaces(PR precision) {
    using AX = FloatApproximation<PR>;
    const AX zero_value(0u, precision);

    // Two distinct 1x1 blocks. Swapping their invariant graph
    // must expose the nonleading eigenvalue 5 as a leading line.
    Matrix<AX> real(2u, 2u, zero_value);
    real[0u][0u] = AX(2u, precision);
    real[0u][1u] = AX(1u, precision);
    real[1u][1u] = AX(5u, precision);
    const SchurDecompositionData<ApproximateTag, PR> real_factors{
        Matrix<AX>::identity(2u, precision), real, {1u, 1u}};
    const std::vector<SizeType> swap_two{1u, 0u};
    const auto swapped = real_factors.reorder_blocks(swap_two, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        swapped.has_value());
    if (swapped) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            swapped->block_sizes.size() == 2u &&
            std::abs((swapped->quasi_triangular[0u][0u] -
                AX(5u, precision)).get_d()) < 1e-8);
        check_schur_invariant_prefix(real, *swapped, 1u, precision);
        check_schur_invariant_prefix(real, *swapped, 2u, precision);
    }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        std::abs((real_factors.quasi_triangular[0u][0u] -
            AX(2u, precision)).get_d()) < 1e-8);

    const std::vector<SizeType> duplicate{0u, 0u};
    const std::vector<SizeType> wrong_size{0u};
    const std::vector<SizeType> out_of_range{0u, 2u};
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !real_factors.reorder_blocks(duplicate, precision) &&
        !real_factors.reorder_blocks(wrong_size, precision) &&
        !real_factors.reorder_blocks(out_of_range, precision));

    // Equal diagonal eigenvalues with a nonzero coupling give a
    // SINGULAR Sylvester equation. Never certify a block swap.
    Matrix<AX> repeated(real);
    repeated[1u][1u] = AX(2u, precision);
    const SchurDecompositionData<ApproximateTag, PR> repeated_factors{
        Matrix<AX>::identity(2u, precision), repeated, {1u, 1u}};
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !repeated_factors.reorder_blocks(swap_two, precision));

    // Interleaved block types 1x1, 2x2 and 1x1. Reorder to
    // last scalar / conjugate pair / first scalar, traversing
    // all combinations of adjacent 1x1 and 2x2 exchanges.
    Matrix<AX> mixed(4u, 4u, zero_value);
    mixed[0u][0u] = AX(7u, precision);
    mixed[1u][2u] = AX(-3, precision);
    mixed[2u][1u] = AX(3u, precision);
    mixed[3u][3u] = AX(9u, precision);
    mixed[0u][1u] = AX(2u, precision);
    mixed[0u][2u] = AX(-1, precision);
    mixed[0u][3u] = AX(3u, precision);
    mixed[1u][3u] = AX(1u, precision);
    mixed[2u][3u] = AX(-2, precision);
    const SchurDecompositionData<ApproximateTag, PR> mixed_factors{
        Matrix<AX>::identity(4u, precision), mixed, {1u, 2u, 1u}};
    const std::vector<SizeType> reversed{2u, 1u, 0u};
    const auto mixed_reordered =
        mixed_factors.reorder_blocks(reversed, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        mixed_reordered.has_value());
    if (mixed_reordered) {
        const auto& data = *mixed_reordered;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            data.block_sizes.size() == 3u &&
            data.block_sizes[0u] == 1u &&
            data.block_sizes[1u] == 2u &&
            data.block_sizes[2u] == 1u &&
            std::abs((data.quasi_triangular[0u][0u] -
                AX(9u, precision)).get_d()) < 1e-8 &&
            std::abs((data.quasi_triangular[3u][3u] -
                AX(7u, precision)).get_d()) < 1e-8);
        for (SizeType k = 1u; k <= 3u; ++k) {
            check_schur_invariant_prefix(mixed, data, k, precision);
        }
    }

    // Two different conjugate pairs with a nonzero coupling.
    // Neither 2x2 block may be split by the swap.
    Matrix<AX> pairs(4u, 4u, zero_value);
    pairs[0u][1u] = AX(-1, precision);
    pairs[1u][0u] = AX(1u, precision);
    pairs[2u][3u] = AX(-3, precision);
    pairs[3u][2u] = AX(3u, precision);
    pairs[0u][2u] = AX(1u, precision);
    pairs[0u][3u] = AX(-2, precision);
    pairs[1u][2u] = AX(3u, precision);
    pairs[1u][3u] = AX(1u, precision);
    const SchurDecompositionData<ApproximateTag, PR> pair_factors{
        Matrix<AX>::identity(4u, precision), pairs, {2u, 2u}};
    const auto pair_reordered =
        pair_factors.reorder_blocks(swap_two, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        pair_reordered.has_value());
    if (pair_reordered) {
        const auto& data = *pair_reordered;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            data.block_sizes.size() == 2u &&
            data.block_sizes[0u] == 2u && data.block_sizes[1u] == 2u);
        for (SizeType i = 0u; i != 2u; ++i) {
            const SizeType k = 2u * i;
            const AX disc =
                sqr(data.quasi_triangular[k][k] -
                    data.quasi_triangular[k + 1u][k + 1u]) +
                AX(4u, precision) * data.quasi_triangular[k][k + 1u] *
                    data.quasi_triangular[k + 1u][k];
            const double imaginary = std::sqrt(-disc.get_d()) / 2.0;
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                disc.get_d() < 0.0 &&
                std::abs(imaginary - (i == 0u ? 3.0 : 1.0)) < 1e-8);
        }
        check_schur_invariant_prefix(pairs, data, 1u, precision);
        check_schur_invariant_prefix(pairs, data, 2u, precision);
    }

    // Actual Francis-computed Schur factors, not only synthetic
    // already-triangular matrices: reorder a 3x3 with spectrum
    // {9,-3i,+3i} and verify the first invariant subspace.
    const Int entries[3][3] = {
        {4, 5, 0}, {3, 4, -4}, {-4, 0, 1}
    };
    Matrix<AX> coupled(3u, 3u, zero_value);
    for (SizeType row = 0u; row != 3u; ++row) {
        for (SizeType col = 0u; col != 3u; ++col) {
            coupled[row][col] = AX(entries[row][col], precision);
        }
    }
    const auto schur =
        SchurDecomposition<ApproximateTag, PR>(precision)(coupled);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        schur.status == SchurStatus::SUCCESS && schur.data.has_value());
    if (schur.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            schur.data->block_sizes.size() == 2u);
        const auto reordered =
            schur.data->reorder_blocks(swap_two, precision);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            reordered.has_value());
        if (reordered) {
            check_schur_invariant_prefix(
                coupled, *reordered, 1u, precision);
            check_schur_invariant_prefix(
                coupled, *reordered, 2u, precision);
        }
    }
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
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_francis_two_pairs(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_francis_two_pairs(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_leading_invariant_subspaces(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_leading_invariant_subspaces(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_reordered_schur_subspaces(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_reordered_schur_subspaces(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
