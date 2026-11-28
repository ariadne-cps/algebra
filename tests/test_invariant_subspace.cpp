/***************************************************************************
 *            tests/test_invariant_subspace.cpp
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

/** @file tests/test_invariant_subspace.cpp
 *  @brief Exact interval-family coordinate invariant-subspace certificates.
 */
#include "algebra/invariant_subspace.hpp"
#include "algebra_testing_context.hpp"

using namespace Ariadne;

namespace {

//! \brief Check the full interval invariance identity, not a midpoint.
template<class PR>
Void check_certificate(const Matrix<FloatBounds<PR>>& A,
                       const CertifiedCoordinateInvariantSubspaceData<PR>&
                           certificate,
                       PR precision) {
    using X = FloatBounds<PR>;
    const SizeType n = A.row_size();
    const SizeType k = certificate.coordinates.size();
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        k > 0u && certificate.basis.row_size() == n &&
        certificate.basis.column_size() == k &&
        certificate.restriction.row_size() == k &&
        certificate.restriction.column_size() == k);
    // U is a point coordinate basis, so the interval product A*U
    // equals U*B entrywise for the entire matrix family.
    const Matrix<X> action = A * certificate.basis;
    const Matrix<X> represented =
        certificate.basis * certificate.restriction;
    const X exact_zero(0u, precision);
    const X exact_one(1u, precision);
    for (SizeType row = 0u; row != n; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            const X lhs = action[row][col];
            const X rhs = represented[row][col];
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                lhs.lower_raw() <= rhs.upper_raw() &&
                rhs.lower_raw() <= lhs.upper_raw());
            const X expected_basis =
                certificate.coordinates[col] == row ? exact_one : exact_zero;
            const X observed = certificate.basis[row][col];
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                observed.lower_raw() == expected_basis.lower_raw() &&
                observed.upper_raw() == expected_basis.upper_raw());
        }
    }
    for (SizeType row = 0u; row != k; ++row) {
        for (SizeType col = 0u; col != k; ++col) {
            const X original =
                A[certificate.coordinates[row]][certificate.coordinates[col]];
            const X restricted = certificate.restriction[row][col];
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                original.lower_raw() == restricted.lower_raw() &&
                original.upper_raw() == restricted.upper_raw());
        }
    }
}

template<class PR>
Void test_coordinate_certificates(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);

    // A reducible 4x4 interval family. Coordinates 0 and 2 form
    // a NONCONSECUTIVE invariant subspace because all entries
    // from exterior rows 1,3 into columns 0,2 are exact zeros.
    // Its restriction is an uncertain upper triangular matrix.
    Matrix<X> A(4u, 4u, zero);
    A[0u][0u] = X(1, 2, precision);
    A[0u][2u] = X(1, 3, precision);
    A[2u][2u] = X(4, 5, precision);
    A[1u][1u] = X(6u, precision);
    A[3u][3u] = X(7u, precision);
    A[0u][1u] = one;
    A[2u][1u] = one;
    A[1u][3u] = one;

    const std::vector<SizeType> selected{2u, 0u};
    const auto certificate = certified_coordinate_invariant_subspace(
        A, selected, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        certificate.has_value());
    if (certificate) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            certificate->coordinates.size() == 2u &&
            certificate->coordinates[0u] == 2u &&
            certificate->coordinates[1u] == 0u);
        check_certificate(A, *certificate, precision);
    }

    // Merely containing zero is insufficient for a uniform proof.
    auto uncertain = A;
    uncertain[1u][2u] = X(-1, 1, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_coordinate_invariant_subspace(
            uncertain, selected, precision));
    // A certain NONZERO cross-block entry also destroys invariance.
    uncertain[1u][2u] = one;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_coordinate_invariant_subspace(
            uncertain, selected, precision));

    // Other couplings may be uncertain: they only influence the
    // restriction and do not invalidate the structural certificate.
    A[2u][0u] = X(-1, 1, precision);
    const auto coupled = certified_coordinate_invariant_subspace(
        A, selected, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        coupled.has_value());
    if (coupled) { check_certificate(A, *coupled, precision); }

    const std::vector<SizeType> duplicate{0u, 0u};
    const std::vector<SizeType> nonexistent{0u, 4u};
    const std::vector<SizeType> none;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_coordinate_invariant_subspace(A, duplicate, precision) &&
        !certified_coordinate_invariant_subspace(A, nonexistent, precision) &&
        !certified_coordinate_invariant_subspace(A, none, precision));

    // All coordinates always span an invariant space for finite A.
    const std::vector<SizeType> everything{3u, 1u, 0u, 2u};
    const auto whole = certified_coordinate_invariant_subspace(
        A, everything, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        whole.has_value());
    if (whole) { check_certificate(A, *whole, precision); }

    // Defective Jordan 3x3: span(e0,e1) is invariant even though
    // no complete eigenvector basis exists.
    Matrix<X> jordan(3u, 3u, zero);
    for (SizeType i = 0u; i != 3u; ++i) {
        jordan[i][i] = X(2u, precision);
    }
    jordan[0u][1u] = one;
    jordan[1u][2u] = one;
    const auto defective = certified_coordinate_invariant_subspace(
        jordan, {0u, 1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        defective.has_value());
    if (defective) { check_certificate(jordan, *defective, precision); }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_coordinate_invariant_subspace(
            jordan, {1u, 2u}, precision));

    // Rotation has a nonreal conjugate pair but its full real plane
    // remains invariant with an exact real basis.
    Matrix<X> rotation(2u, 2u, zero);
    rotation[0u][1u] = X(-1, precision);
    rotation[1u][0u] = one;
    const auto full_rotation = certified_coordinate_invariant_subspace(
        rotation, {0u, 1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        full_rotation.has_value());
    if (full_rotation) {
        check_certificate(rotation, *full_rotation, precision);
    }
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_coordinate_invariant_subspace(
            rotation, {0u}, precision));

    // Nonfinite inputs cannot be certified even when their
    // problematic coefficient lies outside the requested blocks.
    // Test zero order, too: there is no nonempty coordinate basis.
    Matrix<X> empty(0u, 0u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_coordinate_invariant_subspace(
            empty, {0u}, precision));
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("invariant_subspace");
    suite.add("exact interval-family coordinate invariant subspaces",
        [](Ariadne::Testing::Context& ctx) {
            AlgebraTesting::CaseScope capture(ctx);
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_coordinate_certificates(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_coordinate_certificates(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
