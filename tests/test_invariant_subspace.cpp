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

//! \brief Test exact fixed invariant graphs for entire interval families.
//! \details Use the graph Riccati identity, NOT midpoint residuals.
//! Only singleton-zero interval enclosures can prove SUCCESS.
template<class PR>
Void test_invariant_graph_certificates(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);
    const X two_value(2u, precision);
    const X three_value(3u, precision);

    // Noncoordinate invariant line span((1,1)) for the point
    // triangular matrix [2 1; 0 3], eigenvalue 3.
    Matrix<X> A(2u, 2u, zero);
    A[0u][0u] = two_value;
    A[0u][1u] = one;
    A[1u][1u] = three_value;
    Matrix<X> graph(1u, 1u, one);
    const auto line = certified_invariant_graph(A, graph, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        line.has_value());
    if (line) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            line->basis.row_size() == 2u &&
            line->basis.column_size() == 1u &&
            line->restriction.row_size() == 1u &&
            line->restriction.column_size() == 1u);
        const Matrix<X> action = A * line->basis;
        const Matrix<X> represented =
            line->basis * line->restriction;
        for (SizeType i = 0u; i != 2u; ++i) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                action[i][0u].lower_raw() <=
                    represented[i][0u].upper_raw() &&
                represented[i][0u].lower_raw() <=
                    action[i][0u].upper_raw());
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            line->restriction[0u][0u].lower_raw() ==
                three_value.lower_raw() &&
            line->restriction[0u][0u].upper_raw() ==
                three_value.upper_raw());
    }

    // Exact zero residual is essential. Perturbing the exterior
    // coupling, or admitting an uncertain graph coefficient,
    // must not be passed off as an interval-family certificate.
    Matrix<X> perturbed(A);
    perturbed[1u][0u] = one;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph(perturbed, graph, precision));
    Matrix<X> uncertain_graph(1u, 1u, X(0, 1, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph(A, uncertain_graph, precision));
    const Matrix<X> wrong_graph(2u, 1u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph(A, wrong_graph, precision));

    // A noncoordinate TWO-dimensional invariant graph even for a
    // defective Jordan matrix. The restriction remains defective:
    // a certificate of an invariant plane is NOT a full eigenbasis.
    // A=[2 1 0;0 2 0;0 0 2], X=[0 1] and U=[I2;X].
    Matrix<X> jordan(3u, 3u, zero);
    for (SizeType i = 0u; i != 3u; ++i) {
        jordan[i][i] = two_value;
    }
    jordan[0u][1u] = one;
    Matrix<X> plane_graph(1u, 2u, zero);
    plane_graph[0u][1u] = one;
    const auto plane = certified_invariant_graph(
        jordan, plane_graph, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        plane.has_value());
    if (plane) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            plane->basis.row_size() == 3u &&
            plane->basis.column_size() == 2u &&
            plane->restriction.row_size() == 2u &&
            plane->restriction.column_size() == 2u &&
            plane->basis[2u][1u].lower_raw() == one.lower_raw() &&
            plane->restriction[0u][1u].lower_raw() ==
                one.lower_raw());
        const Matrix<X> action = jordan * plane->basis;
        const Matrix<X> represented =
            plane->basis * plane->restriction;
        for (SizeType i = 0u; i != 3u; ++i) {
            for (SizeType j = 0u; j != 2u; ++j) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    action[i][j].lower_raw() <=
                        represented[i][j].upper_raw() &&
                    represented[i][j].lower_raw() <=
                        action[i][j].upper_raw());
            }
        }
    }

    // This uncertainty is INTERNAL to the invariant plane's
    // restriction and has no effect on the graph Riccati equation:
    // A[0][2] affects only the first row of B and X=[0,1]
    // selects its second row in X*B.
    Matrix<X> family(jordan);
    family[0u][2u] = X(1, 2, precision);
    const auto uniform = certified_invariant_graph(
        family, plane_graph, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        uniform.has_value());
    if (uniform) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            uniform->restriction[0u][1u].lower_raw() <=
                X(2u, precision).lower_raw() &&
            uniform->restriction[0u][1u].upper_raw() >=
                X(3u, precision).upper_raw());
    }
    // But an uncertain cross-graph block with a genuine effect
    // cannot be certified merely because its interval includes 0.
    family[2u][0u] = X(-1, 1, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph(family, plane_graph, precision));

    Matrix<X> empty(0u, 0u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph(empty, graph, precision));
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
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_invariant_graph_certificates(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_invariant_graph_certificates(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
