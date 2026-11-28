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

//! \brief Check a Krawczyk-certified graph family, not merely midpoint.
template<class PR>
Void test_certified_invariant_graph_near(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);
    const X two_value(2u, precision);
    const X three_value(3u, precision);
    const X four_value(4u, precision);
    const X seventh(7u, precision);
    const X radius = X(-1, 1, precision) / X(16u, precision);
    const X small = X(-1, 1, precision) / X(1024u, precision);

    // A noncoordinate invariant line for a point matrix, with
    // the nonlinear Riccati equation F(x)=x-x^2, isolated at x=1.
    Matrix<X> A(2u, 2u, zero);
    A[0u][0u] = two_value;
    A[0u][1u] = one;
    A[1u][1u] = three_value;
    Matrix<X> center(1u, 1u, one);
    Matrix<X> box(1u, 1u, one + radius);
    Matrix<X> inverse_jacobian(1u, 1u, X(-1, precision));
    const auto exact = certified_invariant_graph_near(
        A, center, box, inverse_jacobian, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        exact.has_value());
    if (exact) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            exact->basis_enclosure.row_size() == 2u &&
            exact->basis_enclosure.column_size() == 1u &&
            exact->graph_enclosure.row_size() == 1u &&
            exact->graph_enclosure.column_size() == 1u &&
            exact->graph_enclosure[0u][0u].lower_raw() <
                one.lower_raw() &&
            one.upper_raw() <
                exact->graph_enclosure[0u][0u].upper_raw() &&
            exact->restriction_enclosure.row_size() == 1u &&
            exact->restriction_enclosure.column_size() == 1u);
    }

    // A whole input family changes A21 independently. The graph
    // is NOT constant across the family, yet Krawczyk proves an
    // individually existing and locally unique x(A) in the box.
    A[1u][0u] = small;
    const auto family = certified_invariant_graph_near(
        A, center, box, inverse_jacobian, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        family.has_value());
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph(A, center, precision));

    // A false numerical seed must never be accepted without
    // strictly interior Krawczyk inclusion; no result on a box
    // that misses the original nearby root.
    Matrix<X> wrong_center(1u, 1u, X(4u, precision));
    Matrix<X> wrong_box(1u, 1u, X(4u, precision) + radius);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_near(
            A, wrong_center, wrong_box, inverse_jacobian, precision));
    Matrix<X> uncertain_center(1u, 1u, one + radius);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_near(
            A, uncertain_center, box, inverse_jacobian, precision));
    Matrix<X> flat_box(1u, 1u, one);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_near(
            A, center, flat_box, inverse_jacobian, precision));

    // Two coupled unknown graph coefficients. The graph is
    // X=[1,1], A11=diag(2,3), A12=[1,0]^T,
    // A21=[-4,-3], A22=7. The Jacobian is nonsingular
    // but not diagonal: its offdiagonal coupling is essential.
    Matrix<X> coupled(3u, 3u, zero);
    coupled[0u][0u] = two_value;
    coupled[1u][1u] = three_value;
    coupled[0u][2u] = one;
    coupled[2u][0u] = X(-4, precision);
    coupled[2u][1u] = X(-3, precision);
    coupled[2u][2u] = seventh;
    Matrix<X> pair_center(1u, 2u, one);
    Matrix<X> pair_box(1u, 2u, one + radius);
    Matrix<X> pair_preconditioner(2u, 2u, zero);
    pair_preconditioner[0u][0u] = one / four_value;
    pair_preconditioner[1u][1u] = one / four_value;
    const auto two_variables = certified_invariant_graph_near(
        coupled, pair_center, pair_box, pair_preconditioner,
        precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        two_variables.has_value());
    if (two_variables) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            two_variables->basis_enclosure.row_size() == 3u &&
            two_variables->basis_enclosure.column_size() == 2u &&
            two_variables->restriction_enclosure.row_size() == 2u &&
            two_variables->restriction_enclosure.column_size() == 2u);
    }
    coupled[2u][0u] += small;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        certified_invariant_graph_near(
            coupled, pair_center, pair_box, pair_preconditioner,
            precision).has_value());

    // Four unknown graph entries (m=k=2), with a genuinely
    // coupled 4-by-4 Jacobian. This tests row-major indexing
    // rather than only the scalar or one-row graph cases.
    // At X=ones(2,2), B=[[3,1],[0,4]] and F_A(X)=0.
    Matrix<X> four_unknowns(4u, 4u, zero);
    four_unknowns[0u][0u] = two_value;
    four_unknowns[1u][1u] = four_value;
    four_unknowns[0u][2u] = one;
    four_unknowns[2u][0u] = X(-6, precision);
    four_unknowns[2u][1u] = X(-4, precision);
    four_unknowns[3u][0u] = X(-8, precision);
    four_unknowns[3u][1u] = X(-6, precision);
    four_unknowns[2u][2u] = X(9u, precision);
    four_unknowns[3u][3u] = X(11u, precision);
    Matrix<X> square_center(2u, 2u, one);
    Matrix<X> square_box(2u, 2u, one + radius);
    Matrix<X> square_preconditioner(4u, 4u, zero);
    for (SizeType i = 0u; i != 4u; ++i) {
        square_preconditioner[i][i] = one / X(8u, precision);
    }
    const auto four_variables = certified_invariant_graph_near(
        four_unknowns, square_center, square_box,
        square_preconditioner, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        four_variables.has_value());
    if (four_variables) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            four_variables->graph_enclosure.row_size() == 2u &&
            four_variables->graph_enclosure.column_size() == 2u &&
            four_variables->basis_enclosure.row_size() == 4u &&
            four_variables->basis_enclosure.column_size() == 2u &&
            four_variables->restriction_enclosure.row_size() == 2u &&
            four_variables->restriction_enclosure.column_size() == 2u);
    }
    four_unknowns[2u][0u] += small;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        certified_invariant_graph_near(
            four_unknowns, square_center, square_box,
            square_preconditioner, precision).has_value());

    // C=0 makes the Krawczyk defect I, which has norm exactly
    // one. A center with zero residual MUST NOT be accepted
    // without the independent strict contraction predicate.
    const Matrix<X> zero_preconditioner(4u, 4u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_near(
            four_unknowns, square_center, square_box,
            zero_preconditioner, precision));

    // Wrong dimensions: the preconditioner is in graph-variable
    // space, not in the n-by-n original matrix space.
    Matrix<X> wrong_preconditioner(3u, 3u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_near(
            coupled, pair_center, pair_box, wrong_preconditioner,
            precision));
}

//! \brief From numerical Schur graphs to rigorous family certification.
//! \details These tests keep approximate candidates separate from the
//! Krawczyk interval proof (Stewart-Sun, Bai-Demmel, Krawczyk).
template<class PR>
Void test_schur_graph_candidates(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    const AX azero(0u, precision);
    const AX aone(1u, precision);
    const VX vzero(0u, precision);
    const VX vone(1u, precision);

    // Select the NONLEADING 1x1 Schur block with eigenvalue 3.
    // Its eigenvector for [[2,1],[0,3]] is (1,1), yielding
    // a noncoordinate graph X=[1] over the first coordinate.
    Matrix<AX> A(2u, 2u, azero);
    A[0u][0u] = AX(2u, precision);
    A[0u][1u] = aone;
    A[1u][1u] = AX(3u, precision);
    const auto candidate = approximate_invariant_graph_from_schur(
        A, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        candidate.has_value());
    if (candidate) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            candidate->row_size() == 1u &&
            candidate->column_size() == 1u &&
            std::abs(((*candidate)[0u][0u] - aone).get_d()) < 1e-8);
    }

    // The same center is generated from the midpoint of a
    // genuinely uncertain matrix family. A21 varies freely
    // around 0, so its invariant graph MUST vary with A.
    // The certificate must be obtained on the ORIGINAL
    // interval family rather than the midpoint matrix.
    Matrix<VX> family(2u, 2u, vzero);
    family[0u][0u] = VX(2u, precision);
    family[0u][1u] = vone;
    family[1u][1u] = VX(3u, precision);
    family[1u][0u] =
        VX(-1, 1, precision) / VX(1024u, precision);
    Matrix<VX> box(1u, 1u, vone +
        VX(-1, 1, precision) / VX(16u, precision));
    const Matrix<VX> preconditioner(1u, 1u, VX(-1, precision));
    const auto certified = certified_invariant_graph_from_schur(
        family, {1u}, box, preconditioner, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        certified.has_value());
    if (certified) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            certified->graph_enclosure.row_size() == 1u &&
            certified->graph_enclosure.column_size() == 1u &&
            certified->graph_enclosure[0u][0u].lower_raw() <
                vone.lower_raw() &&
            vone.upper_raw() <
                certified->graph_enclosure[0u][0u].upper_raw());
    }

    // A candidate alone cannot certify an uncertain family.
    // A box that excludes the true graph for some input A
    // must fail the strict inclusion test.
    const Matrix<VX> narrow_box(
        1u, 1u, vone + VX(-1, 1, precision) / VX(8192u, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(
            family, {1u}, narrow_box, preconditioner, precision));

    // Invalid permutations, missing spectral blocks and a
    // singular graph chart are numerical INCONCLUSIVE outcomes.
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !approximate_invariant_graph_from_schur(
            A, {0u, 1u}, precision) &&
        !approximate_invariant_graph_from_schur(
            A, {2u}, precision));

    // The repeated-eigenvalue Jordan case cannot swap its
    // coincident 1x1 Schur blocks via a singular Sylvester
    // equation, so it cannot produce this candidate.
    Matrix<AX> jordan(A);
    jordan[1u][1u] = AX(2u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !approximate_invariant_graph_from_schur(
            jordan, {1u}, precision));

    // A pure complex pair occupies one indivisible real
    // 2x2 Schur block. There is no proper real invariant
    // graph subspace to return from such a selection.
    Matrix<AX> rotation(2u, 2u, azero);
    rotation[0u][1u] = AX(-1, precision);
    rotation[1u][0u] = aone;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !approximate_invariant_graph_from_schur(
            rotation, {0u}, precision));
}

//! \brief Fully automatic Schur candidate and Krawczyk box for a cluster.
//! \details The numerical Schur candidate MUST pass the interval
//! inclusion/contractivity proof on each original input family.
template<class PR>
Void test_automatic_schur_graph_certificate(PR precision) {
    using VX = FloatBounds<PR>;
    const VX zero(0u, precision);
    const VX one(1u, precision);
    const VX small = VX(-1, 1, precision) / VX(1024u, precision);

    // A=[2 1;0 3] has a non-coordinate eigenline span((1,1)).
    // The second original Schur block is moved to the front.
    Matrix<VX> A(2u, 2u, zero);
    A[0u][0u] = VX(2u, precision);
    A[0u][1u] = one;
    A[1u][1u] = VX(3u, precision);
    const auto exact = certified_invariant_graph_from_schur(
        A, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        exact.has_value());
    if (exact) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            exact->graph_enclosure.row_size() == 1u &&
            exact->graph_enclosure.column_size() == 1u &&
            exact->graph_enclosure[0u][0u].lower_raw() <
                one.lower_raw() &&
            exact->graph_enclosure[0u][0u].upper_raw() >
                one.upper_raw());
    }

    // An actual interval family with varying A21 has different
    // exact invariant graph X(A) for different point matrices.
    // An approximate midpoint eigenvector is NOT a certificate:
    // the automatic routine must prove existence for the family.
    A[1u][0u] = small;
    const auto family = certified_invariant_graph_from_schur(
        A, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        family.has_value());
    if (family) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            family->graph_enclosure[0u][0u].lower_raw() <
                one.lower_raw() &&
            family->graph_enclosure[0u][0u].upper_raw() >
                one.upper_raw() &&
            family->basis_enclosure.row_size() == 2u);
    }
    const auto lower_root = certified_invariant_graph_from_schur(
        A, {0u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        lower_root.has_value());
    if (lower_root) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            lower_root->graph_enclosure[0u][0u].lower_raw() < zero.lower_raw() &&
            lower_root->graph_enclosure[0u][0u].upper_raw() > zero.upper_raw());
    }

    // Repeat with two graph coordinates (m=2, k=1), testing the
    // row-major Riccati Jacobian and numerical inverse.
    Matrix<VX> three_by_three(3u, 3u, zero);
    three_by_three[0u][0u] = VX(2u, precision);
    three_by_three[0u][1u] = one;
    three_by_three[1u][1u] = VX(3u, precision);
    three_by_three[2u][2u] = VX(5u, precision);
    const auto higher = certified_invariant_graph_from_schur(
        three_by_three, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        higher.has_value());
    if (higher) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            higher->graph_enclosure.row_size() == 2u &&
            higher->graph_enclosure.column_size() == 1u &&
            higher->basis_enclosure.row_size() == 3u &&
            higher->restriction_enclosure.row_size() == 1u);
    }
    three_by_three[1u][0u] = small;
    const auto higher_family = certified_invariant_graph_from_schur(
        three_by_three, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        higher_family.has_value());

    // The interval family now contains point matrices for
    // which a real invariant line with eigenvalue near 3
    // does not exist. The proof must REFUSE the family.
    A[1u][0u] = VX(-1, 1, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(A, {1u}, precision));

    // Zero-order, absent blocks, and a full Schur selection
    // cannot yield a proper real graph chart.
    const Matrix<VX> empty(0u, 0u, zero);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(empty, {0u}, precision) &&
        !certified_invariant_graph_from_schur(A, {}, precision) &&
        !certified_invariant_graph_from_schur(A, {0u, 1u}, precision));

    // A defective Jordan block is not a simple isolated graph
    // chart for its repeated second scalar Schur block.
    Matrix<VX> jordan(2u, 2u, zero);
    jordan[0u][0u] = VX(2u, precision);
    jordan[0u][1u] = one;
    jordan[1u][1u] = VX(2u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(
            jordan, {1u}, precision));

    // The conjugate pair of a real rotation is a single
    // indivisible 2x2 Schur block, not a real invariant line.
    Matrix<VX> rotation(2u, 2u, zero);
    rotation[0u][1u] = VX(-1, precision);
    rotation[1u][0u] = one;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(
            rotation, {0u}, precision));
}

//! \brief Certify the REAL invariant plane of a conjugate Schur pair.
//! \details Although the eigenvalues are nonreal, the conjugate pair
//! spans a real invariant graph. Schur block selection uses the
//! Bai-Demmel/Wilkinson real Schur theorem; the verification is
//! the original-input interval Riccati/Krawczyk theorem (1969).
template<class PR>
Void test_automatic_schur_conjugate_plane(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);
    // 2 is a scalar Schur root; the trailing 2x2 rotation has
    // +/-i. The nonleading invariant plane has the exact
    // graph X=[1,0], so U=[e0+e2,e1] and A*U=U*[0 1;-1 0].
    Matrix<X> matrix(3u, 3u, zero);
    matrix[0u][0u] = X(2u, precision);
    matrix[0u][1u] = one;
    matrix[0u][2u] = X(-2, precision);
    matrix[1u][2u] = X(-1, precision);
    matrix[2u][1u] = one;

    const auto real_plane = certified_invariant_graph_from_schur(
        matrix, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        real_plane.has_value());
    if (real_plane) {
        const auto& graph = real_plane->graph_enclosure;
        const auto& basis = real_plane->basis_enclosure;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            graph.row_size() == 1u && graph.column_size() == 2u &&
            basis.row_size() == 3u && basis.column_size() == 2u &&
            graph[0u][0u].lower_raw() < one.lower_raw() &&
            graph[0u][0u].upper_raw() > one.upper_raw() &&
            graph[0u][1u].lower_raw() < zero.lower_raw() &&
            graph[0u][1u].upper_raw() > zero.upper_raw());
    }

    // Perturb an UPPER coupling independently in the original
    // interval family. The real invariant plane changes with
    // A; a Schur midpoint basis alone cannot certify it.
    matrix[0u][1u] +=
        X(-1, 1, precision) / X(1048576u, precision);
    const auto family = certified_invariant_graph_from_schur(
        matrix, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        family.has_value());
    if (family) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            family->graph_enclosure.row_size() == 1u &&
            family->graph_enclosure.column_size() == 2u &&
            family->graph_enclosure[0u][0u].lower_raw() < one.lower_raw() &&
            family->graph_enclosure[0u][0u].upper_raw() > one.upper_raw());
    }

    // Selecting the scalar block describes another subspace,
    // not the requested conjugate plane. Reject a whole-space
    // selection as an inappropriate graph chart.
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(
            matrix, {0u, 1u}, precision));
}

//! \brief Distinguish isolated near-Jordan roots from an uncertain
//! family containing point matrices with NONREAL eigenvalues.
//! \details Krawczyk must establish a uniformly existing real
//! graph, not merely converge on the midpoint (Stewart-Sun 1990).
template<class PR>
Void test_automatic_schur_nearly_defective(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);
    const X two_value(2u, precision);

    // The two eigenvalues of [[2,1],[0,2+delta]] are separated
    // by delta=1/64. The nonleading real eigenvector has graph
    // x=delta; a sufficiently thin input family preserves a
    // locally unique real invariant graph.
    const X separation = one / X(64u, precision);
    Matrix<X> matrix(2u, 2u, zero);
    matrix[0u][0u] = two_value;
    matrix[0u][1u] = one;
    matrix[1u][1u] = two_value + separation;
    const auto isolated = certified_invariant_graph_from_schur(
        matrix, {1u}, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        isolated.has_value());
    if (isolated) {
        const X& graph = isolated->graph_enclosure[0u][0u];
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            graph.lower_raw() < separation.lower_raw() &&
            graph.upper_raw() > separation.upper_raw());
    }
    matrix[1u][0u] =
        X(-1, 1, precision) / X(4194304u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        certified_invariant_graph_from_schur(
            matrix, {1u}, precision).has_value());

    // With delta=1/1024 and an independent A21 interval
    // [-delta^2,delta^2], the characteristic discriminant
    // delta^2+4*A21 is negative for some point matrices.
    // Those members admit NO invariant real line, so an
    // enclosure common to the whole family cannot be proved.
    const X narrow_separation = one / X(1024u, precision);
    matrix[1u][1u] = two_value + narrow_separation;
    matrix[1u][0u] =
        X(-1, 1, precision) / X(1048576u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(
            matrix, {1u}, precision));

    // Exactly coalesced Jordan roots remain INCONCLUSIVE for
    // the nonleading Schur-block swap, not SUCCESS by tolerance.
    matrix[1u][1u] = two_value;
    matrix[1u][0u] = zero;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !certified_invariant_graph_from_schur(
            matrix, {1u}, precision));
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
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_certified_invariant_graph_near(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_certified_invariant_graph_near(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_schur_graph_candidates(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_schur_graph_candidates(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_automatic_schur_graph_certificate(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_automatic_schur_graph_certificate(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_automatic_schur_conjugate_plane(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_automatic_schur_conjugate_plane(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_automatic_schur_nearly_defective(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_automatic_schur_nearly_defective(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
