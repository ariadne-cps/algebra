/***************************************************************************
 *            test_linear_independence.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/numeric.hpp"
#include "algebra/matrix.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {

using X=FloatDPBounds;
using M=Matrix<X>;
using Status=LinearIndependenceStatus;

Void test_empty_matrices_and_dimension_limits()
{
    const X zero(0u,dp);

    M empty(0u,0u,zero);
    ARIADNE_TEST_ASSERT(validated_independent_rows(empty)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(empty)==Status::INDEPENDENT);

    M no_rows(0u,3u,zero);
    ARIADNE_TEST_ASSERT(validated_independent_rows(no_rows)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(no_rows)==Status::DEPENDENT);

    M no_columns(3u,0u,zero);
    ARIADNE_TEST_ASSERT(validated_independent_rows(no_columns)==Status::DEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(no_columns)==Status::INDEPENDENT);

    // For a nonempty input, the number of vectors cannot exceed dimension.
    M short_wide(1u,3u,zero);
    short_wide[0u][0u]=X(1u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(short_wide)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(short_wide)==Status::DEPENDENT);

    M tall_thin(3u,1u,zero);
    tall_thin[0u][0u]=X(1u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(tall_thin)==Status::DEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(tall_thin)==Status::INDEPENDENT);
}

Void test_full_rank_and_rank_deficient()
{
    const X zero(0u,dp);

    // Distinct rows and columns; the Gram matrices are strictly positive.
    M nonsingular(2u,2u,zero);
    nonsingular[0u][0u]=X(1u,dp);
    nonsingular[0u][1u]=X(2u,dp);
    nonsingular[1u][0u]=X(3u,dp);
    nonsingular[1u][1u]=X(4u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(nonsingular)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(nonsingular)==Status::INDEPENDENT);

    // Exact rank-one matrix: the second Schur complement vanishes.
    M rank_one(2u,2u,zero);
    rank_one[0u][0u]=X(1u,dp);
    rank_one[0u][1u]=X(2u,dp);
    rank_one[1u][0u]=X(2u,dp);
    rank_one[1u][1u]=X(4u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(rank_one)==Status::DEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(rank_one)==Status::DEPENDENT);

    // The zero Gram matrix proves dependence without numerical tolerance.
    M zero_matrix(1u,1u,zero);
    ARIADNE_TEST_ASSERT(validated_independent_rows(zero_matrix)==Status::DEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(zero_matrix)==Status::DEPENDENT);

    // A negative off-diagonal Gram entry must not break SPD certification.
    M orthogonal(2u,2u,zero);
    orthogonal[0u][0u]=X(1u,dp);
    orthogonal[0u][1u]=X(1u,dp);
    orthogonal[1u][0u]=X(-1,dp);
    orthogonal[1u][1u]=X(1u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(orthogonal)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(orthogonal)==Status::INDEPENDENT);
}

Void test_interval_families()
{
    const X zero(0u,dp);

    // All members are nonsingular: diagonal coefficients exclude zero.
    M full_rank_family(2u,2u,zero);
    full_rank_family[0u][0u]=X(2,3,dp);
    full_rank_family[1u][1u]=X(4,5,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(full_rank_family)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(full_rank_family)==Status::INDEPENDENT);

    // Some members have full rank and some are singular.
    // Independence cannot be inferred from the midpoint matrix alone.
    M unknown(1u,1u,zero);
    unknown[0u][0u]=X(-1,1,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(unknown)==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(validated_independent_columns(unknown)==Status::INDETERMINATE);

    // A nonnegative diagonal Gram interval touching zero is also unknown.
    M boundary(1u,1u,zero);
    boundary[0u][0u]=X(0,2,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(boundary)==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(validated_independent_columns(boundary)==Status::INDETERMINATE);

    // A positive interval matrix with a possibly rank-deficient member.
    M mixed(2u,2u,zero);
    mixed[0u][0u]=X(1u,dp);
    mixed[0u][1u]=X(1u,dp);
    mixed[1u][0u]=X(1u,dp);
    mixed[1u][1u]=X(1,2,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(mixed)==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(validated_independent_columns(mixed)==Status::INDETERMINATE);
}

Void test_rectangular_full_rank()
{
    const X zero(0u,dp);

    M independent_rows(2u,3u,zero);
    independent_rows[0u][0u]=X(2u,dp);
    independent_rows[0u][2u]=X(1u,dp);
    independent_rows[1u][1u]=X(3u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(independent_rows)==Status::INDEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(independent_rows)==Status::DEPENDENT);

    M independent_columns(3u,2u,zero);
    independent_columns[0u][0u]=X(2u,dp);
    independent_columns[2u][0u]=X(1u,dp);
    independent_columns[1u][1u]=X(3u,dp);
    ARIADNE_TEST_ASSERT(validated_independent_rows(independent_columns)==Status::DEPENDENT);
    ARIADNE_TEST_ASSERT(validated_independent_columns(independent_columns)==Status::INDEPENDENT);
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_empty_matrices_and_dimension_limits());
    ARIADNE_TEST_CALL(test_full_rank_and_rank_deficient());
    ARIADNE_TEST_CALL(test_interval_families());
    ARIADNE_TEST_CALL(test_rectangular_full_rank());
    return ARIADNE_TEST_FAILURES;
}
