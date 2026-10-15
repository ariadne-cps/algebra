/***************************************************************************
 *            test_cholesky_decomposition.cpp
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
#include "algebra/cholesky_decomposition.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {

using X=FloatDPBounds;
using S=SymmetricMatrix<X>;
using Status=PositiveDefinitenessStatus;

// Check an enclosure, not equality of interval matrices: for every
// admissible real A there is a real Cholesky factor within L.
Void check_reconstruction(const S& A, const Matrix<X>& L)
{
    ARIADNE_TEST_EQUALS(L.row_size(),A.size());
    ARIADNE_TEST_EQUALS(L.column_size(),A.size());

    Matrix<X> B=L*transpose(L);
    for(SizeType i=0; i!=A.size(); ++i) {
        for(SizeType j=0; j!=A.size(); ++j) {
            ARIADNE_TEST_ASSERT(definitely(B[i][j].lower()<=A[i][j].lower()));
            ARIADNE_TEST_ASSERT(definitely(B[i][j].upper()>=A[i][j].upper()));
            if(j>i) {
                ARIADNE_TEST_EQUALS(L[i][j],X(0u,dp));
            }
        }
    }
}

Void test_empty_and_scalar_matrices()
{
    S empty(0u,X(0u,dp));
    auto empty_result=validated_cholesky_decomposition(empty);
    ARIADNE_TEST_ASSERT(empty_result.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(empty_result.lower.has_value());
    check_reconstruction(empty,*empty_result.lower);

    S exact_square(1u,X(0u,dp));
    exact_square[0u][0u]=X(4u,dp);
    auto exact_result=validated_cholesky_decomposition(exact_square);
    ARIADNE_TEST_ASSERT(exact_result.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(exact_result.lower.has_value());
    ARIADNE_TEST_EQUALS((*exact_result.lower)[0u][0u],X(2u,dp));
    check_reconstruction(exact_square,*exact_result.lower);

    // The non-dyadic square root must be enclosed, not rounded to a point.
    S irrational(1u,X(0u,dp));
    irrational[0u][0u]=X(2u,dp);
    auto irrational_result=validated_cholesky_decomposition(irrational);
    ARIADNE_TEST_ASSERT(irrational_result.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(irrational_result.lower.has_value());
    check_reconstruction(irrational,*irrational_result.lower);

    S negative(1u,X(0u,dp));
    negative[0u][0u]=X(-1,dp);
    auto negative_result=validated_cholesky_decomposition(negative);
    ARIADNE_TEST_ASSERT(negative_result.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_result.lower.has_value());

    S zero(1u,X(0u,dp));
    auto zero_result=validated_cholesky_decomposition(zero);
    ARIADNE_TEST_ASSERT(zero_result.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!zero_result.lower.has_value());

    S uncertain(1u,X(0u,dp));
    uncertain[0u][0u]=X(-1,1,dp);
    auto uncertain_result=validated_cholesky_decomposition(uncertain);
    ARIADNE_TEST_ASSERT(uncertain_result.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_result.lower.has_value());

    S boundary(1u,X(0u,dp));
    boundary[0u][0u]=X(0,2,dp);
    auto boundary_result=validated_cholesky_decomposition(boundary);
    ARIADNE_TEST_ASSERT(boundary_result.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!boundary_result.lower.has_value());
}

Void test_nonsingular_and_degenerate_schur_complements()
{
    S positive(2u,X(0u,dp));
    positive[0u][0u]=X(2u,dp);
    positive[0u][1u]=X(1u,dp);
    positive[1u][1u]=X(2u,dp);
    auto positive_result=validated_cholesky_decomposition(positive);
    ARIADNE_TEST_ASSERT(positive_result.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(positive_result.lower.has_value());
    check_reconstruction(positive,*positive_result.lower);

    S indefinite(2u,X(0u,dp));
    indefinite[0u][0u]=X(1u,dp);
    indefinite[0u][1u]=X(2u,dp);
    indefinite[1u][1u]=X(1u,dp);
    auto indefinite_result=validated_cholesky_decomposition(indefinite);
    ARIADNE_TEST_ASSERT(indefinite_result.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!indefinite_result.lower.has_value());

    S semidefinite(2u,X(0u,dp));
    semidefinite[0u][0u]=X(1u,dp);
    semidefinite[0u][1u]=X(1u,dp);
    semidefinite[1u][1u]=X(1u,dp);
    auto semidefinite_result=validated_cholesky_decomposition(semidefinite);
    ARIADNE_TEST_ASSERT(semidefinite_result.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!semidefinite_result.lower.has_value());

    S uncertain(2u,X(0u,dp));
    uncertain[0u][0u]=X(1u,dp);
    uncertain[0u][1u]=X(1u,dp);
    uncertain[1u][1u]=X(0,2,dp);
    auto uncertain_result=validated_cholesky_decomposition(uncertain);
    ARIADNE_TEST_ASSERT(uncertain_result.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_result.lower.has_value());
}

Void test_interval_family()
{
    S family(3u,X(0u,dp));
    family[0u][0u]=X(4,5,dp);
    family[1u][1u]=X(4,5,dp);
    family[2u][2u]=X(4,5,dp);
    family[0u][1u]=X(-1,1,dp);
    family[0u][2u]=X(-1,1,dp);
    family[1u][2u]=X(-1,1,dp);
    auto family_result=validated_cholesky_decomposition(family);
    ARIADNE_TEST_ASSERT(family_result.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(family_result.lower.has_value());
    check_reconstruction(family,*family_result.lower);

    // Zero crossing after eliminating a certified positive first pivot.
    S uncertain(2u,X(0u,dp));
    uncertain[0u][0u]=X(1u,dp);
    uncertain[0u][1u]=X(1u,dp);
    uncertain[1u][1u]=X(1,2,dp);
    auto uncertain_result=validated_cholesky_decomposition(uncertain);
    ARIADNE_TEST_ASSERT(uncertain_result.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_result.lower.has_value());
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_empty_and_scalar_matrices());
    ARIADNE_TEST_CALL(test_nonsingular_and_degenerate_schur_complements());
    ARIADNE_TEST_CALL(test_interval_family());
    return ARIADNE_TEST_FAILURES;
}
