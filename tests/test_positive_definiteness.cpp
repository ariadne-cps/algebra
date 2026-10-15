/***************************************************************************
 *            test_positive_definiteness.cpp
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
#include "algebra/positive_definiteness.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {

using X=FloatDPBounds;
using S=SymmetricMatrix<X>;
using Status=PositiveDefinitenessStatus;

Void test_empty_and_scalar_matrices()
{
    S empty(0u,X(0u,dp));
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(empty)==Status::POSITIVE_DEFINITE);

    S positive(1u,X(0u,dp));
    positive[0u][0u]=X(1u,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(positive)==Status::POSITIVE_DEFINITE);

    S negative(1u,X(0u,dp));
    negative[0u][0u]=X(-2,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(negative)==Status::NOT_POSITIVE_DEFINITE);

    S singular(1u,X(0u,dp));
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(singular)==Status::NOT_POSITIVE_DEFINITE);

    S uncertain(1u,X(0u,dp));
    uncertain[0u][0u]=X(-1,1,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(uncertain)==Status::INDETERMINATE);

    // A nonnegative interval containing zero does not prove strict positivity.
    S boundary(1u,X(0u,dp));
    boundary[0u][0u]=X(0,2,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(boundary)==Status::INDETERMINATE);
}

Void test_validated_schur_complements()
{
    // A strictly diagonally dominant, symmetric positive definite matrix.
    S positive(2u,X(0u,dp));
    positive[0u][0u]=X(2u,dp);
    positive[0u][1u]=X(1u,dp);
    positive[1u][1u]=X(2u,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(positive)==Status::POSITIVE_DEFINITE);

    // A negative second Schur pivot proves failure even with positive diagonal.
    S indefinite(2u,X(0u,dp));
    indefinite[0u][0u]=X(1u,dp);
    indefinite[0u][1u]=X(2u,dp);
    indefinite[1u][1u]=X(1u,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(indefinite)==Status::NOT_POSITIVE_DEFINITE);

    // A rank-one matrix has a zero exact second Schur pivot.
    S semidefinite(2u,X(0u,dp));
    semidefinite[0u][0u]=X(1u,dp);
    semidefinite[0u][1u]=X(1u,dp);
    semidefinite[1u][1u]=X(1u,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(semidefinite)==Status::NOT_POSITIVE_DEFINITE);

    // Uncertainty in the last Schur pivot cannot be collapsed to one sign.
    S uncertain(2u,X(0u,dp));
    uncertain[0u][0u]=X(1u,dp);
    uncertain[0u][1u]=X(1u,dp);
    uncertain[1u][1u]=X(0,2,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(uncertain)==Status::INDETERMINATE);

    // A negative first pivot rejects positive definiteness before elimination.
    S negative(2u,X(0u,dp));
    negative[0u][0u]=X(-1,dp);
    negative[1u][1u]=X(2u,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(negative)==Status::NOT_POSITIVE_DEFINITE);

    // Uncertainty in an earlier pivot must not be bypassed using later rows.
    S uncertain_first(2u,X(0u,dp));
    uncertain_first[0u][0u]=X(-1,1,dp);
    uncertain_first[1u][1u]=X(3u,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(uncertain_first)==Status::INDETERMINATE);
}

Void test_interval_matrix_family()
{
    S positive_family(3u,X(0u,dp));
    positive_family[0u][0u]=X(4,5,dp);
    positive_family[1u][1u]=X(4,5,dp);
    positive_family[2u][2u]=X(4,5,dp);
    positive_family[0u][1u]=X(-1,1,dp);
    positive_family[0u][2u]=X(-1,1,dp);
    positive_family[1u][2u]=X(-1,1,dp);
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(positive_family)==Status::POSITIVE_DEFINITE);

    S uncertain_family(2u,X(0u,dp));
    uncertain_family[0u][0u]=X(1u,dp);
    uncertain_family[0u][1u]=X(1u,dp);
    uncertain_family[1u][1u]=X(1,2,dp);
    // The family includes a singular member; interval evaluation may
    // either certify non-positive-definiteness only when universal or
    // return indeterminate. It must not certify positive definiteness.
    ARIADNE_TEST_ASSERT(validated_positive_definiteness(uncertain_family)==Status::INDETERMINATE);
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_empty_and_scalar_matrices());
    ARIADNE_TEST_CALL(test_validated_schur_complements());
    ARIADNE_TEST_CALL(test_interval_matrix_family());
    return ARIADNE_TEST_FAILURES;
}
