/***************************************************************************
 *            test_positive_definite_solve.cpp
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
#include "algebra/positive_definite_solve.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {

using X=FloatDPBounds;
using S=SymmetricMatrix<X>;
using Status=PositiveDefinitenessStatus;

// Compare the stored exact floating-point endpoints: comparisons on
// LowerBound/UpperBound objects are not ordinary exact endpoint comparisons.
Void check_contains(const X& enclosure, const X& required)
{
    ARIADNE_TEST_ASSERT(enclosure.lower_raw()<=required.lower_raw());
    ARIADNE_TEST_ASSERT(enclosure.upper_raw()>=required.upper_raw());
}

// Check the algebraic residual as an interval inclusion, not equality of
// interval matrices. Also check known solutions independently of the solve.
Void check_matrix_solution(const S& A, const Matrix<X>& B,
                           const Matrix<X>& enclosure)
{
    ARIADNE_TEST_EQUALS(enclosure.row_size(),B.row_size());
    ARIADNE_TEST_EQUALS(enclosure.column_size(),B.column_size());
    const Matrix<X> dense_A=static_cast<Matrix<X>>(A);
    const Matrix<X> residual=dense_A*enclosure;
    for(SizeType i=0; i!=B.row_size(); ++i) {
        for(SizeType j=0; j!=B.column_size(); ++j) {
            check_contains(residual[i][j],B[i][j]);
        }
    }
}

Void check_vector_solution(const S& A, const Vector<X>& b,
                           const Vector<X>& enclosure)
{
    ARIADNE_TEST_EQUALS(enclosure.size(),b.size());
    const Matrix<X> dense_A=static_cast<Matrix<X>>(A);
    const Vector<X> residual=dense_A*enclosure;
    for(SizeType i=0; i!=b.size(); ++i) {
        check_contains(residual[i],b[i]);
    }
}

Void test_empty_and_zero_rhs()
{
    const X zero(0u,dp);
    const S empty(0u,zero);
    const Matrix<X> no_rhs(0u,0u,zero);
    const auto empty_matrix_solution=validated_positive_definite_solve(empty,no_rhs);
    ARIADNE_TEST_ASSERT(empty_matrix_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(empty_matrix_solution.solution.has_value());
    check_matrix_solution(empty,no_rhs,*empty_matrix_solution.solution);

    // Empty equation system with two separate empty right-hand sides.
    const Matrix<X> two_empty_rhs(0u,2u,zero);
    const auto two_empty_solution=validated_positive_definite_solve(empty,two_empty_rhs);
    ARIADNE_TEST_ASSERT(two_empty_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(two_empty_solution.solution.has_value());
    check_matrix_solution(empty,two_empty_rhs,*two_empty_solution.solution);

    const Vector<X> empty_rhs(0u,zero);
    const auto empty_vector_solution=validated_positive_definite_solve(empty,empty_rhs);
    ARIADNE_TEST_ASSERT(empty_vector_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(empty_vector_solution.solution.has_value());
    check_vector_solution(empty,empty_rhs,*empty_vector_solution.solution);

    const S positive(2u,X(0u,dp));
    // A non-SPD matrix remains non-SPD even when there are no right-hand sides.
    const Matrix<X> absent_rhs(2u,0u,zero);
    const auto absent_rhs_solution=validated_positive_definite_solve(positive,absent_rhs);
    ARIADNE_TEST_ASSERT(absent_rhs_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!absent_rhs_solution.solution.has_value());

    S identity(2u,zero);
    identity[0u][0u]=X(1u,dp);
    identity[1u][1u]=X(1u,dp);
    const auto zero_column_solution=validated_positive_definite_solve(identity,absent_rhs);
    ARIADNE_TEST_ASSERT(zero_column_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(zero_column_solution.solution.has_value());
    check_matrix_solution(identity,absent_rhs,*zero_column_solution.solution);
}

Void test_scalar_and_dense_spd()
{
    const X zero(0u,dp);

    S scalar(1u,zero);
    scalar[0u][0u]=X(4u,dp);
    const Vector<X> b1({X(8u,dp)});
    const auto scalar_solution=validated_positive_definite_solve(scalar,b1);
    ARIADNE_TEST_ASSERT(scalar_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(scalar_solution.solution.has_value());
    check_contains((*scalar_solution.solution)[0u],X(2u,dp));
    check_vector_solution(scalar,b1,*scalar_solution.solution);

    // The Cholesky factor is irrational; solution must still enclose 1.
    S irrational(1u,zero);
    irrational[0u][0u]=X(2u,dp);
    const Vector<X> b2({X(2u,dp)});
    const auto irrational_solution=validated_positive_definite_solve(irrational,b2);
    ARIADNE_TEST_ASSERT(irrational_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(irrational_solution.solution.has_value());
    check_contains((*irrational_solution.solution)[0u],X(1u,dp));
    check_vector_solution(irrational,b2,*irrational_solution.solution);

    // Two nontrivial columns; the exact solutions are (1,1) and (2,0).
    S matrix(2u,zero);
    matrix[0u][0u]=X(2u,dp);
    matrix[0u][1u]=X(1u,dp);
    matrix[1u][1u]=X(2u,dp);
    Matrix<X> right_hand_sides(2u,2u,zero);
    right_hand_sides[0u][0u]=X(3u,dp);
    right_hand_sides[1u][0u]=X(3u,dp);
    right_hand_sides[0u][1u]=X(4u,dp);
    right_hand_sides[1u][1u]=X(2u,dp);

    const auto matrix_solution=validated_positive_definite_solve(matrix,right_hand_sides);
    ARIADNE_TEST_ASSERT(matrix_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(matrix_solution.solution.has_value());
    check_matrix_solution(matrix,right_hand_sides,*matrix_solution.solution);
    check_contains((*matrix_solution.solution)[0u][0u],X(1u,dp));
    check_contains((*matrix_solution.solution)[1u][0u],X(1u,dp));
    check_contains((*matrix_solution.solution)[0u][1u],X(2u,dp));
    check_contains((*matrix_solution.solution)[1u][1u],zero);

    const Vector<X> b({X(3u,dp),X(3u,dp)});
    const auto vector_solution=validated_positive_definite_solve(matrix,b);
    ARIADNE_TEST_ASSERT(vector_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(vector_solution.solution.has_value());
    check_vector_solution(matrix,b,*vector_solution.solution);
    check_contains((*vector_solution.solution)[0u],X(1u,dp));
    check_contains((*vector_solution.solution)[1u],X(1u,dp));
}

Void test_three_by_three_and_interval_family()
{
    const X zero(0u,dp);
    S matrix(3u,zero);
    for(SizeType i=0; i!=3u; ++i) {
        matrix[i][i]=X(4u,dp);
    }
    matrix[0u][1u]=X(1u,dp);
    matrix[0u][2u]=X(1u,dp);
    matrix[1u][2u]=X(1u,dp);

    const Vector<X> rhs({X(9u,dp),X(12u,dp),X(15u,dp)});
    const auto dense_solution=validated_positive_definite_solve(matrix,rhs);
    ARIADNE_TEST_ASSERT(dense_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(dense_solution.solution.has_value());
    check_vector_solution(matrix,rhs,*dense_solution.solution);
    for(SizeType i=0; i!=3u; ++i) {
        check_contains((*dense_solution.solution)[i],X(static_cast<Int>(i+1u),dp));
    }

    S interval_matrix(2u,zero);
    interval_matrix[0u][0u]=X(4,5,dp);
    interval_matrix[0u][1u]=X(-1,1,dp);
    interval_matrix[1u][1u]=X(4,5,dp);
    Vector<X> interval_rhs(2u,zero);
    interval_rhs[0u]=X(4,5,dp);
    interval_rhs[1u]=X(8,10,dp);

    const auto interval_solution=validated_positive_definite_solve(interval_matrix,interval_rhs);
    ARIADNE_TEST_ASSERT(interval_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(interval_solution.solution.has_value());
    check_vector_solution(interval_matrix,interval_rhs,*interval_solution.solution);

    // Choose one exact member of the family: A=diag(4,4), b=(4,8).
    check_contains((*interval_solution.solution)[0u],X(1u,dp));
    check_contains((*interval_solution.solution)[1u],X(2u,dp));
}

Void test_failure_and_indeterminate()
{
    const X zero(0u,dp);
    S negative(1u,zero);
    negative[0u][0u]=X(-1,dp);
    const Vector<X> b1({X(1u,dp)});
    const auto negative_vector_solution=validated_positive_definite_solve(negative,b1);
    ARIADNE_TEST_ASSERT(negative_vector_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_vector_solution.solution.has_value());

    const Matrix<X> B1(1u,1u,X(1u,dp));
    const auto negative_matrix_solution=validated_positive_definite_solve(negative,B1);
    ARIADNE_TEST_ASSERT(negative_matrix_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_matrix_solution.solution.has_value());

    S uncertain(1u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    const auto uncertain_vector_solution=validated_positive_definite_solve(uncertain,b1);
    ARIADNE_TEST_ASSERT(uncertain_vector_solution.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_vector_solution.solution.has_value());

    const auto uncertain_matrix_solution=validated_positive_definite_solve(uncertain,B1);
    ARIADNE_TEST_ASSERT(uncertain_matrix_solution.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_matrix_solution.solution.has_value());

    // A singular positive-semidefinite matrix cannot use an SPD solve.
    S semidefinite(2u,zero);
    semidefinite[0u][0u]=X(1u,dp);
    semidefinite[0u][1u]=X(1u,dp);
    semidefinite[1u][1u]=X(1u,dp);
    const Vector<X> b2({X(2u,dp),X(2u,dp)});
    const auto semidefinite_solution=validated_positive_definite_solve(semidefinite,b2);
    ARIADNE_TEST_ASSERT(semidefinite_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!semidefinite_solution.solution.has_value());

    const Matrix<X> bad_matrix_rhs(2u,1u,zero);
    const Vector<X> bad_vector_rhs(2u,zero);
    ARIADNE_TEST_FAIL(validated_positive_definite_solve(negative,bad_matrix_rhs));
    ARIADNE_TEST_FAIL(validated_positive_definite_solve(negative,bad_vector_rhs));
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_empty_and_zero_rhs());
    ARIADNE_TEST_CALL(test_scalar_and_dense_spd());
    ARIADNE_TEST_CALL(test_three_by_three_and_interval_family());
    ARIADNE_TEST_CALL(test_failure_and_indeterminate());
    return ARIADNE_TEST_FAILURES;
}
