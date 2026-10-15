/***************************************************************************
 *            test_exact_linear_algebra.cpp
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

#include "macro/testing.hpp"
#include "algebra/exact_linear_algebra.hpp"

using namespace Ariadne;

namespace {

Void test_empty_and_zero_null_spaces()
{
    const Rational z(0);
    {
        Matrix<Rational> A(0u,0u,z);
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,0u);
        ARIADNE_TEST_EQUALS(result.basis.row_size(),0u);
        ARIADNE_TEST_EQUALS(result.basis.column_size(),0u);
    }
    {
        Matrix<Rational> A(0u,3u,z);
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,0u);
        ARIADNE_TEST_EQUAL(result.basis,Matrix<Rational>::identity(3u,A.element_characteristics()));
    }
    {
        Matrix<Rational> A(2u,0u,z);
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,0u);
        ARIADNE_TEST_EQUALS(result.basis.row_size(),0u);
        ARIADNE_TEST_EQUALS(result.basis.column_size(),0u);
    }
    {
        Matrix<Rational> A(2u,3u,z);
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,0u);
        ARIADNE_TEST_EQUAL(result.basis,Matrix<Rational>::identity(3u,A.element_characteristics()));
    }
}

Void test_full_rank_null_spaces()
{
    {
        Matrix<Rational> A({{Rational(1),Rational(0)},{Rational(0),Rational(1)}});
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,2u);
        ARIADNE_TEST_EQUALS(result.basis.row_size(),2u);
        ARIADNE_TEST_EQUALS(result.basis.column_size(),0u);
    }
    {
        // A pivot in the second row forces a row exchange.
        Matrix<Rational> A({{Rational(0),Rational(2)},{Rational(1),Rational(3)}});
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,2u);
        ARIADNE_TEST_EQUALS(result.basis.column_size(),0u);
    }
    {
        // Elimination of a nonzero coefficient below and above a pivot.
        Matrix<Rational> A({{Rational(1),Rational(2)},{Rational(3),Rational(4)}});
        auto result=exact_null_space(A);
        ARIADNE_TEST_EQUALS(result.rank,2u);
        ARIADNE_TEST_EQUALS(result.basis.column_size(),0u);
    }
}

Void test_rank_deficient_null_spaces()
{
    const Rational z(0);
    {
        // First column is free; second column is a pivot; last is free.
        Matrix<Rational> A({{Rational(0),Rational(2),Rational(4)},
                            {Rational(0),Rational(4),Rational(8)}});
        auto result=exact_null_space(A);
        Matrix<Rational> expected({{Rational(1),Rational(0)},
                                   {Rational(0),Rational(-2)},
                                   {Rational(0),Rational(1)}});
        ARIADNE_TEST_EQUALS(result.rank,1u);
        ARIADNE_TEST_EQUAL(result.basis,expected);
        ARIADNE_TEST_EQUAL(A*result.basis,Matrix<Rational>(2u,2u,z));
    }
    {
        // Two independent rows with a one-dimensional exact null space.
        Matrix<Rational> A({{Rational(0),Rational(2),Rational(4)},
                            {Rational(1),Rational(3),Rational(5)}});
        auto result=exact_null_space(A);
        Matrix<Rational> expected({{Rational(1)},{Rational(-2)},{Rational(1)}});
        ARIADNE_TEST_EQUALS(result.rank,2u);
        ARIADNE_TEST_EQUAL(result.basis,expected);
        ARIADNE_TEST_EQUAL(A*result.basis,Matrix<Rational>(2u,1u,z));
    }
    {
        // Exact fractions, rather than numerical zero thresholds.
        Matrix<Rational> A({{Rational(2),Rational(3)}});
        auto result=exact_null_space(A);
        Matrix<Rational> expected({{Rational(-3)/Rational(2)},{Rational(1)}});
        ARIADNE_TEST_EQUALS(result.rank,1u);
        ARIADNE_TEST_EQUAL(result.basis,expected);
        ARIADNE_TEST_EQUAL(A*result.basis,Matrix<Rational>(1u,1u,z));
    }
}

Void test_exact_symmetric_definiteness()
{
    using D=ExactSymmetricDefiniteness;
    const Rational z(0);
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(Matrix<Rational>(0u,0u,z))==D::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(Matrix<Rational>({{Rational(2)}}))==D::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(Matrix<Rational>({{Rational(0)}}))==D::POSITIVE_SEMIDEFINITE);
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(Matrix<Rational>({{Rational(-1)}}))==D::INDEFINITE);

    Matrix<Rational> positive({{Rational(2),Rational(1)},{Rational(1),Rational(2)}});
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(positive)==D::POSITIVE_DEFINITE);

    Matrix<Rational> positive_semidefinite({{Rational(1),Rational(1)},{Rational(1),Rational(1)}});
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(positive_semidefinite)==D::POSITIVE_SEMIDEFINITE);

    Matrix<Rational> zero_then_positive({{Rational(0),Rational(0)},{Rational(0),Rational(3)}});
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(zero_then_positive)==D::POSITIVE_SEMIDEFINITE);

    Matrix<Rational> negative_schur({{Rational(1),Rational(2)},{Rational(2),Rational(1)}});
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(negative_schur)==D::INDEFINITE);

    Matrix<Rational> zero_pivot_nonzero_row({{Rational(0),Rational(1)},{Rational(1),Rational(0)}});
    ARIADNE_TEST_ASSERT(exact_symmetric_definiteness(zero_pivot_nonzero_row)==D::INDEFINITE);

    Matrix<Rational> nonsquare(2u,3u,z);
    ARIADNE_TEST_FAIL(exact_symmetric_definiteness(nonsquare));

    Matrix<Rational> nonsymmetric({{Rational(1),Rational(2)},{Rational(3),Rational(1)}});
    ARIADNE_TEST_FAIL(exact_symmetric_definiteness(nonsymmetric));
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_empty_and_zero_null_spaces());
    ARIADNE_TEST_CALL(test_full_rank_null_spaces());
    ARIADNE_TEST_CALL(test_rank_deficient_null_spaces());
    ARIADNE_TEST_CALL(test_exact_symmetric_definiteness());
    return ARIADNE_TEST_FAILURES;
}
