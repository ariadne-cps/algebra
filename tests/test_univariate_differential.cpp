/***************************************************************************
 *            test_univariate_differential.cpp
 *
 *  Copyright  2026  Pieter Collins
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
#include "algebra/differential.hpp"
#include "algebra/univariate_differential.tpl.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

namespace {

using X=FloatDPApproximation;
using D=UnivariateDifferential<X>;

Void test_data_access()
{
    Array<X> coefficients(InitializerList<ExactDouble>{1.0_x,2.0_x,3.0_x},dp);
    D differential(coefficients);

    ARIADNE_TEST_EQUALS(differential.value(),X(1u,dp));
    ARIADNE_TEST_EQUALS(differential.gradient(),X(2u,dp));
    ARIADNE_TEST_EQUALS(differential.half_hessian(),X(3u,dp));
    ARIADNE_TEST_EQUALS(differential.hessian(),X(6u,dp));

    D degree_zero(0u,X(1u,dp));
    ARIADNE_TEST_FAIL((void)degree_zero.gradient());

    D degree_one(1u,X(1u,dp));
    ARIADNE_TEST_FAIL((void)degree_one.half_hessian());
    ARIADNE_TEST_FAIL((void)degree_one.hessian());
}

Void test_derivative()
{
    auto x=D::variable(4u,X(1u,dp));
    auto y=exp(x);
    auto dy=derivative(y);

    ARIADNE_TEST_EQUAL(dy.degree(),3u);
    for(DegreeType i=0; i<=dy.degree(); ++i) {
        ARIADNE_TEST_EQUAL(dy[i],(i+1u)*y[i+1u]);
    }

    auto constant=D::constant(0u,X(2u,dp));
    auto constant_derivative=derivative(constant);
    ARIADNE_TEST_EQUAL(constant_derivative.degree(),0u);
    ARIADNE_TEST_EQUALS(constant_derivative[0],0);
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_data_access());
    ARIADNE_TEST_CALL(test_derivative());

    return ARIADNE_TEST_FAILURES;
}
