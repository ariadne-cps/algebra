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
    X one(1u,dp);
    X two(2u,dp);
    X three(3u,dp);

    Array<X> coefficients(InitializerList<ExactDouble>{1.0_x,2.0_x,3.0_x},dp);
    D differential(coefficients);

    ARIADNE_TEST_EQUALS(differential.value(),one);
    ARIADNE_TEST_EQUALS(differential.gradient(),two);
    ARIADNE_TEST_EQUALS(differential.half_hessian(),three);
    ARIADNE_TEST_EQUALS(differential.hessian(),X(6u,dp));
    ARIADNE_TEST_EQUALS(differential.argument_size(),1u);

    ARIADNE_TEST_EQUALS(differential.array()[1u],two);
    D const& const_differential=differential;
    ARIADNE_TEST_EQUALS(const_differential.array()[2u],three);
    ARIADNE_TEST_EXECUTE(differential.array()[2u]=X(4u,dp));
    ARIADNE_TEST_EQUALS(differential[2u],X(4u,dp));

    D from_list(2u,InitializerList<X>{one,two,three});
    ARIADNE_TEST_EQUALS(from_list[0u],one);
    ARIADNE_TEST_EQUALS(from_list[1u],two);
    ARIADNE_TEST_EQUALS(from_list[2u],three);
    ARIADNE_TEST_FAIL((void)D(2u,InitializerList<X>{one,two}));

    D degree_zero(0u,one);
    ARIADNE_TEST_FAIL((void)degree_zero.gradient());
    auto degree_zero_variable=D::variable(0u,one);
    ARIADNE_TEST_EQUALS(degree_zero_variable[0u],one);

    D degree_one(1u,one);
    ARIADNE_TEST_FAIL((void)degree_one.half_hessian());
    ARIADNE_TEST_FAIL((void)degree_one.hessian());
}

Void test_inplace_and_calculus()
{
    X one(1u,dp);
    X two(2u,dp);
    X three(3u,dp);

    D differential(2u,InitializerList<X>{one,two,three});
    ARIADNE_TEST_EXECUTE(differential+=one);
    ARIADNE_TEST_EQUALS(differential[0u],two);
    ARIADNE_TEST_EXECUTE(differential*=two);
    ARIADNE_TEST_EQUALS(differential[0u],X(4u,dp));
    ARIADNE_TEST_EQUALS(differential[1u],X(4u,dp));
    ARIADNE_TEST_EQUALS(differential[2u],X(6u,dp));

    D source(2u,InitializerList<X>{one,two,three});
    auto primitive=antiderivative(source);
    ARIADNE_TEST_EQUAL(primitive.degree(),3u);
    ARIADNE_TEST_EQUALS(primitive[0u],X(0u,dp));
    ARIADNE_TEST_EQUALS(primitive[1u],one);
    ARIADNE_TEST_EQUALS(primitive[2u],one);
    ARIADNE_TEST_EQUALS(primitive[3u],one);

    auto shifted_primitive=antiderivative(source,X(4u,dp));
    ARIADNE_TEST_EQUALS(shifted_primitive[0u],X(4u,dp));
    ARIADNE_TEST_EQUALS(shifted_primitive[1u],one);
    ARIADNE_TEST_EQUALS(shifted_primitive[2u],one);
    ARIADNE_TEST_EQUALS(shifted_primitive[3u],one);

    Series<X> exponential_series(Exp(),X(0u,dp));
    D variable=D::variable(2u,X(0u,dp));
    auto composition=compose(exponential_series,variable);
    ARIADNE_TEST_EQUALS(composition[0u],one);
    ARIADNE_TEST_EQUALS(composition[1u],one);
    ARIADNE_TEST_EQUALS(composition[2u],one/two);
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
    ARIADNE_TEST_CALL(test_inplace_and_calculus());
    ARIADNE_TEST_CALL(test_derivative());

    return ARIADNE_TEST_FAILURES;
}
