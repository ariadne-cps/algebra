/***************************************************************************
 *            test_univariate_differential.cpp
 *
 *  Copyright  2026  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "algebra/differential.hpp"
#include "algebra/univariate_differential.tpl.hpp"
#include "algebra_testing_context.hpp"
#include "numeric/numeric.hpp"

using namespace Ariadne;

namespace {

using X = FloatDPApproximation;
using D = UnivariateDifferential<X>;

Void test_data_access() {
    X one(1u, dp);
    X coefficient_two(2u, dp);
    X three(3u, dp);

    Array<X> coefficients(InitializerList<ExactDouble>{1.0_x, 2.0_x, 3.0_x}, dp);
    D differential(coefficients);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential.value(), one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential.gradient(), coefficient_two);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential.half_hessian(), three);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential.hessian(), X(6u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential.argument_size(), 1u);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential.array()[1u], coefficient_two);
    const D& const_differential = differential;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), const_differential.array()[2u], three);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), differential.array()[2u] = X(4u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential[2u], X(4u, dp));

    D from_list(2u, InitializerList<X>{one, coefficient_two, three});
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), from_list[0u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), from_list[1u], coefficient_two);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), from_list[2u], three);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)D(2u, InitializerList<X>{one, coefficient_two}));

    D degree_zero(0u, one);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)degree_zero.gradient());
    auto degree_zero_variable = D::variable(0u, one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), degree_zero_variable[0u], one);

    D degree_one(1u, one);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)degree_one.half_hessian());
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)degree_one.hessian());
}

Void test_inplace_and_calculus() {
    X one(1u, dp);
    X coefficient_two(2u, dp);
    X three(3u, dp);

    D differential(2u, InitializerList<X>{one, coefficient_two, three});
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), differential += one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential[0u], coefficient_two);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), differential *= coefficient_two);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential[0u], X(4u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential[1u], X(4u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), differential[2u], X(6u, dp));

    D source(2u, InitializerList<X>{one, coefficient_two, three});
    auto primitive = antiderivative(source);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), primitive.degree(), 3u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), primitive[0u], X(0u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), primitive[1u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), primitive[2u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), primitive[3u], one);

    auto shifted_primitive = antiderivative(source, X(4u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), shifted_primitive[0u], X(4u, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), shifted_primitive[1u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), shifted_primitive[2u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), shifted_primitive[3u], one);

    Series<X> exponential_series(Exp(), X(0u, dp));
    D variable = D::variable(2u, X(0u, dp));
    auto composition = compose(exponential_series, variable);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), composition[0u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), composition[1u], one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), composition[2u], one / coefficient_two);
}

Void test_derivative() {
    auto x = D::variable(4u, X(1u, dp));
    auto y = exp(x);
    auto dy = derivative(y);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), dy.degree(), 3u);
    for (DegreeType i = 0; i <= dy.degree(); ++i) {
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), dy[i], (i + 1u) * y[i + 1u]);
    }

    auto constant = D::constant(0u, X(2u, dp));
    auto constant_derivative = derivative(constant);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), constant_derivative.degree(), 0u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), constant_derivative[0], 0);
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("univariate_differential");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_data_access());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_inplace_and_calculus());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_derivative());
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
