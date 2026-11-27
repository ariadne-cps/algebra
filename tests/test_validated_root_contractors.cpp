/***************************************************************************
 *            tests/test_validated_root_contractors.cpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/test_validated_root_contractors.cpp
 *  @brief Test interval Newton and Krawczyk image semantics without equation-solving.
 */

#include "algebra/validated_root_contractors.hpp"
#include "algebra_testing_context.hpp"
#include "numeric/numeric.hpp"

using namespace Ariadne;

namespace {

//! \brief Require interval inclusion of a known real integer.
template<class PR>
bool encloses_integer(const FloatBounds<PR>& value, Int integer, PR precision) {
    const FloatBounds<PR> exact(integer, precision);
    return value.lower_raw() <= exact.lower_raw() &&
           value.upper_raw() >= exact.upper_raw();
}

//! \brief Check both image maps on a coupled two-variable linear problem.
template<class PR>
Void test_linear_problem(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);
    const X two_value(2u, precision);
    const X three(3u, precision);

    // f(x) = [2*x0+x1-4, x0+2*x1-5]; exact root (1,2).
    Matrix<X> jacobian(2u, 2u, zero);
    jacobian[0u][0u] = two_value;
    jacobian[0u][1u] = one;
    jacobian[1u][0u] = one;
    jacobian[1u][1u] = two_value;

    Vector<X> center(2u, zero);
    Vector<X> residual(2u, zero);
    residual[0u] = X(-4, precision);
    residual[1u] = X(-5, precision);

    Vector<X> domain(2u, X(-1, 4, precision));
    Matrix<X> preconditioner(2u, 2u, zero);
    preconditioner[0u][0u] = two_value / three;
    preconditioner[0u][1u] = -one / three;
    preconditioner[1u][0u] = -one / three;
    preconditioner[1u][1u] = two_value / three;

    const Vector<X> newton = interval_newton_image(center, residual, jacobian);
    const Vector<X> krawczyk = krawczyk_image(center, residual, jacobian,
                                              domain, preconditioner);

    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     encloses_integer(newton[0u], 1, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     encloses_integer(newton[1u], 2, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     encloses_integer(krawczyk[0u], 1, precision));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     encloses_integer(krawczyk[1u], 2, precision));

    // The Krawczyk map, unlike interval Newton, accepts an arbitrary
    // preconditioner. C=0 cannot prove a root: its image is the input box.
    Matrix<X> ineffective(2u, 2u, zero);
    const Vector<X> trivial_image = krawczyk_image(center, residual, jacobian,
                                                  domain, ineffective);
    for (SizeType i = 0u; i != 2u; ++i) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            trivial_image[i].lower_raw() <= domain[i].lower_raw() &&
            trivial_image[i].upper_raw() >= domain[i].upper_raw());
    }
}

//! \brief An interval-Jacobian example with independently known irrational root.
template<class PR>
Void test_nonlinear_problem(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);
    const X one(1u, precision);
    const X two_value(2u, precision);
    const X four(4u, precision);

    // f(x)=x^2-2 with X=[1,2], x0=3/2, f(x0)=1/4 and J(X)=[2,4].
    // Every output must contain sqrt(2), but the image by itself does not
    // assert existence or uniqueness.
    Vector<X> center(1u, X(3u, precision) / two_value);
    Vector<X> residual(1u, one / four);
    Vector<X> box(1u, X(1, 2, precision));
    Matrix<X> jacobian(1u, 1u, zero);
    jacobian[0u][0u] = X(2, 4, precision);
    Matrix<X> preconditioner(1u, 1u, zero);
    preconditioner[0u][0u] = one / X(3u, precision);

    const auto newton = interval_newton_image(center, residual, jacobian);
    const auto krawczyk = krawczyk_image(center, residual, jacobian,
                                          box, preconditioner);
    const X radical = sqrt(two_value);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        newton[0u].lower_raw() <= radical.lower_raw() &&
        newton[0u].upper_raw() >= radical.upper_raw());
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        krawczyk[0u].lower_raw() <= radical.lower_raw() &&
        krawczyk[0u].upper_raw() >= radical.upper_raw());
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("validated_root_contractors");
    suite.add("linear and nonlinear interval images", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_linear_problem(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_linear_problem(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_nonlinear_problem(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_nonlinear_problem(MultiplePrecision(80u)));
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
