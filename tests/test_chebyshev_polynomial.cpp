/***************************************************************************
 *            test_chebyshev_polynomial.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <iomanip>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/differential.hpp"
#include "algebra/polynomial.hpp"

#include "algebra/chebyshev_polynomial.hpp"
#include "algebra/chebyshev_polynomial.tpl.hpp"


#include "algebra_testing_context.hpp"

using std::cout; using std::cerr; using std::endl;
using namespace Ariadne;


extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatMP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatMP>::output_places;

template<class I, class X> decltype(auto) sup_norm(Polynomial<I,X> const& p) {
    auto r=mag(p.expansion().zero_coefficient());
    for (auto term : p) { r += mag(term.coefficient()); }
    return r;
}

template<class X> class TestChebyshevPolynomial
{
    typedef PrecisionType<X> PR;
    PR pr;
  public:
    TestChebyshevPolynomial(PrecisionType<X> prec);
    Void test();
  private:
    Void test_univariate();
    Void test_multivariate();
};

template<class X> TestChebyshevPolynomial<X>::TestChebyshevPolynomial(PrecisionType<X> prec)
    : pr(prec)
{
}

template<class X> Void TestChebyshevPolynomial<X>::test()
{
    Float<PR>::set_output_places(18);
    FloatBounds<PR>::set_output_places(18);
    FloatApproximation<PR>::set_output_places(8);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), pow2(0u), 1);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), pow2(3u), 8);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), powm1(0u), 1);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), powm1(1u), -1);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_univariate());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_multivariate());
}

namespace {

template<class X>
constexpr bool check_concept()
{
    return requires(X c, UnivariateChebyshevPolynomial<X> x) {
        x=+x; x=-x; x=x+x; x=x-x; x=x*x;
        x=x+c; x=c+x; x=x-c; x=c-x; x=x*c; x=c*x; x=x/c;
        x=pow(x,2u);
    };
}

static_assert(check_concept<FloatDPApproximation>());
static_assert(check_concept<FloatMPApproximation>());

} // namespace


template<class X> Void TestChebyshevPolynomial<X>::test_univariate() {
    UnivariateChebyshevPolynomial<X> one = UnivariateChebyshevPolynomial<X>::constant(1,pr);
    UnivariateChebyshevPolynomial<X> x = UnivariateChebyshevPolynomial<X>::coordinate(pr);

    std::function<X(DegreeType)> coefficients=[this](DegreeType i) {
        return X(static_cast<Nat>(i)+1u,pr);
    };
    UnivariateChebyshevPolynomial<X> generated(2u,coefficients);
    UnivariateChebyshevPolynomial<X> const& cgenerated=generated;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), cgenerated[0u], X(1u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), cgenerated[2u], X(3u,pr));

    UnivariateChebyshevPolynomial<X> empty(pr);
    auto x_constant=UnivariateChebyshevPolynomial<X>::constant(X(2u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), x_constant[0u], X(2u,pr));
    auto empty_plus_constant=UnivariateChebyshevPolynomial<X>::apply(Add(),empty,X(2u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), empty_plus_constant[0u], X(2u,pr));
    auto one_plus_constant=UnivariateChebyshevPolynomial<X>::apply(Add(),one,X(2u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), one_plus_constant[0u], X(3u,pr));
    auto x_plus_constant=UnivariateChebyshevPolynomial<X>::apply(Add(),x,X(2u,pr));
    auto left_tail=UnivariateChebyshevPolynomial<X>::apply(Add(),UnivariateChebyshevPolynomial<X>::basis(2u,pr),one);
    auto right_tail=UnivariateChebyshevPolynomial<X>::apply(Add(),one,UnivariateChebyshevPolynomial<X>::basis(2u,pr));
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "left_tail", left_tail);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "right_tail", right_tail);

    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "x", x);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "-x", -x);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "x*x", x*x);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "x*x*x", x*x*x);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "x*x*x*x", x*x*x*x);

    X y(-0.75_dy,pr);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), evaluate(empty,y), X(0u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), x_plus_constant(y), y+X(2u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (+x)(y), x(y));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), x(y), (y));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x)(y), (y*y));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x*x)(y), (y*y*y));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x*x*x)(y), (y*y*y*y));

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), UnivariateChebyshevPolynomial<X>::basis(0u,pr), one);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), UnivariateChebyshevPolynomial<X>::basis(1u,pr), x);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), UnivariateChebyshevPolynomial<X>::basis(2u,pr), x*x*2-1);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), UnivariateChebyshevPolynomial<X>::basis(3u,pr), x*x*x*4-x*3);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), UnivariateChebyshevPolynomial<X>::basis(4u,pr), x*x*x*x*8-x*x*8+1);

    UnivariatePolynomial<X> p=UnivariatePolynomial<X>::coordinate(pr);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "p", p);

    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "pnorm(x^4-x^2)=" << sup_norm(p*p*p*p-p*p) << "\n");
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sup_norm(x*x), 1.0_dy);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sup_norm(pow(x,4u)), 1.0_dy);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sup_norm(pow(x,4)-pow(x,2)), 0.25_dy);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sup_norm(x*x*x*x-x*x), 0.25_dy);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sup_norm(3*x*x*x*x-2*x*x-1), 1.75_dy);
    ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), sup_norm(3*x*x*x*x-2*x*x-1), <=,  sup_norm(3*p*p*p*p-2*p*p-1));
//    std::cout << "pnorm(3x^4-2x^2-1)=sup_norm("<<(X(3,pr)*p*p*p*p-X(2,pr)*p*p-X(1,pr)) << ")=" << sup_norm(X(3,pr)*p*p*p*p-X(2,pr)*p*p-X(1,pr)) << "\n";
    // Blank diagnostic omitted.
}


template<class X> Void TestChebyshevPolynomial<X>::test_multivariate() {
    MultivariateChebyshevPolynomial<X> one = MultivariateChebyshevPolynomial<X>::constant(2u,1,pr);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), one.zero_coefficient(), X(0u,pr));
    MultivariateChebyshevPolynomial<X> x = MultivariateChebyshevPolynomial<X>::coordinate(2u,0u,pr);
    MultivariateChebyshevPolynomial<X> y = MultivariateChebyshevPolynomial<X>::coordinate(2u,1u,pr);
    auto multivariate_constant=MultivariateChebyshevPolynomial<X>::constant(2u,X(2u,pr));
    auto multivariate_basis=MultivariateChebyshevPolynomial<X>::basis(2u,0u,2u,pr);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "multivariate_constant", multivariate_constant);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "multivariate_basis", multivariate_basis);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "-x", -x);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), multivariate_constant.sup_norm());

    MultivariateChebyshevPolynomial<X> multivariate_empty(2u,pr);
    auto multivariate_empty_plus=MultivariateChebyshevPolynomial<X>::apply(Add(),multivariate_empty,X(2u,pr));
    auto multivariate_one_plus=MultivariateChebyshevPolynomial<X>::apply(Add(),one,X(2u,pr));
    auto multivariate_x_plus=MultivariateChebyshevPolynomial<X>::apply(Add(),x,X(2u,pr));
    auto multivariate_scaled=MultivariateChebyshevPolynomial<X>::apply(Mul(),x,X(2u,pr));
    auto multivariate_empty_scaled=MultivariateChebyshevPolynomial<X>::apply(Mul(),multivariate_empty,X(2u,pr));
    auto multivariate_empty_product=MultivariateChebyshevPolynomial<X>::apply(Mul(),multivariate_empty,x);
    auto multivariate_zero=MultivariateChebyshevPolynomial<X>::constant(2u,X(0u,pr));
    auto multivariate_zero_product=MultivariateChebyshevPolynomial<X>::apply(Mul(),multivariate_zero,x);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "multivariate_empty_plus", multivariate_empty_plus);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "multivariate_one_plus", multivariate_one_plus);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "multivariate_empty_scaled", multivariate_empty_scaled);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "multivariate_zero_product", multivariate_zero_product);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), multivariate_empty_product.sup_norm(), mag(X(0u,pr)));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), multivariate_zero_product.sup_norm(), mag(X(0u,pr)));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultivariateChebyshevPolynomial<X>::apply(Add(),x,MultivariateChebyshevPolynomial<X>::coordinate(3u,0u,pr)));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultivariateChebyshevPolynomial<X>::apply(Mul(),x,MultivariateChebyshevPolynomial<X>::coordinate(3u,0u,pr)));
    auto xy=MultivariateChebyshevPolynomial<X>::apply(Add(),x,y);
    auto yx=MultivariateChebyshevPolynomial<X>::apply(Add(),y,x);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "xy", xy);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "yx", yx);

    Vector<X> v({0.5,-0.75},pr);
    Vector<X> wrong_dimension(1u,X(0u,pr));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), x(wrong_dimension));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), multivariate_empty(v), X(0u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), multivariate_empty_product(v), X(0u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), x(v), v[0]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), multivariate_x_plus(v), v[0]+X(2u,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), multivariate_scaled(v), X(2u,pr)*v[0]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (+x)(v), x(v));
    auto xx=MultivariateChebyshevPolynomial<X>::apply(Mul(),x,x);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), xx.sup_norm(), mag(X(1u,pr)));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), xx(v), v[0]*v[0]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x*x)(v), v[0]*v[0]*v[0]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x*x*x)(v), v[0]*v[0]*v[0]*v[0]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), pow(x,3)(v), pow(v[0],3));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), pow(x,4)(v), pow(v[0],4));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x*x*x-x*x)(v), v[0]*v[0]*v[0]*v[0]-v[0]*v[0]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (pow(x,4)-pow(x,2))(v), pow(v[0],4)-pow(v[0],2));

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), y(v), v[1]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (y*y)(v), v[1]*v[1]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (y*y*y)(v), v[1]*v[1]*v[1]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (y*y*y*y)(v), v[1]*v[1]*v[1]*v[1]);

    assert(v.size()==2);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), ((x+y)*(x-y))(v), v[0]*v[0]-v[1]*v[1]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*y)(v), v[0]*v[1]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*x*y)(v), v[0]*v[0]*v[1]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), (x*y*y)(v), v[0]*v[1]*v[1]);
}



int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("chebyshev_polynomial");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

    MultiplePrecision mp(128);
    TestChebyshevPolynomial<FloatDPApproximation>(dp).test();
    TestChebyshevPolynomial<FloatMPApproximation>(mp).test();

    

    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}



