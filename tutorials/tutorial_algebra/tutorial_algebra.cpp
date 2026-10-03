#include <iostream>

#include "ariadne-algebra.hpp"

using namespace Ariadne;

namespace {

void section(const char* title) {
    std::cout << "\n== " << title << " ==\n";
}

} // namespace

int main() {
    DoublePrecision dp;

    section("Linear algebra");

    Vector<FloatDP> v({1.0_x,2.0_x},dp);
    Matrix<FloatDP> A({{1.0_x,2.0_x},{3.0_x,4.0_x}},dp);

    std::cout << "v = " << v << "\n";
    std::cout << "A = " << A << "\n";
    std::cout << "transpose(A) = " << transpose(A) << "\n";

    // Exact floating-point operands can produce validated results: the
    // matrix-vector product encloses the rounding error in FloatDPBounds.
    Vector<FloatDPBounds> Av=A*v;
    std::cout << "A*v = " << Av << "\n";

    section("Multi-indices and sparse expansions");

    // A MultiIndex stores the exponents of a multivariate monomial.
    // alpha=(2,1) represents x_0^2*x_1 and has total degree 3.
    MultiIndex alpha({2u,1u});
    MultiIndex e0=MultiIndex::unit(2u,0u);

    std::cout << "alpha = " << alpha << "\n";
    std::cout << "degree(alpha) = " << alpha.degree() << "\n";
    std::cout << "unit multi-index e0 = " << e0 << "\n";

    // Expansion<MultiIndex,X> is a sparse collection of monomial
    // coefficients. Here:
    //
    //   p(x0,x1) = 1 + 2*x0 - x1 + 0.5*x0*x1.
    Expansion<MultiIndex,FloatDPApproximation> p(
        {{{0u,0u},1.0_x},
         {{1u,0u},2.0_x},
         {{0u,1u},-1.0_x},
         {{1u,1u},0.5_x}},
        dp);

    p.graded_sort();
    std::cout << "p = " << p << "\n";
    std::cout << "number of terms = " << p.number_of_terms() << "\n";
    std::cout << "coefficient of x0 = " << p[MultiIndex({1u,0u})] << "\n";

    section("Power-basis polynomials");

    // Polynomial stores coefficients in the usual monomial basis.
    // Unlike Differential, it has no expansion centre: it is a global
    // algebraic polynomial representation.
    using ApproximationPolynomial=MultivariatePolynomial<FloatDPApproximation>;
    auto polynomial_variables=ApproximationPolynomial::variables(2u,dp);
    auto px=polynomial_variables[0];
    auto py=polynomial_variables[1];

    FloatDPApproximation one_a(1,dp);
    FloatDPApproximation two_a(2,dp);
    FloatDPApproximation half_a(0.5_x,dp);

    auto polynomial=ApproximationPolynomial::constant(2u,one_a);
    polynomial += two_a*px;
    polynomial -= py;
    polynomial += half_a*px*py;

    Vector<FloatDPApproximation> polynomial_point({0.25_x,-0.5_x},dp);

    std::cout << "p = 1 + 2*x0 - x1 + 0.5*x0*x1 = " << polynomial << "\n";
    std::cout << "p" << polynomial_point << " = "
              << evaluate(polynomial,polynomial_point) << "\n";
    std::cout << "d/dx0 p = " << derivative(polynomial,0u) << "\n";

    section("Chebyshev polynomials");

    // Chebyshev polynomials represent the same kind of algebraic object
    // in the T_k basis rather than the monomial power basis.
    using UnivariateChebyshev=UnivariateChebyshevPolynomial<FloatDPApproximation>;
    auto chebyshev_x=UnivariateChebyshev::coordinate(dp);
    auto T2=UnivariateChebyshev::basis(2u,dp);
    auto T3=UnivariateChebyshev::basis(3u,dp);
    FloatDPApproximation quarter(0.25_x,dp);

    std::cout << "T2 = " << T2 << "\n";
    std::cout << "2*x^2-1 = " << chebyshev_x*chebyshev_x*2-1 << "\n";
    std::cout << "T3(0.25) = " << T3(quarter) << "\n";

    using MultivariateChebyshev=MultivariateChebyshevPolynomial<FloatDPApproximation>;
    auto chebyshev_x0=MultivariateChebyshev::coordinate(2u,0u,dp);
    auto chebyshev_x1=MultivariateChebyshev::coordinate(2u,1u,dp);

    std::cout << "(T1(x0)*T1(x1))" << polynomial_point << " = "
              << (chebyshev_x0*chebyshev_x1)(polynomial_point) << "\n";

    section("Power series");

    // Series computes Taylor coefficients lazily.  For exp around zero
    // the coefficients are 1, 1, 1/2, 1/6, ...
    FloatDPApproximation zero(0,dp);
    Series<FloatDPApproximation> exp_series(Exp(),zero);

    std::cout << "exp series around 0 = " << exp_series << "\n";
    std::cout << "coefficients through degree 5 = "
              << exp_series.coefficients(5u) << "\n";

    section("Multivariate differential algebra");

    // Differential<X> represents a truncated multivariate jet.  The
    // variables below are centred at (1, 1/2) and retain terms through
    // total degree three.
    Vector<FloatDPBounds> centre({1.0_x,0.5_x},dp);
    auto vars=Differential<FloatDPBounds>::variables(3u,centre);
    auto x=vars[0];
    auto y=vars[1];

    auto f=sqr(x)+x*y+sin(y);

    std::cout << "f = x^2 + x*y + sin(y) at " << centre << "\n";
    std::cout << "jet(f) = " << f << "\n";
    std::cout << "value(f) = " << f.value() << "\n";
    std::cout << "gradient(f) = " << f.gradient() << "\n";
    std::cout << "d/dx f = " << derivative(f,0u) << "\n";
    std::cout << "d/dy f = " << derivative(f,1u) << "\n";

    section("Univariate differential algebra");

    // The univariate form is lighter when there is only one argument.
    // It stores Taylor coefficients through the requested degree while
    // exposing value, first derivative and second derivative directly.
    FloatDPBounds one(1,dp);
    auto t=UnivariateDifferential<FloatDPBounds>::variable(4u,one);
    auto g=exp(t);

    std::cout << "g = exp(t) at t=1: " << g << "\n";
    std::cout << "value(g) = " << g.value() << "\n";
    std::cout << "g'(1) = " << g.gradient() << "\n";
    std::cout << "g''(1) = " << g.hessian() << "\n";
    std::cout << "differential of g' = " << derivative(g) << "\n";

    return 0;
}
