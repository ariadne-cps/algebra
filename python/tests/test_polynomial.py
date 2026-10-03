#!/usr/bin/python3

from pyariadne import *


def test_multivariate_polynomial():
    dp = DoublePrecision()
    Polynomial = FloatDPApproximationMultivariatePolynomial

    assert MultivariatePolynomial[FloatDPApproximation] == Polynomial

    variables = Polynomial.variables(2, dp)
    x = variables[0]
    y = variables[1]

    one = FloatDPApproximation(1, dp)
    p = Polynomial.constant(2, one) + x*x + y

    assert p.argument_size() == 2
    assert p.degree() == 2
    assert p.number_of_terms() == 3

    point = FloatDPApproximationVector([exact(0.25), exact(-0.5)], dp)
    value = p(point)
    assert isinstance(value, FloatDPApproximation)

    dpdx = derivative(p, 0)
    assert dpdx.degree() == 1


def test_univariate_chebyshev_polynomial():
    dp = DoublePrecision()
    Chebyshev = FloatDPApproximationUnivariateChebyshevPolynomial

    x = Chebyshev.coordinate(dp)
    T2 = Chebyshev.basis(2, dp)
    T3 = Chebyshev.basis(3, dp)

    assert T2.degree() == 2
    assert T3.degree() == 3

    point = FloatDPApproximation(exact(0.25), dp)
    assert isinstance(x(point), FloatDPApproximation)
    assert isinstance(T3(point), FloatDPApproximation)


def test_multivariate_chebyshev_polynomial():
    dp = DoublePrecision()
    Chebyshev = FloatDPApproximationMultivariateChebyshevPolynomial

    x = Chebyshev.coordinate(2, 0, dp)
    y = Chebyshev.coordinate(2, 1, dp)

    assert x.argument_size() == 2
    assert y.argument_size() == 2

    point = FloatDPApproximationVector([exact(0.25), exact(-0.5)], dp)
    value = (x*y)(point)
    assert isinstance(value, FloatDPApproximation)
