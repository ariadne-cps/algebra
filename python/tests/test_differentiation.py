#!/usr/bin/python3

from pyariadne import *


def test_differential():
    dp = DoublePrecision()
    one = FloatDPBounds(1, dp)
    x = FloatDPBoundsDifferential.variable(1, 3, one, 0)
    y = sqr(x) + one

    assert y.degree() == 3
    assert definitely(y.value() == FloatDPBounds(2, dp))
    assert y.gradient().size() == 1

    dy = derivative(y, 0)
    assert dy.degree() == 2


def test_univariate_differential():
    dp = DoublePrecision()
    x = FloatDPBoundsUnivariateDifferential.variable(3, FloatDPBounds(1, dp))
    assert x.argument_size() == 1
    y = exp(x)
    assert y.degree() == 3
    expected = exp(FloatDPBounds(1, dp))
    assert definitely(y.value() > FloatDPBounds(2, dp))
    assert definitely(y.value() < FloatDPBounds(3, dp))
    assert possibly(y.value() == expected)

    dy = derivative(y)
    assert dy.degree() == 2
