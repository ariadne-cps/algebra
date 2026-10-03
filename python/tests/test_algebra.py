#!/usr/bin/python3

from pyariadne import *


def test_multi_index():
    alpha = MultiIndex((2, 1))
    assert len(alpha) == 2
    assert alpha[0] == 2
    assert alpha[1] == 1
    assert alpha.degree() == 3

    e0 = MultiIndex.unit(2, 0)
    assert e0.degree() == 1
    assert e0[0] == 1
    assert e0[1] == 0


def test_expansion():
    dp = DoublePrecision()
    p = FloatDPApproximationExpansion({
        (0, 0): 1,
        (1, 0): 2,
        (0, 1): -1,
        (1, 1): 0.5,
    }, dp)

    p.graded_sort()
    assert p.argument_size() == 2
    assert p.number_of_terms() == 4
    assert isinstance(p[MultiIndex((1, 0))], FloatDPApproximation)


def test_series():
    dp = DoublePrecision()
    s = FloatDPApproximationSeries.exp(FloatDPApproximation(0, dp))
    coefficients = s.coefficients(5)

    assert len(coefficients) == 6
    assert all(isinstance(c, FloatDPApproximation) for c in coefficients)
