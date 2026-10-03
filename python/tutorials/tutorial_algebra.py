#!/usr/bin/python3

from pyariadne import *


def section(title):
    print(f"\n== {title} ==")


def show(label, value):
    print(f"{label}: {value}")


def tutorial_algebra():
    dp = DoublePrecision()

    section("Linear algebra")

    v = FloatDPVector([1, 2], dp)
    A = FloatDPMatrix([[1, 2], [3, 4]], dp)

    show("v", v)
    show("A", A)
    show("transpose(A)", transpose(A))

    # As in C++, exact floating-point operands may produce a validated
    # result which encloses rounding error.
    show("A*v", A*v)

    section("Multi-indices and sparse expansions")

    # alpha=(2,1) represents x_0^2*x_1 and has total degree three.
    alpha = MultiIndex((2, 1))
    e0 = MultiIndex.unit(2, 0)

    show("alpha", alpha)
    show("degree(alpha)", alpha.degree())
    show("unit multi-index e0", e0)

    # p(x0,x1) = 1 + 2*x0 - x1 + 0.5*x0*x1.
    p = FloatDPApproximationExpansion({
        (0, 0): 1,
        (1, 0): 2,
        (0, 1): -1,
        (1, 1): 0.5,
    }, dp)

    p.graded_sort()
    show("p", p)
    show("number of terms", p.number_of_terms())
    show("coefficient of x0", p[MultiIndex((1, 0))])

    section("Power series")

    # Series computes Taylor coefficients lazily.
    zero = FloatDPApproximation(0, dp)
    exp_series = FloatDPApproximationSeries.exp(zero)

    show("exp series around 0", exp_series)
    show("coefficients through degree 5", exp_series.coefficients(5))

    section("Multivariate differential algebra")

    centre = FloatDPBoundsVector([exact(1), exact(0.5)], dp)
    variables = FloatDPBoundsDifferential.variables(3, centre)
    x = variables[0]
    y = variables[1]

    f = sqr(x) + x*y + sin(y)

    show(f"f = x^2 + x*y + sin(y) at {centre}", f)
    show("jet(f)", f)
    show("value(f)", f.value())
    show("gradient(f)", f.gradient())
    show("d/dx f", derivative(f, 0))
    show("d/dy f", derivative(f, 1))

    section("Univariate differential algebra")

    one = FloatDPBounds(1, dp)
    t = FloatDPBoundsUnivariateDifferential.variable(4, one)
    g = exp(t)

    show("g = exp(t) at t=1", g)
    show("value(g)", g.value())
    show("g'(1)", g.gradient())
    show("g''(1)", g.hessian())
    show("differential of g'", derivative(g))


if __name__ == "__main__":
    tutorial_algebra()
