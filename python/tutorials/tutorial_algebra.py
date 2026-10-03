#!/usr/bin/python3

from pyariadne import *


def show(label, value):
    print(f"{label}: {value}")


def tutorial_algebra():
    dp = DoublePrecision()

    v = FloatDPVector([1, 2], dp)
    A = FloatDPMatrix([[1, 2], [3, 4]], dp)
    show("v", v)
    show("A", A)
    show("A*v", A*v)

    one = FloatDPBounds(1, dp)
    x = FloatDPBoundsDifferential.variable(1, 3, one, 0)
    y = sqr(x) + one
    show("y = x^2+1 at x=1", y)
    show("value", y.value())
    show("gradient", y.gradient())
    show("derivative", derivative(y, 0))


if __name__ == "__main__":
    tutorial_algebra()
