#include <iostream>

#include "ariadne-algebra.hpp"

using namespace Ariadne;

int main() {
    DoublePrecision dp;

    Vector<FloatDP> v({1.0_x,2.0_x},dp);
    Matrix<FloatDP> A({{1.0_x,2.0_x},{3.0_x,4.0_x}},dp);
    std::cout << "v = " << v << "\n";
    std::cout << "A = " << A << "\n";
    Vector<FloatDPBounds> Av=A*v;
    std::cout << "A*v = " << Av << "\n";

    FloatDPBounds one(1,dp);
    auto x=Differential<FloatDPBounds>::variable(1u,3u,one,0u);
    auto y=sqr(x)+one;
    std::cout << "y = x^2+1 at x=1: " << y << "\n";
    std::cout << "value = " << y.value() << "\n";
    std::cout << "gradient = " << y.gradient() << "\n";

    return 0;
}
