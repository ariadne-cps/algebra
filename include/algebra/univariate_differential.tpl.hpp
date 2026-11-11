/***************************************************************************
 *            algebra/univariate_differential.tpl.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/univariate_differential.tpl.hpp
 *  @brief Template implementations for univariate differentials.
 */

#include "algebra/univariate_differential.hpp"

namespace Ariadne {

template<class X>
UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, const X& c)
    : _ary(d + 1u, nul(c)) { _ary[0] = c; }

template<class X>
UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, InitializerList<X> lst)
    : _ary(d + 1u, nul(*lst.begin())) {
    ARIADNE_PRECONDITION(lst.size() == d + 1u);
    std::copy(lst.begin(), lst.end(), _ary.begin());
}

template<class X>
UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, const Series<X>& s)
    : _ary(d + 1u, nul(s[0])) {
    for (DegreeType i = 0; i <= d; ++i) {
        this->_ary[i] = s[i];
    }
}

template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::constant(DegreeType d, const X& c) {
    return UnivariateDifferential(d, c);
}

template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::variable(DegreeType d, const X& c) {
    UnivariateDifferential r(d, c);
    if (d >= 1) {
        r[1] = 1;
    }
    return r;
}

template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::create_zero() const {
    return UnivariateDifferential<X>(this->degree(), this->zero_coefficient());
}

template<class X>
UnivariateDifferential<X>
    UnivariateDifferential<X>::create_constant(const NumericType& c) const {
    return UnivariateDifferential<X>(this->degree(), c);
}

template<class X>
UnivariateDifferential<X>
    UnivariateDifferential<X>::create_variable(const NumericType& c) const {
    return UnivariateDifferential<X>::variable(this->degree(), c);
}

template<class X>
UnivariateDifferential<X>&
    UnivariateDifferential<X>::operator=(const NumericType& c) {
    X zero = nul(c);
    for (SizeType i = 0; i != this->_ary.size(); ++i) {
        this->_ary[i] = zero;
    }
    this->_ary[0u] = c;
    return *this;
}

template<class X>
Void UnivariateDifferential<X>::clear() {
    X zero = this->zero_coefficient();
    for (SizeType i = 0; i != this->_ary.size(); ++i) {
        this->_ary[i] = zero;
    }
}

// Compute the truncated reciprocal from f*g=1, coefficient by coefficient.
template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::apply(Rec) const {
    UnivariateDifferential<X> result(this->degree(), this->zero_coefficient());
    result[0u] = rec(this->_ary[0u]);
    for (DegreeType n = 1u; n <= this->degree(); ++n) {
        X sum = this->zero_coefficient();
        for (DegreeType k = 1u; k <= n; ++k) {
            sum += this->_ary[k] * result[n - k];
        }
        result[n] = -result[0u] * sum;
    }
    return result;
}

template<class X>
Differential<X> UnivariateDifferential<X>::_compose(
    const Series<X>& f, const Differential<X>& dx) {
    return Differential<X>::_compose(f, dx);
}

template<class X>
X UnivariateDifferential<X>::zero_coefficient() const {
    return nul(this->_ary[0]);
}

template<class X>
SizeType UnivariateDifferential<X>::argument_size() const {
    return 1u;
}

template<class X>
DegreeType UnivariateDifferential<X>::degree() const {
    return static_cast<DegreeType>(this->_ary.size() - 1u);
}

template<class X>
const Array<X>& UnivariateDifferential<X>::array() const {
    return this->_ary;
}

template<class X>
Array<X>& UnivariateDifferential<X>::array() {
    return this->_ary;
}

template<class X>
const X& UnivariateDifferential<X>::operator[](SizeType k) const {
    return this->_ary[k];
}

template<class X>
X& UnivariateDifferential<X>::operator[](SizeType k) {
    return this->_ary[k];
}

template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::_compose(const Series<X>& x, const UnivariateDifferential<X>& y) {
    return compose(UnivariateDifferential(y.degree(), x), y);
}

template<class X>
UnivariateDifferential<X> UnivariateDifferential<X>::_compose(const UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    DegreeType d = std::min(x.degree(), y.degree());
    X y0 = y[0];
    X z = Ariadne::create_zero(y0);
    const_cast<X&>(y[0]) = z;
    UnivariateDifferential<X> r(d, z);
    r[0] = x[d];
    for (DegreeType n = 1; n <= d; ++n) {
        r = r * y;
        r[0] += x[static_cast<SizeType>(d - n)];
    }
    const_cast<X&>(y[0]) = y0;
    return r;
}

template<class X>
UnivariateDifferential<X>
    UnivariateDifferential<X>::_derivative(const UnivariateDifferential<X>& x) {
    DegreeType n = x.degree();
    UnivariateDifferential<X> r(n == 0u ? 0u : n - 1u, x.zero_coefficient());
    if (n == 0u) {
        r[0] = x[0] * 0;
        return r;
    }
    for (DegreeType i = 0; i < n; ++i) {
        r[i] = (i + 1u) * x[i + 1u];
    }
    return r;
}

template<class X>
UnivariateDifferential<X>
    UnivariateDifferential<X>::_antiderivative(const UnivariateDifferential<X>& x) {
    DegreeType n = x.degree();
    UnivariateDifferential<X> r(n + 1u, nul(x[0]));
    for (DegreeType i = 0; i <= n; ++i) {
        r[i + 1u] = x[i] / (i + 1u);
    }
    return r;
}

template<class X>
UnivariateDifferential<X>
    UnivariateDifferential<X>::_antiderivative(const UnivariateDifferential<X>& x, const X& c) {
    DegreeType n = x.degree();
    UnivariateDifferential<X> r(n + 1u, Ariadne::create_zero(x[0]));
    r[0] = c;
    for (DegreeType i = 0; i <= n; ++i) {
        r[i + 1u] = x[i] / (i + 1u);
    }
    return r;
}

template<class X>
OutputStream& UnivariateDifferential<X>::_write(OutputStream& os) const {
    const UnivariateDifferential<X>& x = *this;
    os << "D<" << x.degree() << ">";
    for (DegreeType i = 0; i <= x.degree(); ++i) {
        os << (i == 0 ? '[' : ',') << x[i];
    }
    return os << "]";
}

/*
template<class X> OutputStream& UnivariateDifferential<X>::_write(OutputStream& os) const {
    os << this->_ary[0];
    for(DegreeType i=1; i<=this->degree(); ++i) {
        os <<" ";
        if(decide(this->_ary[i]>=X(0))) { os << "+"; }
        os << this->_ary[i] << "*dx^" << (uint)i;
    }
    return os;
}
*/

} // namespace Ariadne
