/***************************************************************************
 *            algebra/chebyshev_polynomial.hpp
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
 *  GNU General Public License for more details.
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

/** @file algebra/chebyshev_polynomial.hpp
 *  @brief Chebyshev polynomials and their algebraic operations.
 */

#ifndef ARIADNE_CHEBYSHEV_POLYNOMIAL_HPP
#define ARIADNE_CHEBYSHEV_POLYNOMIAL_HPP

#include <functional>
#include <map>

#include "algebra/covector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/operations.hpp"
#include "algebra/vector.hpp"
#include "utility/array.hpp"
#include "utility/macros.hpp"
#include "utility/pointer.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

struct IndexLess;

template<class X>
class UnivariateChebyshevPolynomial;
template<class X>
class MultivariateChebyshevPolynomial;

template<class X>
class UnivariateChebyshevPolynomial
    : public DispatchAlgebraOperations<UnivariateChebyshevPolynomial<X>, X> {
    typedef typename X::Paradigm P;
    typedef typename X::PrecisionType PR;
    SortedExpansion<UniIndex, X, IndexLess> _terms;

  public:
    template<class XX>
    using ChebyshevPolynomial = UnivariateChebyshevPolynomial<XX>;

    typedef P Paradigm;
    typedef PR PrecisionType;

    typedef X NumericType;

    UnivariateChebyshevPolynomial(PR pr) : _terms(SizeOne(), pr) {}
    UnivariateChebyshevPolynomial(DegreeType d, const std::function<X(DegreeType)>& g) : _terms(SizeOne(), get_characteristics(g(0))) {
        for (DegreeType i = 0; i <= d; ++i) {
            _terms.append(i, g(i));
        }
    }

    static ChebyshevPolynomial<X> constant(Number<P> y, PR pr);
    static ChebyshevPolynomial<X> constant(const X& c);
    static ChebyshevPolynomial<X> coordinate(PR pr);
    static ChebyshevPolynomial<X> basis(SizeType k, PR pr);

    UnivariateChebyshevPolynomial<X> create_constant(Number<P> y);

    DegreeType degree() const { return _terms.back().index(); }
    PR precision() const { return _terms.zero_coefficient().precision(); }
    const X& zero_coefficient() const { return _terms.zero_coefficient(); }
    const X& operator[](SizeType i) const { return _terms[UniIndex(static_cast<DegreeType>(i))]; }

    MagType<X> sup_norm() const;
    friend MagType<X> sup_norm(ChebyshevPolynomial<X> cm) { return cm.sup_norm(); }

    friend decltype(auto) operator==(const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2) {
        return cm1._terms == cm2._terms;
    }

  private: // FIXME: Put these concrete-generic operations in proper place
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator+(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p + xc;
    }
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator-(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p - xc;
    }
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator*(const Y& c, ChebyshevPolynomial<X> p) {
        X xc = p.zero_coefficient();
        xc = c;
        return xc * p;
    }
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator*(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p * xc;
    }
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator/(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p / xc;
    }

  public:
    friend OutputStream& operator<<(OutputStream& os, const ChebyshevPolynomial<X>& cm) { return cm._write(os); }

    X operator()(const X& x) const;
    template<class Y>
    friend decltype(auto) evaluate(const ChebyshevPolynomial<X>& p, const Y& y) { return _evaluate(p, y); }

  private:

  public:
    static ChebyshevPolynomial<X> apply(Pos, ChebyshevPolynomial<X> cm);
    static ChebyshevPolynomial<X> apply(Neg, ChebyshevPolynomial<X> cm);
    static ChebyshevPolynomial<X> apply(Add, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    static ChebyshevPolynomial<X> apply(Sub, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    static ChebyshevPolynomial<X> apply(Mul, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    static ChebyshevPolynomial<X> apply(Add, ChebyshevPolynomial<X> cm1, const X& s2);
    static ChebyshevPolynomial<X> apply(Mul, ChebyshevPolynomial<X> cm1, const X& s2);
    static ChebyshevPolynomial<X> apply(Sqr, ChebyshevPolynomial<X> cm);

    OutputStream& _write(OutputStream& os) const;

  public:
    template<class Y>
    static auto _evaluate(const ChebyshevPolynomial<X>& f, const Y& v) -> ArithmeticType<X, Y>;
};

template<class X>
struct AlgebraOperations<UnivariateChebyshevPolynomial<X>, X> : public UnivariateChebyshevPolynomial<X> {};

template<class X>
template<class Y>
auto UnivariateChebyshevPolynomial<X>::_evaluate(const UnivariateChebyshevPolynomial<X>& f, const Y& x) -> ArithmeticType<X, Y> {
    typedef ArithmeticType<X, Y> R;
    R r = x * 0;
    if (f._terms.empty()) {
        return r;
    }
    R p = x * 0 + 1;
    R c(x);
    R n(x);
    Y tx = x * 2;
    SizeType k = 1;
    auto iter = f._terms.begin();
    if (iter->index() == 0) {
        r = iter->coefficient();
        ++iter;
    }
    while (iter != f._terms.end()) {
        while (k < iter->index()) {
            n = tx * c - p;
            p = c;
            c = n;
            ++k;
        }
        r += c * iter->coefficient();
        ++iter;
    }
    return r;
}

template<class X>
class MultivariateChebyshevPolynomial
    : public DispatchAlgebraOperations<MultivariateChebyshevPolynomial<X>, X> {
    typedef typename X::Paradigm P;
    typedef typename X::PrecisionType PR;
    Expansion<MultiIndex, X> _terms;

  public:
    template<class XX>
    using ChebyshevPolynomial = MultivariateChebyshevPolynomial<XX>;

    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef X NumericType;
    typedef typename Expansion<MultiIndex, X>::ConstIterator ConstIterator;

    MultivariateChebyshevPolynomial(SizeType as, PR pr) : _terms(as, pr) {}

    static ChebyshevPolynomial<X> constant(SizeType as, Number<P> y, PR pr);
    static ChebyshevPolynomial<X> constant(SizeType as, const X& c);
    static ChebyshevPolynomial<X> coordinate(SizeType as, SizeType i, PR pr);
    static ChebyshevPolynomial<X> basis(SizeType as, SizeType i, SizeType k, PR pr);

    ChebyshevPolynomial<X> create_constant(Number<P> y);

    SizeType argument_size() const { return this->_terms.argument_size(); }
    PR precision() const { return this->_terms.zero_coefficient().precision(); }
    const X& zero_coefficient() const { return _terms.zero_coefficient(); }

    MagType<X> sup_norm() const;
    friend MagType<X> sup_norm(const ChebyshevPolynomial<X>& cm) { return cm.sup_norm(); }

    friend ChebyshevPolynomial<X> operator*(const Number<P>& s1, ChebyshevPolynomial<X> cm2) { return X(s1, cm2.precision()) * cm2; }
    friend OutputStream& operator<<(OutputStream& os, const ChebyshevPolynomial<X>& cm) { return cm._write(os); }

    X operator()(const Vector<X>& x) const;
    friend X evaluate(const ChebyshevPolynomial<X>& f, const Vector<X>& x) { return f(x); }

  private:

  public:
    static ChebyshevPolynomial<X> apply(Pos, ChebyshevPolynomial<X> cm);
    static ChebyshevPolynomial<X> apply(Neg, ChebyshevPolynomial<X> cm);
    static ChebyshevPolynomial<X> apply(Sqr, ChebyshevPolynomial<X> cm);
    static ChebyshevPolynomial<X> apply(Add, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    static ChebyshevPolynomial<X> apply(Sub, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    static ChebyshevPolynomial<X> apply(Mul, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    static ChebyshevPolynomial<X> apply(Add, ChebyshevPolynomial<X> cm1, const Scalar<X>& s2);
    static ChebyshevPolynomial<X> apply(Mul, ChebyshevPolynomial<X> cm1, const Scalar<X>& s2);
    OutputStream& _write(OutputStream& os) const;
};

template<class X>
struct AlgebraOperations<MultivariateChebyshevPolynomial<X>, X> : public MultivariateChebyshevPolynomial<X> {};

} // namespace Ariadne

#endif // ARIADNE_CHEBYSHEV_POLYNOMIAL_HPP
