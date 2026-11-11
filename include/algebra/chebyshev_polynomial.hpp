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

/** @brief Univariate polynomial in the Chebyshev basis. */
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

    /** @brief Construct a \c UnivariateChebyshevPolynomial from pr. */
    UnivariateChebyshevPolynomial(PR pr) : _terms(SizeOne(), pr) {}
    /** @brief Construct a \c UnivariateChebyshevPolynomial from d, g. */
    UnivariateChebyshevPolynomial(DegreeType d, const std::function<X(DegreeType)>& g) : _terms(SizeOne(), get_characteristics(g(0))) {
        for (DegreeType i = 0; i <= d; ++i) {
            _terms.append(i, g(i));
        }
    }

    /** @brief Construct a constant with the supplied value. */
    static ChebyshevPolynomial<X> constant(Number<P> y, PR pr);
    /** @brief Construct a constant with the supplied value. */
    static ChebyshevPolynomial<X> constant(const X& c);
    /** @brief Construct or return the selected coordinate. */
    static ChebyshevPolynomial<X> coordinate(PR pr);
    /** @brief Return the selected basis element. */
    static ChebyshevPolynomial<X> basis(SizeType k, PR pr);

    /** @brief Construct a constant object with the supplied value. */
    UnivariateChebyshevPolynomial<X> create_constant(Number<P> y);

    /** @brief Return the polynomial or expansion degree. */
    DegreeType degree() const { return _terms.back().index(); }
    /** @brief Return the working precision. */
    PR precision() const { return _terms.zero_coefficient().precision(); }
    /** @brief Return a zero with the coefficient characteristics. */
    const X& zero_coefficient() const { return _terms.zero_coefficient(); }
    /** @brief Return the selected element or component. */
    const X& operator[](SizeType i) const { return _terms[UniIndex(static_cast<DegreeType>(i))]; }

    /** @brief Return the supremum norm. */
    MagType<X> sup_norm() const;
    /** @brief Return the supremum norm. */
    friend MagType<X> sup_norm(ChebyshevPolynomial<X> cm) { return cm.sup_norm(); }

    /** @brief Test equality of the operands. */
    friend decltype(auto) operator==(const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2) {
        return cm1._terms == cm2._terms;
    }

  private: // FIXME: Put these concrete-generic operations in proper place
    /** @brief Return the sum or unary positive value. */
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator+(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p + xc;
    }
    /** @brief Return the difference or additive inverse. */
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator-(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p - xc;
    }
    /** @brief Return the product. */
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator*(const Y& c, ChebyshevPolynomial<X> p) {
        X xc = p.zero_coefficient();
        xc = c;
        return xc * p;
    }
    /** @brief Return the product. */
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator*(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p * xc;
    }
    /** @brief Return the quotient. */
    template<AssignableTo<X> Y>
    friend ChebyshevPolynomial<X> operator/(ChebyshevPolynomial<X> p, const Y& c) {
        X xc = p.zero_coefficient();
        xc = c;
        return p / xc;
    }

  public:
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ChebyshevPolynomial<X>& cm) { return cm._write(os); }

    /** @brief Apply the callable to the supplied arguments. */
    X operator()(const X& x) const;
    /** @brief Evaluate at the supplied argument. */
    template<class Y>
    friend decltype(auto) evaluate(const ChebyshevPolynomial<X>& p, const Y& y) { return _evaluate(p, y); }

  private:

  public:
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Pos, ChebyshevPolynomial<X> cm);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Neg, ChebyshevPolynomial<X> cm);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Add, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Sub, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Mul, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Add, ChebyshevPolynomial<X> cm1, const X& s2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Mul, ChebyshevPolynomial<X> cm1, const X& s2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Sqr, ChebyshevPolynomial<X> cm);

    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;

  public:
    /** @brief Evaluate at the supplied argument. */
    template<class Y>
    static auto _evaluate(const ChebyshevPolynomial<X>& f, const Y& v) -> ArithmeticType<X, Y>;
};

/** @brief Operations supported by algebra objects. */
template<class X>
struct AlgebraOperations<UnivariateChebyshevPolynomial<X>, X> : public UnivariateChebyshevPolynomial<X> {};

template<class X>
/** @brief Evaluate at the supplied argument. */
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

/** @brief Multivariate polynomial in the tensor-product Chebyshev basis. */
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

    /** @brief Construct a \c MultivariateChebyshevPolynomial from as, pr. */
    MultivariateChebyshevPolynomial(SizeType as, PR pr) : _terms(as, pr) {}

    /** @brief Construct a constant with the supplied value. */
    static ChebyshevPolynomial<X> constant(SizeType as, Number<P> y, PR pr);
    /** @brief Construct a constant with the supplied value. */
    static ChebyshevPolynomial<X> constant(SizeType as, const X& c);
    /** @brief Construct or return the selected coordinate. */
    static ChebyshevPolynomial<X> coordinate(SizeType as, SizeType i, PR pr);
    /** @brief Return the selected basis element. */
    static ChebyshevPolynomial<X> basis(SizeType as, SizeType i, SizeType k, PR pr);

    /** @brief Construct a constant object with the supplied value. */
    ChebyshevPolynomial<X> create_constant(Number<P> y);

    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return this->_terms.argument_size(); }
    /** @brief Return the working precision. */
    PR precision() const { return this->_terms.zero_coefficient().precision(); }
    /** @brief Return a zero with the coefficient characteristics. */
    const X& zero_coefficient() const { return _terms.zero_coefficient(); }

    /** @brief Return the supremum norm. */
    MagType<X> sup_norm() const;
    /** @brief Return the supremum norm. */
    friend MagType<X> sup_norm(const ChebyshevPolynomial<X>& cm) { return cm.sup_norm(); }

    /** @brief Return the product. */
    friend ChebyshevPolynomial<X> operator*(const Number<P>& s1, ChebyshevPolynomial<X> cm2) { return X(s1, cm2.precision()) * cm2; }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ChebyshevPolynomial<X>& cm) { return cm._write(os); }

    /** @brief Apply the callable to the supplied arguments. */
    X operator()(const Vector<X>& x) const;
    /** @brief Evaluate at the supplied argument. */
    friend X evaluate(const ChebyshevPolynomial<X>& f, const Vector<X>& x) { return f(x); }

  private:

  public:
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Pos, ChebyshevPolynomial<X> cm);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Neg, ChebyshevPolynomial<X> cm);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Sqr, ChebyshevPolynomial<X> cm);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Add, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Sub, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Mul, const ChebyshevPolynomial<X>& cm1, const ChebyshevPolynomial<X>& cm2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Add, ChebyshevPolynomial<X> cm1, const Scalar<X>& s2);
    /** @brief Apply the supplied operation to the operands. */
    static ChebyshevPolynomial<X> apply(Mul, ChebyshevPolynomial<X> cm1, const Scalar<X>& s2);
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
};

/** @brief Operations supported by algebra objects. */
template<class X>
struct AlgebraOperations<MultivariateChebyshevPolynomial<X>, X> : public MultivariateChebyshevPolynomial<X> {};

} // namespace Ariadne

#endif // ARIADNE_CHEBYSHEV_POLYNOMIAL_HPP
