/***************************************************************************
 *            algebra/univariate_differential.hpp
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

/** @file algebra/univariate_differential.hpp
 *  @brief Univariate differential representations and operations.
 */

#ifndef ARIADNE_UNIVARIATE_DIFFERENTIAL_HPP
#define ARIADNE_UNIVARIATE_DIFFERENTIAL_HPP

#include "algebra/series.hpp"
#include "numeric/numeric.hpp"
#include "utility/array.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

template<class X>
class Series;
template<class X>
class UnivariateDifferential;

//! \ingroup DifferentiationSubModule
//! \brief Arbitrary-order derivatives with respect to a single argument.
template<class X>
class UnivariateDifferential
    : public DispatchTranscendentalAlgebraOperations<UnivariateDifferential<X>, X>,
      public DispatchLatticeAlgebraOperations<UnivariateDifferential<X>, X>,
      public ProvideConcreteGenericArithmeticOperations<UnivariateDifferential<X>> {
    Array<X> _ary;

  public:
    typedef X NumericType;
    typedef X ValueType;
    typedef typename X::Paradigm Paradigm;
    typedef UnivariateDifferential<X> SelfType;
    typedef typename Array<X>::Iterator Iterator;
    typedef typename Array<X>::ConstIterator ConstIterator;

    /** @brief Disallow this construction to preserve the ownership contract. */
    UnivariateDifferential() = delete;
    /** @brief Construct a \c UnivariateDifferential from d, c. */
    explicit UnivariateDifferential(DegreeType d, const NumericType& c);
    /** @brief Construct a \c UnivariateDifferential from d, prs. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    explicit UnivariateDifferential(DegreeType d, PRS... prs) : UnivariateDifferential(d, X(0u, prs...)) {}
    /** @brief Construct a \c UnivariateDifferential from d, e. */
    UnivariateDifferential(DegreeType d, InitializerList<X> e);
    /** @brief Construct a \c UnivariateDifferential from ary. */
    UnivariateDifferential(Array<X> ary) : _ary(ary) {}
    /** @brief Construct a \c UnivariateDifferential from deg, lst, prs. */
    template<class... PRS>
        requires Constructible<X, ExactDouble, PRS...>
    explicit UnivariateDifferential(DegreeType deg, InitializerList<ExactDouble> lst, PRS... prs);
    /** @brief Construct a \c UnivariateDifferential from d, s. */
    UnivariateDifferential(DegreeType d, const Series<X>& s); // explicit

    /** @brief Construct a \c UnivariateDifferential from op, d, c. */
    template<class OP>
    UnivariateDifferential(OP op, DegreeType d, const X& c);

    /** @brief Construct a constant with the supplied value. */
    static SelfType constant(DegreeType d, const NumericType& c);
    /** @brief Construct the selected independent variable. */
    static SelfType variable(DegreeType d, const NumericType& c);

    /** @brief Construct a constant with the supplied value. */
    template<class Y, class... PRS>
        requires Constructible<X, Y, PRS...>
    static SelfType constant(DegreeType deg, const Y& c, const PRS&... prs) {
        return constant(deg, X(c, prs...));
    }
    /** @brief Construct the selected independent variable. */
    template<class Y, class... PRS>
        requires Constructible<X, Y, PRS...>
    static SelfType variable(DegreeType deg, const Y& v, const PRS&... prs) {
        return variable(deg, X(v, prs...));
    }

    /** @brief Construct a zero object with compatible characteristics. */
    SelfType create_zero() const;
    /** @brief Construct a constant object with the supplied value. */
    SelfType create_constant(const NumericType& c) const;
    /** @brief Construct the selected independent variable with compatible characteristics. */
    SelfType create_variable(const NumericType& c) const;

    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const;
    /** @brief Return the polynomial or expansion degree. */
    DegreeType degree() const;
    /** @brief Return a zero with the coefficient characteristics. */
    X zero_coefficient() const;
    /** @brief Return the underlying array. */
    const Array<X>& array() const;
    /** @brief Return the underlying array. */
    Array<X>& array();
    /** @brief Return the selected element or component. */
    const X& operator[](SizeType k) const;
    /** @brief Return the selected element or component. */
    X& operator[](SizeType k);

    /** @brief Assign the supplied value and return this object. */
    SelfType& operator=(const NumericType& c);
    /** @brief Assign the supplied value and return this object. */
    template<AssignableTo<X> W>
    SelfType& operator=(const W& c) {
        X xc = nul(this->value());
        xc = c;
        return (*this) = xc;
    }

    /** @brief Apply the supplied operation to the operands. */
    SelfType apply(Rec) const;

    /** @brief Return the stored value. */
    const X& value() const;
    /** @brief Return the vector of first partial derivatives. */
    const X& gradient() const;
    /** @brief Return the matrix of second partial derivatives. */
    const X hessian() const;
    /** @brief Return half the matrix of second partial derivatives. */
    const X& half_hessian() const;

    /** @brief Remove all elements. */
    Void clear();
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;

    /** @brief Compose the outer object with the inner object. */
    friend UnivariateDifferential<X> compose(const Series<X>& f, const UnivariateDifferential<X>& dx) {
        return UnivariateDifferential<X>::_compose(f, dx);
    }
    /** @brief Compose the outer object with the inner object. */
    friend UnivariateDifferential<X> compose(const UnivariateDifferential<X>& f, const UnivariateDifferential<X>& dx) {
        return UnivariateDifferential<X>::_compose(f, dx);
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend UnivariateDifferential<X> derivative(const UnivariateDifferential<X>& dx) {
        return UnivariateDifferential<X>::_derivative(dx);
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    friend UnivariateDifferential<X> antiderivative(const UnivariateDifferential<X>& dx) {
        return UnivariateDifferential<X>::_antiderivative(dx);
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    friend UnivariateDifferential<X> antiderivative(const UnivariateDifferential<X>& dx, const X& c) {
        return UnivariateDifferential<X>::_antiderivative(dx, c);
    }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const UnivariateDifferential<X>& dx) {
        return dx._write(os);
    }

  public:
    /** @brief Compose the outer object with the inner object. */
    static Differential<X> _compose(const Series<X>& f, const Differential<X>& dx);

  private:
    /** @brief Compose the outer object with the inner object. */
    static UnivariateDifferential<X> _compose(const Series<X>& f, const UnivariateDifferential<X>& dx);
    /** @brief Compose the outer object with the inner object. */
    static UnivariateDifferential<X> _compose(const UnivariateDifferential<X>& f, const UnivariateDifferential<X>& dx);
    /** @brief Return the derivative with respect to the selected argument. */
    static UnivariateDifferential<X> _derivative(const UnivariateDifferential<X>& dx);
    /** @brief Return an antiderivative with respect to the selected argument. */
    static UnivariateDifferential<X> _antiderivative(const UnivariateDifferential<X>& dx);
    /** @brief Return an antiderivative with respect to the selected argument. */
    static UnivariateDifferential<X> _antiderivative(const UnivariateDifferential<X>& dx, const X& c);
};

template<class X>
template<class OP>
UnivariateDifferential<X>::UnivariateDifferential(OP op, DegreeType d, const X& c)
    : UnivariateDifferential(d, Series<X>(op, c)) {}

template<class X>
inline const X& UnivariateDifferential<X>::value() const { return _ary[0]; }
template<class X>
inline const X& UnivariateDifferential<X>::gradient() const {
    ARIADNE_PRECONDITION(this->degree() >= 1);
    return _ary[1];
}
template<class X>
inline const X& UnivariateDifferential<X>::half_hessian() const {
    ARIADNE_PRECONDITION(this->degree() >= 2);
    return _ary[2];
}
template<class X>
inline const X UnivariateDifferential<X>::hessian() const {
    ARIADNE_PRECONDITION(this->degree() >= 2);
    return _ary[2] * 2;
}

template<class X>
template<class... PRS>
    requires Constructible<X, ExactDouble, PRS...>
UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, InitializerList<ExactDouble> lst, PRS... prs)
    : _ary(d + 1u, X(0, prs...)) {
    auto iter = lst.begin();
    for (SizeType i = 0; iter != lst.end(); ++i, ++iter) {
        _ary[i] = X(*iter, prs...);
    }
}

/** @brief Construct a zero object with compatible characteristics. */
template<class X>
UnivariateDifferential<X> create_zero(const UnivariateDifferential<X>& c) {
    return UnivariateDifferential<X>(c.degree(), create_zero(c[0]));
}

/** @brief Operations supported by algebra objects. */
template<class X>
struct AlgebraOperations<UnivariateDifferential<X>, X> {
    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Pos, UnivariateDifferential<X> x) {
        return x;
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Neg, UnivariateDifferential<X> x) {
        for (DegreeType i = 0; i <= x.degree(); ++i) {
            x[i] = -x[i];
        }
        return x;
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Sqr, const UnivariateDifferential<X>& dx) {
        return apply(Mul(), dx, dx);
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Add, const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
        UnivariateDifferential<X> r(std::min(x1.degree(), x2.degree()), nul(x1[0] + x2[0]));
        for (DegreeType i = 0; i <= r.degree(); ++i) {
            r[i] = x1[i] + x2[i];
        }
        return r;
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Sub, const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
        UnivariateDifferential<X> r(std::min(x1.degree(), x2.degree()), nul(x1[0] - x2[0]));
        for (DegreeType i = 0; i <= r.degree(); ++i) {
            r[i] = x1[i] - x2[i];
        }
        return r;
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Mul, const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
        UnivariateDifferential<X> r(std::min(x1.degree(), x2.degree()), nul(x1[0] * x2[0]));
        for (DegreeType i1 = 0; i1 <= r.degree(); ++i1) {
            for (DegreeType i2 = 0; i2 <= r.degree() - i1; ++i2) {
                r[static_cast<DegreeType>(i1 + i2)] += x1[i1] * x2[i2];
            }
        }
        return r;
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Div, const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
        return x1 * rec(x2);
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Add, UnivariateDifferential<X> x, const X& c) {
        x[0] += c;
        return x;
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Mul, UnivariateDifferential<X> x, const X& c) {
        for (DegreeType i = 0; i <= x.degree(); ++i) {
            x[i] *= c;
        }
        return x;
    }

    /** @brief Apply the supplied operation to the operands. */
    template<class OP>
    static UnivariateDifferential<X> apply(OP op, const UnivariateDifferential<X>& dx) {
        return compose(UnivariateDifferential<X>(op, dx.degree(), dx[0]), dx);
    }

    /** @brief Apply the supplied operation to the operands. */
    static UnivariateDifferential<X> apply(Pow, UnivariateDifferential<X> x, Int n) {
        return generic_pow(x, n);
    }
};

} // namespace Ariadne

#endif // ARIADNE_UNIVARIATE_DIFFERENTIAL_HPP

