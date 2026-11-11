/***************************************************************************
 *            algebra/covector.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/** @file algebra/covector.hpp
 *  @brief Covectors and their linear algebra operations.
 */

#ifndef ARIADNE_COVECTOR_HPP
#define ARIADNE_COVECTOR_HPP

#include "algebra/vector.hpp"
#include "utility/container.hpp"
#include "utility/metaprogramming.hpp"

namespace Ariadne {

template<class X>
class Vector;
template<class X>
class Covector;

/** @brief Declare covector operations for the supported operand types. */
struct DeclareCovectorOperations {
    /** @brief Return the sum or unary positive value. */
    template<class X>
    friend Covector<X> operator+(const Covector<X>& v);
    /** @brief Return the difference or additive inverse. */
    template<class X>
    friend Covector<NegationType<X>> operator-(const Covector<X>& v);
    /** @brief Return the sum or unary positive value. */
    template<class X1, class X2>
    friend Covector<SumType<X1, X2>> operator+(const Covector<X1>& v1, const Covector<X2>& v2);
    /** @brief Return the difference or additive inverse. */
    template<class X1, class X2>
    friend Covector<DifferenceType<X1, X2>> operator-(const Covector<X1>& v1, const Covector<X2>& v2);
    /** @brief Return the product. */
    template<class X1, class X2>
    friend Covector<ProductType<Scalar<X1>, X2>> operator*(X1 const& x1, const Covector<X2>& v2);
    /** @brief Return the product. */
    template<class X1, class X2>
    friend Covector<ProductType<X1, Scalar<X2>>> operator*(const Covector<X1>& v1, X2 const& x2);
    /** @brief Return the quotient. */
    template<class X1, class X2>
    friend Covector<QuotientType<X1, Scalar<X2>>> operator/(const Covector<X1>& v1, X2 const& x2);
    /** @brief Add the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceSumType<X1, X2>>& operator+=(Covector<X1>& u1, const Covector<X2>& u2);
    /** @brief Subtract the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceDifferenceType<X1, X2>>& operator-=(Covector<X1>& u1, const Covector<X2>& u2);
    /** @brief Multiply by the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceProductType<X1, X2>>& operator*=(Covector<X1>& u1, const X2& s2);
    /** @brief Divide by the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceQuotientType<X1, X2>>& operator/=(Covector<X1>& u1, const X2& s2);
    /** @brief Return the product. */
    template<class X1, class X2>
    friend ArithmeticType<X1, X2> operator*(const Covector<X1>& u1, const Vector<X2>& v2);
    /** @brief Test equality of the operands. */
    template<class X1, class X2>
    friend EqualsType<X1, X2> operator==(const Covector<X1>& u1, const Covector<X2>& u2);
    /** @brief Write the value to the output stream. */
    template<class X>
    friend OutputStream& operator<<(OutputStream& os, const Covector<X>& u);
};

/** @brief Lazy expression representing covector operations. */
template<class U>
struct CovectorExpression : public DeclareCovectorOperations {
    /** @brief Apply the callable to the supplied arguments. */
    const U& operator()() const { return static_cast<const U&>(*this); }
};
/** @brief Container interface for covector elements. */
template<class U>
struct CovectorContainer : public CovectorExpression<U> {};

//! \ingroup LinearAlgebraSubModule
//! \brief A covector (dual to a vector, or a row vector) over some type \a X.
template<class X>
class Covector
    : public CovectorContainer<Covector<X>> {
  private:

  public:
    UniformArray<X> _ary;

  public:
    typedef X ScalarType;
    /** @brief Construct a \c Covector from lst. */
    template<ConvertibleTo<X> XX>
    Covector(const InitializerList<XX>& lst) : _ary(lst) {}
    /** @brief Construct a \c Covector from lst. */
    template<ExplicitlyConvertibleTo<X> XX>
    explicit Covector(const InitializerList<XX>& lst) : _ary(lst._ary) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<ConvertibleTo<X> XX>
    Covector(const Covector<XX>& u) : _ary(u._ary) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<ExplicitlyConvertibleTo<X> XX>
    explicit Covector(const Covector<XX>& u) : _ary(u._ary) {}
    /** @brief Construct a \c Covector from n, prs. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    explicit Covector(SizeType n, PRS... prs) : _ary(n, X(0u, prs...)) {}
    /** @brief Construct a \c Covector from u, prs. */
    template<class XX, class... PRS>
        requires Constructible<X, XX, PRS...>
    explicit Covector(const Covector<XX>& u, PRS... prs) : _ary(u.array(), prs...) {}
    /** @brief Construct a \c Covector with the specified characteristics. */
    explicit Covector() : _ary() {}
    /** @brief Construct a \c Covector from n. */
    explicit Covector(SizeType n) : _ary(n) {}
    /** @brief Construct a \c Covector from n, x. */
    explicit Covector(SizeType n, const X& x) : _ary(n, x) {}
    /** @brief Construct a \c Covector from ary. */
    explicit Covector(Array<X> ary) : _ary(std::move(ary)) {}
    /** @brief Construct the multiplicative unit. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    static Covector<X> unit(SizeType n, SizeType j, PRS... prs) {
        Covector<X> r(n, prs...);
        r[j] = 1;
        return r;
    }
    /** @brief Construct a \c Covector from n, g, prs. */
    template<class G, class... PRS>
        requires InvocableReturning<X, G, SizeType> and Constructible<X, PRS...>
    explicit Covector(SizeType n, const G& g, PRS... prs) : _ary(n, g, prs...) {}
    /** @brief Construct a \c Covector from n, g, prs. */
    template<class G, class... PRS>
        requires InvocableReturning<X, G, SizeType> and Constructible<X, PRS...>
    explicit Covector(SizeType n, const G& g, Tuple<PRS...> prs) : _ary(n, g, prs) {}

    /** @brief Return the number of elements. */
    SizeType size() const { return _ary.size(); }
    /** @brief Change the number of elements. */
    Void resize(SizeType n) { _ary.resize(n); }
    /** @brief Return the selected element or component. */
    const X& operator[](SizeType j) const { return _ary[j]; }
    /** @brief Return the selected element or component. */
    ElementReference<X> operator[](SizeType j) { return _ary[j]; }
    /** @brief Return the element at the specified index. */
    const X& at(SizeType j) const {
        ARIADNE_PRECONDITION_MSG(j < this->size(), *this << "[" << j << "]");
        return _ary[j];
    }
    /** @brief Return the element at the specified index. */
    ElementReference<X>& at(SizeType j) {
        ARIADNE_PRECONDITION_MSG(j < this->size(), *this << "[" << j << "]");
        return _ary[j];
    }
    /** @brief Return a zero with the element characteristics. */
    X zero_element() const {
        ARIADNE_DEBUG_ASSERT(not _ary.empty());
        return create_zero(_ary[0]);
    }
    /** @brief Return the underlying array. */
    const UniformArray<X>& array() const { return _ary; }

    /** @brief Copy-construct the object from the supplied instance. */
    template<class CVE>
        requires Convertible<typename CVE::ScalarType, X>
    Covector(const CovectorExpression<CVE>& cve) : _ary(cve().size(), cve().zero_element()) {
        for (SizeType i = 0; i != this->size(); ++i) {
            this->_ary[i] = cve()[i];
        }
    }
    /** @brief Assign the supplied value and return this object. */
    template<class CVE>
        requires Assignable<typename CVE::ScalarType, X>
    const Covector<X>& operator=(const CovectorExpression<CVE>& cve) {
        this->resize(cve.size());
        for (SizeType i = 0; i != this->size(); ++i) {
            this->_ary[i] = cve()[i];
        }
    }
};

/** @brief Return the transposed matrix. */
template<class X>
inline const Covector<X>& transpose(const Vector<X>& v) {
    return reinterpret_cast<const Covector<X>&>(v);
}

/** @brief Return the transposed matrix. */
template<class X>
inline const Vector<X>& transpose(const Covector<X>& u) {
    return reinterpret_cast<const Vector<X>&>(u);
}

/** @brief View of a contiguous range of covector components. */
template<class U>
struct CovectorRange
    : public CovectorContainer<CovectorRange<U>> {
    typedef typename U::ScalarType ScalarType;
    const U& _u;
    Range _rng;
    /** @brief Construct a \c CovectorRange from u, rng. */
    CovectorRange(const U& u, Range rng) : _u(u), _rng(rng) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _rng.size(); }
    /** @brief Return a zero with the element characteristics. */
    ScalarType zero_element() const { return _u.zero_element(); }
    /** @brief Return the selected element or component. */
    ScalarType operator[](SizeType i) const { return _u[i + _rng.start()]; }
};

/** @brief Project onto the selected coordinates. */
template<class V>
inline CovectorRange<V> project(const CovectorExpression<V>& v, Range rng) {
    return CovectorRange<V>(v(), rng);
}

/** @brief Provide covector operations for the supported operand types. */
class ProvideCovectorOperations {
    /** @brief Return the sum or unary positive value. */
    template<class X>
    friend Covector<X> operator+(const Covector<X>& u) {
        Covector<X> r(u.size(), +u.zero_element());
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = +u[i];
        }
        return r;
    }

    /** @brief Return the difference or additive inverse. */
    template<class X>
    friend Covector<NegationType<X>> operator-(const Covector<X>& u) {
        Covector<NegationType<X>> r(u.size(), -u.zero_element());
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = -u[i];
        }
        return r;
    }

    /** @brief Return the sum or unary positive value. */
    template<class X1, class X2>
    friend Covector<SumType<X1, X2>> operator+(const Covector<X1>& u1, const Covector<X2>& u2) {
        ARIADNE_PRECONDITION(u1.size() == u2.size());
        Covector<SumType<X1, X2>> r(u1.size(), u1.zero_element() + u2.zero_element());
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = u1[i] + u2[i];
        }
        return r;
    }

    /** @brief Return the difference or additive inverse. */
    template<class X1, class X2>
    friend Covector<DifferenceType<X1, X2>> operator-(const Covector<X1>& u1, const Covector<X2>& u2) {
        ARIADNE_PRECONDITION(u1.size() == u2.size());
        Covector<DifferenceType<X1, X2>> r(u1.size(), u1.zero_element() - u2.zero_element());
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = u1[i] - u2[i];
        }
        return r;
    }

    /** @brief Return the product. */
    template<class X1, class X2>
    friend Covector<ProductType<Scalar<X1>, X2>> operator*(X1 const& s1, const Covector<X2>& u2) {
        Covector<ProductType<X1, X2>> r(u2.size(), s1 * u2.zero_element());
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = s1 * u2[i];
        }
        return r;
    }

    /** @brief Return the product. */
    template<class X1, class X2>
    friend Covector<ProductType<X1, Scalar<X2>>> operator*(const Covector<X1>& u1, X2 const& s2) {
        Covector<ProductType<X1, X2>> r(u1.size(), u1.zero_element() * s2);
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = u1[i] * s2;
        }
        return r;
    }

    /** @brief Return the quotient. */
    template<class X1, class X2>
    friend Covector<QuotientType<X1, Scalar<X2>>> operator/(const Covector<X1>& u1, X2 const& s2) {
        Covector<QuotientType<X1, X2>> r(u1.size(), u1.zero_element() * s2);
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = u1[i] / s2;
        }
        return r;
    }

    /** @brief Add the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceSumType<X1, X2>>& operator+=(Covector<X1>& u1, const Covector<X2>& u2) {
        ARIADNE_PRECONDITION(u1.size() == u2.size());
        for (SizeType i = 0; i != u1.size(); ++i) {
            u1[i] += u2[i];
        }
        return u1;
    }

    /** @brief Subtract the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceDifferenceType<X1, X2>>& operator-=(Covector<X1>& u1, const Covector<X2>& u2) {
        ARIADNE_PRECONDITION(u1.size() == u2.size());
        for (SizeType i = 0; i != u1.size(); ++i) {
            u1[i] -= u2[i];
        }
        return u1;
    }

    /** @brief Multiply by the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceProductType<X1, X2>>& operator*=(Covector<X1>& u1, X2 const& s2) {
        for (SizeType i = 0; i != u1.size(); ++i) {
            u1[i] *= s2;
        }
        return u1;
    }

    /** @brief Divide by the operand in place. */
    template<class X1, class X2>
    friend Covector<InplaceQuotientType<X1, X2>>& operator/=(Covector<X1>& u1, X2 const& s2) {
        for (SizeType i = 0; i != u1.size(); ++i) {
            u1[i] /= s2;
        }
        return u1;
    }

    /** @brief Return the product. */
    template<class X1, class X2>
    friend ArithmeticType<X1, X2> operator*(const Covector<X1>& u1, const Vector<X2>& v2) {
        ARIADNE_PRECONDITION(u1.size() == v2.size());
        ArithmeticType<X1, X2> r = u1.zero_element() * v2.zero_element();
        for (SizeType i = 0; i != v2.size(); ++i) {
            r += u1[i] * v2[i];
        }
        return r;
    }

    /** @brief Test equality of the operands. */
    template<class X1, class X2>
    friend EqualsType<X1, X2> operator==(const Covector<X1>& u1, const Covector<X2>& u2) {
        if (u1.size() != u2.size()) {
            return false;
        }
        decltype(u1[0] == u2[0]) r = true;
        for (SizeType i = 0; i != u1.size(); ++i) {
            r = r && (u1[i] == u2[i]);
        }
        return r;
    }

    /** @brief Write the value to the output stream. */
    template<class X>
    friend OutputStream& operator<<(OutputStream& os, const Covector<X>& u) {
        if (u.size() == 0) {
            os << "{";
        }
        for (SizeType i = 0; i != u.size(); ++i) {
            os << (i == 0u ? "{" : ",") << u[i];
        }
        return os << "}";
    }
};

/** @brief Join the operands along the complementary matrix dimension. */
template<class X>
Covector<X> cojoin(const Covector<X>& u1, const Covector<X>& u2) {
    Array<X> ra(u1.size() + u2.size(), Uninitialised());
    X* rp = ra.begin();
    for (const X* u1p = u1.array().begin(); u1p != u1.array().end(); ++rp, ++u1p) {
        new (rp) X(*u1p);
    }
    for (const X* u2p = u2.array().begin(); u2p != u2.array().end(); ++rp, ++u2p) {
        new (rp) X(*u2p);
    }
    assert(rp == ra.end());
    return Covector<X>(std::move(ra));
}

/** @brief Join the operands along the complementary matrix dimension. */
template<class X>
Covector<X> cojoin(const Covector<X>& u1, const typename Covector<X>::ScalarType& s2) {
    Array<X> ra(u1.size() + 1u, Uninitialised());
    X* rp = ra.begin();
    for (const X* u1p = u1.array().begin(); u1p != u1.array().end(); ++rp, ++u1p) {
        new (rp) X(*u1p);
    }
    new (rp) X(s2);
    return Covector<X>(std::move(ra));
}

/** @brief Join the operands along the complementary matrix dimension. */
template<class X>
Covector<X> cojoin(const typename Covector<X>::ScalarType& s1, const Covector<X>& u2) {
    Array<X> ra(1u + u2.size(), Uninitialised());
    X* rp = ra.begin();
    new (rp) X(s1);
    ++rp;
    for (const X* u2p = u2.array().begin(); u2p != u2.array().end(); ++rp, ++u2p) {
        new (rp) X(*u2p);
    }
    return Covector<X>(std::move(ra));
}

} // namespace Ariadne

#endif
