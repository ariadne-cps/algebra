/***************************************************************************
 *            algebra/vector-crtp.hpp
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

/** @file algebra/vector-crtp.hpp
 *  @brief CRTP helpers for vector implementations.
 */

#ifndef ARIADNE_VECTOR_HPP
#define ARIADNE_VECTOR_HPP

#include "utility/container.hpp"
#include "utility/module.hpp"

namespace Ariadne {

/************ Vector *********************************************************/

template<class X>
using NegationType = decltype(-declval<X>());
template<class X1, class X2>
using SumType = decltype(declval<X1>() + declval<X2>());
template<class X1, class X2>
using DifferenceType = decltype(declval<X1>() - declval<X2>());
template<class X1, class X2>
using ProductType = decltype(declval<X1>() * declval<X2>());
template<class X1, class X2>
using QuotientType = decltype(declval<X1>() / declval<X2>());

/** @brief Lazy expression representing vector operations. */
template<class V>
struct VectorExpression : Object<V> {};
/** @brief Common representation base for vector objects. */
template<class V>
struct VectorObject : VectorExpression<V> {};

template<class X>
struct Vector : VectorObject<Vector<X>> {
    X _zero;
    Array<X> _ary;

  public:
    typedef X ScalarType;
    /** @brief Construct a \c Vector with the specified characteristics. */
    Vector() : _zero(), _ary() {}
    /** @brief Construct a \c Vector from n. */
    Vector(SizeType n) : _zero(), _ary(n) {}
    /** @brief Construct a \c Vector from n, x. */
    Vector(SizeType n, const X& x) : _zero(x * 0), _ary(n, x) {}
    /** @brief Construct a \c Vector from ary. */
    Vector(Array<X> ary) : _zero(*ary.begin() * 0), _ary(std::move(ary)) {}
    /** @brief Construct a \c Vector from lst. */
    Vector(InitializerList<X> lst) : _zero(X(0)), _ary(lst.begin(), lst.end()) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class V>
    Vector(const VectorExpression<V>& ve) : _zero(ve.upcast().zero_element()), _ary(ve.upcast().size(), _zero) {
        const V& v = ve.upcast();
        for (SizeType i = 0; i != v.size(); ++i) {
            this->_ary[i] = v[i];
        }
    }
    /** @brief Assign the supplied value and return this object. */
    template<class V>
    Vector<X>& operator=(const VectorExpression<V>& ve) {
        const V& v = ve.upcast();
        this->_ary.resize(v.size());
        for (SizeType i = 0; i != v.size(); ++i) {
            this->_ary[i] = v[i];
        }
        return *this;
    }
    /** @brief Return the number of elements. */
    SizeType size() const { return _ary.size(); }
    /** @brief Return the selected element or component. */
    const X& operator[](SizeType i) const { return _ary.at(i); }
    /** @brief Return the selected element or component. */
    X& operator[](SizeType i) { return _ary.at(i); }
    /** @brief Return the element at the specified index. */
    const X& at(SizeType i) const { return _ary.at(i); }
    /** @brief Return the element at the specified index. */
    X& at(SizeType i) { return _ary.at(i); }
    /** @brief Return a zero with the element characteristics. */
    X zero_element() const { return _zero; }
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
};
template<class X>
OutputStream& Vector<X>::_write(OutputStream& os) const {
    if (size() == 0) {
        os << "{";
    }
    for (SizeType i = 0; i != size(); ++i) {
        os << (i == 0u ? "{" : ",") << this->_ary[i];
    }
    return os << "}";
}
/** @brief Write the value to the output stream. */
template<class X>
OutputStream& operator<<(OutputStream& os, const Vector<X>& v) {
    return v._write(os);
}

/** @brief Write the value to the output stream. */
template<class V>
OutputStream& operator<<(OutputStream& os, const VectorExpression<V>& ve) {
    const V& v = ve.upcast();
    typedef decltype(v[0]) X;
    return os << Vector<X>(v);
}

/** @brief Test equality of the operands. */
template<class V1, class V2>
auto operator==(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2)
    -> decltype(ve1.upcast()[0] == ve2.upcast()[0]) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    decltype(v1[0] == v2[0]) r = true;
    for (SizeType i = 0; i != v1.size(); ++i) {
        r = r && (v1[i] == v2[i]);
    }
    return r;
}

/** @brief Test inequality of the operands. */
template<class V1, class V2>
auto operator!=(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2)
    -> decltype(!(ve1 == ve2)) {
    return !(ve1 == ve2);
}

/** @brief Return the sum or unary positive value. */
template<class V>
inline const V& operator+(const VectorExpression<V>& v) {
    return v.upcast();
}

/** @brief Lazy expression for the additive inverse of a vector. */
template<class V>
struct VectorNegation : VectorExpression<VectorNegation<V>> {
    typedef NegationType<typename V::ScalarType> ScalarType;
    const V& _v;
    /** @brief Construct a \c VectorNegation from v. */
    VectorNegation(const V& v) : _v(v) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _v.size(); }
    /** @brief Return a zero with the element characteristics. */
    ScalarType zero_element() const { return -_v.zero_element(); }
    /** @brief Return the selected element or component. */
    ScalarType operator[](SizeType i) const { return -_v[i]; }
};
/** @brief Return the difference or additive inverse. */
template<class V>
inline VectorNegation<V> operator-(const VectorExpression<V>& v) {
    return VectorNegation<V>(v.upcast());
}

/** @brief Lazy expression for the sum of vectors. */
template<class V1, class V2>
struct VectorSum : VectorExpression<VectorSum<V1, V2>> {
    typedef SumType<typename V1::ScalarType, typename V2::ScalarType> ScalarType;
    const V1& _v1;
    const V2& _v2;
    /** @brief Construct a \c VectorSum from v1, v2. */
    VectorSum(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _v1.size(); }
    /** @brief Return a zero with the element characteristics. */
    ScalarType zero_element() const { return _v1.zero_element() + _v2.zero_element(); }
    /** @brief Return the selected element or component. */
    ScalarType operator[](SizeType i) const { return _v1[i] + _v2[i]; }
};
/** @brief Return the sum or unary positive value. */
template<class V1, class V2>
inline VectorSum<V1, V2> operator+(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    ARIADNE_PRECONDITION(v1.size() == v2.size());
    return VectorSum<V1, V2>(v1, v2);
}

/** @brief Add the operand in place. */
template<class X, class XX>
inline Vector<X>& operator+=(const Vector<X>& v1, const Vector<XX>& v2) {
    ARIADNE_PRECONDITION(v1.size() == v2.size());
    for (SizeType i = 0; i != v1.size(); ++i) {
        v1[i] += v2[i];
    }
    return v1;
}

/** @brief Lazy expression for the difference of vectors. */
template<class V1, class V2>
struct VectorDifference : VectorExpression<VectorDifference<V1, V2>> {
    typedef DifferenceType<typename V1::ScalarType, typename V2::ScalarType> ScalarType;
    const V1& _v1;
    const V2& _v2;
    /** @brief Construct a \c VectorDifference from v1, v2. */
    VectorDifference(const V1& v1, const V2& v2) : _v1(v1), _v2(v2) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _v1.size(); }
    /** @brief Return a zero with the element characteristics. */
    ScalarType zero_element() const { return _v1.zero_element() - _v2.zero_element(); }
    /** @brief Return the selected element or component. */
    ScalarType operator[](SizeType i) const { return _v1[i] - _v2[i]; }
};
/** @brief Return the difference or additive inverse. */
template<class V1, class V2>
inline VectorDifference<V1, V2> operator-(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    ARIADNE_PRECONDITION(v1.size() == v2.size());
    return VectorDifference<V1, V2>(v1, v2);
}

/** @brief Subtract the operand in place. */
template<class X, class XX>
inline Vector<X>& operator-=(const Vector<X>& v1, const Vector<XX>& v2) {
    ARIADNE_PRECONDITION(v1.size() == v2.size());
    for (SizeType i = 0; i != v1.size(); ++i) {
        v1[i] -= v2[i];
    }
    return v1;
}

/** @brief Lazy expression for a vector multiplied by a scalar. */
template<class V1, class X2>
struct VectorScalarProduct : VectorExpression<VectorScalarProduct<V1, X2>> {
    typedef ProductType<typename V1::ScalarType, X2> ScalarType;
    const V1& _v1;
    const X2& _x2;
    /** @brief Construct a \c VectorScalarProduct from v1, x2. */
    VectorScalarProduct(const V1& v1, const X2& x2) : _v1(v1), _x2(x2) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _v1.size(); }
    /** @brief Return a zero with the element characteristics. */
    ScalarType zero_element() const { return _v1.zero_element() * _x2; }
    /** @brief Return the selected element or component. */
    ScalarType operator[](SizeType i) const { return _v1[i] * _x2; }
};
/** @brief Return the product. */
template<class X1, class V2>
inline VectorScalarProduct<V2, X1> operator*(const ScalarObject<X1>& x1, const VectorExpression<V2>& v2) {
    return VectorScalarProduct<V2, X1>(v2.upcast(), x1.upcast());
}
/** @brief Return the product. */
template<class V1, class X2>
inline VectorScalarProduct<V1, X2> operator*(const VectorExpression<V1>& v1, const ScalarObject<X2>& x2) {
    return VectorScalarProduct<V1, X2>(v1.upcast(), x2.upcast());
}

/** @brief Multiply by the operand in place. */
template<class X>
inline Vector<X>& operator*=(const Vector<X>& v, const typename Vector<X>::ScalarType& s) {
    for (SizeType i = 0; i != v.size(); ++i) {
        v[i] *= s;
    }
    return v;
}

/** @brief Lazy expression for a vector divided by a scalar. */
template<class V1, class X2>
struct VectorScalarQuotient : VectorExpression<VectorScalarQuotient<V1, X2>> {
    typedef QuotientType<typename V1::ScalarType, X2> ScalarType;
    const V1& _v1;
    const X2& _x2;
    /** @brief Construct a \c VectorScalarQuotient from v1, x2. */
    VectorScalarQuotient(const V1& v1, const X2& x2) : _v1(v1), _x2(x2) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _v1.size(); }
    /** @brief Return a zero with the element characteristics. */
    ScalarType zero_element() const { return _v1.zero_element() / _x2; }
    /** @brief Return the selected element or component. */
    ScalarType operator[](SizeType i) const { return _v1[i] / _x2; }
};
/** @brief Return the quotient. */
template<class V1, class X2>
inline VectorScalarQuotient<V1, X2> operator/(const VectorExpression<V1>& v1, const ScalarObject<X2>& x2) {
    return VectorScalarQuotient<V1, X2>(v1.upcast(), x2.upcast());
}

template<class V>
using ScalarType = typename V::ScalarType;

/** @brief Return the norm. */
template<class X>
inline auto norm(const Vector<X>& v) -> decltype(abs(declval<X>())) {
    decltype(abs(declval<X>())) r = 0u;
    for (SizeType i = 0; i != v.size(); ++i) {
        r = max(r, abs(v[i]));
    }
    return r;
}

/** @brief Return the supremum norm. */
template<class X>
inline auto sup_norm(const Vector<X>& v) -> decltype(mag(declval<X>())) {
    decltype(mag(declval<X>())) r = 0u;
    for (SizeType i = 0; i != v.size(); ++i) {
        r = max(r, mag(v[i]));
    }
    return r;
}

/** @brief Return the scalar product of the operands. */
template<class X1, class X2>
inline auto dot(const Vector<X1>& v1, const Vector<X2>& v2) -> decltype(v1[0] * v2[0]) {
    ARIADNE_PRECONDITION(v1.size() == v2.size());
    decltype(declval<X1>() * declval<X2>() + declval<X1>() * declval<X2>()) r = 0u;
    for (SizeType i = 0; i != v1.size(); ++i) {
        r += v1[i] * v2[i];
    }
    return r;
}

/** @brief Return the scalar product of the operands. */
template<class V1, class V2>
inline auto dot(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2)
    -> decltype(ve1.upcast()[0] * ve2.upcast()[0]) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    ARIADNE_PRECONDITION(v1.size() == v2.size());
    decltype(v1[0] * v2[0]) r = 0u;
    for (SizeType i = 0; i != v1.size(); ++i) {
        r += v1[i] * v2[i];
    }
    return r;
}

/*
template<class V1, class V2> inline auto join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2)
    -> Vector<decltype(max(ve1.upcast()[0],ve2.upcast()[0]))>
{
    const V1& v1=ve1.upcast(); const V2& v2=ve2.upcast();
    typedef ScalarType<V1> X1; typedef ScalarType<V2> X2;
    typedef decltype(max(declval<X1>(),declval<X2>())) X0;
    Vector<X0> r(v1.size()+v2.size());
    for(SizeType i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(SizeType i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    return r;
}
*/

/** @brief Join the operands in their natural order. */
template<class V1, class V2>
    requires Same<ScalarType<V1>, ScalarType<V2>>
Vector<ScalarType<V1>>
    join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    if (v1.size() == 0) {
        return v2;
    }
    if (v2.size() == 0) {
        return v1;
    }
    SizeType n1 = v1.size();
    SizeType n2 = v2.size();
    Vector<ScalarType<V1>> r(n1 + n2, v1[0]);
    for (SizeType i = 0; i != v1.size(); ++i) {
        r[i] = v1[i];
    }
    for (SizeType i = 0; i != v2.size(); ++i) {
        r[v1.size() + i] = v2[i];
    }
    return r;
}

/** @brief Join the operands in their natural order. */
template<class V1, class V2, class V3>
    requires AllSame<ScalarType<V1>, ScalarType<V2>, ScalarType<V3>>
Vector<ScalarType<V1>> join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2, const VectorExpression<V3>& ve3) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    const V3& v3 = ve3.upcast();
    //    typedef ScalarType<V1> X1; typedef ScalarType<V2> X2; typedef ScalarType<V3> X3;
    typedef ScalarType<V1> X0;
    Vector<X0> r(v1.size() + v2.size() + v3.size());
    for (SizeType i = 0; i != v1.size(); ++i) {
        r[i] = v1[i];
    }
    for (SizeType i = 0; i != v2.size(); ++i) {
        r[v1.size() + i] = v2[i];
    }
    for (SizeType i = 0; i != v3.size(); ++i) {
        r[v1.size() + v2.size() + i] = v3[i];
    }
    return r;
}

/** @brief Join the operands in their natural order. */
template<class V1, class X2>
    requires Same<ScalarType<V1>, X2>
Vector<ScalarType<V1>> join(const VectorExpression<V1>& ve1, const ScalarObject<X2>& xe2) {
    const V1& v1 = ve1.upcast();
    const X2& x2 = xe2.upcast();
    // typedef ScalarType<V1> X1;
    typedef X2 X0;
    Vector<X0> r(v1.size() + 1u);
    for (SizeType i = 0; i != v1.size(); ++i) {
        r[i] = v1[i];
    }
    r[v1.size()] = x2;
    return r;
}

/** @brief Join the operands in their natural order. */
template<class V1, class V2, class X3>
    requires AllSame<ScalarType<V1>, ScalarType<V2>, X3>
Vector<ScalarType<V1>> join(const VectorExpression<V1>& ve1, const VectorExpression<V2>& ve2, const ScalarObject<X3>& xe3) {
    const V1& v1 = ve1.upcast();
    const V2& v2 = ve2.upcast();
    const X3& x3 = xe3.upcast();
    // typedef ScalarType<V1> X1; typedef ScalarType<V2> X2;
    typedef X3 X0;
    Vector<X0> r(v1.size() + v2.size() + 1u);
    for (SizeType i = 0; i != v1.size(); ++i) {
        r[i] = v1[i];
    }
    for (SizeType i = 0; i != v2.size(); ++i) {
        r[v1.size() + i] = v2[i];
    }
    r[v1.size() + v2.size()] = x3;
    return r;
}

template<class X>
struct Covector {
    X _zero;
    Array<X> _ary;

  public:
    /** @brief Construct a \c Covector from n. */
    Covector(SizeType n) : _zero(), _ary(n) {}
    /** @brief Construct a \c Covector from n, x. */
    Covector(SizeType n, const X& x) : _zero(x * 0), _ary(n, x) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _ary.size(); }
    /** @brief Return the selected element or component. */
    const X& operator[](SizeType i) const { return _ary.at(i); }
    /** @brief Return the selected element or component. */
    X& operator[](SizeType i) { return _ary.at(i); }
    /** @brief Return a zero with the element characteristics. */
    X zero_element() const { return _zero; }
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
};
template<class X>
OutputStream& Covector<X>::_write(OutputStream& os) const {
    if (size() == 0) {
        os << "{";
    }
    for (SizeType i = 0; i != size(); ++i) {
        os << (i == 0u ? "{" : ",") << this->_ary[i];
    }
    return os << "}";
}
/** @brief Write the value to the output stream. */
template<class X>
OutputStream& operator<<(OutputStream& os, const Covector<X>& v) {
    return v._write(os);
}

} // namespace Ariadne

#endif
