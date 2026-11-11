/***************************************************************************
 *            algebra/algebra_concepts.hpp
 *
 *  Copyright  2023  Pieter Collins
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

/** @file algebra/algebra_concepts.hpp
 *  @brief C++ concepts for algebraic structures and vector spaces.
 */

#ifndef ARIADNE_ALGEBRA_CONCEPTS_HPP
#define ARIADNE_ALGEBRA_CONCEPTS_HPP

#include "numeric/concepts.hpp"
#include "utility/metaprogramming.hpp"

namespace Ariadne {

//! \brief Concept for a vector over a ring or field \a X.
template<class V, class X>
concept AVectorOver = Ring<X> and requires(const V& v, X x, typename V::IndexType i) {
    { declval<V>()[declval<typename V::IndexType>] = declval<X>() };
    { v[i] } -> ConvertibleTo<X>;
    { +v } -> ConvertibleTo<V>;
    { -v } -> ConvertibleTo<V>;
    { v + v } -> ConvertibleTo<V>;
    { v - v } -> ConvertibleTo<V>;
    { x * v } -> ConvertibleTo<V>;
    { v * x } -> ConvertibleTo<V>;
    { v / x } -> ConvertibleTo<V>;
};

//! \brief Concept for a unital algebra over a field \a X.
template<class A, class X>
concept AnAlgebraOver = Field<X> and Assignable<A, X> and requires(const A& a, const X& x, Nat m) {
    { operator+(a) } -> SameAsOrConvertibleTo<A>;
    { operator-(a) } -> SameAsOrConvertibleTo<A>;
    { operator+(a, a) } -> SameAsOrConvertibleTo<A>;
    { operator-(a, a) } -> SameAsOrConvertibleTo<A>;
    { operator*(a, a) } -> SameAsOrConvertibleTo<A>;
    { operator+(x, a) } -> SameAsOrConvertibleTo<A>;
    { operator-(x, a) } -> SameAsOrConvertibleTo<A>;
    { operator*(x, a) } -> SameAsOrConvertibleTo<A>;
    { operator+(a, x) } -> SameAsOrConvertibleTo<A>;
    { operator-(a, x) } -> SameAsOrConvertibleTo<A>;
    { operator*(a, x) } -> SameAsOrConvertibleTo<A>;
    { operator/(a, x) } -> SameAsOrConvertibleTo<A>;
    { operator+=(declval<A&>(), a) } -> SameAs<A&>;
    { operator-=(declval<A&>(), a) } -> SameAs<A&>;
    { operator+=(declval<A&>(), x) } -> SameAs<A&>;
    { operator-=(declval<A&>(), x) } -> SameAs<A&>;
    { operator*=(declval<A&>(), x) } -> SameAs<A&>;
    { operator/=(declval<A&>(), x) } -> SameAs<A&>;

    { nul(a) } -> SameAsOrConvertibleTo<A>;
    { pos(a) } -> SameAsOrConvertibleTo<A>;
    { neg(a) } -> SameAsOrConvertibleTo<A>;
    { sqr(a) } -> SameAsOrConvertibleTo<A>;
    { hlf(a) } -> SameAsOrConvertibleTo<A>;
    { add(a, a) } -> SameAsOrConvertibleTo<A>;
    { sub(a, a) } -> SameAsOrConvertibleTo<A>;
    { mul(a, a) } -> SameAsOrConvertibleTo<A>;
    { add(a, x) } -> SameAsOrConvertibleTo<A>;
    { mul(a, x) } -> SameAsOrConvertibleTo<A>;
    { pow(a, m) } -> SameAsOrConvertibleTo<A>;
};

template<class A, class X>
concept ATranscendentalAlgebraOver = AnAlgebraOver<A, X> and Transcendental<A>;

template<class A, class X>
concept AnElementaryAlgebraOver = AnAlgebraOver<A, X> and Transcendental<A> and requires(A a, X x) {
    { max(a, a) } -> SameAsOrConvertibleTo<A>;
    { min(a, a) } -> SameAsOrConvertibleTo<A>;
    { max(a, x) } -> SameAsOrConvertibleTo<A>;
    { min(a, x) } -> SameAsOrConvertibleTo<A>;
    { max(x, a) } -> SameAsOrConvertibleTo<A>;
    { min(x, a) } -> SameAsOrConvertibleTo<A>;
};

//! \brief Concept for a unital algebra over a field \a X supporting inplace operations.
template<class A, class X>
concept AnInplaceAlgebraOver = Field<X> and requires(A& r, const X& c, const A& a) {
    { r.iadd(c) };
    { r.imul(c) };
    { r.isma(c, a) };
    { r.ifma(a, a) };
};

//! \brief Concept for a unital algebra.
template<class A>
concept AnAlgebra = requires { typename A::NumericType; } and AnAlgebraOver<A, typename A::NumericType>;

/** @brief Reference type specifying the requirements of algebra. */
template<class X>
class AlgebraArchetype {
    using ALG = AlgebraArchetype<X>;

  public:
    typedef X NumericType;
    /** @brief Assign the supplied value and return this object. */
    ALG& operator=(const X&);
    /** @brief Return the sum or unary positive value. */
    friend ALG operator+(ALG const&);
    /** @brief Return the difference or additive inverse. */
    friend ALG operator-(ALG const&);
    /** @brief Return the sum or unary positive value. */
    friend ALG operator+(ALG const&, ALG const&);
    /** @brief Return the difference or additive inverse. */
    friend ALG operator-(ALG const&, ALG const&);
    /** @brief Return the product. */
    friend ALG operator*(ALG const&, ALG const&);
    /** @brief Return the sum or unary positive value. */
    friend ALG operator+(const X&, ALG const&);
    /** @brief Return the difference or additive inverse. */
    friend ALG operator-(const X&, ALG const&);
    /** @brief Return the product. */
    friend ALG operator*(const X&, ALG const&);
    /** @brief Return the sum or unary positive value. */
    friend ALG operator+(ALG const&, const X&);
    /** @brief Return the difference or additive inverse. */
    friend ALG operator-(ALG const&, const X&);
    /** @brief Return the product. */
    friend ALG operator*(ALG const&, const X&);
    /** @brief Return the quotient. */
    friend ALG operator/(ALG const&, const X&);

    /** @brief Add the operand in place. */
    friend ALG& operator+=(ALG&, ALG const&);
    /** @brief Subtract the operand in place. */
    friend ALG& operator-=(ALG&, ALG const&);
    /** @brief Add the operand in place. */
    friend ALG& operator+=(ALG&, const X&);
    /** @brief Subtract the operand in place. */
    friend ALG& operator-=(ALG&, const X&);
    /** @brief Multiply by the operand in place. */
    friend ALG& operator*=(ALG&, const X&);
    /** @brief Divide by the operand in place. */
    friend ALG& operator/=(ALG&, const X&);

    /** @brief Return the additive zero with the operand characteristics. */
    friend ALG nul(ALG const&);
    /** @brief Return the unary positive value. */
    friend ALG pos(ALG const&);
    /** @brief Return the additive inverse. */
    friend ALG neg(ALG const&);
    /** @brief Return the square. */
    friend ALG sqr(ALG const&);
    /** @brief Divide the value by two. */
    friend ALG hlf(ALG const&);
    /** @brief Return the sum of the operands. */
    friend ALG add(ALG const&, ALG const&);
    /** @brief Return the difference of the operands. */
    friend ALG sub(ALG const&, ALG const&);
    /** @brief Return the product of the operands. */
    friend ALG mul(ALG const&, ALG const&);
    /** @brief Return the sum of the operands. */
    friend ALG add(const X&, ALG const&);
    /** @brief Return the difference of the operands. */
    friend ALG sub(const X&, ALG const&);
    /** @brief Return the product of the operands. */
    friend ALG mul(const X&, ALG const&);
    /** @brief Return the sum of the operands. */
    friend ALG add(ALG const&, const X&);
    /** @brief Return the difference of the operands. */
    friend ALG sub(ALG const&, const X&);
    /** @brief Return the product of the operands. */
    friend ALG mul(ALG const&, const X&);
    /** @brief Return the quotient of the operands. */
    friend ALG div(ALG const&, const X&);
    /** @brief Raise the value to the specified power. */
    friend ALG pow(ALG const&, const Nat&);

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const AlgebraArchetype<X>&);
};

/** @brief Reference type specifying the requirements of inplace algebra. */
template<class X>
class InplaceAlgebraArchetype {
    using ALG = InplaceAlgebraArchetype<X>;

  public:
    typedef X NumericType;
    /** @brief Assign the supplied value and return this object. */
    ALG& operator=(const X&);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    ALG create() const;
    /** @brief Construct a constant object with the supplied value. */
    ALG create_constant(const X&) const;
    /** @brief Add the operand in place. */
    ALG& iadd(const X&);
    /** @brief Multiply by the operand in place. */
    ALG& imul(const X&);
    /** @brief Accumulate the scaled value in place. */
    ALG& isma(const X&, ALG const&);
    /** @brief Accumulate the product in place. */
    ALG& ifma(ALG const&, ALG const&);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, ALG const& a);
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_CONCEPTS_HPP */
