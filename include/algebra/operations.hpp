/***************************************************************************
 *            algebra/operations.hpp
 *
 *  Copyright  2010-20  Pieter Collins
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

/** @file algebra/operations.hpp
 *  @brief Generic operations shared across algebra types.
 */

#ifndef ARIADNE_ALGEBRA_OPERATIONS_HPP
#define ARIADNE_ALGEBRA_OPERATIONS_HPP

#include "algebra/series.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

template<class A, class X = typename A::NumericType>
struct AlgebraOperations;

/** @brief Operations supported by graded algebra objects. */
template<class A>
class GradedAlgebraOperations {
    typedef typename A::NumericType X;

  private:
    /** @brief Compose the outer object with the inner object. */
    static A _compose(const Series<X>& f, const A& a) { return A::_compose(f, a); }

  public:
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Hlf, const A& a) {
        return a / 2;
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Sqr, const A& a) {
        return a * a;
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Div, const A& a1, const A& a2) {
        return mul(a1, rec(a2));
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Pow, const A& a, Nat m) {
        A s = a;
        A r = s.create_constant(1);
        while (m) {
            if (m % 2) {
                r *= s;
            }
            s = sqr(s);
            m /= 2;
        }
        return r;
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Pow, const A& a, Int n) {
        if (n >= 0) {
            return pow(a, Nat(n));
        } else {
            return rec(pow(a, Nat(-n)));
        }
    }

    /** @brief Apply the supplied operation to the operands. */
    template<class OP>
    static A apply(OP op, const A& a) { return _compose(Series<X>(op, a.value()), a); }
};

/** @brief Operations supported by normed algebra objects. */
template<class A>
class NormedAlgebraOperations {
    typedef typename A::NumericType X;

  public:
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Hlf, const A& a) {
        return a / 2;
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Sqr, const A& a) {
        return a * a;
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Pow, const A& a, Nat m) {
        A s = a;
        A r = s.create_constant(1);
        while (m) {
            if (m % 2) {
                r *= s;
            }
            s = sqr(s);
            m /= 2;
        }
        return r;
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Pow, const A& a, Int n) {
        if (n >= 0) {
            return pow(a, Nat(n));
        } else {
            return rec(pow(a, Nat(-n)));
        }
    }
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Div, const A& a1, const A& a2) {
        return mul(a1, rec(a2));
    }

    /** @brief Apply the supplied operation to the operands. */
    static A apply(Rec, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Sqrt, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Exp, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Log, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Sin, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Cos, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Tan, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Tanh, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Asin, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Acos, const A& a);
    /** @brief Apply the supplied operation to the operands. */
    static A apply(Atan, const A& a);
};

/** @brief Declare vector operators for the supported operand types. */
template<class V, class X>
struct DeclareVectorOperators {
    /** @brief Return the sum or unary positive value. */
    friend V operator+(const V& v);
    /** @brief Return the difference or additive inverse. */
    friend V operator-(const V& v);
    /** @brief Return the sum or unary positive value. */
    friend V operator+(const V& v1, const V& v2);
    /** @brief Return the difference or additive inverse. */
    friend V operator-(const V& v1, const V& v2);
    /** @brief Return the product. */
    friend V operator*(const V& v1, const X& x2);
    /** @brief Return the quotient. */
    friend V operator/(const V& v1, const X& x2);
    /** @brief Return the product. */
    friend V operator*(const X& x1, const V& v2);
    /** @brief Add the operand in place. */
    friend V& operator+=(V& v1, const V& v2);
    /** @brief Subtract the operand in place. */
    friend V& operator-=(V& v1, const V& v2);
    /** @brief Multiply by the operand in place. */
    friend V& operator*=(V& v1, const X& x2);
    /** @brief Divide by the operand in place. */
    friend V& operator/=(V& v1, const X& x2);
};

/** @brief Declare vector algebra operators for the supported operand types. */
template<class VA, class A, class VX, class X>
struct DeclareVectorAlgebraOperators
    : DeclareVectorOperators<VA, A> {
    /** @brief Return the sum or unary positive value. */
    friend VA operator+(VA const& va1, VX const& vx2);
    /** @brief Return the sum or unary positive value. */
    friend VA operator+(VX const& vx1, VA const& va2);
    /** @brief Return the difference or additive inverse. */
    friend VA operator-(VA const& va1, VX const& vx2);
    /** @brief Return the difference or additive inverse. */
    friend VA operator-(VX const& vx1, VA const& va2);
    /** @brief Return the product. */
    friend VA operator*(VA const& va1, const X& x2);
    /** @brief Return the quotient. */
    friend VA operator/(VA const& va1, const X& x2);
    /** @brief Return the product. */
    friend VA operator*(const X& x1, VA const& va2);
    /** @brief Add the operand in place. */
    friend VA& operator+=(VA& va1, VX const& vx2);
    /** @brief Subtract the operand in place. */
    friend VA& operator-=(VA& va1, VX const& vx2);
    /** @brief Multiply by the operand in place. */
    friend VA& operator*=(VA& va1, const X& x2);
    /** @brief Divide by the operand in place. */
    friend VA& operator/=(VA& va1, const X& x2);
};

/** @brief Declare algebra operators for the supported operand types. */
template<class A, class X>
struct DeclareAlgebraOperators
    : DeclareVectorOperators<A, X> {
    /** @brief Return the product. */
    friend A operator*(const A& a1, const A& a2);
    /** @brief Multiply by the operand in place. */
    friend A& operator*=(A& a1, const A& a2);
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const A& a1, const X& x2);
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const A& a1, const X& x2);
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const X& x1, const A& a2);
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const X& x1, const A& a2);
    /** @brief Add the operand in place. */
    friend A& operator+=(A& a1, const X& x2);
    /** @brief Subtract the operand in place. */
    friend A& operator-=(A& a1, const X& x2);
};

/** @brief Declare algebra operations for the supported operand types. */
template<class A, class X>
struct DeclareAlgebraOperations
    : DeclareAlgebraOperators<A, X> {
    /** @brief Return the additive zero with the operand characteristics. */
    friend A nul(const A& a);
    /** @brief Return the unary positive value. */
    friend A pos(const A& a);
    /** @brief Return the additive inverse. */
    friend A neg(const A& a);
    /** @brief Divide the value by two. */
    friend A hlf(const A& a);
    /** @brief Return the square. */
    friend A sqr(const A& a);
    /** @brief Return the sum of the operands. */
    friend A add(const A& a1, const A& a2);
    /** @brief Return the difference of the operands. */
    friend A sub(const A& a1, const A& a2);
    /** @brief Return the product of the operands. */
    friend A mul(const A& a1, const A& a2);

    /** @brief Return the sum of the operands. */
    friend A add(const A& a1, const X& x2);
    /** @brief Return the difference of the operands. */
    friend A sub(const A& a1, const X& x2);
    /** @brief Return the product of the operands. */
    friend A mul(const A& a1, const X& x2);
    /** @brief Return the quotient of the operands. */
    friend A div(const A& a1, const X& x2);
    /** @brief Return the sum of the operands. */
    friend A add(const X& x1, const A& a2);
    /** @brief Return the difference of the operands. */
    friend A sub(const X& x1, const A& a2);
    /** @brief Return the product of the operands. */
    friend A mul(const X& x1, const A& a2);
    /** @brief Raise the value to the specified power. */
    friend A pow(const A& a, Nat m);
    /** @brief Raise the value to the specified power. */
    friend A pow(const A& a, Int n);
};

/** @brief Declare division algebra operations for the supported operand types. */
template<class A, class X>
struct DeclareDivisionAlgebraOperations : DeclareAlgebraOperations<A, X> {
    /** @brief Return the quotient. */
    friend A operator/(const A& a1, const A& a2);
    /** @brief Divide by the operand in place. */
    friend A& operator/=(A& a1, const A& a2);
    /** @brief Return the quotient. */
    friend A operator/(const X& x1, const A& a2);

    /** @brief Return the multiplicative inverse. */
    friend A rec(const A& a);
    /** @brief Return the quotient of the operands. */
    friend A div(const A& a1, const A& a2);
    /** @brief Return the quotient of the operands. */
    friend A div(const X& x1, const A& a2);
    /** @brief Raise the value to the specified power. */
    friend A pow(const A& a, Int n);
};

/** @brief Declare transcendental algebra operations for the supported operand types. */
template<class A, class X>
struct DeclareTranscendentalAlgebraOperations : DeclareDivisionAlgebraOperations<A, X> {
    /** @brief Return the square root. */
    friend A sqrt(const A& a);
    /** @brief Return the exponential. */
    friend A exp(const A& a);
    /** @brief Return the natural logarithm. */
    friend A log(const A& a);
    /** @brief Return the sine. */
    friend A sin(const A& a);
    /** @brief Return the cosine. */
    friend A cos(const A& a);
    /** @brief Return the tangent. */
    friend A tan(const A& a);
    /** @brief Return the hyperbolic tangent. */
    friend A tanh(const A& a);
    /** @brief Return the inverse sine. */
    friend A asin(const A& a);
    /** @brief Return the inverse cosine. */
    friend A acos(const A& a);
    /** @brief Return the inverse tangent. */
    friend A atan(const A& a);
};

/** @brief Declare lattice algebra operations for the supported operand types. */
template<class A, class X>
struct DeclareLatticeAlgebraOperations {
    /** @brief Return the maximum of the operands. */
    friend A max(const A& a1, const A& a2);
    /** @brief Return the minimum of the operands. */
    friend A min(const A& a1, const A& a2);
    /** @brief Return the absolute value. */
    friend A abs(const A& a);
};

/** @brief Declare elementary algebra operations for the supported operand types. */
template<class A, class X>
struct DeclareElementaryAlgebraOperations
    : DeclareTranscendentalAlgebraOperations<A, X>,
      DeclareLatticeAlgebraOperations<A, X> {
};

/** @brief Dispatch mixed algebra number operations for the supported operand types. */
template<class A, class Y>
struct DispatchMixedAlgebraNumberOperations;

/** @brief Dispatch concrete generic algebra number operators for the supported operand types. */
template<class A, class X, class Y>
struct DispatchConcreteGenericAlgebraNumberOperators {
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    friend X create(const Y& y, const A& a) { return X(y, a.precision()); }
    /** @brief Add the operand in place. */
    friend A& operator+=(A& a1, const Y& y2) {
        X x2 = create(y2, a1);
        return a1 += x2;
    }
    /** @brief Subtract the operand in place. */
    friend A& operator-=(A& a1, const Y& y2) {
        X x2 = create(y2, a1);
        return a1 -= x2;
    }
    /** @brief Multiply by the operand in place. */
    friend A& operator*=(A& a1, const Y& y2) {
        X x2 = create(y2, a1);
        return a1 *= x2;
    }
    /** @brief Divide by the operand in place. */
    friend A& operator/=(A& a1, const Y& y2) {
        X x2 = create(y2, a1);
        return a1 /= x2;
    }
    /** @brief Return the sum or unary positive value. */
    friend A operator+(A a1, const Y& y2) {
        a1 += y2;
        return a1;
    }
    /** @brief Return the difference or additive inverse. */
    friend A operator-(A a1, const Y& y2) {
        a1 -= y2;
        return a1;
    }
    /** @brief Return the product. */
    friend A operator*(A a1, const Y& y2) {
        a1 *= y2;
        return a1;
    }
    /** @brief Return the quotient. */
    friend A operator/(A a1, const Y& y2) {
        a1 /= y2;
        return a1;
    }
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const Y& y1, A a2) { return a2 + y1; }
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const Y& y1, A a2) { return (-a2) + y1; }
    /** @brief Return the product. */
    friend A operator*(const Y& y1, A a2) { return a2 * y1; }
    /** @brief Return the quotient. */
    friend A operator/(const Y& y1, A a2) { return rec(a2) * y1; }
};

/** @brief Dispatch concrete generic algebra number operations for the supported operand types. */
template<class A, class X, class Y>
struct DispatchConcreteGenericAlgebraNumberOperations : DispatchConcreteGenericAlgebraNumberOperators<A, X, Y> {
    /** @brief Return the sum of the operands. */
    friend A add(A a1, const Y& y2) {
        a1 += y2;
        return a1;
    }
    /** @brief Return the difference of the operands. */
    friend A sub(A a1, const Y& y2) {
        a1 -= y2;
        return a1;
    }
    /** @brief Return the product of the operands. */
    friend A mul(A a1, const Y& y2) {
        a1 *= y2;
        return a1;
    }
    /** @brief Return the quotient of the operands. */
    friend A div(A a1, const Y& y2) {
        a1 /= y2;
        return a1;
    }
    /** @brief Return the sum of the operands. */
    friend A add(const Y& y1, A a2) { return a2 + y1; }
    /** @brief Return the difference of the operands. */
    friend A sub(const Y& y1, A a2) { return (-a2) + y1; }
    /** @brief Return the product of the operands. */
    friend A mul(const Y& y1, A a2) { return a2 * y1; }
    /** @brief Return the quotient of the operands. */
    friend A div(const Y& y1, A a2) { return rec(a2) * y1; }
};

/** @brief Dispatch algebra operators for the supported operand types. */
template<class A, class X>
struct DispatchAlgebraOperators {
    typedef AlgebraOperations<A, X> OperationsType;

  public:
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const A& a) { return pos(a); }
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const A& a) { return neg(a); }
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const A& a1, const A& a2) { return add(a1, a2); }
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const A& a1, const A& a2) { return sub(a1, a2); }
    /** @brief Return the product. */
    friend A operator*(const A& a1, const A& a2) { return mul(a1, a2); }
    /** @brief Add the operand in place. */
    friend A& operator+=(A& a1, const A& a2) { return a1 = add(a1, a2); }
    /** @brief Subtract the operand in place. */
    friend A& operator-=(A& a1, const A& a2) { return a1 = sub(a1, a2); }
    /** @brief Multiply by the operand in place. */
    friend A& operator*=(A& a1, const A& a2) { return a1 = mul(a1, a2); }
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const A& a1, const X& x2) { return add(a1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const A& a1, const X& x2) { return sub(a1, x2); }
    /** @brief Return the product. */
    friend A operator*(const A& a1, const X& x2) { return mul(a1, x2); }
    /** @brief Return the quotient. */
    friend A operator/(const A& a1, const X& x2) { return div(a1, x2); }
    /** @brief Add the operand in place. */
    friend A& operator+=(A& a1, const X& x2) { return a1 = add(a1, x2); }
    /** @brief Subtract the operand in place. */
    friend A& operator-=(A& a1, const X& x2) { return a1 = sub(a1, x2); }
    /** @brief Multiply by the operand in place. */
    friend A& operator*=(A& a1, const X& x2) { return a1 = mul(a1, x2); }
    /** @brief Divide by the operand in place. */
    friend A& operator/=(A& a1, const X& x2) { return a1 = div(a1, x2); }
    /** @brief Return the sum or unary positive value. */
    friend A operator+(const X& x1, const A& a2) { return add(x1, a2); }
    /** @brief Return the difference or additive inverse. */
    friend A operator-(const X& x1, const A& a2) { return sub(x1, a2); }
    /** @brief Return the product. */
    friend A operator*(const X& x1, const A& a2) { return mul(x1, a2); }
};

/** @brief Dispatch algebra operations for the supported operand types. */
template<class A, class X>
struct DispatchAlgebraOperations
    : DispatchAlgebraOperators<A, X> {
    typedef AlgebraOperations<A, X> OperationsType;

  public:
    /** @brief Return the additive zero with the operand characteristics. */
    friend A nul(const A& a) { return OperationsType::apply(Nul(), a); }
    /** @brief Return the unary positive value. */
    friend A pos(const A& a) { return OperationsType::apply(Pos(), a); }
    /** @brief Return the additive inverse. */
    friend A neg(const A& a) { return OperationsType::apply(Neg(), a); }
    /** @brief Divide the value by two. */
    friend A hlf(const A& a) { return OperationsType::apply(Hlf(), a); }
    /** @brief Return the square. */
    friend A sqr(const A& a) { return OperationsType::apply(Sqr(), a); }
    /** @brief Return the sum of the operands. */
    friend A add(const A& a1, const A& a2) { return OperationsType::apply(Add(), a1, a2); }
    /** @brief Return the difference of the operands. */
    friend A sub(const A& a1, const A& a2) { return OperationsType::apply(Sub(), a1, a2); }
    /** @brief Return the product of the operands. */
    friend A mul(const A& a1, const A& a2) { return OperationsType::apply(Mul(), a1, a2); }

    /** @brief Return the sum of the operands. */
    friend A add(const A& a1, const X& x2) { return OperationsType::apply(Add(), a1, x2); }
    /** @brief Return the difference of the operands. */
    friend A sub(const A& a1, const X& x2) { return OperationsType::apply(Add(), a1, neg(x2)); }
    /** @brief Return the product of the operands. */
    friend A mul(const A& a1, const X& x2) { return OperationsType::apply(Mul(), a1, x2); }
    /** @brief Return the quotient of the operands. */
    friend A div(const A& a1, const X& x2) { return OperationsType::apply(Mul(), a1, rec(x2)); }
    /** @brief Return the sum of the operands. */
    friend A add(const X& x1, const A& a2) { return OperationsType::apply(Add(), a2, x1); }
    /** @brief Return the difference of the operands. */
    friend A sub(const X& x1, const A& a2) { return OperationsType::apply(Add(), neg(a2), x1); }
    /** @brief Return the product of the operands. */
    friend A mul(const X& x1, const A& a2) { return OperationsType::apply(Mul(), a2, x1); }
    /** @brief Raise the value to the specified power. */
    friend A pow(const A& a, Nat m) { return GradedAlgebraOperations<A>::apply(Pow(), a, m); }
};

/** @brief Dispatch transcendental algebra operations for the supported operand types. */
template<class A, class X>
struct DispatchTranscendentalAlgebraOperations : DispatchAlgebraOperations<A, X> {
    typedef AlgebraOperations<A, X> OperationsType;

  public:
    // BUG: Cannot declare these due to bug in Clang
    // friend A operator/(A const& a1, X const& x2);
    // friend A div(A const& a1, X const& x2);
    // friend A pow(A const& a, Nat m);

    /** @brief Return the quotient. */
    friend A operator/(const A& a1, const A& a2) { return div(a1, a2); }
    /** @brief Divide by the operand in place. */
    friend A& operator/=(A& a1, const A& a2) { return a1 = div(a1, a2); }
    /** @brief Return the quotient. */
    friend A operator/(const X& x1, const A& a2) { return div(x1, a2); }
    /** @brief Return the quotient of the operands. */
    friend A div(const A& a1, const A& a2) { return OperationsType::apply(Div(), a1, a2); }

    /** @brief Return the multiplicative inverse. */
    friend A rec(const A& a) { return OperationsType::apply(Rec(), a); }
    /** @brief Return the quotient of the operands. */
    friend A div(const X& x1, const A& a2) { return x1 * rec(a2); }
    /** @brief Raise the value to the specified power. */
    friend A pow(const A& a, Int n) { return OperationsType::apply(Pow(), a, n); }

    /** @brief Return the square root. */
    friend A sqrt(const A& a) { return OperationsType::apply(Sqrt(), a); }
    /** @brief Return the exponential. */
    friend A exp(const A& a) { return OperationsType::apply(Exp(), a); }
    /** @brief Return the natural logarithm. */
    friend A log(const A& a) { return OperationsType::apply(Log(), a); }
    /** @brief Return the sine. */
    friend A sin(const A& a) { return OperationsType::apply(Sin(), a); }
    /** @brief Return the cosine. */
    friend A cos(const A& a) { return OperationsType::apply(Cos(), a); }
    /** @brief Return the tangent. */
    friend A tan(const A& a) { return OperationsType::apply(Tan(), a); }
    /** @brief Return the hyperbolic tangent. */
    friend A tanh(const A& a) { return OperationsType::apply(Tanh(), a); }
    /** @brief Return the inverse sine. */
    friend A asin(const A& a) { return OperationsType::apply(Asin(), a); }
    /** @brief Return the inverse cosine. */
    friend A acos(const A& a) { return OperationsType::apply(Acos(), a); }
    /** @brief Return the inverse tangent. */
    friend A atan(const A& a) { return OperationsType::apply(Atan(), a); }
};

/** @brief Dispatch lattice algebra operations for the supported operand types. */
template<class A, class X>
struct DispatchLatticeAlgebraOperations {
    typedef AlgebraOperations<A, X> OperationsType;
    /** @brief Return the maximum of the operands. */
    friend A max(const A& a1, const A& a2) { return OperationsType::apply(Max(), a1, a2); }
    /** @brief Return the minimum of the operands. */
    friend A min(const A& a1, const A& a2) { return OperationsType::apply(Min(), a1, a2); }
    /** @brief Return the absolute value. */
    friend A abs(const A& a) { return OperationsType::apply(Abs(), a); }

    /** @brief Return the maximum of the operands. */
    friend A max(const A& a1, const X& x2) { return OperationsType::apply(Max(), a1, x2); }
    /** @brief Return the minimum of the operands. */
    friend A min(const A& a1, const X& x2) { return OperationsType::apply(Min(), a1, x2); }
    /** @brief Return the maximum of the operands. */
    friend A max(const X& x1, const A& a2) { return OperationsType::apply(Max(), x1, a2); }
    /** @brief Return the minimum of the operands. */
    friend A min(const X& x1, const A& a2) { return OperationsType::apply(Min(), x1, a2); }
};

/** @brief Dispatch elementary algebra operations for the supported operand types. */
template<class A, class X>
struct DispatchElementaryAlgebraOperations
    : DispatchTranscendentalAlgebraOperations<A, X>,
      DispatchLatticeAlgebraOperations<A, X> {
};

/** @brief Factorial value helper for series and differential coefficients. */
struct Factorial {
    Nat _n;
    /** @brief Construct a \c Factorial from n. */
    Factorial(Nat n) : _n(n) {}
    /** @brief Convert to float dp bounds. */
    operator FloatDPBounds() {
        FloatDPBounds r(1, dp);
        for (Nat i = 1; i <= _n; ++i) {
            r *= i;
        }
        return r;
    }
    /** @brief Return the multiplicative inverse. */
    friend FloatDPBounds rec(Factorial x) { return rec(FloatDPBounds(x)); }
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_OPERATIONS_HPP */
