/***************************************************************************
 *            algebra/algebra_wrapper.hpp
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

/** @file algebra/algebra_wrapper.hpp
 *  @brief Wrappers adapting concrete types to the algebra interfaces.
 */

#ifndef ARIADNE_ALGEBRA_WRAPPER_HPP
#define ARIADNE_ALGEBRA_WRAPPER_HPP

#include "algebra/algebra_interface.hpp"
#include "algebra/algebra_mixin.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

/** @brief Adapter exposing algebra. */
template<class A, class X = NumericType<A>>
class AlgebraWrapper;

/** @brief Adapter exposing algebra. */
template<class A, class X>
    requires AnAlgebraOver<A, X>
class AlgebraWrapper<A, X>
    : public virtual AlgebraInterface<X>,
      public A {
    /** @brief Return a dynamically allocated moved value. */
    static AlgebraInterface<X>* heap_move(A&& a) { return new AlgebraWrapper<A, X>(std::forward<A>(a)); }

  public:
    typedef X NumericType;

    /** @brief Move-construct the object from the supplied instance. */
    AlgebraWrapper(A&& a) : A(std::move(a)) {}
    /** @brief Construct a \c AlgebraWrapper from a. */
    AlgebraWrapper(const A& a) : A(a) {}
    virtual AlgebraInterface<X>* _copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual AlgebraInterface<X>* _create_zero() const override { return heap_move(nul(*this)); }
    virtual AlgebraInterface<X>* _create_constant(const X& c) const override {
        A a = nul(*this);
        a = c;
        return heap_move(std::move(a));
    }
    virtual AlgebraInterface<X>* _create_copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual Void _iadd(const X& c) override {
        static_cast<A&>(*this) += c;
    }
    virtual Void _imul(const X& c) override {
        static_cast<A&>(*this) *= c;
    }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& a) override {
        static_cast<A&>(*this) += c * dynamic_cast<const A&>(a);
    }
    virtual Void _ifma(const AlgebraInterface<X>& a1, const AlgebraInterface<X>& a2) override {
        static_cast<A&>(*this) += dynamic_cast<const A&>(a1) * dynamic_cast<const A&>(a2);
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(UnaryRingOperator op) const override {
        return heap_move(op(static_cast<const A&>(*this)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(BinaryRingOperator op, const AlgebraInterface<X>& a) const override {
        return heap_move(op(static_cast<const A&>(*this), dynamic_cast<const AlgebraWrapper<A, X>&>(a)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(BinaryFieldOperator op, const X& c) const override {
        return heap_move(op(static_cast<const A&>(*this), c));
    }
    /** @brief Apply the operation with reversed operand order. */
    virtual AlgebraInterface<X>* _rapply(BinaryRingOperator op, const X& c) const override {
        return heap_move(op(c, static_cast<const A&>(*this)));
    };
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(GradedRingOperator op, Nat m) const override {
        return heap_move(op(static_cast<const A&>(*this), m));
    };
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const override {
        os << static_cast<const A&>(*this);
        return os;
    }
};

template<class A, class X>
    requires AnInplaceAlgebraOver<A, X>
class AlgebraWrapper<A, X>
    : public virtual AlgebraInterface<X>,
      public A {
    /** @brief Return a dynamically allocated moved value. */
    static AlgebraInterface<X>* heap_move(A&& a) { return new AlgebraWrapper<A, X>(std::forward<A>(a)); }

  public:
    typedef X NumericType;

    /** @brief Move-construct the object from the supplied instance. */
    AlgebraWrapper(A&& a) : A(std::move(a)) {}
    /** @brief Construct a \c AlgebraWrapper from a. */
    AlgebraWrapper(const A& a) : A(a) {}
    virtual AlgebraInterface<X>* _copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual AlgebraInterface<X>* _create_zero() const override { return heap_move(this->A::create()); }
    virtual AlgebraInterface<X>* _create_constant(const X& c) const override { return heap_move(this->A::create_constant(c)); }
    virtual AlgebraInterface<X>* _create_copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual Void _iadd(const X& c) override { static_cast<A*>(this)->A::iadd(c); }
    virtual Void _imul(const X& c) override { static_cast<A*>(this)->A::imul(c); }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) override {
        static_cast<A*>(this)->A::isma(c, dynamic_cast<const A&>(x));
    }
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2) override {
        static_cast<A*>(this)->A::ifma(dynamic_cast<const A&>(x1), dynamic_cast<const A&>(x2));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(Neg) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(BinaryRingOperator, const AlgebraInterface<X>&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(BinaryFieldOperator, const X&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Apply the operation with reversed operand order. */
    virtual AlgebraInterface<X>* _rapply(BinaryRingOperator, const X&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(GradedRingOperator, Nat) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const override {
        os << static_cast<const A&>(*this);
        return os;
    }
};

/** @brief Adapter exposing normed algebra. */
template<class A, class X = NumericType<A>>
class NormedAlgebraWrapper
    : public virtual NormedAlgebraInterface<X>,
      public AlgebraWrapper<A, X> {
    using typename AlgebraInterface<X>::ErrorType;
    using AlgebraWrapper<A, X>::AlgebraWrapper;
    /** @brief Construct a value with the supplied error bound. */
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType r) const { return new A(static_cast<const A&>(*this).A::create_ball(r)); }
    /** @brief Construct a constant object with the supplied value. */
    virtual NormedAlgebraInterface<X>* _create_constant(X c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual NormedAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual NormedAlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
};

/** @brief Adapter exposing graded algebra. */
template<class A, class X = NumericType<A>>
class GradedAlgebraWrapper
    : public virtual GradedAlgebraInterface<X>,
      public AlgebraWrapper<A, X> {
    using AlgebraWrapper<A, X>::AlgebraWrapper;
    virtual GradedAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual GradedAlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    /** @brief Apply the supplied operation to the operands. */
    virtual GradedAlgebraInterface<X>* _apply(const Series<X>& f) const { return new A(compose(f, static_cast<const A&>(*this))); }
};

/** @brief Adapter exposing symbolic algebra. */
template<class A, class X = NumericType<A>>
class SymbolicAlgebraWrapper
    : public virtual SymbolicAlgebraInterface<X>,
      public AlgebraWrapper<A, X> {
    using AlgebraWrapper<A, X>::AlgebraWrapper;
    using SymbolicAlgebraInterface<X>::_apply;
    using AlgebraWrapper<A, X>::_apply;
    virtual SymbolicAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual SymbolicAlgebraInterface<X>* _create_copy() const { return new SymbolicAlgebraWrapper<A>(*this); }
    /** @brief Apply the supplied operation to the operands. */
    virtual SymbolicAlgebraInterface<X>* _apply(UnaryElementaryOperator op) { return new A(op, static_cast<const A&>(*this)); }
};

/** @brief Adapter exposing transcendental algebra. */
template<class A, class X>
class TranscendentalAlgebraWrapper
    : public TranscendentalAlgebraMixin<TranscendentalAlgebraWrapper<A, X>, X>,
      public A {
  public:
    using A::A;
    using TranscendentalAlgebraMixin<TranscendentalAlgebraWrapper<A, X>, X>::_apply;
    using TranscendentalAlgebraMixin<TranscendentalAlgebraWrapper<A, X>, X>::_rapply;
    /** @brief Construct a \c TranscendentalAlgebraWrapper from a. */
    TranscendentalAlgebraWrapper(const A& a) : A(a) {}
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const TranscendentalAlgebraWrapper<A, X>& a) { return os << static_cast<const A&>(a); }
};

/** @brief Adapter exposing elementary algebra. */
template<class A, class X>
class ElementaryAlgebraWrapper
    : public ElementaryAlgebraMixin<ElementaryAlgebraWrapper<A, X>, X>,
      public A {
  public:
    using A::A;
    using ElementaryAlgebraMixin<ElementaryAlgebraWrapper<A, X>, X>::_apply;
    using ElementaryAlgebraMixin<ElementaryAlgebraWrapper<A, X>, X>::_rapply;
    /** @brief Construct a \c ElementaryAlgebraWrapper from a. */
    ElementaryAlgebraWrapper(const A& a) : A(a) {}
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ElementaryAlgebraWrapper<A, X>& a) { return os << static_cast<const A&>(a); }
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_WRAPPER_HPP */
