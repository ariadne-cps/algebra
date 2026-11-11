/***************************************************************************
 *            algebra/algebra_mixin.hpp
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

/** @file algebra/algebra_mixin.hpp
 *  @brief Mixin operations shared by algebra implementations.
 */

#ifndef ARIADNE_ALGEBRA_MIXIN_HPP
#define ARIADNE_ALGEBRA_MIXIN_HPP

#include "algebra/algebra_concepts.hpp"
#include "algebra/algebra_interface.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

/** @brief Reusable implementation of algebra. */
template<class A, class X>
class AlgebraMixin
    : public virtual AlgebraInterface<X> {
    static_assert(AnInplaceAlgebraOver<A, X> or AnAlgebraOver<A, X>);
    typedef X NumericType;

  private:
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static AlgebraInterface<X>* _eval(OP op, const AlgebraMixin<A, X>& am1, const AlgebraInterface<X>& ai2) {
        const AlgebraMixin<A, X>* amp2 = dynamic_cast<const AlgebraMixin<A, X>*>(&ai2);
        assert(amp2);
        const A& a1 = static_cast<const A&>(am1);
        const A& a2 = static_cast<const A&>(*amp2);
        return new A(op(a1, a2));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static AlgebraInterface<X>* _eval(OP op, const AlgebraMixin<A, X>& am, const X& c) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a, c));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static AlgebraInterface<X>* _eval(OP op, const X& c, const AlgebraMixin<A, X>& am) {
        const A& a = static_cast<const A&>(am);
        return new A(op(c, a));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static AlgebraInterface<X>* _eval(OP op, const AlgebraMixin<A, X>& am) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static AlgebraInterface<X>* _eval(OP op, const AlgebraMixin<A, X>& am, Nat m) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a, m));
    }

  public:
    virtual AlgebraInterface<X>* _copy() const { return new A(static_cast<const A&>(*this)); }
    virtual AlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create_zero()); }
    virtual AlgebraInterface<X>* _create_constant(const X& c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual AlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    virtual Void _iadd(const X& c) {
        if constexpr (AnInplaceAlgebraOver<A, X>) {
            static_cast<A*>(this)->A::iadd(c);
        } else {
            static_cast<A&>(*this) += c;
        }
    }
    virtual Void _imul(const X& c) {
        if constexpr (AnInplaceAlgebraOver<A, X>) {
            static_cast<A*>(this)->A::imul(c);
        } else {
            static_cast<A&>(*this) *= c;
        }
    }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) {
        if constexpr (AnInplaceAlgebraOver<A, X>) {
            static_cast<A*>(this)->A::isma(c, dynamic_cast<const A&>(x));
        } else {
            static_cast<A&>(*this) += c * dynamic_cast<const A&>(x);
        }
    }
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2) {
        if constexpr (AnInplaceAlgebraOver<A, X>) {
            static_cast<A*>(this)->A::isma(dynamic_cast<const A&>(x1), dynamic_cast<const A&>(x2));
        } else {
            static_cast<A&>(*this) += dynamic_cast<const A&>(x1) * dynamic_cast<const A&>(x2);
        }
    }

    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(UnaryRingOperator op) const { return _eval(op, *this); }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(BinaryRingOperator op, const AlgebraInterface<X>& other) const { return _eval(op, *this, other); }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(BinaryFieldOperator op, const X& cnst) const { return _eval(op, *this, cnst); }
    /** @brief Apply the operation with reversed operand order. */
    virtual AlgebraInterface<X>* _rapply(BinaryRingOperator op, const X& cnst) const { return _eval(op, cnst, *this); }
    /** @brief Apply the supplied operation to the operands. */
    virtual AlgebraInterface<X>* _apply(GradedRingOperator op, Nat m) const { return _eval(op, *this, m); }

    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const { return os << *static_cast<const A*>(this); }
};

/** @brief Reusable implementation of transcendental algebra. */
template<class A, class X>
class TranscendentalAlgebraMixin
    : public virtual TranscendentalAlgebraInterface<X> {
    typedef X NumericType;

  private:
    /** @brief Convert to the concrete representation required by dispatch. */
    static const A& _cast(const TranscendentalAlgebraMixin<A, X>& am) {
        return static_cast<const A&>(am);
    }
    /** @brief Convert to the concrete representation required by dispatch. */
    static const A& _cast(const TranscendentalAlgebraInterface<X>& ai) {
        return static_cast<const A&>(dynamic_cast<const TranscendentalAlgebraMixin<A, X>&>(ai));
    }
    /** @brief Return a dynamically allocated moved value. */
    static A* _heap_move(A&& a) { return new A(std::move(a)); }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static TranscendentalAlgebraInterface<X>* _eval(OP op, const TranscendentalAlgebraMixin<A, X>& am1, const TranscendentalAlgebraInterface<X>& ai2) {
        const TranscendentalAlgebraMixin<A, X>* amp2 = dynamic_cast<const TranscendentalAlgebraMixin<A, X>*>(&ai2);
        assert(amp2);
        const A& a1 = static_cast<const A&>(am1);
        const A& a2 = static_cast<const A&>(*amp2);
        return new A(op(a1, a2));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static TranscendentalAlgebraInterface<X>* _eval(OP op, const TranscendentalAlgebraMixin<A, X>& am, const X& c) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a, c));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static TranscendentalAlgebraInterface<X>* _eval(OP op, const X& c, const TranscendentalAlgebraMixin<A, X>& am) {
        const A& a = static_cast<const A&>(am);
        return new A(op(c, a));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static TranscendentalAlgebraInterface<X>* _eval(OP op, const TranscendentalAlgebraMixin<A, X>& am) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static TranscendentalAlgebraInterface<X>* _eval(OP op, const TranscendentalAlgebraMixin<A, X>& am, Nat m) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a, m));
    }

  public:
    /** @brief Construct a zero object with compatible characteristics. */
    virtual TranscendentalAlgebraInterface<X>* _create_zero() const override {
        return new A(static_cast<const A&>(*this).A::create_zero());
    }
    /** @brief Construct a constant object with the supplied value. */
    virtual TranscendentalAlgebraInterface<X>* _create_constant(const X& c) const override {
        return new A(static_cast<const A&>(*this).A::create_constant(c));
    }
    /** @brief Construct a compatible copy of the supplied object. */
    virtual TranscendentalAlgebraInterface<X>* _create_copy() const override {
        return new A(static_cast<const A&>(*this));
    }

    /** @brief Apply the supplied operation to the operands. */
    virtual TranscendentalAlgebraInterface<X>* _apply(BinaryFieldOperator op, const TranscendentalAlgebraInterface<X>& other) const override {
        return _heap_move(op(_cast(*this), _cast(other)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual TranscendentalAlgebraInterface<X>* _apply(UnaryTranscendentalOperator op) const override {
        return _heap_move(op(_cast(*this)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual TranscendentalAlgebraInterface<X>* _apply(BinaryFieldOperator op, const X& cnst) const override {
        return _heap_move(op(_cast(*this), cnst));
    }
    /** @brief Apply the operation with reversed operand order. */
    virtual TranscendentalAlgebraInterface<X>* _rapply(BinaryFieldOperator op, const X& cnst) const override {
        return _heap_move(op(cnst, _cast(*this)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual TranscendentalAlgebraInterface<X>* _apply(GradedFieldOperator op, Int n) const override {
        return _heap_move(op(_cast(*this), n));
    }

    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const override { return os << *static_cast<const A*>(this); }
};

/** @brief Reusable implementation of elementary algebra. */
template<class A, class X>
class ElementaryAlgebraMixin
    : public virtual ElementaryAlgebraInterface<X> {
    typedef X NumericType;

  private:
    /** @brief Convert to the concrete representation required by dispatch. */
    static const A& _cast(const ElementaryAlgebraMixin<A, X>& am) {
        return static_cast<const A&>(am);
    }
    /** @brief Convert to the concrete representation required by dispatch. */
    static const A& _cast(const ElementaryAlgebraInterface<X>& ai) {
        return static_cast<const A&>(dynamic_cast<const ElementaryAlgebraMixin<A, X>&>(ai));
    }
    /** @brief Return a dynamically allocated moved value. */
    static A* _heap_move(A&& a) { return new A(std::move(a)); }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static ElementaryAlgebraInterface<X>* _eval(OP op, const ElementaryAlgebraMixin<A, X>& am1, const ElementaryAlgebraInterface<X>& ai2) {
        const ElementaryAlgebraMixin<A, X>* amp2 = dynamic_cast<const ElementaryAlgebraMixin<A, X>*>(&ai2);
        assert(amp2);
        const A& a1 = static_cast<const A&>(am1);
        const A& a2 = static_cast<const A&>(*amp2);
        return new A(op(a1, a2));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static ElementaryAlgebraInterface<X>* _eval(OP op, const ElementaryAlgebraMixin<A, X>& am, const X& c) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a, c));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static ElementaryAlgebraInterface<X>* _eval(OP op, const X& c, const ElementaryAlgebraMixin<A, X>& am) {
        const A& a = static_cast<const A&>(am);
        return new A(op(c, a));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static ElementaryAlgebraInterface<X>* _eval(OP op, const ElementaryAlgebraMixin<A, X>& am) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a));
    }
    /** @brief Evaluate the expansion at the supplied argument. */
    template<class OP>
    static ElementaryAlgebraInterface<X>* _eval(OP op, const ElementaryAlgebraMixin<A, X>& am, Nat m) {
        const A& a = static_cast<const A&>(am);
        return new A(op(a, m));
    }

  public:
    /** @brief Construct a zero object with compatible characteristics. */
    virtual ElementaryAlgebraInterface<X>* _create_zero() const override {
        return new A(static_cast<const A&>(*this).A::create_zero());
    }
    /** @brief Construct a constant object with the supplied value. */
    virtual ElementaryAlgebraInterface<X>* _create_constant(const X& c) const override {
        return new A(static_cast<const A&>(*this).A::create_constant(c));
    }
    /** @brief Construct a compatible copy of the supplied object. */
    virtual ElementaryAlgebraInterface<X>* _create_copy() const override {
        return new A(static_cast<const A&>(*this));
    }

    /** @brief Apply the supplied operation to the operands. */
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryElementaryOperator op, const ElementaryAlgebraInterface<X>& other) const override {
        return _heap_move(op(_cast(*this), _cast(other)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual ElementaryAlgebraInterface<X>* _apply(UnaryElementaryOperator op) const override {
        return _heap_move(op(_cast(*this)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryElementaryOperator op, const X& cnst) const override {
        return _heap_move(op(_cast(*this), cnst));
    }
    /** @brief Apply the operation with reversed operand order. */
    virtual ElementaryAlgebraInterface<X>* _rapply(BinaryElementaryOperator op, const X& cnst) const override {
        return _heap_move(op(cnst, _cast(*this)));
    }
    /** @brief Apply the supplied operation to the operands. */
    virtual ElementaryAlgebraInterface<X>* _apply(GradedElementaryOperator op, Int n) const override {
        return _heap_move(op(_cast(*this), n));
    }

    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const override { return os << *static_cast<const A*>(this); }
};

/** @brief Reusable implementation of normed algebra. */
template<class A, class X>
class NormedAlgebraMixin
    : public virtual NormedAlgebraInterface<X>,
      public AlgebraMixin<A, X> {
    using typename NormedAlgebraInterface<X>::ErrorType;

    /** @brief Construct a value with the supplied error bound. */
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType r) const { return new A(static_cast<const A&>(*this).A::create_ball(r)); }
    /** @brief Construct a constant object with the supplied value. */
    virtual NormedAlgebraInterface<X>* _create_constant(X c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual NormedAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual NormedAlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
};

/** @brief Reusable implementation of graded algebra. */
template<class A, class X>
class GradedAlgebraMixin
    : public virtual GradedAlgebraInterface<X>,
      public AlgebraMixin<A, X> {
    virtual GradedAlgebraMixin<A, X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual GradedAlgebraMixin<A, X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    /** @brief Apply the supplied operation to the operands. */
    virtual GradedAlgebraMixin<A, X>* _apply(const Series<X>& f) const { return new A(compose(f, static_cast<const A&>(*this))); }
};

/** @brief Reusable implementation of symbolic algebra. */
template<class A, class X>
class SymbolicAlgebraMixin
    : public virtual SymbolicAlgebraInterface<X>,
      public AlgebraMixin<A, X> {
    using SymbolicAlgebraInterface<X>::_apply;
    using AlgebraMixin<A, X>::_apply;
    virtual SymbolicAlgebraMixin<A, X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    virtual SymbolicAlgebraMixin<A, X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual SymbolicAlgebraMixin<A, X>* _create_constant(const X& c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    /** @brief Apply the supplied operation to the operands. */
    virtual SymbolicAlgebraMixin<A, X>* _apply(UnaryElementaryOperator op) { return new A(op, static_cast<const A&>(*this)); }
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_MIXIN_HPP */
