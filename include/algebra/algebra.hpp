/***************************************************************************
 *            algebra/algebra.hpp
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

/** @file algebra/algebra.hpp
 *  @brief Algebraic structures and normed algebra operations.
 */
#ifndef ARIADNE_ALGEBRA_HPP
#define ARIADNE_ALGEBRA_HPP

#include <iosfwd>
#include <iostream>

#include "algebra/algebra_interface.hpp"
#include "algebra/algebra_wrapper.hpp"
#include "algebra/operations.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.hpp"
#include "utility/exceptions.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/pointer.hpp"

namespace Ariadne {

/** @brief Operations supported by algebra objects. */
template<class A, class X>
struct AlgebraOperations;

/** @brief Operations supported by algebra objects. */
template<class X>
struct AlgebraOperations<Algebra<X>, X> {
    /** @brief Apply the supplied operation to the operands. */
    static Algebra<X> apply(UnaryRingOperator op, Algebra<X> a) { return Algebra<X>(a.managed_pointer()->_apply(op)); }
    /** @brief Apply the supplied operation to the operands. */
    static Algebra<X> apply(BinaryRingOperator op, Algebra<X> a1, Algebra<X> a2) { return Algebra<X>(a1.managed_pointer()->_apply(op, *a2.managed_pointer())); }
    /** @brief Apply the supplied operation to the operands. */
    static Algebra<X> apply(BinaryFieldOperator op, Algebra<X> a1, const X& c2) { return Algebra<X>(a1.managed_pointer()->_apply(op, c2)); }
    /** @brief Apply the supplied operation to the operands. */
    static Algebra<X> apply(BinaryRingOperator op, const X& c1, Algebra<X> a2) { return Algebra<X>(a2.managed_pointer()->_rapply(op, c1)); }
    /** @brief Apply the supplied operation to the operands. */
    static Algebra<X> apply(GradedRingOperator op, Algebra<X> a1, Nat m2) { return Algebra<X>(a1.managed_pointer()->_apply(op, m2)); }
};

/** @brief Operations supported by algebra objects. */
template<class X>
struct AlgebraOperations<TranscendentalAlgebra<X>, X> {
    /** @brief Convert to the common interface or value type. */
    static const TranscendentalAlgebraInterface<X>* _upcast(const AlgebraInterface<X>* ap) {
        return dynamic_cast<const TranscendentalAlgebraInterface<X>*>(ap);
    }
    /** @brief Convert to the common interface or value type. */
    static const TranscendentalAlgebraInterface<X>* _upcast(const TranscendentalAlgebraInterface<X>* ap) {
        return ap;
    }
    /** @brief Apply the supplied operation to the operands. */
    static TranscendentalAlgebra<X> apply(BinaryFieldOperator op, TranscendentalAlgebra<X> a1, TranscendentalAlgebra<X> a2) {
        return TranscendentalAlgebra<X>(_upcast(a1.raw_pointer())->_apply(op, a2.reference()));
    }
    /** @brief Apply the supplied operation to the operands. */
    static TranscendentalAlgebra<X> apply(BinaryFieldOperator op, TranscendentalAlgebra<X> a1, const X& c2) {
        return TranscendentalAlgebra<X>(_upcast(a1.raw_pointer())->_apply(op, c2));
    }
    /** @brief Apply the supplied operation to the operands. */
    static TranscendentalAlgebra<X> apply(BinaryFieldOperator op, const X& c1, TranscendentalAlgebra<X> a2) {
        return TranscendentalAlgebra<X>(_upcast(a2.raw_pointer())->_rapply(op, c1));
    }
    /** @brief Apply the supplied operation to the operands. */
    static TranscendentalAlgebra<X> apply(UnaryTranscendentalOperator op, TranscendentalAlgebra<X> a) {
        return TranscendentalAlgebra<X>(_upcast(a.raw_pointer())->_apply(op));
    }
    /** @brief Apply the supplied operation to the operands. */
    static TranscendentalAlgebra<X> apply(Pow op, TranscendentalAlgebra<X> a1, Int n2) {
        return TranscendentalAlgebra<X>(_upcast(a1.raw_pointer())->_apply(op, n2));
    }
};

/** @brief Operations supported by algebra objects. */
template<class X>
struct AlgebraOperations<ElementaryAlgebra<X>, X> {
    /** @brief Convert to the common interface or value type. */
    static const ElementaryAlgebraInterface<X>* _upcast(const AlgebraInterface<X>* ap) {
        return dynamic_cast<const ElementaryAlgebraInterface<X>*>(ap);
    }
    /** @brief Convert to the common interface or value type. */
    static const ElementaryAlgebraInterface<X>* _upcast(const ElementaryAlgebraInterface<X>* ap) {
        return ap;
    }
    /** @brief Apply the supplied operation to the operands. */
    static ElementaryAlgebra<X> apply(BinaryElementaryOperator op, ElementaryAlgebra<X> a1, ElementaryAlgebra<X> a2) {
        return ElementaryAlgebra<X>(_upcast(a1.raw_pointer())->_apply(op, a2.reference()));
    }
    /** @brief Apply the supplied operation to the operands. */
    static ElementaryAlgebra<X> apply(BinaryElementaryOperator op, ElementaryAlgebra<X> a1, const X& c2) {
        return ElementaryAlgebra<X>(_upcast(a1.raw_pointer())->_apply(op, c2));
    }
    /** @brief Apply the supplied operation to the operands. */
    static ElementaryAlgebra<X> apply(BinaryElementaryOperator op, const X& c1, ElementaryAlgebra<X> a2) {
        return ElementaryAlgebra<X>(_upcast(a2.raw_pointer())->_rapply(op, c1));
    }
    /** @brief Apply the supplied operation to the operands. */
    static ElementaryAlgebra<X> apply(UnaryElementaryOperator op, ElementaryAlgebra<X> a) {
        return ElementaryAlgebra<X>(_upcast(a.raw_pointer())->_apply(op));
    }
    /** @brief Apply the supplied operation to the operands. */
    static ElementaryAlgebra<X> apply(GradedElementaryOperator op, ElementaryAlgebra<X> a1, Int n2) {
        return ElementaryAlgebra<X>(_upcast(a1.raw_pointer())->_apply(op, n2));
    }
};

//! \brief Generic class for elements of unital algebras.
template<class X>
class Algebra
    : public Handle<const AlgebraInterface<X>>,
      public DispatchAlgebraOperations<Algebra<X>, X> {
  public:
    typedef AlgebraInterface<X> Interface;
    typedef X ScalarType;
    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;

  public:
    using Handle<const Interface>::Handle;
    /** @brief Construct a \c Algebra with the specified characteristics. */
    template<AnAlgebraOver<X> A>
    Algebra(const A&);
    /** @brief Extract the requested concrete representation. */
    template<class A>
        requires BaseOf<AlgebraInterface<X>, A>
    A extract() const;
    /** @brief Extract the requested concrete representation. */
    template<class A>
        requires(!BaseOf<AlgebraInterface<X>, A>)
    A extract() const;
    /** @brief Copy-construct the object from the supplied instance. */
    Algebra(const Algebra<X>& a) : Handle<const AlgebraInterface<X>>(a) {}
    /** @brief Assign the supplied value and return this object. */
    Algebra<X>& operator=(const Algebra<X>& a) {
        (*this).Handle<const AlgebraInterface<X>>::operator=(std::shared_ptr<const AlgebraInterface<X>>(a.managed_pointer()->_create_copy()));
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    Algebra<X>& operator=(const X& c) { return *this = this->create_constant(c); }
    /** @brief Convert to const algebra interface< x > &. */
    operator const AlgebraInterface<X>&() const { return *this->managed_pointer(); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    Algebra<X> create() const { return Algebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    Algebra<X> clone() const { return Algebra<X>(this->managed_pointer()->_create_copy()); }
    /** @brief Construct a zero object with compatible characteristics. */
    Algebra<X> create_zero() const { return Algebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Construct a constant object with the supplied value. */
    Algebra<X> create_constant(const X& c) const { return Algebra<X>(this->managed_pointer()->_create_constant(c)); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Algebra<X>& a) { return a.managed_pointer()->_write(os); }

  private:
    /** @brief Construct algebra from the supplied arguments. */
    template<class A>
        requires BaseOf<AlgebraInterface<X>, A>
    static AlgebraInterface<X>* make_algebra(const A& a) { return new A(a); }
    /** @brief Construct algebra from the supplied arguments. */
    template<class A>
        requires(!BaseOf<AlgebraInterface<X>, A>)
    static AlgebraInterface<X>* make_algebra(const A& a) { return new AlgebraWrapper<A, X>(a); }
};

//! \brief Generic class for elements of unital algebras.
template<class X>
class TranscendentalAlgebra
    : public Handle<const TranscendentalAlgebraInterface<X>>,
      public DispatchTranscendentalAlgebraOperations<TranscendentalAlgebra<X>, X> {
  public:
    typedef TranscendentalAlgebraInterface<X> Interface;
    typedef X ScalarType;
    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;

  public:
    using Handle<const Interface>::Handle;
    /** @brief Construct a \c TranscendentalAlgebra with the specified characteristics. */
    TranscendentalAlgebra() : TranscendentalAlgebra(nullptr) {}
    /** @brief Construct a \c TranscendentalAlgebra from a. */
    template<ATranscendentalAlgebraOver<X> A>
    explicit TranscendentalAlgebra(const A& a);
    /** @brief Extract the requested concrete representation. */
    template<class A>
        requires BaseOf<TranscendentalAlgebraInterface<X>, A>
    A extract() const;
    /** @brief Extract the requested concrete representation. */
    template<class A>
        requires(!BaseOf<TranscendentalAlgebraInterface<X>, A>)
    A extract() const;
    /** @brief Assign the supplied value and return this object. */
    TranscendentalAlgebra<X>& operator=(const X& c) { return *this = this->create_constant(c); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    TranscendentalAlgebra<X> create() const { return TranscendentalAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    TranscendentalAlgebra<X> clone() const { return TranscendentalAlgebra<X>(this->managed_pointer()->_create_copy()); }
    /** @brief Construct a zero object with compatible characteristics. */
    TranscendentalAlgebra<X> create_zero() const { return TranscendentalAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Construct a constant object with the supplied value. */
    TranscendentalAlgebra<X> create_constant(const X& c) const { return TranscendentalAlgebra<X>(this->managed_pointer()->_create_constant(c)); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const TranscendentalAlgebra<X>& a) { return a.managed_pointer()->_write(os); }

  private:
    /** @brief Construct elementary algebra from the supplied arguments. */
    template<class A>
        requires BaseOf<TranscendentalAlgebraInterface<X>, A>
    static TranscendentalAlgebraInterface<X>* make_elementary_algebra(const A& a) { return new A(a); }
    /** @brief Construct elementary algebra from the supplied arguments. */
    template<class A>
        requires(!BaseOf<TranscendentalAlgebraInterface<X>, A>)
    static TranscendentalAlgebraInterface<X>* make_elementary_algebra(const A& a) { return new TranscendentalAlgebraWrapper<A, X>(a); }
};

//! \brief Generic class for elements of unital algebras.
template<class X>
class ElementaryAlgebra
    : public Handle<const ElementaryAlgebraInterface<X>>,
      public DispatchElementaryAlgebraOperations<ElementaryAlgebra<X>, X> {
  public:
    typedef ElementaryAlgebraInterface<X> Interface;
    typedef X ScalarType;
    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;

  public:
    using Handle<const Interface>::Handle;
    /** @brief Construct a \c ElementaryAlgebra with the specified characteristics. */
    ElementaryAlgebra() : ElementaryAlgebra(nullptr) {}
    /** @brief Construct a \c ElementaryAlgebra from a. */
    template<AnElementaryAlgebraOver<X> A>
    explicit ElementaryAlgebra(const A& a);
    /** @brief Extract the requested concrete representation. */
    template<class A>
        requires BaseOf<ElementaryAlgebraInterface<X>, A>
    A extract() const;
    /** @brief Extract the requested concrete representation. */
    template<class A>
        requires(!BaseOf<ElementaryAlgebraInterface<X>, A>)
    A extract() const;
    /** @brief Assign the supplied value and return this object. */
    ElementaryAlgebra<X>& operator=(const X& c) { return *this = this->create_constant(c); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    ElementaryAlgebra<X> create() const { return ElementaryAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    ElementaryAlgebra<X> clone() const { return ElementaryAlgebra<X>(this->managed_pointer()->_create_copy()); }
    /** @brief Construct a zero object with compatible characteristics. */
    ElementaryAlgebra<X> create_zero() const { return ElementaryAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Construct a constant object with the supplied value. */
    ElementaryAlgebra<X> create_constant(const X& c) const { return ElementaryAlgebra<X>(this->managed_pointer()->_create_constant(c)); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ElementaryAlgebra<X>& a) { return a.managed_pointer()->_write(os); }

  private:
    /** @brief Construct elementary algebra from the supplied arguments. */
    template<class A>
        requires BaseOf<ElementaryAlgebraInterface<X>, A>
    static ElementaryAlgebraInterface<X>* make_elementary_algebra(const A& a) { return new A(a); }
    /** @brief Construct elementary algebra from the supplied arguments. */
    template<class A>
        requires(!BaseOf<ElementaryAlgebraInterface<X>, A>)
    static ElementaryAlgebraInterface<X>* make_elementary_algebra(const A& a) { return new ElementaryAlgebraWrapper<A, X>(a); }
};

/** @brief Polymorphic algebra element with a norm and error-aware operations. */
template<class X>
class NormedAlgebra
    : public Handle<const NormedAlgebraInterface<X>> {
    static_assert(not Same<X, Real>);

  public:
    typedef NormedAlgebraInterface<X> Interface;
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;

    using Handle<const Interface>::Handle;
    /** @brief Construct a \c NormedAlgebra from a. */
    NormedAlgebra(const Algebra<X>& a) : Handle<const Interface>(std::dynamic_pointer_cast(a.managed_pointer())) {}
    /** @brief Convert to algebra< x >. */
    operator Algebra<X>() const { return this->managed_pointer()->_create_copy(); }
    /** @brief Assign the supplied value and return this object. */
    NormedAlgebra<X>& operator=(Int c) {
        *this = this->create();
        this->iadd(c);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    NormedAlgebra<X>& operator=(const X& c) {
        *this = this->create();
        this->iadd(c);
        return *this;
    }
    //! \brief The norm of the element.
    typename AlgebraTraits<X>::NormType norm() const { return this->managed_pointer()->norm(); }
    //! \brief A value \a c minimising |a-c|.
    typename AlgebraTraits<X>::ValueType average() const { return this->managed_pointer()->average(); }
    //! \brief The tolerance used to determine the truncation error when applying an analytic function.
    FloatDP tolerance() const { return this->managed_pointer()->tolerance(); }
    //! \brief A value \c r such that \c |a-c1|<=r.
    typename AlgebraTraits<X>::NormType radius() const { return this->managed_pointer()->radius(); }
    /** @brief Construct a zero object with compatible characteristics. */
    NormedAlgebra<X> create_zero() const { return NormedAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Construct a constant object with the supplied value. */
    NormedAlgebra<X> create_constant(X c) const { return NormedAlgebra<X>(this->managed_pointer()->_create_constant(c)); }
    /** @brief Construct a value with the supplied error bound. */
    NormedAlgebra<X> create_ball(FloatDPError r) const { return NormedAlgebra<X>(this->managed_pointer()->_create_ball(r)); }
    /** @brief Construct a value with the supplied error bound. */
    NormedAlgebra<X> create_ball(FloatDPApproximation r) const {
        return create_ball(FloatDPError(r.raw()));
    }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    NormedAlgebra<X> create() const { return NormedAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    NormedAlgebra<X> clone() const { return NormedAlgebra<X>(this->managed_pointer()->_create_copy()); }
    /** @brief Remove all elements. */
    Void clear() { this->imul(0); }
    /** @brief Add the operand in place. */
    Void iadd(const X& c) { this->managed_pointer()->_iadd(c); }
    /** @brief Multiply by the operand in place. */
    Void imul(const X& c) { this->managed_pointer()->_imul(c); }
    /** @brief Accumulate the scaled value in place. */
    Void isma(const X& c, const NormedAlgebra<X>& x) { this->managed_pointer()->_isma(c, *x.managed_pointer()); }
    /** @brief Accumulate the product in place. */
    Void ifma(const NormedAlgebra<X>& x1, const NormedAlgebra<X>& x2) { this->managed_pointer()->_ifma(*x1.managed_pointer(), *x2.managed_pointer()); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const NormedAlgebra<X>& a) { return a.managed_pointer()->_write(os); }
};

//! \brief Generic class for elements of unital algebras.
template<class X>
class GradedAlgebra
    : public Handle<const GradedAlgebraInterface<X>> {
  public:
    typedef GradedAlgebraInterface<X> Interface;
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;

    using Handle<const Interface>::Handle;
    /** @brief Assign the supplied value and return this object. */
    GradedAlgebra<X>& operator=(Int c) {
        *this = this->create();
        this->iadd(c);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    GradedAlgebra<X>& operator=(const X& c) {
        *this = this->create();
        this->iadd(c);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    GradedAlgebra<X>& operator=(const Algebra<X>& a) {
        this->managed_pointer() = std::shared_ptr<AlgebraInterface<X>>(a.managed_pointer()->_create_copy());
        return *this;
    }
    /** @brief Convert to algebra< x >. */
    operator Algebra<X>() const { return this->managed_pointer()->_create_copy(); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    GradedAlgebra<X> create() const { return GradedAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    GradedAlgebra<X> clone() const { return GradedAlgebra<X>(this->managed_pointer()->_create_copy()); }
    //! \brief The degree of the algebra.
    Nat degree() const { return this->managed_pointer()->degree(); }
    //! \brief The value in the null grade.
    const X& value() const { return this->managed_pointer()->value(); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const GradedAlgebra<X>& a) { return a.managed_pointer()->_write(os); }

  public:
    /** @brief Add the operand in place. */
    Void iadd(const X& c) { this->managed_pointer()->_iadd(c); }
    /** @brief Multiply by the operand in place. */
    Void imul(const X& c) { this->managed_pointer()->_imul(c); }
    /** @brief Accumulate the scaled value in place. */
    Void isma(const X& c, const GradedAlgebra<X>& x) { this->managed_pointer()->_isma(c, *x.managed_pointer()); }
    /** @brief Accumulate the product in place. */
    Void ifma(const GradedAlgebra<X>& x1, const GradedAlgebra<X>& x2) { this->managed_pointer()->_ifma(*x1.managed_pointer(), *x2.managed_pointer()); }
};

//! \brief Generic class for elements of unital algebras.
template<class X>
class SymbolicAlgebra
    : public Handle<const SymbolicAlgebraInterface<X>> {
  public:
    typedef SymbolicAlgebraInterface<X> Interface;
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;

    using Handle<const Interface>::Handle;
    //! \brief Create the representation of the operator \a op applied to \a a.
    SymbolicAlgebra(UnaryElementaryOperator op, const SymbolicAlgebra<X>& a) : Handle<const Interface>(a.managed_pointer()->_apply(op)) {}
    /** @brief Assign the supplied value and return this object. */
    SymbolicAlgebra<X>& operator=(Int c) {
        *this = this->create();
        this->iadd(c);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    SymbolicAlgebra<X>& operator=(const X& c) {
        *this = this->create();
        this->iadd(c);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    SymbolicAlgebra<X>& operator=(const Algebra<X>& a) {
        this->managed_pointer() = std::shared_ptr<AlgebraInterface<X>>(a.managed_pointer()->_create_copy());
        return *this;
    }
    /** @brief Convert to algebra< x >. */
    operator Algebra<X>() const { return this->managed_pointer()->_create_copy(); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    SymbolicAlgebra<X> create() const { return SymbolicAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Construct a zero object with compatible characteristics. */
    SymbolicAlgebra<X> create_zero() const { return SymbolicAlgebra<X>(this->managed_pointer()->_create_zero()); }
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    SymbolicAlgebra<X> clone() const { return SymbolicAlgebra<X>(this->managed_pointer()->_create_copy()); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const SymbolicAlgebra<X>& a) { return a.managed_pointer()->_write(os); }

  public:
    /** @brief Add the operand in place. */
    Void iadd(const X& c) { this->managed_pointer()->_iadd(c); }
    /** @brief Multiply by the operand in place. */
    Void imul(const X& c) { this->managed_pointer()->_imul(c); }
    /** @brief Accumulate the scaled value in place. */
    Void isma(const X& c, const SymbolicAlgebra<X>& x) { this->managed_pointer()->_isma(c, *x.managed_pointer()); }
    /** @brief Accumulate the product in place. */
    Void ifma(const SymbolicAlgebra<X>& x1, const SymbolicAlgebra<X>& x2) { this->managed_pointer()->_ifma(*x1.managed_pointer(), *x2.managed_pointer()); }
};

template<class X>
template<AnAlgebraOver<X> A>
Algebra<X>::Algebra(const A& a)
    : Handle<const AlgebraInterface<X>>(make_algebra(a)) {}

template<class X>
template<class A>
    requires BaseOf<AlgebraInterface<X>, A>
A Algebra<X>::extract() const {
    auto ap = dynamic_cast<const A*>(this->raw_pointer());
    ARIADNE_PRECONDITION(ap);
    return *ap;
}

template<class X>
template<class A>
    requires(!BaseOf<AlgebraInterface<X>, A>)
A Algebra<X>::extract() const {
    auto awp = dynamic_cast<const AlgebraWrapper<A, X>*>(this->raw_pointer());
    ARIADNE_PRECONDITION(awp);
    return *awp;
}

template<class X>
template<ATranscendentalAlgebraOver<X> A>
TranscendentalAlgebra<X>::TranscendentalAlgebra(const A& a)
    : Handle<const TranscendentalAlgebraInterface<X>>(make_elementary_algebra(a)) {}

template<class X>
template<class A>
    requires BaseOf<TranscendentalAlgebraInterface<X>, A>
A TranscendentalAlgebra<X>::extract() const {
    auto ap = dynamic_cast<const A*>(this->raw_pointer());
    ARIADNE_PRECONDITION(ap);
    return *ap;
}

template<class X>
template<class A>
    requires(!BaseOf<TranscendentalAlgebraInterface<X>, A>)
A TranscendentalAlgebra<X>::extract() const {
    auto awp = dynamic_cast<const TranscendentalAlgebraWrapper<A, X>*>(this->raw_pointer());
    ARIADNE_PRECONDITION(awp);
    return *awp;
}

template<class X>
template<AnElementaryAlgebraOver<X> A>
ElementaryAlgebra<X>::ElementaryAlgebra(const A& a)
    : Handle<const ElementaryAlgebraInterface<X>>(make_elementary_algebra(a)) {}

template<class X>
template<class A>
    requires BaseOf<ElementaryAlgebraInterface<X>, A>
A ElementaryAlgebra<X>::extract() const {
    auto ap = dynamic_cast<const A*>(this->raw_pointer());
    ARIADNE_PRECONDITION(ap);
    return *ap;
}

template<class X>
template<class A>
    requires(!BaseOf<ElementaryAlgebraInterface<X>, A>)
A ElementaryAlgebra<X>::extract() const {
    auto awp = dynamic_cast<const ElementaryAlgebraWrapper<A, X>*>(this->raw_pointer());
    ARIADNE_PRECONDITION(awp);
    return *awp;
}

/** @brief Construction characteristics obtained from a prototype coefficient. */
template<class T>
class PrototypeCharacteristics {
    T _prototype;

  public:
    /** @brief Construct a \c PrototypeCharacteristics from prototype. */
    PrototypeCharacteristics(T prototype) : _prototype(prototype.create_zero()) {}
    /** @brief Convert to const t &. */
    operator const T&() const { return _prototype; }
    /** @brief Test equality of the operands. */
    friend Bool operator==(const PrototypeCharacteristics<T>&, const PrototypeCharacteristics<T>&) {
        ARIADNE_NOT_IMPLEMENTED;
    }
};

/** @brief Return the characteristics required to construct compatible values. */
template<class X>
PrototypeCharacteristics<ElementaryAlgebra<X>> characteristics(const ElementaryAlgebra<X>& a) { return PrototypeCharacteristics(a); }

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_HPP */
