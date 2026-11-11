/***************************************************************************
 *            algebra/expansion.inl.hpp
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

/** @file algebra/expansion.inl.hpp
 *  @brief Inline operations and helpers for sparse expansions.
 */
#ifndef ARIADNE_EXPANSION_INL_HPP
#define ARIADNE_EXPANSION_INL_HPP

#include "expansion.hpp"
#include "multi_index.hpp"
#include "utility/iterator.hpp"
#include "utility/macros.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

/************ Expansion ******************************************************/

/** @brief Associated type definitions for reference. */
template<class T>
struct ReferenceTypedef {
    typedef T& Type;
};
template<class T>
using ReferenceType = typename ReferenceTypedef<T>::Type;
template<class T>
using ConstReferenceType = typename ReferenceTypedef<const T>::Type;

/** @brief Associated type definitions for pointer. */
template<class T>
struct PointerTypedef {
    typedef T* Type;
};
template<class T>
using PointerType = typename PointerTypedef<T>::Type;
template<class T>
using ConstPointerType = typename PointerTypedef<const T>::Type;

/** @brief Associated type definitions for reference. */
template<>
struct ReferenceTypedef<MultiIndex> {
    typedef MultiIndexReference Type;
};
/** @brief Associated type definitions for reference. */
template<>
struct ReferenceTypedef<const MultiIndex> {
    typedef MultiIndexConstReference Type;
};

/** @brief Associated type definitions for pointer. */
template<>
struct PointerTypedef<MultiIndex> {
    typedef MultiIndexPointer Type;
};
/** @brief Associated type definitions for pointer. */
template<>
struct PointerTypedef<const MultiIndex> {
    typedef MultiIndexConstPointer Type;
};

/** @brief Comparison policy for expansion indices. */
template<class CMP>
struct IndexComparison {
    /** @brief Apply the callable to the supplied arguments. */
    template<class M1, class M2>
    bool operator()(M1 const& m1, M2 const& m2) const {
        CMP cmp;
        return cmp(m1.index(), m2.index());
    }
    /** @brief Apply the callable to the supplied arguments. */
    template<class M>
    bool operator()(const M& m, const typename M::IndexType& a) const {
        CMP cmp;
        return cmp(m.index(), a);
    }
};

/** @brief Ordering predicate using index order. */
struct IndexLess : public IndexComparison<Less> {};

/** @brief Ordering predicate using graded index order. */
struct GradedIndexLess : public IndexComparison<GradedLess> {};

/** @brief Ordering predicate using lexicographic index order. */
struct LexicographicIndexLess : public IndexComparison<LexicographicLess> {};

/** @brief Ordering predicate using reverse lexicographic index order. */
struct ReverseLexicographicIndexLess : public IndexComparison<ReverseLexicographicLess> {};

/** @brief Ordering predicate using coefficient order. */
struct CoefficientLess {
    /** @brief Apply the callable to the supplied arguments. */
    template<class M1, class M2>
    bool operator()(M1 const& m1, M2 const& m2) const {
        return m1.coefficient() < m2.coefficient();
    }
};

/** @brief Predicate selecting expansion entries with zero coefficients. */
struct CoefficientIsZero {
    /** @brief Apply the callable to the supplied arguments. */
    template<class M>
    bool operator()(const M& m) const {
        return decide(m.coefficient() == 0);
    }
};

/** @brief Owning index-coefficient pair from an expansion. */
template<class I, class X>
class ExpansionValue {
    I _a;
    X _c;

  public:
    typedef I IndexType;
    typedef X CoefficientType;

  public:
    /** @brief Construct a \c ExpansionValue from a, c. */
    ExpansionValue(const IndexType& a, const CoefficientType& c) : _a(a), _c(c) {}
    /** @brief Return the selected index. */
    IndexType& index() { return _a; }
    /** @brief Return the selected index. */
    const IndexType& index() const { return _a; }
    /** @brief Return the selected coefficient. */
    CoefficientType& coefficient() { return _c; }
    /** @brief Return the selected coefficient. */
    const CoefficientType& coefficient() const { return _c; }
    /** @brief Test equality of the operands. */
    decltype(auto) operator==(const ExpansionValue<I, X>& other) const { return this->_a == other._a && this->_c == other._c; }
    /** @brief Test inequality of the operands. */
    decltype(auto) operator!=(const ExpansionValue<I, X>& other) const { return !(*this == other); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ExpansionValue<I, X>& m) {
        return os << m.index() << ":" << m.coefficient();
    }
};

template<class I, class X>
class ExpansionReference;
template<class I, class X>
class ExpansionConstReference;

/** @brief Reference proxy into expansion storage. */
template<class I, class X>
class ExpansionReference {
  public:
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename Expansion<I, X>::IndexReference IndexReference;
    typedef typename Expansion<I, X>::CoefficientReference CoefficientReference;
    typedef typename Expansion<I, X>::IndexConstReference IndexConstReference;
    typedef typename Expansion<I, X>::CoefficientConstReference CoefficientConstReference;

  private:
    IndexReference _a;
    CoefficientReference _c;

  public:
    /** @brief Construct a \c ExpansionReference from a, c. */
    ExpansionReference(IndexReference a, CoefficientReference c) : _a(a), _c(c) {}
    /** @brief Return the selected index. */
    IndexReference index() { return _a; }
    /** @brief Return the selected index. */
    IndexConstReference index() const { return _a; }
    /** @brief Return the selected coefficient. */
    CoefficientReference coefficient() { return _c; }
    /** @brief Return the selected coefficient. */
    CoefficientConstReference coefficient() const { return _c; }
    /** @brief Convert to expansion value< i, x >. */
    operator ExpansionValue<I, X>() const { return ExpansionValue<I, X>(_a, _c); }
    /** @brief Assign the supplied value and return this object. */
    ExpansionReference<I, X>& operator=(const ExpansionValue<I, X>&);
    /** @brief Assign the supplied value and return this object. */
    ExpansionReference<I, X>& operator=(const ExpansionReference<I, X>&);
    /** @brief Assign the supplied value and return this object. */
    ExpansionReference<I, X>& operator=(const ExpansionConstReference<I, X>&);
    /** @brief Exchange the contents of the operands. */
    friend Void swap(ExpansionReference<I, X> m1, ExpansionReference<I, X> m2) {
        swap(m1.index(), m2.index());
        std::swap(m1.coefficient(), m2.coefficient());
    }
};

/** @brief Reference proxy into expansion const storage. */
template<class I, class X>
class ExpansionConstReference {
  public:
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename Expansion<I, X>::IndexConstReference IndexConstReference;
    typedef typename Expansion<I, X>::CoefficientConstReference CoefficientConstReference;

  private:
    IndexConstReference _a;
    CoefficientConstReference _c;

  public:
    /** @brief Construct a \c ExpansionConstReference from a, c. */
    ExpansionConstReference(IndexConstReference a, CoefficientConstReference c) : _a(a), _c(c) {}
    /** @brief Return the selected index. */
    IndexConstReference index() const { return _a; }
    /** @brief Return the selected coefficient. */
    CoefficientConstReference coefficient() const { return _c; }
    /** @brief Construct a \c ExpansionConstReference from other. */
    ExpansionConstReference(const ExpansionValue<I, X>& other) : ExpansionConstReference(other.index(), other.coefficient()) {}
    /** @brief Construct a \c ExpansionConstReference from other. */
    ExpansionConstReference(const ExpansionReference<I, X>& other) : ExpansionConstReference(other.index(), other.coefficient()) {}
    /** @brief Convert to expansion value< i, x >. */
    operator ExpansionValue<I, X>() const { return ExpansionValue<I, X>(_a, _c); }
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const ExpansionConstReference<I, X>& ac1, const ExpansionConstReference<I, X>& ac2) {
        if constexpr (Convertible<EqualityType<X>, Bool>) {
            return ac1.index() == ac2.index() && ac1.coefficient() == ac2.coefficient();
        } else {
            return ac1.index() == ac2.index() && same(ac1.coefficient(), ac2.coefficient());
        }
    }
};

template<class I, class X>
ExpansionReference<I, X>& ExpansionReference<I, X>::operator=(const ExpansionValue<I, X>& other) {
    this->index() = other.index();
    this->coefficient() = other.coefficient();
    return *this;
}

template<class I, class X>
ExpansionReference<I, X>& ExpansionReference<I, X>::operator=(const ExpansionReference<I, X>& other) {
    this->index() = other.index();
    this->coefficient() = other.coefficient();
    return *this;
}

template<class I, class X>
ExpansionReference<I, X>& ExpansionReference<I, X>::operator=(const ExpansionConstReference<I, X>& other) {
    this->index() = other.index();
    this->coefficient() = other.coefficient();
    return *this;
}

/** @brief Test equality of the operands. */
template<class I, class X>
decltype(auto) operator==(const ExpansionConstReference<I, X>& m1, const ExpansionConstReference<I, X>& m2) {
    return m1.index() == m2.index() && (m1.coefficient() == m2.coefficient());
}

/** @brief Test inequality of the operands. */
template<class I, class X>
inline decltype(auto) operator!=(const ExpansionConstReference<I, X>& m1, const ExpansionConstReference<I, X>& m2) {
    return !(m1 == m2);
}

/** @brief Write the value to the output stream. */
template<class I, class X>
OutputStream& operator<<(OutputStream& os, const ExpansionReference<I, X>& m) {
    return os << ExpansionValue<I, X>(m);
}

/** @brief Write the value to the output stream. */
template<class I, class X>
OutputStream& operator<<(OutputStream& os, const ExpansionConstReference<I, X>& m) {
    return os << ExpansionValue<I, X>(m);
}

/** @brief Pointer proxy into expansion storage. */
template<class I, class X>
class ExpansionPointer {
    ExpansionReference<I, X> _r;

  public:
    /** @brief Construct a \c ExpansionPointer from ap, cp. */
    ExpansionPointer(typename UniformList<I>::Pointer ap, typename UniformList<X>::Pointer cp) : _r(*ap, *cp) {}
    /** @brief Construct a \c ExpansionPointer from p. */
    ExpansionPointer(ExpansionReference<I, X>* p) : _r(*p) {}
    /** @brief Return the value referenced by the iterator or pointer. */
    ExpansionReference<I, X>& operator*() { return _r; }
    /** @brief Return a pointer to the referenced object. */
    ExpansionReference<I, X>* operator->() { return &_r; }
};

/** @brief Pointer proxy into expansion const storage. */
template<class I, class X>
class ExpansionConstPointer {
    ExpansionConstReference<I, X> _r;

  public:
    /** @brief Construct a \c ExpansionConstPointer from ap, cp. */
    ExpansionConstPointer(typename UniformList<I>::ConstPointer ap, typename UniformList<X>::ConstPointer cp) : _r(*ap, *cp) {}
    /** @brief Construct a \c ExpansionConstPointer from p. */
    ExpansionConstPointer(ExpansionConstReference<I, X>* p) : _r(*p) {}
    /** @brief Return the value referenced by the iterator or pointer. */
    ExpansionConstReference<I, X>& operator*() { return _r; }
    /** @brief Return a pointer to the referenced object. */
    ExpansionConstReference<I, X>* operator->() { return &_r; }
};

/** @brief Iterator over expansion entries. */
template<class I, class X>
class ExpansionIterator
    : public IteratorFacade<ExpansionIterator<I, X>, ExpansionValue<I, X>, RandomAccessTraversalTag, ExpansionReference<I, X>> {
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename UniformList<I>::Iterator IndexIterator;
    typedef typename UniformList<X>::Iterator CoefficientIterator;
    friend class ExpansionConstIterator<I, X>;

  private:

  public:
    IndexIterator _ap;
    CoefficientIterator _cp;

  public:
    /** @brief Construct a \c ExpansionIterator from ap, cp. */
    ExpansionIterator(IndexIterator ap, CoefficientIterator cp) : _ap(ap), _cp(cp) {}
    /** @brief Test equality with the supplied iterator or value. */
    template<class XX>
    Bool equal(const ExpansionIterator<I, XX>& other) const { return this->_cp == other._cp; }
    /** @brief Test equality with the supplied iterator or value. */
    template<class XX>
    Bool equal(const ExpansionConstIterator<I, XX>& other) const { return this->_cp == other._cp; }
    /** @brief Move the iterator by the specified offset. */
    Void advance(PointerDifferenceType k) {
        _ap += k;
        _cp += k;
    }
    /** @brief Return the signed distance to the supplied iterator. */
    template<class XX>
    PointerDifferenceType distance_to(const ExpansionIterator<I, XX>& other) const { return other._cp - this->_cp; }
    /** @brief Return the value at the current iterator position. */
    ExpansionReference<I, X> dereference() const { return ExpansionReference<I, X>(*_ap, *_cp); }
    /** @brief Return a pointer to the referenced object. */
    ExpansionPointer<I, X> operator->() { return ExpansionPointer<I, X>(_ap.operator->(), _cp.operator->()); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ExpansionIterator<I, X>&) {
        return os << "{ExpansionIterator}";
    }
};

/** @brief Iterator over expansion const entries. */
template<class I, class X>
class ExpansionConstIterator
    : public IteratorFacade<ExpansionConstIterator<I, X>, ExpansionValue<I, X>, RandomAccessTraversalTag, ExpansionConstReference<I, X>> {
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename UniformList<I>::ConstIterator IndexConstIterator;
    typedef typename UniformList<X>::ConstIterator CoefficientConstIterator;
    friend class ExpansionIterator<I, X>;
    IndexConstIterator _ap;
    CoefficientConstIterator _cp;

  public:
    /** @brief Construct a \c ExpansionConstIterator from ap, cp. */
    ExpansionConstIterator(IndexConstIterator ap, CoefficientConstIterator cp) : _ap(ap), _cp(cp) {}
    /** @brief Construct a \c ExpansionConstIterator from other. */
    ExpansionConstIterator(const ExpansionIterator<I, X>& other) : ExpansionConstIterator(other._ap, other._cp) {}
    /** @brief Test equality with the supplied iterator or value. */
    template<class XX>
    Bool equal(const ExpansionIterator<I, XX>& other) const { return this->_cp == other._cp; }
    /** @brief Test equality with the supplied iterator or value. */
    template<class XX>
    Bool equal(const ExpansionConstIterator<I, XX>& other) const { return this->_cp == other._cp; }
    /** @brief Move the iterator by the specified offset. */
    Void advance(PointerDifferenceType k) {
        _ap += k;
        _cp += k;
    }
    /** @brief Return the signed distance to the supplied iterator. */
    template<class XX>
    PointerDifferenceType distance_to(const ExpansionConstIterator<I, XX>& other) const { return other._cp - this->_cp; }
    /** @brief Return the value at the current iterator position. */
    ExpansionConstReference<I, X> dereference() const { return ExpansionConstReference<I, X>(*_ap, *_cp); }
    /** @brief Return a pointer to the referenced object. */
    ExpansionConstPointer<I, X> operator->() { return ExpansionConstPointer<I, X>(_ap.operator->(), _cp.operator->()); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ExpansionConstIterator<I, X>&) {
        return os << "{ExpansionConstIterator}";
    }
};

/** @brief Reference proxy into expansion value storage. */
template<class I, class X>
class ExpansionValueReference {
    Expansion<I, X>& _e;
    const I& _a;

  public:
    /** @brief Construct a \c ExpansionValueReference from e, a. */
    ExpansionValueReference(Expansion<I, X>& e, const I& a) : _e(e), _a(a) {}
    /** @brief Convert to const x &. */
    operator const X&() const { return _e.get(_a); }
    // operator X& () { return _e.at(_a); }
    /** @brief Assign the supplied value and return this object. */
    ExpansionValueReference<I, X>& operator=(const X& x) {
        _e.at(_a) = x;
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    template<class Y>
        requires Assignable<X, Y>
    ExpansionValueReference<I, X>& operator=(const Y& y) {
        _e.at(_a) = y;
        return *this;
    }
    /** @brief Add the operand in place. */
    ExpansionValueReference<I, X>& operator+=(const X& c) {
        _e.at(_a) += c;
        return *this;
    }
    /** @brief Subtract the operand in place. */
    ExpansionValueReference<I, X>& operator-=(const X& c) {
        _e.at(_a) -= c;
        return *this;
    }
    /** @brief Multiply by the operand in place. */
    ExpansionValueReference<I, X>& operator*=(const X& c) {
        _e.at(_a) *= c;
        return *this;
    }
    /** @brief Divide by the operand in place. */
    ExpansionValueReference<I, X>& operator/=(const X& c) {
        _e.at(_a) /= c;
        return *this;
    }
};

template<class I, class X>
template<class Y, class... PRS>
inline Void Expansion<I, X>::_fill(const Expansion<I, Y>& other, PRS... prs) {
    MultiIndex a(other.argument_size());
    X x = this->zero_coefficient();
    for (auto iter = other.begin();
         iter != other.end(); ++iter) {
        a = iter->index();
        x = X(iter->coefficient(), prs...);
        if (decide(x != 0)) {
            this->append(a, x);
        }
    }
}

} // namespace Ariadne

#endif
