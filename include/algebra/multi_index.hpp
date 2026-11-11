/***************************************************************************
 *            algebra/multi_index.hpp
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

/*! \brief \file algebra/multi_index.hpp
 *  \brief An index specifying the degree of differentiation.
 */

#ifndef ARIADNE_MULTI_INDEX_HPP
#define ARIADNE_MULTI_INDEX_HPP

#include <cassert>
#include <initializer_list>
#include <iostream>

#include "numeric/numeric.hpp"
#include "utility/array.hpp"
#include "utility/container.hpp"
#include "utility/macros.hpp"
#include "utility/uniform_list.hpp"

namespace Ariadne {

/** @brief Single nonnegative exponent identifying a univariate monomial. */
class UniIndex {
    DegreeType _a;

  public:
    /** @brief Construct a \c UniIndex from a. */
    UniIndex(DegreeType a = 0u) : _a(a) {}
    /** @brief Increment the value or iterator. */
    UniIndex& operator++() {
        ++_a;
        return *this;
    }
    /** @brief Decrement the value or iterator. */
    UniIndex& operator--() {
        --_a;
        return *this;
    }
    /** @brief Return the sum or unary positive value. */
    UniIndex operator+(const UniIndex& other) const { return UniIndex(this->_a + other._a); }
    /** @brief Return the difference or additive inverse. */
    short int operator-(const UniIndex& other) const { return (short int)this->_a - (short int)other._a; }
    /** @brief Add the operand in place. */
    UniIndex& operator+=(const UniIndex& other) {
        this->_a += other._a;
        return *this;
    }
    /** @brief Return the sum or unary positive value. */
    template<BuiltinIntegral D>
    UniIndex operator+(D d) const { return UniIndex(static_cast<DegreeType>(this->_a + d)); }
    /** @brief Return the selected element or component. */
    DegreeType& operator[](IndexZero) { return _a; }
    /** @brief Return the selected element or component. */
    const DegreeType& operator[](IndexZero) const { return _a; }
    /** @brief Return the number of elements. */
    SizeOne size() const { return SizeOne(); }
    /** @brief Return the polynomial or expansion degree. */
    DegreeType degree() const { return _a; }
    /** @brief Convert to degree type. */
    operator DegreeType() const { return _a; }
};

/** @brief Storage for a multivariate exponent vector and its total degree. */
struct MultiIndexData {
    friend class MultiIndexList;

  public:
    typedef DegreeType IndexType;
    typedef DegreeType& Reference;
    typedef const DegreeType& ConstReference;

  protected:
    /** @brief Destroy the object and release its owned resources. */
    ~MultiIndexData();
    /** @brief Construct a \c MultiIndexData from n, p. */
    explicit MultiIndexData(SizeType n, IndexType* p);

  public:
    /** @brief Return the number of elements. */
    SizeType size() const;
    /** @brief Return the number of variables. */
    SizeType number_of_variables() const;
    /** @brief Return the polynomial or expansion degree. */
    DegreeType degree() const;
    /** @brief Return the selected element or component. */
    const IndexType& operator[](SizeType i) const;
    /** @brief Return the selected element or component. */
    IndexType& operator[](SizeType i);

    /** @brief Return the requested element or stored value. */
    DegreeType get(SizeType i) const;
    /** @brief Replace the selected value. */
    Void set(SizeType i, DegreeType n);

    /** @brief Assign the supplied values to the existing storage. */
    Void assign(const MultiIndexData& a);

  public:
    /** @brief Test equality of the operands. */
    friend Bool operator==(const MultiIndexData& a1, const MultiIndexData& a2);
    /** @brief Test inequality of the operands. */
    friend Bool operator!=(const MultiIndexData& a1, const MultiIndexData& a2);

  public:
    /** @brief Compare indices by total degree and then by the graded ordering. */
    friend Bool graded_less(const MultiIndexData& a1, const MultiIndexData& a2);
    /** @brief Compare indices in lexicographic order. */
    friend Bool lexicographic_less(const MultiIndexData& a1, const MultiIndexData& a2);
    /** @brief Compare indices in reverse lexicographic order. */
    friend Bool reverse_lexicographic_less(const MultiIndexData& a1, const MultiIndexData& a2);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream&, const MultiIndexData&);

  public:
    /** @brief Return an iterator to the first element. */
    IndexType* begin();
    /** @brief Return the past-the-end iterator. */
    IndexType* end();
    /** @brief Return an iterator to the first element. */
    const IndexType* begin() const;
    /** @brief Return the past-the-end iterator. */
    const IndexType* end() const;

  protected:

  public:
    SizeType _n;
    IndexType* _p;
};

class MultiIndex
    : public MultiIndexData {
  public:
    /** @brief Destroy the object and release its owned resources. */
    ~MultiIndex();
    /** @brief Construct a \c MultiIndex with the specified characteristics. */
    explicit MultiIndex();
    /** @brief Construct a \c MultiIndex from nv. */
    explicit MultiIndex(SizeType nv);
    /** @brief Construct a \c MultiIndex from nv, ary. */
    explicit MultiIndex(SizeType nv, const DegreeType* ary);
    /** @brief Construct a \c MultiIndex from lst. */
    MultiIndex(InitializerList<DegreeType> lst);

    /** @brief Copy-construct the object from the supplied instance. */
    MultiIndex(const MultiIndex& a);
    /** @brief Assign the supplied value and return this object. */
    MultiIndex& operator=(const MultiIndex& a);

    /** @brief Return the additive zero. */
    static MultiIndex zero(SizeType nv);
    /** @brief Construct the multiplicative unit. */
    static MultiIndex unit(SizeType nv, SizeType j);

    /** @brief Change the number of elements. */
    Void resize(SizeType n);
    /** @brief Remove all elements. */
    Void clear();

    /** @brief Increment the value or iterator. */
    MultiIndex& operator++();

    /** @brief Add the operand in place. */
    MultiIndex& operator+=(const MultiIndex& a);
    /** @brief Subtract the operand in place. */
    MultiIndex& operator-=(const MultiIndex& a);
    /** @brief Multiply by the operand in place. */
    MultiIndex& operator*=(const DegreeType& a);
    /** @brief Return the sum or unary positive value. */
    friend MultiIndex operator+(MultiIndex a1, const MultiIndex& a2);
    /** @brief Return the difference or additive inverse. */
    friend MultiIndex operator-(MultiIndex a1, const MultiIndex& a2);
    /** @brief Return the product. */
    friend MultiIndex operator*(MultiIndex a, IndexType s);
    /** @brief Return the product. */
    friend MultiIndex operator*(DegreeType s, MultiIndex a);

    /** @brief Exchange the contents of the operands. */
    friend Void swap(MultiIndex& a1, MultiIndex& a2);

    //    SizeType position() const;
    //    SizeType factorial() const;
    //    SizeType number() const;
};

/** @brief Ordering predicate using lexicographic order. */
struct LexicographicLess {
    /** @brief Apply the callable to the supplied arguments. */
    Bool operator()(const DegreeType& a1, const DegreeType& a2) const;
    /** @brief Apply the callable to the supplied arguments. */
    Bool operator()(const MultiIndex& a1, const MultiIndex& a2) const;
};

/** @brief Ordering predicate using reverse lexicographic order. */
struct ReverseLexicographicLess {
    /** @brief Apply the callable to the supplied arguments. */
    Bool operator()(const DegreeType& a1, const DegreeType& a2) const;
    /** @brief Apply the callable to the supplied arguments. */
    Bool operator()(const MultiIndex& a1, const MultiIndex& a2) const;
};

/** @brief Ordering predicate using graded order. */
struct GradedLess {
    /** @brief Apply the callable to the supplied arguments. */
    Bool operator()(const DegreeType& a1, const DegreeType& a2) const;
    /** @brief Apply the callable to the supplied arguments. */
    Bool operator()(const MultiIndex& a1, const MultiIndex& a2) const;
};

static const SizeType DEFAULT_CAPACITY = 4u;

struct MultiIndexReference : public MultiIndexData {
    /** @brief Construct a \c MultiIndexReference from n, p. */
    MultiIndexReference(SizeType n, IndexType* p);
    /** @brief Construct a \c MultiIndexReference with the specified characteristics. */
    MultiIndexReference(const MultiIndexData&);
    /** @brief Copy-construct the object from the supplied instance. */
    MultiIndexReference(const MultiIndexReference&) = default;
    /** @brief Assign the supplied value and return this object. */
    MultiIndexReference& operator=(const MultiIndexReference&);
    /** @brief Assign the supplied value and return this object. */
    MultiIndexReference& operator=(const MultiIndexData&);
    /** @brief Convert to multi index &. */
    operator MultiIndex&();
    /** @brief Convert to const multi index &. */
    operator const MultiIndex&() const;
    /** @brief Add the operand in place. */
    MultiIndexReference& operator+=(const MultiIndexData&);
    /** @brief Exchange the contents of the operands. */
    friend Void swap(MultiIndex&, MultiIndex&);
    //    friend Void swap(MultiIndexReference, MultiIndexReference);
};

/** @brief Reference proxy into multi index const storage. */
struct MultiIndexConstReference : public MultiIndexData {
    /** @brief Construct a \c MultiIndexConstReference from n, p. */
    MultiIndexConstReference(SizeType n, const DegreeType* p);
    /** @brief Copy-construct the object from the supplied instance. */
    MultiIndexConstReference(const MultiIndexConstReference&) = default;
    /** @brief Disallow assignment to preserve the ownership contract. */
    MultiIndexConstReference& operator=(const MultiIndexConstReference&) = delete;
    /** @brief Construct a \c MultiIndexConstReference with the specified characteristics. */
    MultiIndexConstReference(const MultiIndexReference&);
    /** @brief Convert to const multi index &. */
    operator const MultiIndex&() const;
    /** @brief Return the sum or unary positive value. */
    friend MultiIndex operator+(MultiIndex a1, const MultiIndex& a2);
};

/** @brief Pointer proxy into multi index storage. */
struct MultiIndexPointer {
    MultiIndexReference _r;

  public:
    /** @brief Construct a \c MultiIndexPointer from _p. */
    MultiIndexPointer(MultiIndexReference* _p) : _r(*_p) {}
    /** @brief Return the value referenced by the iterator or pointer. */
    MultiIndexReference& operator*() { return _r; }
    /** @brief Return a pointer to the referenced object. */
    MultiIndexReference* operator->() { return &_r; }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const MultiIndexPointer&);
};
/** @brief Pointer proxy into multi index const storage. */
struct MultiIndexConstPointer {
    MultiIndexConstReference _r;

  public:
    /** @brief Construct a \c MultiIndexConstPointer from _p. */
    MultiIndexConstPointer(MultiIndexConstReference* _p) : _r(*_p) {}
    /** @brief Return the value referenced by the iterator or pointer. */
    MultiIndexConstReference& operator*() const { return const_cast<MultiIndexConstReference&>(_r); }
    /** @brief Return a pointer to the referenced object. */
    MultiIndexConstReference* operator->() const { return const_cast<MultiIndexConstReference*>(&_r); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const MultiIndexConstPointer&);
};

class MultiIndexListIterator;
class MultiIndexListConstIterator;

class MultiIndexList {
    SizeType _capacity;
    SizeType _size;
    SizeType _argument_size;
    DegreeType* _indices;

  public:
    typedef MultiIndex ValueType;
    typedef MultiIndexReference Reference;
    typedef MultiIndexConstReference ConstReference;
    typedef MultiIndexPointer Pointer;
    typedef MultiIndexConstPointer ConstPointer;
    typedef MultiIndexListIterator Iterator;
    typedef MultiIndexListConstIterator ConstIterator;

  protected:

  public:
    /** @brief Construct a \c MultiIndexList from as. */
    explicit MultiIndexList(SizeType as);

  public:
    /** @brief Destroy the object and release its owned resources. */
    ~MultiIndexList();
    //    MultiIndexList(SizeType as, SizeType cap);
    /** @brief Construct a \c MultiIndexList from lst. */
    MultiIndexList(const InitializerList<InitializerList<DegreeType>>& lst);
    /** @brief Construct a \c MultiIndexList from lst. */
    MultiIndexList(const InitializerList<MultiIndex>& lst);
    /** @brief Construct a \c MultiIndexList from n, a. */
    MultiIndexList(SizeType n, const MultiIndex& a);
    /** @brief Copy-construct the object from the supplied instance. */
    MultiIndexList(const MultiIndexList&);
    /** @brief Move-construct the object from the supplied instance. */
    MultiIndexList(MultiIndexList&&);
    /** @brief Assign the supplied value and return this object. */
    MultiIndexList& operator=(const MultiIndexList&);
    /** @brief Assign the supplied value and return this object. */
    MultiIndexList& operator=(MultiIndexList&&);
    /** @brief Change the number of elements. */
    Void resize(SizeType n);
    /** @brief Reserve storage for at least the requested number of elements. */
    Void reserve(SizeType n);
    /** @brief Return the number of elements. */
    SizeType size() const;
    /** @brief Return the allocated element capacity. */
    SizeType capacity() const;
    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const;
    /** @brief Append the supplied element or sequence. */
    Void append(const MultiIndexData& a);
    /** @brief Append the sum of the supplied expansion terms. */
    Void append_sum(const MultiIndexData& a1, const MultiIndexData& a2);
    /** @brief Return the selected element or component. */
    Reference operator[](SizeType i);
    /** @brief Return the selected element or component. */
    ConstReference operator[](SizeType i) const;
    /** @brief Return the first element. */
    Reference front();
    /** @brief Return the first element. */
    ConstReference front() const;
    /** @brief Return the last element. */
    Reference back();
    /** @brief Return the last element. */
    ConstReference back() const;
    /** @brief Return an iterator to the first element. */
    Iterator begin();
    /** @brief Return the past-the-end iterator. */
    Iterator end();
    /** @brief Return an iterator to the first element. */
    ConstIterator begin() const;
    /** @brief Return the past-the-end iterator. */
    ConstIterator end() const;
    /** @brief Erase the selected entry or range. */
    Iterator erase(Iterator pos);
    /** @brief Remove all elements. */
    Void clear();
    /** @brief Test equality of the operands. */
    friend Bool operator==(const MultiIndexList& lst1, const MultiIndexList& lst2);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const MultiIndexList& lst);
};

class MultiIndexListIterator {
    friend class MultiIndexListConstIterator;
    MultiIndexReference _r;

  public:
    /** @brief Construct a \c MultiIndexListIterator from n, p. */
    explicit MultiIndexListIterator(SizeType n, DegreeType* p);
    /** @brief Copy-construct the object from the supplied instance. */
    MultiIndexListIterator(const MultiIndexListIterator&);
    /** @brief Assign the supplied value and return this object. */
    MultiIndexListIterator& operator=(const MultiIndexListIterator&);
    /** @brief Test equality of the operands. */
    Bool operator==(const MultiIndexListIterator& other) const;
    /** @brief Test inequality of the operands. */
    Bool operator!=(const MultiIndexListIterator& other) const;
    /** @brief Increment the value or iterator. */
    MultiIndexListIterator& operator++();
    /** @brief Decrement the value or iterator. */
    MultiIndexListIterator& operator--();
    /** @brief Add the operand in place. */
    MultiIndexListIterator& operator+=(PointerDifferenceType k);
    /** @brief Return the sum or unary positive value. */
    MultiIndexListIterator operator+(PointerDifferenceType k) const;
    /** @brief Return the value referenced by the iterator or pointer. */
    MultiIndexReference operator*() const;
    /** @brief Return a pointer to the referenced object. */
    MultiIndexPointer operator->() const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const MultiIndexListIterator&);
};

class MultiIndexListConstIterator {
    MultiIndexConstReference _r;

  public:
    /** @brief Construct a \c MultiIndexListConstIterator from n, p. */
    explicit MultiIndexListConstIterator(SizeType n, const DegreeType* p);
    /** @brief Construct a \c MultiIndexListConstIterator from iter. */
    MultiIndexListConstIterator(const MultiIndexListIterator& iter);
    /** @brief Copy-construct the object from the supplied instance. */
    MultiIndexListConstIterator(const MultiIndexListConstIterator& iter);
    /** @brief Assign the supplied value and return this object. */
    MultiIndexListConstIterator& operator=(const MultiIndexListConstIterator& iter);
    /** @brief Test equality of the operands. */
    Bool operator==(const MultiIndexListConstIterator& other);
    /** @brief Test inequality of the operands. */
    Bool operator!=(const MultiIndexListConstIterator& other) const;
    /** @brief Increment the value or iterator. */
    MultiIndexListConstIterator& operator++();
    /** @brief Decrement the value or iterator. */
    MultiIndexListConstIterator& operator--();
    /** @brief Add the operand in place. */
    MultiIndexListConstIterator& operator+=(PointerDifferenceType k);
    /** @brief Return the sum or unary positive value. */
    MultiIndexListConstIterator operator+(PointerDifferenceType k) const;
    /** @brief Return the value referenced by the iterator or pointer. */
    MultiIndexConstReference operator*() const;
    /** @brief Return a pointer to the referenced object. */
    MultiIndexConstPointer operator->() const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const MultiIndexListConstIterator&);
};

/** @brief Sequence whose elements share construction characteristics. */
template<>
class UniformList<MultiIndex>
    : public MultiIndexList {
    /** @brief Move-construct the object from the supplied instance. */
    using MultiIndexList::MultiIndexList;
};

} // namespace Ariadne

#include "multi_index.inl.hpp"

#endif /* ARIADNE_MULTI_INDEX_HPP */
