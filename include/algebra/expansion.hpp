/***************************************************************************
 *            algebra/expansion.hpp
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

/** @file algebra/expansion.hpp
 *  @brief Sparse expansions in multivariate monomials.
 */

#ifndef ARIADNE_EXPANSION_HPP
#define ARIADNE_EXPANSION_HPP

#include "algebra/linear_algebra.decl.hpp"
#include "multi_index.hpp"
#include "utility/iterator.hpp"
#include "utility/macros.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

/************ Expansion ******************************************************/

/** @brief Exchange the contents of the operands. */
inline Void swap(UniIndex& a1, UniIndex& a2) { std::swap(a1, a2); }
// template<class T> inline Void swap(T& t1, T& t2) { std::swap(t1,t2); }

/** @brief Return the degree of the supplied object. */
inline DegreeType degree_of(const MultiIndex& a) { return a.degree(); }
/** @brief Return the degree of the supplied object. */
inline DegreeType degree_of(const UniIndex& a) { return a; }

/** @brief Return the dimension or element count of the supplied object. */
inline SizeType size_of(const MultiIndex& a) { return a.size(); }
/** @brief Return the dimension or element count of the supplied object. */
inline SizeOne size_of(const UniIndex&) { return SizeOne(); }

/** @brief Return the argument dimension of the supplied object. */
inline SizeType argument_size_of(const UniformList<MultiIndex>& as) { return as.argument_size(); }
/** @brief Return the argument dimension of the supplied object. */
inline SizeOne argument_size_of(const UniformList<UniIndex>&) { return SizeOne(); }

template<class T>
using UniformReference = typename UniformList<T>::Reference;
template<class T>
using UniformConstReference = typename UniformList<T>::ConstReference;
template<class T>
using UniformPointer = typename UniformList<T>::Pointer;
template<class T>
using UniformConstPointer = typename UniformList<T>::ConstPointer;

template<class I, class X>
class ExpansionValue;
template<class I, class X>
class ExpansionReference;
template<class I, class X>
class ExpansionConstReference;
template<class I, class X>
class ExpansionPointer;
template<class I, class X>
class ExpansionConstPointer;
template<class I, class X>
class ExpansionIterator;
template<class I, class X>
class ExpansionConstIterator;
template<class I, class X>
class ExpansionValueReference;

template<class CMP>
struct IndexComparison;
struct GradedIndexLess;
struct LexicographicIndexLess;
struct ReverseLexicographicIndexLess;
struct CoefficientLess;
struct CoefficientIsZero;

/** @brief Compile-time traits for index. */
template<class I>
struct IndexTraits;

/** @brief Compile-time traits for index. */
template<>
struct IndexTraits<UniIndex> {
    typedef DegreeType InitializerType;
    typedef SizeOne SizeOfType;
    typedef IndexZero IndexIntoType;
    typedef String NameType;
    template<class Y>
    using Argument = Scalar<Y>;
};

/** @brief Compile-time traits for index. */
template<>
struct IndexTraits<MultiIndex> {
    typedef InitializerList<DegreeType> InitializerType;
    typedef SizeType SizeOfType;
    typedef SizeType IndexIntoType;
    typedef Array<String> NameType;
    template<class Y>
    using Argument = Vector<Y>;
};

template<class I, class X>
using ArgumentOfType = typename IndexTraits<I>::template Argumentq<X>;

template<class I, class X>
class Expansion;

//! \brief A sparse list of index-value pairs.
template<class I, class X>
class Expansion {
    using S = typename IndexTraits<I>::SizeOfType;
    using V = typename IndexTraits<I>::IndexIntoType;

  public:
    using IndexInitializerType = typename IndexTraits<I>::InitializerType;
    template<class Y>
    using Argument = typename IndexTraits<I>::template Argument<Y>;

    //    static const SizeType DEFAULT_CAPACITY=16;

  public:
    UniformList<I> _indices;
    UniformList<X> _coefficients;

  public:
    typedef S ArgumentSizeType;
    typedef V VariableIndexType;
    typedef I IndexType;
    typedef X CoefficientType;

    typedef ExpansionValue<I, X> ValueType;
    typedef ExpansionReference<I, X> Reference;
    typedef ExpansionConstReference<I, X> ConstReference;
    typedef ExpansionPointer<I, X> Pointer;
    typedef ExpansionConstPointer<I, X> ConstPointer;
    typedef ExpansionIterator<I, X> Iterator;
    typedef ExpansionConstIterator<I, X> ConstIterator;

    typedef typename UniformList<I>::Reference IndexReference;
    typedef typename UniformList<X>::Reference CoefficientReference;
    typedef typename UniformList<I>::ConstReference IndexConstReference;
    typedef typename UniformList<X>::ConstReference CoefficientConstReference;

    /** @brief Destroy the object and release its owned resources. */
    ~Expansion();
    /** @brief Construct a \c Expansion from as, prs. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    explicit Expansion(ArgumentSizeType as, PRS... prs) : Expansion(as, get_characteristics(X(0u, prs...))) {} // DEPRECATED
    /** @brief Construct a \c Expansion from as, prs, cap. */
    explicit Expansion(ArgumentSizeType as, const CharacteristicsType<X>& prs, SizeType cap = DEFAULT_CAPACITY);
    /** @brief Construct a \c Expansion from lst. */
    Expansion(InitializerList<Pair<IndexInitializerType, X>> lst);
    /** @brief Construct a \c Expansion from as, prs, cap. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    explicit Expansion(ArgumentSizeType as, PRS... prs, SizeType cap = DEFAULT_CAPACITY);
    /** @brief Construct a \c Expansion from lst, prs. */
    template<class... PRS>
        requires Constructible<X, ExactDouble, PRS...>
    Expansion(InitializerList<Pair<IndexInitializerType, ExactDouble>> lst, PRS... prs);
    /** @brief Construct a \c Expansion from prs. */
    template<class Y, class... PRS>
        requires Constructible<X, Y, PRS...>
    explicit Expansion(const Expansion<I, Y>&, PRS... prs);
    /** @brief Copy-construct the object from the supplied instance. */
    Expansion(const Expansion<I, X>&);
    /** @brief Assign the supplied value and return this object. */
    Expansion<I, X>& operator=(const Expansion<I, X>&);
    /** @brief Move-construct the object from the supplied instance. */
    Expansion(Expansion<I, X>&&);
    /** @brief Assign the supplied value and return this object. */
    Expansion<I, X>& operator=(Expansion<I, X>&&);
    /** @brief Exchange the contents of the operands. */
    Void swap(Expansion<I, X>&);

    /** @brief Test whether the operands have the same representation. */
    Bool same_as(const Expansion<I, X>& e) const;

    /** @brief Test whether the container is empty. */
    Bool empty() const;

    /** @brief Return the characteristics required to construct compatible values. */
    Pair<ArgumentSizeType, CharacteristicsType<X>> characteristics() const;
    /** @brief Return the common characteristics of the coefficients. */
    CharacteristicsType<X> coefficient_characteristics() const;
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const;
    /** @brief Return the number of terms. */
    SizeType number_of_terms() const;
    /** @brief Return the number of nonzero stored terms. */
    SizeType number_of_nonzeros() const; // DEPRECATED
    /** @brief Return the number of elements. */
    SizeType size() const;
    /** @brief Return the allocated element capacity. */
    SizeType capacity() const;
    /** @brief Change the number of elements. */
    Void resize(SizeType sz);
    /** @brief Reserve storage for at least the requested number of elements. */
    Void reserve(SizeType cap);

    /** @brief Return a zero with the coefficient characteristics. */
    const CoefficientType& zero_coefficient() const;
    /** @brief Return the selected element or component. */
    const CoefficientType& operator[](const IndexType& a) const;
    /** @brief Return the selected element or component. */
    ExpansionValueReference<I, X> operator[](const IndexType& a);

    /** @brief Replace the selected value. */
    Void set(const IndexType& a, const CoefficientType& c);
    /** @brief Return the element at the specified index. */
    CoefficientType& at(const IndexType& a);
    /** @brief Return the requested element or stored value. */
    const CoefficientType& get(const IndexType& a) const;

    /** @brief Return the first element. */
    Reference front();
    /** @brief Return the last element. */
    Reference back();
    /** @brief Return the first element. */
    ConstReference front() const;
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

    /** @brief Insert the supplied element or sequence at the beginning. */
    Void prepend(const IndexType& a, const CoefficientType& x);
    /** @brief Append the supplied element or sequence. */
    Void append(const IndexType& a, const CoefficientType& x);
    /** @brief Append the sum of the supplied expansion terms. */
    Void append_sum(const IndexType& a1, const IndexType& a2, const CoefficientType& x);

    /** @brief Append the supplied element or sequence. */
    template<class Y>
        requires Assignable<X, Y>
    Void append(const IndexType& a, const Y& y) {
        X x = this->zero_coefficient();
        x = y;
        this->append(a, x);
    }

    /** @brief Find the entry matching the supplied key. */
    Iterator find(const IndexType& a);
    /** @brief Find the entry matching the supplied key. */
    ConstIterator find(const IndexType& a) const;

    /** @brief Insert the supplied entry. */
    Iterator insert(Iterator pos, const IndexType& a, const CoefficientType& x);
    /** @brief Erase the selected entry or range. */
    Iterator erase(Iterator pos);

    /** @brief Remove all elements. */
    Void clear();
    /** @brief Remove entries with zero coefficients. */
    Void remove_zeros();
    /** @brief Combine entries with equal indices. */
    Void combine_terms();
    /** @brief Check the predicate at the supplied effort or validate the supplied data. */
    Void check() const;

    // template<class CMP> Void sort(CMP cmp);
    /** @brief Sort expansion terms by their indices. */
    Void index_sort(ReverseLexicographicLess cmp);
    /** @brief Sort expansion terms by their indices. */
    Void index_sort(GradedLess cmp);
    /** @brief Sort the stored entries. */
    Void sort(ReverseLexicographicIndexLess cmp);
    /** @brief Sort the stored entries. */
    Void sort(GradedIndexLess cmp);

    /** @brief Sort indices in reverse lexicographic order. */
    Void reverse_lexicographic_sort();
    /** @brief Sort expansion terms by graded index order. */
    Void graded_sort();

    /** @brief Test whether the entries are sorted. */
    Bool is_sorted(ReverseLexicographicIndexLess cmp);
    /** @brief Test whether the entries are sorted. */
    Bool is_sorted(GradedIndexLess cmp);

    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os, const typename IndexTraits<I>::NameType& vars) const;

  public:
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Expansion<I, X>& self) { return self._write(os); }
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const Expansion<I, X>& e1, const Expansion<I, X>& e2) { return e1.same_as(e2); }
    /** @brief Embed in the larger argument space. */
    friend Expansion<MultiIndex, CoefficientType> embed(SizeType as1, const Expansion<IndexType, CoefficientType>& e2, SizeType as3) { return _embed(as1, e2, as3); }

  private:
    /** @brief Embed in the larger argument space. */
    static Expansion<MultiIndex, CoefficientType> _embed(SizeType as1, const Expansion<IndexType, CoefficientType>& e2, SizeType as3);
    /** @brief Fill the stored entries from the supplied source. */
    template<class Y, class... PRS>
    Void _fill(const Expansion<I, Y>& other, PRS... prs);
};

/** @brief Sparse index-coefficient expansion maintained in sorted order. */
template<class I, class X, class CMP>
class SortedExpansion : public Expansion<I, X> {
  public:
    typedef I IndexType;
    typedef X CoefficientType;
    typedef typename Expansion<I, X>::Iterator Iterator;
    typedef typename Expansion<I, X>::ConstIterator ConstIterator;

  public:
    using Expansion<I, X>::Expansion;
    /** @brief Construct a \c SortedExpansion from e. */
    SortedExpansion(Expansion<I, X> e);
    /** @brief Sort the stored entries. */
    Void sort();
    /** @brief Insert the supplied entry. */
    Void insert(const IndexType& a, const CoefficientType& c);
    /** @brief Replace the selected value. */
    Void set(const IndexType& a, const CoefficientType& c);
    /** @brief Return the element at the specified index. */
    CoefficientType& at(const IndexType& a);
    /** @brief Return the requested element or stored value. */
    const CoefficientType& get(const IndexType& a) const;
    //    Iterator find(const IndexType& a);
    //    ConstIterator find(const IndexType& a) const;
    /** @brief Check the predicate at the supplied effort or validate the supplied data. */
    Void check() const; // Check the expansion is sorted and has unique terms

    /** @brief Test equality of the operands. */
    EqualityType<X> operator==(const SortedExpansion<I, X, CMP>& other) const;
    /** @brief Test inequality of the operands. */
    EqualityType<X> operator!=(const SortedExpansion<I, X, CMP>& other) const;
};

template<class I, class X>
template<class... PRS>
    requires Constructible<X, Nat, PRS...>
Expansion<I, X>::Expansion(ArgumentSizeType as, PRS... prs, SizeType cap)
    : Expansion(as, X(0u, prs...), cap) {
}

template<class I, class X>
template<class... PRS>
    requires Constructible<X, ExactDouble, PRS...>
Expansion<I, X>::Expansion(InitializerList<Pair<IndexInitializerType, ExactDouble>> lst, PRS... prs)
    : Expansion(ArgumentSizeType(), make_characteristics<X>(prs...)) {
    ARIADNE_PRECONDITION(lst.size() != 0);

    _indices = UniformList<I>(0u, I(lst.begin()->first) * 0);
    _coefficients = UniformList<X>(prs...);

    SizeType cap = std::max(DEFAULT_CAPACITY, lst.size());
    _indices.reserve(cap);
    _coefficients.reserve(cap);

    for (auto iter = lst.begin();
         iter != lst.end(); ++iter) {
        I a = iter->first;
        X x(iter->second, prs...);
        if (decide(x != 0)) {
            this->append(a, x);
        }
    }
}

template<class I, class X>
template<class Y, class... PRS>
    requires Constructible<X, Y, PRS...>
Expansion<I, X>::Expansion(const Expansion<I, Y>& other, PRS... prs)
    : Expansion(other.argument_size(), get_characteristics(X(other.zero_coefficient(), prs...))) {
    this->_fill(other, prs...);
}

} // namespace Ariadne

#include "expansion.inl.hpp"

#endif
