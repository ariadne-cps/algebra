/***************************************************************************
 *            test_expansion.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <vector>
#include <utility>
#include <initializer_list>
#include <sstream>
#include "numeric/numeric.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "interval/interval.hpp"
#include "algebra/expansion.tpl.hpp"

#include "algebra_testing_context.hpp"

using namespace std;
using namespace Ariadne;

namespace Ariadne {
DP characteristics(FloatDP flt);
MultiplePrecision characteristics(FloatMP flt);
}

template<class K, class V> struct MapValue {
    typedef K key_type;
    K first; V second;
    MapValue(const K& k, const V& v) : first(k), second(v) { }
    const K& index() const { return first; }
    const V& coefficient() const { return second; }
};


template<class F> class TestExpansion
{
    typedef PrecisionType<F> PR;
    typedef MultiIndex MI;
    typedef Expansion<MI,F> ExpansionType;

    typedef typename Expansion<MI,F>::ValueType ExpansionValueType;
    typedef typename Expansion<MI,F>::Reference ExpansionReference;
    typedef typename Expansion<MI,F>::ConstReference ExpansionConstReference;
    typedef typename Expansion<MI,F>::Pointer ExpansionExpansionPointer;
    typedef typename Expansion<MI,F>::ConstPointer ExpansionExpansionConstPointer;
    typedef typename Expansion<MI,F>::Iterator ExpansionIterator;
    typedef typename Expansion<MI,F>::ConstIterator ExpansionConstIterator;

    typedef typename Expansion<MI,F>::IndexReference ExpansionIndexReference;
    typedef typename Expansion<MI,F>::IndexConstReference ExpansionIndexConstReference;
    typedef typename Expansion<MI,F>::CoefficientReference ExpansionCoefficientReference;
    typedef typename Expansion<MI,F>::CoefficientConstReference ExpansionCoefficientConstReference;
  public:
    PR prec; F zero;
    GradedLess graded_less;
    LexicographicLess lexicographic_less;
    ReverseLexicographicLess reverse_lexicographic_less;
  public:
    TestExpansion(F const& z);
    Void test();
  private:
    Void test_working();
    Void test_helpers();
    Void test_data_access();
    Void test_equality();
    Void test_sort();
    Void test_cleanup();
    Void test_constructors();
    //Void test_indexing();
    Void test_find();
    Void test_embed();
};


template<class F> TestExpansion<F>::TestExpansion(F const& z)
    : prec(z.precision()), zero(z)
{
}

template<class F> Void TestExpansion<F>::test()
{
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_working());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_helpers());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_data_access());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_equality());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_sort());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_cleanup());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_constructors());
    //ARIADNE_TEST_CALL(test_indexing());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_find());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_embed());
}


template<class F> Void TestExpansion<F>::test_working()
{
    MultiIndex helper_multi({1u,2u});
    UniIndex helper_uni(3u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), degree_of(helper_multi), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), degree_of(helper_uni), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), size_of(helper_multi), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), static_cast<SizeType>(size_of(helper_uni)), 1u);

    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "zero", zero);
    MultiIndexList as(3u);
    ExpansionType e(3u,prec);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.argument_size(), 3u);
    // Append values
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({0,0,0},2.0_x));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.begin()->index().number_of_variables(), 3);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({1,0,0},3.0_x));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({0,1,0},5.0_x));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({1,0,1},7.0_x));
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.find({0,0,0}), e.begin());
    //assert(false);

    // Regression test to ensure that indices and coefficients always have same capacity
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), e._indices.capacity()==e._coefficients.capacity());
    ExpansionType ce=e;
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), ce._indices.capacity()==ce._coefficients.capacity());

}


template<class F> Void TestExpansion<F>::test_helpers()
{
    MultiIndex zero_index({0u,0u});
    MultiIndex first_index({1u,0u});
    MultiIndex second_index({0u,1u});
    F coefficient_one(1,prec);
    F coefficient_two(2,prec);
    F coefficient_three(3,prec);

    ExpansionType checked(2u,prec);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), checked.check());

    ExpansionType self_copy(2u,prec);
    self_copy.append(zero_index,coefficient_one);
    ExpansionType& self_copy_alias=self_copy;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), self_copy=self_copy_alias);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), self_copy.size(), 1u);

    ExpansionType self_move(2u,prec);
    self_move.append(first_index,coefficient_two);
    ExpansionType* self_move_alias=&self_move;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), self_move=std::move(*self_move_alias));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), self_move.size(), 1u);
    auto checked_characteristics=checked.characteristics();
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), checked_characteristics.first, 2u);

    ExpansionType broken(2u,prec);
    broken._indices.append(zero_index);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), broken.check());

    ExpansionType left(2u,prec);
    ExpansionType right(2u,prec);
    left.append(zero_index,coefficient_one);
    right.append(first_index,coefficient_two);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), left.swap(right));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), left.front().index(), first_index);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), left.front().coefficient(), coefficient_two);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), left.back().index(), first_index);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), left.back().coefficient(), coefficient_two);
    ExpansionType const& const_left=left;
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), const_left.front().index(), first_index);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), const_left.back().coefficient(), coefficient_two);

    ExpansionType resized(2u,prec);
    SizeType initial_capacity=resized.capacity();
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), resized.reserve(initial_capacity));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), resized.reserve(initial_capacity+1u));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), resized.capacity()>=initial_capacity+1u);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), resized.resize(2u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), resized.size(), 2u);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), resized.resize(1u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), resized.size(), 1u);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), resized.resize(3u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), resized.size(), 3u);

    ExpansionType different_argument_size(3u,prec);
    different_argument_size.append(MultiIndex({0u,0u,0u}),coefficient_one);
    ExpansionType one_term_two_variables(2u,prec);
    one_term_two_variables.append(zero_index,coefficient_one);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !same(one_term_two_variables,different_argument_size));

    ExpansionType values(2u,prec);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), values.set(zero_index,coefficient_one));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), values.set(zero_index,coefficient_two));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), values.get(zero_index), coefficient_two);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), values.get(first_index), zero);

    ExpansionType spare(2u,prec);
    spare.reserve(4u);
    spare.append(zero_index,coefficient_one);
    spare.append(second_index,coefficient_three);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), spare.append_sum(zero_index,first_index,coefficient_two));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), spare.back().index(), first_index);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), spare.back().coefficient(), coefficient_two);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), spare.erase(spare.end()-1));
    auto inserted=spare.insert(spare.begin()+1,first_index,coefficient_two);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), inserted->index(), first_index);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), spare.erase(inserted));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), spare.size(), 2u);

    ExpansionType tight(2u,characteristics(coefficient_one),1u);
    tight.append(zero_index,coefficient_one);
    auto tight_inserted=tight.insert(tight.begin(),first_index,coefficient_two);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), tight_inserted->index(), first_index);

    ExpansionType tight_sum(2u,characteristics(coefficient_one),1u);
    tight_sum.append(zero_index,coefficient_one);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), tight_sum.append_sum(zero_index,first_index,coefficient_two));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), tight_sum.size(), 2u);

    ExpansionType write_cases(2u,prec);
    write_cases.append(MultiIndex({1u,1u}),coefficient_one);
    write_cases.append(zero_index,coefficient_one);
    std::ostringstream write_stream;
    Array<String> write_names({"x","y"});
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), write_cases._write(write_stream,write_names));

    ExpansionType reverse_sortable(2u,prec);
    reverse_sortable.append(zero_index,coefficient_one);
    reverse_sortable.append(first_index,coefficient_two);
    reverse_sortable.append(second_index,coefficient_three);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), reverse_sortable.sort(ReverseLexicographicIndexLess()));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), reverse_sortable.is_sorted(ReverseLexicographicIndexLess()));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), reverse_sortable.reverse_lexicographic_sort());
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), reverse_sortable.is_sorted(ReverseLexicographicIndexLess()));

    if constexpr (Same<F,RoundedFloatDP>) {
        SortedExpansion<MI,F,ReverseLexicographicIndexLess> bad_reverse(2u,prec);
        bad_reverse.append(zero_index,coefficient_one);
        bad_reverse.append(first_index,coefficient_two);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), bad_reverse.check());

        ExpansionType sortable(2u,prec);
        sortable.append(first_index,coefficient_two);
        sortable.append(zero_index,coefficient_one);
        SortedExpansion<MI,F,GradedIndexLess> sorted(std::move(sortable));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), sorted.check());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), sorted.insert(second_index,coefficient_three));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), sorted.insert(second_index,coefficient_three));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), sorted.set(first_index,coefficient_three));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), sorted.get(first_index), coefficient_three);
        MultiIndex missing_index({2u,0u});
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), sorted.get(missing_index), zero);

        SortedExpansion<MI,F,GradedIndexLess> bad_sorted(2u,prec);
        bad_sorted.append(first_index,coefficient_two);
        bad_sorted.append(zero_index,coefficient_one);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), bad_sorted.check());

        ExpansionType first_only_base(2u,prec);
        first_only_base.append(first_index,zero);
        SortedExpansion<MI,F,GradedIndexLess> first_only(std::move(first_only_base));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), first_only.get(zero_index), zero);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), first_only.insert(zero_index,coefficient_one));

        ExpansionType zero_only_base(2u,prec);
        zero_only_base.append(zero_index,zero);
        SortedExpansion<MI,F,GradedIndexLess> zero_only(std::move(zero_only_base));

        ExpansionType first_zero_base(2u,prec);
        first_zero_base.append(first_index,zero);
        SortedExpansion<MI,F,GradedIndexLess> first_zero(std::move(first_zero_base));

        ExpansionType zero_first_base(2u,prec);
        zero_first_base.append(zero_index,zero);
        zero_first_base.append(first_index,zero);
        SortedExpansion<MI,F,GradedIndexLess> zero_first(std::move(zero_first_base));

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(zero_only==first_zero));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(first_zero==zero_only));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(zero_first==zero_only));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(zero_only==zero_first));

        Expansion<UniIndex,F> univariate_embed(SizeOne(),prec);
        univariate_embed.append(UniIndex(0u),coefficient_one);
        univariate_embed.append(UniIndex(2u),coefficient_two);
        auto embedded=embed(1u,univariate_embed,1u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), embedded.argument_size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), embedded[MultiIndex({0u,0u,0u})], coefficient_one);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), embedded[MultiIndex({0u,2u,0u})], coefficient_two);
    }

    std::ostringstream os;
    Array<String> names({"x","y"});
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), values._write(os,names));
    Array<String> wrong_names({"x"});
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)values._write(os,wrong_names));
}


namespace {

template<class F>
constexpr bool check_concept()
{
    using E=Expansion<MultiIndex,F>;
    return requires(F x, PrecisionType<F> prec, SizeType as, E e, E const ce, MultiIndex a) {
        E(as,prec);
        e=E(ce);
        e.reserve(2u);
        e.set(a,x);
        e.prepend(a,x);
        e.append(a,x);
        e.append_sum(a,a,x);
        e.clear();
        e.index_sort(GradedLess());
        e.sort(GradedIndexLess());
        e.sort(ReverseLexicographicIndexLess());
        x=ce[a];
        ce.number_of_terms();
        ce.argument_size();
        e.erase(e.begin());
        ce.check();
    };
}

template<class F>
constexpr bool check_iterator_concept()
{
    using E=Expansion<MultiIndex,F>;
    using Iterator=typename E::Iterator;
    using ConstIterator=typename E::ConstIterator;
    using Reference=typename E::Reference;
    using ConstReference=typename E::ConstReference;
    return requires(E e, E const ce, MultiIndex a, Iterator iter, ConstIterator citer,
                    Reference ref, ConstReference cref) {
        iter=e.begin(); iter=e.end(); iter=e.find(a);
        citer=e.begin(); citer=e.end(); citer=e.find(a);
        citer=ce.begin(); citer=ce.end(); citer=ce.find(a);
        *iter; *citer;
        ++iter; --iter; ++citer; --citer;
        iter==iter; iter!=iter; citer==citer; citer!=citer;
        citer==iter; citer!=iter; iter==citer; iter!=citer;
        ref=cref;
    };
}

static_assert(check_concept<RoundedFloatDP>());
static_assert(check_concept<Rounded<FloatMP>>());
static_assert(check_iterator_concept<RoundedFloatDP>());
static_assert(check_iterator_concept<Rounded<FloatMP>>());

} // namespace


// Test dereferencing of iterators
template<class F> Void TestExpansion<F>::test_data_access()
{
    ExpansionType e(3,prec);
    ExpansionType const& ce=e;
    MultiIndex a(3);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e[a]=11);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e[a], 11);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), ce[a], 11);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.at(a), 11);
    e.clear();

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.at(a)=17);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e[a], 17);
    e.clear();

    SortedExpansion<MultiIndex,F,GradedIndexLess> se(3,prec);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), se[a]=11);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), se[a], 11);
    se.clear();


    // Append values
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({0,0,0},F(2.0_x,prec)));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({1,0,0},F(3.0_x,prec)));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({0,1,0},F(5.0_x,prec)));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.append({1,0,1},F(7.0_x,prec)));

    // Test ExpansionIterator difference
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e.begin()", e.begin());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e.end()", e.end());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.begin()+static_cast<PointerDifferenceType>(e.number_of_terms()), e.end());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.end()-static_cast<PointerDifferenceType>(e.number_of_terms()), e.begin());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.end()-e.begin(), static_cast<PointerDifferenceType>(e.number_of_terms()));

    // Test derefencing of iterators
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e.begin()", e.begin());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "*e.begin()", *e.begin());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.begin()->index(), MultiIndex({0,0,0}));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.begin()->coefficient(), 2.0_x);

    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "ce.begin()", ce.begin());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "*ce.begin()", *ce.begin());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), ce.begin()->index(), MultiIndex({0,0,0}));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), ce.begin()->coefficient(), 2.0_x);


    // The behaviour of iterators is rather odd and not what might be expected
    // A MultiIndex reference assigned to by iter->index() changes its value
    // when the ExpansionIterator is incremented, but a F reference does not.
    // This behaviour should be changed in future versions if technologically
    // feasible.
    {
        ExpansionIterator iter=e.begin();
        ExpansionIndexReference aref=iter->index();
        ExpansionCoefficientConstReference xref=iter->coefficient();
        MultiIndex a1=iter->index();
        F x1=iter->coefficient();
        ++iter;
        MultiIndex a2=iter->index();
        F x2=iter->coefficient();
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), a1==aref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), x1==xref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), a2!=aref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), x2!=xref);
    }

    {
        ExpansionConstIterator citer=e.begin();
        ExpansionIndexConstReference caref=citer->index();
        ExpansionCoefficientConstReference cxref=citer->coefficient();
        MultiIndex ca1=citer->index();
        F cx1=citer->coefficient();
        ++citer;
        MultiIndex ca2=citer->index();
        F cx2=citer->coefficient();
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), ca1==caref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), cx1==cxref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), ca2!=caref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), cx2!=cxref);
    }

    {
        ExpansionIterator iter=e.begin()+2;
        ExpansionPointer mptr=iter.operator->();
        ExpansionReference mref=iter.operator*();
        ExpansionValueType m1val=*iter;
        ++iter;
        ExpansionValueType m2val=*iter;
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), m1val==mref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), m1val==*mptr);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), m2val!=mref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), m2val!=*mptr);
    }

    {
        ExpansionConstIterator iter=e.end()-1;
        ExpansionConstPointer mptr=iter.operator->();
        ExpansionConstReference mref=iter.operator*();
        ExpansionValueType m1val=*iter;
        ++iter;
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), m1val==mref);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), m1val==*mptr);
    }

    {
        ExpansionConstIterator src=ce.begin()+1;
        ExpansionIterator dst=e.begin();
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), *dst=*src);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), dst->index(), src->index());
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), dst->coefficient(), src->coefficient());
    }

    // Test finding of values of iterators
    ExpansionIterator iter=e.begin();
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), iter=e.find(MultiIndex({1,0,0})));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), iter->index(), MultiIndex({1,0,0}));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), iter->coefficient(), 3.0_x);



    // Test hand-coded swap of values
    ExpansionIterator iter1(e.begin()+1);
    ExpansionIterator iter2(e.begin()+3);

    // Perform swap
    ExpansionValueType tmp(*iter2);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), tmp.index()==MultiIndex({1,0,1}));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), tmp.coefficient()==7.0_x);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), *iter2=*iter1);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), iter2->index()==MultiIndex({1,0,0}));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), iter2->coefficient()==3.0_x);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), *iter1=tmp);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), iter1->index()==MultiIndex({1,0,1}));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), iter1->coefficient()==7.0_x);


}

template<class F> Void TestExpansion<F>::test_equality()
{
    MultiIndex a(2);
    MultiIndex b(2); ++b;
    Expansion<MI,F> e1(2,prec),e2(2,prec);
    e1.append(a,1.0_x); e1.append(b,2.0_x);
    e2.append(a,1.0_x); e2.append(b,3.0_x);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE_CTX(AlgebraTesting::context(), same, e1, e2);
    e2.clear(); e2.append(a,1.0_x); e2.append(b,2.0_x);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE_CTX(AlgebraTesting::context(), same, e1, e2);
    e1.clear(); e1.append(b,2.0_x);
    e2.clear(); e2.append(a,0.0_x); e2.append(b,2.0_x);
    if(!same(e1,e2)) { ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "Expansion<MI,F> objects differing by explicit zeros are considered not same."); }
    e1.clear(); e1.append(a,-0.0_x);
    e1.clear(); e1.append(a,+0.0_x);
    if(!same(e1,e2)) { ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "Expansion<MI,F> objects differing by +0 versus -0 coefficients are considered not same."); }
    e1.clear(); e1.append(a,1.0_x); e1.append(b,2.0_x);
    e2.clear(); e2.append(b,2.0_x); e2.append(a,1.0_x);
    if(!same(e1,e2)) { ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "Expansion<MI,F> objects differing by order of set operators are considered not same."); }
}

template<class F> Void TestExpansion<F>::test_sort()
{
    MultiIndexList as1({{0,0},{1,0},{2,0}});
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "as1", as1);
    MultiIndexList as2({{2,0},{0,0},{1,0}});
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "as1", as1);

    GradedIndexLess graded_index_less;
    Expansion<MI,F> e2s({{{0,0},5.0_x},{{1,0},11.0_x}},prec);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e2s", e2s)
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), e2s.is_sorted(graded_index_less));
    Expansion<MI,F> e2u({{{1,0},2.0_x},{{0,0},5.0_x}},prec);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e2u", e2u)
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not e2u.is_sorted(graded_index_less));

//    ARIADNE_TEST_CONSTRUCT( Expansion<MI,F>, e1, ({{{0,1},1.0_x},{{1,1},1.0_x},{{2,1},1.0_x},{{0,1},11.0_x},{{0,1},11.0_x},{{1,1},11.0_x}},prec) );
    Expansion<MI,F> e1({{{0,1},1.0_x},{{1,1},1.0_x},{{2,1},1.0_x},{{0,1},11.0_x},{{0,1},11.0_x},{{1,1},11.0_x}},prec);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e1", e1);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not e1.is_sorted(graded_index_less));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e1.sort(graded_index_less));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), e1.is_sorted(graded_index_less));
}

template<class F> Void TestExpansion<F>::test_cleanup()
{
    // Test to see if the cleanup/sort operations work.
    // Since these are used in the constructors, we can't use the main constructors to test this
    MultiIndex a(3);
    MultiIndex b(3); ++b;
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "a", a);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "b", b);

    Expansion<MI,F> e(3,prec);
    F c(1,prec);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    for(Nat i=0; i!=2; ++i) {
        c=hlf(c);
        if(i%2) { e.append(a,c); ++b; ++b; a=b; ++b; } else { e.append(b,c);}
        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    }

    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.index_sort(graded_less));
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.combine_terms());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), e.remove_zeros());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);

}

template<class F> Void TestExpansion<F>::test_constructors()
{
    // Empty initialiser list causes failure
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), ExpansionType e0({}));

    InitializerList<Pair<typename ExpansionType::IndexInitializerType,ExactDouble>> empty_exact;
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), ExpansionType(empty_exact,prec));

    // Empty expansion
    ExpansionType e1(3,prec);
    // Expansion with all entries; useful for checking ordering of indices
    //ARIADNE_TEST_CONSTRUCT(Expansion<MI,F>,e2,(3,4, {1., 2.,3.,4., 5.,6.,7.,8.,9.,10.,
    //    11.,12.,13.,14.,15.,16.,17.,18.,19.,20., 21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35}));

    // Dense expansion
    ExpansionType e3({{{2,0,0},5.0_x},{{1,1,0},2.0_x},{{1,0,1},0.0_x},{{0,2,0},0.0_x},{{0,1,1},3.0_x},{{0,2,0},0.0_x}}, prec);
    //ARIADNE_TEST_CONSTRUCT(Expansion<MI,F>,e3,(3,2, {0., 0.,0.,0., 5.,2.,0.,0.,3.,0.}));
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e3", e3);
    ARIADNE_TESTING_EXPECT_NOT_EQUAL_CTX(AlgebraTesting::context(), e3.find(MultiIndex({2,0,0})), e3.end());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e3[MultiIndex({2,0,0})], 5.0_x);

    // Sparse expansion with unordered indiced
    Expansion<MI,F> pp3({{{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x}}, prec);
    ExpansionType p3({{{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x}}, prec);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), p3[MultiIndex({1,2})], 5.0_x);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), p3[MultiIndex({0,0})], 2.0_x);

    // Unordered indices
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE_CTX(AlgebraTesting::context(), same, ExpansionType({{{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x}}, prec), ExpansionType({{{0,0},2.0_x}, {{1,0},3.0_x}, {{0,1},11.0_x}, {{3,0},7.0_x}, {{1,2},5.0_x}}, prec));
    // Repeated indices; do not sum in expansion class
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE_CTX(AlgebraTesting::context(), same, ExpansionType({{{1,0},2.0_x}, {{1,0},3.0_x}, {{0,2},7.0_x}}, prec), ExpansionType({{{1,0},5.0_x}, {{0,2},7.0_x}}, prec));

    // Regression tests for expansions with only preces
    ExpansionType pr2({{{},0.0_x}},prec);
    ExpansionType pr1({{{3,0,0},0.0_x}},prec);

    // Regression tests for higher-order expansions
    ExpansionType ho1({{{0,1,0,0,0},2.0_x}, {{0,1,0,0,1},3.0_x}, {{2,0,1,0,0},5.0_x}, {{0,0,0,0,0},7.0_x}},prec);
}


template<class F> Void TestExpansion<F>::test_find()
{
    ExpansionType e({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} }, prec);
    MultiIndex a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e", e);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "e.find(a)-e.begin()", e.find(a)-e.begin());
    ARIADNE_TESTING_EXPECT_NOT_EQUAL_CTX(AlgebraTesting::context(), e.find(a), e.end());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.find(a)->index(), a);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.find(a)->coefficient(), 5.0_x);
    a[1]=1;
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), e.find(a), e.end());

}

template<class F> Void TestExpansion<F>::test_embed()
{
    ExpansionType e({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },prec);
    ARIADNE_TESTING_EXPECT_SAME_CTX(AlgebraTesting::context(), embed(0,e,2), ExpansionType({ {{1,2,0,0},5.0_x}, {{0,0,0,0},2.0_x}, {{1,0,0,0},3.0_x}, {{3,0,0,0},7.0_x}, {{0,1,0,0},11.0_x} },prec));
    ARIADNE_TESTING_EXPECT_SAME_CTX(AlgebraTesting::context(), embed(1,e,0), ExpansionType({ {{0,1,2},5.0_x}, {{0,0,0},2.0_x}, {{0,1,0},3.0_x}, {{0,3,0},7.0_x}, {{0,0,1},11.0_x} },prec));
    ARIADNE_TESTING_EXPECT_SAME_CTX(AlgebraTesting::context(), embed(1,e,2), ExpansionType({ {{0,1,2,0,0},5.0_x}, {{0,0,0,0,0},2.0_x}, {{0,1,0,0,0},3.0_x}, {{0,3,0,0,0},7.0_x}, {{0,0,1,0,0},11.0_x} },prec));

}

namespace {

// A strict weak ordering that groups different monomials of equal degree.
struct DegreeIndexLess {
    template<class T> Bool operator()(T const& lhs, T const& rhs) const {
        return lhs.index().degree()<rhs.index().degree();
    }
};

Void test_sorted_equality_edge_cases() {
    using E=SortedExpansion<MultiIndex,RoundedFloatDP,GradedIndexLess>;
    E two_variables(2u,dp);
    E three_variables(3u,dp);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(two_variables==three_variables));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(three_variables==two_variables));

    two_variables.insert(MultiIndex({1u,0u}),RoundedFloatDP(2,dp));
    E original=two_variables;
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, two_variables.insert(MultiIndex({1u,0u}),RoundedFloatDP(3,dp)));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), same(two_variables,original));

    using Reverse=SortedExpansion<MultiIndex,RoundedFloatDP,ReverseLexicographicIndexLess>;

    Reverse reverse_two_variables(2u,dp);
    Reverse reverse_three_variables(3u,dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(reverse_two_variables==reverse_three_variables));

    Reverse reverse_left_only(2u,dp);
    reverse_left_only.insert(MultiIndex({1u,0u}),RoundedFloatDP(0,dp));
    Reverse reverse_empty(2u,dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(reverse_left_only==reverse_empty));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(reverse_empty==reverse_left_only));

    Reverse reverse_first(2u,dp);
    Reverse reverse_second(2u,dp);
    reverse_first.insert(MultiIndex({1u,0u}),RoundedFloatDP(0,dp));
    reverse_second.insert(MultiIndex({0u,1u}),RoundedFloatDP(0,dp));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(reverse_first==reverse_second));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), (void)(reverse_second==reverse_first));

    using Grouped=SortedExpansion<MultiIndex,RoundedFloatDP,DegreeIndexLess>;
    Grouped grouped_empty(2u,dp), grouped_three_variables(3u,dp);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(grouped_empty==grouped_three_variables));

    // Unequal degrees must also compare correctly with this ordering,
    // including unmatched zero and nonzero terms on either side.
    Grouped lower_degree(2u,dp), higher_degree(2u,dp);
    lower_degree.append(MultiIndex({1u,0u}),RoundedFloatDP(0,dp));
    higher_degree.append(MultiIndex({2u,0u}),RoundedFloatDP(0,dp));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lower_degree==higher_degree);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), higher_degree==lower_degree);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lower_degree==grouped_empty);
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), grouped_empty==lower_degree);
    lower_degree.front().coefficient()=1;
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(lower_degree==higher_degree));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(higher_degree==lower_degree));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(lower_degree==grouped_empty));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(grouped_empty==lower_degree));

    for(Int left_value=0; left_value!=2; ++left_value) {
        for(Int right_value=0; right_value!=2; ++right_value) {
            Grouped left(2u,dp), right(2u,dp);
            left.append(MultiIndex({1u,0u}),RoundedFloatDP(left_value,dp));
            right.append(MultiIndex({0u,1u}),RoundedFloatDP(right_value,dp));
            // A following common term verifies that comparison advances past
            // the distinct indices equivalent under DegreeIndexLess.
            left.append(MultiIndex({2u,0u}),RoundedFloatDP(3,dp));
            right.append(MultiIndex({2u,0u}),RoundedFloatDP(3,dp));
            Bool expected=(left_value==0 && right_value==0);
            ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), decide(left==right), expected);
            ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), decide(right==left), expected);
            right.back().coefficient()=4;
            ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(left==right));
        }
    }
}

Void test_boolean_sorted_equality_short_circuit() {
    // Built-in double equality is Boolean: unlike multi-valued equality,
    // each conjunction must preserve the short-circuit behavior.
    using E=SortedExpansion<MultiIndex,double,DegreeIndexLess>;
    using Term=std::pair<MultiIndex,double>;
    auto make_degree_entry=[](DegreeType first, DegreeType second, double coefficient) -> Term {
        return {MultiIndex({first,second}),coefficient};
    };
    auto build_sorted=[](std::initializer_list<Term> entries) {
        E expansion(2u,std::tuple<>());
        for(auto const& entry : entries) { expansion.append(entry.first,entry.second); }
        return expansion;
    };
    auto check_both_orders=[](E const& lhs, E const& rhs, bool expected) {
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lhs==rhs, expected);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), rhs==lhs, expected);
    };

    E empty_expansion=build_sorted({});
    E different_dimension(3u,std::tuple<>());
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not (empty_expansion==different_dimension));
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not (different_dimension==empty_expansion));

    E matched=build_sorted({make_degree_entry(0u,0u,1.0),
                            make_degree_entry(2u,0u,4.0)});
    check_both_orders(matched,matched,true);
    E different_coefficient=build_sorted({make_degree_entry(0u,0u,2.0),
                                          make_degree_entry(2u,0u,4.0)});
    check_both_orders(matched,different_coefficient,false);

    // Both unmatched-term branches: zero preserves equality, nonzero rejects it.
    E left_unmatched=build_sorted({make_degree_entry(1u,0u,0.0),
                                   make_degree_entry(3u,0u,5.0)});
    E right_unmatched=build_sorted({make_degree_entry(2u,0u,0.0),
                                    make_degree_entry(3u,0u,5.0)});
    check_both_orders(left_unmatched,right_unmatched,true);
    left_unmatched.front().coefficient()=1.0;
    check_both_orders(left_unmatched,right_unmatched,false);
    left_unmatched.front().coefficient()=0.0;
    right_unmatched.front().coefficient()=1.0;
    check_both_orders(left_unmatched,right_unmatched,false);

    // Once equality is false, later comparisons must not restore it.
    E prior_mismatch_left=build_sorted({make_degree_entry(0u,0u,1.0),
                                        make_degree_entry(1u,0u,0.0),
                                        make_degree_entry(3u,0u,1.0)});
    E prior_mismatch_right=build_sorted({make_degree_entry(0u,0u,2.0),
                                         make_degree_entry(2u,0u,0.0),
                                         make_degree_entry(3u,0u,1.0)});
    check_both_orders(prior_mismatch_left,prior_mismatch_right,false);

    // DegreeIndexLess considers distinct indices of equal total degree
    // equivalent under the ordering, but not necessarily as monomials.
    for(double lhs_coefficient : {0.0,1.0}) {
        for(double rhs_coefficient : {0.0,1.0}) {
            E lhs=build_sorted({make_degree_entry(1u,0u,lhs_coefficient)});
            E rhs=build_sorted({make_degree_entry(0u,1u,rhs_coefficient)});
            check_both_orders(lhs,rhs,lhs_coefficient==0.0 && rhs_coefficient==0.0);
        }
    }
    E prior_mismatch_grouped_left=build_sorted({make_degree_entry(0u,0u,1.0),
                                                 make_degree_entry(1u,0u,0.0)});
    E prior_mismatch_grouped_right=build_sorted({make_degree_entry(0u,0u,2.0),
                                                  make_degree_entry(0u,1u,0.0)});
    check_both_orders(prior_mismatch_grouped_left,prior_mismatch_grouped_right,false);

    // Both tail loops: compare a zero/nonzero tail against an empty operand.
    E zero_tail=build_sorted({make_degree_entry(1u,0u,0.0)});
    E nonzero_tail=build_sorted({make_degree_entry(1u,0u,1.0)});
    check_both_orders(empty_expansion,zero_tail,true);
    check_both_orders(empty_expansion,nonzero_tail,false);

    // The tail must also be visited after a preceding mismatching term.
    E prior_mismatch_tail_left=build_sorted({make_degree_entry(0u,0u,1.0),
                                              make_degree_entry(1u,0u,0.0)});
    E prior_mismatch_tail_right=build_sorted({make_degree_entry(0u,0u,2.0)});
    check_both_orders(prior_mismatch_tail_left,prior_mismatch_tail_right,false);
}

Void test_exact_coefficient_helpers() {
    Expansion<MultiIndex,double> plain(2u,std::tuple<>());
    ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), plain.coefficient_characteristics()==std::tuple<>());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), plain.characteristics().first, 2u);

    // Raw floats do not have inplace addition: combining them must fail
    // without losing either of the duplicate terms.
    static_assert(!CanInplaceAdd<FloatDP,FloatDP>);
    Expansion<MultiIndex,FloatDP> raw(2u,dp);
    raw.append(MultiIndex({1u,0u}),FloatDP(2,dp));
    raw.append(MultiIndex({1u,0u}),FloatDP(3,dp));
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, raw.combine_terms());
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), raw.size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), raw.front().coefficient(), FloatDP(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), raw.back().coefficient(), FloatDP(3,dp));
}

template<class X, class PR>
Void test_additional_expansion_instantiation(PR const& pr)
{
    X coefficient_zero(0u,pr);
    X coefficient_one(1u,pr);

    Expansion<MultiIndex,X> init_case({
        {{0u,0u},coefficient_zero},
        {{1u,0u},coefficient_one}
    });
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), init_case.size(), 1u);

    Expansion<MultiIndex,X> source(2u,pr);
    source.append(MultiIndex({1u,0u}),coefficient_one);

    Expansion<MultiIndex,X> copy_target(2u,pr);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), copy_target=source);
    Expansion<MultiIndex,X>& copy_alias=copy_target;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), copy_target=copy_alias);

    Expansion<MultiIndex,X> move_target(2u,pr);
    Expansion<MultiIndex,X> move_source(source);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), move_target=std::move(move_source));
    Expansion<MultiIndex,X>* move_alias=&move_target;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), move_target=std::move(*move_alias));
}

template<class X, class PR>
Void test_univariate_expansion_assignment_instantiation(PR const& pr)
{
    X coefficient_one(1u,pr);

    Expansion<UniIndex,X> source(SizeOne(),pr);
    source.append(UniIndex(1u),coefficient_one);

    Expansion<UniIndex,X> copy_target(SizeOne(),pr);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), copy_target=source);
    Expansion<UniIndex,X>& copy_alias=copy_target;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), copy_target=copy_alias);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), copy_target.size(), 1u);

    Expansion<UniIndex,X> move_target(SizeOne(),pr);
    Expansion<UniIndex,X> move_source(source);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), move_target=std::move(move_source));
    Expansion<UniIndex,X>* move_alias=&move_target;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), move_target=std::move(*move_alias));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), move_target.size(), 1u);
}

template<class X, class PR>
Void test_graded_check_instantiation(PR const& pr)
{
    using E=SortedExpansion<MultiIndex,X,GradedIndexLess>;
    X one(1u,pr);

    E sorted(2u,pr);
    sorted.append(MultiIndex({0u,0u}),one);
    sorted.append(MultiIndex({1u,0u}),one);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), sorted.check());

    E unsorted(2u,pr);
    unsorted.append(MultiIndex({1u,0u}),one);
    unsorted.append(MultiIndex({0u,0u}),one);
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, unsorted.check());
}

template<class X, class PR>
Void test_univariate_sorted_equality_ordering(PR const& pr)
{
    // A sparse comparison must visit unmatched terms in both directions,
    // and then resume at a matching term without losing a previous result.
    using E=SortedExpansion<UniIndex,X,IndexLess>;
    X zero(0u,pr), one(1u,pr), coefficient_two(2u,pr);

    auto check_symmetric = [](E const& lhs, E const& rhs) {
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), AriadneTesting::truth(lhs==rhs)==AriadneTesting::truth(rhs==lhs));
    };

    E empty(SizeOne(),pr);
    E left_only(SizeOne(),pr), right_only(SizeOne(),pr);
    left_only.append(UniIndex(0u),zero);
    right_only.append(UniIndex(2u),zero);
    check_symmetric(left_only,empty);
    check_symmetric(empty,right_only);
    check_symmetric(left_only,right_only);

    E left(SizeOne(),pr), right(SizeOne(),pr);
    left.append(UniIndex(0u),zero);
    left.append(UniIndex(2u),one);
    right.append(UniIndex(1u),zero);
    right.append(UniIndex(2u),one);
    check_symmetric(left,right);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), left.size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), right.size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), left.back().index(), UniIndex(2u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), right.back().index(), UniIndex(2u));

    right.back().coefficient()=coefficient_two;
    check_symmetric(left,right);

    E equal_left(SizeOne(),pr), equal_right(SizeOne(),pr);
    equal_left.append(UniIndex(1u),one);
    equal_right.append(UniIndex(1u),one);
    check_symmetric(equal_left,equal_right);
    equal_right.front().coefficient()=coefficient_two;
    check_symmetric(equal_left,equal_right);
}

Void test_unit_coefficient_output() {
    // Built-in doubles have decidable equality, so unit coefficients can be
    // omitted while unit constants must still be printed.
    Expansion<MultiIndex,double> multi(2u,std::tuple<>());
    multi.append(MultiIndex({1u,1u}),1.0);
    multi.append(MultiIndex({0u,2u}),-1.0);
    multi.append(MultiIndex({0u,0u}),1.0);
    multi.append(MultiIndex({1u,0u}),2.0);
    std::ostringstream multi_stream;
    multi._write(multi_stream,Array<String>({"x","y"}));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), multi_stream.str(), String(" x*y -y^2 +1 +2*x"));

    Expansion<UniIndex,double> uni(SizeOne{},std::tuple<>{});
    uni.append(UniIndex(1u),1.0);
    uni.append(UniIndex(2u),-1.0);
    uni.append(UniIndex(0u),1.0);
    uni.append(UniIndex(3u),2.0);
    std::ostringstream uni_stream;
    uni._write(uni_stream,String("t"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), uni_stream.str(), String(" t -t^2 +1 +2*t^3"));
}

Void test_bounds_coefficient_output() {
    // Singleton bounds can prove that a coefficient is exactly one, unlike
    // approximate coefficients. Check omission, signs, constants and zero.
    Expansion<MultiIndex,FloatDPBounds> expansion(2u,dp);
    Array<String> names({"x","y"});
    std::ostringstream empty_stream;
    expansion._write(empty_stream,names);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), empty_stream.str(), String("0"));

    FloatDPBounds one(1,dp), coefficient_two(2,dp);
    expansion.append(MultiIndex({1u,1u}),one);
    expansion.append(MultiIndex({0u,2u}),FloatDPBounds(-1,dp));
    expansion.append(MultiIndex({0u,0u}),one);
    expansion.append(MultiIndex({1u,0u}),coefficient_two);
    std::ostringstream stream, expected;
    expansion._write(stream,names);
    expected << " x*y -y^2 +" << one << " +" << coefficient_two << "*x";
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), stream.str(), expected.str());
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("expansion");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_sorted_equality_edge_cases());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_boolean_sorted_equality_short_circuit());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_exact_coefficient_helpers());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_unit_coefficient_output());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_bounds_coefficient_output());
    {
        Expansion<MultiIndex,ExactDouble> exact_expansion({{{0u,0u},ExactDouble(1.0)},{{0u,0u},ExactDouble(2.0)}});
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), exact_expansion.coefficient_characteristics());
    }

    {
        Expansion<MultiIndex,double> double_expansion({{{1u,0u},-2.0},{{0u,0u},1.0}});
        std::ostringstream double_stream;
        Array<String> double_names({"x","y"});
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), double_expansion._write(double_stream,double_names));
    }

    {
        RoundedFloatDP coefficient_zero(0,double_precision);
        RoundedFloatDP coefficient_one(1,double_precision);
        RoundedFloatDP coefficient_two(2,double_precision);
        Expansion<UniIndex,RoundedFloatDP> univariate({{DegreeType(0u),coefficient_one},{DegreeType(1u),coefficient_zero},{DegreeType(2u),coefficient_two}});
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), univariate.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), univariate[UniIndex(2u)], coefficient_two);
        univariate.append(UniIndex(1u),coefficient_one);
        std::ostringstream univariate_stream;
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), univariate._write(univariate_stream,String("x")));
    }

    FloatDP raw_dp(0,double_precision);
    FloatMP raw_mp(0,MultiplePrecision(128_bits));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), characteristics(raw_dp));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), characteristics(raw_mp));

    {
        Expansion<MultiIndex,FloatDPBounds> source(2u,double_precision);
        source.append({0u,0u},FloatDPBounds(1,double_precision));
        source.append({1u,0u},FloatDPBounds(0,double_precision));
        Expansion<MultiIndex,FloatDPApproximation> approximation(source);
        Expansion<MultiIndex,FloatDPUpperInterval> upper(source);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), approximation.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), upper.size(), 1u);
    }

    {
        auto mp=MultiplePrecision(128_bits);
        Expansion<MultiIndex,FloatMPBounds> source(2u,mp);
        source.append({0u,0u},FloatMPBounds(1,mp));
        source.append({1u,0u},FloatMPBounds(0,mp));
        Expansion<MultiIndex,FloatMPApproximation> approximation(source);
        Expansion<MultiIndex,FloatMPUpperInterval> upper(source);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), approximation.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), upper.size(), 1u);
    }

    test_additional_expansion_instantiation<RoundedFloatDP>(double_precision);
    test_additional_expansion_instantiation<FloatDPApproximation>(double_precision);
    test_additional_expansion_instantiation<FloatDPBounds>(double_precision);
    test_additional_expansion_instantiation<FloatDPUpperInterval>(double_precision);

    test_univariate_expansion_assignment_instantiation<RoundedFloatDP>(double_precision);
    test_univariate_expansion_assignment_instantiation<FloatDPApproximation>(double_precision);
    test_univariate_expansion_assignment_instantiation<FloatDPBounds>(double_precision);

    auto multiple_precision=MultiplePrecision(128_bits);
    test_additional_expansion_instantiation<RoundedFloatMP>(multiple_precision);
    test_additional_expansion_instantiation<FloatMPApproximation>(multiple_precision);
    test_additional_expansion_instantiation<FloatMPBounds>(multiple_precision);
    test_additional_expansion_instantiation<FloatMPUpperInterval>(multiple_precision);

    test_univariate_expansion_assignment_instantiation<FloatMPApproximation>(multiple_precision);
    test_univariate_expansion_assignment_instantiation<FloatMPBounds>(multiple_precision);

    test_univariate_sorted_equality_ordering<FloatDPApproximation>(double_precision);
    test_univariate_sorted_equality_ordering<FloatMPApproximation>(multiple_precision);

    test_graded_check_instantiation<FloatDPApproximation>(double_precision);
    test_graded_check_instantiation<FloatDPBounds>(double_precision);
    test_graded_check_instantiation<FloatDPUpperInterval>(double_precision);
    test_graded_check_instantiation<FloatMPApproximation>(multiple_precision);
    test_graded_check_instantiation<FloatMPBounds>(multiple_precision);
    test_graded_check_instantiation<FloatMPUpperInterval>(multiple_precision);

    {
        Expansion<MultiIndex,FloatDPUpperInterval> dp_upper_reserve(2u,double_precision);
        SizeType dp_upper_capacity=dp_upper_reserve.capacity();
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), dp_upper_reserve.reserve(dp_upper_capacity+1u));

        Expansion<MultiIndex,FloatMPBounds> mp_bounds_reserve(2u,multiple_precision);
        SizeType mp_bounds_capacity=mp_bounds_reserve.capacity();
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), mp_bounds_reserve.reserve(mp_bounds_capacity+1u));

        Expansion<UniIndex,FloatMPApproximation> uni_mp_reserve(SizeOne(),multiple_precision);
        SizeType uni_mp_capacity=uni_mp_reserve.capacity();
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), uni_mp_reserve.reserve(uni_mp_capacity));
    }

    RoundedFloatDP zero_dp{0,double_precision};
    RoundedFloatMP zero_mp{0,multiple_precision};
    TestExpansion<RoundedFloatDP>(zero_dp).test();
    TestExpansion<Rounded<FloatMP>>(zero_mp).test();
    

    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
