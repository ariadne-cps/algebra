/***************************************************************************
 *            test_multi_index.cpp
 *
 *  Copyright  2007-20  Pieter Collins
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
#include <fstream>

#include "algebra/multi_index.hpp"

#include "utility/test.hpp"

using namespace std;
using namespace Ariadne;


class TestMultiIndex
{
  public:
    TestMultiIndex() { }

    Void test() {
        ARIADNE_TEST_CALL(test_constructor());
        ARIADNE_TEST_CALL(test_uni_index());
        ARIADNE_TEST_CALL(test_inline_helpers());
        ARIADNE_TEST_CALL(test_preconditions());
        ARIADNE_TEST_CALL(test_comparison());
        ARIADNE_TEST_CALL(test_addition());
        ARIADNE_TEST_CALL(test_increment());
        ARIADNE_TEST_CALL(test_list());
    }

    Void test_constructor() {
        DegreeType a1[3]={2,0,3};
        MultiIndex i1(3u,a1);
        ARIADNE_TEST_EQUAL((Int)i1.number_of_variables(),3);
        ARIADNE_TEST_EQUAL((Int)i1.degree(),5);
        ARIADNE_TEST_EQUAL((Int)i1[0],2);
        ARIADNE_TEST_EQUAL((Int)i1[1],0);
        ARIADNE_TEST_EQUAL((Int)i1[2],3);

        DegreeType a2[5]={2,0,3,0,5};
        MultiIndex i2(5,a2);
        ARIADNE_TEST_EQUAL((Int)i2.number_of_variables(),5);
        ARIADNE_TEST_EQUAL((Int)i2.degree(),10);
        ARIADNE_TEST_EQUAL((Int)i2[3],0);
        ARIADNE_TEST_EQUAL((Int)i2[4],5);

        MultiIndex i3=MultiIndex::unit(9u,2u);
        ARIADNE_TEST_EQUAL((Int)i3.number_of_variables(),9);
        ARIADNE_TEST_EQUAL((Int)i3.degree(),1);
        ARIADNE_TEST_EQUAL((Int)i3[2],1);
        ARIADNE_TEST_EQUAL((Int)i3[8],0);
    }

    Void test_uni_index() {
        UniIndex a(2u);
        ARIADNE_TEST_EQUALS(a.degree(),DegreeType(2u));
        ARIADNE_TEST_EXECUTE(a.size());
        ARIADNE_TEST_EQUALS(a[IndexZero()],DegreeType(2u));

        a[IndexZero()]=3u;
        const UniIndex& ca=a;
        ARIADNE_TEST_EQUALS(ca[IndexZero()],DegreeType(3u));

        ARIADNE_TEST_EXECUTE(++a);
        ARIADNE_TEST_EQUALS(a.degree(),DegreeType(4u));
        ARIADNE_TEST_EXECUTE(--a);
        ARIADNE_TEST_EQUALS(a.degree(),DegreeType(3u));

        UniIndex b(2u);
        ARIADNE_TEST_EQUALS(a-b,static_cast<short int>(1));
        ARIADNE_TEST_EXECUTE(a+=b);
        ARIADNE_TEST_EQUALS(a.degree(),DegreeType(5u));
    }

    Void test_inline_helpers() {
        MultiIndex a({1u,2u});
        MultiIndex b({3u,4u});

        ARIADNE_TEST_ASSERT(a.begin()!=a.end());
        const MultiIndex& ca=a;
        ARIADNE_TEST_ASSERT(ca.begin()!=ca.end());

        MultiIndex assigned(2u);
        ARIADNE_TEST_EXECUTE(assigned.assign(a));
        ARIADNE_TEST_EQUALS(assigned,a);

        ARIADNE_TEST_EQUALS(DegreeType(2u)*a,MultiIndex({2u,4u}));

        GradedLess graded_less_object;
        LexicographicLess lexicographic_less_object;
        ARIADNE_TEST_ASSERT(graded_less_object(DegreeType(1u),DegreeType(2u)));
        ARIADNE_TEST_ASSERT(lexicographic_less_object(DegreeType(1u),DegreeType(2u)));
        ARIADNE_TEST_ASSERT(lexicographic_less_object(a,b));
    }

    Void test_preconditions() {
        MultiIndex a({1u,2u});
        const MultiIndex& ca=a;
        ARIADNE_TEST_FAIL(a[a.size()]);
        ARIADNE_TEST_FAIL(ca[ca.size()]);
        ARIADNE_TEST_FAIL(ca.get(ca.size()));
        ARIADNE_TEST_FAIL(a.set(a.size(),0u));

        MultiIndex short_index({1u});
        ARIADNE_TEST_ASSERT(!(a==short_index));
        ARIADNE_TEST_FAIL(a.assign(short_index));
        ARIADNE_TEST_FAIL(graded_less(a,short_index));
        ARIADNE_TEST_FAIL(lexicographic_less(a,short_index));
        ARIADNE_TEST_FAIL(reverse_lexicographic_less(a,short_index));

        MultiIndex long_index({2u,0u});
        ARIADNE_TEST_FAIL(long_index+=short_index);
        ARIADNE_TEST_FAIL(long_index-=short_index);
        ARIADNE_TEST_FAIL(MultiIndex({0u,0u})-=MultiIndex({1u,0u}));
        ARIADNE_TEST_FAIL(swap(long_index,short_index));

        MultiIndexReference long_reference(long_index);
        MultiIndexReference short_reference(short_index);
        ARIADNE_TEST_FAIL(long_reference=short_index);
        ARIADNE_TEST_FAIL(long_reference+=short_index);
        ARIADNE_TEST_FAIL(swap(long_reference,short_reference));

        MultiIndex same_size({3u,4u});
        MultiIndex& same_alias=same_size;
        ARIADNE_TEST_EXECUTE(same_size=same_alias);
        ARIADNE_TEST_EXECUTE(same_size.resize(same_size.size()));
        ARIADNE_TEST_EXECUTE(same_size.resize(3u));
        ARIADNE_TEST_EQUALS(same_size.size(),3u);
        ARIADNE_TEST_EXECUTE(same_size.clear());
    }

    Void test_comparison() {
        DegreeType a1p[3]={2,0,3};
        DegreeType a2p[3]={1,0,4};
        DegreeType a3p[3]={3,0,1};
        ARIADNE_TEST_CONSTRUCT(MultiIndex,a1,(3,a1p));
        ARIADNE_TEST_CONSTRUCT(MultiIndex,a2,(3,a2p));
        ARIADNE_TEST_CONSTRUCT(MultiIndex,a3,(3,a3p));

        ARIADNE_TEST_ASSERT(graded_less(a1,a2));
        ARIADNE_TEST_ASSERT(graded_less(a3,a1));
        ARIADNE_TEST_ASSERT(graded_less(a3,a2));
        ARIADNE_TEST_ASSERT(!graded_less(a2,a1));

        ARIADNE_TEST_ASSERT(lexicographic_less(a2,a1));
        ARIADNE_TEST_ASSERT(lexicographic_less(a1,a3));
        ARIADNE_TEST_ASSERT(lexicographic_less(a2,a3));
        ARIADNE_TEST_ASSERT(!lexicographic_less(a1,a1));

        ARIADNE_TEST_ASSERT(reverse_lexicographic_less(a2,a1));
        ARIADNE_TEST_ASSERT(reverse_lexicographic_less(a1,a3));
        ARIADNE_TEST_ASSERT(reverse_lexicographic_less(a2,a3));

        ARIADNE_TEST_EXECUTE(swap(a1,a2));
        ARIADNE_TEST_EQUALS(a1,MultiIndex(3,a2p));
        ARIADNE_TEST_EQUALS(a2,MultiIndex(3,a1p));

        MultiIndexReference a1r(a1);
        MultiIndexReference a2r(a2);
        ARIADNE_TEST_EXECUTE(swap(a1r,a2r));
        ARIADNE_TEST_EQUALS(a1,MultiIndex(3,a1p));

    }

    Void test_addition() {
        ARIADNE_TEST_CONSTRUCT(MultiIndex,a1,({2,0,3}));
        ARIADNE_TEST_CONSTRUCT(MultiIndexReference,a1r,(a1));
        ARIADNE_TEST_CONSTRUCT(MultiIndex,a2,({1,1,0}));
        ARIADNE_TEST_CONSTRUCT(MultiIndex,a1pa2,({3,1,3}));
        ARIADNE_TEST_EQUALS(a1+a2,a1pa2);
        ARIADNE_TEST_EQUALS(a1pa2-a1,a2);
        ARIADNE_TEST_EQUALS(a1*2u,MultiIndex({4,0,6}));
        ARIADNE_TEST_EQUALS(a1+=a2,a1pa2);

        ARIADNE_TEST_EQUALS(a1r,a1);
        ARIADNE_TEST_EQUALS(a1r+=a2,a1pa2+a2);
    }

    Void test_increment() {
        MultiIndex empty;
        ARIADNE_TEST_FAIL(++empty);

        MultiIndex u(1);
        for(SizeType n=0; n!=8; ++n) {
            ARIADNE_TEST_EQUAL(u.degree(),n);
            ++u;
        }

        MultiIndex a(4);
        SizeType n=0;
        while(a.degree()<=5) {
            MultiIndex b=a; ++a; ++n;
            ARIADNE_TEST_BINARY_PREDICATE(graded_less,b,a);
            MultiIndex::IndexType d=0; for(SizeType i=0; i!=a.size(); ++i) { d+=a[i]; }
            ARIADNE_TEST_EQUAL(a.degree(),d);
        }
        ARIADNE_ASSERT_EQUAL(n,126);
    }

    Void test_list() {
        MultiIndex az={0,0,0};
        MultiIndex a0={2,0,3};
        MultiIndex a1={1,2,2};
        MultiIndex a2={1,0,4};
        MultiIndex a3={0,0,0};
        MultiIndex a4={1,0,1};

        ARIADNE_TEST_CONSTRUCT(MultiIndexList,lst,(0u,MultiIndex(3u)));
        ARIADNE_TEST_EQUALS(lst.capacity(),4u);
        ARIADNE_TEST_EQUALS(lst.size(),0u);
        ARIADNE_TEST_EQUALS(lst.argument_size(),3u);
        ARIADNE_TEST_EXECUTE(lst.append(a0));
        ARIADNE_TEST_PRINT(*lst.begin());
        ARIADNE_TEST_EQUALS(lst.begin()->size(),3u);
        const MultiIndexList& const_lst=lst;
        ARIADNE_TEST_EQUALS(const_lst.begin()->size(),3u);
        ARIADNE_TEST_EQUALS(*lst.begin(),a0);

        ARIADNE_TEST_EXECUTE(lst.append(a1));
        ARIADNE_TEST_EQUALS(*lst.begin(),a0);
        ARIADNE_TEST_EQUALS(*++lst.begin(),a1);

        const MultiIndexList& clst=lst;
        auto cit=clst.begin();
        auto cend=clst.end();
        ARIADNE_TEST_ASSERT(cit.operator!=(cend));
        auto cit_copy=cit;
        ARIADNE_TEST_ASSERT(cit_copy==cit);
        ARIADNE_TEST_EXECUTE(++cit_copy);
        ARIADNE_TEST_EXECUTE(--cit_copy);
        ARIADNE_TEST_EQUALS(*(cit+1),a1);
        ARIADNE_TEST_EQUALS(clst.front(),a0);
        ARIADNE_TEST_EQUALS(clst.back(),a1);

        ARIADNE_TEST_CONSTRUCT(MultiIndexList,lstc,(lst));
        ARIADNE_TEST_EQUAL(lstc,lst);

        ARIADNE_TEST_CONSTRUCT(MultiIndexList,lstm,(std::move(lst)));
        ARIADNE_TEST_EQUAL(lstm,lstc);

        ARIADNE_TEST_ASSIGN(lst,lstc);
        ARIADNE_TEST_EQUAL(lst,lstc);
        MultiIndex& lst_alias=lst[0u];
        (void)lst_alias;
        ARIADNE_TEST_EXECUTE(lst=lst);

        ARIADNE_TEST_ASSIGN(lstm,std::move(lstc));
        MultiIndexList* lstm_alias=&lstm;
        ARIADNE_TEST_EXECUTE(lstm=std::move(*lstm_alias));

        ARIADNE_TEST_PRINT(lst);
        ARIADNE_TEST_PRINT(lstm);
        ARIADNE_TEST_EQUAL(lstm,lst);
        ARIADNE_TEST_EXECUTE(lst.append(a2));
        ARIADNE_TEST_EXECUTE(lst.append(a3));
        ARIADNE_TEST_EQUALS(lst.back(),a3);
        ARIADNE_TEST_EQUALS(lst.capacity(),4u);
        ARIADNE_TEST_PRINT(&lst.front()[0]);
        ARIADNE_TEST_EXECUTE(lst.append(a4));
        ARIADNE_TEST_PRINT(&lst.front()[0]);
        ARIADNE_TEST_EQUALS(lst.capacity(),8u);
        ARIADNE_TEST_EQUALS(lst.size(),5u);
        ARIADNE_TEST_EQUALS(*--lst.end(),a4);
        ARIADNE_TEST_EQUALS(lst.back(),a4);

        ARIADNE_TEST_EXECUTE(lst.resize(3u));
        ARIADNE_TEST_EQUALS(lst.size(),3u);
        ARIADNE_TEST_EQUALS(lst.capacity(),8u);

        ARIADNE_TEST_EXECUTE(lst.resize(11u));
        ARIADNE_TEST_ASSERT(lst.capacity()==16u or lst.capacity()==11u);
        ARIADNE_TEST_EQUALS(lst.size(),11u);
        ARIADNE_TEST_EQUALS(lst[2],a2);
        ARIADNE_TEST_EQUALS(lst[3],az);
        ARIADNE_TEST_EQUALS(lst.back(),az);
        ARIADNE_TEST_FAIL(lst[lst.size()]);
        const MultiIndexList& const_resized_lst=lst;
        ARIADNE_TEST_FAIL(const_resized_lst[const_resized_lst.size()]);
        ARIADNE_TEST_PRINT(lst);

        ARIADNE_TEST_EXECUTE(*(lst.begin()+3)=a3);
        ARIADNE_TEST_EQUALS(lst[3],a3);
        ARIADNE_TEST_EQUALS(*(lst.begin()+3),a3);

        MultiIndexList summed(3u);
        ARIADNE_TEST_EXECUTE(summed.append_sum(a0,a1));
        ARIADNE_TEST_EQUALS(summed.front(),a0+a1);

        MultiIndexList repeated(2u,a0);
        ARIADNE_TEST_EQUALS(repeated.size(),2u);
        ARIADNE_TEST_EQUALS(repeated[0],a0);
        ARIADNE_TEST_EQUALS(repeated[1],a0);

        MultiIndexList from_multi({a0,a1,a2});
        ARIADNE_TEST_EQUALS(from_multi.size(),3u);
        ARIADNE_TEST_EQUALS(from_multi[1],a1);

        ARIADNE_TEST_FAIL(MultiIndexList(InitializerList<InitializerList<DegreeType>>{}));
        ARIADNE_TEST_FAIL(MultiIndexList(InitializerList<MultiIndex>{}));
        ARIADNE_TEST_FAIL(MultiIndexList(InitializerList<InitializerList<DegreeType>>{{1u,2u},{3u}}));
        ARIADNE_TEST_FAIL(MultiIndexList(InitializerList<MultiIndex>{MultiIndex({1u,2u}),MultiIndex({3u})}));

        MultiIndexList erased({a0,a1,a2});
        auto erased_end=erased.erase(erased.begin()+1);
        ARIADNE_TEST_EQUALS(erased.size(),2u);
        ARIADNE_TEST_EQUALS(erased[1],a2);
        ARIADNE_TEST_ASSERT(erased_end==erased.end());

        MultiIndexList& erased_alias=erased;
        ARIADNE_TEST_EXECUTE(erased=erased_alias);
        ARIADNE_TEST_EQUALS(erased.size(),2u);

        MultiIndexList size_mismatch({a0});
        ARIADNE_TEST_ASSERT(not (erased==size_mismatch));

        MultiIndexList argument_mismatch({MultiIndex({1u,2u})});
        MultiIndexList argument_mismatch_other({MultiIndex({1u,2u,0u})});
        ARIADNE_TEST_ASSERT(not (argument_mismatch==argument_mismatch_other));

        MultiIndexList value_mismatch({MultiIndex({1u,2u})});
        MultiIndexList value_mismatch_other({MultiIndex({1u,3u})});
        ARIADNE_TEST_ASSERT(not (value_mismatch==value_mismatch_other));

    }


};

Int main() {
  TestMultiIndex().test();
  return ARIADNE_TEST_FAILURES;
}

