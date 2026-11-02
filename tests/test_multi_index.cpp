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
#include <sstream>

#include "algebra/multi_index.hpp"

#include "algebra_testing_context.hpp"

using namespace std;
using namespace Ariadne;


class TestMultiIndex
{
  public:
    TestMultiIndex() { }

    Void test() {
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_constructor());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_uni_index());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_inline_helpers());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_preconditions());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_comparison());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_addition());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_increment());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_list());
    }

    Void test_constructor() {
        DegreeType a1[3]={2,0,3};
        MultiIndex i1(3u,a1);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i1.number_of_variables(), 3);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i1.degree(), 5);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i1[0], 2);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i1[1], 0);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i1[2], 3);

        DegreeType a2[5]={2,0,3,0,5};
        MultiIndex i2(5,a2);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i2.number_of_variables(), 5);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i2.degree(), 10);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i2[3], 0);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i2[4], 5);

        MultiIndex i3=MultiIndex::unit(9u,2u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i3.number_of_variables(), 9);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i3.degree(), 1);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i3[2], 1);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), (Int)i3[8], 0);
    }

    Void test_uni_index() {
        UniIndex a(2u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a.degree(), DegreeType(2u));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), a.size());
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a[IndexZero()], DegreeType(2u));

        a[IndexZero()]=3u;
        const UniIndex& ca=a;
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), ca[IndexZero()], DegreeType(3u));

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), ++a);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a.degree(), DegreeType(4u));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), --a);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a.degree(), DegreeType(3u));

        UniIndex b(2u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a-b, static_cast<short int>(1));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), a+=b);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a.degree(), DegreeType(5u));
    }

    Void test_inline_helpers() {
        MultiIndex a({1u,2u});
        MultiIndex b({3u,4u});

        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), a.begin()!=a.end());
        const MultiIndex& ca=a;
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), ca.begin()!=ca.end());

        MultiIndex assigned(2u);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), assigned.assign(a));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), assigned, a);

        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), DegreeType(2u)*a, MultiIndex({2u,4u}));

        GradedLess graded_less_object;
        LexicographicLess lexicographic_less_object;
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), graded_less_object(DegreeType(1u),DegreeType(2u)));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lexicographic_less_object(DegreeType(1u),DegreeType(2u)));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lexicographic_less_object(a,b));

        // Pointer proxy streams expose the address of the referenced index storage.
        MultiIndexReference mutable_reference(a);
        MultiIndexConstReference const_reference(mutable_reference);
        MultiIndexPointer mutable_pointer(&mutable_reference);
        MultiIndexConstPointer const_pointer(&const_reference);
        std::ostringstream mutable_stream, const_stream, expected_mutable, expected_const;
        mutable_stream << mutable_pointer;
        const_stream << const_pointer;
        expected_mutable << a.begin();
        expected_const << static_cast<MultiIndex const&>(a).begin();
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), mutable_stream.str(), expected_mutable.str());
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), const_stream.str(), expected_const.str());
    }

    Void test_preconditions() {
        MultiIndex a({1u,2u});
        const MultiIndex& ca=a;
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), a[a.size()]);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), ca[ca.size()]);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), ca.get(ca.size()));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), a.set(a.size(),0u));

        MultiIndex short_index({1u});
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !(a==short_index));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), a.assign(short_index));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), graded_less(a,short_index));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), lexicographic_less(a,short_index));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), reverse_lexicographic_less(a,short_index));

        MultiIndex long_index({2u,0u});
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), long_index+=short_index);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), long_index-=short_index);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultiIndex({0u,0u})-=MultiIndex({1u,0u}));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), swap(long_index,short_index));

        MultiIndexReference long_reference(long_index);
        MultiIndexReference short_reference(short_index);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), long_reference=short_index);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), long_reference+=short_index);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), swap(long_reference,short_reference));

        MultiIndex same_size({3u,4u});
        MultiIndex& same_alias=same_size;
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), same_size=same_alias);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), same_size.resize(same_size.size()));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), same_size.resize(3u));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), same_size.size(), 3u);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), same_size.clear());
    }

    Void test_comparison() {
        DegreeType a1p[3]={2,0,3};
        DegreeType a2p[3]={1,0,4};
        DegreeType a3p[3]={3,0,1};
        MultiIndex a1(3,a1p);
        MultiIndex a2(3,a2p);
        MultiIndex a3(3,a3p);

        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), graded_less(a1,a2));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), graded_less(a3,a1));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), graded_less(a3,a2));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !graded_less(a2,a1));

        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lexicographic_less(a2,a1));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lexicographic_less(a1,a3));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lexicographic_less(a2,a3));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), !lexicographic_less(a1,a1));

        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), reverse_lexicographic_less(a2,a1));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), reverse_lexicographic_less(a1,a3));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), reverse_lexicographic_less(a2,a3));

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), swap(a1,a2));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1, MultiIndex(3,a2p));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a2, MultiIndex(3,a1p));

        MultiIndexReference a1r(a1);
        MultiIndexReference a2r(a2);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), swap(a1r,a2r));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1, MultiIndex(3,a1p));

    }

    Void test_addition() {
        MultiIndex a1({2,0,3});
        MultiIndexReference a1r(a1);
        MultiIndex a2({1,1,0});
        MultiIndex a1pa2({3,1,3});
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1+a2, a1pa2);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1pa2-a1, a2);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1*2u, MultiIndex({4,0,6}));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1+=a2, a1pa2);

        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1r, a1);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a1r+=a2, a1pa2+a2);
    }

    Void test_increment() {
        MultiIndex empty;
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), ++empty);

        MultiIndex u(1);
        for(SizeType n=0; n!=8; ++n) {
            ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), u.degree(), n);
            ++u;
        }

        MultiIndex a(4);
        SizeType n=0;
        while(a.degree()<=5) {
            MultiIndex b=a; ++a; ++n;
            ARIADNE_TESTING_EXPECT_BINARY_PREDICATE_CTX(AlgebraTesting::context(), graded_less, b, a);
            MultiIndex::IndexType d=0; for(SizeType i=0; i!=a.size(); ++i) { d+=a[i]; }
            ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), a.degree(), d);
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

        MultiIndexList lst(0u,MultiIndex(3u));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.capacity(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.argument_size(), 3u);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.append(a0));
        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "*lst.begin()", *lst.begin());
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.begin()->size(), 3u);
        const MultiIndexList& const_lst=lst;
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), const_lst.begin()->size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), *lst.begin(), a0);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.append(a1));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), *lst.begin(), a0);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), *++lst.begin(), a1);

        const MultiIndexList& clst=lst;
        auto cit=clst.begin();
        auto cend=clst.end();
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), cit.operator!=(cend));
        auto cit_copy=cit;
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), cit_copy==cit);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), ++cit_copy);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), --cit_copy);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), *(cit+1), a1);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), clst.front(), a0);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), clst.back(), a1);

        MultiIndexList lstc(lst);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lstc, lst);

        MultiIndexList lstm(std::move(lst));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lstm, lstc);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst = lstc);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst, lstc);
        MultiIndex& lst_alias=lst[0u];
        (void)lst_alias;
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst=lst);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lstm = std::move(lstc));
        MultiIndexList* lstm_alias=&lstm;
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lstm=std::move(*lstm_alias));

        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "lst", lst);
        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "lstm", lstm);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lstm, lst);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.append(a2));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.append(a3));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.back(), a3);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.capacity(), 4u);
        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "&lst.front()[0]", &lst.front()[0]);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.append(a4));
        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "&lst.front()[0]", &lst.front()[0]);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.capacity(), 8u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.size(), 5u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), *--lst.end(), a4);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.back(), a4);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.resize(3u));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.capacity(), 8u);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lst.resize(11u));
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), lst.capacity()==16u or lst.capacity()==11u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.size(), 11u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst[2], a2);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst[3], az);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst.back(), az);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), lst[lst.size()]);
        const MultiIndexList& const_resized_lst=lst;
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), const_resized_lst[const_resized_lst.size()]);
        ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "lst", lst);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), *(lst.begin()+3)=a3);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), lst[3], a3);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), *(lst.begin()+3), a3);

        MultiIndexList summed(3u);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), summed.append_sum(a0,a1));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), summed.front(), a0+a1);

        MultiIndexList repeated(2u,a0);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), repeated.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), repeated[0], a0);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), repeated[1], a0);

        MultiIndexList from_multi({a0,a1,a2});
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), from_multi.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), from_multi[1], a1);

        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultiIndexList(InitializerList<InitializerList<DegreeType>>{}));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultiIndexList(InitializerList<MultiIndex>{}));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultiIndexList(InitializerList<InitializerList<DegreeType>>{{1u,2u},{3u}}));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), MultiIndexList(InitializerList<MultiIndex>{MultiIndex({1u,2u}),MultiIndex({3u})}));

        MultiIndexList erased({a0,a1,a2});
        auto erased_end=erased.erase(erased.begin()+1);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), erased.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), erased[1], a2);
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), erased_end==erased.end());

        MultiIndexList& erased_alias=erased;
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), erased=erased_alias);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(AlgebraTesting::context(), erased.size(), 2u);

        MultiIndexList size_mismatch({a0});
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not (erased==size_mismatch));

        MultiIndexList argument_mismatch({MultiIndex({1u,2u})});
        MultiIndexList argument_mismatch_other({MultiIndex({1u,2u,0u})});
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not (argument_mismatch==argument_mismatch_other));

        MultiIndexList value_mismatch({MultiIndex({1u,2u})});
        MultiIndexList value_mismatch_other({MultiIndex({1u,3u})});
        ARIADNE_TESTING_PROVE_CTX(AlgebraTesting::context(), not (value_mismatch==value_mismatch_other));

    }


};

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("multi_index");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

  TestMultiIndex().test();
  

    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}

