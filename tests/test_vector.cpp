/***************************************************************************
 *            test_vector.cpp
 *
 *  Copyright  2006-20  Pieter Collins, Alberto Casagrande
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

#include <iostream>
#include <fstream>
#include <cassert>

#include "numeric/numeric.hpp"
#include "numeric/floats.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "interval/interval.hpp"

#include "algebra_testing_context.hpp"

using namespace std;
using namespace Ariadne;

class TestVector {
  public:
    Void test();
  private:
    Void test_range();
    Void test_slice();
    Void test_constructors();
    Void test_comparisons();
    Void test_arithmetic();
    Void test_misc();
};

Void
TestVector::test()
{
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_range());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_slice());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_constructors());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_comparisons());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_arithmetic());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_misc());
}

namespace {

constexpr bool check_concept()
{
    return requires(FloatDPApproximation ax, FloatDPBounds ix, FloatDP ex,
                    Vector<FloatDPApproximation> av, Vector<FloatDPBounds> iv, Vector<FloatDP> ev) {
        iv=Vector<FloatDPBounds>(ev);
        av=av+av;
        iv=ev+ev; iv=ev+iv; iv=iv+ev; iv=iv+iv;
        av=av-av; iv=ev-ev; iv=ev-iv; iv=iv-ev; iv=iv-iv;
        av=ax*av; iv=ex*ev; iv=ex*iv; iv=ix*ev; iv=ix*iv;
        av=av*ax; iv=ev*ex; iv=ev*ix; iv=iv*ex; iv=iv*ix;
        av=av/ax; iv=ev/ex; iv=ev/ix; iv=iv/ex; iv=iv/ix;
    };
}

static_assert(check_concept());
static_assert(Same<
    decltype(cast_exact(declval<Vector<FloatDPApproximation> const&>())),
    Vector<FloatDP>>);

} // namespace


Void
TestVector::test_range()
{
    Range prefix_range=range(4u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), prefix_range.start(), 0u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), prefix_range.stop(), 4u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), prefix_range.size(), 4u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), prefix_range.stride(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), prefix_range[2u], 2u);

    Range offset_range=range(2u,5u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), offset_range.start(), 2u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), offset_range.stop(), 5u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), offset_range.size(), 3u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), offset_range[1u], 3u);

    RangeIterator range_iterator(2u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), *range_iterator, 2u);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), range_iterator!=RangeIterator(3u));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), ++range_iterator);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), *range_iterator, 3u);

    auto first=begin(offset_range);
    auto last=end(offset_range);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), *first, 2u);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), first!=last);
}

Void
TestVector::test_slice()
{
    Slice sample_slice=Ariadne::slice(3u,2u,4u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice.size(), 3u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice.start(), 2u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice.stride(), 4u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice.stop(), 14u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice[0u], 2u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice[1u], 6u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sample_slice[2u], 10u);
}

Void
TestVector::test_constructors()
{
    Vector<RoundedFloatDP> v0(3,dp);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v0.size()==3);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v0[0]==0.0_x);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v0[1]==0.0_x);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v0[2]==0.0_x);

    Vector<RoundedFloatDP> v1({3.25_x,-0.75_x,0.0_x,1.375_x},dp);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v1.size()==4);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v1[0]==3.25_x);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v1[1]==-0.75_x);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v1[2]==0.0_x);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), v1[3]==1.375_x);

    Vector<FloatDPApproximation> va({3.25,-0.75,0.0,1.375},dp);
    Vector<FloatDPBounds> vb({3.25_x,-0.75_x,0.0_x,1.375_x},dp);
    Vector<FloatDP> vx({3.25_x,-0.75_x,0.0_x,1.375_x},dp);
}


Void
TestVector::test_comparisons()
{
    Vector<RoundedFloatDP> v(2,dp); v[0]=1.25_x; v[1]=-1.375_x;
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "v", v);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), v, Vector<RoundedFloatDP>({1.25_x,-1.375_x},dp));
    ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), v, !=,  Vector<RoundedFloatDP>({1.25_x,-1.625_x},dp));
}


Void
TestVector::test_arithmetic()
{
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), + Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp), Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), - Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp), Vector<RoundedFloatDP>({-2.0_x,+3.0_x,-5.0_x},dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp) + Vector<RoundedFloatDP>({1.25_x,2.75_x,-3.5_x},dp), Vector<RoundedFloatDP>({3.25_x,-0.25_x,1.5_x},dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp) - Vector<RoundedFloatDP>({1.25_x,2.75_x,-3.5_x},dp), Vector<RoundedFloatDP>({0.75_x,-5.75_x,8.5_x},dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), RoundedFloatDP(4.0_x,dp) * Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp), Vector<RoundedFloatDP>({8.0_x,-12.0_x,20.0_x},dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp) * RoundedFloatDP(4.0_x,dp), Vector<RoundedFloatDP>({8.0_x,-12.0_x,20.0_x},dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), Vector<RoundedFloatDP>({2.0_x,-3.0_x,5.0_x},dp) / RoundedFloatDP(4.0_x,dp), Vector<RoundedFloatDP>({0.5_x,-0.75_x,1.25_x},dp));

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sup_norm(Vector<RoundedFloatDP>({2.0_x,-3.0_x,1.0_x},dp)), RoundedFloatDP(3.0_x,dp));
}


Void
TestVector::test_misc()
{
    DoublePrecision pr;

    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), Vector<RoundedFloatDP>::unit(2u,2u,dp));

    Vector<RoundedFloatDP> invalid_at({1.0_x,2.0_x},dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), invalid_at.at(0u));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), invalid_at.at(2u));
    Vector<RoundedFloatDP> const& cinvalid_at=invalid_at;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), cinvalid_at.at(0u));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), cinvalid_at.at(2u));

    Vector<RoundedFloatDP> short_v({1.0_x},dp);
    Vector<RoundedFloatDP> long_v({1.0_x,2.0_x},dp);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), short_v+long_v);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), short_v-long_v);
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), dot(short_v,long_v));
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), dot(long_v,long_v));
    Array<FloatDPApproximation> vary(InitializerList<ExactDouble>{-4.0_x,3.0_x,1.0_x},pr);
    FloatDPApproximation x={1.5_x,pr};

    Vector<FloatDPApproximation> v0;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), cout << Vector<RoundedFloatDP>());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), cout << Vector<FloatDP>());
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v0.size()=" << v0.size());
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v0=" << flush; cout << v0);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "Constructor Vector<X>(SizeType, const X*) is unsafe and has been removed.");
    Array<FloatDPApproximation> a1(vary.begin(),vary.end());
    Vector<FloatDPApproximation> v1(a1);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v1=" << v1);
    Vector<FloatDPApproximation> v2=Vector<FloatDPApproximation>({2.375_x,4.25_x,-1.25_x},pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v2=" << v2);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "norm(v1)=" << norm(v1) << "  norm(v2)=" << norm(v2));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), norm(v1).raw(), 4);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), norm(v2).raw(), 4.25_x);

    Vector<FloatDPApproximation> v3(1,pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v3=" << v3);
    Vector<FloatDPApproximation> v4=v2;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v4=" << v4);
    Vector<FloatDPApproximation> v5({-4.0_x,3.0_x,1.0_x},pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "v5=" << v5);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), v1, v5);
    // Blank diagnostic omitted.

    Vector<FloatDPApproximation> vf0;
    v1=Vector<FloatDPApproximation>({0.25_x,-1.5_x},pr);
    v2=Vector<FloatDPApproximation>({-0.5_x,2.25_x},pr);
    vf0=-v1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), vf0 << " = -" << v1);
    vf0=Vector<FloatDPApproximation>(v1)+v2;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), vf0 << " = " << v1 << " + " << v2);
    vf0=Vector<FloatDPApproximation>(v1)-v2;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), vf0 << " = " << v1 << " - " << v2);
    vf0=x*Vector<FloatDPApproximation>(v2);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), vf0 << " = " << x << " * " << v2);
    vf0=Vector<FloatDPApproximation>(v1)*x;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), vf0 << " = " << v1 << " * " << x);
    vf0=Vector<FloatDPApproximation>(v1)/x;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), vf0 << " = " << v1 << " / " << x);
    // Blank diagnostic omitted.

    Vector< FloatDPBounds > iv1=Vector<FloatDPBounds>({{0.984375_x,1.015625_x},{2.25_x,2.375_x},{4.0_x,4.375_x},{-0.03125_x,0.015625_x}},pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv1=" << iv1);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "norm(iv1)=" << norm(iv1));
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "norm(iv1).upper()=" << norm(iv1).upper());

    Vector< FloatDPBounds > iv2=Vector<FloatDPBounds>({{-1.0_x,1.0_x},{-1.0_x,1.0_x}},pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv2=" << iv2);
    Vector< FloatDPBounds > iv3(3,pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv3=" << iv3);
    iv3=Vector<FloatDPBounds>({{4.25_x,4.25_x},{2.375_x,2.375_x}},pr);
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv3=" << iv3);
    FloatDPBounds ix=FloatDPBounds(-2,1,pr);

    Vector< FloatDPBounds > iv0;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv0=" << iv0);
    iv1=iv0;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv1=" << iv1);
    iv1=iv2;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv1=" << iv1);
    // Blank diagnostic omitted.

    FloatDPBounds ix2=iv2[0];
    FloatDPBounds ix3=iv3[0];
    FloatDPBounds ix1=ix2+ix3;
    ix1=ix2+ix3;

    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), "iv2=" << iv2 << ", iv3=" << iv3);
    iv1=iv2+iv3;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv1 << " = " << iv2 << " + " << iv3);
    iv1=iv2-iv3;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv1 << " = " << iv2 << " - " << iv3);
    iv1=ix*iv3;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv1 << " = " << ix << " * " << iv3);
    iv1=iv2*ix;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv1 << " = " << iv2 << " * " << ix);
    ix=FloatDPBounds(1,2,dp);
    iv1=iv2/ix;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv1 << " = " << iv2 << " / " << ix);
    // Blank diagnostic omitted.

    Vector<FloatDP> ev1=cast_exact(v1);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), ev1.size(), v1.size());
    for(SizeType i=0; i!=v1.size(); ++i) {
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), ev1[i], v1[i].raw());
    }
    FloatDP ex(x.raw());
    iv0=iv1+ev1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << iv1 << " + " << ev1);
    iv0=ev1+iv1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << ev1 << " + " << iv1);
    iv0=iv1-ev1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << iv1 << " - " << ev1);
    iv0=ev1-iv1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << ev1 << " - " << iv1);
    iv0=ex*iv1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << ex << " * " << iv1);
    iv0=ix*ev1;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << ix << " * " << ev1);
    iv0=iv1*ex;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << iv1 << " * " << ex);
    iv0=ev1*ix;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << ev1 << " * " << ix);
    iv0=iv1/ex;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << iv1 << " / " << ex);
    iv0=ev1/ix;
    ARIADNE_TESTING_NOTE_CTX(AlgebraTesting::context(), iv0 << " = " << ev1 << " / " << ix);

    iv0=ev1;
    iv0/=ix;
    iv0=Vector<FloatDPBounds>({2,1},pr);
    iv1=Vector<FloatDPBounds>({0,1},pr);
    /*
      ARIADNE_TEST_ASSERT( (iv0+=Vector<FloatDPBounds>("[0,1]")) == Vector<FloatDPBounds>("[2,2]") );
      ARIADNE_TEST_ASSERT( (iv0-=Vector<FloatDPBounds>("[0,1]")) == Vector<FloatDPBounds>("[2,1]") );
      ARIADNE_TEST_ASSERT( (iv0*=2) == Vector<FloatDPBounds>("[4,2]") );
      ARIADNE_TEST_ASSERT( (iv0/=4) == Vector<FloatDPBounds>("[1,0.5_x]") );
    */

    {
        RoundedFloatDP s(1,dp);
        Vector<RoundedFloatDP> v({2.0_x},dp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), join(s,v), Vector<RoundedFloatDP>({1.0_x,2.0_x},dp));
    }
    {
        FloatDPBounds s(1,dp);
        Vector<FloatDPBounds> v({2.0_x},dp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), join(s,v), Vector<FloatDPBounds>({1.0_x,2.0_x},dp));
    }
    {
        FloatDPApproximation s(1,dp);
        Vector<FloatDPApproximation> v({2.0_x},dp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), join(s,v), Vector<FloatDPApproximation>({1.0_x,2.0_x},dp));
    }
    {
        FloatDPUpperInterval s({1.0_x,1.0_x});
        Vector<FloatDPUpperInterval> v({FloatDPUpperInterval({2.0_x,2.0_x})});
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), join(s,v).size(), 2u);
    }
    {
        MultiplePrecision mp(128);
        FloatMPBounds s(1,mp);
        Vector<FloatMPBounds> v({2.0_x},mp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), join(s,v).size(), 2u);
    }
    {
        MultiplePrecision mp(128);
        FloatMPApproximation s(1,mp);
        Vector<FloatMPApproximation> v({2.0_x},mp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), join(s,v).size(), 2u);
    }

    Covector<int> ci({1,2});
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), characteristics(ci));
    Covector<FloatDPApproximation> ca(1u,dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), characteristics(ca));

    {
        Vector<FloatDPBounds> v({1.0_x,2.0_x},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), dot(v,v));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), two_norm(v));
    }
    {
        MultiplePrecision mp(128);
        Vector<FloatMPApproximation> v({1.0_x,2.0_x},mp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), dot(v,v));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), two_norm(v));
    }
    {
        MultiplePrecision mp(128);
        Vector<FloatMPBounds> v({1.0_x,2.0_x},mp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), dot(v,v));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), two_norm(v));
    }

    /*
      cout << "test_vector_slice" << endl;
      v1=Vector<FloatDPApproximation>("[-1.25_x,0.75_x,-0.5_x,-4.25_x,2.375_x]");
      cout << v1 << endl;
      VectorSlice<FloatDPApproximation> vs1(2,v1.begin()+2,2);
      cout << vs1 << endl;
      VectorSlice<FloatDPApproximation> vs2(2,v1.begin(),3);
      cout << vs2 << endl;

      iv1=vs1+vs2;
      cout << iv1 << endl;
      iv1=vs1-vs2;
      cout << iv1 << endl;
      iv1=x*vs1;
      cout << iv1 << endl;
      iv1=vs1*x;
      cout << iv1 << endl;
      iv1=vs1/x;
      cout << iv1 << endl;
    */

}


int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("vector");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

    TestVector().test();

    

    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
