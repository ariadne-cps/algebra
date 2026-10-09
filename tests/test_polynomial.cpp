/***************************************************************************
 *            test_polynomial.cpp
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
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/algebra.hpp"
#include "algebra/polynomial.hpp"
#include "algebra/polynomial.tpl.hpp"

#include "utility/test.hpp"

using namespace std;
using namespace Ariadne;

class TestPolynomial
{
    typedef MultiIndex MI;
    typedef MultivariatePolynomial<RoundedFloatDP> P;
  public:
    Void test();
  private:
    Void test_cleanup();
    Void test_constructors();
    Void test_indexing();
    Void test_modifiers();
    Void test_arithmetic();
    Void test_partial_evaluate();
    Void test_evaluate_horner();
    Void test_variables();
    Void test_find();
};


Void TestPolynomial::test()
{
    ARIADNE_TEST_CALL(test_cleanup())
    ARIADNE_TEST_CALL(test_constructors())
    ARIADNE_TEST_CALL(test_indexing())
    ARIADNE_TEST_CALL(test_modifiers())
    ARIADNE_TEST_CALL(test_arithmetic())
    ARIADNE_TEST_CALL(test_partial_evaluate())
    ARIADNE_TEST_CALL(test_evaluate_horner())
    ARIADNE_TEST_CALL(test_variables())

    ARIADNE_TEST_CALL(test_find())
}


namespace {

constexpr bool check_concept()
{
    using P=MultivariatePolynomial<RoundedFloatDP>;
    return requires(RoundedFloatDP x, Vector<RoundedFloatDP> v, MultiIndex a, P p, P const cp, Vector<P> pv) {
        p=P(3u,dp);
        p=P(cp);
        p=x;
        p.reserve(2u);
        p.insert(a,x);
        x=cp[a];
        p[a]=1.0_x;
        cp.argument_size();
        p.erase(p.begin());
        p.check();
        p.cleanup();
        p.clear();
        evaluate(p,v);
        compose(p,pv);
    };
}

static_assert(check_concept());

} // namespace

Void TestPolynomial::test_cleanup()
{
/*
    {
        MultiIndex a(3);
        MultiIndex b(3); ++b;
        std::vector<ValueType> v;

        for(Nat i=0; i!=20; ++i) {
            if(i%2) { v.push_back(ValueType(a,1/(1.+i)));  ++b; ++b; a=b; ++b; } else { v.push_back(ValueType(b,1/(1.+i)));}
        }
        std::sort(v.begin(),v.end());
    }
    {
        std::vector<Int> v;
        const Int element_size=3;
        const Int index_size=1;
        const Int data_size=2;
        for(Nat i=0; i!=20; ++i) {
            Int word= ( i%2 ? (3*i-2)/2 : (3*i+2)/2 );
            double value=1/(2.+i);
            v.resize(v.size()+3);
            v[v.size()-3]=word;
            reinterpret_cast<double&>(v[v.size()-2])=value;
        }
        for(Nat i=0; i!=v.size()/3; ++i) {
            Int word; double value;
            word=v[3*i];
            value=reinterpret_cast<double&>(v[3*i+1]);
        }

        typedef Expansion<ApproximateDouble>::Iterator Iterator;
        Iterator iter1(3,&*v.begin());
        Iterator iter2(3,&*v.end());
        std::sort(iter1,iter2);

        for(Nat i=0; i!=v.size()/3; ++i) {
            Int word; double value;
            word=v[3*i];
            value=reinterpret_cast<double&>(v[3*i+1]);
        }

    }
*/

    // Test to see if the cleanup/sort operations work.
    // Since these are used in the constructors, we can't use the main constructors to test this
    MultiIndex a(3);
    MultiIndex b(3); ++b;
    MultivariatePolynomial<RoundedFloatDP> p(3,dp);
    RoundedFloatDP c(1,dp);
    for(Nat i=0; i!=2; ++i) {
        if(i%2) { p.expansion().append(a,c); ++b; ++b; a=b; ++b; } else { p.expansion().append(b,c);}
        c=hlf(c);
    }
    ARIADNE_TEST_PRINT(p)
    ARIADNE_TEST_EXECUTE(p.cleanup())
    ARIADNE_TEST_PRINT(p)

    MultivariatePolynomial<FloatDP> raw_p(3,dp);
    ARIADNE_TEST_EXECUTE(raw_p.cleanup())

}

Void TestPolynomial::test_constructors()
{
    // Empty polynomial
    ARIADNE_TEST_CONSTRUCT(UnivariatePolynomial<RoundedFloatDP>,q1,(SizeOne(),dp))
    // Dense polynomial
    ARIADNE_TEST_CONSTRUCT(UnivariatePolynomial<RoundedFloatDP>,q2,({ {0,0.0_x}, {1,0.0_x},{2,5.0_x},{3,2.0_x} },dp))
    ARIADNE_TEST_EQUAL(q2[1],0.0_x)
    ARIADNE_TEST_EQUAL(q2[2],5.0_x)
    ARIADNE_TEST_EQUAL(q2[3],2.0_x)

    // Empty polynomial
    ARIADNE_TEST_CONSTRUCT(MultivariatePolynomial<RoundedFloatDP>,p1,(3,dp))
    // Dense polynomial
    ARIADNE_TEST_CONSTRUCT(MultivariatePolynomial<RoundedFloatDP>,p2,({ {{0,0,0},0.0_x}, {{1,0,0},0.0_x},{{0,1,0},0.0_x},{{0,0,1},0.0_x}, {{2,0,0},5.0_x},{{1,1,0},2.0_x},{{1,0,1},0.0_x},{{0,2,0},0.0_x},{{0,1,2},3.0_x},{{0,0,2},0.0_x} },dp))
    ARIADNE_TEST_EQUAL(p2[MultiIndex({2,0,0})],5.0_x)
    // Sparse polynomial with unordered indiced
    ARIADNE_TEST_CONSTRUCT(MultivariatePolynomial<RoundedFloatDP>,p3,({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },dp))
    ARIADNE_TEST_EQUAL(p3[MultiIndex({1,2})],5.0_x)
    ARIADNE_TEST_EQUAL(p3[MultiIndex({0,0})],2.0_x)

    // Unordered indices
    ARIADNE_TEST_EQUAL(MultivariatePolynomial<RoundedFloatDP>({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },dp), MultivariatePolynomial<RoundedFloatDP>({ {{0,0},2.0_x}, {{1,0},3.0_x}, {{0,1},11.0_x}, {{3,0},7.0_x}, {{1,2},5.0_x} },dp))
    // Repeated indices
    ARIADNE_TEST_EQUAL(MultivariatePolynomial<RoundedFloatDP>({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{1,0},7.0_x}, {{1,2},11.0_x} },dp), MultivariatePolynomial<RoundedFloatDP>({ {{0,0},2.0_x}, {{1,0},10.0_x}, {{1,2},16.0_x} },dp))

    // Exercise both cleanup outcomes in a Polynomial coefficient type closed under addition.
    UnivariatePolynomial<FloatDPBounds> exact_cancellation(SizeOne(),dp);
    exact_cancellation.expansion().append(UniIndex(1u),FloatDPBounds(1,dp));
    exact_cancellation.expansion().append(UniIndex(1u),FloatDPBounds(-1,dp));
    exact_cancellation.expansion().append(UniIndex(2u),FloatDPBounds(2,dp));
    ARIADNE_TEST_EXECUTE(exact_cancellation.cleanup())
    ARIADNE_TEST_EQUALS(exact_cancellation.number_of_terms(),1u)
    ARIADNE_TEST_EQUALS(exact_cancellation[UniIndex(2u)],FloatDPBounds(2,dp))

}

Void TestPolynomial::test_indexing()
{
    MultivariatePolynomial<RoundedFloatDP> p({ {{0,0,0},2.0_x},  {{1,0,0},3.0_x}, {{1,0,1},5.0_x}, {{2,1,0},7.0_x} },dp);
    const MultivariatePolynomial<RoundedFloatDP>& pc=p;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,0,0})],3.0_x)

    p[MultiIndex({1,0,0})]-=0.5_x;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,0,0})],2.5_x)

    p[MultiIndex({1,1,0})]=11.0_x;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,1,0})],11.0_x)

    MultivariatePolynomial<RoundedFloatDP> q(3,dp);
    q[MultiIndex({0,0,0})]=2.0_x;
    q[MultiIndex({0,1,0})]=3.0_x;
    ARIADNE_TEST_EQUALS(q.number_of_terms(),2)
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,0,0})],2.0_x)
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,1,0})],3.0_x)
    q[MultiIndex({1,0,0})]=5.0_x;
    ARIADNE_TEST_EQUALS(q[MultiIndex({1,0,0})],5.0_x)
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,1,0})],3.0_x)
    q[MultiIndex({0,0,1})]=7.0_x;
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,0,1})],7.0_x)

    // Test insert at beginning
    p.clear();
    ARIADNE_TEST_PRINT(p.expansion())
    p[MultiIndex({0,1,0})]=2.0_x;
    p[MultiIndex({0,0,1})]=3.0_x;
    p[MultiIndex({2,1,0})]=5.0_x;
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EQUALS(pc[MultiIndex({0,0,0})],0.0_x)
    ARIADNE_TEST_EQUALS(p.number_of_terms(),3)
    ARIADNE_TEST_PRINT(p.expansion())
    ARIADNE_TEST_EXECUTE(p[MultiIndex({0,0,0})]=7)
    ARIADNE_TEST_PRINT(p.expansion())
    ARIADNE_TEST_EQUALS(p.number_of_terms(),4)
    p.expansion().graded_sort();
    ARIADNE_TEST_PRINT(p.expansion())

    MultivariatePolynomial<RoundedFloatDP>::ConstIterator iter=p.begin();
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({0,0,0}))
    ARIADNE_TEST_EQUALS(iter->coefficient(),7.0_x)
    ++iter;
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({0,1,0}))
    ARIADNE_TEST_EQUALS(iter->coefficient(),2.0_x)
    ++iter;
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({0,0,1}))
    ARIADNE_TEST_EQUALS(iter->coefficient(),3.0_x)
    ++iter;
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({2,1,0}))
    ARIADNE_TEST_EQUALS(iter->coefficient(),5.0_x)
}

Void TestPolynomial::test_modifiers()
{
    RoundedFloatDP coefficient_zero(0,dp);
    RoundedFloatDP coefficient_one(1,dp);
    RoundedFloatDP coefficient_two(2,dp);
    RoundedFloatDP coefficient_three(3,dp);

    P direct_coefficients({
        {{0u,0u},coefficient_one},
        {{1u,0u},coefficient_two}
    });
    ARIADNE_TEST_EQUALS(direct_coefficients.number_of_terms(),2u)

    auto compatible_zero=direct_coefficients.create_zero();
    ARIADNE_TEST_EQUALS(compatible_zero.argument_size(),2u)
    ARIADNE_TEST_EQUALS(compatible_zero.number_of_terms(),0u)

    P const& const_direct=direct_coefficients;
    auto const_found=const_direct.find(MultiIndex({1u,0u}));
    ARIADNE_TEST_COMPARE(const_found,!=,const_direct.end())
    ARIADNE_TEST_EQUALS(const_found->coefficient(),coefficient_two)

    P edited(2u,dp);
    ARIADNE_TEST_EXECUTE(edited.reserve(4u))
    ARIADNE_TEST_EXECUTE(edited.insert(MultiIndex({0u,0u}),coefficient_one))
    ARIADNE_TEST_EXECUTE(edited.insert(MultiIndex({1u,0u}),coefficient_two))
    ARIADNE_TEST_EQUALS(edited.number_of_terms(),2u)
    auto erase_iter=edited.find(MultiIndex({0u,0u}));
    ARIADNE_TEST_COMPARE(erase_iter,!=,edited.end())
    ARIADNE_TEST_EXECUTE(edited.erase(erase_iter))
    ARIADNE_TEST_EQUALS(edited.number_of_terms(),1u)

    P derivative({
        {{0u,0u},coefficient_three},
        {{1u,0u},coefficient_two},
        {{2u,0u},coefficient_one}
    });
    ARIADNE_TEST_EXECUTE(derivative.differentiate(0u))
    ARIADNE_TEST_EQUALS(derivative[MultiIndex({0u,0u})],coefficient_two)
    ARIADNE_TEST_EQUALS(derivative[MultiIndex({1u,0u})],coefficient_two)

    P antiderivative({
        {{0u,0u},coefficient_two},
        {{1u,0u},coefficient_two}
    });
    ARIADNE_TEST_EXECUTE(antiderivative.antidifferentiate(0u))
    ARIADNE_TEST_EQUALS(antiderivative[MultiIndex({1u,0u})],coefficient_two)
    ARIADNE_TEST_EQUALS(antiderivative[MultiIndex({2u,0u})],coefficient_one)

    P truncated({
        {{0u,0u},coefficient_one},
        {{1u,0u},coefficient_two},
        {{2u,0u},coefficient_three}
    });
    truncated.expansion().append(MultiIndex({0u,1u}),coefficient_zero);
    ARIADNE_TEST_EXECUTE(truncated.truncate(1u))
    ARIADNE_TEST_EQUALS(truncated.number_of_terms(),2u)
    ARIADNE_TEST_EQUALS(truncated[MultiIndex({0u,0u})],coefficient_one)
    ARIADNE_TEST_EQUALS(truncated[MultiIndex({1u,0u})],coefficient_two)
    ARIADNE_TEST_EQUALS(truncated[MultiIndex({2u,0u})],coefficient_zero)
}

Void TestPolynomial::test_arithmetic()
{
    RoundedFloatDP coefficient_zero(0,dp);
    RoundedFloatDP coefficient_one(1,dp);
    RoundedFloatDP coefficient_two(2,dp);
    RoundedFloatDP coefficient_three(3,dp);

    P unary_source({
        {{0u,0u},coefficient_one},
        {{1u,0u},coefficient_two}
    });
    P negative_expected({
        {{0u,0u},-coefficient_one},
        {{1u,0u},-coefficient_two}
    });
    ARIADNE_TEST_EXECUTE(unary_source.check())
    ARIADNE_TEST_EQUAL(+unary_source,unary_source)
    ARIADNE_TEST_EQUAL(-unary_source,negative_expected)

    P scaled=unary_source;
    ARIADNE_TEST_EXECUTE(scaled*=coefficient_two)
    ARIADNE_TEST_EQUALS(scaled[MultiIndex({0u,0u})],coefficient_two)
    ARIADNE_TEST_EXECUTE(scaled*=coefficient_zero)
    ARIADNE_TEST_EQUALS(scaled.number_of_terms(),0u)

    P subtract_left({
        {{0u,0u},coefficient_one},
        {{1u,0u},coefficient_two},
        {{3u,0u},coefficient_three}
    });
    P subtract_right({
        {{0u,0u},coefficient_three},
        {{2u,0u},coefficient_two},
        {{3u,0u},coefficient_one}
    });
    ARIADNE_TEST_EXECUTE(subtract_left-subtract_right)
    ARIADNE_TEST_EXECUTE(subtract_right-subtract_left)

    P lone_term({{{2u,0u},coefficient_two}});
    P empty_two_variables(2u,dp);
    ARIADNE_TEST_EXECUTE(lone_term-empty_two_variables)
    ARIADNE_TEST_EXECUTE(empty_two_variables-lone_term)

    P mismatched_arguments(3u,dp);
    ARIADNE_TEST_FAIL((void)(subtract_left+mismatched_arguments))
    ARIADNE_TEST_FAIL((void)(subtract_left-mismatched_arguments))
    ARIADNE_TEST_FAIL((void)(subtract_left*mismatched_arguments))

    MultivariateMonomial<RoundedFloatDP> monomial(MultiIndex({1u,1u}),coefficient_two);
    P monomial_product({
        {{0u,0u},coefficient_one},
        {{1u,0u},coefficient_three}
    });
    ARIADNE_TEST_EXECUTE(monomial_product*=monomial)
    ARIADNE_TEST_EQUALS(monomial_product[MultiIndex({1u,1u})],coefficient_two)
    ARIADNE_TEST_EQUALS(monomial_product[MultiIndex({2u,1u})],RoundedFloatDP(6,dp))

    MultivariateMonomial<RoundedFloatDP> zero_monomial(MultiIndex({1u,0u}),coefficient_zero);
    ARIADNE_TEST_EXECUTE(monomial_product*=zero_monomial)
    ARIADNE_TEST_EQUALS(monomial_product.number_of_terms(),0u)

    ARIADNE_TEST_EQUAL(P(3,dp)+P(3,dp),P(3,dp))
    ARIADNE_TEST_EQUAL(P(3,dp)+P({ {{2,1,0},2.0_x} },dp),P({ {{2,1,0},2.0_x} },dp))
    ARIADNE_TEST_EQUAL(P(3,dp)+P({ {{2,1,0},2.0_x}, {{0,1,0},3.0_x}, {{1,1,0},5.0_x} },dp), P({ {{0,1,0},3.0_x}, {{1,1,0},5.0_x}, {{2,1,0},2.0_x} },dp))

    MultivariatePolynomial<RoundedFloatDP> x0(3,dp); x0[MultiIndex({1,0,0})]=1.0_x;
    MultivariatePolynomial<RoundedFloatDP> x1(3,dp); x1[MultiIndex({0,1,0})]=1.0_x;
    MultivariatePolynomial<RoundedFloatDP> x2=MultivariatePolynomial<RoundedFloatDP>::coordinate(3,2,dp);
    RoundedFloatDP zero(0,dp);
    UnivariatePolynomial<RoundedFloatDP> y=UnivariatePolynomial<RoundedFloatDP>::coordinate(SizeOne(),IndexZero(),dp);
    ARIADNE_TEST_EXECUTE(y=UnivariatePolynomial<RoundedFloatDP>::coordinate(SizeOne(),IndexZero(),zero))
    y=UnivariatePolynomial<RoundedFloatDP>::coordinate(dp);

    RoundedFloatDP w(3,dp);
    Vector<RoundedFloatDP> v({3,5,2},dp);
    UnivariatePolynomial<FloatDP> raw_y=UnivariatePolynomial<FloatDP>::coordinate(SizeOne(),IndexZero(),dp);
    MultivariatePolynomial<FloatDP> raw_x0=MultivariatePolynomial<FloatDP>::coordinate(3u,0u,dp);
    FloatDPBounds bw(3,dp);
    Vector<FloatDPBounds> bv({3.0_x,5.0_x,2.0_x},dp);

    ARIADNE_TEST_EQUALS(evaluate(2*x0*x0-1,v),2*v[0]*v[0]-1)
    ARIADNE_TEST_EQUALS(evaluate(2*x1*x1-1,v),2*v[1]*v[1]-1)
    ARIADNE_TEST_EXECUTE(evaluate(raw_y,bw))
    Expansion<UniIndex,FloatDP> raw_q(SizeOne(),dp);
    raw_q.append(UniIndex(2u),FloatDP(1,dp));
    raw_q.append(UniIndex(0u),FloatDP(3,dp));
    ARIADNE_TEST_EQUALS(horner_evaluate(raw_q,bw),bw*bw+FloatDPBounds(3,dp))
    ARIADNE_TEST_EXECUTE(evaluate(raw_x0,bv))
    /* Failing with UnivariatePolynomial
    ARIADNE_TEST_EQUALS(evaluate(2*y*y-1,w),2*w*w-1);
    ARIADNE_TEST_EQUALS(evaluate(8*y*y*(y*y-1)+1,w),8*w*w*(w*w-1)+1);
    ARIADNE_TEST_EQUALS(compose(2*y*y-1,x0),2*x0*x0-1);
     */
}

Void TestPolynomial::test_partial_evaluate()
{
    RoundedFloatDP coefficient_zero(0,dp);
    RoundedFloatDP coefficient_one(1,dp);
    RoundedFloatDP coefficient_two(2,dp);

    P source({
        {{0u,0u,0u},RoundedFloatDP(2,dp)},
        {{0u,1u,0u},RoundedFloatDP(3,dp)},
        {{0u,2u,0u},RoundedFloatDP(5,dp)},
        {{1u,3u,0u},RoundedFloatDP(7,dp)},
        {{0u,0u,1u},RoundedFloatDP(11,dp)}
    });

    auto at_zero=partial_evaluate(source,1u,coefficient_zero);
    ARIADNE_TEST_EQUALS(at_zero.argument_size(),2u)
    ARIADNE_TEST_EQUALS(at_zero[MultiIndex({0u,0u})],RoundedFloatDP(2,dp))
    ARIADNE_TEST_EQUALS(at_zero[MultiIndex({0u,1u})],RoundedFloatDP(11,dp))

    auto at_one=partial_evaluate(source,1u,coefficient_one);
    ARIADNE_TEST_EQUALS(at_one[MultiIndex({0u,0u})],RoundedFloatDP(10,dp))
    ARIADNE_TEST_EQUALS(at_one[MultiIndex({1u,0u})],RoundedFloatDP(7,dp))
    ARIADNE_TEST_EQUALS(at_one[MultiIndex({0u,1u})],RoundedFloatDP(11,dp))

    auto at_two=partial_evaluate(source,1u,coefficient_two);
    ARIADNE_TEST_EQUALS(at_two[MultiIndex({0u,0u})],RoundedFloatDP(28,dp))
    ARIADNE_TEST_EQUALS(at_two[MultiIndex({1u,0u})],RoundedFloatDP(56,dp))
    ARIADNE_TEST_EQUALS(at_two[MultiIndex({0u,1u})],RoundedFloatDP(11,dp))

    P linear({
        {{0u,0u},coefficient_one},
        {{0u,1u},coefficient_two}
    });
    auto linear_at_two=partial_evaluate(linear,1u,coefficient_two);
    ARIADNE_TEST_EQUALS(linear_at_two[MultiIndex({0u})],RoundedFloatDP(5,dp))

    // Degree-zero polynomials have no entry for the first power of c.
    P constant=P::constant(2u,coefficient_two);
    P empty(2u,dp);
    for(SizeType k=0; k!=2u; ++k) {
        auto constant_at_two=partial_evaluate(constant,k,coefficient_two);
        ARIADNE_TEST_EQUALS(constant_at_two.argument_size(),1u)
        ARIADNE_TEST_EQUALS(constant_at_two.number_of_terms(),1u)
        ARIADNE_TEST_EQUALS(constant_at_two.value(),coefficient_two)
        auto empty_at_two=partial_evaluate(empty,k,coefficient_two);
        ARIADNE_TEST_EQUALS(empty_at_two.argument_size(),1u)
        ARIADNE_TEST_EQUALS(empty_at_two.number_of_terms(),0u)
    }

    // Unsupported univariate partial evaluation must report an error without
    // terminating the process or modifying the source polynomial.
    auto univariate=UnivariatePolynomial<RoundedFloatDP>::coordinate(dp);
    ARIADNE_TEST_THROWS(partial_evaluate(univariate,0u,coefficient_two),std::logic_error)
    ARIADNE_TEST_EQUALS(univariate.degree(),1u)
    ARIADNE_TEST_EQUALS(univariate[UniIndex(1u)],coefficient_one)

    ARIADNE_TEST_FAIL((void)partial_evaluate(source,3u,coefficient_one))
}

Void TestPolynomial::test_evaluate_horner()
{
    using X=RoundedFloatDP;
    using HornerPolynomial=MultivariatePolynomial<X>;

    Vector<X> v({2.0_x,3.0_x,5.0_x},dp);

    HornerPolynomial empty(3u,dp);
    ARIADNE_TEST_EQUALS(evaluate(empty,v),X(0u,dp))

    HornerPolynomial p({ {{2,0,2},1.0_x}, {{1,1,2},2.0_x}, {{0,2,1},3.0_x}, {{1,0,1},4.0_x},
          {{0,1,0},5.0_x}, {{2,0,0},6.0_x}, {{0,0,0},7.0_x} },dp);
    X expected =
        X(1u,dp)*v[0]*v[0]*v[2]*v[2] +
        X(2u,dp)*v[0]*v[1]*v[2]*v[2] +
        X(3u,dp)*v[1]*v[1]*v[2] +
        X(4u,dp)*v[0]*v[2] +
        X(5u,dp)*v[1] +
        X(6u,dp)*v[0]*v[0] +
        X(7u,dp);
    ARIADNE_TEST_EQUALS(evaluate(p,v),expected)

    Vector<X> short_v({2.0_x,3.0_x},dp);
    ARIADNE_TEST_FAIL((void)evaluate(p,short_v))

    Expansion<MultiIndex,X> register_expansion(3u,dp);
    register_expansion.append(MultiIndex({2u,1u,2u}),X(1u,dp));
    register_expansion.append(MultiIndex({1u,1u,2u}),X(2u,dp));
    register_expansion.append(MultiIndex({0u,0u,1u}),X(3u,dp));
    X register_expected =
        v[0]*v[0]*v[1]*v[2]*v[2] +
        X(2u,dp)*v[0]*v[1]*v[2]*v[2] +
        X(3u,dp)*v[2];
    ARIADNE_TEST_EQUALS(horner_evaluate(register_expansion,v),register_expected)

    // With a single term, all powers are applied in the final Horner loop.
    Expansion<MultiIndex,X> monomial(3u,dp);
    monomial.append(MultiIndex({2u,1u,2u}),X(3u,dp));
    ARIADNE_TEST_EQUALS(horner_evaluate(monomial,v),X(900u,dp))

    Expansion<MultiIndex,X> unordered(2u,dp);
    unordered.append(MultiIndex({1u,0u}),X(1u,dp));
    unordered.append(MultiIndex({0u,1u}),X(1u,dp));
    Vector<X> uv({2.0_x,3.0_x},dp);
    ARIADNE_TEST_FAIL((void)horner_evaluate(unordered,uv))

    Expansion<UniIndex,X> univariate_empty(SizeOne(),dp);
    ARIADNE_TEST_EQUALS(horner_evaluate(univariate_empty,X(2u,dp)),X(0u,dp))

    Expansion<UniIndex,X> univariate_ordered(SizeOne(),dp);
    univariate_ordered.append(UniIndex(2u),X(1u,dp));
    univariate_ordered.append(UniIndex(0u),X(3u,dp));
    ARIADNE_TEST_EQUALS(horner_evaluate(univariate_ordered,X(2u,dp)),X(7u,dp))

    Expansion<UniIndex,X> univariate_unordered(SizeOne(),dp);
    univariate_unordered.append(UniIndex(0u),X(1u,dp));
    univariate_unordered.append(UniIndex(1u),X(1u,dp));
    ARIADNE_TEST_FAIL((void)horner_evaluate(univariate_unordered,X(2u,dp)))

    UnivariatePolynomial<X> u=UnivariatePolynomial<X>::coordinate(SizeOne(),IndexZero(),dp);
    auto univariate_polynomial_zero=horner_evaluate(univariate_empty,u);
    ARIADNE_TEST_EQUALS(univariate_polynomial_zero.number_of_terms(),0u)
    auto univariate_polynomial_value=horner_evaluate(univariate_ordered,u);
    ARIADNE_TEST_EQUAL(univariate_polynomial_value,u*u+X(3u,dp))
    ARIADNE_TEST_EQUALS(univariate_polynomial_value.degree(),2u)

    HornerPolynomial m0=HornerPolynomial::coordinate(2u,0u,dp);
    auto multivariate_polynomial_zero=horner_evaluate(univariate_empty,m0);
    ARIADNE_TEST_EQUALS(multivariate_polynomial_zero.number_of_terms(),0u)
    auto multivariate_polynomial_value=horner_evaluate(univariate_ordered,m0);
    ARIADNE_TEST_EQUAL(multivariate_polynomial_value,m0*m0+X(3u,dp))

    Vector<UnivariatePolynomial<X>> uq({u,u});
    MultivariatePolynomial<X> q2({ {{1,0},1.0_x}, {{0,1},1.0_x}, {{0,0},1.0_x} },dp);
    ARIADNE_TEST_EXECUTE(compose(q2,uq))

    HornerPolynomial m1=HornerPolynomial::coordinate(2u,1u,dp);
    Vector<HornerPolynomial> mq({m0,m1});
    ARIADNE_TEST_EXECUTE(compose(q2,mq))
}

Void TestPolynomial::test_variables()
{
    ARIADNE_TEST_EXECUTE(UnivariatePolynomial<RoundedFloatDP>::coordinates(SizeOne(),dp))
    ARIADNE_TEST_FAIL((void)MultivariatePolynomial<RoundedFloatDP>::coordinate(2u,2u,dp))

    Vector< MultivariatePolynomial<RoundedFloatDP> > x=MultivariatePolynomial<RoundedFloatDP>::variables(3,dp);
    Array< Vector<RoundedFloatDP> > e=Vector<RoundedFloatDP>::basis(2,dp);

    MultivariatePolynomial<RoundedFloatDP> p1=x[1]*3.0_x;
    MultivariatePolynomial<RoundedFloatDP> p2=p1+x[0]; p2=x[1]*3; p2=x[0]+0;
    MultivariatePolynomial<RoundedFloatDP> p3=x[0]*p2; p3=x[0]*(x[1]*3.0_x+x[0]);
    MultivariatePolynomial<RoundedFloatDP> p4=x[1]*x[2];
    MultivariatePolynomial<RoundedFloatDP> p5=p3+p4;
    MultivariatePolynomial<RoundedFloatDP> p=x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2];

    ARIADNE_TEST_EQUAL(x[0], MultivariatePolynomial<RoundedFloatDP>({ {{1,0,0},1.0_x} },dp))
    ARIADNE_TEST_EQUAL(x[1], MultivariatePolynomial<RoundedFloatDP>({ {{0,1,0},1.0_x} },dp))
    ARIADNE_TEST_EQUAL(x[2], MultivariatePolynomial<RoundedFloatDP>({ {{0,0,1},1.0_x} },dp))
    ARIADNE_TEST_EQUAL(x[0]+x[1], MultivariatePolynomial<RoundedFloatDP>({ {{1,0,0},1.0_x}, {{0,1,0},1.0_x} },dp))
    ARIADNE_TEST_EQUAL(x[0]*x[1], MultivariatePolynomial<RoundedFloatDP>({ {{1,1,0},1.0_x} },dp))
    ARIADNE_TEST_EVALUATE(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2])
    ARIADNE_TEST_EQUAL((x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]), MultivariatePolynomial<RoundedFloatDP>({ {{1,1,0},3.0_x}, {{2,0,0},1.0_x}, {{0,1,1},1.0_x} },dp))
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[1], MultivariatePolynomial<RoundedFloatDP>({ {{1,1,0},3.0_x}, {{2,0,0},1.0_x}, {{0,1,1},1.0_x} },dp))
    ARIADNE_TEST_PRINT((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[0])
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[0], MultivariatePolynomial<RoundedFloatDP>(3,dp))
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[0], MultivariatePolynomial<RoundedFloatDP>({ {{3,0,0},0.0_x} },dp))
}

Void TestPolynomial::test_find()
{
    MultivariatePolynomial<RoundedFloatDP> p({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },dp);
    MultiIndex a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TEST_PRINT(p)
    ARIADNE_TEST_PRINT(p.find(a)-p.begin())
    ARIADNE_TEST_COMPARE(p.find(a),!=,p.end())
    ARIADNE_TEST_EQUAL(p.find(a)->index(),a)
    ARIADNE_TEST_EQUAL(p.find(a)->coefficient(),5.0_x)
    a[1]=1;
    ARIADNE_TEST_EQUAL(p.find(a),p.end())
}

Int main() {
    TestPolynomial().test();
    return ARIADNE_TEST_FAILURES;
}
