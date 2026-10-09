/***************************************************************************
 *            test_differential.cpp
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>


#include "numeric/numeric.hpp"
#include "interval/interval.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/differential.hpp"

#include "algebra/expansion.tpl.hpp"

#include "utility/test.hpp"


using namespace Ariadne;
using namespace std;

template<class R, class A, class P>
Void henon(R& r, const A& x, const P& p)
{
    r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
    r[1]=x[0];
}

template<class DF>
Vector<DF>
henon(const Vector<DF>& x, const Vector<typename DF::NumericType>& p)
{
    Vector<DF> r(2,2,x.degree(),p.zero_element()); henon(r,x,p); return r;
}



template<class DF>
class TestDifferential {
    typedef typename DF::ValueType X;
    typedef typename X::PrecisionType PR;
    typedef X ScalarType;
    typedef typename DF::SeriesType SeriesType;
    typedef DF DifferentialType;
    typedef Vector<DF> DifferentialVectorType;
  private:
    PR pr;
    X c1;
    DifferentialType x1,x2,x3;
  public:
    TestDifferential()
        : c1(pr), x1(2,4,pr), x2(2,4,pr), x3(1,4,pr)
    {
        c1=ScalarType(3.0_x,pr);
        x1=DifferentialType(2,4,{{{0,0},2.0_x},{{1,0},1.0_x},{{2,0},0.5_x}},pr);
        x2=DifferentialType(2,4,{{{0,0},3.0_x},{{1,0},1.0_x},{{2,0},0.25_x}},pr);
        x3=DifferentialType(1,4,{{{0},2.0_x},{{1},1.0_x},{{2},0.125_x}},pr);

        ARIADNE_TEST_PRINT(x1);
        ARIADNE_TEST_PRINT(x2);

        ARIADNE_TEST_CALL(test_construct());
        ARIADNE_TEST_CALL(test_degree());
        ARIADNE_TEST_CALL(test_neg());
        ARIADNE_TEST_CALL(test_add());
        ARIADNE_TEST_CALL(test_sub());
        ARIADNE_TEST_CALL(test_mul());
        ARIADNE_TEST_CALL(test_div());
        ARIADNE_TEST_CALL(test_rec());
        ARIADNE_TEST_CALL(test_pow());
        ARIADNE_TEST_CALL(test_compose());
        ARIADNE_TEST_CALL(test_gradient());
        ARIADNE_TEST_CALL(test_hessian());
        ARIADNE_TEST_CALL(test_tpl_helpers());
    }

    Void test_degree() {
        ARIADNE_TEST_ASSERT(x1.degree()==4);
    }

    Void test_construct() {
        ARIADNE_TEST_CONSTRUCT(DifferentialType,x,(2,4, { {{0,0},2.0_x}, {{1,0},1.0_x}, {{2,0},0.5_x} }, pr));
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_ASSERT(x.expansion().is_sorted(GradedIndexLess()));
    }

    Void test_neg() {
        DifferentialType nx1(2,4, { {{0,0},-2.0_x}, {{1,0},-1.0_x}, {{2,0},-0.5_x} }, pr);
        ARIADNE_TEST_EQUALS(-x1,nx1);
    }

    Void test_add() {
        DifferentialType x1px2(2,4, { {{0,0},5.0_x}, {{1,0},2.0_x}, {{2,0},0.75_x} }, pr);
        ARIADNE_TEST_EQUALS(x1+x2,x1px2);
        ARIADNE_TEST_EVALUATE(x1+c1);
        ARIADNE_TEST_EVALUATE(c1+x1);
        //assert((x1+x2)==DifferentialType("[3,2,0,0]"));
    }

    Void test_sub() {
        DifferentialType x1mx2(2,4, { {{0,0},-1.0_x}, {{1,0},0.0_x}, {{2,0},0.25_x} }, pr);
        ARIADNE_TEST_EQUALS(x1-x2,x1mx2);
        ARIADNE_TEST_EVALUATE(x1-c1);
        ARIADNE_TEST_EVALUATE(c1-x1);
        //assert((x1-x2)==DifferentialType("[-1,0,0,0]"));
    }

    Void test_mul() {
        DifferentialType y1(2,2,{{{0,0},1.0_x}, {{1,0},2.0_x}, {{0,1},3.0_x}, {{2,0},4.0_x}, {{1,1},5.0_x}, {{0,2},6.0_x}},pr);
        DifferentialType y2(2,2,{{{0,0},2.0_x}, {{1,0},3.0_x}, {{0,1},5.0_x}, {{2,0},7.0_x}, {{1,1},11.0_x}, {{0,2},13.0_x}},pr);
        DifferentialType y1my2(2,2,{{{0,0},2.0_x}, {{1,0},7.0_x}, {{0,1},11.0_x}, {{2,0},21.0_x}, {{1,1},40.0_x}, {{0,2},40.0_x}},pr);
        DifferentialType cmy2(2,2,{{{0,0},10.0_x}, {{1,0},15.0_x}, {{0,1},25.0_x}, {{2,0},35.0_x}, {{1,1},55.0_x}, {{0,2},65.0_x}},pr);
        X c={5,pr};
        ARIADNE_TEST_EQUAL(y1*y2,y1my2);
        ARIADNE_TEST_EQUAL(c*y2,cmy2);
        ARIADNE_TEST_EQUAL(y2*c,cmy2);
    }

    Void test_div() {
        ARIADNE_TEST_CALL(test_rec());

        DifferentialType y2(2,2,{{{0,0},4.0_x}, {{1,0},3.0_x}, {{0,1},5.0_x}, {{2,0},7.0_x}, {{1,1},11.0_x}, {{0,2},13.0_x}},pr);

        ARIADNE_TEST_PRINT(x1);
        ARIADNE_TEST_PRINT(y2);
        ARIADNE_TEST_PRINT(c1);
        ARIADNE_TEST_EVALUATE(x1/c1);
        ARIADNE_TEST_EVALUATE(c1/x1);
        ARIADNE_TEST_EQUAL(((x1/c1)*c1),x1);
        ARIADNE_TEST_EQUAL((c1/x1)*x1,DifferentialType::constant(2,4,c1));
        ARIADNE_TEST_EVALUATE(x1/y2);
        ARIADNE_TEST_EQUAL((x1/y2)*y2,x1);
        /*
          DifferentialType x3("[2,3,4]");
          DifferentialType x4("[1,0,0]");
          cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
          assert((x3/x4)==x3);
          cout << x4 << "/" << x3 << " = " << x4/x3 << std::endl;
          assert((x4/x3)==DifferentialType("[0.5_x,-0.75_x,1.25_x]"));
          cout << 1 << "/" << x2 << " = " << 1/x2 << std::endl;
          assert((1/x2)==DifferentialType("[0.5_x,-0.25_x,0.25_x,-0.375_x]"));
          cout << x1 << "/" << x2 << " = " << x1/x2 << std::endl;
          assert((x1/x2)==DifferentialType("[0.5_x,0.25_x,-0.25_x,0.375_x]"));
        */
    }

    Void test_rec() {
        ARIADNE_TEST_CONSTRUCT(DifferentialType,y1,(2,2,{{{0,0},1.0_x}, {{1,0},2.0_x}, {{0,1},3.0_x}, {{2,0},4.0_x}, {{1,1},5.0_x}, {{0,2},6.0_x}},pr));
        ARIADNE_TEST_PRINT(Series<X>(Rec(),2*y1.value()));
        auto drec=UnivariateDifferential<X>(Rec(),y1.degree(),y1.value());
        ARIADNE_TEST_PRINT(drec);
        ARIADNE_TEST_PRINT(compose(drec,y1));

        ARIADNE_TEST_PRINT(DifferentialType::_compose(Series<X>(Rec(),y1.value()),y1));
        ARIADNE_TEST_PRINT(rec(y1));
        ARIADNE_TEST_EQUAL(rec(rec(y1)),y1);

    }

    Void test_pow() {
        cout << x2 << "^5 = " << pow(x2,5) << std::endl;
        //    assert(pow(x2,5)==DifferentialType("[32,80,160,240]"));
    }

    Void test_compose() {
        //double ax[10] = { 3.0_x, 1.0_x, 0.0_x, 0.5_x, 0.0_x, 0.0_x, 0.0_x, 0.0_x, 0.0_x, 0.0_x };
        ARIADNE_TEST_CONSTRUCT(DifferentialType,x,(2,3,{{{0,0},3.0_x},{{1,0},1.0_x},{{0,1},2.0_x},{{2,0},1.0_x},{{1,1},0.5_x},{{0,2},2.0_x}},pr));
        ARIADNE_TEST_CONSTRUCT(SeriesType,y,(3,{1.0_x,2.0_x,-3.0_x,5.0_x},pr));
        ARIADNE_TEST_CONSTRUCT(SeriesType,id,(3,{3.0_x,1.0_x,0.0_x},pr));
        ARIADNE_TEST_CONSTRUCT(DifferentialType,r,(2,3,{{{0,0},1.0_x},{{1,0},2.0_x},{{0,1},4.0_x},{{2,0},-1.0_x},{{1,1},-11.0_x},{{0,2},-8.0_x},{{3,0},-1.0_x},{{2,1},15.0_x},{{1,2},42.0_x},{{0,3},16.0_x}},pr));
        ARIADNE_TEST_EQUAL(compose(y,x),r);
        ARIADNE_TEST_EQUAL(compose(id,x),x);
    }

    Void test_gradient() {
        // Regression test based on errors in Henon evaluation.
        Vector<FloatDPApproximation> x({0.875_x,-0.125_x},pr);
        Vector< Differential<FloatDPApproximation> > dx=Differential<FloatDPApproximation>::variables(1u,x);
        Differential<FloatDPApproximation> dfx=1.5_x-dx[0]*dx[0]-0.25_x*dx[1];
        ARIADNE_TEST_PRINT(dfx);
        Covector<FloatDPApproximation> g = dfx.gradient();
        ARIADNE_TEST_PRINT(g);
        ARIADNE_TEST_EQUALS(g[0],-1.75_x);
        ARIADNE_TEST_EQUALS(g[1],-0.25_x);
    }

    Void test_hessian() {
        // Test Hessian matrix of
        FloatDPApproximation a00={1.5_x,pr}; FloatDPApproximation a01={2.5_x,pr}; FloatDPApproximation a11={3.5_x,pr};
        ExactDouble y0=0.875_x; ExactDouble y1=-1.25_x;
        Vector<FloatDPApproximation> x({y0,y1},pr);
        Vector< Differential<FloatDPApproximation> > dx=Differential<FloatDPApproximation>::variables(2u,x);
        Differential<FloatDPApproximation> dfx=a00*dx[0]*dx[0]+a01*dx[0]*dx[1]+a11*dx[1]*dx[1];
        ARIADNE_TEST_PRINT(dfx);
        Matrix<FloatDPApproximation> H = dfx.hessian();
        ARIADNE_TEST_PRINT(H);
        ARIADNE_TEST_EQUAL(H[0][1],H[1][0]);
        ARIADNE_TEST_EQUALS(H[0][0],a00*2);
        ARIADNE_TEST_EQUALS(H[0][1],a01);
        ARIADNE_TEST_EQUALS(H[1][1],a11*2);
    }

    Void test_tpl_helpers() {
        ARIADNE_TEST_EXECUTE(+x1);
        ARIADNE_TEST_ASSERT(x1!=x3);
        ARIADNE_TEST_EXECUTE(x1.coefficient_characteristics());
        auto created=x1.create();
        ARIADNE_TEST_EQUALS(created.argument_size(),x1.argument_size());
        ARIADNE_TEST_EQUALS(created.degree(),x1.degree());

        DifferentialType indexed=x1;
        ARIADNE_TEST_EXECUTE(indexed[0u]=X(7u,pr));
        DifferentialType const& cindexed=indexed;
        ARIADNE_TEST_EXECUTE(cindexed[0u]);
        ARIADNE_TEST_EXECUTE(indexed.set_gradient(1u,X(8u,pr)));

        Covector<X> g(2u,pr);
        g[0u]=X(2u,pr); g[1u]=X(-3,pr);
        auto affine1=DifferentialType::affine(2u,2u,X(4u,pr),g);
        auto affine2=DifferentialType::affine(2u,X(4u,pr),g);
        ARIADNE_TEST_PRINT(affine1);
        ARIADNE_TEST_PRINT(affine2);
        Covector<X> short_g(1u,pr);
        ARIADNE_TEST_FAIL(DifferentialType::affine(2u,2u,X(4u,pr),short_g));

        Vector<X> values({1.0_x,2.0_x},pr);
        Matrix<X> G({{1.0_x,0.0_x},{0.0_x,2.0_x}},pr);
        ARIADNE_TEST_EXECUTE(DifferentialType::affine(2u,values,G));
        ARIADNE_TEST_EXECUTE(DifferentialType::identity(2u,X(1u,pr)));
        ARIADNE_TEST_EXECUTE(DifferentialType::identity(2u,values));

        ARIADNE_TEST_FAIL(DifferentialType::constants(1u,2u,2u,values));
        ARIADNE_TEST_FAIL(DifferentialType::variables(1u,2u,2u,values));
        ARIADNE_TEST_FAIL(DifferentialType::variables(2u,3u,2u,values));

        DifferentialType low=DifferentialType::constant(1u,1u,X(-2,pr));
        DifferentialType high=DifferentialType::constant(1u,1u,X(3u,pr));
        DifferentialType zero=DifferentialType::constant(1u,1u,X(0u,pr));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(Min(),low,high));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(Min(),high,low));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(Max(),low,high));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(Max(),high,low));
        ARIADNE_TEST_FAIL(AlgebraOperations<DifferentialType>::apply(Min(),low,low));
        ARIADNE_TEST_FAIL(AlgebraOperations<DifferentialType>::apply(Max(),high,high));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(Abs(),low));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(Abs(),high));
        ARIADNE_TEST_FAIL(AlgebraOperations<DifferentialType>::apply(Abs(),zero));
        DifferentialType wrong_as=DifferentialType::constant(2u,1u,X(1u,pr));
        ARIADNE_TEST_FAIL(AlgebraOperations<DifferentialType>::apply(Min(),low,wrong_as));
        ARIADNE_TEST_FAIL(AlgebraOperations<DifferentialType>::apply(Max(),high,wrong_as));

        DifferentialType explicit_values(2u,2u,{{{0u,0u},1.0_x},{{1u,0u},2.0_x}},pr);
        ARIADNE_TEST_PRINT(explicit_values);
        ARIADNE_TEST_FAIL((void)DifferentialType(2u,1u,{{{0u},1.0_x}},pr));

        DifferentialType degree_one=DifferentialType::constant(1u,1u,X(1u,pr));
        ARIADNE_TEST_FAIL((void)degree_one.half_hessian());

        DifferentialVectorType compose_one(1u,1u,1u,pr);
        DifferentialType compose_wrong_size=DifferentialType::constant(2u,1u,X(1u,pr));
        ARIADNE_TEST_FAIL((void)compose(compose_wrong_size,compose_one));
        DifferentialType compose_degree_zero=DifferentialType::constant(1u,0u,X(1u,pr));
        ARIADNE_TEST_FAIL((void)compose(compose_degree_zero,compose_one));

        Array<DifferentialType> empty_differential_array;
        ARIADNE_TEST_FAIL(delete new DifferentialVectorType(empty_differential_array));

        Vector<X> one_value(1u,pr);
        Vector<X> two_values(2u,pr);
        Vector<X> one_parameter(1u,pr);

        DifferentialVectorType solve_bad_dimension(2u,1u,1u,pr);
        ARIADNE_TEST_FAIL((void)solve(solve_bad_dimension,two_values));
        DifferentialVectorType solve_bad_value_size(1u,1u,1u,pr);
        ARIADNE_TEST_FAIL((void)solve(solve_bad_value_size,two_values));

        DifferentialVectorType autonomous_bad_result(1u,2u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(autonomous_bad_result,two_values));
        DifferentialVectorType autonomous_bad_initial(1u,1u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(autonomous_bad_initial,two_values));

        DifferentialVectorType timed_bad_result(1u,2u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(timed_bad_result,two_values,X(0u,pr)));
        DifferentialVectorType timed_bad_arguments(1u,1u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(timed_bad_arguments,one_value,X(0u,pr)));

        DifferentialVectorType parameter_bad_result(1u,2u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(parameter_bad_result,two_values,one_parameter));
        DifferentialVectorType parameter_bad_arguments(1u,1u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(parameter_bad_arguments,one_value,one_parameter));

        DifferentialVectorType timed_parameter_bad_result(1u,3u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(timed_parameter_bad_result,two_values,X(0u,pr),one_parameter));
        DifferentialVectorType timed_parameter_bad_arguments(1u,2u,1u,pr);
        ARIADNE_TEST_FAIL((void)flow(timed_parameter_bad_arguments,one_value,X(0u,pr),one_parameter));

        DifferentialVectorType internal_df_first(1u,2u,1u,pr);
        DifferentialVectorType internal_dx_first(2u,2u,1u,pr);
        DifferentialVectorType internal_dt_first(1u,2u,1u,pr);
        ARIADNE_TEST_FAIL((void)DifferentialVectorType::_flow(internal_df_first,internal_dx_first,internal_dt_first));

        DifferentialVectorType internal_df_second(1u,1u,1u,pr);
        DifferentialVectorType internal_dx_second(1u,2u,1u,pr);
        DifferentialVectorType internal_dt_second(1u,2u,1u,pr);
        ARIADNE_TEST_FAIL((void)DifferentialVectorType::_flow(internal_df_second,internal_dx_second,internal_dt_second));

        DifferentialVectorType internal_df_third(1u,2u,1u,pr);
        DifferentialVectorType internal_dx_third(1u,2u,1u,pr);
        DifferentialVectorType internal_dt_third(1u,3u,1u,pr);
        ARIADNE_TEST_FAIL((void)DifferentialVectorType::_flow(internal_df_third,internal_dx_third,internal_dt_third));

        DifferentialVectorType internal_df_fourth(1u,2u,1u,pr);
        DifferentialVectorType internal_dx_fourth(1u,1u,1u,pr);
        DifferentialVectorType internal_dt_fourth(1u,1u,1u,pr);
        ARIADNE_TEST_FAIL((void)DifferentialVectorType::_flow(internal_df_fourth,internal_dx_fourth,internal_dt_fourth));
    }



};



template<class DF>
class TestDifferentialVector {
    typedef typename DF::ValueType X;
    typedef typename X::PrecisionType PrecisionType;
    typedef X ScalarType;
    typedef Vector<X> VectorType;
    typedef Series<X> SeriesType;
    typedef DF DifferentialType;
    typedef Vector<DF> DifferentialVectorType;
  private:
    PrecisionType pr;
    DifferentialVectorType x1,x2,x3;
  public:
    TestDifferentialVector()
        : x1(1,2,4,pr), x2(1,2,4,pr), x3(1,1,4,pr)
    {
        x1=DifferentialVectorType(1,2,4,{{ {{0,0},2.0_x},{{1,0},1.0_x},{{2,0},0.5_x} }},pr);
        x2=DifferentialVectorType(1,2,4,{{ {{0,0},3.0_x},{{1,0},1.0_x},{{2,0},0.25_x} }},pr);
        x3=DifferentialVectorType(1,2,4,{{ {{0,0},2.0_x},{{1,0},1.0_x},{{2,0},0.125_x} }},pr);

        ARIADNE_TEST_CALL(test_degree());
        ARIADNE_TEST_CALL(test_add());
        ARIADNE_TEST_CALL(test_sub());
        ARIADNE_TEST_CALL(test_mul());
        ARIADNE_TEST_CALL(test_div());
        ARIADNE_TEST_CALL(test_jacobian());
        ARIADNE_TEST_CALL(test_differentiate());
        ARIADNE_TEST_CALL(test_compose());
        ARIADNE_TEST_CALL(test_autonomous_flow());
        ARIADNE_TEST_CALL(test_time_dependent_flow());
        ARIADNE_TEST_CALL(test_linear_flow());
        ARIADNE_TEST_CALL(test_solve());
        ARIADNE_TEST_CALL(test_mapping());
        ARIADNE_TEST_CALL(test_header_helpers());
    }

    Void test_degree() {
        ARIADNE_TEST_EQUAL(x1.degree(),4);

        // Regression test to check setting of degree for null vector
        DifferentialVectorType dx(0u,2u,3u,pr);

        if(dx.argument_size()!=2u) {
            ARIADNE_TEST_WARN("Vector<Differential<X>>::argument_size() returns 0 if the vector of differentials has no elements.");
        }
        if(dx.degree()!=3u) {
            ARIADNE_TEST_WARN("Vector<Differential<X>>::degree() returns 0 if the vector of differentials has no elements.");
        }

    }

    Void test_add() {
        cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
        //assert((x1+x2)==DifferentialVectorType("[3,2,0,0]"));
    }

    Void test_sub() {
        cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
        //assert((x1-x2)==DifferentialVectorType("[-1,0,0,0]"));
    }

    Void test_mul() {
        X c={2,pr};
        cout << x1 << "*" << c << " = " << x1*c << std::endl;
        cout << c << "*" << x1 << " = " << c*x1 << std::endl;
        //assert((x1*x2)==DifferentialVectorType("[2,3,2,0]"));
    }

    Void test_div() {
        X c={2,pr};
        cout << x1 << "/" << c << " = " << x1/c << std::endl;
    }

    Void test_jacobian() {
        DifferentialVectorType dv(0u,2u,3u,pr);
        Matrix<FloatDPApproximation> J=dv.jacobian();
        ARIADNE_TEST_EQUALS(J.row_size(),dv.result_size());
        ARIADNE_TEST_EQUALS(J.column_size(),dv.argument_size());
    }

    Void test_differentiate() {
        DifferentialType x(2,2,{{{0,0},2.0_x}, {{1,0},12.0_x}, {{0,1},7.0_x}, {{2,0},33.0_x}, {{1,1},24.0_x}, {{0,2},13.0_x}},pr);
        DifferentialType y(2,3,{{{0,0},1.0_x}, {{1,0},2}, {{0,1},3}, {{2,0},6}, {{1,1},7}, {{0,2},9},{{3,0},11},{{2,1},12},{{1,2},13},{{0,3},14}},pr);
        DifferentialType z(2,4,{ {{1,0},1.0_x}, {{2,0},1},{{1,1},3}, {{3,0},2.0_x},{{2,1},3.5_x},{{1,2},9.0_x}, {{4,0},2.75_x},{{3,1},4.0_x},{{2,2},6.5_x},{{1,3},14.0_x} },pr);
        ARIADNE_TEST_EQUAL(derivative(y,0),x);
        ARIADNE_TEST_EQUAL(antiderivative(y,0),z);
        ARIADNE_TEST_PRINT(y);
        ARIADNE_TEST_PRINT(antiderivative(y,0));
        ARIADNE_TEST_PRINT(derivative(antiderivative(y,0),0));
        ARIADNE_TEST_PRINT(derivative(antiderivative(y,0),0)-y);
        ARIADNE_TEST_EQUAL(derivative(antiderivative(y,0),0),y);
    }

    Void test_compose() {
        DifferentialVectorType x(1,2,3,{ {{{0,0},3.0_x}, {{1,0},1.0_x}, {{1,1},0.125_x}, {{0,2},0.25_x}} },pr);
        DifferentialType y(1,3,{{{0},1.0_x},{{1},-1.0_x},{{2},0.5_x},{{3},-0.25_x}},pr);
        DifferentialType z(2,3,{{{0,0},1.0_x},{{1,0},-1.0_x},{{2,0},0.5_x},{{1,1},-0.125_x},{{0,2},-0.25_x},{{3,0},-0.25_x},{{2,1},0.125_x},{{1,2},0.25_x}},pr);
        DifferentialVectorType id(1,1,3,{ {{{0},3.0_x},{{1},1.0_x}} },pr);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_PRINT(y);
        ARIADNE_TEST_PRINT(z);
        ARIADNE_TEST_EQUAL(compose(y,x),z);
        ARIADNE_TEST_EQUAL(compose(id,x),x);
    }

    Void test_solve() {
        SizeType m=2;
        SizeType n=4;
        DegreeType deg=3;
        X z(pr);
        Vector<X> x0({0,0},pr);
        Vector<X> y0({0,0},pr);
        Vector<X> xy0=join(x0,y0);
        DifferentialVectorType xy=DifferentialType::variables(n,n,deg,xy0);
        DifferentialVectorType x={xy[0],xy[1]};
        DifferentialVectorType y={xy[2],xy[3]};
        DifferentialVectorType f={x[0]+3*y[0]-2*y[1]+x[1]*y[0]*y[1], x[0]-x[1]*x[1]+y[0]+2*y[1]-x[0]*x[1]+x[0]*y[0]*y[0]/2};
        DifferentialVectorType h=solve(f,y0);
        ARIADNE_TEST_EQUAL(h.result_size(),f.result_size());
        ARIADNE_TEST_EQUAL(h.argument_size(),f.argument_size()-h.result_size());
        ARIADNE_TEST_EQUAL(h.degree(),f.degree());

        x=DifferentialType::variables(n-m,n-m,deg,x0);
        DifferentialVectorType zero(m,n-m,deg,z);

        ARIADNE_TEST_EQUAL(compose(f,join(x,h)),zero);
    }

    Void test_autonomous_flow() {
        SizeType n=2;
        DegreeType deg=4;
        X z(pr);
        Vector<X> x0={z,z};
        DifferentialVectorType x=DifferentialType::variables(n,n,deg,x0);
        DifferentialVectorType f={5+2*x[0]-3*x[1]+x[0]*x[1], 2+x[1]-x[0]*x[1]+x[0]*x[0]*x[0]/2};
        DifferentialVectorType phi=flow(f,x0);
        ARIADNE_TEST_EQUAL(phi.result_size(),f.result_size());
        ARIADNE_TEST_EQUAL(phi.argument_size(),f.argument_size()+1u);
        ARIADNE_TEST_EQUAL(phi.degree(),f.degree()+1u);
        ARIADNE_TEST_EQUAL(derivative(phi,n),compose(f,phi));

        n=1;
        x0={z+1}; auto t0=z;
        x=DifferentialType::variables(n,n,deg,x0);
        f={x[0]};
        phi=flow(f,x0);
        x={DifferentialType::variable(n+1,deg+1,x0[0],0)};
        DifferentialType t=DifferentialType::variable(n+1,deg+1,t0,1);
        ARIADNE_TEST_EQUAL(phi[0],x[0]*exp(t));
    }

    Void test_time_dependent_flow() {
        SizeType n=2;
        DegreeType deg=4;
        X z(pr);
        Vector<X> x0={z,z}; X t0=z;
        DifferentialVectorType xt=DifferentialType::variables(n+1,n+1,deg,join(x0,t0));
        DifferentialVectorType x=project(xt,Range(0,n));
        DifferentialType t=xt[n];

        DifferentialVectorType f={5+2*x[0]-3*x[1]+x[0]*x[1]+t, 2+x[1]-x[0]*x[1]*t+x[0]*x[0]*x[0]/2};
        DifferentialVectorType phi=flow(f,x0,t0);
        ARIADNE_TEST_EQUAL(phi.result_size(),f.result_size());
        ARIADNE_TEST_EQUAL(phi.argument_size(),f.argument_size());
        ARIADNE_TEST_EQUAL(phi.degree(),f.degree()+1u);
        t=DifferentialType::variable(n+1,deg+1,t0,n);
        ARIADNE_TEST_EQUAL(derivative(phi,n),compose(f,join(phi,t)));
    }

    Void test_linear_flow() {
        SizeType n=1;
        DegreeType deg=4;
        X z(pr);
        Vector<X> x0={z+1}; X t0=z;
        DifferentialVectorType xt=DifferentialType::variables(n+1,n+1,deg,join(x0,t0));
        DifferentialVectorType x=project(xt,Range(0,n));
        DifferentialType t=xt[n];
        DifferentialVectorType f={2*x[0]+12*t};
        DifferentialVectorType phi=flow(f,x0,t0);
        xt=DifferentialType::variables(n+1,n+1,deg+1,join(x0,t0));
        x=static_cast<DifferentialVectorType>(project(xt,Range(0,n)));
        t=xt[n];
        ARIADNE_TEST_EQUAL(phi[0],(x[0]+3)*exp(2*t)-6*t-3);
    }

    Void test_header_helpers() {
        DifferentialType delegated(2u,2u,pr);
        ARIADNE_TEST_EQUALS(delegated.argument_size(),2u);
        ARIADNE_TEST_EQUALS(delegated.degree(),2u);

        DifferentialType positive=DifferentialType::constant(1u,2u,X(4u,pr));
        ARIADNE_TEST_EXECUTE(AlgebraOperations<DifferentialType>::apply(UnaryElementaryOperator(Sqrt()),positive));

        DifferentialVectorType v(2u,2u,2u,pr);
        v[0]=DifferentialType::variable(2u,2u,X(1u,pr),0u);
        v[1]=DifferentialType::variable(2u,2u,X(2u,pr),1u);
        DifferentialVectorType const& cv=v;
        auto vr=cv[Range(0u,1u)];
        ARIADNE_TEST_EQUALS(vr.size(),1u);
        auto chrs=cv.element_characteristics();
        ARIADNE_TEST_EQUALS(chrs.argument_size(),2u);
        ARIADNE_TEST_EQUAL(cv.get(0u),cv[0u]);
        ARIADNE_TEST_EXECUTE(v.set(0u,cv[1u]));
        ARIADNE_TEST_FAIL(v.set(v.size(),cv[0u]));
        DifferentialType wrong=DifferentialType::constant(1u,2u,X(0u,pr));
        ARIADNE_TEST_FAIL(v.set(0u,wrong));

        DifferentialType wrong_argument_size=DifferentialType::constant(1u,2u,X(0u,pr));
        DifferentialType wrong_degree=DifferentialType::constant(2u,1u,X(0u,pr));
        ARIADNE_TEST_FAIL(v[0u]=wrong_argument_size);
        ARIADNE_TEST_FAIL(v[0u]=wrong_degree);
        ARIADNE_TEST_FAIL((DifferentialType(2u,1u,{{{0u},1.0_x}},pr)));

        auto generator=[this](SizeType i) {
            return DifferentialType::constant(1u,1u,X(static_cast<Nat>(i),pr));
        };
        DifferentialVectorType generated(1u,generator);
        ARIADNE_TEST_EQUALS(generated.size(),1u);
        ARIADNE_TEST_FAIL(DifferentialVectorType(0u,generator));

        auto empty_range=cv[Range(0u,0u)];
        ARIADNE_TEST_FAIL((void)DifferentialVectorType(empty_range));
        ARIADNE_TEST_FAIL(v=empty_range);

        DifferentialVectorType source=v;
        DifferentialVectorType const& csource=source;
        auto source_range=csource[Range(0u,1u)];
        ARIADNE_TEST_EXECUTE(v=source_range);
        ARIADNE_TEST_EQUALS(v.size(),1u);

        X z(pr);
        DifferentialVectorType df={DifferentialType::constant(2u,1u,z)};
        DifferentialVectorType dx0={DifferentialType::variable(2u,1u,z,0u)};
        DifferentialVectorType da={DifferentialType::variable(2u,1u,z,1u)};
        ARIADNE_TEST_EXECUTE(v.flow(df,dx0,da));

        auto vv=v.value();
        ARIADNE_TEST_EXECUTE(v.set_value(vv));
        Vector<X> wrong_values(2u,pr);
        ARIADNE_TEST_FAIL(v.set_value(wrong_values));

        DifferentialVectorType compose_source(1u,1u,2u,pr);
        compose_source[0u]=DifferentialType::variable(1u,2u,X(1u,pr),0u);
        ARIADNE_TEST_EXECUTE(DifferentialType::_compose(compose_source,compose_source));

        DifferentialVectorType lie_source(1u,1u,2u,pr);
        lie_source[0u]=DifferentialType::variable(1u,2u,X(1u,pr),0u);
        ARIADNE_TEST_EXECUTE(lie_derivative(lie_source,lie_source));

        Vector<X> x0({0.0_x},pr);
        Vector<X> a({0.0_x},pr);
        DifferentialVectorType autonomous_param_df={DifferentialType::constant(2u,1u,X(0u,pr))};
        ARIADNE_TEST_EXECUTE(flow(autonomous_param_df,x0,a));

        X t0(0u,pr);
        DifferentialVectorType timed_param_df={DifferentialType::constant(3u,1u,X(0u,pr))};
        ARIADNE_TEST_EXECUTE(flow(timed_param_df,x0,t0,a));
    }

    Void test_mapping() {
        DifferentialVectorType x(2,2,2,pr);
        x[0][MultiIndex::unit(2,0)]=1; x[1][MultiIndex::unit(2,1)]=1;
        cout << "x=" << x << endl;
        Vector<FloatDPApproximation> p(2,dp); p[0]=1.5_x; p[1]=0.375_x;
        DifferentialVectorType hxp(2,2,2,{ {{{0,0},1.5_x}, {{0,1},-0.375_x}, {{2,0},-1.0_x}}, {{{1,0},1.0_x}} },pr);
        ARIADNE_TEST_EQUAL(henon(x,p),hxp);
    }
};

template<class X, class PR>
Void test_differential_tpl_instantiation_branches(PR const& pr) {
    using D=Differential<X>;
    X coefficient_zero(0u,pr);
    X coefficient_one(1u,pr);
    X coefficient_two(2u,pr);

    Expansion<MultiIndex,X> expansion(1u,pr);
    expansion.append(MultiIndex({0u}),coefficient_one);
    expansion.append(MultiIndex({2u}),coefficient_two);
    ARIADNE_TEST_EXECUTE((void)D(expansion,1u));

    D empty(1u,2u,coefficient_zero);
    D constant=D::constant(1u,2u,coefficient_one);
    D high(1u,2u,coefficient_zero);
    high[MultiIndex({2u})]=coefficient_one;
    D low(1u,1u,coefficient_zero);
    low[MultiIndex({0u})]=coefficient_one;

    ARIADNE_TEST_EXECUTE((void)(empty==constant));
    ARIADNE_TEST_EXECUTE((void)(constant==high));
    ARIADNE_TEST_EXECUTE((void)(high==constant));

    D const_missing(1u,2u,coefficient_zero);
    D const& const_ref=const_missing;
    ARIADNE_TEST_EXECUTE((void)const_ref[MultiIndex({1u})]);

    D linear_only(1u,2u,coefficient_zero);
    linear_only.set_gradient(0u,coefficient_one);
    ARIADNE_TEST_EXECUTE(linear_only.half_hessian());

    D cubic_only(1u,3u,coefficient_zero);
    cubic_only[MultiIndex({3u})]=coefficient_one;
    ARIADNE_TEST_EXECUTE(cubic_only.half_hessian());

    D mixed_gap(3u,2u,coefficient_zero);
    mixed_gap[MultiIndex({1u,0u,1u})]=coefficient_one;
    ARIADNE_TEST_EXECUTE(mixed_gap.half_hessian());

    ARIADNE_TEST_EXECUTE(AlgebraOperations<D>::apply(Add(),high,low));
    ARIADNE_TEST_EXECUTE(AlgebraOperations<D>::apply(Add(),low,high));
    ARIADNE_TEST_EXECUTE(AlgebraOperations<D>::apply(Sub(),high,low));
    ARIADNE_TEST_EXECUTE(AlgebraOperations<D>::apply(Sub(),low,high));
    ARIADNE_TEST_EXECUTE(AlgebraOperations<D>::apply(Mul(),high,low));
    ARIADNE_TEST_EXECUTE(AlgebraOperations<D>::apply(Mul(),low,high));

    D degree_zero=D::constant(1u,0u,coefficient_one);
    ARIADNE_TEST_EXECUTE(derivative(degree_zero,0u));

    D valid_check=D::constant(1u,1u,coefficient_one);
    ARIADNE_TEST_EXECUTE(valid_check.check());
    D empty_check(1u,1u,coefficient_zero);
    ARIADNE_TEST_EXECUTE(empty_check.check());
    D invalid_check(1u,1u,coefficient_zero);
    invalid_check.expansion().append(MultiIndex({2u}),coefficient_one);
    ARIADNE_TEST_FAIL(invalid_check.check());

    MultiIndex wrong_index({0u,0u});
    ARIADNE_TEST_FAIL((void)constant[wrong_index]);
    D const& constant_ref=constant;
    ARIADNE_TEST_FAIL((void)constant_ref[wrong_index]);

    D wrong_as=D::constant(2u,1u,coefficient_one);
    ARIADNE_TEST_FAIL(AlgebraOperations<D>::apply(Add(),low,wrong_as));
    ARIADNE_TEST_FAIL(AlgebraOperations<D>::apply(Sub(),low,wrong_as));
    ARIADNE_TEST_FAIL(AlgebraOperations<D>::apply(Mul(),low,wrong_as));
}

template<class X, class PR>
Void test_differential_tpl_compose_branch(PR const& pr) {
    X coefficient_one(1u,pr);
    Differential<X> linear_only(1u,2u,X(0u,pr));
    linear_only.set_gradient(0u,coefficient_one);
    UnivariateDifferential<X> series(2u,pr);
    series[0u]=coefficient_one;
    series[1u]=coefficient_one;
    ARIADNE_TEST_EXECUTE(compose(series,linear_only));
}

Int main() {
    {
        using X=FloatDPApproximation;
        auto x=UnivariateDifferential<X>::variable(4u,X(1,dp));
        auto y=exp(x);
        auto dy=derivative(y);

        ARIADNE_TEST_EQUAL(dy.degree(),3u);
        for(DegreeType i=0; i<=dy.degree(); ++i) {
            ARIADNE_TEST_EQUAL(dy[i],(i+1u)*y[i+1u]);
        }

        auto c=UnivariateDifferential<X>::constant(0u,X(2,dp));
        auto dc=derivative(c);
        ARIADNE_TEST_EQUAL(dc.degree(),0u);
        ARIADNE_TEST_EQUALS(dc[0],0);
    }

    MultiplePrecision mp(128);
    test_differential_tpl_instantiation_branches<RoundedFloatDP>(dp);
    test_differential_tpl_instantiation_branches<FloatDPApproximation>(dp);
    test_differential_tpl_instantiation_branches<FloatDPBounds>(dp);
    test_differential_tpl_instantiation_branches<FloatDPUpperInterval>(dp);
    test_differential_tpl_instantiation_branches<FloatMPApproximation>(mp);
    test_differential_tpl_instantiation_branches<FloatMPBounds>(mp);
    test_differential_tpl_instantiation_branches<FloatMPUpperInterval>(mp);

    test_differential_tpl_compose_branch<RoundedFloatDP>(dp);
    test_differential_tpl_compose_branch<FloatDPApproximation>(dp);
    test_differential_tpl_compose_branch<FloatDPBounds>(dp);
    test_differential_tpl_compose_branch<FloatDPUpperInterval>(dp);
    test_differential_tpl_compose_branch<FloatMPApproximation>(mp);
    test_differential_tpl_compose_branch<FloatMPBounds>(mp);

    TestDifferential< Differential<FloatDPApproximation> > tf;
    TestDifferentialVector< Differential<FloatDPApproximation> > tfv;
    return ARIADNE_TEST_FAILURES;
}
