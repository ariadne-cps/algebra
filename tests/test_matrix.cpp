/***************************************************************************
 *            test_matrix.cpp
 *
 *  Copyright  2006-20  Pieter Collins, Alberto Casagrande
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

#define NO_CBLAS

#include <iostream>
#include <fstream>
#include <sstream>


#include "algebra_testing_context.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/matrix.tpl.hpp"
#include "algebra/symmetric_matrix.hpp"
#include "interval/interval.hpp"
#include "algebra/covector.hpp"

namespace Ariadne {
typedef Matrix<FloatDPApproximation> FloatApproximationMatrix;
}

using namespace std;
using namespace Ariadne;


class TestMatrix {
    DoublePrecision pr;
  public:
    Void test();
  private:
    Void test_project();
    Void test_data_access();
    Void test_factorisations();
    Void test_misc();
};

Void
TestMatrix::test()
{
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_project());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_data_access());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_factorisations());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_misc());
}

namespace {

constexpr bool check_concept()
{
    return requires(FloatDPApproximation fx, FloatDPBounds ix, FloatDP ex,
                    Vector<FloatDPApproximation> fv, Vector<FloatDPBounds> iv, Vector<FloatDP> ev,
                    Matrix<FloatDPApproximation> fA, Matrix<FloatDPBounds> iA, Matrix<FloatDP> eA) {
        fv=fv+fv; iv=ev+ev; iv=ev+iv; iv=iv+ev; iv=iv+iv;
        fv=fv-fv; iv=ev-ev; iv=ev-iv; iv=iv-ev; iv=iv-iv;
        fv=fx*fv; iv=ex*ev; iv=ex*iv; iv=ix*ev; iv=ix*iv;
        fv=fv*fx; iv=ev*ex; iv=ev*ix; iv=iv*ex; iv=iv*ix;
        fv=fv/fx; iv=ev/ex; iv=ev/ix; iv=iv/ex; iv=iv/ix;
        fA=fA+fA; iA=eA+eA; iA=eA+iA; iA=iA+eA; iA=iA+iA;
        fA=fA-fA; iA=eA-eA; iA=eA-iA; iA=iA-eA; iA=iA-iA;
        fA=fx*fA; iA=ex*eA; iA=ex*iA; iA=ix*eA; iA=ix*iA;
        fA=fA*fx; iA=eA*ex; iA=eA*ix; iA=iA*ex; iA=iA*ix;
        fv=fA*fv; iv=eA*ev; iv=eA*iv; iv=iA*ev; iv=iA*iv;
        fA=fA*fA; iA=eA*eA; iA=eA*iA; iA=iA*eA; iA=iA*iA;
    };
}

static_assert(check_concept());

} // namespace

Void
TestMatrix::test_project()
{
    typedef FloatDP X;
    Matrix<X> A({{11,12,13,14,15},{21,22,23,24,25},{31,32,33,34,35}},pr);
    Matrix<X> B({{120,130,140},{220,230,240}},pr);

    Vector<X> v = A[range(0,2)][1];
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), v[1], A[1][1]);
    Covector<X> u = A[2][range(1,4)];
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), u[2], A[2][3]);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A[range(0,2)][range(1,4)].row_size(), 2);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A[range(0,2)][range(1,4)].column_size(), 3);
    MatrixRange<Matrix<X>> AR(A[range(0,2)][range(1,4)]);
    Matrix<X> ACR(A[range(0,2)][range(1,4)]);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), AR[1][2], A[1][3]);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), ACR[1][2], A[1][3]);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), AR=B);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A", A);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A[range(0,2)][range(1,4)], B);

    const Matrix<X> original=A;
    Matrix<X> wrong_rows(1u,3u,pr);
    Matrix<X> wrong_columns(2u,2u,pr);
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, AR=wrong_rows);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A, original);
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, AR=wrong_columns);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A, original);
//    ARIADNE_TEST_EQUALS(AR,B);
//    ARIADNE_TEST_EQUALS(ACR,B);

}


Void
TestMatrix::test_data_access()
{
    Matrix<FloatDP> A({{1,2,3},{4,5,6}},pr);
    const Matrix<FloatDP>& C=A;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A.at(1u,2u), FloatDP(6,pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), C.at(1u,2u), FloatDP(6,pr));

    // Keep the other index valid to exercise each part of the bounds check.
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, A.at(A.row_size(),0u));
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, A.at(0u,A.column_size()));
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, C.at(C.row_size(),0u));
    ARIADNE_TESTING_EXPECT_THROWS_CTX(AlgebraTesting::context(), std::runtime_error, C.at(0u,C.column_size()));
}


Void
TestMatrix::test_factorisations()
{
    {
        PivotMatrix P=PivotMatrix::identity(2u);
        PivotMatrix const& CP=P;
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), CP.size(), 2u);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), CP.pivot(0u), 0u);
        std::ostringstream oss;
        oss << P;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), oss.str().find("PivotMatrix(")==0u);
    }

    {
        Matrix<FloatDP> A({{1,0},{0,1}},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), inverse(A));
    }

    {
        Matrix<FloatDPBounds> A({{1,0},{0,1}},dp);
        Matrix<FloatDPBounds> B({{2,0},{0,3}},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), inverse(A));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), solve(A,B));

        Matrix<FloatDPBounds> singular({{1,1},{1,1}},dp);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), inverse(singular));
    }

    {
        Matrix<Rational> A({{Rational(1),Rational(2)},{Rational(3),Rational(4)}});
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), midpoint(A), A);
    }

    // Set all diagonal elements to 0.
    auto set_diagonal_to_zero = [](FloatApproximationMatrix A){ for (SizeType i=0; i!=min(A.row_size(),A.column_size()); ++i) { A[i][i]=0; } return A; };

    {
        // LU test_factorisation
        auto A = FloatApproximationMatrix({{1.0_x,2.0_x,3.0_x},{4.0_x,5.0_x,6.0_x},{7.0_x,8.0_x,10.0_x}},pr);
        auto PLU = triangular_decomposition(A); auto P=std::get<0>(PLU); auto L=std::get<1>(PLU); auto U=std::get<2>(PLU);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(L*U-P*A), <,  1e-12);
    }

    {
        // QR decomposition
        auto I = FloatApproximationMatrix::identity(3,dp);
        auto A = FloatApproximationMatrix({{1.0_x,2.0_x,3.0_x},{4.0_x,5.0_x,6.0_x},{7.0_x,8.0_x,10.0_x}},pr);
        auto QR = gram_schmidt_orthogonalisation(A); auto Q=std::get<0>(QR); auto R=std::get<1>(QR);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(Q*R-A), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(transpose(Q)*Q-I), <,  1e-12);
    }

    {
        // Orthogonal decomposition
        auto A1 = FloatApproximationMatrix({{1,7,4},{4,-2,5},{-5,-5,-3}},pr);
        auto OR1 = orthogonal_decomposition(A1); auto O1=std::get<0>(OR1); auto R1=std::get<1>(OR1);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(O1*R1-A1), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(set_diagonal_to_zero(transpose(O1)*O1)), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(R1), <=,  1.0+1e-12);

        auto A2 = FloatApproximationMatrix({{9,-5,6},{-7,-1,-1},{8,-5,5},{-5,5,2}},pr);
        auto OR2 = orthogonal_decomposition(A2); auto O2=std::get<0>(OR2); auto R2=std::get<1>(OR2);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A2.column_size(), O2.column_size());
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(O2*R2-A2), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(set_diagonal_to_zero(transpose(O2)*O2)), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(R2), <=,  1.0+1e-12);

        auto A3 = FloatApproximationMatrix({{-2,-7,1,5},{-6,-5,-7,7},{3,2,0,8}},dp);
        auto OR3 = orthogonal_decomposition(A3); auto O3=std::get<0>(OR3); auto R3=std::get<1>(OR3);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A3.row_size(), O3.column_size());
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(O3*R3-A3), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(set_diagonal_to_zero(transpose(O3)*O3)), <,  1e-12);
        ARIADNE_ALGEBRA_DECIDED_REL_CTX(AlgebraTesting::context(), norm(R3), <=,  1.0+1e-12);
    }

}

Void
TestMatrix::test_misc()
{
    // The initializer-list constructor must reuse Matrix's shape checks
    // and reject nonsymmetric input without changing stored values.
    {
        SymmetricMatrix<Rational> symmetric(
            {{Rational(1),Rational(2)},{Rational(2),Rational(3)}});
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(
            AlgebraTesting::context(), symmetric[0u][1u], Rational(2));
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(
            AlgebraTesting::context(), symmetric[1u][0u], Rational(2));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(
            AlgebraTesting::context(),
            (SymmetricMatrix<Rational>(
                {{Rational(1),Rational(2)},{Rational(3),Rational(4)}})));
    }

    {
        InitializerList<InitializerList<Rational>> empty_matrix{};
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (Matrix<Rational>(empty_matrix)));
        InitializerList<InitializerList<Rational>> empty_row{{}};
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (Matrix<Rational>(empty_row)));
        InitializerList<InitializerList<Rational>> uneven_rows{{Rational(1),Rational(2)},{Rational(3)}};
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (Matrix<Rational>(uneven_rows)));

        Matrix<FloatDPBounds> nonsquare_bounds({{1.0_x,0.0_x}},pr);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), log_norm(nonsquare_bounds));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), lu_inverse(nonsquare_bounds));

        Matrix<FloatDPBounds> singular_zero({{0.0_x,0.0_x},{0.0_x,0.0_x}},pr);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), lu_inverse(singular_zero));

        Matrix<FloatDPBounds> square_bounds({{1.0_x,0.0_x},{0.0_x,1.0_x}},pr);
        Matrix<FloatDPBounds> bounds_rhs_ok({{1.0_x},{2.0_x}},pr);
        Matrix<FloatDPBounds> bounds_rhs_bad({{1.0_x}},pr);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), gs_solve(nonsquare_bounds,bounds_rhs_ok));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), gs_solve(square_bounds,bounds_rhs_bad));

        Matrix<FloatDPApproximation> nonsquare_approx({{1.0_x,0.0_x}},pr);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), triangular_decomposition(nonsquare_approx));

        Matrix<FloatDPApproximation> identity_approx=Matrix<FloatDPApproximation>::identity(2u,pr);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), triangular_decomposition(identity_approx));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), orthogonal_decomposition(identity_approx,true));
    }

    {
        Matrix<FloatDPBounds> A({{2.0_x,0.0_x},{0.0_x,3.0_x}},pr);
        Vector<FloatDPBounds> b({1.0_x,2.0_x},pr);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), sup_norm(A));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), log_norm(A));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), lu_solve(A,b));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), gs_solve(A,b));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), gs_inverse(A));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), midpoint(A));
    }

    {
        Matrix<FloatDPApproximation> A({{2.0_x,0.0_x},{0.0_x,3.0_x}},pr);
        Matrix<FloatDPApproximation> B=Matrix<FloatDPApproximation>::identity(2u,pr);
        Vector<FloatDPApproximation> b({1.0_x,2.0_x},pr);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), solve(A,B));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), solve(A,b));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), row_norms(A));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), orthogonal_decomposition(A,false));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), orthogonal_decomposition(A,true));
    }

    {
        Matrix<RoundedFloatDP> A({{1.0_x,1.0_x},{2.0_x,2.0_x}},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), normalise_rows(A));
    }

    {
        FloatDPUpperInterval z({0.0_x,0.0_x});
        Matrix<FloatDPUpperInterval> A(2u,2u,z);
        A[0][0]=FloatDPUpperInterval({1.0_x,1.0_x});
        A[1][1]=FloatDPUpperInterval({2.0_x,2.0_x});
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), cast_singleton(A));
    }

    Covector<FloatDPApproximation> empty_covector(0u,pr);
    std::ostringstream empty_covector_stream;
    empty_covector_stream << empty_covector;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), empty_covector_stream.str(), String("{}"));
#ifndef NDEBUG
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), empty_covector.zero_element());
#endif

    Covector<FloatDP> empty_float_covector(0u,dp);
    std::ostringstream empty_float_covector_stream;
    empty_float_covector_stream << empty_float_covector;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), empty_float_covector_stream.str(), String("{}"));

    Array<FloatDPApproximation> Aary(InitializerList<Dbl>{-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0},pr);
    Array<FloatDPBounds> iAary(InitializerList<ExactDouble>{-1.0_x,3.0_x, -1.0_x,1.0_x},pr);
    FloatDPApproximation* Aptr=Aary.begin();

    Matrix<FloatDPApproximation> A0(0,0,pr);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A0", A0);
    Matrix<FloatDPApproximation> A1(3,2,pr);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A1", A1);
    Matrix<FloatDPApproximation> A2(3,3,Aary.begin());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A2", A2);
    Matrix<FloatDPApproximation> A3({{-1.0,3.0,1.0}, {-1.0,1.0,2.0}, {2.0,1.0,1.0}},pr);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A3", A3);

    {
        Matrix<FloatDPApproximation> R(2u,3u,pr);
        R.resize(3u,2u);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), R.row_size(), 3u);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), R.column_size(), 2u);
        R.resize(4u,2u);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), R.row_size(), 4u);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), R.column_size(), 2u);
    }

    {
        Matrix<FloatDPApproximation> empty_row(1u,0u,pr);
        std::ostringstream oss;
        oss << empty_row;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !oss.str().empty());
    }

    {
        Matrix<FloatDPApproximation> const C({{1.0,2.0},{3.0,4.0}},pr);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), C[range(0,1)]);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), C.begin());
    }

    {
        auto rc=characteristics(Rational(0));
        Matrix<Rational> Z=Matrix<Rational>::zero(2u,3u,rc);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), Z.row_size(), 2u);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), Z.column_size(), 3u);
    }

    {
        Matrix<RoundedFloatDP> M({{1.0_x,2.0_x}},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), M*2);
    }

    {
        Matrix<FloatDPBounds> I=Matrix<FloatDPBounds>::identity(2u,dp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), I.row_size(), 2u);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), Matrix<FloatDPBounds>::identity<DoublePrecision>(2u,dp));
        MultiplePrecision mp(128);
        Matrix<FloatMPBounds> J=Matrix<FloatMPBounds>::identity(2u,mp);
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), J.column_size(), 2u);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), Matrix<FloatMPBounds>::identity<MultiplePrecision>(2u,mp));
    }

    {
        Matrix<FloatDPBounds> A({{1.0_x},{2.0_x}},dp);
        MatrixColumn<Matrix<FloatDPBounds>> c(A,0u);
        Vector<FloatDPBounds> v=c;
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), v.size(), 2u);
        Vector<FloatDPBounds> w({3.0_x,4.0_x},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), c=w);
    }

    {
        MultiplePrecision mp(128);
        Matrix<FloatMPApproximation> A({{1.0_x},{2.0_x}},mp);
        MatrixColumn<Matrix<FloatMPApproximation>> c(A,0u);
        Vector<FloatMPApproximation> v=c;
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), v.size(), 2u);
        Vector<FloatMPApproximation> w({3.0_x,4.0_x},mp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), c=w);
    }

    {
        MultiplePrecision mp(128);
        Matrix<FloatMPBounds> A({{1.0_x},{2.0_x}},mp);
        MatrixColumn<Matrix<FloatMPBounds>> c(A,0u);
        Vector<FloatMPBounds> v=c;
        ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), v.size(), 2u);
        Vector<FloatMPBounds> w({3.0_x,4.0_x},mp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), c=w);
    }

    {
        Matrix<RoundedFloatDP> A({{1.0_x,2.0_x}},dp);
        Matrix<RoundedFloatDP> row_mismatch({{1.0_x,2.0_x},{3.0_x,4.0_x}},dp);
        Matrix<RoundedFloatDP> col_mismatch({{1.0_x}},dp);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A+row_mismatch);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A+col_mismatch);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A-row_mismatch);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A-col_mismatch);
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A*Matrix<RoundedFloatDP>({{1.0_x}},dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A*Vector<RoundedFloatDP>({1.0_x},dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), A*transpose(Matrix<RoundedFloatDP>({{1.0_x}},dp)));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), transpose(A)*Matrix<RoundedFloatDP>({{1.0_x},{2.0_x}},dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), transpose(A)*Vector<RoundedFloatDP>({1.0_x,2.0_x},dp));

        Matrix<RoundedFloatDP> B({{3.0_x,4.0_x}},dp);
        Matrix<RoundedFloatDP> C({{1.0_x},{2.0_x}},dp);
        Vector<RoundedFloatDP> v({1.0_x,2.0_x},dp);
        Vector<RoundedFloatDP> one_v({1.0_x},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), A+B);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), A-B);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), A*C);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), A*v);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), A*transpose(B));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), transpose(A)*Matrix<RoundedFloatDP>({{1.0_x}},dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), transpose(A)*one_v);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(A==B));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(A==row_mismatch));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(A==col_mismatch));

        Matrix<FloatDP> exact_float({{1.0_x,2.0_x}},dp);
        Matrix<FloatDP> exact_float_row_mismatch({{1.0_x,2.0_x},{3.0_x,4.0_x}},dp);
        Matrix<FloatDP> exact_float_col_mismatch({{1.0_x}},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(exact_float==exact_float_row_mismatch));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(exact_float==exact_float_col_mismatch));

        Matrix<FloatDPApproximation> approx_matrix({{1.0_x,2.0_x}},dp);
        Matrix<FloatDPApproximation> approx_row_mismatch({{1.0_x,2.0_x},{3.0_x,4.0_x}},dp);
        Matrix<FloatDPApproximation> approx_col_mismatch({{1.0_x}},dp);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(approx_matrix==approx_row_mismatch));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(approx_matrix==approx_col_mismatch));

        Matrix<RoundedFloatDP> R(2u,2u,dp);
        R.resize(1u,4u);
        R.resize(3u,2u);
        R.set(0u,0u,RoundedFloatDP(1.0_x,dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), R.get(0u,0u));
        Matrix<RoundedFloatDP> const& CR=R;
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), CR.begin());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), CR[range(0,1)]);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), R[range(0,1)]);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), CR.element_characteristics());
        std::ostringstream rounded_stream;
        rounded_stream << R;

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), Matrix<RoundedFloatDP>::identity(2u,dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), Matrix<RoundedFloatDP>::identity<DoublePrecision>(2u,dp));

        Matrix<FloatDPBounds> bounds_resize(2u,2u,dp);
        bounds_resize.resize(1u,4u);
        bounds_resize.resize(3u,2u);

        Matrix<FloatDP> float_resize(2u,2u,dp);
        float_resize.resize(1u,4u);
        float_resize.resize(3u,2u);

        std::ostringstream empty_rounded_stream;
        empty_rounded_stream << Matrix<RoundedFloatDP>(0u,0u,dp);
        std::ostringstream empty_bounds_stream;
        empty_bounds_stream << Matrix<FloatDPBounds>(0u,0u,dp);
        std::ostringstream empty_float_stream;
        empty_float_stream << Matrix<FloatDP>(0u,0u,dp);

        Matrix<Rational> exact_a({{Rational(1),Rational(2)}});
        Matrix<Rational> exact_row_mismatch({{Rational(1),Rational(2)},{Rational(3),Rational(4)}});
        Matrix<Rational> exact_col_mismatch({{Rational(1)}});
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !(exact_a==exact_row_mismatch));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !(exact_a==exact_col_mismatch));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), static_cast<void>(exact_a==exact_a));
    }

    {
        InitializerList<InitializerList<ExactDouble>> bad_exact{{1.0_x,2.0_x},{3.0_x}};
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), Matrix<RoundedFloatDP>(bad_exact,dp));
        InitializerList<InitializerList<Dbl>> bad_dbl{{1.0,2.0},{3.0}};
        ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), Matrix<FloatDPApproximation>(bad_dbl,dp));
    }

    for(SizeType i=0; i!=A2.row_size(); ++i) {
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A2[i][j], A2.get(i,j));
            ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A2[i][j], Aptr[i*A2.column_size()+j]);
        }
    }

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A2, A3);

    A1[0][0]=1.0;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A1[0][0], 1.0);
    A1.set(0,1,3.0);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A1[0][1], 3.0);

    A1.at(1,0)=2.0;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), A1[1][0], 2.0);

    A0=Matrix<FloatDPApproximation>::zero(2,3,dp);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A0", A0);
    A1=Matrix<FloatDPApproximation>::identity(4,dp);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "A1", A1);

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), +FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr), FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), -FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr), FloatApproximationMatrix({{-1.0_x,-2.0_x},{-3.0_x,-4.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr)+FloatApproximationMatrix({{5.0_x,7.0_x},{8.0_x,6.0_x}},pr), FloatApproximationMatrix({{6.0_x,9.0_x},{11.0_x,10.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr)-FloatApproximationMatrix({{5.0_x,7.0_x},{8.0_x,6.0_x}},pr), FloatApproximationMatrix({{-4.0_x,-5.0_x},{-5.0_x,-2.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr)*FloatApproximationMatrix({{5.0_x,7.0_x},{8.0_x,6.0_x}},pr), FloatApproximationMatrix({{21.0_x,19.0_x},{47.0_x,45.0_x}},pr));

    // Transpose operations
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), FloatApproximationMatrix(transpose(FloatApproximationMatrix({{1.0_x,2.0_x,3.0_x},{4.0_x,5.0_x,6.0_x}},pr))), FloatApproximationMatrix({{1.0_x,4.0_x},{2.0_x,5.0_x},{3.0_x,6.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), transpose(FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr))*FloatApproximationMatrix({{5.0_x,7.0_x},{8.0_x,6.0_x}},pr), FloatApproximationMatrix({{1.0_x,3.0_x},{2.0_x,4.0_x}},pr)*FloatApproximationMatrix({{5.0_x,7.0_x},{8.0_x,6.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr)*transpose(FloatApproximationMatrix({{5.0_x,7.0_x},{8.0_x,6.0_x}},pr)), FloatApproximationMatrix({{1.0_x,2.0_x},{3.0_x,4.0_x}},pr)*FloatApproximationMatrix({{5.0_x,8.0_x},{7.0_x,6.0_x}},pr));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), transpose(FloatApproximationMatrix({{1.0_x,2.0_x,3.0_x},{4.0_x,5.0_x,6.0_x}},pr))*FloatDPApproximationVector({5.0_x,7.0_x},pr), FloatApproximationMatrix({{1.0_x,4.0_x},{2.0_x,5.0_x},{3.0_x,6.0_x}},pr)*FloatDPApproximationVector({5.0_x,7.0_x},pr));

}


int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("matrix");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

    TestMatrix().test();
    

    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}

