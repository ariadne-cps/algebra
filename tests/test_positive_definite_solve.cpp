/***************************************************************************
 *            test_positive_definite_solve.cpp
 *
 *  Copyright  2026  Luca Geretti
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

#include "numeric/numeric.hpp"
#include "algebra/symmetric_matrix.hpp"
#include "macro/testing.hpp"

#include <iostream>
#include <source_location>
#include <string>
#include <utility>

using namespace Ariadne;

namespace {

using X=FloatDPBounds;
using S=SymmetricMatrix<X>;
using Status=PositiveDefinitenessStatus;

// A compact inclusion audit: print one meaningful summary per check
// group and show the index and both actual intervals only on failure.
// Compare exact stored endpoints rather than logical LowerBound objects.
class InclusionAudit {
    std::string _description;
    SizeType _checked=0u;
    bool _ok=true;
  public:
    explicit InclusionAudit(std::string description) : _description(std::move(description)) { }

    Void add(const X& enclosure, const X& required, const std::string& component) {
        ++_checked;
        const bool lower_ok=enclosure.lower_raw()<=required.lower_raw();
        const bool upper_ok=enclosure.upper_raw()>=required.upper_raw();
        if(!lower_ok || !upper_ok) {
            _ok=false;
            std::cerr << "INCLUSION FAILED: " << _description
                      << " [" << component << "]\\n"
                      << "  required: " << required << "\n"
                      << "  enclosure: " << enclosure << "\n"
                      << "  lower bound OK: " << lower_ok
                      << ", upper bound OK: " << upper_ok << std::endl;
        }
    }

    bool valid() const { return _ok; }
};

// Known exact-point solutions are checked independently of residuals.
// The source line identifies which expected primary/dual component failed.
Void check_contains(const X& enclosure, const X& required,
                    std::source_location caller=std::source_location::current())
{
    InclusionAudit audit("known solution, test line "+std::to_string(caller.line()));
    audit.add(enclosure,required,"expected exact value");
    const bool known_exact_value_enclosed=audit.valid();
    ARIADNE_TEST_ASSERT(known_exact_value_enclosed);
}

// Check the algebraic residual as an interval inclusion, not equality of
// interval matrices. One aggregate report replaces per-endpoint spam.
Void check_matrix_solution(const S& A, const Matrix<X>& B,
                           const Matrix<X>& enclosure)
{
    ARIADNE_TEST_EQUALS(enclosure.row_size(),B.row_size());
    ARIADNE_TEST_EQUALS(enclosure.column_size(),B.column_size());
    const Matrix<X> dense_A=static_cast<Matrix<X>>(A);
    const Matrix<X> residual=dense_A*enclosure;
    InclusionAudit audit("SPD matrix residual: B contained in A*X");
    for(SizeType i=0; i!=B.row_size(); ++i) {
        for(SizeType j=0; j!=B.column_size(); ++j) {
            audit.add(residual[i][j],B[i][j],
                      "row "+std::to_string(i)+", column "+std::to_string(j));
        }
    }
    const bool matrix_residual_enclosed=audit.valid();
    ARIADNE_TEST_ASSERT(matrix_residual_enclosed);
}

Void check_vector_solution(const S& A, const Vector<X>& b,
                           const Vector<X>& enclosure)
{
    ARIADNE_TEST_EQUALS(enclosure.size(),b.size());
    const Matrix<X> dense_A=static_cast<Matrix<X>>(A);
    const Vector<X> residual=dense_A*enclosure;
    InclusionAudit audit("SPD vector residual: b contained in A*x");
    for(SizeType i=0; i!=b.size(); ++i) {
        audit.add(residual[i],b[i],"row "+std::to_string(i));
    }
    const bool vector_residual_enclosed=audit.valid();
    ARIADNE_TEST_ASSERT(vector_residual_enclosed);
}

Void check_saddle_point_residuals(
    const S& H, const Matrix<X>& A, const Vector<X>& f, const Vector<X>& g,
    const Vector<X>& x, const Vector<X>& second_block)
{
    ARIADNE_TEST_EQUALS(x.size(),H.size());
    ARIADNE_TEST_EQUALS(second_block.size(),A.row_size());
    const Matrix<X> dense_H=static_cast<Matrix<X>>(H);
    const Vector<X> stationary=dense_H*x+transpose(A)*second_block;
    const Vector<X> coupled=A*x;
    InclusionAudit primary_audit("saddle-point stationarity: f contained in H*x+A^T*second_block");
    for(SizeType i=0; i!=f.size(); ++i) {
        primary_audit.add(stationary[i],f[i],"row "+std::to_string(i));
    }
    const bool saddle_point_stationarity_enclosed=primary_audit.valid();
    ARIADNE_TEST_ASSERT(saddle_point_stationarity_enclosed);

    InclusionAudit equality_audit("saddle-point equalities: g contained in A*x");
    for(SizeType i=0; i!=g.size(); ++i) {
        equality_audit.add(coupled[i],g[i],"row "+std::to_string(i));
    }
    const bool saddle_point_equalities_enclosed=equality_audit.valid();
    ARIADNE_TEST_ASSERT(saddle_point_equalities_enclosed);
}

// Check both primary and secondary variable components against a selected exact
// member of an interval saddle-point family, with a meaningful scenario label.
Void check_known_saddle_point_solution(
    const ValidatedSaddlePointResult& solution,
    const Vector<X>& expected_primary,
    const Vector<X>& expected_secondary,
    const std::string& scenario)
{
    ARIADNE_TEST_ASSERT(solution.primary.has_value());
    ARIADNE_TEST_ASSERT(solution.secondary.has_value());
    if(!solution.primary.has_value() || !solution.secondary.has_value()) {
        return;
    }
    ARIADNE_TEST_EQUALS(solution.primary->size(),expected_primary.size());
    ARIADNE_TEST_EQUALS(solution.secondary->size(),expected_secondary.size());
    InclusionAudit audit("saddle-point known solution: "+scenario);
    for(SizeType i=0; i!=expected_primary.size(); ++i) {
        audit.add((*solution.primary)[i],expected_primary[i],
                  "primary["+std::to_string(i)+"]");
    }
    for(SizeType i=0; i!=expected_secondary.size(); ++i) {
        audit.add((*solution.secondary)[i],expected_secondary[i],
                  "secondary variable["+std::to_string(i)+"]");
    }
    const bool saddle_point_reference_solution_enclosed=audit.valid();
    ARIADNE_TEST_ASSERT(saddle_point_reference_solution_enclosed);
}

Void test_validated_saddle_point_solve()
{
    const X zero(0u,dp);
    const Vector<X> empty_vector(0u,zero);

    // Both the zero-dimensional saddle-point system and the uncoupled
    // case have unique solutions under positive definiteness.
    const S empty_leading_block(0u,zero);
    const Matrix<X> empty_coupling_rows(0u,0u,zero);
    const auto empty_solution=validated_saddle_point_solve(
        empty_leading_block,empty_coupling_rows,empty_vector,empty_vector);
    ARIADNE_TEST_ASSERT(empty_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(empty_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(empty_solution.secondary.has_value());
    check_saddle_point_residuals(empty_leading_block,empty_coupling_rows,empty_vector,
                        empty_vector,empty_solution.primary.value(),empty_solution.secondary.value());

    S scalar_leading_block(1u,zero);
    scalar_leading_block[0u][0u]=X(2u,dp);
    const Matrix<X> no_coupling_rows(0u,1u,zero);
    const Vector<X> uncoupled_rhs({X(6u,dp)});
    const auto uncoupled_solution=validated_saddle_point_solve(
        scalar_leading_block,no_coupling_rows,uncoupled_rhs,empty_vector);
    ARIADNE_TEST_ASSERT(uncoupled_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(uncoupled_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(uncoupled_solution.secondary.has_value());
    check_known_saddle_point_solution(uncoupled_solution,Vector<X>({X(3u,dp)}),
                             empty_vector,"uncoupled SPD system");
    check_saddle_point_residuals(scalar_leading_block,no_coupling_rows,uncoupled_rhs,
                        empty_vector,uncoupled_solution.primary.value(),
                        uncoupled_solution.secondary.value());

    // One equality: H=2I, A=(1,1), f=0, g=2.
    // Exact solution x=(1,1), second_block=-2.
    S leading_block(2u,zero);
    leading_block[0u][0u]=X(2u,dp);
    leading_block[1u][1u]=X(2u,dp);
    Matrix<X> coupling_rows(1u,2u,zero);
    coupling_rows[0u][0u]=X(1u,dp);
    coupling_rows[0u][1u]=X(1u,dp);
    const Vector<X> f(2u,zero);
    const Vector<X> g({X(2u,dp)});
    const auto single_solution=validated_saddle_point_solve(leading_block,coupling_rows,f,g);
    ARIADNE_TEST_ASSERT(single_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(single_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(single_solution.secondary.has_value());
    check_known_saddle_point_solution(single_solution,
                             Vector<X>({X(1u,dp),X(1u,dp)}),
                             Vector<X>({X(-2,dp)}),"one equality");
    check_saddle_point_residuals(leading_block,coupling_rows,f,g,
                        single_solution.primary.value(),single_solution.secondary.value());

    // A non-diagonal leading block and a nonzero dual variable.
    S coupled_leading_block(2u,zero);
    coupled_leading_block[0u][0u]=X(2u,dp);
    coupled_leading_block[0u][1u]=X(1u,dp);
    coupled_leading_block[1u][1u]=X(2u,dp);
    Matrix<X> first_coordinate(1u,2u,zero);
    first_coordinate[0u][0u]=X(1u,dp);
    const Vector<X> coupled_f({X(11u,dp),X(8u,dp)});
    const Vector<X> coupled_g({X(2u,dp)});
    const auto coupled_solution=validated_saddle_point_solve(
        coupled_leading_block,first_coordinate,coupled_f,coupled_g);
    ARIADNE_TEST_ASSERT(coupled_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(coupled_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(coupled_solution.secondary.has_value());
    check_known_saddle_point_solution(coupled_solution,
                             Vector<X>({X(2u,dp),X(3u,dp)}),
                             Vector<X>({X(4u,dp)}),"coupled leading block");
    check_saddle_point_residuals(coupled_leading_block,first_coordinate,coupled_f,coupled_g,
                        coupled_solution.primary.value(),coupled_solution.secondary.value());

    // Two independent coupling_rows: exact x=(1,2,3), second_block=(4,5).
    S h3(3u,zero);
    for(SizeType i=0; i!=3u; ++i) {
        h3[i][i]=X(1u,dp);
    }
    Matrix<X> a2(2u,3u,zero);
    a2[0u][0u]=X(1u,dp);
    a2[0u][2u]=X(1u,dp);
    a2[1u][1u]=X(1u,dp);
    a2[1u][2u]=X(1u,dp);
    const Vector<X> f3({X(5u,dp),X(7u,dp),X(12u,dp)});
    const Vector<X> g2({X(4u,dp),X(5u,dp)});
    const auto multi_solution=validated_saddle_point_solve(h3,a2,f3,g2);
    ARIADNE_TEST_ASSERT(multi_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(multi_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(multi_solution.secondary.has_value());
    check_saddle_point_residuals(h3,a2,f3,g2,
                        multi_solution.primary.value(),multi_solution.secondary.value());
    check_known_saddle_point_solution(multi_solution,
                             Vector<X>({X(1u,dp),X(2u,dp),X(3u,dp)}),
                             Vector<X>({X(4u,dp),X(5u,dp)}),"two equalities");

    // A family of admissible saddle-point systems with a known exact member.
    S interval_leading_block(2u,zero);
    interval_leading_block[0u][0u]=X(4,5,dp);
    interval_leading_block[1u][1u]=X(4,5,dp);
    Matrix<X> interval_coupling_rows(1u,2u,zero);
    interval_coupling_rows[0u][0u]=X(1,2,dp);
    const Vector<X> interval_f({X(4,5,dp),X(8,10,dp)});
    const Vector<X> interval_g({X(1,2,dp)});
    const auto interval_solution=validated_saddle_point_solve(
        interval_leading_block,interval_coupling_rows,interval_f,interval_g);
    ARIADNE_TEST_ASSERT(interval_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(interval_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(interval_solution.secondary.has_value());
    check_saddle_point_residuals(interval_leading_block,interval_coupling_rows,interval_f,interval_g,
                        interval_solution.primary.value(),interval_solution.secondary.value());
    check_known_saddle_point_solution(interval_solution,
                             Vector<X>({X(1u,dp),X(2u,dp)}),
                             Vector<X>({zero}),"interval saddle-point family");

    // A provably non-SPD leading block is a failure of this solver's precondition.
    S negative(1u,zero);
    negative[0u][0u]=X(-1,dp);
    Matrix<X> scalar_coupling_rows(1u,1u,X(1u,dp));
    const Vector<X> scalar_rhs({X(1u,dp)});
    const auto negative_solution=validated_saddle_point_solve(
        negative,scalar_coupling_rows,scalar_rhs,scalar_rhs);
    ARIADNE_TEST_ASSERT(
        negative_solution.status==ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(!negative_solution.secondary.has_value());

    S uncertain_leading_block(1u,zero);
    uncertain_leading_block[0u][0u]=X(-1,1,dp);
    const auto uncertain_h_solution=validated_saddle_point_solve(
        uncertain_leading_block,scalar_coupling_rows,scalar_rhs,scalar_rhs);
    ARIADNE_TEST_ASSERT(uncertain_h_solution.status==ValidatedSaddlePointStatus::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_h_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(!uncertain_h_solution.secondary.has_value());

    // More coupling_row rows than unknowns is universally rank deficient.
    Matrix<X> overcoupled(2u,1u,X(1u,dp));
    const Vector<X> two_coupling_rows({X(1u,dp),X(1u,dp)});
    const auto dependent_solution=validated_saddle_point_solve(
        scalar_leading_block,overcoupled,scalar_rhs,two_coupling_rows);
    ARIADNE_TEST_ASSERT(
        dependent_solution.status==ValidatedSaddlePointStatus::DEPENDENT_COUPLING_ROWS);
    ARIADNE_TEST_ASSERT(!dependent_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(!dependent_solution.secondary.has_value());

    // Repeated coupling_rows are dependent even without dimension excess.
    Matrix<X> repeated(2u,2u,zero);
    repeated[0u][0u]=X(1u,dp);
    repeated[1u][0u]=X(1u,dp);
    const auto repeated_solution=validated_saddle_point_solve(
        leading_block,repeated,f,two_coupling_rows);
    ARIADNE_TEST_ASSERT(
        repeated_solution.status==ValidatedSaddlePointStatus::DEPENDENT_COUPLING_ROWS);
    ARIADNE_TEST_ASSERT(!repeated_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(!repeated_solution.secondary.has_value());

    // Interval dependence of the same A entry is unresolved. Never use
    // a midpoint-rank decision, and do not infer dependence from an
    // interval Schur matrix which overestimates correlations.
    Matrix<X> ambiguous(1u,1u,X(-1,1,dp));
    const auto ambiguous_solution=validated_saddle_point_solve(
        scalar_leading_block,ambiguous,scalar_rhs,scalar_rhs);
    ARIADNE_TEST_ASSERT(ambiguous_solution.status==ValidatedSaddlePointStatus::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!ambiguous_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(!ambiguous_solution.secondary.has_value());

    Matrix<X> bad_columns(1u,2u,zero);
    const Vector<X> bad_f(2u,zero);
    const Vector<X> bad_g(2u,zero);
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(scalar_leading_block,bad_columns,scalar_rhs,scalar_rhs));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(scalar_leading_block,scalar_coupling_rows,bad_f,scalar_rhs));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(scalar_leading_block,scalar_coupling_rows,scalar_rhs,bad_g));
}


Void test_validated_saddle_point_matrix_rhs()
{
    const X zero(0u,dp);
    S H(2u,zero);
    H[0u][0u]=X(2u,dp);
    H[1u][1u]=X(2u,dp);
    Matrix<X> A(1u,2u,zero);
    A[0u][0u]=X(1u,dp);
    A[0u][1u]=X(1u,dp);

    // First column: x=(1,1), second_block=-2; second: x=(2,3), second_block=4.
    Matrix<X> F(2u,2u,zero);
    F[0u][1u]=X(8u,dp);
    F[1u][1u]=X(10u,dp);
    Matrix<X> G(1u,2u,zero);
    G[0u][0u]=X(2u,dp);
    G[0u][1u]=X(5u,dp);

    const auto multiple=validated_saddle_point_solve(H,A,F,G);
    ARIADNE_TEST_ASSERT(multiple.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(multiple.primary.has_value());
    ARIADNE_TEST_ASSERT(multiple.secondary.has_value());
    if(multiple.primary.has_value() && multiple.secondary.has_value()) {
        const Matrix<X>& X_solution=*multiple.primary;
        const Matrix<X>& Secondary=*multiple.secondary;
        ARIADNE_TEST_EQUALS(X_solution.row_size(),2u);
        ARIADNE_TEST_EQUALS(X_solution.column_size(),2u);
        ARIADNE_TEST_EQUALS(Secondary.row_size(),1u);
        ARIADNE_TEST_EQUALS(Secondary.column_size(),2u);
        const Matrix<X> residual=static_cast<Matrix<X>>(H)*X_solution+
                                 transpose(A)*Secondary;
        const Matrix<X> equality=A*X_solution;
        InclusionAudit audit("saddle-point two-right-hand-side residual");
        for(SizeType i=0; i!=2u; ++i) {
            for(SizeType col=0; col!=2u; ++col) {
                audit.add(residual[i][col],F[i][col],
                          "stationarity["+std::to_string(i)+","+std::to_string(col)+"]");
            }
        }
        for(SizeType col=0; col!=2u; ++col) {
            audit.add(equality[0u][col],G[0u][col],
                      "equality[0,"+std::to_string(col)+"]");
        }
        audit.add(X_solution[0u][0u],X(1u,dp),"primary[0,0]");
        audit.add(X_solution[1u][0u],X(1u,dp),"primary[1,0]");
        audit.add(Secondary[0u][0u],X(-2,dp),"second_block[0,0]");
        audit.add(X_solution[0u][1u],X(2u,dp),"primary[0,1]");
        audit.add(X_solution[1u][1u],X(3u,dp),"primary[1,1]");
        audit.add(Secondary[0u][1u],X(4u,dp),"second_block[0,1]");
        const bool both_saddle_point_columns_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(both_saddle_point_columns_enclosed);
    }

    // A zero-column right-hand side must still certify the shared
    // saddle-point operator and produce correctly shaped, empty solutions.
    Matrix<X> empty_F(2u,0u,zero);
    Matrix<X> empty_G(1u,0u,zero);
    const auto no_rhs=validated_saddle_point_solve(H,A,empty_F,empty_G);
    ARIADNE_TEST_ASSERT(no_rhs.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(no_rhs.primary.has_value());
    ARIADNE_TEST_ASSERT(no_rhs.secondary.has_value());
    if(no_rhs.primary.has_value() && no_rhs.secondary.has_value()) {
        ARIADNE_TEST_EQUALS(no_rhs.primary->row_size(),2u);
        ARIADNE_TEST_EQUALS(no_rhs.primary->column_size(),0u);
        ARIADNE_TEST_EQUALS(no_rhs.secondary->row_size(),1u);
        ARIADNE_TEST_EQUALS(no_rhs.secondary->column_size(),0u);
    }

    const Matrix<X> wrong_F(1u,2u,zero);
    const Matrix<X> wrong_G_rows(2u,2u,zero);
    const Matrix<X> wrong_G_cols(1u,1u,zero);
    const Matrix<X> wrong_A(1u,3u,zero);
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,wrong_A,F,G));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,wrong_F,G));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,F,wrong_G_rows));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,F,wrong_G_cols));

    // Verify failure classification for the matrix overload as well.
    S nonspd(2u,zero);
    nonspd[0u][0u]=X(-1,dp);
    nonspd[1u][1u]=X(1u,dp);
    const auto rejected_h=validated_saddle_point_solve(nonspd,A,F,G);
    ARIADNE_TEST_ASSERT(rejected_h.status==ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!rejected_h.primary.has_value());
    ARIADNE_TEST_ASSERT(!rejected_h.secondary.has_value());

    S uncertain(2u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    uncertain[1u][1u]=X(1u,dp);
    const auto rejected_uncertain=validated_saddle_point_solve(uncertain,A,F,G);
    ARIADNE_TEST_ASSERT(rejected_uncertain.status==ValidatedSaddlePointStatus::INDETERMINATE);

    Matrix<X> dependent_A(2u,2u,zero);
    dependent_A[0u][0u]=X(1u,dp);
    dependent_A[1u][0u]=X(1u,dp);
    const Matrix<X> G2(2u,2u,zero);
    const auto rejected_rank=validated_saddle_point_solve(H,dependent_A,F,G2);
    ARIADNE_TEST_ASSERT(rejected_rank.status==ValidatedSaddlePointStatus::DEPENDENT_COUPLING_ROWS);

    Matrix<X> ambiguous_A(1u,2u,zero);
    ambiguous_A[0u][0u]=X(-1,1,dp);
    const auto rejected_ambiguous=validated_saddle_point_solve(H,ambiguous_A,F,G);
    ARIADNE_TEST_ASSERT(rejected_ambiguous.status==ValidatedSaddlePointStatus::INDETERMINATE);

    // Zero-dimensional operator with multiple columns.
    S h0(0u,zero);
    const Matrix<X> a0(0u,0u,zero);
    const Matrix<X> f0(0u,2u,zero);
    const Matrix<X> g0(0u,2u,zero);
    const auto zero_system=validated_saddle_point_solve(h0,a0,f0,g0);
    ARIADNE_TEST_ASSERT(zero_system.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(zero_system.primary.has_value());
    ARIADNE_TEST_ASSERT(zero_system.secondary.has_value());
    if(zero_system.primary.has_value() && zero_system.secondary.has_value()) {
        ARIADNE_TEST_EQUALS(zero_system.primary->column_size(),2u);
        ARIADNE_TEST_EQUALS(zero_system.secondary->column_size(),2u);
    }

    // Uncoupled multi-RHS: the same SPD solve is used with m=0.
    const Matrix<X> no_A(0u,2u,zero);
    const Matrix<X> no_G(0u,2u,zero);
    const auto uncoupled=validated_saddle_point_solve(H,no_A,F,no_G);
    ARIADNE_TEST_ASSERT(uncoupled.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(uncoupled.primary.has_value());
    ARIADNE_TEST_ASSERT(uncoupled.secondary.has_value());
    if(uncoupled.primary.has_value() && uncoupled.secondary.has_value()) {
        ARIADNE_TEST_EQUALS(uncoupled.primary->column_size(),2u);
        ARIADNE_TEST_EQUALS(uncoupled.secondary->row_size(),0u);
    }

    // Interval family with a certified unique solution for each member.
    S interval_H(2u,zero);
    interval_H[0u][0u]=X(4,5,dp);
    interval_H[1u][1u]=X(4,5,dp);
    Matrix<X> interval_A(1u,2u,zero);
    interval_A[0u][0u]=X(1,2,dp);
    Matrix<X> interval_F(2u,2u,zero);
    interval_F[0u][0u]=X(4,5,dp);
    interval_F[1u][0u]=X(8,10,dp);
    interval_F[0u][1u]=X(8,10,dp);
    interval_F[1u][1u]=X(4,5,dp);
    Matrix<X> interval_G(1u,2u,zero);
    interval_G[0u][0u]=X(1,2,dp);
    interval_G[0u][1u]=X(1,2,dp);
    const auto interval_system=validated_saddle_point_solve(
        interval_H,interval_A,interval_F,interval_G);
    ARIADNE_TEST_ASSERT(interval_system.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(interval_system.primary.has_value());
    ARIADNE_TEST_ASSERT(interval_system.secondary.has_value());
    if(interval_system.primary.has_value() && interval_system.secondary.has_value()) {
        InclusionAudit audit("saddle-point interval family, two right-hand sides");
        const Matrix<X> exact_X=static_cast<Matrix<X>>(interval_H)*(*interval_system.primary)+
                                 transpose(interval_A)*(*interval_system.secondary);
        const Matrix<X> exact_G=interval_A*(*interval_system.primary);
        for(SizeType i=0; i!=2u; ++i) {
            for(SizeType col=0; col!=2u; ++col) {
                audit.add(exact_X[i][col],interval_F[i][col],
                          "stationarity["+std::to_string(i)+","+std::to_string(col)+"]");
            }
        }
        for(SizeType col=0; col!=2u; ++col) {
            audit.add(exact_G[0u][col],interval_G[0u][col],
                      "equality[0,"+std::to_string(col)+"]");
        }
        const bool family_rhs_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(family_rhs_enclosed);
    }
}


// A semidefinite leading block can still define a strictly convex quadratic
// on the affine feasible space when H is positive on ker(A). Test the
// augmented route without weakening any of the original SPD tests.
Void test_saddle_point_augmented_leading_block()
{
    const X zero(0u,dp);

    // H=diag(0,2) is positive semidefinite but not SPD.
    // A fixes the first coordinate. The exact solution is
    // x=(2,3), second_block=5, for f=(5,6) and g=2.
    S semidefinite(2u,zero);
    semidefinite[1u][1u]=X(2u,dp);
    Matrix<X> A(1u,2u,zero);
    A[0u][0u]=X(1u,dp);
    const Vector<X> f({X(5u,dp),X(6u,dp)});
    const Vector<X> g({X(2u,dp)});
    const auto coupled=validated_saddle_point_solve(semidefinite,A,f,g);
    ARIADNE_TEST_ASSERT(coupled.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(coupled.primary.has_value());
    ARIADNE_TEST_ASSERT(coupled.secondary.has_value());
    if(coupled.primary.has_value() && coupled.secondary.has_value()) {
        check_known_saddle_point_solution(coupled,
                                 Vector<X>({X(2u,dp),X(3u,dp)}),
                                 Vector<X>({X(5u,dp)}),
                                 "positive semidefinite H, coupled nullspace");
        check_saddle_point_residuals(semidefinite,A,f,g,
                            *coupled.primary,*coupled.secondary);
    }

    // The original H=[-1,1] is not certifiably SPD, but A=[2]
    // makes H+A^T*A=[3,5] certified SPD for every interval member.
    // The member H=0 has x=1 and second_block=2.
    S uncertain(1u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    Matrix<X> A_scalar(1u,1u,X(2u,dp));
    const Vector<X> f_scalar({X(4u,dp)});
    const Vector<X> g_scalar({X(2u,dp)});
    const auto uncertain_coupled=validated_saddle_point_solve(
        uncertain,A_scalar,f_scalar,g_scalar);
    ARIADNE_TEST_ASSERT(uncertain_coupled.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(uncertain_coupled.primary.has_value());
    ARIADNE_TEST_ASSERT(uncertain_coupled.secondary.has_value());
    if(uncertain_coupled.primary.has_value()
       && uncertain_coupled.secondary.has_value()) {
        check_known_saddle_point_solution(uncertain_coupled,
                                 Vector<X>({X(1u,dp)}),
                                 Vector<X>({X(2u,dp)}),
                                 "interval H with coupled negative curvature");
        check_saddle_point_residuals(uncertain,A_scalar,f_scalar,g_scalar,
                            *uncertain_coupled.primary,
                            *uncertain_coupled.secondary);
    }

    // Even H=-1 is admissible when all negative-curvature directions
    // are fixed by A=[2]. The secondary variable in the original saddle-point system
    // is unchanged by the augmentation.
    S negative_scalar(1u,zero);
    negative_scalar[0u][0u]=X(-1,dp);
    const Vector<X> negative_f({X(3u,dp)});
    const auto negative_coupled=validated_saddle_point_solve(
        negative_scalar,A_scalar,negative_f,g_scalar);
    ARIADNE_TEST_ASSERT(negative_coupled.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(negative_coupled.primary.has_value());
    ARIADNE_TEST_ASSERT(negative_coupled.secondary.has_value());
    if(negative_coupled.primary.has_value()
       && negative_coupled.secondary.has_value()) {
        check_known_saddle_point_solution(negative_coupled,
                                 Vector<X>({X(1u,dp)}),
                                 Vector<X>({X(2u,dp)}),
                                 "negative leading block, coupled degree of freedom");
        check_saddle_point_residuals(negative_scalar,A_scalar,negative_f,g_scalar,
                            *negative_coupled.primary,
                            *negative_coupled.secondary);
    }

    // H=0 and full-rank A=I: the coupling_rows determine every primary
    // variable and second_block=F. Check two right-hand sides at once.
    S zero_leading_block(2u,zero);
    Matrix<X> identity(2u,2u,zero);
    identity[0u][0u]=X(1u,dp);
    identity[1u][1u]=X(1u,dp);
    Matrix<X> F(2u,2u,zero);
    F[0u][0u]=X(3u,dp);
    F[1u][0u]=X(5u,dp);
    F[0u][1u]=X(7u,dp);
    F[1u][1u]=X(11u,dp);
    Matrix<X> G(2u,2u,zero);
    G[0u][0u]=X(1u,dp);
    G[1u][0u]=X(2u,dp);
    G[0u][1u]=X(4u,dp);
    G[1u][1u]=X(6u,dp);
    const auto zero_full_rank=validated_saddle_point_solve(
        zero_leading_block,identity,F,G);
    ARIADNE_TEST_ASSERT(zero_full_rank.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(zero_full_rank.primary.has_value());
    ARIADNE_TEST_ASSERT(zero_full_rank.secondary.has_value());
    if(zero_full_rank.primary.has_value()
       && zero_full_rank.secondary.has_value()) {
        InclusionAudit audit("zero leading block, full-rank coupling_rows, two RHS");
        const Matrix<X>& x=*zero_full_rank.primary;
        const Matrix<X>& second_block=*zero_full_rank.secondary;
        const Matrix<X> stationarity=static_cast<Matrix<X>>(zero_leading_block)*x+
                                     transpose(identity)*second_block;
        const Matrix<X> equalities=identity*x;
        for(SizeType i=0; i!=2u; ++i) {
            for(SizeType col=0; col!=2u; ++col) {
                audit.add(x[i][col],G[i][col],
                          "primary["+std::to_string(i)+","+std::to_string(col)+"]");
                audit.add(second_block[i][col],F[i][col],
                          "secondary variable["+std::to_string(i)+","+std::to_string(col)+"]");
                audit.add(stationarity[i][col],F[i][col],
                          "stationarity["+std::to_string(i)+","+std::to_string(col)+"]");
                audit.add(equalities[i][col],G[i][col],
                          "equality["+std::to_string(i)+","+std::to_string(col)+"]");
            }
        }
        const bool coupled_zero_leading_block_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(coupled_zero_leading_block_enclosed);
    }

    // A zero-curvature direction not coupled by A cannot be
    // certified with H+A^T*A; this is not an infeasibility result.
    Matrix<X> missing_direction(1u,2u,zero);
    missing_direction[0u][1u]=X(1u,dp);
    const auto incomplete=validated_saddle_point_solve(
        semidefinite,missing_direction,f,g);
    ARIADNE_TEST_ASSERT(
        incomplete.status==ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!incomplete.primary.has_value());
    ARIADNE_TEST_ASSERT(!incomplete.secondary.has_value());

    // No coupling_rows: a semidefinite leading block remains uncertifiable by
    // either direct or augmented SPD tests.
    const Matrix<X> no_coupling_rows(0u,2u,zero);
    const Vector<X> empty_g(0u,zero);
    const auto uncoupled=validated_saddle_point_solve(
        semidefinite,no_coupling_rows,f,empty_g);
    ARIADNE_TEST_ASSERT(
        uncoupled.status==ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!uncoupled.primary.has_value());
    ARIADNE_TEST_ASSERT(!uncoupled.secondary.has_value());
}


Void test_generalized_saddle_point()
{
    const X zero(0u,dp);
    S H(2u,zero);
    H[0u][0u]=X(2u,dp);
    H[1u][1u]=X(2u,dp);
    Matrix<X> A(1u,2u,zero);
    A[0u][0u]=X(1u,dp);
    A[0u][1u]=X(1u,dp);
    S C(1u,zero);
    C[0u][0u]=X(2u,dp);
    const Vector<X> f(2u,zero);
    const Vector<X> g({X(6u,dp)});
    const auto one=validated_saddle_point_solve(H,A,C,f,g);
    ARIADNE_TEST_ASSERT(one.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(one.primary.has_value());
    ARIADNE_TEST_ASSERT(one.secondary.has_value());
    if(one.primary.has_value() && one.secondary.has_value()) {
        check_known_saddle_point_solution(one,Vector<X>({X(1u,dp),X(1u,dp)}),
                                          Vector<X>({X(-2,dp)}),"nonzero lower block");
        InclusionAudit audit("generalized saddle-point vector residual");
        const Vector<X> upper=static_cast<Matrix<X>>(H)*(*one.primary)+
                              transpose(A)*(*one.secondary);
        const Vector<X> lower=A*(*one.primary)-static_cast<Matrix<X>>(C)*(*one.secondary);
        for(SizeType i=0; i!=2u; ++i) {
            audit.add(upper[i],f[i],"upper["+std::to_string(i)+"]");
        }
        audit.add(lower[0u],g[0u],"lower[0]");
        const bool generalized_vector_residual_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(generalized_vector_residual_enclosed);
    }

    Matrix<X> F(2u,2u,zero);
    F[0u][1u]=X(8u,dp);
    F[1u][1u]=X(10u,dp);
    Matrix<X> G(1u,2u,zero);
    G[0u][0u]=X(6u,dp);
    G[0u][1u]=X(-3,dp);
    const auto many=validated_saddle_point_solve(H,A,C,F,G);
    ARIADNE_TEST_ASSERT(many.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(many.primary.has_value());
    ARIADNE_TEST_ASSERT(many.secondary.has_value());
    if(many.primary.has_value() && many.secondary.has_value()) {
        const Matrix<X> upper=static_cast<Matrix<X>>(H)*(*many.primary)+
                              transpose(A)*(*many.secondary);
        const Matrix<X> lower=A*(*many.primary)-
                              static_cast<Matrix<X>>(C)*(*many.secondary);
        InclusionAudit audit("generalized saddle-point two RHS");
        for(SizeType i=0; i!=2u; ++i) {
            for(SizeType j=0; j!=2u; ++j) {
                audit.add(upper[i][j],F[i][j],"upper");
            }
        }
        for(SizeType j=0; j!=2u; ++j) {
            audit.add(lower[0u][j],G[0u][j],"lower");
        }
        audit.add((*many.primary)[0u][0u],X(1u,dp),"x00");
        audit.add((*many.primary)[1u][0u],X(1u,dp),"x10");
        audit.add((*many.secondary)[0u][0u],X(-2,dp),"y00");
        audit.add((*many.primary)[0u][1u],X(2u,dp),"x01");
        audit.add((*many.primary)[1u][1u],X(3u,dp),"x11");
        audit.add((*many.secondary)[0u][1u],X(4u,dp),"y01");
        const bool generalized_matrix_residual_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(generalized_matrix_residual_enclosed);
    }

    // Dependent coupling rows are allowed when C regularizes the lower block.
    S scalar(1u,zero);
    scalar[0u][0u]=X(2u,dp);
    Matrix<X> repeated(2u,1u,X(1u,dp));
    S identity(2u,zero);
    identity[0u][0u]=X(1u,dp);
    identity[1u][1u]=X(1u,dp);
    const auto dependent=validated_saddle_point_solve(
        scalar,repeated,identity,Vector<X>({X(15u,dp)}),
        Vector<X>({X(-1,dp),X(-2,dp)}));
    ARIADNE_TEST_ASSERT(dependent.status==ValidatedSaddlePointStatus::SOLVED);
    if(dependent.primary.has_value() && dependent.secondary.has_value()) {
        check_known_saddle_point_solution(dependent,
            Vector<X>({X(3u,dp)}),Vector<X>({X(4u,dp),X(5u,dp)}),
            "dependent coupling, nonzero lower block");
    }

    Matrix<X> no_coupling(2u,1u,zero);
    const auto decoupled=validated_saddle_point_solve(
        scalar,no_coupling,identity,Vector<X>({X(6u,dp)}),
        Vector<X>({X(-4,dp),X(-5,dp)}));
    ARIADNE_TEST_ASSERT(decoupled.status==ValidatedSaddlePointStatus::SOLVED);
    if(decoupled.primary.has_value() && decoupled.secondary.has_value()) {
        check_known_saddle_point_solution(decoupled,
            Vector<X>({X(3u,dp)}),Vector<X>({X(4u,dp),X(5u,dp)}),
            "zero coupling, nonzero lower block");
    }

    // S can be SPD even when C is indefinite.
    Matrix<X> strong(1u,1u,X(2u,dp));
    S negative_C(1u,zero);
    negative_C[0u][0u]=X(-1,dp);
    const auto indefinite=validated_saddle_point_solve(
        scalar,strong,negative_C,Vector<X>({X(6u,dp)}),
        Vector<X>({X(4u,dp)}));
    ARIADNE_TEST_ASSERT(indefinite.status==ValidatedSaddlePointStatus::SOLVED);
    if(indefinite.primary.has_value() && indefinite.secondary.has_value()) {
        check_known_saddle_point_solution(indefinite,
            Vector<X>({X(1u,dp)}),Vector<X>({X(2u,dp)}),
            "indefinite lower-right block");
    }

    // A positive interval C gives certified enclosures for a family.
    S interval_C(1u,zero);
    interval_C[0u][0u]=X(1,2,dp);
    const auto interval=validated_saddle_point_solve(
        scalar,Matrix<X>(1u,1u,X(1u,dp)),interval_C,
        Vector<X>({X(4u,dp)}),Vector<X>({X(-3,-1,dp)}));
    ARIADNE_TEST_ASSERT(interval.status==ValidatedSaddlePointStatus::SOLVED);
    if(interval.primary.has_value() && interval.secondary.has_value()) {
        check_known_saddle_point_solution(interval,
            Vector<X>({X(1u,dp)}),Vector<X>({X(2u,dp)}),
            "interval lower-right block");
    }

    // C=0 calls the original certified augmented route.
    S zero_H(1u,zero);
    S zero_C(1u,zero);
    const Matrix<X> coupling(1u,1u,X(1u,dp));
    const auto augmented=validated_saddle_point_solve(
        zero_H,coupling,zero_C,Vector<X>({X(4u,dp)}),
        Vector<X>({X(2u,dp)}));
    ARIADNE_TEST_ASSERT(augmented.status==ValidatedSaddlePointStatus::SOLVED);
    if(augmented.primary.has_value() && augmented.secondary.has_value()) {
        check_known_saddle_point_solution(augmented,
            Vector<X>({X(2u,dp)}),Vector<X>({X(4u,dp)}),
            "explicit C=0");
    }

    S empty_H(0u,zero);
    Matrix<X> empty_A(1u,0u,zero);
    S lower_only(1u,zero);
    lower_only[0u][0u]=X(2u,dp);
    const auto no_leading=validated_saddle_point_solve(
        empty_H,empty_A,lower_only,Vector<X>(0u,zero),
        Vector<X>({X(-4,dp)}));
    ARIADNE_TEST_ASSERT(no_leading.status==ValidatedSaddlePointStatus::SOLVED);
    if(no_leading.primary.has_value() && no_leading.secondary.has_value()) {
        check_known_saddle_point_solution(no_leading,
            Vector<X>(0u,zero),Vector<X>({X(2u,dp)}),
            "empty leading block");
    }

    const auto no_rhs=validated_saddle_point_solve(
        H,A,C,Matrix<X>(2u,0u,zero),Matrix<X>(1u,0u,zero));
    ARIADNE_TEST_ASSERT(no_rhs.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(no_rhs.primary.has_value());
    ARIADNE_TEST_ASSERT(no_rhs.secondary.has_value());
    if(no_rhs.primary.has_value() && no_rhs.secondary.has_value()) {
        ARIADNE_TEST_EQUALS(no_rhs.primary->column_size(),0u);
        ARIADNE_TEST_EQUALS(no_rhs.secondary->column_size(),0u);
    }
}

Void test_generalized_saddle_point_failures()
{
    const X zero(0u,dp);
    S H(1u,zero);
    H[0u][0u]=X(2u,dp);
    const Matrix<X> A(1u,1u,X(1u,dp));
    S C(1u,zero);
    C[0u][0u]=X(1u,dp);
    const Vector<X> f({X(2u,dp)});
    const Vector<X> g({X(1u,dp)});
    const Matrix<X> F(1u,1u,X(2u,dp));
    const Matrix<X> G(1u,1u,X(1u,dp));
    S bad_C(2u,zero);
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,bad_C,f,g));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,bad_C,F,G));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,Matrix<X>(1u,2u,zero),C,f,g));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,C,Vector<X>(2u,zero),g));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,C,f,Vector<X>(2u,zero)));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,Matrix<X>(1u,2u,zero),C,F,G));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,C,Matrix<X>(2u,1u,zero),G));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,C,F,Matrix<X>(2u,1u,zero)));
    ARIADNE_TEST_FAIL(validated_saddle_point_solve(H,A,C,F,Matrix<X>(1u,2u,zero)));

    S negative(1u,zero);
    negative[0u][0u]=X(-1,dp);
    const auto no_augmentation=validated_saddle_point_solve(negative,A,C,f,g);
    ARIADNE_TEST_ASSERT(no_augmentation.status==
        ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);

    S uncertain(1u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    const auto unknown=validated_saddle_point_solve(uncertain,A,C,f,g);
    ARIADNE_TEST_ASSERT(unknown.status==ValidatedSaddlePointStatus::INDETERMINATE);

    // An indefinite Schur complement is inconclusive even if A has dependent rows.
    Matrix<X> dependent(2u,1u,X(1u,dp));
    S negative_C(2u,zero);
    negative_C[0u][0u]=X(-1,dp);
    const auto schur_failure=validated_saddle_point_solve(
        H,dependent,negative_C,F,Matrix<X>(2u,1u,zero));
    ARIADNE_TEST_ASSERT(schur_failure.status==ValidatedSaddlePointStatus::INDETERMINATE);

    // The zero lower endpoint does NOT make C identically zero:
    // [0,1] exercises the upper-endpoint comparison in the structural
    // zero test. At C=1 the system (H=2,A=1,f=2,g=1) has x=1,y=0.
    S one_sided_C(1u,zero);
    one_sided_C[0u][0u]=X(0,1,dp);
    const auto one_sided_solution=validated_saddle_point_solve(
        H,A,one_sided_C,f,g);
    ARIADNE_TEST_ASSERT(one_sided_solution.status==ValidatedSaddlePointStatus::SOLVED);
    ARIADNE_TEST_ASSERT(one_sided_solution.primary.has_value());
    ARIADNE_TEST_ASSERT(one_sided_solution.secondary.has_value());
    if(one_sided_solution.primary.has_value() && one_sided_solution.secondary.has_value()) {
        check_known_saddle_point_solution(
            one_sided_solution,Vector<X>({X(1u,dp)}),
            Vector<X>({zero}),"C=[0,1], same lower endpoint as zero");
        InclusionAudit audit("generalized saddle-point with C=[0,1]");
        const Vector<X> upper=static_cast<Matrix<X>>(H)*(*one_sided_solution.primary)+
                              transpose(A)*(*one_sided_solution.secondary);
        const Vector<X> lower=A*(*one_sided_solution.primary)-
                              static_cast<Matrix<X>>(one_sided_C)*(*one_sided_solution.secondary);
        audit.add(upper[0u],f[0u],"upper[0]");
        audit.add(lower[0u],g[0u],"lower[0]");
        const bool one_sided_interval_residual_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(one_sided_interval_residual_enclosed);
    }

    // If H cannot be certified SPD, a C=[0,1] interval still
    // forbids the C=0 augmentation: some represented C are nonzero.
    const auto one_sided_rejected=validated_saddle_point_solve(
        negative,A,one_sided_C,f,g);
    ARIADNE_TEST_ASSERT(one_sided_rejected.status==
        ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!one_sided_rejected.primary.has_value());
    ARIADNE_TEST_ASSERT(!one_sided_rejected.secondary.has_value());

    // A zero-containing interval is not a structural zero: no fallback.
    S uncertain_C(1u,zero);
    uncertain_C[0u][0u]=X(-1,1,dp);
    const auto not_structurally_zero=validated_saddle_point_solve(
        negative,A,uncertain_C,f,g);
    ARIADNE_TEST_ASSERT(not_structurally_zero.status==
        ValidatedSaddlePointStatus::LEADING_BLOCK_NOT_POSITIVE_DEFINITE);
}

Void test_empty_and_zero_rhs()
{
    const X zero(0u,dp);
    const S empty(0u,zero);
    const Matrix<X> no_rhs(0u,0u,zero);
    const auto empty_matrix_solution=validated_positive_definite_solve(empty,no_rhs);
    ARIADNE_TEST_ASSERT(empty_matrix_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(empty_matrix_solution.solution.has_value());
    check_matrix_solution(empty,no_rhs,*empty_matrix_solution.solution);

    // Empty equation system with two separate empty right-hand sides.
    const Matrix<X> two_empty_rhs(0u,2u,zero);
    const auto two_empty_solution=validated_positive_definite_solve(empty,two_empty_rhs);
    ARIADNE_TEST_ASSERT(two_empty_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(two_empty_solution.solution.has_value());
    check_matrix_solution(empty,two_empty_rhs,*two_empty_solution.solution);

    const Vector<X> empty_rhs(0u,zero);
    const auto empty_vector_solution=validated_positive_definite_solve(empty,empty_rhs);
    ARIADNE_TEST_ASSERT(empty_vector_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(empty_vector_solution.solution.has_value());
    check_vector_solution(empty,empty_rhs,*empty_vector_solution.solution);

    const S positive(2u,X(0u,dp));
    // A non-SPD matrix remains non-SPD even when there are no right-hand sides.
    const Matrix<X> absent_rhs(2u,0u,zero);
    const auto absent_rhs_solution=validated_positive_definite_solve(positive,absent_rhs);
    ARIADNE_TEST_ASSERT(absent_rhs_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!absent_rhs_solution.solution.has_value());

    S identity(2u,zero);
    identity[0u][0u]=X(1u,dp);
    identity[1u][1u]=X(1u,dp);
    const auto zero_column_solution=validated_positive_definite_solve(identity,absent_rhs);
    ARIADNE_TEST_ASSERT(zero_column_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(zero_column_solution.solution.has_value());
    check_matrix_solution(identity,absent_rhs,*zero_column_solution.solution);
}

Void test_scalar_and_dense_spd()
{
    const X zero(0u,dp);

    S scalar(1u,zero);
    scalar[0u][0u]=X(4u,dp);
    const Vector<X> b1({X(8u,dp)});
    const auto scalar_solution=validated_positive_definite_solve(scalar,b1);
    ARIADNE_TEST_ASSERT(scalar_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(scalar_solution.solution.has_value());
    check_contains((*scalar_solution.solution)[0u],X(2u,dp));
    check_vector_solution(scalar,b1,*scalar_solution.solution);

    // The Cholesky factor is irrational; solution must still enclose 1.
    S irrational(1u,zero);
    irrational[0u][0u]=X(2u,dp);
    const Vector<X> b2({X(2u,dp)});
    const auto irrational_solution=validated_positive_definite_solve(irrational,b2);
    ARIADNE_TEST_ASSERT(irrational_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(irrational_solution.solution.has_value());
    check_contains((*irrational_solution.solution)[0u],X(1u,dp));
    check_vector_solution(irrational,b2,*irrational_solution.solution);

    // Two nontrivial columns; the exact solutions are (1,1) and (2,0).
    S matrix(2u,zero);
    matrix[0u][0u]=X(2u,dp);
    matrix[0u][1u]=X(1u,dp);
    matrix[1u][1u]=X(2u,dp);
    Matrix<X> right_hand_sides(2u,2u,zero);
    right_hand_sides[0u][0u]=X(3u,dp);
    right_hand_sides[1u][0u]=X(3u,dp);
    right_hand_sides[0u][1u]=X(4u,dp);
    right_hand_sides[1u][1u]=X(2u,dp);

    const auto matrix_solution=validated_positive_definite_solve(matrix,right_hand_sides);
    ARIADNE_TEST_ASSERT(matrix_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(matrix_solution.solution.has_value());
    check_matrix_solution(matrix,right_hand_sides,*matrix_solution.solution);
    check_contains((*matrix_solution.solution)[0u][0u],X(1u,dp));
    check_contains((*matrix_solution.solution)[1u][0u],X(1u,dp));
    check_contains((*matrix_solution.solution)[0u][1u],X(2u,dp));
    check_contains((*matrix_solution.solution)[1u][1u],zero);

    const Vector<X> b({X(3u,dp),X(3u,dp)});
    const auto vector_solution=validated_positive_definite_solve(matrix,b);
    ARIADNE_TEST_ASSERT(vector_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(vector_solution.solution.has_value());
    check_vector_solution(matrix,b,*vector_solution.solution);
    check_contains((*vector_solution.solution)[0u],X(1u,dp));
    check_contains((*vector_solution.solution)[1u],X(1u,dp));
}

Void test_three_by_three_and_interval_family()
{
    const X zero(0u,dp);
    S matrix(3u,zero);
    for(SizeType i=0; i!=3u; ++i) {
        matrix[i][i]=X(4u,dp);
    }
    matrix[0u][1u]=X(1u,dp);
    matrix[0u][2u]=X(1u,dp);
    matrix[1u][2u]=X(1u,dp);

    const Vector<X> rhs({X(9u,dp),X(12u,dp),X(15u,dp)});
    const auto dense_solution=validated_positive_definite_solve(matrix,rhs);
    ARIADNE_TEST_ASSERT(dense_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(dense_solution.solution.has_value());
    check_vector_solution(matrix,rhs,*dense_solution.solution);
    for(SizeType i=0; i!=3u; ++i) {
        check_contains((*dense_solution.solution)[i],X(static_cast<Int>(i+1u),dp));
    }

    S interval_matrix(2u,zero);
    interval_matrix[0u][0u]=X(4,5,dp);
    interval_matrix[0u][1u]=X(-1,1,dp);
    interval_matrix[1u][1u]=X(4,5,dp);
    Vector<X> interval_rhs(2u,zero);
    interval_rhs[0u]=X(4,5,dp);
    interval_rhs[1u]=X(8,10,dp);

    const auto interval_solution=validated_positive_definite_solve(interval_matrix,interval_rhs);
    ARIADNE_TEST_ASSERT(interval_solution.status==Status::POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(interval_solution.solution.has_value());
    check_vector_solution(interval_matrix,interval_rhs,*interval_solution.solution);

    // Choose one exact member of the family: A=diag(4,4), b=(4,8).
    check_contains((*interval_solution.solution)[0u],X(1u,dp));
    check_contains((*interval_solution.solution)[1u],X(2u,dp));
}

Void test_failure_and_indeterminate()
{
    const X zero(0u,dp);
    S negative(1u,zero);
    negative[0u][0u]=X(-1,dp);
    const Vector<X> b1({X(1u,dp)});
    const auto negative_vector_solution=validated_positive_definite_solve(negative,b1);
    ARIADNE_TEST_ASSERT(negative_vector_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_vector_solution.solution.has_value());

    const Matrix<X> B1(1u,1u,X(1u,dp));
    const auto negative_matrix_solution=validated_positive_definite_solve(negative,B1);
    ARIADNE_TEST_ASSERT(negative_matrix_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_matrix_solution.solution.has_value());

    S uncertain(1u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    const auto uncertain_vector_solution=validated_positive_definite_solve(uncertain,b1);
    ARIADNE_TEST_ASSERT(uncertain_vector_solution.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_vector_solution.solution.has_value());

    const auto uncertain_matrix_solution=validated_positive_definite_solve(uncertain,B1);
    ARIADNE_TEST_ASSERT(uncertain_matrix_solution.status==Status::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_matrix_solution.solution.has_value());

    // A singular positive-semidefinite matrix cannot use an SPD solve.
    S semidefinite(2u,zero);
    semidefinite[0u][0u]=X(1u,dp);
    semidefinite[0u][1u]=X(1u,dp);
    semidefinite[1u][1u]=X(1u,dp);
    const Vector<X> b2({X(2u,dp),X(2u,dp)});
    const auto semidefinite_solution=validated_positive_definite_solve(semidefinite,b2);
    ARIADNE_TEST_ASSERT(semidefinite_solution.status==Status::NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!semidefinite_solution.solution.has_value());

    const Matrix<X> bad_matrix_rhs(2u,1u,zero);
    const Vector<X> bad_vector_rhs(2u,zero);
    ARIADNE_TEST_FAIL(validated_positive_definite_solve(negative,bad_matrix_rhs));
    ARIADNE_TEST_FAIL(validated_positive_definite_solve(negative,bad_vector_rhs));
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_empty_and_zero_rhs());
    ARIADNE_TEST_CALL(test_scalar_and_dense_spd());
    ARIADNE_TEST_CALL(test_three_by_three_and_interval_family());
    ARIADNE_TEST_CALL(test_failure_and_indeterminate());
    ARIADNE_TEST_CALL(test_validated_saddle_point_solve());
    ARIADNE_TEST_CALL(test_validated_saddle_point_matrix_rhs());
    ARIADNE_TEST_CALL(test_saddle_point_augmented_leading_block());
    ARIADNE_TEST_CALL(test_generalized_saddle_point());
    ARIADNE_TEST_CALL(test_generalized_saddle_point_failures());
    return ARIADNE_TEST_FAILURES;
}
