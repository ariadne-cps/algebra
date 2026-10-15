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
// The source line identifies which expected primal/dual component failed.
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

Void check_kkt_residuals(
    const S& H, const Matrix<X>& A, const Vector<X>& f, const Vector<X>& g,
    const Vector<X>& x, const Vector<X>& lambda)
{
    ARIADNE_TEST_EQUALS(x.size(),H.size());
    ARIADNE_TEST_EQUALS(lambda.size(),A.row_size());
    const Matrix<X> dense_H=static_cast<Matrix<X>>(H);
    const Vector<X> stationary=dense_H*x+transpose(A)*lambda;
    const Vector<X> constrained=A*x;
    InclusionAudit primal_audit("KKT stationarity: f contained in H*x+A^T*lambda");
    for(SizeType i=0; i!=f.size(); ++i) {
        primal_audit.add(stationary[i],f[i],"row "+std::to_string(i));
    }
    const bool kkt_stationarity_enclosed=primal_audit.valid();
    ARIADNE_TEST_ASSERT(kkt_stationarity_enclosed);

    InclusionAudit equality_audit("KKT equalities: g contained in A*x");
    for(SizeType i=0; i!=g.size(); ++i) {
        equality_audit.add(constrained[i],g[i],"row "+std::to_string(i));
    }
    const bool kkt_equalities_enclosed=equality_audit.valid();
    ARIADNE_TEST_ASSERT(kkt_equalities_enclosed);
}

// Check both primal and multiplier components against a selected exact
// member of an interval KKT family, with a meaningful scenario label.
Void check_known_kkt_solution(
    const ValidatedKktResult& solution,
    const Vector<X>& expected_primal,
    const Vector<X>& expected_multipliers,
    const std::string& scenario)
{
    ARIADNE_TEST_ASSERT(solution.primal.has_value());
    ARIADNE_TEST_ASSERT(solution.multipliers.has_value());
    if(!solution.primal.has_value() || !solution.multipliers.has_value()) {
        return;
    }
    ARIADNE_TEST_EQUALS(solution.primal->size(),expected_primal.size());
    ARIADNE_TEST_EQUALS(solution.multipliers->size(),expected_multipliers.size());
    InclusionAudit audit("KKT known solution: "+scenario);
    for(SizeType i=0; i!=expected_primal.size(); ++i) {
        audit.add((*solution.primal)[i],expected_primal[i],
                  "primal["+std::to_string(i)+"]");
    }
    for(SizeType i=0; i!=expected_multipliers.size(); ++i) {
        audit.add((*solution.multipliers)[i],expected_multipliers[i],
                  "multiplier["+std::to_string(i)+"]");
    }
    const bool kkt_reference_solution_enclosed=audit.valid();
    ARIADNE_TEST_ASSERT(kkt_reference_solution_enclosed);
}

Void test_validated_kkt_solve()
{
    const X zero(0u,dp);
    const Vector<X> empty_vector(0u,zero);

    // Both the zero-dimensional KKT system and the unconstrained
    // case have unique solutions under positive definiteness.
    const S empty_hessian(0u,zero);
    const Matrix<X> empty_constraints(0u,0u,zero);
    const auto empty_solution=validated_kkt_solve(
        empty_hessian,empty_constraints,empty_vector,empty_vector);
    ARIADNE_TEST_ASSERT(empty_solution.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(empty_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(empty_solution.multipliers.has_value());
    check_kkt_residuals(empty_hessian,empty_constraints,empty_vector,
                        empty_vector,empty_solution.primal.value(),empty_solution.multipliers.value());

    S scalar_hessian(1u,zero);
    scalar_hessian[0u][0u]=X(2u,dp);
    const Matrix<X> no_constraints(0u,1u,zero);
    const Vector<X> unconstrained_rhs({X(6u,dp)});
    const auto unconstrained_solution=validated_kkt_solve(
        scalar_hessian,no_constraints,unconstrained_rhs,empty_vector);
    ARIADNE_TEST_ASSERT(unconstrained_solution.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(unconstrained_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(unconstrained_solution.multipliers.has_value());
    check_known_kkt_solution(unconstrained_solution,Vector<X>({X(3u,dp)}),
                             empty_vector,"unconstrained SPD system");
    check_kkt_residuals(scalar_hessian,no_constraints,unconstrained_rhs,
                        empty_vector,unconstrained_solution.primal.value(),
                        unconstrained_solution.multipliers.value());

    // One equality: H=2I, A=(1,1), f=0, g=2.
    // Exact solution x=(1,1), lambda=-2.
    S hessian(2u,zero);
    hessian[0u][0u]=X(2u,dp);
    hessian[1u][1u]=X(2u,dp);
    Matrix<X> constraints(1u,2u,zero);
    constraints[0u][0u]=X(1u,dp);
    constraints[0u][1u]=X(1u,dp);
    const Vector<X> f(2u,zero);
    const Vector<X> g({X(2u,dp)});
    const auto single_solution=validated_kkt_solve(hessian,constraints,f,g);
    ARIADNE_TEST_ASSERT(single_solution.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(single_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(single_solution.multipliers.has_value());
    check_known_kkt_solution(single_solution,
                             Vector<X>({X(1u,dp),X(1u,dp)}),
                             Vector<X>({X(-2,dp)}),"one equality");
    check_kkt_residuals(hessian,constraints,f,g,
                        single_solution.primal.value(),single_solution.multipliers.value());

    // A non-diagonal Hessian and a nonzero dual variable.
    S coupled_hessian(2u,zero);
    coupled_hessian[0u][0u]=X(2u,dp);
    coupled_hessian[0u][1u]=X(1u,dp);
    coupled_hessian[1u][1u]=X(2u,dp);
    Matrix<X> first_coordinate(1u,2u,zero);
    first_coordinate[0u][0u]=X(1u,dp);
    const Vector<X> coupled_f({X(11u,dp),X(8u,dp)});
    const Vector<X> coupled_g({X(2u,dp)});
    const auto coupled_solution=validated_kkt_solve(
        coupled_hessian,first_coordinate,coupled_f,coupled_g);
    ARIADNE_TEST_ASSERT(coupled_solution.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(coupled_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(coupled_solution.multipliers.has_value());
    check_known_kkt_solution(coupled_solution,
                             Vector<X>({X(2u,dp),X(3u,dp)}),
                             Vector<X>({X(4u,dp)}),"coupled Hessian");
    check_kkt_residuals(coupled_hessian,first_coordinate,coupled_f,coupled_g,
                        coupled_solution.primal.value(),coupled_solution.multipliers.value());

    // Two independent constraints: exact x=(1,2,3), lambda=(4,5).
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
    const auto multi_solution=validated_kkt_solve(h3,a2,f3,g2);
    ARIADNE_TEST_ASSERT(multi_solution.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(multi_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(multi_solution.multipliers.has_value());
    check_kkt_residuals(h3,a2,f3,g2,
                        multi_solution.primal.value(),multi_solution.multipliers.value());
    check_known_kkt_solution(multi_solution,
                             Vector<X>({X(1u,dp),X(2u,dp),X(3u,dp)}),
                             Vector<X>({X(4u,dp),X(5u,dp)}),"two equalities");

    // A family of admissible KKT systems with a known exact member.
    S interval_hessian(2u,zero);
    interval_hessian[0u][0u]=X(4,5,dp);
    interval_hessian[1u][1u]=X(4,5,dp);
    Matrix<X> interval_constraints(1u,2u,zero);
    interval_constraints[0u][0u]=X(1,2,dp);
    const Vector<X> interval_f({X(4,5,dp),X(8,10,dp)});
    const Vector<X> interval_g({X(1,2,dp)});
    const auto interval_solution=validated_kkt_solve(
        interval_hessian,interval_constraints,interval_f,interval_g);
    ARIADNE_TEST_ASSERT(interval_solution.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(interval_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(interval_solution.multipliers.has_value());
    check_kkt_residuals(interval_hessian,interval_constraints,interval_f,interval_g,
                        interval_solution.primal.value(),interval_solution.multipliers.value());
    check_known_kkt_solution(interval_solution,
                             Vector<X>({X(1u,dp),X(2u,dp)}),
                             Vector<X>({zero}),"interval KKT family");

    // A provably non-SPD Hessian is a failure of this solver's precondition.
    S negative(1u,zero);
    negative[0u][0u]=X(-1,dp);
    Matrix<X> scalar_constraints(1u,1u,X(1u,dp));
    const Vector<X> scalar_rhs({X(1u,dp)});
    const auto negative_solution=validated_kkt_solve(
        negative,scalar_constraints,scalar_rhs,scalar_rhs);
    ARIADNE_TEST_ASSERT(
        negative_solution.status==ValidatedKktStatus::HESSIAN_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!negative_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(!negative_solution.multipliers.has_value());

    S uncertain_hessian(1u,zero);
    uncertain_hessian[0u][0u]=X(-1,1,dp);
    const auto uncertain_h_solution=validated_kkt_solve(
        uncertain_hessian,scalar_constraints,scalar_rhs,scalar_rhs);
    ARIADNE_TEST_ASSERT(uncertain_h_solution.status==ValidatedKktStatus::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!uncertain_h_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(!uncertain_h_solution.multipliers.has_value());

    // More constraint rows than unknowns is universally rank deficient.
    Matrix<X> overconstrained(2u,1u,X(1u,dp));
    const Vector<X> two_constraints({X(1u,dp),X(1u,dp)});
    const auto dependent_solution=validated_kkt_solve(
        scalar_hessian,overconstrained,scalar_rhs,two_constraints);
    ARIADNE_TEST_ASSERT(
        dependent_solution.status==ValidatedKktStatus::DEPENDENT_CONSTRAINTS);
    ARIADNE_TEST_ASSERT(!dependent_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(!dependent_solution.multipliers.has_value());

    // Repeated constraints are dependent even without dimension excess.
    Matrix<X> repeated(2u,2u,zero);
    repeated[0u][0u]=X(1u,dp);
    repeated[1u][0u]=X(1u,dp);
    const auto repeated_solution=validated_kkt_solve(
        hessian,repeated,f,two_constraints);
    ARIADNE_TEST_ASSERT(
        repeated_solution.status==ValidatedKktStatus::DEPENDENT_CONSTRAINTS);
    ARIADNE_TEST_ASSERT(!repeated_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(!repeated_solution.multipliers.has_value());

    // Interval dependence of the same A entry is unresolved. Never use
    // a midpoint-rank decision, and do not infer dependence from an
    // interval Schur matrix which overestimates correlations.
    Matrix<X> ambiguous(1u,1u,X(-1,1,dp));
    const auto ambiguous_solution=validated_kkt_solve(
        scalar_hessian,ambiguous,scalar_rhs,scalar_rhs);
    ARIADNE_TEST_ASSERT(ambiguous_solution.status==ValidatedKktStatus::INDETERMINATE);
    ARIADNE_TEST_ASSERT(!ambiguous_solution.primal.has_value());
    ARIADNE_TEST_ASSERT(!ambiguous_solution.multipliers.has_value());

    Matrix<X> bad_columns(1u,2u,zero);
    const Vector<X> bad_f(2u,zero);
    const Vector<X> bad_g(2u,zero);
    ARIADNE_TEST_FAIL(validated_kkt_solve(scalar_hessian,bad_columns,scalar_rhs,scalar_rhs));
    ARIADNE_TEST_FAIL(validated_kkt_solve(scalar_hessian,scalar_constraints,bad_f,scalar_rhs));
    ARIADNE_TEST_FAIL(validated_kkt_solve(scalar_hessian,scalar_constraints,scalar_rhs,bad_g));
}


Void test_validated_kkt_matrix_rhs()
{
    const X zero(0u,dp);
    S H(2u,zero);
    H[0u][0u]=X(2u,dp);
    H[1u][1u]=X(2u,dp);
    Matrix<X> A(1u,2u,zero);
    A[0u][0u]=X(1u,dp);
    A[0u][1u]=X(1u,dp);

    // First column: x=(1,1), lambda=-2; second: x=(2,3), lambda=4.
    Matrix<X> F(2u,2u,zero);
    F[0u][1u]=X(8u,dp);
    F[1u][1u]=X(10u,dp);
    Matrix<X> G(1u,2u,zero);
    G[0u][0u]=X(2u,dp);
    G[0u][1u]=X(5u,dp);

    const auto multiple=validated_kkt_solve(H,A,F,G);
    ARIADNE_TEST_ASSERT(multiple.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(multiple.primal.has_value());
    ARIADNE_TEST_ASSERT(multiple.multipliers.has_value());
    if(multiple.primal.has_value() && multiple.multipliers.has_value()) {
        const Matrix<X>& X_solution=*multiple.primal;
        const Matrix<X>& Lambda=*multiple.multipliers;
        ARIADNE_TEST_EQUALS(X_solution.row_size(),2u);
        ARIADNE_TEST_EQUALS(X_solution.column_size(),2u);
        ARIADNE_TEST_EQUALS(Lambda.row_size(),1u);
        ARIADNE_TEST_EQUALS(Lambda.column_size(),2u);
        const Matrix<X> residual=static_cast<Matrix<X>>(H)*X_solution+
                                 transpose(A)*Lambda;
        const Matrix<X> equality=A*X_solution;
        InclusionAudit audit("KKT two-right-hand-side residual");
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
        audit.add(X_solution[0u][0u],X(1u,dp),"primal[0,0]");
        audit.add(X_solution[1u][0u],X(1u,dp),"primal[1,0]");
        audit.add(Lambda[0u][0u],X(-2,dp),"lambda[0,0]");
        audit.add(X_solution[0u][1u],X(2u,dp),"primal[0,1]");
        audit.add(X_solution[1u][1u],X(3u,dp),"primal[1,1]");
        audit.add(Lambda[0u][1u],X(4u,dp),"lambda[0,1]");
        const bool both_kkt_columns_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(both_kkt_columns_enclosed);
    }

    // A zero-column right-hand side must still certify the shared
    // KKT operator and produce correctly shaped, empty solutions.
    Matrix<X> empty_F(2u,0u,zero);
    Matrix<X> empty_G(1u,0u,zero);
    const auto no_rhs=validated_kkt_solve(H,A,empty_F,empty_G);
    ARIADNE_TEST_ASSERT(no_rhs.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(no_rhs.primal.has_value());
    ARIADNE_TEST_ASSERT(no_rhs.multipliers.has_value());
    if(no_rhs.primal.has_value() && no_rhs.multipliers.has_value()) {
        ARIADNE_TEST_EQUALS(no_rhs.primal->row_size(),2u);
        ARIADNE_TEST_EQUALS(no_rhs.primal->column_size(),0u);
        ARIADNE_TEST_EQUALS(no_rhs.multipliers->row_size(),1u);
        ARIADNE_TEST_EQUALS(no_rhs.multipliers->column_size(),0u);
    }

    const Matrix<X> wrong_F(1u,2u,zero);
    const Matrix<X> wrong_G_rows(2u,2u,zero);
    const Matrix<X> wrong_G_cols(1u,1u,zero);
    const Matrix<X> wrong_A(1u,3u,zero);
    ARIADNE_TEST_FAIL(validated_kkt_solve(H,wrong_A,F,G));
    ARIADNE_TEST_FAIL(validated_kkt_solve(H,A,wrong_F,G));
    ARIADNE_TEST_FAIL(validated_kkt_solve(H,A,F,wrong_G_rows));
    ARIADNE_TEST_FAIL(validated_kkt_solve(H,A,F,wrong_G_cols));

    // Verify failure classification for the matrix overload as well.
    S nonspd(2u,zero);
    nonspd[0u][0u]=X(-1,dp);
    nonspd[1u][1u]=X(1u,dp);
    const auto rejected_h=validated_kkt_solve(nonspd,A,F,G);
    ARIADNE_TEST_ASSERT(rejected_h.status==ValidatedKktStatus::HESSIAN_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!rejected_h.primal.has_value());
    ARIADNE_TEST_ASSERT(!rejected_h.multipliers.has_value());

    S uncertain(2u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    uncertain[1u][1u]=X(1u,dp);
    const auto rejected_uncertain=validated_kkt_solve(uncertain,A,F,G);
    ARIADNE_TEST_ASSERT(rejected_uncertain.status==ValidatedKktStatus::INDETERMINATE);

    Matrix<X> dependent_A(2u,2u,zero);
    dependent_A[0u][0u]=X(1u,dp);
    dependent_A[1u][0u]=X(1u,dp);
    const Matrix<X> G2(2u,2u,zero);
    const auto rejected_rank=validated_kkt_solve(H,dependent_A,F,G2);
    ARIADNE_TEST_ASSERT(rejected_rank.status==ValidatedKktStatus::DEPENDENT_CONSTRAINTS);

    Matrix<X> ambiguous_A(1u,2u,zero);
    ambiguous_A[0u][0u]=X(-1,1,dp);
    const auto rejected_ambiguous=validated_kkt_solve(H,ambiguous_A,F,G);
    ARIADNE_TEST_ASSERT(rejected_ambiguous.status==ValidatedKktStatus::INDETERMINATE);

    // Zero-dimensional operator with multiple columns.
    S h0(0u,zero);
    const Matrix<X> a0(0u,0u,zero);
    const Matrix<X> f0(0u,2u,zero);
    const Matrix<X> g0(0u,2u,zero);
    const auto zero_system=validated_kkt_solve(h0,a0,f0,g0);
    ARIADNE_TEST_ASSERT(zero_system.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(zero_system.primal.has_value());
    ARIADNE_TEST_ASSERT(zero_system.multipliers.has_value());
    if(zero_system.primal.has_value() && zero_system.multipliers.has_value()) {
        ARIADNE_TEST_EQUALS(zero_system.primal->column_size(),2u);
        ARIADNE_TEST_EQUALS(zero_system.multipliers->column_size(),2u);
    }

    // Unconstrained multi-RHS: the same SPD solve is used with m=0.
    const Matrix<X> no_A(0u,2u,zero);
    const Matrix<X> no_G(0u,2u,zero);
    const auto unconstrained=validated_kkt_solve(H,no_A,F,no_G);
    ARIADNE_TEST_ASSERT(unconstrained.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(unconstrained.primal.has_value());
    ARIADNE_TEST_ASSERT(unconstrained.multipliers.has_value());
    if(unconstrained.primal.has_value() && unconstrained.multipliers.has_value()) {
        ARIADNE_TEST_EQUALS(unconstrained.primal->column_size(),2u);
        ARIADNE_TEST_EQUALS(unconstrained.multipliers->row_size(),0u);
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
    const auto interval_system=validated_kkt_solve(
        interval_H,interval_A,interval_F,interval_G);
    ARIADNE_TEST_ASSERT(interval_system.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(interval_system.primal.has_value());
    ARIADNE_TEST_ASSERT(interval_system.multipliers.has_value());
    if(interval_system.primal.has_value() && interval_system.multipliers.has_value()) {
        InclusionAudit audit("KKT interval family, two right-hand sides");
        const Matrix<X> exact_X=static_cast<Matrix<X>>(interval_H)*(*interval_system.primal)+
                                 transpose(interval_A)*(*interval_system.multipliers);
        const Matrix<X> exact_G=interval_A*(*interval_system.primal);
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


// A semidefinite Hessian can still define a strictly convex quadratic
// on the affine feasible space when H is positive on ker(A). Test the
// augmented route without weakening any of the original SPD tests.
Void test_kkt_augmented_hessian()
{
    const X zero(0u,dp);

    // H=diag(0,2) is positive semidefinite but not SPD.
    // A fixes the first coordinate. The exact solution is
    // x=(2,3), lambda=5, for f=(5,6) and g=2.
    S semidefinite(2u,zero);
    semidefinite[1u][1u]=X(2u,dp);
    Matrix<X> A(1u,2u,zero);
    A[0u][0u]=X(1u,dp);
    const Vector<X> f({X(5u,dp),X(6u,dp)});
    const Vector<X> g({X(2u,dp)});
    const auto constrained=validated_kkt_solve(semidefinite,A,f,g);
    ARIADNE_TEST_ASSERT(constrained.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(constrained.primal.has_value());
    ARIADNE_TEST_ASSERT(constrained.multipliers.has_value());
    if(constrained.primal.has_value() && constrained.multipliers.has_value()) {
        check_known_kkt_solution(constrained,
                                 Vector<X>({X(2u,dp),X(3u,dp)}),
                                 Vector<X>({X(5u,dp)}),
                                 "positive semidefinite H, constrained nullspace");
        check_kkt_residuals(semidefinite,A,f,g,
                            *constrained.primal,*constrained.multipliers);
    }

    // The original H=[-1,1] is not certifiably SPD, but A=[2]
    // makes H+A^T*A=[3,5] certified SPD for every interval member.
    // The member H=0 has x=1 and lambda=2.
    S uncertain(1u,zero);
    uncertain[0u][0u]=X(-1,1,dp);
    Matrix<X> A_scalar(1u,1u,X(2u,dp));
    const Vector<X> f_scalar({X(4u,dp)});
    const Vector<X> g_scalar({X(2u,dp)});
    const auto uncertain_constrained=validated_kkt_solve(
        uncertain,A_scalar,f_scalar,g_scalar);
    ARIADNE_TEST_ASSERT(uncertain_constrained.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(uncertain_constrained.primal.has_value());
    ARIADNE_TEST_ASSERT(uncertain_constrained.multipliers.has_value());
    if(uncertain_constrained.primal.has_value()
       && uncertain_constrained.multipliers.has_value()) {
        check_known_kkt_solution(uncertain_constrained,
                                 Vector<X>({X(1u,dp)}),
                                 Vector<X>({X(2u,dp)}),
                                 "interval H with constrained negative curvature");
        check_kkt_residuals(uncertain,A_scalar,f_scalar,g_scalar,
                            *uncertain_constrained.primal,
                            *uncertain_constrained.multipliers);
    }

    // Even H=-1 is admissible when all negative-curvature directions
    // are fixed by A=[2]. The multiplier in the original KKT system
    // is unchanged by the augmentation.
    S negative_scalar(1u,zero);
    negative_scalar[0u][0u]=X(-1,dp);
    const Vector<X> negative_f({X(3u,dp)});
    const auto negative_constrained=validated_kkt_solve(
        negative_scalar,A_scalar,negative_f,g_scalar);
    ARIADNE_TEST_ASSERT(negative_constrained.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(negative_constrained.primal.has_value());
    ARIADNE_TEST_ASSERT(negative_constrained.multipliers.has_value());
    if(negative_constrained.primal.has_value()
       && negative_constrained.multipliers.has_value()) {
        check_known_kkt_solution(negative_constrained,
                                 Vector<X>({X(1u,dp)}),
                                 Vector<X>({X(2u,dp)}),
                                 "negative Hessian, constrained degree of freedom");
        check_kkt_residuals(negative_scalar,A_scalar,negative_f,g_scalar,
                            *negative_constrained.primal,
                            *negative_constrained.multipliers);
    }

    // H=0 and full-rank A=I: the constraints determine every primal
    // variable and lambda=F. Check two right-hand sides at once.
    S zero_hessian(2u,zero);
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
    const auto zero_full_rank=validated_kkt_solve(
        zero_hessian,identity,F,G);
    ARIADNE_TEST_ASSERT(zero_full_rank.status==ValidatedKktStatus::SOLVED);
    ARIADNE_TEST_ASSERT(zero_full_rank.primal.has_value());
    ARIADNE_TEST_ASSERT(zero_full_rank.multipliers.has_value());
    if(zero_full_rank.primal.has_value()
       && zero_full_rank.multipliers.has_value()) {
        InclusionAudit audit("zero Hessian, full-rank constraints, two RHS");
        const Matrix<X>& x=*zero_full_rank.primal;
        const Matrix<X>& lambda=*zero_full_rank.multipliers;
        const Matrix<X> stationarity=static_cast<Matrix<X>>(zero_hessian)*x+
                                     transpose(identity)*lambda;
        const Matrix<X> equalities=identity*x;
        for(SizeType i=0; i!=2u; ++i) {
            for(SizeType col=0; col!=2u; ++col) {
                audit.add(x[i][col],G[i][col],
                          "primal["+std::to_string(i)+","+std::to_string(col)+"]");
                audit.add(lambda[i][col],F[i][col],
                          "multiplier["+std::to_string(i)+","+std::to_string(col)+"]");
                audit.add(stationarity[i][col],F[i][col],
                          "stationarity["+std::to_string(i)+","+std::to_string(col)+"]");
                audit.add(equalities[i][col],G[i][col],
                          "equality["+std::to_string(i)+","+std::to_string(col)+"]");
            }
        }
        const bool constrained_zero_hessian_enclosed=audit.valid();
        ARIADNE_TEST_ASSERT(constrained_zero_hessian_enclosed);
    }

    // A zero-curvature direction not constrained by A cannot be
    // certified with H+A^T*A; this is not an infeasibility result.
    Matrix<X> missing_direction(1u,2u,zero);
    missing_direction[0u][1u]=X(1u,dp);
    const auto incomplete=validated_kkt_solve(
        semidefinite,missing_direction,f,g);
    ARIADNE_TEST_ASSERT(
        incomplete.status==ValidatedKktStatus::HESSIAN_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!incomplete.primal.has_value());
    ARIADNE_TEST_ASSERT(!incomplete.multipliers.has_value());

    // No constraints: a semidefinite Hessian remains uncertifiable by
    // either direct or augmented SPD tests.
    const Matrix<X> no_constraints(0u,2u,zero);
    const Vector<X> empty_g(0u,zero);
    const auto unconstrained=validated_kkt_solve(
        semidefinite,no_constraints,f,empty_g);
    ARIADNE_TEST_ASSERT(
        unconstrained.status==ValidatedKktStatus::HESSIAN_NOT_POSITIVE_DEFINITE);
    ARIADNE_TEST_ASSERT(!unconstrained.primal.has_value());
    ARIADNE_TEST_ASSERT(!unconstrained.multipliers.has_value());
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
    ARIADNE_TEST_CALL(test_validated_kkt_solve());
    ARIADNE_TEST_CALL(test_validated_kkt_matrix_rhs());
    ARIADNE_TEST_CALL(test_kkt_augmented_hessian());
    return ARIADNE_TEST_FAILURES;
}
