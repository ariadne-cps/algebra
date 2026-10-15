/***************************************************************************
 *            algebra/matrix.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

#include "numeric/module.hpp"
#include "numeric/numeric.hpp"

#include <utility>

#include "matrix.hpp"
#include "symmetric_matrix.hpp"

#include "utility/exceptions.hpp"
#include "numeric/floats.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/rational.hpp"
#include "numeric/casts.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"

#include "matrix.tpl.hpp"

namespace Ariadne {

OutputStream& operator<<(OutputStream& os, const PivotMatrix& pv) {
    return os << "PivotMatrix(" << static_cast< Matrix<Int> >(pv) << ")";
}

template<> Matrix<FloatDPBounds> inverse<>(const Matrix<FloatDP>& A) {
    return lu_inverse(Matrix<FloatDPBounds>(A));
}

template<> Matrix<FloatDPBounds> inverse<FloatDPBounds>(const Matrix<FloatDPBounds>& A) {
    return lu_inverse(A);
}

template<> Matrix<FloatDPBounds> solve(const Matrix<FloatDPBounds>& A, const Matrix<FloatDPBounds>& B) {
    return gs_solve(A,B);
}


template class Matrix<RoundedFloatDP>;
template Matrix<RoundedFloatDP> inverse(const Matrix<RoundedFloatDP>&);
template Vector<RoundedFloatDP> solve(const Matrix<RoundedFloatDP>&, const Vector<RoundedFloatDP>&);
template Void normalise_rows(Matrix<RoundedFloatDP>&);


template class Matrix<FloatDPApproximation>;
template Matrix<FloatDPApproximation> inverse(const Matrix<FloatDPApproximation>&);
template Matrix<FloatDPApproximation> solve(const Matrix<FloatDPApproximation>&, const Matrix<FloatDPApproximation>&);
template Vector<FloatDPApproximation> solve(const Matrix<FloatDPApproximation>&, const Vector<FloatDPApproximation>&);
template Tuple<PivotMatrix,Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>> triangular_decomposition(Matrix<FloatDPApproximation> const&);
template Tuple<Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>> orthogonal_decomposition(Matrix<FloatDPApproximation> const&);
template Tuple<Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>> gram_schmidt_orthogonalisation(Matrix<FloatDPApproximation> const&);
template Tuple<Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>,PivotMatrix> orthogonal_decomposition(Matrix<FloatDPApproximation> const&, Bool);
template Vector<FloatDPApproximation> row_norms(Matrix<FloatDPApproximation> const&);
template Matrix<FloatDPApproximation> operator*(PivotMatrix, Matrix<FloatDPApproximation>);

template class Matrix<FloatDPBounds>;
// template Matrix<FloatDPBounds> inverse(const Matrix<FloatDPBounds>&);
template Matrix<FloatDPBounds> lu_inverse(const Matrix<FloatDPBounds>&);
template Matrix<FloatDPBounds> gs_inverse(const Matrix<FloatDPBounds>&);
template Matrix<FloatDPBounds> lu_solve(const Matrix<FloatDPBounds>&, const Matrix<FloatDPBounds>&);
template Vector<FloatDPBounds> solve(const Matrix<FloatDPBounds>&, const Vector<FloatDPBounds>&);
template Vector<FloatDPBounds> lu_solve(const Matrix<FloatDPBounds>&, const Vector<FloatDPBounds>&);
template Vector<FloatDPBounds> gs_solve(const Matrix<FloatDPBounds>&, const Vector<FloatDPBounds>&);
template Matrix<FloatDPBounds> gs_solve(const Matrix<FloatDPBounds>&, const Matrix<FloatDPBounds>&);
template Matrix<MidpointType<FloatDPBounds>> midpoint(Matrix<FloatDPBounds> const&);
template Tuple<PivotMatrix,Matrix<FloatDPBounds>,Matrix<FloatDPBounds>> triangular_decomposition(Matrix<FloatDPBounds> const&);
template Tuple<Matrix<FloatDPBounds>,Matrix<FloatDPBounds>> orthogonal_decomposition(Matrix<FloatDPBounds> const&);
template Tuple<Matrix<FloatDPBounds>,Matrix<FloatDPBounds>> gram_schmidt_orthogonalisation(Matrix<FloatDPBounds> const&);
template Matrix<FloatDPBounds> operator*(PivotMatrix, Matrix<FloatDPBounds>);

template class Matrix<FloatDP>;
template Matrix<ExactType<FloatDPApproximation>> cast_exact(const Matrix<FloatDPApproximation>&);


template class Matrix<FloatMPApproximation>;
template Matrix<FloatMPApproximation> inverse(const Matrix<FloatMPApproximation>&);
template Matrix<FloatMPApproximation> solve(const Matrix<FloatMPApproximation>&, const Matrix<FloatMPApproximation>&);
template Vector<FloatMPApproximation> solve(const Matrix<FloatMPApproximation>&, const Vector<FloatMPApproximation>&);
template Tuple<PivotMatrix,Matrix<FloatMPApproximation>,Matrix<FloatMPApproximation>> triangular_decomposition(Matrix<FloatMPApproximation> const&);
template Tuple<Matrix<FloatMPApproximation>,Matrix<FloatMPApproximation>> orthogonal_decomposition(Matrix<FloatMPApproximation> const&);
template Tuple<Matrix<FloatMPApproximation>,Matrix<FloatMPApproximation>> gram_schmidt_orthogonalisation(Matrix<FloatMPApproximation> const&);
template Tuple<Matrix<FloatMPApproximation>,Matrix<FloatMPApproximation>,PivotMatrix> orthogonal_decomposition(Matrix<FloatMPApproximation> const&, Bool);
template Vector<FloatMPApproximation> row_norms(Matrix<FloatMPApproximation> const&);
template Matrix<FloatMPApproximation> operator*(PivotMatrix, Matrix<FloatMPApproximation>);

template class Matrix<FloatMPBounds>;
template Matrix<FloatMPBounds> inverse(const Matrix<FloatMPBounds>&);
template Matrix<FloatMPBounds> lu_inverse(const Matrix<FloatMPBounds>&);
template Matrix<FloatMPBounds> gs_inverse(const Matrix<FloatMPBounds>&);
template Matrix<FloatMPBounds> lu_solve(const Matrix<FloatMPBounds>&, const Matrix<FloatMPBounds>&);
template Vector<FloatMPBounds> solve(const Matrix<FloatMPBounds>&, const Vector<FloatMPBounds>&);
template Vector<FloatMPBounds> lu_solve(const Matrix<FloatMPBounds>&, const Vector<FloatMPBounds>&);
template Vector<FloatMPBounds> gs_solve(const Matrix<FloatMPBounds>&, const Vector<FloatMPBounds>&);
template Matrix<FloatMPBounds> gs_solve(const Matrix<FloatMPBounds>&, const Matrix<FloatMPBounds>&);
template Matrix<MidpointType<FloatMPBounds>> midpoint(Matrix<FloatMPBounds> const&);
template Tuple<PivotMatrix,Matrix<FloatMPBounds>,Matrix<FloatMPBounds>> triangular_decomposition(Matrix<FloatMPBounds> const&);
template Tuple<Matrix<FloatMPBounds>,Matrix<FloatMPBounds>> orthogonal_decomposition(Matrix<FloatMPBounds> const&);
template Tuple<Matrix<FloatMPBounds>,Matrix<FloatMPBounds>> gram_schmidt_orthogonalisation(Matrix<FloatMPBounds> const&);
template Matrix<FloatMPBounds> operator*(PivotMatrix, Matrix<FloatMPBounds>);

template class Matrix<FloatMP>;
template Matrix<ExactType<FloatMPApproximation>> cast_exact(const Matrix<FloatMPApproximation>&);


template class Matrix<Real>;

template PositiveFloatDPUpperBound sup_norm(const Matrix<FloatDPBounds>& A);
template FloatDPUpperBound log_norm(const Matrix<FloatDPBounds>& A);

template class Matrix<Dyadic>;

template class Matrix<Rational>;
template Matrix<Rational> inverse(const Matrix<Rational>&);
template Matrix<Rational> solve(const Matrix<Rational>&, const Matrix<Rational>&);
template Vector<Rational> solve(const Matrix<Rational>&, const Vector<Rational>&);
template<> Matrix<Rational> midpoint(Matrix<Rational> const& A) { return A; }

template class SymmetricMatrix<FloatDPApproximation>;

// Validated matrix operations: positive definiteness, Cholesky,
// SPD solves and interval linear independence.
PositiveDefinitenessStatus validated_positive_definiteness(
    const SymmetricMatrix<FloatDPBounds>& A)
{
    const SizeType n=A.size();
    Matrix<FloatDPBounds> S=static_cast<Matrix<FloatDPBounds>>(A);

    for(SizeType k=0; k!=n; ++k) {
        const FloatDPBounds pivot=S[k][k];

        // A positive Schur pivot permits a validated LDL^T elimination.
        if(definitely(pivot.lower()>0)) {
            for(SizeType i=k+1; i!=n; ++i) {
                for(SizeType j=i; j!=n; ++j) {
                    const FloatDPBounds entry=S[i][j]-S[i][k]*S[k][j]/pivot;
                    S[i][j]=entry;
                    S[j][i]=entry;
                }
            }
            continue;
        }

        // With all earlier pivots positive, a nonpositive Schur pivot
        // excludes positive definiteness for every represented matrix.
        if(definitely(pivot.upper()<=0)) {
            return PositiveDefinitenessStatus::NOT_POSITIVE_DEFINITE;
        }

        return PositiveDefinitenessStatus::INDETERMINATE;
    }

    return PositiveDefinitenessStatus::POSITIVE_DEFINITE;
}


ValidatedCholeskyResult validated_cholesky_decomposition(
    const SymmetricMatrix<FloatDPBounds>& A)
{
    const SizeType n=A.size();
    Matrix<FloatDPBounds> L(n,n,A.zero_element());

    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=0; j!=i; ++j) {
            FloatDPBounds entry=A[i][j];
            for(SizeType k=0; k!=j; ++k) {
                entry-=L[i][k]*L[j][k];
            }
            // The earlier positive diagonal pivot is separated from zero.
            L[i][j]=entry/L[j][j];
        }

        FloatDPBounds pivot=A[i][i];
        for(SizeType k=0; k!=i; ++k) {
            // sqr handles intervals that straddle zero without inventing
            // a negative lower bound for the square.
            pivot-=sqr(L[i][k]);
        }

        if(definitely(pivot.lower()>0)) {
            // Directed double sqrt_rnd may cycle between adjacent values
            // near a representable square. MPFR square roots terminate
            // and round correctly; conversion back to FloatDPBounds is
            // outward rounded, preserving the enclosure guarantee.
            const FloatMPBounds pivot_mp(pivot,MultiplePrecision(64));
            L[i][i]=FloatDPBounds(sqrt(pivot_mp),dp);
            continue;
        }

        if(definitely(pivot.upper()<=0)) {
            return {PositiveDefinitenessStatus::NOT_POSITIVE_DEFINITE,std::nullopt};
        }

        return {PositiveDefinitenessStatus::INDETERMINATE,std::nullopt};
    }

    return {PositiveDefinitenessStatus::POSITIVE_DEFINITE,std::move(L)};
}


ValidatedPositiveDefiniteSolveResult<Matrix<FloatDPBounds>>
validated_positive_definite_solve(
    const SymmetricMatrix<FloatDPBounds>& A,
    const Matrix<FloatDPBounds>& B)
{
    ARIADNE_PRECONDITION(B.row_size()==A.size());

    const auto factorization=validated_cholesky_decomposition(A);
    if(!factorization.lower.has_value()) {
        return {factorization.status,std::nullopt};
    }

    const Matrix<FloatDPBounds>& L=*factorization.lower;
    const SizeType n=A.size();
    const SizeType rhs_count=B.column_size();
    Matrix<FloatDPBounds> X=B;

    // Forward substitution: LY=B. The solution matrix stores Y in place.
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType col=0; col!=rhs_count; ++col) {
            FloatDPBounds entry=X[i][col];
            for(SizeType k=0; k!=i; ++k) {
                entry-=L[i][k]*X[k][col];
            }
            X[i][col]=entry/L[i][i];
        }
    }

    // Back substitution: L^T X=Y, overwriting Y from bottom to top.
    for(SizeType i=n; i!=0;) {
        --i;
        for(SizeType col=0; col!=rhs_count; ++col) {
            FloatDPBounds entry=X[i][col];
            for(SizeType k=i+1; k!=n; ++k) {
                entry-=L[k][i]*X[k][col];
            }
            X[i][col]=entry/L[i][i];
        }
    }

    return {PositiveDefinitenessStatus::POSITIVE_DEFINITE,std::move(X)};
}

ValidatedPositiveDefiniteSolveResult<Vector<FloatDPBounds>>
validated_positive_definite_solve(
    const SymmetricMatrix<FloatDPBounds>& A,
    const Vector<FloatDPBounds>& b)
{
    ARIADNE_PRECONDITION(b.size()==A.size());

    Matrix<FloatDPBounds> B(b.size(),1u,b.zero_element());
    for(SizeType i=0; i!=b.size(); ++i) {
        B[i][0u]=b[i];
    }

    const auto matrix_solution=validated_positive_definite_solve(A,B);
    if(!matrix_solution.solution.has_value()) {
        return {matrix_solution.status,std::nullopt};
    }

    Vector<FloatDPBounds> x(b.size(),b.zero_element());
    for(SizeType i=0; i!=b.size(); ++i) {
        x[i]=(*matrix_solution.solution)[i][0u];
    }
    return {matrix_solution.status,std::move(x)};
}


namespace {

LinearIndependenceStatus classify_linear_independence(
    const Matrix<FloatDPBounds>& A, bool rows)
{
    const SizeType count=rows ? A.row_size() : A.column_size();
    const SizeType length=rows ? A.column_size() : A.row_size();

    // More vectors than coordinates always implies dependence.
    if(count>length) {
        return LinearIndependenceStatus::DEPENDENT;
    }

    SymmetricMatrix<FloatDPBounds> gram(count,A.zero_element());
    for(SizeType i=0; i!=count; ++i) {
        for(SizeType j=i; j!=count; ++j) {
            FloatDPBounds entry=A.zero_element();
            for(SizeType k=0; k!=length; ++k) {
                const FloatDPBounds& ai=A.at(rows ? i : k,rows ? k : i);
                const FloatDPBounds& aj=A.at(rows ? j : k,rows ? k : j);
                entry+=(i==j) ? sqr(ai) : ai*aj;
            }
            gram[i][j]=entry;
        }
    }

    const auto definiteness=validated_positive_definiteness(gram);
    if(definiteness==PositiveDefinitenessStatus::POSITIVE_DEFINITE) {
        return LinearIndependenceStatus::INDEPENDENT;
    }
    if(definiteness==PositiveDefinitenessStatus::NOT_POSITIVE_DEFINITE) {
        // Since every exact Gram matrix is positive semidefinite,
        // this result proves that every represented matrix is rank-deficient.
        return LinearIndependenceStatus::DEPENDENT;
    }
    return LinearIndependenceStatus::INDETERMINATE;
}

} // namespace

LinearIndependenceStatus validated_independent_rows(
    const Matrix<FloatDPBounds>& A)
{
    return classify_linear_independence(A,true);
}

LinearIndependenceStatus validated_independent_columns(
    const Matrix<FloatDPBounds>& A)
{
    return classify_linear_independence(A,false);
}


// Validated strictly convex equality-constrained KKT solve.
// All computed matrices enclose the exact expressions for each common
// selection of the input intervals; correlations need not be retained.
ValidatedKktResult validated_kkt_solve(
    const SymmetricMatrix<FloatDPBounds>& H,
    const Matrix<FloatDPBounds>& A,
    const Vector<FloatDPBounds>& f,
    const Vector<FloatDPBounds>& g)
{
    const SizeType n=H.size();
    const SizeType m=A.row_size();
    ARIADNE_PRECONDITION(A.column_size()==n);
    ARIADNE_PRECONDITION(f.size()==n);
    ARIADNE_PRECONDITION(g.size()==m);

    // A single validated solve encloses H^{-1}*f and all columns of
    // H^{-1}*A^T without forming an interval inverse.
    Matrix<FloatDPBounds> rhs(n,m+1,H.zero_element());
    for(SizeType i=0; i!=n; ++i) {
        rhs[i][0u]=f[i];
        for(SizeType j=0; j!=m; ++j) {
            rhs[i][j+1]=A[j][i];
        }
    }

    const auto h_solve=validated_positive_definite_solve(H,rhs);
    if(!h_solve.solution.has_value()) {
        if(h_solve.status==PositiveDefinitenessStatus::NOT_POSITIVE_DEFINITE) {
            return {ValidatedKktStatus::HESSIAN_NOT_POSITIVE_DEFINITE,std::nullopt,std::nullopt};
        }
        return {ValidatedKktStatus::INDETERMINATE,std::nullopt,std::nullopt};
    }

    // The Schur complement itself can certify full row rank; an
    // independent Gram test is unnecessary on the success path.
    const Matrix<FloatDPBounds>& Y=*h_solve.solution;
    SymmetricMatrix<FloatDPBounds> schur(m,H.zero_element());
    Vector<FloatDPBounds> dual_rhs(m,H.zero_element());

    // Exact S is symmetric; fill each independent interval coefficient
    // once so the representation preserves that structural symmetry.
    for(SizeType i=0; i!=m; ++i) {
        FloatDPBounds bi=-g[i];
        for(SizeType k=0; k!=n; ++k) {
            bi+=A[i][k]*Y[k][0u];
        }
        dual_rhs[i]=bi;

        for(SizeType j=i; j!=m; ++j) {
            FloatDPBounds sij=H.zero_element();
            for(SizeType k=0; k!=n; ++k) {
                sij+=A[i][k]*Y[k][j+1];
            }
            schur[i][j]=sij;
        }
    }

    const auto schur_solve=validated_positive_definite_solve(schur,dual_rhs);
    if(!schur_solve.solution.has_value()) {
        // An inconclusive interval Schur complement is not a proof of
        // rank deficiency: cancellations in A*H^{-1}*A^T are correlated.
        // Only a separate universal rank certificate can prove
        // dependent equality constraints.
        if(validated_independent_rows(A)==LinearIndependenceStatus::DEPENDENT) {
            return {ValidatedKktStatus::DEPENDENT_CONSTRAINTS,std::nullopt,std::nullopt};
        }
        return {ValidatedKktStatus::INDETERMINATE,std::nullopt,std::nullopt};
    }

    const Vector<FloatDPBounds>& multipliers=*schur_solve.solution;
    Vector<FloatDPBounds> primal(n,H.zero_element());
    for(SizeType i=0; i!=n; ++i) {
        FloatDPBounds xi=Y[i][0u];
        for(SizeType j=0; j!=m; ++j) {
            xi-=Y[i][j+1]*multipliers[j];
        }
        primal[i]=xi;
    }

    return {ValidatedKktStatus::SOLVED,std::move(primal),multipliers};
}

} // namespace Ariadne

#include "interval/interval.hpp"

namespace Ariadne {
template class Matrix<FloatDPUpperInterval>;
template Matrix<SingletonType<FloatDPUpperInterval>> cast_singleton(Matrix<FloatDPUpperInterval> const&);
template Matrix<MidpointType<FloatDPUpperInterval>> midpoint(Matrix<FloatDPUpperInterval> const&);
template Matrix<FloatDPUpperInterval> inverse(const Matrix<FloatDPUpperInterval>&);
template Vector<FloatDPUpperInterval> solve(const Matrix<FloatDPUpperInterval>&, const Vector<FloatDPUpperInterval>&);

template class Matrix<FloatMPUpperInterval>;
} // namespace Ariadne

