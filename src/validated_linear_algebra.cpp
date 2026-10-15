/***************************************************************************
 *            validated_linear_algebra.cpp
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

#include "algebra/symmetric_matrix.hpp"

#include <utility>

namespace Ariadne {

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
            L[i][i]=sqrt(pivot);
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

} // namespace Ariadne
