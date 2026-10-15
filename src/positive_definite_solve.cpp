/***************************************************************************
 *            positive_definite_solve.cpp
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

#include "algebra/positive_definite_solve.hpp"

#include <utility>

namespace Ariadne {

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

} // namespace Ariadne
