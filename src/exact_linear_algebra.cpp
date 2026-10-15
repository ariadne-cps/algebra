/***************************************************************************
 *            exact_linear_algebra.cpp
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

#include "algebra/exact_linear_algebra.hpp"

#include <vector>

namespace Ariadne {

ExactNullSpaceResult exact_null_space(const Matrix<Rational>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    const Rational zero(0);
    Matrix<Rational> R=A;
    std::vector<SizeType> pivot_columns;
    SizeType rank=0;

    for(SizeType col=0; col!=n && rank!=m; ++col) {
        SizeType pivot=rank;
        while(pivot!=m && decide(R[pivot][col]==zero)) { ++pivot; }
        if(pivot==m) { continue; }

        if(pivot!=rank) {
            for(SizeType j=0; j!=n; ++j) {
                const Rational value=R[pivot][j];
                R[pivot][j]=R[rank][j];
                R[rank][j]=value;
            }
        }

        const Rational divisor=R[rank][col];
        for(SizeType j=col; j!=n; ++j) { R[rank][j]/=divisor; }

        for(SizeType i=0; i!=m; ++i) {
            if(i==rank) { continue; }
            const Rational factor=R[i][col];
            if(decide(factor==zero)) { continue; }
            for(SizeType j=col; j!=n; ++j) {
                R[i][j]-=factor*R[rank][j];
            }
        }
        pivot_columns.push_back(col);
        ++rank;
    }

    Matrix<Rational> Z(n,n-rank,zero);
    SizeType free_index=0;
    SizeType pivot_index=0;
    for(SizeType col=0; col!=n; ++col) {
        if(pivot_index!=rank && pivot_columns[pivot_index]==col) {
            ++pivot_index;
            continue;
        }
        Z[col][free_index]=1;
        for(SizeType i=0; i!=rank; ++i) {
            Z[pivot_columns[i]][free_index]=-R[i][col];
        }
        ++free_index;
    }

    return ExactNullSpaceResult{Z,rank};
}

ExactSymmetricDefiniteness exact_symmetric_definiteness(const Matrix<Rational>& A)
{
    ARIADNE_PRECONDITION(A.row_size()==A.column_size());
    const SizeType n=A.row_size();
    const Rational zero(0);

    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=i+1; j!=n; ++j) {
            ARIADNE_PRECONDITION_MSG(decide(A[i][j]==A[j][i]),"Matrix must be symmetric");
        }
    }

    Matrix<Rational> S=A;
    bool singular=false;
    for(SizeType k=0; k!=n; ++k) {
        const Rational pivot=S[k][k];
        if(decide(pivot<zero)) { return ExactSymmetricDefiniteness::INDEFINITE; }

        if(decide(pivot==zero)) {
            for(SizeType j=k+1; j!=n; ++j) {
                if(decide(S[k][j]!=zero)) { return ExactSymmetricDefiniteness::INDEFINITE; }
            }
            singular=true;
            continue;
        }

        for(SizeType i=k+1; i!=n; ++i) {
            for(SizeType j=i; j!=n; ++j) {
                const Rational value=S[i][j]-S[i][k]*S[k][j]/pivot;
                S[i][j]=value;
                if(i!=j) { S[j][i]=value; }
            }
        }
    }

    return singular ? ExactSymmetricDefiniteness::POSITIVE_SEMIDEFINITE
                    : ExactSymmetricDefiniteness::POSITIVE_DEFINITE;
}

} // namespace Ariadne
