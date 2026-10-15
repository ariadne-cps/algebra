/***************************************************************************
 *            cholesky_decomposition.cpp
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

#include "algebra/cholesky_decomposition.hpp"

namespace Ariadne {

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

} // namespace Ariadne
