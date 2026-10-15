/***************************************************************************
 *            positive_definiteness.cpp
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

#include "algebra/positive_definiteness.hpp"

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

} // namespace Ariadne
