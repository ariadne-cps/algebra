/***************************************************************************
 *            linear_independence.cpp
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

#include "algebra/linear_independence.hpp"

namespace Ariadne {

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
