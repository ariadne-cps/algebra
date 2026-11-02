/***************************************************************************
 *            test_symmetric_matrix.cpp
 *
 *  Copyright  2026  Pieter Collins
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

#include <sstream>

// Ensure symmetric_matrix.hpp is self-contained and works as the first include.
#include "algebra/symmetric_matrix.hpp"
#include "numeric/numeric.hpp"
#include "algebra/matrix.hpp"

#include "algebra_testing_context.hpp"

using namespace Ariadne;

namespace {

using X=FloatDPApproximation;
using S=SymmetricMatrix<X>;

Void test_construction_and_access()
{
    S matrix(3u,X(0u,dp));

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.size(), 3u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.row_size(), 3u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.column_size(), 3u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.zero_element(), X(0u,dp));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), matrix.set(0u,2u,X(3u,dp)));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.get(0u,2u), X(3u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.get(2u,0u), X(3u,dp));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), matrix.at(1u,2u)=X(5u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.at(2u,1u), X(5u,dp));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), matrix[0u][1u]=X(7u,dp));
    S const& const_matrix=matrix;
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), const_matrix[1u][0u], X(7u,dp));

    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), *matrix.begin(), matrix.at(0u,0u));

    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)matrix.at(3u,0u));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)matrix.at(0u,3u));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), matrix.resize(4u));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.size(), 4u);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), matrix.set(3u,3u,X(11u,dp)));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.get(3u,3u), X(11u,dp));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), matrix.resize(2u));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), matrix.size(), 2u);
}

Void test_dense_conversion()
{
    Matrix<X> dense_symmetric({
        {X(1u,dp),X(2u,dp)},
        {X(2u,dp),X(4u,dp)}
    });

    S symmetric(dense_symmetric);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), symmetric.at(0u,0u), X(1u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), symmetric.at(0u,1u), X(2u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), symmetric.at(1u,1u), X(4u,dp));

    Matrix<X> dense_again=static_cast<Matrix<X>>(symmetric);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), dense_again[0u][0u], X(1u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), dense_again[0u][1u], X(2u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), dense_again[1u][0u], X(2u,dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), dense_again[1u][1u], X(4u,dp));

    std::ostringstream stream;
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), stream << symmetric);

    Matrix<X> dense_nonsymmetric({
        {X(1u,dp),X(2u,dp)},
        {X(3u,dp),X(4u,dp)}
    });
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)S(dense_nonsymmetric));

    Matrix<X> rectangular(2u,3u,X(0u,dp));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)S(rectangular));
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("symmetric_matrix");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_construction_and_access());
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_dense_conversion());

    

    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
