/***************************************************************************
 *            tests/test_spectrum.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Algebra.
 *
 *  Ariadne Algebra is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Algebra is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/test_spectrum.cpp
 *  @brief Independently test complex spectra and defective eigenvalues.
 */
#include <cmath>

#include "algebra/spectrum.hpp"
#include "algebra_testing_context.hpp"

using namespace Ariadne;

namespace {

template<class PR>
bool encloses(const FloatBounds<PR>& x, Int value, PR precision) {
    const FloatBounds<PR> exact(value, precision);
    return x.lower_raw() <= exact.lower_raw() &&
           x.upper_raw() >= exact.upper_raw();
}

template<class PR>
Void test_spectrum(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    using Status = SpectralStatus;
    const AX azero(0u, precision);
    const VX vzero(0u, precision);
    Spectrum<ApproximateTag, PR> approximate(precision);
    Spectrum<ValidatedTag, PR> validated(precision);

    // Real matrices with nonreal conjugate eigenvalues are not
    // real-diagonalizable but do have a complete spectrum.
    Matrix<AX> rotation(2u, 2u, azero);
    Matrix<VX> rotation_v(2u, 2u, vzero);
    rotation[0u][1u] = AX(-1, precision);
    rotation[1u][0u] = AX(1u, precision);
    rotation_v[0u][1u] = VX(-1, precision);
    rotation_v[1u][0u] = VX(1u, precision);
    const auto a = approximate(rotation);
    const auto v = validated(rotation_v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        a.status == Status::SUCCESS && a.data.has_value());
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        v.status == Status::SUCCESS && v.data.has_value());
    if (a.data && v.data) {
        const Int signs[2] = {-1, 1};
        for (SizeType k = 0u; k != 2u; ++k) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(a.data->eigenvalues[k].imaginary_part().get_d() -
                    static_cast<double>(signs[k])) < 1e-10);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(v.data->eigenvalues[k].real_part(), 0, precision) &&
                encloses(v.data->eigenvalues[k].imaginary_part(),
                         signs[k], precision));
        }
    }
    // Certified complex conjugate pair across an interval family.
    rotation_v[0u][1u] -= VX(0, 1, precision) / VX(1024u, precision);
    const auto family = validated(rotation_v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        family.status == Status::SUCCESS && family.data.has_value());
    if (family.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            encloses(family.data->eigenvalues[0u].imaginary_part(), -1, precision));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            encloses(family.data->eigenvalues[1u].imaginary_part(), 1, precision));
    }

    // Defective 2x2 Jordan block: two eigenvalues, only one eigenvector.
    Matrix<AX> jordan(2u, 2u, azero);
    Matrix<VX> jordan_v(2u, 2u, vzero);
    jordan[0u][0u] = AX(2u, precision);
    jordan[1u][1u] = AX(2u, precision);
    jordan[0u][1u] = AX(1u, precision);
    jordan_v[0u][0u] = VX(2u, precision);
    jordan_v[1u][1u] = VX(2u, precision);
    jordan_v[0u][1u] = VX(1u, precision);
    const auto ja = approximate(jordan);
    const auto jv = validated(jordan_v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        ja.status == Status::SUCCESS && ja.data.has_value());
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        jv.status == Status::SUCCESS && jv.data.has_value());
    if (ja.data && jv.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            ja.data->eigenvalues.size() == 2u &&
            jv.data->eigenvalues.size() == 2u);
        for (SizeType k = 0u; k != 2u; ++k) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(ja.data->eigenvalues[k].real_part().get_d()-2.0) < 1e-10);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(jv.data->eigenvalues[k].real_part(), 2, precision) &&
                encloses(jv.data->eigenvalues[k].imaginary_part(), 0, precision));
        }
    }

    // Higher-order triangular families also include defective matrices.
    Matrix<AX> triangular(3u, 3u, azero);
    Matrix<VX> triangular_v(3u, 3u, vzero);
    for (SizeType i = 0u; i != 3u; ++i) {
        const unsigned diagonal = i == 2u ? 5u : 2u;
        triangular[i][i] = AX(diagonal, precision);
        triangular_v[i][i] = VX(diagonal, precision);
    }
    triangular[0u][1u] = AX(1u, precision);
    triangular_v[0u][1u] = VX(1u, precision);
    triangular[1u][2u] = AX(1u, precision);
    triangular_v[1u][2u] = VX(1u, precision);
    const auto tri = approximate(triangular);
    const auto tri_v = validated(triangular_v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        tri.status == Status::SUCCESS && tri.data.has_value());
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        tri_v.status == Status::SUCCESS && tri_v.data.has_value());
    if (tri.data && tri_v.data) {
        const Int expected[3] = {2, 2, 5};
        for (SizeType i = 0u; i != 3u; ++i) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(tri_v.data->eigenvalues[i].real_part(),
                         expected[i], precision));
        }
    }
    // Diagonal intervals may cross; coordinate ordering remains valid.
    triangular_v[1u][1u] += VX(0, 1, precision) / VX(64u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(triangular_v).status == Status::SUCCESS);

    // Dense 3x3 complex Schur support has not yet been implemented.
    triangular[2u][0u] = AX(1u, precision);
    triangular_v[2u][0u] = VX(1u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        approximate(triangular).status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(triangular_v).status == Status::INCONCLUSIVE);
}


//! \brief Complex spectra of coupled, interleaved 2x2 structural blocks.
template<class PR>
Void test_structural_complex_spectrum(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    const AX azero(0u, precision);
    const VX vzero(0u, precision);
    Spectrum<ApproximateTag, PR> approximate(precision);
    Spectrum<ValidatedTag, PR> validated(precision);

    // The directed cycle {0,2} represents a real rotation with
    // spectrum {-i,+i}; the scalar components {1},{3} have 3,5.
    // Nonzero entries from the cyclic block into both scalar blocks
    // prevent the matrix from being block diagonal or triangular in
    // its stored order; nevertheless its spectrum is complete.
    const Int entries[4][4] = {
        {0, 1, -1, 0},
        {0, 3, 0, 1},
        {1, 1, 0, 0},
        {0, 0, 0, 5}
    };
    Matrix<AX> a(4u, 4u, azero);
    Matrix<VX> v(4u, 4u, vzero);
    for (SizeType i = 0u; i != 4u; ++i) {
        for (SizeType j = 0u; j != 4u; ++j) {
            a[i][j] = AX(entries[i][j], precision);
            v[i][j] = VX(entries[i][j], precision);
        }
    }
    const auto point = approximate(a);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        point.status == SpectralStatus::SUCCESS && point.data.has_value());
    const auto check = [&](const Matrix<VX>& family) {
        const auto certified = validated(family);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            certified.status == SpectralStatus::SUCCESS &&
            certified.data.has_value());
        if (certified.data) {
            const auto& eigenvalues = certified.data->eigenvalues;
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                eigenvalues.size() == 4u);
            const Int imaginary[2] = {-1, 1};
            for (SizeType k = 0u; k != 2u; ++k) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    encloses(eigenvalues[k].real_part(), 0, precision) &&
                    encloses(eigenvalues[k].imaginary_part(),
                             imaginary[k], precision));
            }
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(eigenvalues[2u].real_part(), 3, precision) &&
                encloses(eigenvalues[2u].imaginary_part(), 0, precision) &&
                encloses(eigenvalues[3u].real_part(), 5, precision));
        }
    };
    check(v);
    if (point.data) {
        const auto& eigenvalues = point.data->eigenvalues;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            eigenvalues.size() == 4u);
        const Int imaginary[2] = {-1, 1};
        for (SizeType k = 0u; k != 2u; ++k) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(eigenvalues[k].real_part().get_d()) < 1e-12 &&
                std::abs(eigenvalues[k].imaginary_part().get_d() -
                    static_cast<double>(imaginary[k])) < 1e-12);
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            std::abs(eigenvalues[2u].real_part().get_d() - 3.0) < 1e-12 &&
            std::abs(eigenvalues[3u].real_part().get_d() - 5.0) < 1e-12);
    }

    // Independent family uncertainty in a cross-block edge and scalar
    // diagonal leaves the conjugate block and full root multiplicity sound.
    v[0u][1u] += VX(0, 1, precision) / VX(1024u, precision);
    v[1u][1u] += VX(0, 1, precision) / VX(1024u, precision);
    check(v);

    // The closing edge makes all 4 vertices mutually reachable; the
    // present exact structural method must refuse to assert a result.
    a[3u][0u] = AX(1u, precision);
    v[3u][0u] = VX(1u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        approximate(a).status == SpectralStatus::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(v).status == SpectralStatus::INCONCLUSIVE);
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("spectrum");
    suite.add("complex, defective, and triangular spectra",
        [](Ariadne::Testing::Context& ctx) {
            AlgebraTesting::CaseScope capture(ctx);
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_spectrum(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_spectrum(MultiplePrecision(80u)));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_structural_complex_spectrum(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_structural_complex_spectrum(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
