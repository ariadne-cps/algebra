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

    // Closing the 3x3 cycle produces a dense, non-triangular real
    // matrix with characteristic polynomial
    // lambda^3 - 9*lambda^2 + 24*lambda - 21.
    // Its independently computed roots are one real and one
    // nonreal conjugate pair. The new approximate Schur backend
    // must report the complete spectrum, without eigenvectors.
    triangular[2u][0u] = AX(1u, precision);
    triangular_v[2u][0u] = VX(1u, precision);
    const auto dense = approximate(triangular);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        dense.status == Status::SUCCESS && dense.data.has_value());
    if (dense.data) {
        const double expected_real = 5.10380340273554;
        const double expected_pair_real = 1.94809829863223;
        const double expected_pair_imag = 0.565235851677171;
        SizeType real_count = 0u;
        SizeType negative_count = 0u;
        SizeType positive_count = 0u;
        for (const auto& eigenvalue: dense.data->eigenvalues) {
            const double real_part = eigenvalue.real_part().get_d();
            const double imaginary_part = eigenvalue.imaginary_part().get_d();
            if (std::abs(real_part - expected_real) < 1e-8 &&
                std::abs(imaginary_part) < 1e-8) {
                ++real_count;
            } else if (std::abs(real_part - expected_pair_real) < 1e-8 &&
                       std::abs(imaginary_part + expected_pair_imag) < 1e-8) {
                ++negative_count;
            } else if (std::abs(real_part - expected_pair_real) < 1e-8 &&
                       std::abs(imaginary_part - expected_pair_imag) < 1e-8) {
                ++positive_count;
            }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            dense.data->eigenvalues.size() == 3u &&
            real_count == 1u && negative_count == 1u && positive_count == 1u);
    }
    // The interval family has no general complex certification yet.
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

    // A closing edge makes all four vertices strongly connected.
    // It is no longer a structural 1x1/2x2 case, but the new
    // approximate Schur backend must still compute the spectrum.
    // Characteristic polynomial:
    // lambda^4-8*lambda^3+16*lambda^2-9*lambda+16.
    a[3u][0u] = AX(1u, precision);
    v[3u][0u] = VX(1u, precision);
    const auto dense = approximate(a);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        dense.status == SpectralStatus::SUCCESS && dense.data.has_value());
    if (dense.data) {
        const double expected_real[2] = {2.90371299040682, 5.07346898592983};
        const double pair_real = 0.0114090118316761;
        const double pair_imag = 1.04208853101622;
        SizeType real_count = 0u, positive_count = 0u, negative_count = 0u;
        for (const auto& eigenvalue: dense.data->eigenvalues) {
            const double re = eigenvalue.real_part().get_d();
            const double im = eigenvalue.imaginary_part().get_d();
            if (std::abs(im) < 1e-8 &&
                (std::abs(re - expected_real[0u]) < 1e-8 ||
                 std::abs(re - expected_real[1u]) < 1e-8)) {
                ++real_count;
            } else if (std::abs(re - pair_real) < 1e-8 &&
                       std::abs(im - pair_imag) < 1e-8) {
                ++positive_count;
            } else if (std::abs(re - pair_real) < 1e-8 &&
                       std::abs(im + pair_imag) < 1e-8) {
                ++negative_count;
            }
        }
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            dense.data->eigenvalues.size() == 4u &&
            real_count == 2u && negative_count == 1u && positive_count == 1u);
    }
    // The full interval-family complex spectrum is not certified.
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(v).status == SpectralStatus::INCONCLUSIVE);
}

//! \brief Mixed SCCs: an irreducible 3x3, a 2x2, and a scalar.
template<class PR>
Void test_mixed_schur_components(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    Matrix<AX> a(6u, 6u, AX(0u, precision));
    Matrix<VX> v(6u, 6u, VX(0u, precision));

    // Interleave coordinates so the matrix is not triangular in its
    // stored order. SCC {0,2,5} has characteristic polynomial
    // (lambda-9)(lambda^2+9); SCC {1,4} is a rotation with
    // eigenvalues {-i,+i}; SCC {3} has eigenvalue 7.
    // Directed cross-block edges {0,2,5}->{1,4}->{3} do not
    // alter the determinant of lambda*I-A.
    const Int coefficients[6][6] = {
        { 4,  2,  5,  0,  0,  0},
        { 0,  0,  0,  3, -1,  0},
        { 3,  0,  4, -2,  0, -4},
        { 0,  0,  0,  7,  0,  0},
        { 0,  1,  0,  0,  0,  0},
        {-4,  0,  0,  0, -1,  1}
    };
    for (SizeType i = 0u; i != 6u; ++i) {
        for (SizeType j = 0u; j != 6u; ++j) {
            a[i][j] = AX(coefficients[i][j], precision);
            v[i][j] = VX(coefficients[i][j], precision);
        }
    }

    const auto approximate = Spectrum<ApproximateTag, PR>(precision)(a);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        approximate.status == SpectralStatus::SUCCESS &&
        approximate.data.has_value());
    if (approximate.data) {
        const auto& eigenvalues = approximate.data->eigenvalues;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            eigenvalues.size() == 6u);
        const double expected_real[2] = {7.0, 9.0};
        const double expected_imag[4] = {-3.0, -1.0, 1.0, 3.0};
        SizeType real_hits[2] = {0u, 0u};
        SizeType imaginary_hits[4] = {0u, 0u, 0u, 0u};
        for (const auto& eigenvalue: eigenvalues) {
            const double re = eigenvalue.real_part().get_d();
            const double im = eigenvalue.imaginary_part().get_d();
            for (SizeType k = 0u; k != 2u; ++k) {
                if (std::abs(re - expected_real[k]) < 1e-8 &&
                    std::abs(im) < 1e-8) {
                    ++real_hits[k];
                }
            }
            for (SizeType k = 0u; k != 4u; ++k) {
                if (std::abs(re) < 1e-8 &&
                    std::abs(im - expected_imag[k]) < 1e-8) {
                    ++imaginary_hits[k];
                }
            }
        }
        for (SizeType k = 0u; k != 2u; ++k) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                real_hits[k] == 1u);
        }
        for (SizeType k = 0u; k != 4u; ++k) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                imaginary_hits[k] == 1u);
        }
    }

    // A numerical Schur factor for the 3x3 block is not a rigorous
    // enclosure for the corresponding interval matrix family.
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        Spectrum<ValidatedTag, PR>(precision)(v).status ==
            SpectralStatus::INCONCLUSIVE);

    // An uncertain cross-block coefficient may include zero but
    // must never create a reverse directed edge in SCC analysis.
    v[0u][1u] += VX(-1, 1, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        Spectrum<ValidatedTag, PR>(precision)(v).status ==
            SpectralStatus::INCONCLUSIVE);
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
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_mixed_schur_components(dp));
            ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(),
                test_mixed_schur_components(MultiplePrecision(80u)));
        });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
