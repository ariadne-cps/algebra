/***************************************************************************
 *            tests/test_eigen_decomposition.cpp
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

/** @file tests/test_eigen_decomposition.cpp
 *  @brief Test real approximate and validated eigendecomposition contracts.
 */

#include <cmath>
#include <vector>

#include "algebra/eigen_decomposition.hpp"
#include "algebra_testing_context.hpp"
#include "numeric/numeric.hpp"

using namespace Ariadne;

namespace {

using Status = EigenDecompositionStatus;

// Check inclusion against exact integers using stored directed endpoints.
// A small floating point residual alone is not a validated certificate.
template<class PR>
bool encloses(const FloatBounds<PR>& interval, Int integer, PR pr) {
    const FloatBounds<PR> exact(integer, pr);
    return interval.lower_raw() <= exact.lower_raw()
        && interval.upper_raw() >= exact.upper_raw();
}

template<class PR>
Void test_real_examples(PR pr) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;

    Matrix<AX> a(2u, 2u, AX(0u, pr));
    a[0u][0u] = AX(2u, pr);
    a[0u][1u] = AX(1u, pr);
    a[1u][0u] = AX(1u, pr);
    a[1u][1u] = AX(2u, pr);

    EigenDecomposition<ApproximateTag, PR> approximate(pr);
    const auto approximation = approximate(a);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), approximation.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), approximation.data.has_value());
    if (approximation.data) {
        const auto& values = approximation.data->eigenvalues;
        const auto& vectors = approximation.data->eigenvectors;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(values[0u].get_d() - 1.0) < 1e-12);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(values[1u].get_d() - 3.0) < 1e-12);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(vectors[1u][0u].get_d() + 1.0) < 1e-12);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), std::abs(vectors[1u][1u].get_d() - 1.0) < 1e-12);
        const auto D = approximation.data->diagonal_matrix();
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), D.size() == 2u);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            D.at(0u, 1u).raw() == AX(0u, pr).raw());
        const Matrix<AX> left = a * vectors;
        const Matrix<AX> right = vectors * D;
        for (SizeType row = 0u; row != 2u; ++row) {
            for (SizeType col = 0u; col != 2u; ++col) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs((left[row][col] - right[row][col]).get_d()) < 1e-12);
            }
        }
    }

    Matrix<VX> v(2u, 2u, VX(0u, pr));
    v[0u][0u] = VX(2u, pr);
    v[0u][1u] = VX(1u, pr);
    v[1u][0u] = VX(1u, pr);
    v[1u][1u] = VX(2u, pr);

    EigenDecomposition<ValidatedTag, PR> validated(pr);
    const auto certificate = validated(v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), certificate.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), certificate.data.has_value());
    if (certificate.data) {
        const auto& eigen = *certificate.data;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(eigen.eigenvalues[0u], 1, pr));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(eigen.eigenvalues[1u], 3, pr));
        const auto D = eigen.diagonal_matrix();
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), D.size() == 2u);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            D.at(0u, 1u).lower_raw() == VX(0u, pr).lower_raw() &&
            D.at(0u, 1u).upper_raw() == VX(0u, pr).upper_raw());
        const Matrix<VX> product = v * eigen.eigenvectors;
        const Matrix<VX> expected = eigen.eigenvectors * D;
        for (SizeType row = 0u; row != 2u; ++row) {
            for (SizeType col = 0u; col != 2u; ++col) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    product[row][col].lower_raw() <= expected[row][col].upper_raw()
                    && product[row][col].upper_raw() >= expected[row][col].lower_raw());
            }
        }
    }

    // Nonsymmetric, real distinct spectrum; b!=0 branch.
    v[0u][0u] = VX(2u, pr);
    v[0u][1u] = VX(1u, pr);
    v[1u][0u] = VX(0u, pr);
    v[1u][1u] = VX(3u, pr);
    const auto triangular = validated(v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), triangular.status == Status::SUCCESS);
    if (triangular.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(triangular.data->eigenvalues[0u], 2, pr));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(triangular.data->eigenvalues[1u], 3, pr));
    }

    // b=0, c!=0 branch and diagonal branch (including a zero eigenvalue).
    v[0u][0u] = VX(1u, pr);
    v[0u][1u] = VX(0u, pr);
    v[1u][0u] = VX(1u, pr);
    v[1u][1u] = VX(2u, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), validated(v).status == Status::SUCCESS);
    v[1u][0u] = VX(0u, pr);
    v[0u][0u] = VX(0u, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), validated(v).status == Status::SUCCESS);

    // Avoid losing the small eigenvalue to subtraction cancellation.
    // Both examples have one eigenvalue very close to +1 or -1, while
    // the other has magnitude 1e16.
    Matrix<AX> ill_conditioned(2u, 2u, AX(0u, pr));
    ill_conditioned[0u][0u] = AX(1e16, pr);
    ill_conditioned[0u][1u] = AX(1u, pr);
    ill_conditioned[1u][0u] = AX(1u, pr);
    ill_conditioned[1u][1u] = AX(1u, pr);
    const auto stable_positive = approximate(ill_conditioned);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), stable_positive.status == Status::SUCCESS);
    if (stable_positive.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         std::abs(stable_positive.data->eigenvalues[0u].get_d() - 1.0) < 1e-6);
    }

    ill_conditioned[0u][0u] = AX(-1e16, pr);
    ill_conditioned[1u][1u] = AX(-1, pr);
    const auto stable_negative = approximate(ill_conditioned);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), stable_negative.status == Status::SUCCESS);
    if (stable_negative.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         std::abs(stable_negative.data->eigenvalues[1u].get_d() + 1.0) < 1e-6);
    }

    // Dense symmetric 3x3 Jacobi eigenpairs: exactly 2, 2 +/- sqrt(2).
    // Test the known spectrum independently and check M*V=V*D.
    Matrix<AX> symmetric(3u, 3u, AX(0u, pr));
    symmetric[0u][0u] = AX(2u, pr);
    symmetric[1u][1u] = AX(2u, pr);
    symmetric[2u][2u] = AX(2u, pr);
    symmetric[0u][1u] = AX(1u, pr);
    symmetric[1u][0u] = AX(1u, pr);
    symmetric[1u][2u] = AX(1u, pr);
    symmetric[2u][1u] = AX(1u, pr);
    const auto three = approximate(symmetric);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), three.status == Status::SUCCESS);
    if (three.data) {
        const double expected[3] = {2.0 - std::sqrt(2.0), 2.0, 2.0 + std::sqrt(2.0)};
        for (SizeType col = 0u; col != 3u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             std::abs(three.data->eigenvalues[col].get_d() - expected[col]) < 1e-10);
            for (SizeType row = 0u; row != 3u; ++row) {
                AX action(0u, pr);
                for (SizeType k = 0u; k != 3u; ++k) {
                    action += symmetric[row][k] * three.data->eigenvectors[k][col];
                }
                const AX residual = action - three.data->eigenvalues[col] * three.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                                 std::abs(residual.get_d()) < 1e-10);
            }
        }
    }

    // A dense symmetric 4x4 matrix with an independently constructed
    // orthogonal eigenbasis and known spectrum {1,2,3,4}.
    Matrix<AX> four_by_four(4u, 4u, AX(0u, pr));
    const double entries_4[4][4] = {
        {2.5, -0.5, -1.0, 0.0},
        {-0.5, 2.5, 0.0, -1.0},
        {-1.0, 0.0, 2.5, -0.5},
        {0.0, -1.0, -0.5, 2.5}
    };
    for (SizeType row = 0u; row != 4u; ++row) {
        for (SizeType col = 0u; col != 4u; ++col) {
            four_by_four[row][col] = AX(entries_4[row][col], pr);
        }
    }
    const auto four = approximate(four_by_four);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), four.status == Status::SUCCESS);
    if (four.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             std::abs(four.data->eigenvalues[col].get_d() - static_cast<double>(col + 1u)) < 1e-10);
        }
    }

    // This now enters the general nonsymmetric QR path. Its independent
    // characteristic roots are {0,2,4}, not the former Jacobi spectrum.
    symmetric[0u][1u] = AX(3u, pr);
    EigenDecompositionDiagnostics qr_diagnostics;
    const auto asymmetric = approximate(symmetric, &qr_diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     asymmetric.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        qr_diagnostics.strategy == EigenDecompositionStrategy::SHIFTED_QR);
    if (asymmetric.data) {
        for (SizeType col = 0u; col != 3u; ++col) {
            const double exact = 2.0 * static_cast<double>(col);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(asymmetric.data->eigenvalues[col].get_d() - exact) < 1e-8);
        }
    }
    Matrix<VX> symmetric_bounds(3u, 3u, VX(0u, pr));
    symmetric_bounds[0u][0u] = VX(2u, pr);
    symmetric_bounds[1u][1u] = VX(2u, pr);
    symmetric_bounds[2u][2u] = VX(2u, pr);
    symmetric_bounds[0u][1u] = VX(1u, pr);
    symmetric_bounds[1u][0u] = VX(1u, pr);
    symmetric_bounds[1u][2u] = VX(1u, pr);
    symmetric_bounds[2u][1u] = VX(1u, pr);
    // Gershgorin disks of the original symmetric tridiagonal matrix
    // overlap, but the Jacobi candidate basis permits certification
    // after a verified change of coordinates.
    const auto preconditioned_symmetric = validated(symmetric_bounds);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     preconditioned_symmetric.status == Status::SUCCESS);
    if (preconditioned_symmetric.data) {
        const VX radical = sqrt(VX(2u, pr));
        const VX expected[3] = {VX(2u, pr) - radical,
                                VX(2u, pr),
                                VX(2u, pr) + radical};
        for (SizeType col = 0u; col != 3u; ++col) {
            const VX& enclosed = preconditioned_symmetric.data->eigenvalues[col];
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                enclosed.lower_raw() <= expected[col].lower_raw() &&
                enclosed.upper_raw() >= expected[col].upper_raw());
        }
    }

    // Small nonsymmetric perturbation: the same verification remains
    // rigorous for the actual matrix, not merely for the symmetric seed.
    symmetric_bounds[0u][1u] = VX(9u, pr) / VX(8u, pr);
    const auto preconditioned_nonsymmetric = validated(symmetric_bounds);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     preconditioned_nonsymmetric.status == Status::SUCCESS);

    // Wide uncertain entries are not certified by a narrow midpoint
    // approximation. The method must return INCONCLUSIVE.
    symmetric_bounds[0u][1u] = VX(-2, 2, pr);
    const auto unseparated = validated(symmetric_bounds);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     unseparated.status == Status::INCONCLUSIVE);


    // A symmetric dense interval family with pairwise disjoint Gershgorin
    // disks. Certify real eigenvalue enclosures and normalized eigenvectors
    // of every member, not merely the midpoint matrix.
    Matrix<VX> dense(3u, 3u, VX(0u, pr));
    dense[0u][0u] = VX(1, 2, pr);
    dense[1u][1u] = VX(4u, pr);
    dense[2u][2u] = VX(7u, pr);
    const VX coupling = VX(1u, pr) / VX(16u, pr);
    dense[0u][1u] = coupling;
    dense[1u][0u] = coupling;
    dense[1u][2u] = coupling;
    dense[2u][1u] = coupling;
    const auto dense_certificate = validated(dense);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     dense_certificate.status == Status::SUCCESS);
    if (dense_certificate.data) {
        const auto& eigen = *dense_certificate.data;
        for (SizeType col = 0u; col != 3u; ++col) {
            // Sorted eigenvalue intervals are strictly disjoint.
            if (col > 0u) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    eigen.eigenvalues[col - 1u].upper_raw() <
                    eigen.eigenvalues[col].lower_raw());
            }
        }
        const auto applied = dense * eigen.eigenvectors;
        for (SizeType row = 0u; row != 3u; ++row) {
            for (SizeType col = 0u; col != 3u; ++col) {
                const VX expected = eigen.eigenvalues[col] * eigen.eigenvectors[row][col];
                // The two independently evaluated expressions contain
                // the same exact residual even if their intervals differ.
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    applied[row][col].lower_raw() <= expected.upper_raw() &&
                    applied[row][col].upper_raw() >= expected.lower_raw());
            }
        }
    }

    // The same rigorous method also applies to a nonsymmetric dense
    // matrix with an isolated real spectrum.
    Matrix<VX> nonsymmetric(4u, 4u, VX(0u, pr));
    for (SizeType row = 0u; row != 4u; ++row) {
        nonsymmetric[row][row] = VX(static_cast<unsigned>(2u * row + 1u), pr);
    }
    nonsymmetric[0u][1u] = VX(1u, pr) / VX(8u, pr);
    nonsymmetric[1u][0u] = -coupling;
    nonsymmetric[1u][2u] = coupling;
    nonsymmetric[2u][1u] = coupling;
    nonsymmetric[2u][3u] = VX(1u, pr) / VX(8u, pr);
    nonsymmetric[3u][2u] = -coupling;
    const auto general = validated(nonsymmetric);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     general.status == Status::SUCCESS);
    if (general.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            const VX& value = general.data->eigenvalues[col];
            const VX integer(static_cast<unsigned>(2u * col + 1u), pr);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                value.lower_raw() < integer.lower_raw() &&
                value.upper_raw() > integer.upper_raw());
        }
    }

    // Loss of Gershgorin separation is inconclusive rather than evidence
    // of complex eigenvalues or a defective matrix.
    nonsymmetric[0u][0u] = VX(1u, pr);
    nonsymmetric[1u][1u] = VX(1u, pr);
    const auto overlapping_disks = validated(nonsymmetric);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     overlapping_disks.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     !overlapping_disks.data.has_value());

    // A 1x1 matrix is decomposable even with an interval entry.
    Matrix<VX> scalar(1u, 1u, VX(0u, pr));
    scalar[0u][0u] = VX(1, 2, pr);
    const auto one_by_one = validated(scalar);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), one_by_one.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), one_by_one.data.has_value());
    if (one_by_one.data) {
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            one_by_one.data->eigenvalues[0u].lower_raw() <= VX(1u, pr).lower_raw());
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            one_by_one.data->eigenvalues[0u].upper_raw() >= VX(2u, pr).upper_raw());
    }

    // Real upper triangular matrix of order four. Its eigenvalues are
    // the diagonal entries, independent of the off-diagonal coefficients.
    Matrix<VX> upper(4u, 4u, VX(0u, pr));
    upper[0u][0u] = VX(4u, pr);
    upper[1u][1u] = VX(1u, pr);
    upper[2u][2u] = VX(3u, pr);
    upper[3u][3u] = VX(2u, pr);
    upper[0u][1u] = VX(2u, pr);
    upper[0u][2u] = VX(1u, pr);
    upper[1u][2u] = VX(-1, pr);
    upper[2u][3u] = VX(1u, pr);
    const auto upper_result = validated(upper);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), upper_result.status == Status::SUCCESS);
    if (upper_result.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             encloses(upper_result.data->eigenvalues[col], static_cast<Int>(col + 1u), pr));
        }
        const auto lhs = upper * upper_result.data->eigenvectors;
        for (SizeType row = 0u; row != 4u; ++row) {
            for (SizeType col = 0u; col != 4u; ++col) {
                const VX rhs = upper_result.data->eigenvalues[col] * upper_result.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                                 lhs[row][col].lower_raw() <= rhs.upper_raw() &&
                                                 lhs[row][col].upper_raw() >= rhs.lower_raw());
            }
        }
    }

    // Lower triangular case with the same pairwise distinct spectrum.
    Matrix<VX> lower(4u, 4u, VX(0u, pr));
    for (SizeType row = 0u; row != 4u; ++row) {
        for (SizeType col = 0u; col != 4u; ++col) {
            lower[row][col] = upper[col][row];
        }
    }
    const auto lower_result = validated(lower);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), lower_result.status == Status::SUCCESS);
    if (lower_result.data) {
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                             encloses(lower_result.data->eigenvalues[col], static_cast<Int>(col + 1u), pr));
        }
    }

    // Wide diagonal intervals prevent a consistent eigenvalue ordering.
    Matrix<VX> overlapping(3u, 3u, VX(0u, pr));
    overlapping[0u][0u] = VX(0, 2, pr);
    overlapping[1u][1u] = VX(1, 3, pr);
    overlapping[2u][2u] = VX(4, 5, pr);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     validated(overlapping).status == Status::INCONCLUSIVE);

    // An interval family with strictly separated diagonal values.
    overlapping[0u][0u] = VX(0, 1, pr);
    overlapping[1u][1u] = VX(2, 3, pr);
    const auto separated = validated(overlapping);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), separated.status == Status::SUCCESS);

    // Dimension zero has no individual eigenpairs and is deliberately
    // excluded from this diagonalization API.
    Matrix<VX> empty(0u, 0u, VX(0u, pr));
    const auto empty_result = validated(empty);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), empty_result.status == Status::UNSUPPORTED_DIMENSION);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !empty_result.data.has_value());

    // This reducible symmetric 3x3 matrix has simple real eigenvalues
    // {-1, 0, 1}. Its original Gershgorin disks overlap, but the verified
    // similarity path introduced later must certify the eigenbasis.
    Matrix<VX> decoupled(3u, 3u, VX(0u, pr));
    decoupled[0u][1u] = VX(1u, pr);
    decoupled[1u][0u] = VX(1u, pr);

    const auto decoupled_result = validated(decoupled);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     decoupled_result.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     decoupled_result.data.has_value());
    if (decoupled_result.data) {
        const auto& result = *decoupled_result.data;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         encloses(result.eigenvalues[0u], -1, pr));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         encloses(result.eigenvalues[1u], 0, pr));
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         encloses(result.eigenvalues[2u], 1, pr));
        const auto lhs = decoupled * result.eigenvectors;
        for (SizeType row = 0u; row != 3u; ++row) {
            for (SizeType col = 0u; col != 3u; ++col) {
                const VX rhs = result.eigenvalues[col] * result.eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    lhs[row][col].lower_raw() <= rhs.upper_raw() &&
                    lhs[row][col].upper_raw() >= rhs.lower_raw());
            }
        }
    }

    // Changing only the isolated third diagonal entry to one gives
    // spectrum {-1, 1, 1}. Although the eigenvalue 1 is repeated,
    // the two invariant blocks contribute independent eigenvectors.
    // The block-decomposition algorithm must certify a full real basis.
    decoupled[2u][2u] = VX(1u, pr);
    const auto repeated_spectrum = validated(decoupled);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     repeated_spectrum.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     repeated_spectrum.data.has_value());
    if (repeated_spectrum.data) {
        const auto& result = *repeated_spectrum.data;
        const Int expected[3] = {-1, 1, 1};
        for (SizeType col = 0u; col != 3u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(result.eigenvalues[col], expected[col], pr));
        }
        const auto lhs = decoupled * result.eigenvectors;
        for (SizeType row = 0u; row != 3u; ++row) {
            for (SizeType col = 0u; col != 3u; ++col) {
                const VX rhs = result.eigenvalues[col] * result.eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    lhs[row][col].lower_raw() <= rhs.upper_raw() &&
                    lhs[row][col].upper_raw() >= rhs.lower_raw());
            }
        }
    }

    // A truly defective repeated eigenvalue remains inconclusive even
    // when an independent one-dimensional block is also present.
    Matrix<VX> defective(3u, 3u, VX(0u, pr));
    defective[0u][1u] = VX(1u, pr);
    defective[2u][2u] = VX(1u, pr);
    const auto defective_spectrum = validated(defective);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     defective_spectrum.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     !defective_spectrum.data.has_value());
}

Void test_validated_uncertainty() {
    using X = FloatDPBounds;
    EigenDecomposition<ValidatedTag, DoublePrecision> solver(dp);
    Matrix<X> M(2u, 2u, X(0u, dp));

    // Irrational eigenvalues +/- sqrt(2), checked through exact squares.
    M[0u][1u] = X(1u, dp);
    M[1u][0u] = X(2u, dp);
    auto irrational = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), irrational.status == Status::SUCCESS);
    if (irrational.data) {
        for (SizeType j = 0u; j != 2u; ++j) {
            const X square = sqr(irrational.data->eigenvalues[j]);
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), encloses(square, 2, dp));
        }
    }

    // A nontrivial interval matrix family: every member has two real
    // separated eigenvalues with a certified eigenvector basis.
    M[0u][0u] = X(2, 3, dp);
    M[0u][1u] = X(1u, dp);
    M[1u][0u] = X(1u, dp);
    M[1u][1u] = X(2u, dp);
    const auto family = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), family.status == Status::SUCCESS);

    // A rotation has nonreal eigenvalues; a Jordan block has a double
    // eigenvalue. Neither may claim a real simple eigendecomposition.
    M[0u][0u] = X(0u, dp);
    M[0u][1u] = X(-1, dp);
    M[1u][0u] = X(1u, dp);
    M[1u][1u] = X(0u, dp);
    const auto rotation = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), rotation.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !rotation.data.has_value());

    M[0u][0u] = X(1u, dp);
    M[0u][1u] = X(1u, dp);
    M[1u][0u] = X(0u, dp);
    M[1u][1u] = X(1u, dp);
    const auto jordan = solver(M);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), jordan.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !jordan.data.has_value());
}



/** @brief Exercise repeated real eigenvalues without confusing degeneracy and defect. */
template<class PR>
Void test_repeated_real_spectrum(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    const AX approx_zero(0u, precision);
    const VX validated_zero(0u, precision);
    const VX validated_one(1u, precision);
    EigenDecomposition<ApproximateTag, PR> approximate(precision);
    EigenDecomposition<ValidatedTag, PR> validated(precision);

    // A repeated 2x2 diagonal eigenvalue is not a defective matrix.
    Matrix<AX> diagonal_approx(2u, 2u, approx_zero);
    diagonal_approx[0u][0u] = AX(4u, precision);
    diagonal_approx[1u][1u] = AX(4u, precision);
    const auto two_approx = approximate(diagonal_approx);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     two_approx.status == Status::SUCCESS);
    if (two_approx.data) {
        for (SizeType col = 0u; col != 2u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                two_approx.data->eigenvalues[col].raw() == AX(4u, precision).raw());
        }
    }

    Matrix<VX> diagonal_two(2u, 2u, validated_zero);
    diagonal_two[0u][0u] = VX(4u, precision);
    diagonal_two[1u][1u] = VX(4u, precision);
    const auto two_validated = validated(diagonal_two);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     two_validated.status == Status::SUCCESS);
    if (two_validated.data) {
        for (SizeType col = 0u; col != 2u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(two_validated.data->eigenvalues[col], 4, precision));
        }
    }

    // Pointwise-equal repeated entries are sortable even in ValidatedTag.
    // Another entry can remain a genuine interval if order is unambiguous.
    Matrix<VX> diagonal_four(4u, 4u, validated_zero);
    diagonal_four[0u][0u] = VX(2, 3, precision);
    diagonal_four[1u][1u] = VX(-1, precision);
    diagonal_four[2u][2u] = VX(5u, precision);
    diagonal_four[3u][3u] = VX(-1, precision);
    const auto four = validated(diagonal_four);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     four.status == Status::SUCCESS);
    if (four.data) {
        const auto& data = *four.data;
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            encloses(data.eigenvalues[0u], -1, precision) &&
            encloses(data.eigenvalues[1u], -1, precision) &&
            encloses(data.eigenvalues[2u], 2, precision) &&
            encloses(data.eigenvalues[2u], 3, precision) &&
            encloses(data.eigenvalues[3u], 5, precision));
        const Matrix<VX> product = diagonal_four * data.eigenvectors;
        for (SizeType row = 0u; row != 4u; ++row) {
            for (SizeType col = 0u; col != 4u; ++col) {
                const VX rhs = data.eigenvalues[col] * data.eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    product[row][col].lower_raw() <= rhs.upper_raw() &&
                    product[row][col].upper_raw() >= rhs.lower_raw());
            }
        }
        // The column permutation gives an exact coordinate basis.
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            data.eigenvectors[1u][0u].lower_raw() == validated_one.lower_raw());
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            data.eigenvectors[3u][1u].lower_raw() == validated_one.lower_raw());
    }

    // Identical interval descriptions are not the same real parameter:
    // independent entries may reverse eigenvalue order across the family.
    Matrix<VX> ambiguous(3u, 3u, validated_zero);
    ambiguous[0u][0u] = VX(2, 3, precision);
    ambiguous[1u][1u] = VX(2, 3, precision);
    ambiguous[2u][2u] = VX(5u, precision);
    const auto unorderable = validated(ambiguous);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     unorderable.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     !unorderable.data.has_value());

    // Jacobi's orthogonal eigenvectors remain independent even when
    // the fully connected symmetric matrix 4*I - ones(3,3) has repeated
    // eigenvalues {1,4,4}. This case must NOT be block-decomposable:
    // it preserves independent coverage of the Jacobi path.
    Matrix<AX> symmetric(3u, 3u, approx_zero);
    for (SizeType row = 0u; row != 3u; ++row) {
        for (SizeType col = 0u; col != 3u; ++col) {
            symmetric[row][col] = row == col
                ? AX(3u, precision) : AX(-1, precision);
        }
    }
    EigenDecompositionDiagnostics diagnostics;
    const auto repeated = approximate(symmetric, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     repeated.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy == EigenDecompositionStrategy::JACOBI);
    if (repeated.data) {
        const double expected[3] = {1.0, 4.0, 4.0};
        for (SizeType col = 0u; col != 3u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(repeated.data->eigenvalues[col].get_d() - expected[col]) < 1e-10);
            for (SizeType row = 0u; row != 3u; ++row) {
                AX action(approx_zero);
                for (SizeType k = 0u; k != 3u; ++k) {
                    action += symmetric[row][k] * repeated.data->eigenvectors[k][col];
                }
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs((action - repeated.data->eigenvalues[col] *
                        repeated.data->eigenvectors[row][col]).get_d()) < 1e-10);
            }
        }
    }
}

/** @brief Blockwise diagonalization with noncontiguous invariant components. */
template<class PR>
Void test_block_diagonal_spectrum(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    const AX azero(0u, precision);
    const VX vzero(0u, precision);
    EigenDecomposition<ApproximateTag, PR> approximate(precision);
    EigenDecomposition<ValidatedTag, PR> validated(precision);

    // The two coupled 2x2 blocks occupy indices {0,2} and {1,3},
    // not contiguous intervals. Their spectra {2,4} and {4,6}
    // overlap at 4, but the block eigenvector bases are independent.
    Matrix<AX> a(4u, 4u, azero);
    Matrix<VX> v(4u, 4u, vzero);
    for (SizeType row = 0u; row != 4u; ++row) {
        const unsigned diagonal = (row % 2u == 0u) ? 3u : 5u;
        a[row][row] = AX(diagonal, precision);
        v[row][row] = VX(diagonal, precision);
    }
    a[0u][2u] = AX(1u, precision);
    a[2u][0u] = AX(1u, precision);
    a[1u][3u] = AX(1u, precision);
    a[3u][1u] = AX(1u, precision);
    v[0u][2u] = VX(1u, precision);
    v[2u][0u] = VX(1u, precision);
    v[1u][3u] = VX(1u, precision);
    v[3u][1u] = VX(1u, precision);

    EigenDecompositionDiagnostics diagnostics;
    const auto point = approximate(a, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     point.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy == EigenDecompositionStrategy::BLOCK_DIAGONAL);
    const auto certified = validated(v, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     certified.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy == EigenDecompositionStrategy::BLOCK_DIAGONAL);
    if (point.data && certified.data) {
        // The same diagonal-matrix API must retain the blockwise
        // eigenpair ordering, including the repeated eigenvalue.
        const auto D = certified.data->diagonal_matrix();
        const Matrix<VX> left = v * certified.data->eigenvectors;
        const Matrix<VX> right = certified.data->eigenvectors * D;
        for (SizeType row = 0u; row != 4u; ++row) {
            for (SizeType col = 0u; col != 4u; ++col) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    left[row][col].lower_raw() <= right[row][col].upper_raw() &&
                    left[row][col].upper_raw() >= right[row][col].lower_raw());
            }
        }
        const Int expected[4] = {2, 4, 4, 6};
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                point.data->eigenvalues[col].raw() ==
                    AX(expected[col], precision).raw());
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(certified.data->eigenvalues[col], expected[col], precision));
            for (SizeType row = 0u; row != 4u; ++row) {
                AX action(azero);
                VX rigorous_action(vzero);
                for (SizeType k = 0u; k != 4u; ++k) {
                    action += a[row][k] * point.data->eigenvectors[k][col];
                    rigorous_action += v[row][k] *
                        certified.data->eigenvectors[k][col];
                }
                const AX residual = action -
                    point.data->eigenvalues[col] * point.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs(residual.get_d()) < 1e-11);
                const VX rigorous_rhs = certified.data->eigenvalues[col] *
                    certified.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    rigorous_action.lower_raw() <= rigorous_rhs.upper_raw() &&
                    rigorous_action.upper_raw() >= rigorous_rhs.lower_raw());
            }
        }
    }

    // The second block has spectrum safely above the first, even if its
    // first diagonal coefficient varies independently in a small interval.
    v[1u][1u] = VX(6u, precision) +
                 VX(0, 1, precision) / VX(64u, precision);
    v[3u][3u] = VX(6u, precision);
    const auto family = validated(v);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     family.status == Status::SUCCESS);
    if (family.data) {
        const Int expected[4] = {2, 4, 5, 7};
        for (SizeType col = 0u; col != 4u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(family.data->eigenvalues[col], expected[col], precision));
        }
    }

    // Repeated interval descriptions in distinct blocks do not imply
    // equal real eigenvalues for every possible input matrix.
    v[1u][1u] = VX(5u, precision) +
                 VX(0, 1, precision) / VX(64u, precision);
    v[3u][3u] = VX(5u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(v).status == Status::INCONCLUSIVE);

    // A defective Jordan block is not made diagonalizable by adding
    // another independent diagonalizable block.
    v[0u][0u] = VX(1u, precision);
    v[2u][2u] = VX(1u, precision);
    v[2u][0u] = vzero;
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(v).status == Status::INCONCLUSIVE);
}

/** @brief Repeated triangular eigenvalues that still admit an invertible basis. */
template<class PR>
Void test_repeated_triangular_spectrum(PR precision) {
    using AX = FloatApproximation<PR>;
    using VX = FloatBounds<PR>;
    const AX azero(0u, precision);
    const VX vzero(0u, precision);
    EigenDecomposition<ApproximateTag, PR> approximate(precision);
    EigenDecomposition<ValidatedTag, PR> validated(precision);

    // Non-diagonal upper triangular A, with eigenvalues {1,1,2}.
    // The repeated eigenvalue 1 has two independent eigenvectors:
    // e0 and (0,-1,1); the dependent row equation cancels exactly.
    const Int entries[3][3] = {{1, 1, 1}, {0, 2, 1}, {0, 0, 1}};
    Matrix<AX> upper_approx(3u, 3u, azero);
    Matrix<VX> upper_validated(3u, 3u, vzero);
    for (SizeType row = 0u; row != 3u; ++row) {
        for (SizeType col = 0u; col != 3u; ++col) {
            upper_approx[row][col] = AX(entries[row][col], precision);
            upper_validated[row][col] = VX(entries[row][col], precision);
        }
    }
    const auto aresult = approximate(upper_approx);
    const auto vresult = validated(upper_validated);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     aresult.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     vresult.status == Status::SUCCESS);
    if (aresult.data && vresult.data) {
        const Int expected[3] = {1, 1, 2};
        for (SizeType col = 0u; col != 3u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                aresult.data->eigenvalues[col].raw() ==
                    AX(expected[col], precision).raw());
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                encloses(vresult.data->eigenvalues[col], expected[col], precision));
            for (SizeType row = 0u; row != 3u; ++row) {
                AX lhs(azero);
                for (SizeType k = 0u; k != 3u; ++k) {
                    lhs += upper_approx[row][k] * aresult.data->eigenvectors[k][col];
                }
                const AX residual = lhs -
                    aresult.data->eigenvalues[col] * aresult.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    abs(residual).get_d() < 1e-12);
                VX rigorous_lhs(vzero);
                for (SizeType k = 0u; k != 3u; ++k) {
                    rigorous_lhs += upper_validated[row][k] *
                        vresult.data->eigenvectors[k][col];
                }
                const VX rigorous_rhs = vresult.data->eigenvalues[col] *
                    vresult.data->eigenvectors[row][col];
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    rigorous_lhs.lower_raw() <= rigorous_rhs.upper_raw() &&
                    rigorous_lhs.upper_raw() >= rigorous_rhs.lower_raw());
            }
        }
        // The repeated eigenspace is two-dimensional. Both canonical
        // pivot entries are exact ones in the corresponding columns.
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            vresult.data->eigenvectors[0u][0u].lower_raw() ==
                VX(1u, precision).lower_raw());
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            vresult.data->eigenvectors[2u][1u].lower_raw() ==
                VX(1u, precision).lower_raw());
    }

    // Transposition exercises the distinct lower-triangular recurrence.
    Matrix<VX> lower_validated(3u, 3u, vzero);
    for (SizeType row = 0u; row != 3u; ++row) {
        for (SizeType col = 0u; col != 3u; ++col) {
            lower_validated[row][col] = upper_validated[col][row];
        }
    }
    const auto lower_result = validated(lower_validated);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     lower_result.status == Status::SUCCESS);

    // Different independent interval diagonal entries that happen to use
    // the same bounds do not represent a guaranteed repeated eigenvalue.
    upper_validated[0u][0u] = VX(1, 2, precision);
    upper_validated[2u][2u] = VX(1, 2, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(upper_validated).status == Status::INCONCLUSIVE);

    // A Jordan block is still defective: the repeated-eigenvalue row has
    // a nonzero compatibility obstruction.
    Matrix<VX> jordan(3u, 3u, vzero);
    jordan[0u][0u] = VX(1u, precision);
    jordan[1u][1u] = VX(1u, precision);
    jordan[2u][2u] = VX(2u, precision);
    jordan[0u][1u] = VX(1u, precision);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        validated(jordan).status == Status::INCONCLUSIVE);
}

/** @brief Test dense nonsymmetric approximate simple real eigendecompositions. */
template<class PR>
Void test_nonsymmetric_approximate(PR precision) {
    using X = FloatApproximation<PR>;
    const X zero(0u, precision);
    EigenDecomposition<ApproximateTag, PR> solver(precision);

    // Neither matrix is upper/lower triangular, and both have independent
    // integer spectra by exact characteristic-polynomial factorization.
    const Int three_coefficients[3][3] = {
        {0, 1, 0},
        {-2, 3, 1},
        {3, -3, 6}
    };
    Matrix<X> three(3u, 3u, zero);
    for (SizeType row = 0u; row != 3u; ++row) {
        for (SizeType col = 0u; col != 3u; ++col) {
            three[row][col] = X(three_coefficients[row][col], precision);
        }
    }

    const Int four_coefficients[4][4] = {
        {1, 2, 0, 0},
        {0, 3, 0, 0},
        {-8, 10, 5, 0},
        {-6, 2, 0, 7}
    };
    Matrix<X> four(4u, 4u, zero);
    for (SizeType row = 0u; row != 4u; ++row) {
        for (SizeType col = 0u; col != 4u; ++col) {
            four[row][col] = X(four_coefficients[row][col], precision);
        }
    }

    const auto verify = [&](const Matrix<X>& A,
                            const std::vector<Int>& expected) {
        EigenDecompositionDiagnostics diagnostics;
        const auto answer = solver(A, &diagnostics);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            answer.status == Status::SUCCESS);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            diagnostics.strategy == EigenDecompositionStrategy::SHIFTED_QR);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            answer.data.has_value());
        if (answer.data) {
            const auto& data = *answer.data;
            for (SizeType col = 0u; col != expected.size(); ++col) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    std::abs(data.eigenvalues[col].get_d() -
                             static_cast<double>(expected[col])) < 1e-8);
                for (SizeType row = 0u; row != A.row_size(); ++row) {
                    X action(zero);
                    for (SizeType k = 0u; k != A.row_size(); ++k) {
                        action += A[row][k] * data.eigenvectors[k][col];
                    }
                    const X residual = action -
                        data.eigenvalues[col] * data.eigenvectors[row][col];
                    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                        std::abs(residual.get_d()) < 1e-8);
                }
            }
        }
    };
    verify(three, {1, 3, 5});
    verify(four, {1, 3, 5, 7});

    // The eigenvalues +-i and 2 are not a complete real spectrum.
    Matrix<X> complex(3u, 3u, zero);
    complex[0u][1u] = X(-1, precision);
    complex[1u][0u] = X(1u, precision);
    complex[2u][2u] = X(2u, precision);
    const auto complex_answer = solver(complex);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        complex_answer.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !complex_answer.data.has_value());

    // A repeated eigenvalue does not by itself make a symmetric matrix
    // defective. This matrix is orthogonally diagonalizable with 0,2,2.
    Matrix<X> repeated(3u, 3u, zero);
    repeated[0u][0u] = X(1u, precision);
    repeated[0u][1u] = X(1u, precision);
    repeated[1u][0u] = X(1u, precision);
    repeated[1u][1u] = X(1u, precision);
    repeated[2u][2u] = X(2u, precision);
    const auto repeated_answer = solver(repeated);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        repeated_answer.status == Status::SUCCESS);
    if (repeated_answer.data) {
        const double expected[3] = {0.0, 2.0, 2.0};
        for (SizeType col = 0u; col != 3u; ++col) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                std::abs(repeated_answer.data->eigenvalues[col].get_d() - expected[col]) < 1e-10);
        }
    }

    // A genuine Jordan block has insufficient eigenvectors despite real
    // eigenvalues. It must remain INCONCLUSIVE.
    Matrix<X> defective(3u, 3u, zero);
    defective[0u][0u] = X(2u, precision);
    defective[0u][1u] = X(1u, precision);
    defective[1u][1u] = X(2u, precision);
    defective[2u][2u] = X(3u, precision);
    const auto defective_answer = solver(defective);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        defective_answer.status == Status::INCONCLUSIVE);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        !defective_answer.data.has_value());
}

//! \brief Non-normal dense example with known spectrum that defeats Gershgorin.
template<class PR>
Void test_krawczyk_nonnormal(PR precision) {
    using X = FloatBounds<PR>;
    const X zero(0u, precision);

    // A=S diag(1,3,5,7) S^{-1} with integer unimodular
    // S=[[1,1,0,0],[0,1,0,0],[2,-1,1,0],[1,1,0,1]].
    // The original Gershgorin disks overlap and the symmetric midpoint
    // eigenbasis does not separate the transformed disks. Eigenpair
    // Newton/Krawczyk must certify all four real simple roots.
    const Int coefficients[4][4] = {
        {1, 2, 0, 0},
        {0, 3, 0, 0},
        {-8, 10, 5, 0},
        {-6, 2, 0, 7}
    };
    Matrix<X> A(4u, 4u, zero);
    for (SizeType i = 0u; i != 4u; ++i) {
        for (SizeType j = 0u; j != 4u; ++j) {
            A[i][j] = X(coefficients[i][j], precision);
        }
    }
    EigenDecomposition<ValidatedTag, PR> solver(precision);

    const auto check = [&](const Matrix<X>& matrix) {
        const auto result = solver(matrix);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         result.status == Status::SUCCESS);
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                         result.data.has_value());
        if (result.data) {
            const auto& data = *result.data;
            for (SizeType col = 0u; col != 4u; ++col) {
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                    encloses(data.eigenvalues[col], static_cast<Int>(2u * col + 1u), precision));
                if (col > 0u) {
                    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                        data.eigenvalues[col - 1u].upper_raw() <
                        data.eigenvalues[col].lower_raw());
                }
            }
            const auto action = matrix * data.eigenvectors;
            for (SizeType row = 0u; row != 4u; ++row) {
                for (SizeType col = 0u; col != 4u; ++col) {
                    const X rhs = data.eigenvalues[col] * data.eigenvectors[row][col];
                    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                        action[row][col].lower_raw() <= rhs.upper_raw() &&
                        action[row][col].upper_raw() >= rhs.lower_raw());
                }
            }
        }
    };
    check(A);

    // Uniform certification of a nontrivial interval family must
    // still enclose the unperturbed eigenvalues as one family member.
    A[0u][0u] += X(0, 1, precision) / X(1048576u, precision);
    check(A);
}

//! \brief Optional profiling must not alter eigenpairs or certification.
template<class PR>
Void test_eigen_diagnostics(PR precision) {
    using VX = FloatBounds<PR>;
    using AX = FloatApproximation<PR>;
    const VX zero(0u, precision);
    const VX one(1u, precision);
    Matrix<VX> dense(3u, 3u, zero);
    dense[0u][0u] = one;
    dense[1u][1u] = VX(3u, precision);
    dense[2u][2u] = VX(5u, precision);
    const VX off = one / VX(16u, precision);
    for (SizeType i = 0u; i != 3u; ++i) {
        for (SizeType j = i + 1u; j != 3u; ++j) {
            dense[i][j] = off;
            dense[j][i] = off;
        }
    }
    EigenDecomposition<ValidatedTag, PR> solver(precision);
    const auto baseline = solver(dense);
    EigenDecompositionDiagnostics diagnostics;
    const auto instrumented = solver(dense, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     baseline.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     instrumented.status == baseline.status);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy == EigenDecompositionStrategy::GERSHGORIN);
    if (baseline.data && instrumented.data) {
        for (SizeType i = 0u; i != 3u; ++i) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                baseline.data->eigenvalues[i].lower_raw() ==
                instrumented.data->eigenvalues[i].lower_raw());
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                baseline.data->eigenvalues[i].upper_raw() ==
                instrumented.data->eigenvalues[i].upper_raw());
        }
    }

    Matrix<AX> symmetric(3u, 3u, AX(0u, precision));
    symmetric[0u][0u] = AX(2u, precision);
    symmetric[1u][1u] = AX(2u, precision);
    symmetric[2u][2u] = AX(2u, precision);
    symmetric[0u][1u] = AX(1u, precision);
    symmetric[1u][0u] = AX(1u, precision);
    symmetric[1u][2u] = AX(1u, precision);
    symmetric[2u][1u] = AX(1u, precision);
    EigenDecomposition<ApproximateTag, PR> approximate(precision);
    const auto approx_baseline = approximate(symmetric);
    const auto approx_instrumented = approximate(symmetric, &diagnostics);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     approx_baseline.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     approx_instrumented.status == approx_baseline.status);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
        diagnostics.strategy == EigenDecompositionStrategy::JACOBI);
    if (approx_baseline.data && approx_instrumented.data) {
        for (SizeType i = 0u; i != 3u; ++i) {
            ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                approx_baseline.data->eigenvalues[i].raw() ==
                approx_instrumented.data->eigenvalues[i].raw());
        }
    }
}

//! \brief Regression for an MP-visible but binary64-invisible coupling.
Void test_jacobi_multiprecision() {
    const MultiplePrecision pr(160u);
    using AX = FloatMPApproximation;
    using Raw = typename AX::RawType;
    const AX zero(0u, pr);
    const AX one(1u, pr);
    const AX two_value(2u, pr);
    AX coupling(one);
    for (SizeType k = 0u; k != 75u; ++k) {
        coupling = coupling / two_value;
    }

    Matrix<AX> A(3u, 3u, zero);
    A[0u][0u] = one;
    A[1u][1u] = two_value;
    A[2u][2u] = AX(4u, pr);
    A[0u][1u] = coupling;
    A[1u][0u] = coupling;

    EigenDecomposition<ApproximateTag, MultiplePrecision> solver(pr);
    const auto computed = solver(A);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     computed.status == Status::SUCCESS);
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                     computed.data.has_value());
    if (computed.data) {
        // A binary64-based stopping rule accepts the identity matrix as
        // an eigenbasis; an MP-aware rule must resolve this coupling.
        ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
            abs(computed.data->eigenvectors[1u][0u]).raw() >
            (coupling / two_value).raw());

        const AX tolerance = AX(4096u, pr) * AX(Raw::eps(pr)) *
                             AX(4u, pr) * AX(3u, pr);
        for (SizeType col = 0u; col != 3u; ++col) {
            for (SizeType row = 0u; row != 3u; ++row) {
                AX applied(zero);
                for (SizeType k = 0u; k != 3u; ++k) {
                    applied += A[row][k] * computed.data->eigenvectors[k][col];
                }
                const AX residual = abs(applied -
                    computed.data->eigenvalues[col] * computed.data->eigenvectors[row][col]);
                ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(),
                                                 residual.raw() <= tolerance.raw());
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("eigen_decomposition");
    suite.add("real 1x1 and 2x2 cases", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_real_examples(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_real_examples(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_validated_uncertainty());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_jacobi_multiprecision());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_krawczyk_nonnormal(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_krawczyk_nonnormal(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_eigen_diagnostics(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_eigen_diagnostics(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_nonsymmetric_approximate(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_nonsymmetric_approximate(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_repeated_real_spectrum(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_repeated_real_spectrum(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_repeated_triangular_spectrum(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_repeated_triangular_spectrum(MultiplePrecision(80u)));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_block_diagonal_spectrum(dp));
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_block_diagonal_spectrum(MultiplePrecision(80u)));
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
