/***************************************************************************
 *            algebra/schur_decomposition.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/schur_decomposition.hpp
 *  @brief Real Schur reduction using Householder reflections and shifted QR.
 *
 *  @details Mathematical background and original references are documented
 *  on SchurDecomposition. Complex-conjugate shifts now use Francis's
 *  implicit real double-shift bulge chase.
 */
#ifndef ARIADNE_SCHUR_DECOMPOSITION_HPP
#define ARIADNE_SCHUR_DECOMPOSITION_HPP

#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "algebra/matrix.hpp"

namespace Ariadne {

//! \brief Status of a Schur computation; INCONCLUSIVE is not nonexistence.
enum class SchurStatus { SUCCESS, INCONCLUSIVE, UNSUPPORTED_DIMENSION };

//! \brief A leading real invariant subspace from Schur factors.
//! \details Its n-by-k orthonormal basis U and k-by-k restriction B
//! satisfy A*U approximately equal to U*B. This does not require
//! eigenvectors or diagonalizability, and is NOT a validated certificate.
//! The basis spans only a leading Schur block prefix, not an arbitrary
//! selected eigenvalue cluster.
//! \see SchurDecompositionData::leading_invariant_subspace
template<class P, class PR>
struct SchurInvariantSubspaceData {
    using X = FloatType<P, PR>;
    Matrix<X> basis;       //!< n-by-k orthonormal columns.
    Matrix<X> restriction; //!< k-by-k representation of A on that subspace.
};

//! \brief Real Schur factors approximating A*Q = Q*T.
//! \details Q is approximately orthogonal; T is real upper quasi-triangular.
//! Its diagonal 1x1 blocks hold real eigenvalues, while 2x2 blocks hold
//! complex-conjugate eigenvalue pairs (negative 2x2 discriminant).
//! A real Schur form exists even if A is defective and has no eigenbasis.
//! Complete leading Schur blocks generate real invariant subspaces even
//! for repeated eigenvalues and nonreal complex-conjugate pairs.
//! The factors are floating-point approximations, NOT interval certificates.
//! \see SchurDecomposition for the underlying reductions and references.
template<class P, class PR>
struct SchurDecompositionData {
    using X = FloatType<P, PR>;
    Matrix<X> orthogonal;
    Matrix<X> quasi_triangular;
    std::vector<SizeType> block_sizes;

    //! \brief Extract a leading invariant subspace of the real Schur form.
    //! \param block_count Number of complete, initial diagonal Schur
    //! blocks to include (between one and block_sizes.size()).
    //! \param precision Numeric precision used for the new real matrices.
    //! \return The orthonormal basis U and its restriction B, or no
    //! result if the requested prefix is invalid or cuts a 2x2 block.
    //!
    //! \details Suppose Q=[Q1,Q2] and the real Schur matrix is
    //! \f[ T=\begin{pmatrix}B&C\\0&D\end{pmatrix}. \f]
    //! Since A*Q approximately equals Q*T, its first k columns give
    //! \f[ A Q_1 \approx Q_1 B,\quad Q_1^{\mathsf T}Q_1\approx I. \f]
    //! Hence span(Q1) is an approximate A-invariant REAL subspace.
    //! Unlike an eigenbasis, the restriction B need not be diagonal:
    //! this construction also works for a defective Jordan block
    //! and for real 2x2 conjugate-pair blocks.
    //!
    //! A real Schur 2x2 block must never be split. Only prefixes of
    //! the CURRENT Schur ordering can be extracted; selecting an
    //! arbitrary cluster requires an additional Schur-reordering
    //! algorithm that is not yet implemented. The method checks the
    //! complete block partition and that the lower-left T coupling
    //! at the requested boundary is identically zero. It does not
    //! rerun QR or allocate complex matrices.
    //!
    //! This method is available only for ApproximateTag. The
    //! ValidatedTag Schur backend is not implemented and cannot
    //! certify invariant subspaces merely from floating factors.
    //! See J. H. Wilkinson, The Algebraic Eigenvalue Problem (1965),
    //! Chapters 6 and 8, and SchurDecomposition for the underlying
    //! Schur/Householder/Francis derivations.
    std::optional<SchurInvariantSubspaceData<P, PR>>
    leading_invariant_subspace(SizeType block_count, PR precision) const {
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return std::nullopt;
        } else {
            const SizeType n = orthogonal.row_size();
            if (n == 0u || orthogonal.column_size() != n ||
                quasi_triangular.row_size() != n ||
                quasi_triangular.column_size() != n ||
                block_count == 0u || block_count > block_sizes.size()) {
                return std::nullopt;
            }
            SizeType covered = 0u;
            SizeType dimension = 0u;
            for (SizeType i = 0u; i != block_sizes.size(); ++i) {
                const SizeType width = block_sizes[i];
                if ((width != 1u && width != 2u) ||
                    width > n - covered) {
                    return std::nullopt;
                }
                covered += width;
                if (i + 1u == block_count) {
                    dimension = covered;
                }
            }
            if (covered != n || dimension == 0u) {
                return std::nullopt;
            }
            const X zero(0u, precision);
            // A block Schur prefix has no coupling from its
            // complement into the leading invariant subspace.
            // The Schur solver explicitly zeros deflated entries.
            for (SizeType row = dimension; row != n; ++row) {
                for (SizeType col = 0u; col != dimension; ++col) {
                    const X& entry = quasi_triangular[row][col];
                    if (!is_finite(entry.raw()) ||
                        entry.raw() != zero.raw()) {
                        return std::nullopt;
                    }
                }
            }
            Matrix<X> basis(n, dimension, zero);
            Matrix<X> restriction(dimension, dimension, zero);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != dimension; ++col) {
                    basis[row][col] = orthogonal[row][col];
                }
            }
            for (SizeType row = 0u; row != dimension; ++row) {
                for (SizeType col = 0u; col != dimension; ++col) {
                    restriction[row][col] = quasi_triangular[row][col];
                }
            }
            return SchurInvariantSubspaceData<P, PR>{
                std::move(basis), std::move(restriction)};
        }
    }

    //! \brief Reorder complete real Schur blocks by orthogonal similarity.
    //! \param order A permutation of the ORIGINAL diagonal block indices.
    //! \param precision Precision for the numerical 1x1/2x2 swaps.
    //! \return Independently owned reordered factors, or no value if
    //! any adjacent exchange is ill-conditioned or fails verification.
    //!
    //! \details A desired block order is realized by adjacent swaps.
    //! For adjacent blocks A (p-by-p), B (q-by-q) and coupling C,
    //! \f[ S=\begin{pmatrix}A&C\\0&B\end{pmatrix}, \f]
    //! solve the SMALL real Sylvester equation
    //! \f[ A X-X B=-C. \f]
    //! The graph columns of W=[X;I] span an S-invariant subspace
    //! corresponding to the eigenvalues of B. Orthonormalizing W
    //! and completing its columns gives a local orthogonal Z for
    //! which Z^T*S*Z has block order B,A. Apply Z as a similarity to
    //! the entire Schur matrix and accumulate the same factor in Q.
    //! The Sylvester equation has at most four unknowns, and all
    //! calculations are real even for 2x2 complex-conjugate blocks.
    //!
    //! \par Numerical acceptance
    //! Each adjacent exchange checks scaled Sylvester pivots and
    //! finite entries, orthogonalizes with two Gram-Schmidt passes,
    //! verifies 1x1/2x2 quasi-triangular structure, and zeros
    //! roundoff-only off-block fill. A final check verifies Q^T Q
    //! and preservation of the original Schur factorization.
    //! Nearly coincident block spectra can make Sylvester singular:
    //! this operation then returns no value rather than silently
    //! producing an unreliable invariant subspace. The original
    //! factors remain unchanged on any failure.
    //!
    //! \par Mathematical provenance
    //! See Bai and Demmel, "On Swapping Diagonal Blocks in Real Schur
    //! Form" (LAPACK Working Note 54), and LAPACK DTREXC for real
    //! Schur reordering by orthogonal adjacent block exchanges.
    //! This is a conservative SMALL-block Sylvester implementation,
    //! not LAPACK's more sophisticated canonicalization/conditioning.
    //! https://netlib.org/lapack/lawnspdf/lawn54.pdf
    //! https://netlib.org/lapack/explore-html/
    //!
    //! \see leading_invariant_subspace for extracting the selected
    //! leading subspace AFTER a successful reorder. ValidatedTag
    //! is not supported: floating Schur swaps are not certificates.
    std::optional<SchurDecompositionData<P, PR>>
    reorder_blocks(const std::vector<SizeType>& order, PR precision) const {
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return std::nullopt;
        } else {
            using Raw = typename X::RawType;
            const SizeType n = orthogonal.row_size();
            const SizeType count = block_sizes.size();
            if (n == 0u || orthogonal.column_size() != n ||
                quasi_triangular.row_size() != n ||
                quasi_triangular.column_size() != n ||
                count == 0u || order.size() != count) {
                return std::nullopt;
            }
            const X zero(0u, precision);
            const X one(1u, precision);
            const X eps(Raw::eps(precision));
            const auto finite = [](const X& x) {
                return is_finite(x.raw());
            };
            std::vector<bool> seen(count, false);
            SizeType covered = 0u;
            for (SizeType j = 0u; j != count; ++j) {
                const SizeType width = block_sizes[j];
                if ((width != 1u && width != 2u) || width > n - covered) {
                    return std::nullopt;
                }
                if (width == 2u) {
                    const SizeType k = covered;
                    const X discriminant =
                        sqr(quasi_triangular[k][k] -
                            quasi_triangular[k + 1u][k + 1u]) +
                        X(4u, precision) * quasi_triangular[k][k + 1u] *
                            quasi_triangular[k + 1u][k];
                    if (!finite(discriminant) ||
                        !(discriminant.raw() < zero.raw())) {
                        return std::nullopt;
                    }
                }
                covered += width;
            }
            if (covered != n) { return std::nullopt; }
            for (const SizeType desired: order) {
                if (desired >= count || seen[desired]) {
                    return std::nullopt;
                }
                seen[desired] = true;
            }
            X scale(one);
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    const X& q = orthogonal[row][col];
                    const X& t = quasi_triangular[row][col];
                    if (!finite(q) || !finite(t)) {
                        return std::nullopt;
                    }
                    const X magnitude = abs(t);
                    if (magnitude.raw() > scale.raw()) {
                        scale = magnitude;
                    }
                }
            }
            const X factor = X(static_cast<unsigned>(n * n), precision);
            const X residual_limit =
                X(32768u, precision) * eps * scale * factor;
            const X orthogonal_limit =
                X(32768u, precision) * eps * factor;
            const X structure_limit =
                X(2048u, precision) * eps * scale * factor;
            const X minimum_squared = sqr(X(128u, precision) * eps);
            if (!finite(residual_limit) || !finite(orthogonal_limit) ||
                !finite(structure_limit) || !finite(minimum_squared)) {
                return std::nullopt;
            }

            SchurDecompositionData reordered(*this);
            std::vector<SizeType> current(count);
            for (SizeType i = 0u; i != count; ++i) { current[i] = i; }

            // Swap the blocks at positions left and left+1.
            // This is never used by the ordinary eigensolver, only
            // upon an explicit request to select a Schur subspace.
            const auto exchange = [&](SizeType left) {
                SizeType start = 0u;
                for (SizeType i = 0u; i != left; ++i) {
                    start += reordered.block_sizes[i];
                }
                const SizeType p = reordered.block_sizes[left];
                const SizeType q = reordered.block_sizes[left + 1u];
                const SizeType total = p + q;
                const SizeType unknowns = p * q;
                Matrix<X>& T = reordered.quasi_triangular;
                Matrix<X>& Q = reordered.orthogonal;

                // Row-major Sylvester system for A*X-X*B=-C.
                std::vector<std::vector<X>> system(
                    unknowns, std::vector<X>(unknowns + 1u, zero));
                for (SizeType i = 0u; i != p; ++i) {
                    for (SizeType j = 0u; j != q; ++j) {
                        const SizeType row = i * q + j;
                        for (SizeType k = 0u; k != p; ++k) {
                            system[row][k * q + j] += T[start + i][start + k];
                        }
                        for (SizeType k = 0u; k != q; ++k) {
                            system[row][i * q + k] -=
                                T[start + p + k][start + p + j];
                        }
                        system[row][unknowns] =
                            -T[start + i][start + p + j];
                    }
                }
                X system_scale(one);
                for (SizeType row = 0u; row != unknowns; ++row) {
                    for (SizeType col = 0u; col <= unknowns; ++col) {
                        if (!finite(system[row][col])) { return false; }
                        if (col != unknowns &&
                            abs(system[row][col]).raw() > system_scale.raw()) {
                            system_scale = abs(system[row][col]);
                        }
                    }
                }
                const X pivot_limit =
                    X(128u, precision) * eps * system_scale;
                if (!finite(pivot_limit)) { return false; }

                // Small partial-pivoting Gaussian elimination, at
                // most four unknowns. Reject near-singular exchanges.
                for (SizeType col = 0u; col != unknowns; ++col) {
                    SizeType pivot = col;
                    X biggest(zero);
                    for (SizeType row = col; row != unknowns; ++row) {
                        const X magnitude = abs(system[row][col]);
                        if (!finite(magnitude)) { return false; }
                        if (magnitude.raw() > biggest.raw()) {
                            biggest = magnitude;
                            pivot = row;
                        }
                    }
                    if (biggest.raw() <= pivot_limit.raw()) { return false; }
                    if (pivot != col) {
                        std::swap(system[pivot], system[col]);
                    }
                    for (SizeType row = col + 1u; row != unknowns; ++row) {
                        const X multiplier =
                            system[row][col] / system[col][col];
                        if (!finite(multiplier)) { return false; }
                        for (SizeType j = col + 1u; j <= unknowns; ++j) {
                            system[row][j] -= multiplier * system[col][j];
                            if (!finite(system[row][j])) { return false; }
                        }
                        system[row][col] = zero;
                    }
                }
                std::vector<X> solution(unknowns, zero);
                for (SizeType back = unknowns; back > 0u;) {
                    --back;
                    X value(system[back][unknowns]);
                    for (SizeType j = back + 1u; j != unknowns; ++j) {
                        value -= system[back][j] * solution[j];
                    }
                    solution[back] = value / system[back][back];
                    if (!finite(solution[back])) { return false; }
                }

                // First q columns span graph([X;I]). Complete that
                // orthonormal basis in R^(p+q) with unit vectors.
                std::vector<std::vector<X>> basis;
                const auto orthogonalize = [&](std::vector<X> v,
                                               bool required) {
                    X maximum(zero);
                    for (const X& x: v) {
                        if (!finite(x)) { return false; }
                        if (abs(x).raw() > maximum.raw()) { maximum = abs(x); }
                    }
                    if (!(maximum.raw() > zero.raw())) { return !required; }
                    for (X& x: v) { x /= maximum; }
                    // Reorthogonalized modified Gram-Schmidt,
                    // on at most four-dimensional real vectors.
                    for (SizeType repeat = 0u; repeat != 2u; ++repeat) {
                        for (const auto& u: basis) {
                            X dot(zero);
                            for (SizeType row = 0u; row != total; ++row) {
                                dot += u[row] * v[row];
                            }
                            for (SizeType row = 0u; row != total; ++row) {
                                v[row] -= dot * u[row];
                            }
                        }
                    }
                    X squared(zero);
                    for (const X& x: v) { squared += sqr(x); }
                    if (!finite(squared)) { return false; }
                    if (squared.raw() <= minimum_squared.raw()) {
                        return !required;
                    }
                    const X magnitude = sqrt(squared);
                    if (!finite(magnitude) ||
                        !(magnitude.raw() > zero.raw())) {
                        return false;
                    }
                    for (X& x: v) { x /= magnitude; }
                    basis.push_back(std::move(v));
                    return true;
                };
                for (SizeType col = 0u; col != q; ++col) {
                    std::vector<X> v(total, zero);
                    for (SizeType row = 0u; row != p; ++row) {
                        v[row] = solution[row * q + col];
                    }
                    v[p + col] = one;
                    if (!orthogonalize(std::move(v), true)) { return false; }
                }
                for (SizeType col = 0u; col != total; ++col) {
                    if (basis.size() == total) { break; }
                    std::vector<X> v(total, zero);
                    v[col] = one;
                    if (!orthogonalize(std::move(v), false)) { return false; }
                }
                if (basis.size() != total) { return false; }

                Matrix<X> step = Matrix<X>::identity(n, precision);
                for (SizeType col = 0u; col != total; ++col) {
                    for (SizeType row = 0u; row != total; ++row) {
                        step[start + row][start + col] = basis[col][row];
                    }
                }
                T = transpose(step) * (T * step);
                Q = Q * step;
                std::swap(reordered.block_sizes[left],
                          reordered.block_sizes[left + 1u]);

                // Preserve only 1x1 and conjugate 2x2 blocks. Floating
                // fill outside those diagonal blocks may be discarded
                // ONLY below the conservative structure threshold.
                std::vector<bool> allowed(n, false);
                SizeType offset = 0u;
                for (const SizeType width: reordered.block_sizes) {
                    if (width == 2u) {
                        allowed[offset] = true;
                        const X discriminant =
                            sqr(T[offset][offset] -
                                T[offset + 1u][offset + 1u]) +
                            X(4u, precision) * T[offset][offset + 1u] *
                                T[offset + 1u][offset];
                        if (!finite(discriminant) ||
                            !(discriminant.raw() < zero.raw())) {
                            return false;
                        }
                    }
                    offset += width;
                }
                if (offset != n) { return false; }
                for (SizeType col = 0u; col != n; ++col) {
                    for (SizeType row = 0u; row != n; ++row) {
                        if (!finite(T[row][col]) || !finite(Q[row][col])) {
                            return false;
                        }
                        if (row <= col ||
                            (row == col + 1u && allowed[col])) {
                            continue;
                        }
                        if (abs(T[row][col]).raw() >
                            structure_limit.raw()) {
                            return false;
                        }
                        T[row][col] = zero;
                    }
                }
                return true;
            };

            for (SizeType target = 0u; target != count; ++target) {
                SizeType from = target;
                while (from < count && current[from] != order[target]) {
                    ++from;
                }
                if (from == count) { return std::nullopt; }
                while (from > target) {
                    if (!exchange(from - 1u)) {
                        return std::nullopt;
                    }
                    std::swap(current[from - 1u], current[from]);
                    --from;
                }
            }

            // Reordering must preserve the ORIGINAL Schur-factor
            // reconstruction, not merely produce an orthogonal Q.
            const Matrix<X> old_reconstruction =
                orthogonal * quasi_triangular * transpose(orthogonal);
            const Matrix<X> new_action =
                old_reconstruction * reordered.orthogonal;
            const Matrix<X> new_represented =
                reordered.orthogonal * reordered.quasi_triangular;
            const Matrix<X> new_gram =
                transpose(reordered.orthogonal) * reordered.orthogonal;
            for (SizeType row = 0u; row != n; ++row) {
                for (SizeType col = 0u; col != n; ++col) {
                    const X error = abs(new_action[row][col] -
                                        new_represented[row][col]);
                    const X ortho = abs(new_gram[row][col] -
                        (row == col ? one : zero));
                    if (!finite(error) || !finite(ortho) ||
                        error.raw() > residual_limit.raw() ||
                        ortho.raw() > orthogonal_limit.raw()) {
                        return std::nullopt;
                    }
                }
            }
            return reordered;
        }
    }


    //! \brief Extract a real invariant subspace for selected Schur blocks.
    //! \param selected_blocks Indices of complete diagonal blocks in
    //! the ORIGINAL Schur ordering, without duplicates.
    //! \param precision Arithmetic precision for block exchanges.
    //! \return An orthonormal n-by-k basis and k-by-k restriction,
    //! or no result if selection/reordering is invalid or unstable.
    //!
    //! \details The selected blocks are moved to the leading positions
    //! in the supplied order. Unselected blocks retain their relative
    //! order. The exchange uses the orthogonal Sylvester-based Schur
    //! reordering and then leading_invariant_subspace; it never
    //! separates a conjugate real 2x2 block or assumes that A is
    //! diagonalizable. The original Schur factors are untouched.
    //! See reorder_blocks and LAPACK DTREXC (Bai-Demmel, LAWN 54)
    //! for the mathematical mechanism and its limitations.
    std::optional<SchurInvariantSubspaceData<P, PR>>
    selected_invariant_subspace(
            const std::vector<SizeType>& selected_blocks,
            PR precision) const {
        const SizeType count = block_sizes.size();
        if (selected_blocks.empty() || selected_blocks.size() > count) {
            return std::nullopt;
        }
        std::vector<bool> used(count, false);
        std::vector<SizeType> permutation;
        permutation.reserve(count);
        for (const SizeType index: selected_blocks) {
            if (index >= count || used[index]) {
                return std::nullopt;
            }
            used[index] = true;
            permutation.push_back(index);
        }
        for (SizeType i = 0u; i != count; ++i) {
            if (!used[i]) { permutation.push_back(i); }
        }
        const auto reordered = this->reorder_blocks(permutation, precision);
        if (!reordered) { return std::nullopt; }
        return reordered->leading_invariant_subspace(
            selected_blocks.size(), precision);
    }

    //! \brief Select real invariant subspaces inside a conjugation-closed disk.
    //! \param center Real center of a disk in the complex plane.
    //! \param radius Strictly positive real radius.
    //! \param precision Precision for boundary tests and block swaps.
    //! \return A real invariant subspace for every Schur block strictly
    //! inside the disk, or no value if the selection is empty,
    //! intersects a numerically uncertain boundary, or cannot be
    //! stably reordered.
    //!
    //! \details The eigenvalues in a scalar Schur block are real,
    //! and a nonreal 2x2 real Schur block contains a conjugate pair.
    //! A disk centered ON THE REAL AXIS is invariant under complex
    //! conjugation, so it can never select just one member of a
    //! conjugate pair. For a 2x2 block [a b;c d] with negative
    //! discriminant Delta=(a-d)^2+4bc, the two eigenvalues have
    //! real part (a+d)/2 and squared imaginary magnitude -Delta/4.
    //! Both have the SAME squared distance from the real center r:
    //! \f[ d^2 = ((a+d)/2-r)^2-\Delta/4. \f]
    //! For a scalar block lambda, d^2=(lambda-r)^2.
    //! Classify complete blocks by d^2 compared with radius^2;
    //! if they approach the boundary within a scale/precision
    //! tolerance, report INCONCLUSIVE instead of guessing inclusion.
    //! Selected blocks are then reordered before extraction.
    //!
    //! This is an APPROXIMATE spectral selection, not a verified
    //! eigenvalue-count or contour certificate. In particular the
    //! chosen radius must be separated from the approximate spectrum.
    //! Real-only arithmetic is sufficient, and this method is
    //! explicitly requested rather than invoked by the real fast path.
    //! See Wilkinson (1965), Chapters 6 and 8, and the Schur
    //! reordering literature cited by reorder_blocks.
    std::optional<SchurInvariantSubspaceData<P, PR>>
    invariant_subspace_in_real_disk(
            const X& center, const X& radius, PR precision) const {
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return std::nullopt;
        } else {
            using Raw = typename X::RawType;
            const SizeType n = quasi_triangular.row_size();
            if (n == 0u || quasi_triangular.column_size() != n ||
                orthogonal.row_size() != n ||
                orthogonal.column_size() != n ||
                block_sizes.empty() ||
                !is_finite(center.raw()) ||
                !is_finite(radius.raw())) {
                return std::nullopt;
            }
            const X zero(0u, precision);
            const X one(1u, precision);
            const X two(2u, precision);
            const X four(4u, precision);
            if (!(radius.raw() > zero.raw())) {
                return std::nullopt;
            }
            X scale(one);
            if (abs(center).raw() > scale.raw()) { scale = abs(center); }
            if (radius.raw() > scale.raw()) { scale = radius; }
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != n; ++j) {
                    const X entry = quasi_triangular[i][j];
                    if (!is_finite(entry.raw())) { return std::nullopt; }
                    if (abs(entry).raw() > scale.raw()) {
                        scale = abs(entry);
                    }
                }
            }
            const X eps(Raw::eps(precision));
            const X radius_squared = sqr(radius);
            const X tolerance = X(512u, precision) * eps *
                sqr(scale) * X(static_cast<unsigned>(n), precision);
            if (!is_finite(radius_squared.raw()) ||
                !is_finite(tolerance.raw())) {
                return std::nullopt;
            }
            SizeType offset = 0u;
            std::vector<SizeType> selected_blocks;
            selected_blocks.reserve(block_sizes.size());
            for (SizeType i = 0u; i != block_sizes.size(); ++i) {
                const SizeType width = block_sizes[i];
                if ((width != 1u && width != 2u) ||
                    width > n - offset) {
                    return std::nullopt;
                }
                X distance_squared(zero);
                if (width == 1u) {
                    distance_squared = sqr(
                        quasi_triangular[offset][offset] - center);
                } else {
                    const X a = quasi_triangular[offset][offset];
                    const X b = quasi_triangular[offset][offset + 1u];
                    const X c = quasi_triangular[offset + 1u][offset];
                    const X d = quasi_triangular[offset + 1u][offset + 1u];
                    const X discriminant = sqr(a - d) + four * b * c;
                    if (!is_finite(discriminant.raw()) ||
                        !(discriminant.raw() < zero.raw())) {
                        return std::nullopt;
                    }
                    const X real_part = (a + d) / two;
                    distance_squared = sqr(real_part - center) -
                        discriminant / four;
                }
                if (!is_finite(distance_squared.raw())) {
                    return std::nullopt;
                }
                const X boundary_gap = abs(distance_squared - radius_squared);
                if (!is_finite(boundary_gap.raw()) ||
                    boundary_gap.raw() <= tolerance.raw()) {
                    return std::nullopt;
                }
                if (distance_squared.raw() < radius_squared.raw()) {
                    selected_blocks.push_back(i);
                }
                offset += width;
            }
            if (offset != n || selected_blocks.empty()) {
                return std::nullopt;
            }
            return this->selected_invariant_subspace(
                selected_blocks, precision);
        }
    }

};

//! \brief Status and optional factors on SUCCESS.
template<class P, class PR>
struct SchurDecompositionResult {
    SchurStatus status;
    std::optional<SchurDecompositionData<P, PR>> data;
};

//! \brief Approximate a real Schur decomposition of a real square matrix.
//!
//! \details The real Schur theorem gives an orthogonal Q and a real
//! quasi-upper-triangular T such that
//! \f[ A = Q T Q^{\mathsf T}, \qquad Q^{\mathsf T}Q=I. \f]
//! Scalar diagonal blocks encode real eigenvalues; irreducible 2x2
//! diagonal blocks encode conjugate nonreal eigenvalues. Schur reduction
//! does not assume diagonalizability and does not produce eigenvectors.
//!
//! \par Step 1: Householder reflections and Hessenberg reduction
//! A Householder reflector is an orthogonal, symmetric matrix
//! \f[ U = I - 2vv^{\mathsf T}/(v^{\mathsf T}v). \f]
//! Choosing v from the portion of column k below the diagonal makes
//! the components below the first subdiagonal vanish after applying
//! the similarity T <- U*T*U, while accumulating Q <- Q*U. Repeating
//! this for k=0,...,n-3 yields an upper Hessenberg matrix H, i.e.
//! H(i,j)=0 when i>j+1. Orthogonal similarity preserves eigenvalues
//! and is preferable to a nonorthogonal triangularization because
//! it does not amplify errors through an ill-conditioned basis.
//! The reflection construction derives from Householder's 1958
//! orthogonal triangularization; its use for Hessenberg eigenvalue
//! reduction is described by Wilkinson (1965, Chapter 6).
//!
//! \par Step 2: shifted QR transformations
//! For a real shift mu, form H-mu*I=U*R by Householder QR, then set
//! \f[ H_{\mathrm{new}} = R U+\mu I = U^{\mathsf T} H U. \f]
//! The shift is chosen as the real eigenvalue of the trailing 2x2
//! block closest to its bottom-right diagonal entry (a Wilkinson-type
//! shift). When that 2x2 has a negative discriminant, its eigenvalues
//! are a conjugate pair and one instead forms the REAL quadratic
//! polynomial
//! \f[ p(H)=H^2-\operatorname{tr}(B)H+\det(B)I, \f]
//! where B is the trailing 2x2 block. The implicit Q theorem of
//! Francis (1961-1962) allows its QR similarity to be constructed
//! from the first column of p(H), without forming p(H) explicitly.
//! For upper Hessenberg H, its first three possibly nonzero entries are
//! \f[ x=h_{00}^2+h_{01}h_{10}-s h_{00}+t \f]
//! \f[ y=h_{10}(h_{00}+h_{11}-s) \f]
//! \f[ z=h_{10}h_{21} \f]
//! \f[ s=\operatorname{tr}(B),\qquad t=\det(B) \f]
//! A Householder reflector maps (x,y,z) to the first coordinate,
//! introducing a small out-of-band bulge. Subsequent reflectors on
//! at most three consecutive rows chase the bulge down the active
//! Hessenberg subdiagonal. Applying the full orthogonal similarities
//! and accumulating them in Q implements the implicit real Francis
//! double shift, without explicitly squaring H. Aggressive deflation
//! and multishifts remain future work.
//!
//! \par Step 3: deflation and acceptance criteria
//! A sufficiently small subdiagonal entry is set to zero using a
//! precision/scale-dependent tolerance, splitting off a real 1x1
//! block or a nonreal 2x2 block. A final 2x2 with a real spectrum
//! is converted to upper triangular form by an orthogonal rotation
//! built from a real eigenvector; this also permits a defective
//! repeated real root without requiring a diagonalization.
//! Iterations are bounded. Before
//! SUCCESS, finite factors, the 1x1/2x2 block structure, approximate
//! orthogonality and the ORIGINAL-input residual A*Q-Q*T are checked.
//! The tolerances are numerical heuristics, not proof of a backward
//! error bound in the ValidatedTag sense. ApproximateTag can return
//! INCONCLUSIVE on clusters, overflows and lack of convergence;
//! ValidatedTag currently always returns INCONCLUSIVE.
//!
//! \par Historical sources
//! - A. S. Householder, "Unitary Triangularization of a Nonsymmetric
//!   Matrix", Journal of the ACM 5(4), 339-342 (1958).
//!   https://doi.org/10.1145/320941.320947
//! - J. G. F. Francis, "The QR Transformation", Parts 1 and 2,
//!   The Computer Journal 4(3), 265-271 (1961) and 4(4), 332-345 (1962).
//!   https://doi.org/10.1093/comjnl/4.3.265
//!   https://doi.org/10.1093/comjnl/4.4.332
//! - J. H. Wilkinson, The Algebraic Eigenvalue Problem, Clarendon
//!   Press, Oxford (1965), Chapters 6 and 8 (Hessenberg and QR).
//!
//! \see Spectrum for the complete eigenvalue multiset, including
//! defective matrices. SchurDecompositionData::leading_invariant_subspace
//! extracts leading real invariant subspaces without diagonalization;
//! EigenDecomposition additionally requires a full real or complex basis.
template<class P, class PR>
class SchurDecomposition {
  public:
    using X = FloatType<P, PR>;
    using Result = SchurDecompositionResult<P, PR>;
    //! \brief Select the floating-point precision used by the QR iteration.
    explicit SchurDecomposition(PR pr) : _precision(pr) {}

    //! \brief Attempt approximate real Schur reduction of A.
    //! \param A Real square matrix (non-square input violates precondition).
    //! \return Schur factors only on SUCCESS. Zero order is unsupported;
    //! the validated paradigm is deliberately inconclusive.
    //! \details First uses orthogonal Householder similarities to obtain
    //! Hessenberg form, then performs shifted QR with deflation and
    //! checks A*Q-Q*T and Q^T*Q-I. See the class documentation for
    //! derivation, numerical limitations and historical references.
    Result operator()(const Matrix<X>& A) const {
        ARIADNE_PRECONDITION(A.row_size() == A.column_size());
        const SizeType n = A.row_size();
        if (n == 0u) {
            return {SchurStatus::UNSUPPORTED_DIMENSION, std::nullopt};
        }
        if constexpr (!std::is_same_v<P, ApproximateTag>) {
            return {SchurStatus::INCONCLUSIVE, std::nullopt};
        } else {
            using Raw = typename X::RawType;
            const X zero(0u, _precision);
            const X one(1u, _precision);
            const X two_value(2u, _precision);
            const X eps(Raw::eps(_precision));
            X scale(one);
            if (!this->_finite_matrix(A)) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != n; ++j) {
                    const X mag = abs(A[i][j]);
                    if (mag.raw() > scale.raw()) { scale = mag; }
                }
            }
            const X tolerance = X(128u, _precision) * eps * scale;
            const X residual_limit = X(32768u, _precision) * eps * scale *
                X(static_cast<unsigned>(n * n), _precision);
            const X orthogonal_limit = X(32768u, _precision) * eps *
                X(static_cast<unsigned>(n * n), _precision);
            if (!this->_finite(tolerance) || !this->_finite(residual_limit) ||
                !this->_finite(orthogonal_limit)) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            Matrix<X> T(A);
            Matrix<X> Q = Matrix<X>::identity(n, _precision);

            /** \brief Householder reduction to upper Hessenberg form.
             *  \details Each reflector annihilates the tail of a column
             *  by choosing its sign to avoid cancellation. Apply the
             *  reflector on both sides (orthogonal similarity), and
             *  multiply it into Q. See Householder (1958), Wilkinson (1965).
             */
            for (SizeType k = 0u; k + 2u < n; ++k) {
                X norm(zero);
                for (SizeType i = k + 1u; i != n; ++i) {
                    norm += sqr(T[i][k]);
                }
                if (!this->_finite(norm)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                if (norm.raw() == zero.raw()) { continue; }
                const X mag = sqrt(norm);
                const X alpha = T[k + 1u][k].raw() < zero.raw() ? mag : -mag;
                std::vector<X> v(n - k - 1u, zero);
                v[0u] = T[k + 1u][k] - alpha;
                X vnorm(sqr(v[0u]));
                for (SizeType i = k + 2u; i != n; ++i) {
                    v[i - k - 1u] = T[i][k];
                    vnorm += sqr(v[i - k - 1u]);
                }
                if (!this->_finite(vnorm) || !(vnorm.raw() > zero.raw())) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                const X beta = two_value / vnorm;
                for (SizeType j = k; j != n; ++j) {
                    X dot(zero);
                    for (SizeType i = k + 1u; i != n; ++i) {
                        dot += v[i - k - 1u] * T[i][j];
                    }
                    for (SizeType i = k + 1u; i != n; ++i) {
                        T[i][j] -= beta * dot * v[i - k - 1u];
                    }
                }
                for (SizeType i = 0u; i != n; ++i) {
                    X dot(zero);
                    for (SizeType j = k + 1u; j != n; ++j) {
                        dot += T[i][j] * v[j - k - 1u];
                    }
                    for (SizeType j = k + 1u; j != n; ++j) {
                        T[i][j] -= beta * dot * v[j - k - 1u];
                    }
                    X qdot(zero);
                    for (SizeType j = k + 1u; j != n; ++j) {
                        qdot += Q[i][j] * v[j - k - 1u];
                    }
                    for (SizeType j = k + 1u; j != n; ++j) {
                        Q[i][j] -= beta * qdot * v[j - k - 1u];
                    }
                }
            }
            if (!this->_finite_matrix(T) || !this->_finite_matrix(Q)) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }

            SizeType active = n;
            const SizeType limit = 512u * n * n;
            for (SizeType iter = 0u; active != 0u && iter != limit; ++iter) {
                // The trailing subdiagonal exists only for active >= 2.
                // In particular, never form active-2 for a singleton:
                // SizeType is unsigned and would wrap to its maximum.
                if (active < 2u) {
                    active = 0u;
                    break;
                }
                if (abs(T[active - 1u][active - 2u]).raw() <= tolerance.raw()) {
                    T[active - 1u][active - 2u] = zero;
                    --active;
                    continue;
                }
                const X& a = T[active - 2u][active - 2u];
                const X& b = T[active - 2u][active - 1u];
                const X& c = T[active - 1u][active - 2u];
                const X& d = T[active - 1u][active - 1u];
                const X disc = sqr(a - d) + X(4u, _precision) * b * c;
                if (!this->_finite(disc)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                const bool pair = disc.raw() < zero.raw();
                // At order two, a real Schur block can be finished
                // directly: with nonreal roots keep the 2x2 block;
                // with real roots rotate a real eigenvector to e_1.
                // This also handles a defective repeated real root
                // without demanding two independent eigenvectors.
                if (active == 2u) {
                    if (pair) {
                        active = 0u;
                        break;
                    }
                    const X root = sqrt(disc);
                    const X lambda = (a + d + root) / two_value;
                    // Both vectors are in ker(B-lambda*I) in exact
                    // arithmetic. Choose the larger to reduce loss
                    // through cancellation before normalization.
                    const X first_x = b;
                    const X first_y = lambda - a;
                    const X second_x = lambda - d;
                    const X second_y = c;
                    const X first_norm = sqr(first_x) + sqr(first_y);
                    const X second_norm = sqr(second_x) + sqr(second_y);
                    if (!this->_finite(first_norm) || !this->_finite(second_norm)) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    const bool use_first = first_norm.raw() > second_norm.raw();
                    const X vx = use_first ? first_x : second_x;
                    const X vy = use_first ? first_y : second_y;
                    const X squared = use_first ? first_norm : second_norm;
                    if (!(squared.raw() > zero.raw())) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    const X magnitude = sqrt(squared);
                    const X cs = vx / magnitude;
                    const X sn = vy / magnitude;
                    if (!this->_finite(cs) || !this->_finite(sn)) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    Matrix<X> final_step = Matrix<X>::identity(n, _precision);
                    final_step[0u][0u] = cs;
                    final_step[1u][0u] = sn;
                    final_step[0u][1u] = -sn;
                    final_step[1u][1u] = cs;
                    // Apply the similarity to the full matrix so that
                    // couplings to previously deflated blocks rotate too.
                    T = transpose(final_step) * (T * final_step);
                    Q = Q * final_step;
                    if (!this->_finite_matrix(T) || !this->_finite_matrix(Q) ||
                        abs(T[1u][0u]).raw() > tolerance.raw()) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    T[1u][0u] = zero;
                    active = 0u;
                    break;
                }
                if (pair && active > 2u &&
                    abs(T[active - 2u][active - 3u]).raw() <= tolerance.raw()) {
                    T[active - 2u][active - 3u] = zero;
                    active -= 2u;
                    continue;
                }

                if (pair) {
                    // Francis implicit real double shift: only the
                    // first column of p(T) is required to start bulge
                    // chasing, rather than explicitly forming T^2.
                    if (!this->_francis_double_shift(T, Q, active)) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                } else {
                    // Wilkinson-type single shift: choose the nearer
                    // real root of the trailing 2x2 characteristic
                    // polynomial. See Wilkinson (1965), Chapter 8.
                    const X root = sqrt(disc);
                    const X low = (a + d - root) / two_value;
                    const X high = (a + d + root) / two_value;
                    const X shift = abs(low - d).raw() < abs(high - d).raw()
                        ? low : high;
                    Matrix<X> shifted(active, active, zero);
                    for (SizeType i = 0u; i != active; ++i) {
                        for (SizeType j = 0u; j != active; ++j) {
                            shifted[i][j] = T[i][j] -
                                (i == j ? shift : zero);
                        }
                    }
                    const auto small_Q = this->_qr(shifted);
                    if (!small_Q) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    Matrix<X> step = Matrix<X>::identity(n, _precision);
                    for (SizeType i = 0u; i != active; ++i) {
                        for (SizeType j = 0u; j != active; ++j) {
                            step[i][j] = (*small_Q)[i][j];
                        }
                    }
                    // Full similarity includes previously deflated
                    // trailing Schur blocks and their couplings.
                    T = transpose(step) * (T * step);
                    Q = Q * step;
                }
                if (!this->_finite_matrix(T) || !this->_finite_matrix(Q)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                for (SizeType col = 0u; col + 2u < n; ++col) {
                    for (SizeType row = col + 2u; row != n; ++row) {
                        if (abs(T[row][col]).raw() > tolerance.raw()) {
                            return {SchurStatus::INCONCLUSIVE, std::nullopt};
                        }
                        T[row][col] = zero;
                    }
                }
            }
            if (active != 0u) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            // A Schur partition must exhaust the matrix with disjoint
            // 1x1/2x2 blocks, without probing beyond its final row.
            if (T.row_size() != n || T.column_size() != n) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }
            std::vector<SizeType> blocks;
            SizeType covered = 0u;
            while (covered < n) {
                const SizeType remaining = n - covered;
                if (remaining == 1u) {
                    blocks.push_back(1u);
                    ++covered;
                    break;
                }
                const SizeType i = covered;
                const X first_subdiagonal = abs(T[i + 1u][i]);
                if (!this->_finite(first_subdiagonal)) {
                    return {SchurStatus::INCONCLUSIVE, std::nullopt};
                }
                if (first_subdiagonal.raw() > tolerance.raw()) {
                    if (remaining > 2u) {
                        const X second_subdiagonal = abs(T[i + 2u][i + 1u]);
                        if (!this->_finite(second_subdiagonal) ||
                            second_subdiagonal.raw() > tolerance.raw()) {
                            return {SchurStatus::INCONCLUSIVE, std::nullopt};
                        }
                    }
                    const X disc = sqr(T[i][i] - T[i + 1u][i + 1u]) +
                        X(4u, _precision) * T[i][i + 1u] * T[i + 1u][i];
                    if (!this->_finite(disc) || !(disc.raw() < zero.raw())) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                    blocks.push_back(2u);
                    covered += 2u;
                } else {
                    T[i + 1u][i] = zero;
                    blocks.push_back(1u);
                    ++covered;
                }
            }
            if (covered != n) {
                return {SchurStatus::INCONCLUSIVE, std::nullopt};
            }

            const Matrix<X> A_Q = A * Q;
            const Matrix<X> Q_T = Q * T;
            const Matrix<X> Q_Q = transpose(Q) * Q;
            for (SizeType i = 0u; i != n; ++i) {
                for (SizeType j = 0u; j != n; ++j) {
                    const X error = abs(A_Q[i][j] - Q_T[i][j]);
                    const X orth_error = abs(Q_Q[i][j] -
                        (i == j ? one : zero));
                    if (!this->_finite(error) || !this->_finite(orth_error) ||
                        error.raw() > residual_limit.raw() ||
                        orth_error.raw() > orthogonal_limit.raw()) {
                        return {SchurStatus::INCONCLUSIVE, std::nullopt};
                    }
                }
            }
            return {SchurStatus::SUCCESS, SchurDecompositionData<P, PR>{
                std::move(Q), std::move(T), std::move(blocks)}};
        }
    }

  private:
    bool _finite(const X& value) const {
        if constexpr (std::is_same_v<P, ApproximateTag>) {
            return is_finite(value.raw());
        } else {
            return is_finite(value.lower_raw()) && is_finite(value.upper_raw());
        }
    }

    bool _finite_matrix(const Matrix<X>& A) const {
        for (SizeType i = 0u; i != A.row_size(); ++i) {
            for (SizeType j = 0u; j != A.column_size(); ++j) {
                if (!this->_finite(A[i][j])) { return false; }
            }
        }
        return true;
    }

    //! \brief Perform a real Francis implicit double-shift QR step.
    //! \param T Full Hessenberg/Schur iterate, updated by similarity.
    //! \param Q Accumulated orthogonal transformation, updated in place.
    //! \param active Order of the leading unreduced Hessenberg block.
    //! \return False if a reflector or intermediate result is nonfinite.
    //! \details Write p(H)=H^2-tr(B)H+det(B)I for the trailing 2x2
    //! block B. The implicit Q theorem allows the QR similarity of
    //! p(H) to be obtained from just its first column: exactly three
    //! potentially nonzero entries for upper Hessenberg H. A first
    //! Householder reflector initiates a bulge, then reflectors of
    //! length two or three chase the fill below the subdiagonal.
    //! Every reflector acts on the full matrix T and accumulates into
    //! Q so couplings to deflated trailing blocks remain consistent.
    //! No full quadratic shift matrix is constructed.
    //! See Francis (1961-1962), The QR Transformation, Parts 1 and 2;
    //! Wilkinson (1965), The Algebraic Eigenvalue Problem, Chapter 8.
    bool _francis_double_shift(Matrix<X>& T, Matrix<X>& Q,
                               SizeType active) const {
        const SizeType n = T.row_size();
        if (active < 3u || active > n || T.column_size() != n ||
            Q.row_size() != n || Q.column_size() != n) {
            return false;
        }
        const X zero(0u, _precision);
        const X two_value(2u, _precision);
        const X& a = T[active - 2u][active - 2u];
        const X& b = T[active - 2u][active - 1u];
        const X& c = T[active - 1u][active - 2u];
        const X& d = T[active - 1u][active - 1u];
        const X trace = a + d;
        const X determinant = a * d - b * c;

        // First column of the real quadratic-shift polynomial.
        std::vector<X> first(3u, zero);
        first[0u] = T[0u][0u] * (T[0u][0u] - trace) +
            T[0u][1u] * T[1u][0u] + determinant;
        first[1u] = T[1u][0u] * (T[0u][0u] + T[1u][1u] - trace);
        first[2u] = T[1u][0u] * T[2u][1u];
        for (const auto& value: first) {
            if (!this->_finite(value)) {
                return false;
            }
        }

        for (SizeType k = 0u; k + 1u < active; ++k) {
            const SizeType width = k + 2u < active ? 3u : 2u;
            std::vector<X> v(width, zero);
            for (SizeType i = 0u; i != width; ++i) {
                v[i] = k == 0u ? first[i] : T[k + i][k - 1u];
            }
            X norm_squared(zero);
            for (const auto& value: v) {
                norm_squared += sqr(value);
            }
            if (!this->_finite(norm_squared)) {
                return false;
            }
            if (norm_squared.raw() == zero.raw()) {
                // No reflector needed for a zero bulge; the final
                // Hessenberg, residual and orthogonality checks apply.
                continue;
            }
            const X magnitude = sqrt(norm_squared);
            const X alpha = v[0u].raw() < zero.raw() ? magnitude : -magnitude;
            v[0u] -= alpha;
            X squared(sqr(v[0u]));
            for (SizeType i = 1u; i != width; ++i) {
                squared += sqr(v[i]);
            }
            if (!this->_finite(squared) || !(squared.raw() > zero.raw())) {
                return false;
            }
            const X beta = two_value / squared;

            // Left Householder factor on T.
            for (SizeType col = 0u; col != n; ++col) {
                X dot(zero);
                for (SizeType i = 0u; i != width; ++i) {
                    dot += v[i] * T[k + i][col];
                }
                for (SizeType i = 0u; i != width; ++i) {
                    T[k + i][col] -= beta * dot * v[i];
                }
            }
            // Right factor on T and accumulation of Q.
            for (SizeType row = 0u; row != n; ++row) {
                X dot(zero);
                for (SizeType i = 0u; i != width; ++i) {
                    dot += T[row][k + i] * v[i];
                }
                for (SizeType i = 0u; i != width; ++i) {
                    T[row][k + i] -= beta * dot * v[i];
                }
                X qdot(zero);
                for (SizeType i = 0u; i != width; ++i) {
                    qdot += Q[row][k + i] * v[i];
                }
                for (SizeType i = 0u; i != width; ++i) {
                    Q[row][k + i] -= beta * qdot * v[i];
                }
            }
        }
        return this->_finite_matrix(T) && this->_finite_matrix(Q);
    }

    //! \brief Orthogonal factor of a full Householder QR factorization.
    //! \details Repeated reflections of the form
    //! \f$U=I-2vv^{\mathsf T}/(v^{\mathsf T}v)\f$ reduce the real
    //! single-shift matrix to triangular form R. The accumulated Q satisfies input=Q*R
    //! approximately. A rank-deficient column with zero norm is
    //! skipped, not normalized. The reflector sign reduces cancellation.
    //! This is Householder QR, not classical Gram-Schmidt; see
    //! Householder (1958) and Wilkinson (1965), Chapter 4.
    //! \return Orthogonal Q, or no value if an intermediate value
    //! is not finite or a reflector cannot be formed.
    std::optional<Matrix<X>> _qr(const Matrix<X>& input) const {
        const SizeType n = input.row_size();
        const X zero(0u, _precision);
        const X two_value(2u, _precision);
        Matrix<X> R(input);
        Matrix<X> Q = Matrix<X>::identity(n, _precision);
        for (SizeType k = 0u; k + 1u < n; ++k) {
            X norm(zero);
            for (SizeType i = k; i != n; ++i) {
                norm += sqr(R[i][k]);
            }
            if (!this->_finite(norm)) { return std::nullopt; }
            if (norm.raw() == zero.raw()) { continue; }
            const X mag = sqrt(norm);
            const X alpha = R[k][k].raw() < zero.raw() ? mag : -mag;
            std::vector<X> v(n - k, zero);
            v[0u] = R[k][k] - alpha;
            X vnorm(sqr(v[0u]));
            for (SizeType i = k + 1u; i != n; ++i) {
                v[i - k] = R[i][k];
                vnorm += sqr(v[i - k]);
            }
            if (!this->_finite(vnorm) || !(vnorm.raw() > zero.raw())) {
                return std::nullopt;
            }
            const X beta = two_value / vnorm;
            for (SizeType j = k; j != n; ++j) {
                X dot(zero);
                for (SizeType i = k; i != n; ++i) { dot += v[i - k] * R[i][j]; }
                for (SizeType i = k; i != n; ++i) {
                    R[i][j] -= beta * dot * v[i - k];
                }
            }
            for (SizeType i = 0u; i != n; ++i) {
                X dot(zero);
                for (SizeType j = k; j != n; ++j) { dot += Q[i][j] * v[j - k]; }
                for (SizeType j = k; j != n; ++j) {
                    Q[i][j] -= beta * dot * v[j - k];
                }
            }
        }
        if (!this->_finite_matrix(Q)) { return std::nullopt; }
        return Q;
    }

    PR _precision;
};

} // namespace Ariadne

#endif // ARIADNE_SCHUR_DECOMPOSITION_HPP
