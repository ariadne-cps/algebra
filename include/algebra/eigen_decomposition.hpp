/***************************************************************************
 *            algebra/eigen_decomposition.hpp
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

/** @file algebra/eigen_decomposition.hpp
 *  @brief Approximate and validated real eigenvalue decompositions.
 */

#ifndef ARIADNE_EIGEN_DECOMPOSITION_HPP
#define ARIADNE_EIGEN_DECOMPOSITION_HPP

#include <chrono>
#include <optional>

#include "algebra/diagonal_matrix.hpp"
#include "algebra/matrix.hpp"
#include "foundation/paradigm.hpp"
#include "numeric/float.decl.hpp"

namespace Ariadne {

//! \brief Outcome of a real eigendecomposition attempt.
enum class EigenDecompositionStatus {
    SUCCESS,              //!< All real eigenpairs and an invertible eigenvector basis were obtained.
    INCONCLUSIVE,         //!< The requested real diagonalization was not established.
    UNSUPPORTED_DIMENSION //!< Unsupported structure or dimension for the requested paradigm.
};

//! \brief Algorithm attempted or successfully used by a decomposition.
enum class EigenDecompositionStrategy {
    NONE,          //!< No strategy attempted.
    ANALYTIC,      //!< 1x1/2x2 closed form.
    TRIANGULAR,    //!< Triangular eigenvector substitution.
    BLOCK_DIAGONAL, //!< Eigenbasis assembled from independent blocks.
    JACOBI,        //!< Approximate symmetric Jacobi rotations.
    SHIFTED_QR,    //!< Approximate nonsymmetric real QR and Newton refinement.
    GERSHGORIN,    //!< Direct validated disk separation.
    SIMILARITY,    //!< Verified similarity plus disk separation.
    KRAWCZYK      //!< Verified eigenpair root inclusion.
};

//! \brief Optional per-call performance diagnostics for benchmark analysis.
//! \details Timings are collected only when explicitly requested.
//! The two disjoint counters cover approximate candidate generation and
//! interval verification. Algorithm dispatch and minor bookkeeping are
//! included in neither. Timing incurs overhead and is not a certification.
//! On failure, strategy identifies the last attempted algorithm.
struct EigenDecompositionDiagnostics {
    EigenDecompositionStrategy strategy = EigenDecompositionStrategy::NONE;
    std::chrono::nanoseconds candidate_time{0};
    std::chrono::nanoseconds verification_time{0};
};

//! \brief Eigenvalues (in ascending order) and eigenvectors stored by columns.
//! \details For validated coefficients, each entry encloses the corresponding
//! eigenvalue or normalized eigenvector component for every represented matrix.
//! The closed-form, triangular and direct Gershgorin implementations use
//! a unit pivot. The approximate Jacobi implementation uses norm-one
//! eigenvectors; the validated similarity-preconditioned implementation
//! may use another nonsingular normalization.
template<class P, class PR>
struct EigenDecompositionData {
    using NumericType = FloatType<P, PR>;
    Vector<NumericType> eigenvalues; //!< Real eigenvalue approximations or enclosures.
    Matrix<NumericType> eigenvectors; //!< Corresponding eigenvectors as columns.

    //! \brief Construct the diagonal eigenvalue matrix D in A*V=V*D.
    //! \details D has exactly zero off-diagonal entries; in ValidatedTag
    //! its diagonal entries are eigenvalue enclosures, not necessarily
    //! point eigenvalues. For every represented input matrix, a coherent
    //! choice of eigenpairs exists in the returned enclosures, but this
    //! does not imply literal equality of independently evaluated
    //! interval matrix products.
    DiagonalMatrix<NumericType> diagonal_matrix() const {
        return DiagonalMatrix<NumericType>(eigenvalues);
    }

    //! \brief Attempt to compute the inverse of the eigenvector matrix V.
    //! \details On success, this is the third factor of A=V*D*V^{-1}.
    //! With validated arithmetic it encloses inverses of *every* matrix
    //! represented by the independent interval entries of V, a stronger
    //! requirement than the pointwise existence of a coherent eigenbasis
    //! established by a successful eigendecomposition. The interval
    //! elimination can therefore be inconclusive even when the spectral
    //! decomposition was certified successfully. No inverse is cached.
    //! \return A computed inverse, or no value when elimination cannot
    //! establish a nonzero pivot. In ApproximateTag it is not certified.
    std::optional<Matrix<NumericType>> inverse_eigenvector_matrix() const {
        try {
            return inverse(eigenvectors);
        } catch (const SingularMatrixException&) {
            return std::nullopt;
        }
    }
};

//! \brief Result of an eigendecomposition attempt.
//! \details Data are present only when status is SUCCESS. INCONCLUSIVE never
//! asserts nonexistence of a real eigenbasis.
template<class P, class PR>
struct EigenDecompositionResult {
    EigenDecompositionStatus status; //!< Success, inconclusive, or unsupported dimension.
    std::optional<EigenDecompositionData<P, PR>> data; //!< Present on successful decomposition.
};

//! \brief Compute a real eigendecomposition in the requested paradigm and precision.
//! \details Real diagonal matrices of arbitrary positive dimension support
//! repeated real eigenvalues, including exact repetitions in ValidatedTag.
//! A validated interval diagonal family needs a fixed ascending order:
//! entries must be strictly separated or identical singleton values.
//! Block-diagonal matrices whose off-block entries are exactly zero are
//! decomposed blockwise, including interleaved blocks and repeated
//! eigenvalues in distinct invariant blocks. ValidatedTag requires a
//! consistent ordering of the resulting eigenvalue enclosures.
//! Non-diagonal triangular matrices allow repeated real eigenvalues when
//! their normalized eigenvector equations admit a triangular unit-pivot
//! basis. A defective or interval-undecidable case is INCONCLUSIVE.
//! Dense 2x2 real matrices with simple spectrum are supported by closed form.
//! For dense higher-order matrices, ApproximateTag uses Jacobi rotations
//! for symmetric inputs (including repeated real eigenvalues) and shifted
//! real QR with Newton refinement for nonsymmetric inputs having a
//! resolvable simple real spectrum.
//! ValidatedTag supports matrices with
//! rigorously isolated Gershgorin disks, possibly after a certified
//! similarity using a point-matrix basis seeded by approximate Jacobi
//! rotations, or through a Krawczyk certificate on individually normalized
//! eigenpair equations. This also applies to some nonsymmetric matrices.
//! These are sufficient conditions, not a general-purpose spectral solver.
//! An approximate result is never claimed to be validated.
//! Interval input denotes a family of matrices; a successful validated
//! result encloses eigenpairs of every member, with consistent ascending
//! eigenvalue order. This is not yet a general-purpose eigensolver.
template<class P, class PR>
class EigenDecomposition {
  public:
    using NumericType = FloatType<P, PR>;
    using ResultType = EigenDecompositionResult<P, PR>;

    //! \brief Construct a decomposer using the requested numeric precision.
    explicit EigenDecomposition(PR precision) : _precision(precision) {}

    //! \brief Attempt real diagonalization of the square matrix M.
    //! \details A non-square matrix violates the precondition. Orders
    //! greater than two may be triangular or admit a dense real eigensolver;
    //! ValidatedTag requires a certified real eigenbasis.
    ResultType operator()(const Matrix<NumericType>& M) const;

    //! \brief Run the same eigendecomposition, collecting optional timings.
    //! \details Diagnostics are reset on entry; passing nullptr is equivalent
    //! to the ordinary overload, with no clock reads.
    ResultType operator()(const Matrix<NumericType>& M,
                          EigenDecompositionDiagnostics* diagnostics) const;

  private:
    //! \brief Decide strict positivity using the requested paradigm.
    bool _positive(const NumericType& x) const;
    //! \brief Decide strict negativity using the requested paradigm.
    bool _negative(const NumericType& x) const;
    //! \brief Test whether a number is represented by exact zero.
    bool _zero(const NumericType& x) const;
    //! \brief Certify or numerically test a nonzero number.
    bool _nonzero(const NumericType& x) const;

    //! \brief Approximate a real orthogonal eigenbasis of a symmetric dense matrix.
    //! \details Implemented using Jacobi plane rotations for ApproximateTag.
    //! Other paradigms or inconclusive numerical convergence report failure.
    ResultType _symmetric(const Matrix<NumericType>& A) const;

    //! \brief Approximate a complete simple real spectrum of a dense matrix.
    //! \details Uses shifted real Householder QR followed by normalized
    //! Newton eigenpair refinement; does not enclose or certify eigenpairs.
    //! Non-real, multiple, unresolved or ill-conditioned spectra return
    //! INCONCLUSIVE rather than an incorrect real diagonalization.
    ResultType _nonsymmetric(const Matrix<NumericType>& A) const;

    //! \brief Certify a dense matrix with pairwise disjoint Gershgorin disks.
    //! \details Encloses a real eigenpair for every disk and every matrix in
    //! the interval input; succeeds only after interval linear solves and
    //! strict column-diagonal-dominance verification of the eigenvector basis.
    ResultType _gershgorin(const Matrix<NumericType>& A) const;

    //! \brief Certify dense eigenpairs after an interval-verified similarity.
    //! \details Uses an approximate Jacobi basis as a *candidate* fixed
    //! coordinate change, verifies its inverse by interval elimination,
    //! then invokes the rigorous Gershgorin eigenpair certificate.
    ResultType _preconditioned(const Matrix<NumericType>& A,
                               EigenDecompositionDiagnostics* diagnostics) const;

    //! \brief Certify each real simple eigenpair via a Krawczyk root box.
    //! \details Fixes one eigenvector component to one, performs approximate
    //! Newton refinement only for seeding, preferring a QR eigenbasis of
    //! the actual midpoint and falling back to symmetrized Jacobi, then
    //! verifies pointwise-real eigenpairs using the shared interval
    //! Krawczyk root predicate.
    ResultType _krawczyk(const Matrix<NumericType>& A,
                          EigenDecompositionDiagnostics* diagnostics) const;

    PR _precision;
};

} // namespace Ariadne

#include "algebra/eigen_decomposition.tpl.hpp"

#endif // ARIADNE_EIGEN_DECOMPOSITION_HPP
