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
#include <type_traits>
#include <utility>
#include <variant>

#include "algebra/complex_eigen_decomposition.hpp"
#include "algebra/diagonal_matrix.hpp"
#include "algebra/matrix.hpp"
#include "algebra/schur_decomposition.hpp"
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
    BLOCK_TRIANGULAR, //!< Eigenvectors lifted from coupled triangular blocks.
    JACOBI,        //!< Approximate symmetric Jacobi rotations.
    SHIFTED_QR,    //!< Approximate nonsymmetric real Schur QR and Newton refinement.
    COMPLEX_SCHUR, //!< Approximate complex eigenvectors from real Schur blocks.
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

namespace detail {

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
//! Non-diagonal triangular matrices, also after a simultaneous row/column
//! permutation detected by their acyclic structural nonzero pattern, are
//! solved by triangular substitution. More general structurally block-upper-
//! triangular matrices are decomposed into strongly connected diagonal
//! blocks. Certified interval inverses lift their eigenvectors into the
//! preceding blocks when spectral separation permits it. When a shifted
//! block is singular, spectral-coordinate lifting also permits compatible
//! exactly repeated eigenvalues, provided the coupling has provably zero
//! projection onto each corresponding eigenspace. Exact repeated eigenvalues are
//! accepted when their normalized equations are compatible; a defective
//! or interval-undecidable case is INCONCLUSIVE.
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
//! Nonfinite input coefficients (including NaN and infinity) are
//! inconclusive for all positive orders. SUCCESS additionally requires
//! finite eigenvalues and eigenvector coefficients after every numerical
//! operation, including triangular substitution and block lifting.
//! Interval input denotes a family of matrices; a successful validated
//! result encloses eigenpairs of every member, with consistent ascending
//! eigenvalue order. This is not yet a general-purpose eigensolver.
template<class P, class PR>
class RealEigenDecomposition {
  public:
    using NumericType = FloatType<P, PR>;
    using ResultType = EigenDecompositionResult<P, PR>;

    //! \brief Construct a decomposer using the requested numeric precision.
    explicit RealEigenDecomposition(PR precision) : _precision(precision) {}

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
    //! \brief Execute decomposition strategies before final output checks.
    //! \details The public operator rejects nonfinite results produced by
    //! intermediate numerical overflow, regardless of chosen strategy.
    ResultType _attempt(const Matrix<NumericType>& M,
                        EigenDecompositionDiagnostics* diagnostics) const;

    //! \brief Decide strict positivity using the requested paradigm.
    bool _positive(const NumericType& x) const;
    //! \brief Decide strict negativity using the requested paradigm.
    bool _negative(const NumericType& x) const;
    //! \brief Test whether a number is represented by exact zero.
    bool _zero(const NumericType& x) const;
    //! \brief Certify or numerically test a nonzero number.
    bool _nonzero(const NumericType& x) const;

    //! \brief Approximate a real orthogonal eigenbasis of a symmetric dense matrix.
    //! \details Uses the classical Jacobi eigenvalue method: successive
    //! orthogonal rotations in coordinate planes (p,q) annihilate
    //! off-diagonal entries of a real symmetric matrix. The same rotations
    //! are accumulated into an approximate orthogonal eigenvector basis.
    //! Symmetry is checked; no non-symmetric input is silently symmetrized.
    //! See J. H. Wilkinson, The Algebraic Eigenvalue Problem (1965),
    //! Chapter 5 (symmetric/Hermitian eigenproblems).
    //! Other paradigms or inconclusive numerical convergence report failure.
    ResultType _symmetric(const Matrix<NumericType>& A) const;

    //! \brief Approximate real diagonalization of a dense nonsymmetric matrix.
    //! \details Reuses SchurDecomposition as the numerical eigenvalue
    //! backend: Householder Hessenberg reduction, implicit Francis real
    //! double shifts and real Wilkinson shifts produce Q,T with
    //! \f$A Q \approx Q T\f$. Only a purely 1x1 real Schur block
    //! partition is admissible for this real diagonalization path.
    //! Back-substitution on T produces eigenvector seeds mapped through
    //! Q; Newton refinement then solves the ORIGINAL eigenpair equations
    //! \f$(A-\lambda I)v=0\f$ with a fixed unit-pivot normalization.
    //! An independent invertibility and residual check is required
    //! before SUCCESS. A complex Schur 2x2 block means the matrix
    //! cannot have a complete real eigenbasis; this path returns
    //! INCONCLUSIVE (the public status does not assert impossibility).
    //! Repeated or poorly separated real eigenvalues may also be
    //! inconclusive despite diagonalizability. This is only numerical
    //! ApproximateTag; see SchurDecomposition for the historical and
    //! mathematical sources of Householder, Hessenberg, Francis and
    //! Wilkinson. ValidatedTag retains separate Gershgorin/Krawczyk
    //! certificates and never treats a numerical Schur form as proof.
    ResultType _nonsymmetric(const Matrix<NumericType>& A) const;

    //! \brief Certify a dense matrix using disjoint Gershgorin disks.
    //! \details Gershgorin's circle theorem (1931) localizes eigenvalues
    //! of A within disks centered at diagonal entries, with radii given
    //! by sums of absolute off-diagonal entries. A disjoint disk contains
    //! exactly one eigenvalue with algebraic multiplicity. Here outward
    //! interval arithmetic, certified linear solves, and strict
    //! column-diagonal-dominance establish a real eigenbasis for all
    //! matrices in the input interval family when successful.
    //! See S. Gerschgorin, "Ueber die Abgrenzung der Eigenwerte einer
    //! Matrix", Bulletin de l'Academie des Sciences de l'URSS,
    //! no. 6, 749-754 (1931), https://www.mathnet.ru/eng/im5235.
    ResultType _gershgorin(const Matrix<NumericType>& A) const;

    //! \brief Certify dense eigenpairs after an interval-verified similarity.
    //! \details Tries an approximate basis computed from the actual
    //! midpoint (QR for nonsymmetric matrices, Jacobi for symmetric ones),
    //! then falls back to symmetrized Jacobi when necessary. The fixed
    //! candidate basis is inverted by verified interval elimination;
    //! Gershgorin on the rigorously transformed original matrix
    //! determines whether the eigenpairs can be certified.
    ResultType _preconditioned(const Matrix<NumericType>& A,
                               EigenDecompositionDiagnostics* diagnostics) const;

    //! \brief Certify a real simple eigenpair with an interval Krawczyk operator.
    //! \details Treats the normalized eigenpair equations
    //! \f$F(\lambda,v)=(Av-\lambda v,\ v_p-1)=0\f$ as a square
    //! nonlinear system after fixing a pivot component. Approximate
    //! QR/Newton (or Jacobi) supplies only a candidate center; an
    //! interval Jacobian and approximate inverse enter Krawczyk's
    //! inclusion test. Only a verified inclusion yields SUCCESS.
    //! The interval fixed-point method originates with R. Krawczyk,
    //! "Newton-Algorithmen zur Bestimmung von Nullstellen mit
    //! Fehlerschranken", Computing 4, 187-201 (1969).
    //! Uncertified candidates remain INCONCLUSIVE.
    ResultType _krawczyk(const Matrix<NumericType>& A,
                          EigenDecompositionDiagnostics* diagnostics) const;

    PR _precision;
};

//! \brief End of the specialized real eigenbasis implementation.
//! \details Its proven real and interval methods remain internal.
} // namespace detail

//! \brief A full real or complex diagonalization of a real matrix.
//! \details On SUCCESS, the data holds either a real eigenvector basis
//! (and real eigenvalues) or a complex eigenvector basis. The caller
//! does not choose the scalar field in advance. A real basis is always
//! preferred when it can be obtained; it is never eagerly converted
//! into a complex matrix. ValidatedTag presently returns only a real
//! certificate, because complex interval eigenbasis certification is
//! not implemented. On failure, data is absent.
template<class P, class PR>
struct EigenDecompositionResult {
    using RealData = detail::EigenDecompositionData<P, PR>;
    using ComplexData = detail::ComplexEigenDecompositionData<P, PR>;
    using Data = std::variant<RealData, ComplexData>;

    EigenDecompositionStatus status;
    std::optional<Data> data;

    //! \brief Real factors, or nullptr if the result is complex/absent.
    const RealData* real_data() const {
        return data ? std::get_if<RealData>(&*data) : nullptr;
    }

    //! \brief Complex factors, or nullptr if the result is real/absent.
    const ComplexData* complex_data() const {
        return data ? std::get_if<ComplexData>(&*data) : nullptr;
    }
};

//! \brief Compute a full eigenbasis without requiring the spectral field.
//!
//! \details For a real input A, first attempt the specialized REAL
//! eigenbasis algorithms: diagonal and block structure, symmetric Jacobi,
//! Schur-backed nonsymmetric QR/Newton, and, for ValidatedTag, certified
//! Gershgorin, similarity and Krawczyk. A real success is returned
//! immediately, WITHOUT complex matrix allocations or complex arithmetic.
//! The result representation is a std::variant of real and complex data.
//!
//! \par Complex fallback
//! Only when approximate real diagonalization is inconclusive, a real
//! Schur form is inspected to establish whether it contains a conjugate
//! nonreal 2x2 block. If so, the internal complex eigenvector routine
//! back-substitutes through the real Schur blocks and verifies the
//! original-matrix complex eigenpair residuals and inverse basis.
//! Matrices with repeated real eigenvalues or an unresolved Schur
//! decomposition do not trigger complex matrix allocation. For the
//! uncommon nonreal case, Schur is currently recomputed during the
//! fallback: this is an optimization opportunity but imposes no
//! overhead on successful real decompositions.
//!
//! \par Limitations and mathematical sources
//! EigenDecomposition requires a complete, numerically acceptable
//! eigenvector basis. A defective matrix can have a complete Spectrum
//! and SchurDecomposition but no diagonalization. Inconclusive does
//! not assert mathematical non-diagonalizability; clusters can defeat
//! the sufficient numerical tests. Validated complex eigenbasis
//! certification remains unavailable. See SchurDecomposition for
//! the Householder/Hessenberg/Francis/Wilkinson derivations, and the
//! internal real method for Jacobi, Gershgorin and Krawczyk.
//!
//! \par Performance contract
//! A common REAL success does not call the complex solver, does not
//! allocate ComplexMatrix, and does not perform complex arithmetic.
//! Benchmark DP/MP real-spectrum performance before campaign closure.
template<class P, class PR>
class EigenDecomposition {
  public:
    using NumericType = FloatType<P, PR>;
    using ResultType = EigenDecompositionResult<P, PR>;
    using RealData = typename ResultType::RealData;
    using ComplexData = typename ResultType::ComplexData;
    using Data = typename ResultType::Data;

    explicit EigenDecomposition(PR precision) : _precision(precision) {}

    //! \brief Attempt diagonalization, preferring a real basis.
    ResultType operator()(const Matrix<NumericType>& A) const {
        return (*this)(A, nullptr);
    }

    //! \brief Collect optional real/complex solver diagnostics.
    //! \details A successful real path incurs no complex-search timing
    //! overhead. COMPLEX_SCHUR identifies a complex fallback attempt.
    ResultType operator()(const Matrix<NumericType>& A,
                          EigenDecompositionDiagnostics* diagnostics) const {
        auto real = detail::RealEigenDecomposition<P, PR>(_precision)(
            A, diagnostics);
        if (real.status == EigenDecompositionStatus::SUCCESS) {
            if (!real.data) {
                return {EigenDecompositionStatus::INCONCLUSIVE, std::nullopt};
            }
            return {EigenDecompositionStatus::SUCCESS,
                    Data(std::in_place_type<RealData>, std::move(*real.data))};
        }
        if constexpr (std::is_same_v<P, ApproximateTag>) {
            if (real.status == EigenDecompositionStatus::INCONCLUSIVE) {
                using Clock = std::chrono::steady_clock;
                const auto start = diagnostics
                    ? Clock::now() : Clock::time_point{};
                // Real Schur detection is deliberately performed
                // BEFORE complex allocations or inversion. Merely
                // failing to find a real basis does not justify
                // promoting all coefficients to complex arithmetic.
                const auto schur = SchurDecomposition<P, PR>(_precision)(A);
                bool has_nonreal_block = false;
                if (schur.status == SchurStatus::SUCCESS && schur.data) {
                    for (const SizeType width: schur.data->block_sizes) {
                        if (width == 2u) {
                            has_nonreal_block = true;
                            break;
                        }
                    }
                }
                if (has_nonreal_block) {
                    if (diagnostics) {
                        diagnostics->strategy =
                            EigenDecompositionStrategy::COMPLEX_SCHUR;
                    }
                    auto complex = detail::ComplexEigenDecomposition<P, PR>(
                        _precision)(A);
                    if (diagnostics) {
                        diagnostics->candidate_time +=
                            std::chrono::duration_cast<std::chrono::nanoseconds>(
                                Clock::now() - start);
                    }
                    if (complex.status ==
                        detail::ComplexEigenDecompositionStatus::SUCCESS &&
                        complex.data) {
                        return {EigenDecompositionStatus::SUCCESS,
                                Data(std::in_place_type<ComplexData>,
                                     std::move(*complex.data))};
                    }
                } else if (diagnostics) {
                    diagnostics->candidate_time +=
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            Clock::now() - start);
                }
            }
        }
        return {real.status, std::nullopt};
    }

  private:
    PR _precision;
};

} // namespace Ariadne

#include "algebra/eigen_decomposition.tpl.hpp"

#endif // ARIADNE_EIGEN_DECOMPOSITION_HPP
