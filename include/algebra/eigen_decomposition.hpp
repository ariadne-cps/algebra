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

#include <optional>

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
//! \details Real triangular matrices of arbitrary positive dimension with
//! pairwise strictly separated diagonal entries are supported. Dense matrices
//! of order one or two are supported when the spectrum is real and simple.
//! For dense higher-order matrices, ApproximateTag supports exactly symmetric
//! real matrices by Jacobi rotations; ValidatedTag supports matrices with
//! rigorously isolated Gershgorin disks, possibly after a certified
//! similarity using a point-matrix basis seeded by approximate Jacobi
//! rotations. This also applies to some nonsymmetric matrices. These
//! are sufficient conditions, not a general-purpose spectral certificate.
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
    //! greater than two may be triangular, symmetric in ApproximateTag,
//! or have separated Gershgorin disks in ValidatedTag.
    ResultType operator()(const Matrix<NumericType>& M) const;

  private:
    //! \brief Decide strict positivity using the requested paradigm.
    bool _positive(const NumericType& x) const;
    //! \brief Decide strict negativity using the requested paradigm.
    bool _negative(const NumericType& x) const;
    //! \brief Test whether a number is represented by exact zero.
    bool _zero(const NumericType& x) const;
    //! \brief Certify or numerically test a nonzero number.
    bool _nonzero(const NumericType& x) const;

    //! \brief Approximate real simple eigenpairs of a symmetric dense matrix.
    //! \details Implemented using Jacobi plane rotations for ApproximateTag.
    //! Other paradigms or inconclusive numerical convergence report failure.
    ResultType _symmetric(const Matrix<NumericType>& A) const;

    //! \brief Certify a dense matrix with pairwise disjoint Gershgorin disks.
    //! \details Encloses a real eigenpair for every disk and every matrix in
    //! the interval input; succeeds only after interval linear solves and
    //! strict column-diagonal-dominance verification of the eigenvector basis.
    ResultType _gershgorin(const Matrix<NumericType>& A) const;

    //! \brief Certify dense eigenpairs after an interval-verified similarity.
    //! \details Uses an approximate Jacobi basis as a *candidate* fixed
    //! coordinate change, verifies its inverse by interval elimination,
    //! then invokes the rigorous Gershgorin eigenpair certificate.
    ResultType _preconditioned(const Matrix<NumericType>& A) const;

    PR _precision;
};

} // namespace Ariadne

#include "algebra/eigen_decomposition.tpl.hpp"

#endif // ARIADNE_EIGEN_DECOMPOSITION_HPP
