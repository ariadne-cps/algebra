/***************************************************************************
 *            algebra/symmetric_matrix.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file algebra/symmetric_matrix.hpp
 *  \brief Symmetric matrices
 */

#ifndef ARIADNE_SYMMETRIC_MATRIX_HPP
#define ARIADNE_SYMMETRIC_MATRIX_HPP

#include <initializer_list>
#include <optional>

#include "utility/metaprogramming.hpp"
#include "macro/general.hpp"
// matrix.hpp declares midpoint<Rational>, which requires Rational to be complete.
#include "numeric/numeric.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace Ariadne {

template<class T> using InitializerList = std::initializer_list<T>;

/************ Matrix *********************************************************/

template<class X> class Vector;
template<class X> class Covector;

template<class X> class Matrix;

class SingularMatrixException;

class NonSymmetricMatrixException { };



//! \ingroup LinearAlgebraSubModule
//! \brief Symmetric matrices over some type \a X.
template<class X> class SymmetricMatrix
    : public MatrixContainer<Matrix<X>>
{
    X _zero;
    SizeType _s;
    Array<X> _ary;
  public:
    typedef X ScalarType;
    typedef X ValueType;
  public:

    //!@{
    //! \name Constructors

    //! Destructor
    ~SymmetricMatrix();

    //! Construct a symmetric matrix with \a n rows and \a n columns with values initialised to \a z.
    SymmetricMatrix(SizeType n, X const& z);

    //! Construct a symmetric matrix with \a n rows and \a n columns with values initialised to zero.
    template<class... PRS> requires Constructible<X,Nat,PRS...>
    SymmetricMatrix(SizeType n, PRS... prs) : SymmetricMatrix(n,X(0u,prs...)) { }

    //! Construct a matrix using initializer lists.
    SymmetricMatrix(InitializerList<InitializerList<X>> lst);

    //!@{
    //! \name Static constructors

    //! \brief The zero matrix with \a n rows and \a n columns.
    template<class... PRS> requires Constructible<X,Nat,PRS...>
    static SymmetricMatrix<X> zero(SizeType n, PRS... prs);
    //! \brief The itentity matrix with \a n rows and \a n columns.
    template<class... PRS> requires Constructible<X,Nat,PRS...>
    static SymmetricMatrix<X> identity(SizeType n, PRS... prs);
    //!@}

    //! Explicit conversion of a normal matrix to a symmetric matrix.
    SymmetricMatrix(const Matrix<X>& A);

    //! Conversion of a symmetric matrix to a dense matrix.
    operator Matrix<X> () const;

    //! \brief The number of rows (equivalently, columns) of the matrix.
    SizeType size() const;
     //! \brief The number of rows of the matrix.
    SizeType row_size() const;
    //! \brief The number of columns of the matrix.
    SizeType column_size() const;
    //! \brief Resize to an \a m by \a n matrix.
    Void resize(SizeType n);
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    X& at(SizeType i, SizeType j);
    const X& at(SizeType i, SizeType j) const;
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    const X& get(SizeType i, SizeType j) const ;
    //! \brief Set the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column to \a x.
    Void set(SizeType i, SizeType j, const X& x);
    //! \brief A pointer to the first element of the data storage.
    const X* begin() const;

#ifdef DOXYGEN
    //! \brief C-style subscripting operator.
    X& operator[][](SizeType i, SizeType j);
    //! \brief C-style constant subscripting operator.
    const X& operator[][](SizeType i, SizeType j) const;
#else
    MatrixRow<const SymmetricMatrix<X>> operator[](SizeType i) const;
    MatrixRow<SymmetricMatrix<X>> operator[](SizeType i);
#endif
    //! \brief The zero element of the field/algebra of the matrix.
    X zero_element() const;
  public:
    template<class T> friend OutputStream& operator<<(OutputStream& os, SymmetricMatrix<T>const& A);
    template<class T> friend SymmetricMatrix<T> operator+(SymmetricMatrix<T> const&, SymmetricMatrix<T> const&);
    template<class T> friend SymmetricMatrix<T> operator-(SymmetricMatrix<T> const&, SymmetricMatrix<T> const&);
    template<class T> friend SymmetricMatrix<T> operator*(T const&, SymmetricMatrix<T> const&);
    template<class T> friend SymmetricMatrix<T> operator*(SymmetricMatrix<T> const&, T const&);
    template<class T> friend SymmetricMatrix<T> operator/(SymmetricMatrix<T> const&, T const&);
    template<class T> friend SymmetricMatrix<T> symmetrize(Matrix<T> const& A);
    template<class T> friend SymmetricMatrix<T> outer(Matrix<T> const& A);
    template<class T> friend SymmetricMatrix<T> AxTpTA(Matrix<T> const&, SymmetricMatrix<T> const&);
  private:
    OutputStream& _write(OutputStream& os) const;

    friend Vector<X>& to_vector(SymmetricMatrix<X>& S) { return reinterpret_cast<Vector<X>&>(S._ary); }
    friend Vector<X> const& to_vector(SymmetricMatrix<X> const& S) { return reinterpret_cast<Vector<X>const&>(S._ary); }
  private:
    SizeType _fast_position(SizeType i, SizeType j) const { return i*(2*this->size()-i-1)/2+j; }
    SizeType _position(SizeType i, SizeType j) const {
        ARIADNE_PRECONDITION(i<this->size() && j<this->size());
        if(i<=j) { return this->_fast_position(i,j); } else { return this->_fast_position(j,i); }
    }
};

template<class X> Vector<X>& to_vector(SymmetricMatrix<X>& S);
template<class X> Vector<X> const& to_vector(SymmetricMatrix<X> const& S);




template<class X> SymmetricMatrix<X>::~SymmetricMatrix() {
}

template<class X> SymmetricMatrix<X>::SymmetricMatrix(SizeType n, X const& z)
    : _zero(z), _s(n), _ary(n*(n+1)/2,z)
{
}

template<class X> SymmetricMatrix<X>::SymmetricMatrix(Matrix<X> const& A)
    : _zero(A.zero_element()), _s(A.row_size()), _ary(_s*(_s+1)/2,_zero)
{
    ARIADNE_PRECONDITION(A.row_size()==A.column_size());
    const SizeType n=this->_s;
    for(SizeType i=0; i!=n; ++i) {
        this->_ary[this->_fast_position(i,i)]=A[i][i];
        for(SizeType j=i+1; j!=n; ++j) {
            if(decide(A[i][j]!=A[j][i])) { throw NonSymmetricMatrixException(); }
            this->_ary[this->_fast_position(i,j)]=A[i][j];
        }
    }
}

template<class X> template<class... PRS> requires Constructible<X,Nat,PRS...>
SymmetricMatrix<X> SymmetricMatrix<X>::zero(SizeType n, PRS... prs) {
    return SymmetricMatrix<X>(n,X(0u,prs...));
}

template<class X> template<class... PRS> requires Constructible<X,Nat,PRS...>
SymmetricMatrix<X> SymmetricMatrix<X>::identity(SizeType n, PRS... prs) {
    SymmetricMatrix<X> S(n,X(0u,prs...));
    for(SizeType i=0; i!=n; ++i) {
        S[i][i]=1u;
    }
    return S;
}

template<class X> inline SizeType SymmetricMatrix<X>::size() const {
    return this->_s;
}

template<class X> inline SizeType SymmetricMatrix<X>::row_size() const {
    return this->_s;
}

template<class X> inline SizeType SymmetricMatrix<X>::column_size() const {
    return this->_s;
}

template<class X> inline Void SymmetricMatrix<X>::resize(SizeType n) {
    this->_s = n;
    this->_ary.resize(n*(n+1)/2,this->_zero);
}

template<class X> inline X& SymmetricMatrix<X>::at(SizeType i, SizeType j) {
    return this->_ary[this->_position(i,j)];
}

template<class X> inline const X& SymmetricMatrix<X>::at(SizeType i, SizeType j) const {
    return this->_ary[this->_position(i,j)];
}

template<class X> inline const X& SymmetricMatrix<X>::get(SizeType i, SizeType j) const {
    return this->_ary[this->_position(i,j)];
}

template<class X> inline Void SymmetricMatrix<X>::set(SizeType i, SizeType j, const X& x) {
    this->_ary[this->_position(i,j)]=x;
}

template<class X> inline const X* SymmetricMatrix<X>::begin() const {
    return this->_ary.begin();
}

#ifndef DOXYGEN
template<class X> inline MatrixRow<const SymmetricMatrix<X>> SymmetricMatrix<X>::operator[](SizeType i) const {
    return MatrixRow<const SymmetricMatrix<X>>(*this,i);
}
template<class X> inline MatrixRow<SymmetricMatrix<X>> SymmetricMatrix<X>::operator[](SizeType i) {
    return MatrixRow<SymmetricMatrix<X>>(*this,i);
}
#endif

template<class X> inline X SymmetricMatrix<X>::zero_element() const {
    return _zero;
}

template<class X> inline OutputStream& operator<<(OutputStream& os, SymmetricMatrix<X>const& S) {
    return S._write(os);
}



template<class X> SymmetricMatrix<X>::operator Matrix<X>() const {
    SymmetricMatrix<X> const& S=*this;
    Matrix<X> A(S.row_size(),S.column_size(),S.zero_element());
    for(SizeType i=0; i!=S.row_size(); ++i) {
        A[i][i]=S.at(i,i);
        for(SizeType j=i+1; j!=S.row_size(); ++j) {
            A[i][j]=S.at(i,j);
            A[j][i]=S.at(i,j);
        }
    }
    return A;
}


template<class X> OutputStream& SymmetricMatrix<X>::_write(OutputStream& os) const {
    return os << Matrix<X>(*this);
}

template<class X> SymmetricMatrix<X> operator+(SymmetricMatrix<X> const& S1, SymmetricMatrix<X> const& S2) {
    ARIADNE_PRECONDITION(S1.size()==S2.size());
    const SizeType n=S1.size();
    SymmetricMatrix<X> R(n,S1.zero_element()+S2.zero_element());
    for(SizeType i=0; i!=n*(n+1)/2; ++i) {
        R._ary[i]=S1._ary[i]+S2._ary[i];
    }
    return R;
}

template<class X> SymmetricMatrix<X> operator-(SymmetricMatrix<X> const& S1, SymmetricMatrix<X> const& S2) {
    assert(S1.size()==S2.size());
    const SizeType n=S1.size();
    SymmetricMatrix<X> R(n,S1.zero_element()-S2.zero_element());
    for(SizeType i=0; i!=n*(n+1)/2; ++i) {
        R._ary[i]=S1._ary[i]-S2._ary[i];
    }
    return R;
}

template<class X> SymmetricMatrix<X> operator*(X const& x1, SymmetricMatrix<X> const& S2) {
    const SizeType n=S2.size();
    SymmetricMatrix<X> R(n,x1*S2.zero_element());
    for(SizeType i=0; i!=n*(n+1)/2; ++i) {
        R._ary[i]=x1*S2._ary[i];
    }
    return R;
}

template<class X> SymmetricMatrix<X> operator*(SymmetricMatrix<X> const& S1, X const& x2) {
    const SizeType n=S1.size();
    SymmetricMatrix<X> R(n,S1.zero_element()*x2);
    for(SizeType i=0; i!=n*(n+1)/2; ++i) {
        R._ary[i]=S1._ary[i]*x2;
    }
    return R;
}

template<class X> SymmetricMatrix<X> operator/(SymmetricMatrix<X> const& S1, X const& x2) {
    const SizeType n=S1.size();
    SymmetricMatrix<X> R(n,S1.zero_element()/x2);
    for(SizeType i=0; i!=n*(n+1)/2; ++i) {
        R._ary[i]=S1._ary[i]/x2;
    }
    return R;
}

template<class X> SymmetricMatrix<X> outer(Matrix<X> const& A) {
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    SymmetricMatrix<X> S(n,A.zero_element());
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=i; j!=n; ++j) {
            X& Sij=S.at(i,j);
            Sij=0;
            for(SizeType k=0; k!=m; ++k) {
                Sij+=A[k][i]*A[k][j];
            }
        }
    }
    return S;
}

template<class X> SymmetricMatrix<X> outer(SymmetricMatrix<X> const& A) {
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    SymmetricMatrix<X> S(n,A.zero_element());
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=i; j!=n; ++j) {
            X& Sij=S.at(i,j);
            Sij=0;
            for(SizeType k=0; k!=m; ++k) {
                Sij+=A[k][i]*A[k][j];
            }
        }
    }
    return S;
}

template<class X> SymmetricMatrix<X> symmetrize(Matrix<X> const& A) {
    ARIADNE_PRECONDITION(A.row_size()==A.column_size());
    const SizeType n=A.row_size();
    SymmetricMatrix S(n,A.zero_element());
    for(SizeType i=0; i!=n; ++i) {
        S.at(i,i)=A[i][i];
        for(SizeType j=i+1; j!=n; ++j) {
            S.at(i,j)=hlf(A[i][j]+A[j][i]);
        }
    }
    return S;
}



//!@{
//! \name Validated symmetric linear algebra

//! \ingroup LinearAlgebraSubModule
//! \brief Rigorous three-way classification of positive definiteness.
//! \details For an interval symmetric matrix, the classification refers
//! to *every* real symmetric matrix represented by its entries.
//! Indeterminate is not a claim that some represented matrix is singular.
enum class PositiveDefinitenessStatus {
    POSITIVE_DEFINITE,      //!< Every represented real matrix is positive definite.
    NOT_POSITIVE_DEFINITE,  //!< No represented real matrix is positive definite.
    INDETERMINATE           //!< Neither assertion could be certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Certify positive definiteness using an interval LDL^T elimination.
//! \param A An interval symmetric matrix. Off-diagonal entries use the same
//! interval in both symmetric positions; each interval encloses all
//! possible exact coefficient values.
//! \return POSITIVE_DEFINITE if all successive interval Schur pivots
//! are strictly positive; NOT_POSITIVE_DEFINITE if a pivot is provably
//! nonpositive after all preceding pivots were proved positive;
//! otherwise INDETERMINATE.
//! \details Arithmetic on FloatDPBounds uses directed outward rounding.
//! The conclusion is universally quantified over all represented
//! symmetric real matrices. An undecidable pivot is *never* rounded
//! or treated as zero to select an arbitrary numerical rank.
//! The zero-dimensional matrix is positive definite by convention.
//! This is a sufficient interval test: dependency overestimation can
//! lead to INDETERMINATE even when every represented matrix is positive
//! definite. The method does not certify positive semidefiniteness
//! with zero pivots, and it does not compute eigenvalues.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, LDL^T factorization
//! and symmetric Schur complements.
PositiveDefinitenessStatus validated_positive_definiteness(
    const SymmetricMatrix<FloatDPBounds>& A);

//! \ingroup LinearAlgebraSubModule
//! \brief The result of an attempted validated interval Cholesky decomposition.
//! \details The factor is present if and only if the status is
//! POSITIVE_DEFINITE. Otherwise no complete validated factor is available.
//! A successful lower matrix encloses the lower-triangular Cholesky factor
//! of every real symmetric matrix represented by the interval input.
//! In particular, the result is not a pointwise equality between interval
//! matrices: L*transpose(L) is an interval enclosure of the input family.
struct ValidatedCholeskyResult {
    PositiveDefinitenessStatus status; //!< Rigorous universal classification.
    std::optional<Matrix<FloatDPBounds>> lower; //!< Validated lower factor, if certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Compute a lower-triangular Cholesky enclosure with outward rounding.
//! \param A Symmetric matrix with coefficients represented as FloatDPBounds.
//! \return A certified Cholesky enclosure, a proof of failure of positive
//! definiteness for all represented matrices, or INDETERMINATE.
//! \details For each represented real symmetric matrix a, when the status
//! is POSITIVE_DEFINITE, the unique Cholesky factor l (a=l*l^T,
//! with positive diagonal) belongs componentwise to the returned
//! interval matrix L. The upper triangular portion of L is exactly zero.
//! Directed square roots are computed in binary64 using bounded bisection
//! of IEEE-754 floating-point values, with exact Dyadic square comparisons.
//! They require neither MPFR nor a potentially nonterminating Newton
//! iteration; both returned endpoints rigorously enclose the exact root.
//!
//! The algorithm computes successive Schur pivots
//! d_i=a_ii-sum_{k<i} l_ik^2, and l_ij is formed by the corresponding
//! triangular recurrence. Each interval operation encloses its exact
//! real counterpart. A pivot known strictly positive permits sqrt(d_i);
//! a pivot known nonpositive, with all preceding pivots certified positive,
//! proves that no represented matrix is positive definite; otherwise the
//! method returns INDETERMINATE without dividing by an uncertain pivot.
//!
//! No floating-point tolerance is used. Dependency overestimation may
//! yield INDETERMINATE even for a family of positive-definite matrices.
//! An empty matrix is positive definite by convention.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, Cholesky factorization.
ValidatedCholeskyResult validated_cholesky_decomposition(
    const SymmetricMatrix<FloatDPBounds>& A);

//! \ingroup LinearAlgebraSubModule
//! \brief Result of a validated symmetric positive-definite linear solve.
//! \tparam T Either Matrix<FloatDPBounds> or Vector<FloatDPBounds>.
//! \details A solution is present exactly when positive definiteness was
//! certified. It contains every exact solution for every admissible real
//! symmetric coefficient matrix and every admissible right-hand side.
//! An absent solution is not a claim that the linear system is unsolvable:
//! the positive-definiteness precondition could not be established.
template<class T> struct ValidatedPositiveDefiniteSolveResult {
    PositiveDefinitenessStatus status; //!< Classification of the coefficient matrix.
    std::optional<T> solution;         //!< Validated solution enclosure if certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Enclose the solutions of AX=B using validated Cholesky factorisation.
//! \param A Symmetric interval coefficient matrix.
//! \param B Interval right-hand sides, with exactly A.size() rows.
//! \return POSITIVE_DEFINITE together with an enclosure of every solution
//! if Cholesky can be certified; otherwise the failure or indeterminate
//! status of the factorisation, without a solution.
//! \pre B.row_size()==A.size().
//! \details For every real symmetric a in A and every real b in B,
//! the returned X includes the unique solution x of ax=b componentwise.
//! The algorithm first encloses a=L*L^T via validated Cholesky, then
//! applies interval forward and backward substitution. All arithmetic
//! is outward rounded. Correlations between entries are not retained;
//! the returned interval matrix may be wider than the optimal enclosure.
//! No inverse is constructed, and no numerical rank threshold is used.
//! Zero-dimensional matrices and a zero-column right-hand side are valid.
//! A non-SPD status does not imply that AX=B is inconsistent.
//! \par Reference
//! G. H. Golub and C. F. Van Loan, \em Matrix \em Computations,
//! 4th ed., Johns Hopkins University Press, 2013, Cholesky factorization
//! and triangular substitution.
ValidatedPositiveDefiniteSolveResult<Matrix<FloatDPBounds>>
validated_positive_definite_solve(
    const SymmetricMatrix<FloatDPBounds>& A,
    const Matrix<FloatDPBounds>& B);

//! \ingroup LinearAlgebraSubModule
//! \brief Enclose the solutions of Ax=b via validated Cholesky.
//! \param A Symmetric interval coefficient matrix.
//! \param b Interval right-hand side with exactly A.size() elements.
//! \return A certified interval solution iff the input is certified SPD;
//! otherwise the factorisation status and no solution.
//! \pre b.size()==A.size().
//! \details Equivalent to the matrix right-hand-side overload with
//! one column, including its universal enclosure guarantee.
ValidatedPositiveDefiniteSolveResult<Vector<FloatDPBounds>>
validated_positive_definite_solve(
    const SymmetricMatrix<FloatDPBounds>& A,
    const Vector<FloatDPBounds>& b);

//! \ingroup LinearAlgebraSubModule
//! \brief Certification status for an interval saddle-point linear system.
//! \details The block coefficient matrix has form [H A^T; A 0].
//! SOLVED proves that the unique pair (x,y) is enclosed for every
//! real system represented by the interval data. Other statuses
//! identify why this *sufficient* certificate did not succeed; they
//! do not assert that the block system is inconsistent.
enum class ValidatedSaddlePointStatus {
    SOLVED,                              //!< A unique enclosed solution for every admissible system.
    LEADING_BLOCK_NOT_POSITIVE_DEFINITE, //!< H is universally non-SPD and neither certified route succeeded.
    DEPENDENT_COUPLING_ROWS,             //!< Every represented A has dependent rows.
    INDETERMINATE                        //!< Neither solvability nor a structural obstruction certified.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Validated solution vectors of a symmetric saddle-point system.
//! \details Both optional vectors are present exactly when status==SOLVED.
//! \c primary encloses the first block x; \c secondary encloses the
//! second block y. The names convey no optimization-specific semantics.
struct ValidatedSaddlePointResult {
    ValidatedSaddlePointStatus status; //!< Universal certification outcome.
    std::optional<Vector<FloatDPBounds>> primary; //!< n-vector solving the first block.
    std::optional<Vector<FloatDPBounds>> secondary; //!< m-vector solving the second block.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Solve a symmetric block saddle-point system rigorously.
//! \param H Symmetric n-by-n interval leading block.
//! \param A Interval m-by-n coupling block.
//! \param f Interval n-vector, right-hand side of the first block equation.
//! \param g Interval m-vector, right-hand side of the second block equation.
//! \pre A.column_size()==H.size(), f.size()==H.size(),
//! and g.size()==A.row_size().
//! \return SOLVED and certified componentwise enclosures for (x,y) in
//! H*x+A^T*y=f, A*x=g, if the interval SPD/Schur tests prove uniqueness
//! for every admissible real block matrix and right-hand side.
//! If H is certified non-SPD and no alternative succeeds, returns
//! LEADING_BLOCK_NOT_POSITIVE_DEFINITE; if all represented A have
//! dependent rows, returns DEPENDENT_COUPLING_ROWS; otherwise
//! INDETERMINATE. Neither failure status claims inconsistency.
//! \details First attempts a validated Cholesky solve for H, forming
//! S=A*H^{-1}*A^T. If H cannot be certified SPD and m>0, attempts
//! the algebraically equivalent system using the augmented block
//! H+A^T*A and first right-hand side f+A^T*g. This fixed algebraic
//! augmentation retains the same unknowns (x,y) and can certify
//! systems with semidefinite or indefinite leading blocks.
//! Positive definiteness of the chosen leading block and of its
//! Schur complement certifies uniqueness of the block system.
//! All interval calculations round outward. Correlation overestimation
//! may prevent certification, in which case the routine is inconclusive.
//! With m=0 the routine reduces to a validated SPD solve; n=m=0
//! is admissible. No numerical tolerance or explicit inverse is used.
//! This overload delegates to the multiple-right-hand-side version.
//! \par Reference
//! M. Benzi, G. H. Golub and J. Liesen, Numerical solution of saddle
//! point problems, \em Acta \em Numerica 14 (2005), pp. 1--137.
ValidatedSaddlePointResult validated_saddle_point_solve(
    const SymmetricMatrix<FloatDPBounds>& H,
    const Matrix<FloatDPBounds>& A,
    const Vector<FloatDPBounds>& f,
    const Vector<FloatDPBounds>& g);

//! \ingroup LinearAlgebraSubModule
//! \brief Validated enclosures for multiple saddle-point right-hand sides.
//! \details Both matrices are present exactly when status==SOLVED.
//! Their columns enclose the respective solutions of the upper and
//! lower block equations of the same represented system family.
struct ValidatedSaddlePointMatrixResult {
    ValidatedSaddlePointStatus status; //!< Universal certification outcome.
    std::optional<Matrix<FloatDPBounds>> primary; //!< n-by-p first-block solution.
    std::optional<Matrix<FloatDPBounds>> secondary; //!< m-by-p second-block solution.
};

//! \ingroup LinearAlgebraSubModule
//! \brief Solve multiple right-hand sides of a symmetric saddle-point system.
//! \param H Symmetric n-by-n interval leading block.
//! \param A Interval m-by-n coupling block.
//! \param F Interval n-by-p upper-block right-hand sides.
//! \param G Interval m-by-p lower-block right-hand sides.
//! \pre A.column_size()==H.size(), F.row_size()==H.size(),
//! G.row_size()==A.row_size(), F.column_size()==G.column_size().
//! \return On SOLVED, certified n-by-p and m-by-p solution matrices
//! (X,Y) satisfying H*X+A^T*Y=F, A*X=G for each admissible real
//! selection of H,A,F,G; otherwise the same non-success statuses
//! as the vector overload and no solution matrices.
//! \details Uses one validated Cholesky factorization of H (or of
//! H+A^T*A after the equivalent augmented-system transformation)
//! and one of its symmetric Schur complement, shared among all
//! p right-hand sides. The augmented upper right-hand side is
//! F+A^T*G. Positive definiteness of the two factors is sufficient
//! for universal uniqueness; a failed Schur test alone never proves
//! dependent coupling rows. Supports p=0, m=0 and n=m=0.
//! No numerical rank threshold or explicit inverse is introduced.
//! \par Reference
//! M. Benzi, G. H. Golub and J. Liesen, Numerical solution of saddle
//! point problems, \em Acta \em Numerica 14 (2005), pp. 1--137.
ValidatedSaddlePointMatrixResult validated_saddle_point_solve(
    const SymmetricMatrix<FloatDPBounds>& H,
    const Matrix<FloatDPBounds>& A,
    const Matrix<FloatDPBounds>& F,
    const Matrix<FloatDPBounds>& G);

//!@}

} // namespace Ariadne

#endif
