/***************************************************************************
 *            algebra/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
 *
 ****************************************************************************/

#ifndef ARIADNE_ALGEBRA_DECLARATIONS_HPP
#define ARIADNE_ALGEBRA_DECLARATIONS_HPP

#include <iosfwd>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

#include "foundation/paradigm.hpp"
#include "foundation/logical.decl.hpp"

#include "algebra/linear_algebra.decl.hpp"
#include "algebra/differential.decl.hpp"

namespace Ariadne {

class UniIndex;
class MultiIndex;


template<class X> class Algebra;
template<class X> class ElementaryAlgebra;
template<class X> class Series;

template<class I, class X> class Polynomial;
template<class X> using UnivariatePolynomial = Polynomial<UniIndex,X>;
template<class X> using MultivariatePolynomial = Polynomial<MultiIndex,X>;

template<class X> class UnivariateChebyshevPolynomial;
template<class X> class MultivariateChebyshevPolynomial;

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_DECLARATIONS_HPP */
