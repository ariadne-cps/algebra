/***************************************************************************
 *            algebra/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/declarations.hpp
 *  @brief Forward declarations and common algebra type aliases.
 */
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
