/***************************************************************************
 *            algebra/vector_interface.hpp
 *
 *  Copyright  2010-20  Pieter Collins
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
 *  along with Ariadne Algebra.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file algebra/vector_interface.hpp
 *  @brief Abstract interfaces for vector algebra.
 */

#ifndef ARIADNE_VECTOR_INTERFACE_HPP
#define ARIADNE_VECTOR_INTERFACE_HPP

#include "utility/writable.hpp"

namespace Ariadne {

template<class X, class I = Void>
class VectorInterface;

template<class X, class I>
class VectorInterface {
    typedef I ScalarInterface;
    typedef X ScalarType;
    typedef typename ScalarType::NumericType NumericType;

  public:
    virtual ScalarInterface* _get(SizeType i) const = 0;
    virtual Void _set(SizeType i, const ScalarType& x) const = 0;
    virtual VectorInterface* _add(const VectorInterface*) const = 0;
    virtual VectorInterface* _sub(const VectorInterface*) const = 0;
    virtual VectorInterface* _mul(const ScalarType&) const = 0;
    virtual VectorInterface* _mul(const NumericType&) const = 0;
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_INTERFACE_HPP */
