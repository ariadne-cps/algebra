/***************************************************************************
 *            algebra/tensor.hpp
 *
 *  Copyright  2018-20  Pieter Collins
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

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/** @file algebra/tensor.hpp
 *  @brief Tensor data structures and operations.
 */

#ifndef ARIADNE_TENSOR_HPP
#define ARIADNE_TENSOR_HPP

#include "diagonal_matrix.hpp"
#include "matrix.hpp"
#include "vector.hpp"

namespace Ariadne {

//! \brief A row of a rank-2 tensor.
template<class T>
class TensorRow {
    T _t;
    SizeType _i;

  public:
    /** @brief Construct a \c TensorRow from t, i. */
    TensorRow(T& t, SizeType i) : _t(t), _i(i) {};
    /** @brief Return the selected element or component. */
    const T& operator[](SizeType j) const { return _t[Array<SizeType>({_i, j})]; }
    /** @brief Return the selected element or component. */
    T& operator[](SizeType j) { return _t[Array<SizeType>({_i, j})]; }
};

//! \ingroup LinearAlgebraSubModule
//! \brief A rank-\a N tensor with elements of type \a X.
template<SizeType N, class X>
class Tensor {
    //   TODO: Also provide a version which is not templated on rank.
    Array<SizeType> _ns;
    Array<X> _a;

  public:
    typedef X ValueType;
    /** @brief Construct a \c Tensor from ns, z. */
    Tensor(const Array<SizeType>& ns, const X& z) : _ns(ns), _a(_total_size(), z) {}
    /** @brief Construct a \c Tensor from ns, g. */
    template<class G>
        requires InvocableReturning<X, G, Array<SizeType>>
    Tensor(const Array<SizeType>& ns, const G& g);
    /** @brief Return the matrix rank. */
    constexpr SizeType rank() const { return _ns.size(); }
    /** @brief Return the extents of the dimensions. */
    Array<SizeType> sizes() const { return _ns; }
    /** @brief Return the number of elements. */
    SizeType size(SizeType i) const { return _ns[i]; }
    /** @brief Return the selected element or component. */
    const X& operator[](const Array<SizeType>& is) const { return _a[_index(is)]; }
    /** @brief Return the selected element or component. */
    X& operator[](const Array<SizeType>& is) { return _a[_index(is)]; }

    /** @brief Return the selected element or component. */
    TensorRow<const Tensor<N, X>&> operator[](const SizeType& i) const { return TensorRow<const Tensor<N, X>&>(*this, i); }
    /** @brief Return the selected element or component. */
    TensorRow<Tensor<N, X>&> operator[](const SizeType& i) { return TensorRow<Tensor<N, X>&>(*this, i); }

  private:

  public:
    /** @brief Return the selected index. */
    SizeType _index(Array<SizeType> is) const {
        SizeType k = 0;
        SizeType r = is[k];
        ++k;
        while (k != is.size()) {
            r *= _ns[k];
            r += is[k];
            ++k;
        }
        return r;
    }
    /** @brief Return the total number of storage units required. */
    SizeType _total_size() const {
        SizeType r = 1;
        for (SizeType i = 0; i != N; ++i) {
            r *= _ns[i];
        }
        return r;
    }
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
};

/** @brief Apply the transformation to the supplied elements. */
template<class X, class F>
    requires Invocable<F, X>
auto transform(const Tensor<2, X>& t, const F& f) -> Tensor<2, ResultOf<F(X)>> {
    Tensor<2, ResultOf<F(X)>> r(t.sizes(), f(t[{0, 0}]));
    for (SizeType i0 = 0; i0 != r.size(0); ++i0) {
        for (SizeType i1 = 0; i1 != r.size(1); ++i1) {
            r[{i0, i1}] = f(t[{i0, i1}]);
        }
    }
    return r;
}

template<SizeType N, class X>
template<class G>
    requires InvocableReturning<X, G, Array<SizeType>>
Tensor<N, X>::Tensor(const Array<SizeType>& ns, const G& g)
    : _ns(ns), _a(_total_size(), g({0, 0})) {
    static_assert(N == 2);
    for (SizeType i0 = 0; i0 != ns[0]; ++i0) {
        for (SizeType i1 = 0; i1 != ns[1]; ++i1) {
            _a[i0 * ns[1] + i1] = g({i0, i1});
        }
    }
}

} // namespace Ariadne

#endif // ARIADNE_TENSOR_HPP
