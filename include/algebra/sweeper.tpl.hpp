/***************************************************************************
 *            algebra/sweeper.tpl.hpp
 *
 *  Copyright  2010-20  Pieter Collins
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

#include "algebra/sweeper.hpp"

namespace Ariadne {

template<class F>
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,Float<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,Float<PR>>::Iterator curr=p.begin();

    F::set_rounding_upward();
    FloatError<PR> te(0.0_x,e.precision());
    while(adv!=end) {
        if(this->_discard(adv->index(),adv->coefficient())) {
            te+=cast_positive(abs(adv->coefficient()));
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,FloatBounds<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatBounds<PR>>::Iterator curr=p.begin();

    F::set_rounding_upward();
    FloatError<PR> te(0.0_x,e.precision());
    while(adv!=end) {
        if(this->_discard(adv->index(),mag(adv->coefficient()).raw())) {
            te+=mag(adv->coefficient());
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,FloatApproximation<PR>>& p) const
{
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::Iterator curr=p.begin();
    while(adv!=end) {
        if(this->_discard(adv->index(),adv->coefficient().raw())) {
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,FloatUpperInterval<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::Iterator curr=p.begin();

    F::set_rounding_upward();
    FloatError<PR> te(0.0_x,e.precision());
    while(adv!=end) {
        if(this->_discard(adv->index(),mag(adv->coefficient()).raw())) {
            te+=mag(adv->coefficient());
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

namespace {

template<class I, class X>
decltype(auto) sweeper_radius(Expansion<I,X> const& p)
{
    typedef decltype(mag(declval<X>())) R;
    R r=mag(p.zero_coefficient());
    for (auto term : p) {
        if (term.index().degree() != 0) {
            r += mag(term.coefficient());
        }
    }
    return r;
}

} // namespace

template<class F>
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,Float<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,Float<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=sweeper_radius(p)+e;

    F::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(adv->coefficient(),nrm.raw())) {
            te+=cast_positive(abs(adv->coefficient()));
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,FloatBounds<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatBounds<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=sweeper_radius(p)+e;

    F::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(mag(adv->coefficient()).raw(),nrm.raw())) {
            te+=mag(adv->coefficient());
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,FloatApproximation<PR>>& p) const
{
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::Iterator curr=p.begin();

    FloatApproximation<PR> nrm=sweeper_radius(p);

    while(adv!=end) {
        if(this->_discard(adv->coefficient().raw(),nrm.raw())) {
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,FloatUpperInterval<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=sweeper_radius(p)+e;

    F::set_rounding_upward();
    FloatError<PR> te(0u,e.precision());
    while(adv!=end) {
        if(this->_discard(mag(adv->coefficient()).raw(),nrm.raw())) {
            te+=mag(adv->coefficient());
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

} // namespace Ariadne
