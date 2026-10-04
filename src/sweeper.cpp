/***************************************************************************
 *            algebra/sweeper.cpp
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

#include "numeric/numeric.hpp"
#include "interval/interval.hpp"

#include "algebra/sweeper.hpp"
#include "algebra/sweeper.tpl.hpp"

namespace Ariadne {

template class SweeperBase<FloatDP>;
template class SweeperBase<FloatMP>;

template class RelativeSweeperBase<FloatDP>;
template class RelativeSweeperBase<FloatMP>;

} // namespace Ariadne
