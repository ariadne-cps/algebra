/***************************************************************************
 *            test_sweeper.cpp
 *
 *  Copyright  2026  Pieter Collins
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

#include <cmath>

#include "numeric/numeric.hpp"
#include "interval/interval.hpp"
#include "algebra/expansion.hpp"
#include "algebra/sweeper.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

namespace {

Void test_threshold_sweeper()
{
    MultiplePrecision multiple_precision(128);

    ThresholdSweeper<FloatDP> sweeper_dp(dp,1e-8);
    ThresholdSweeper<FloatMP> sweeper_mp(multiple_precision,std::pow(2.0,-64));

    ARIADNE_TEST_PRINT(sweeper_dp.precision());
    ARIADNE_TEST_PRINT(sweeper_mp.precision());
    ARIADNE_TEST_PRINT(sweeper_dp);
    ARIADNE_TEST_PRINT(sweeper_mp);

    ARIADNE_TEST_FAIL((void)ThresholdSweeper<FloatDP>(dp,-1e-8));

    MultiIndex constant_index({0u});
    ARIADNE_TEST_ASSERT(sweeper_dp.discard(constant_index,FloatDP(0u,dp)));
    ARIADNE_TEST_ASSERT(!sweeper_dp.discard(constant_index,FloatDP(1u,dp)));

    SweeperInterface<FloatDP>* cloned_interface=sweeper_dp._copy();
    ARIADNE_TEST_PRINT(cloned_interface->precision());
    delete cloned_interface;
}

Void test_sweeper_handle()
{
    Sweeper<FloatDP> default_dp_sweeper;
    Sweeper<FloatMP> default_mp_sweeper;
    ARIADNE_TEST_PRINT(default_dp_sweeper.precision());
    ARIADNE_TEST_PRINT(default_mp_sweeper.precision());

    ThresholdSweeper<FloatDP> threshold_sweeper(dp,1e-8);
    Sweeper<FloatDP> threshold_base(threshold_sweeper);
    Sweeper<FloatDP> threshold_copy(threshold_base);
    ARIADNE_TEST_PRINT(threshold_copy.precision());

    Expansion<MultiIndex,FloatDP> sweep_terms(1u,dp);
    sweep_terms.append(MultiIndex({0u}),FloatDP(0u,dp));
    sweep_terms.append(MultiIndex({1u}),FloatDP(1u,dp));
    FloatDPError sweep_error(0u,dp);

    ARIADNE_TEST_EXECUTE(threshold_copy.sweep(sweep_terms,sweep_error));
    ARIADNE_TEST_EQUALS(sweep_terms.size(),1u);
    ARIADNE_TEST_EQUALS(sweep_terms.front().coefficient(),FloatDP(1u,dp));
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_threshold_sweeper());
    ARIADNE_TEST_CALL(test_sweeper_handle());

    return ARIADNE_TEST_FAILURES;
}
