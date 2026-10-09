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
    ARIADNE_TEST_FAIL((void)ThresholdSweeper<FloatDP>(dp,FloatDP(-1,dp)));

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

Void test_sweeper_overloads()
{
    ThresholdSweeper<FloatDP> threshold_sweeper(dp,FloatDP(0.25_x,dp));
    SweeperInterface<FloatDP> const& threshold_interface=threshold_sweeper;

    Expansion<MultiIndex,FloatDP> float_terms(1u,dp);
    float_terms.append(MultiIndex({0u}),FloatDP(0.125_x,dp));
    float_terms.append(MultiIndex({1u}),FloatDP(1.0_x,dp));
    FloatDPError float_error(0u,dp);
    ARIADNE_TEST_EXECUTE(threshold_interface.sweep(float_terms,float_error));
    ARIADNE_TEST_EQUALS(float_terms.size(),1u);
    ARIADNE_TEST_EQUALS(float_terms.front().index(),MultiIndex({1u}));

    Expansion<MultiIndex,FloatDPBounds> bounds_terms(1u,dp);
    bounds_terms.append(MultiIndex({0u}),FloatDPBounds(0.125_x,dp));
    bounds_terms.append(MultiIndex({1u}),FloatDPBounds(1.0_x,dp));

    Expansion<MultiIndex,FloatDPApproximation> approximation_terms(bounds_terms);
    Expansion<MultiIndex,FloatDPUpperInterval> upper_terms(bounds_terms);

    FloatDPError bounds_error(0u,dp);
    ARIADNE_TEST_EXECUTE(threshold_interface.sweep(bounds_terms,bounds_error));
    ARIADNE_TEST_EQUALS(bounds_terms.size(),1u);
    ARIADNE_TEST_EQUALS(bounds_terms.front().index(),MultiIndex({1u}));

    ARIADNE_TEST_EXECUTE(threshold_interface.sweep(approximation_terms));
    ARIADNE_TEST_EQUALS(approximation_terms.size(),1u);
    ARIADNE_TEST_EQUALS(approximation_terms.front().index(),MultiIndex({1u}));

    FloatDPError upper_error(0u,dp);
    ARIADNE_TEST_EXECUTE(threshold_interface.sweep(upper_terms,upper_error));
    ARIADNE_TEST_EQUALS(upper_terms.size(),1u);
    ARIADNE_TEST_EQUALS(upper_terms.front().index(),MultiIndex({1u}));
}

Void test_relative_sweeper()
{
    RelativeThresholdSweeper<FloatDP> relative_sweeper(dp,FloatDP(0.25_x,dp));
    ARIADNE_TEST_EQUALS(relative_sweeper.relative_sweep_threshold(),FloatDP(0.25_x,dp));
    ARIADNE_TEST_FAIL((void)RelativeThresholdSweeper<FloatDP>(dp,FloatDP(0u,dp)));
    ARIADNE_TEST_PRINT(relative_sweeper);

    SweeperInterface<FloatDP> const& relative_interface=relative_sweeper;
    ARIADNE_TEST_PRINT(relative_interface.precision());

    SweeperInterface<FloatDP>* relative_copy=relative_sweeper._copy();
    ARIADNE_TEST_PRINT(relative_copy->precision());
    delete relative_copy;

    Expansion<MultiIndex,FloatDP> float_terms(1u,dp);
    float_terms.append(MultiIndex({0u}),FloatDP(0.125_x,dp));
    float_terms.append(MultiIndex({1u}),FloatDP(1.0_x,dp));
    FloatDPError float_error(0u,dp);
    ARIADNE_TEST_EXECUTE(relative_interface.sweep(float_terms,float_error));
    ARIADNE_TEST_EQUALS(float_terms.size(),1u);
    ARIADNE_TEST_EQUALS(float_terms.front().index(),MultiIndex({1u}));

    Expansion<MultiIndex,FloatDPBounds> bounds_source(1u,dp);
    bounds_source.append(MultiIndex({0u}),FloatDPBounds(0.125_x,dp));
    bounds_source.append(MultiIndex({1u}),FloatDPBounds(1.0_x,dp));

    Expansion<MultiIndex,FloatDPBounds> bounds_terms(bounds_source);
    Expansion<MultiIndex,FloatDPApproximation> approximation_terms(bounds_source);
    Expansion<MultiIndex,FloatDPUpperInterval> upper_terms(bounds_source);

    FloatDPError bounds_error(0u,dp);
    ARIADNE_TEST_EXECUTE(relative_interface.sweep(bounds_terms,bounds_error));
    ARIADNE_TEST_EQUALS(bounds_terms.size(),1u);
    ARIADNE_TEST_EQUALS(bounds_terms.front().index(),MultiIndex({1u}));

    ARIADNE_TEST_EXECUTE(relative_interface.sweep(approximation_terms));
    ARIADNE_TEST_EQUALS(approximation_terms.size(),1u);
    ARIADNE_TEST_EQUALS(approximation_terms.front().index(),MultiIndex({1u}));

    FloatDPError upper_error(0u,dp);
    ARIADNE_TEST_EXECUTE(relative_interface.sweep(upper_terms,upper_error));
    ARIADNE_TEST_EQUALS(upper_terms.size(),1u);
    ARIADNE_TEST_EQUALS(upper_terms.front().index(),MultiIndex({1u}));
}

} // namespace

Int main()
{
    ARIADNE_TEST_CALL(test_threshold_sweeper());
    ARIADNE_TEST_CALL(test_sweeper_handle());
    ARIADNE_TEST_CALL(test_sweeper_overloads());
    ARIADNE_TEST_CALL(test_relative_sweeper());

    return ARIADNE_TEST_FAILURES;
}
