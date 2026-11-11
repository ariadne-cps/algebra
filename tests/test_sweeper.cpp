/***************************************************************************
 *            test_sweeper.cpp
 *
 *  Copyright  2026  Pieter Collins
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

#include <cmath>

#include "algebra/expansion.hpp"
#include "algebra/sweeper.hpp"
#include "algebra_testing_context.hpp"
#include "interval/interval.hpp"
#include "numeric/numeric.hpp"

using namespace Ariadne;

namespace {

Void test_threshold_sweeper() {
    MultiplePrecision multiple_precision(128);

    ThresholdSweeper<FloatDP> sweeper_dp(dp, 1e-8);
    ThresholdSweeper<FloatMP> sweeper_mp(multiple_precision, std::pow(2.0, -64));

    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "sweeper_dp.precision()", sweeper_dp.precision());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "sweeper_mp.precision()", sweeper_mp.precision());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "sweeper_dp", sweeper_dp);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "sweeper_mp", sweeper_mp);

    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)ThresholdSweeper<FloatDP>(dp, -1e-8));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)ThresholdSweeper<FloatDP>(dp, FloatDP(-1, dp)));

    MultiIndex constant_index({0u});
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), sweeper_dp.discard(constant_index, FloatDP(0u, dp)));
    ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(AlgebraTesting::context(), !sweeper_dp.discard(constant_index, FloatDP(1u, dp)));

    SweeperInterface<FloatDP>* cloned_interface = sweeper_dp._copy();
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "cloned_interface->precision()", cloned_interface->precision());
    delete cloned_interface;
}

Void test_sweeper_handle() {
    Sweeper<FloatDP> default_dp_sweeper;
    Sweeper<FloatMP> default_mp_sweeper;
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "default_dp_sweeper.precision()", default_dp_sweeper.precision());
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "default_mp_sweeper.precision()", default_mp_sweeper.precision());

    ThresholdSweeper<FloatDP> threshold_sweeper(dp, 1e-8);
    Sweeper<FloatDP> threshold_base(threshold_sweeper);
    Sweeper<FloatDP> threshold_copy(threshold_base);
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "threshold_copy.precision()", threshold_copy.precision());

    Expansion<MultiIndex, FloatDP> sweep_terms(1u, dp);
    sweep_terms.append(MultiIndex({0u}), FloatDP(0u, dp));
    sweep_terms.append(MultiIndex({1u}), FloatDP(1u, dp));
    FloatDPError sweep_error(0u, dp);

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), threshold_copy.sweep(sweep_terms, sweep_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sweep_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), sweep_terms.front().coefficient(), FloatDP(1u, dp));
}

Void test_sweeper_overloads() {
    ThresholdSweeper<FloatDP> threshold_sweeper(dp, FloatDP(0.25_x, dp));
    const SweeperInterface<FloatDP>& threshold_interface = threshold_sweeper;

    Expansion<MultiIndex, FloatDP> float_terms(1u, dp);
    float_terms.append(MultiIndex({0u}), FloatDP(0.125_x, dp));
    float_terms.append(MultiIndex({1u}), FloatDP(1.0_x, dp));
    FloatDPError float_error(0u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), threshold_interface.sweep(float_terms, float_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), float_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), float_terms.front().index(), MultiIndex({1u}));

    Expansion<MultiIndex, FloatDPBounds> bounds_terms(1u, dp);
    bounds_terms.append(MultiIndex({0u}), FloatDPBounds(0.125_x, dp));
    bounds_terms.append(MultiIndex({1u}), FloatDPBounds(1.0_x, dp));

    Expansion<MultiIndex, FloatDPApproximation> approximation_terms(bounds_terms);
    Expansion<MultiIndex, FloatDPUpperInterval> upper_terms(bounds_terms);

    FloatDPError bounds_error(0u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), threshold_interface.sweep(bounds_terms, bounds_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), bounds_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), bounds_terms.front().index(), MultiIndex({1u}));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), threshold_interface.sweep(approximation_terms));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), approximation_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), approximation_terms.front().index(), MultiIndex({1u}));

    FloatDPError upper_error(0u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), threshold_interface.sweep(upper_terms, upper_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), upper_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), upper_terms.front().index(), MultiIndex({1u}));
}

Void test_relative_sweeper() {
    RelativeThresholdSweeper<FloatDP> relative_sweeper(dp, FloatDP(0.25_x, dp));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), relative_sweeper.relative_sweep_threshold(), FloatDP(0.25_x, dp));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(AlgebraTesting::context(), (void)RelativeThresholdSweeper<FloatDP>(dp, FloatDP(0u, dp)));
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "relative_sweeper", relative_sweeper);

    const SweeperInterface<FloatDP>& relative_interface = relative_sweeper;
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "relative_interface.precision()", relative_interface.precision());

    SweeperInterface<FloatDP>* relative_copy = relative_sweeper._copy();
    ARIADNE_TESTING_INSPECT_CTX(AlgebraTesting::context(), "relative_copy->precision()", relative_copy->precision());
    delete relative_copy;

    Expansion<MultiIndex, FloatDP> float_terms(1u, dp);
    float_terms.append(MultiIndex({0u}), FloatDP(0.125_x, dp));
    float_terms.append(MultiIndex({1u}), FloatDP(1.0_x, dp));
    FloatDPError float_error(0u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), relative_interface.sweep(float_terms, float_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), float_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), float_terms.front().index(), MultiIndex({1u}));

    Expansion<MultiIndex, FloatDPBounds> bounds_source(1u, dp);
    bounds_source.append(MultiIndex({0u}), FloatDPBounds(0.125_x, dp));
    bounds_source.append(MultiIndex({1u}), FloatDPBounds(1.0_x, dp));

    Expansion<MultiIndex, FloatDPBounds> bounds_terms(bounds_source);
    Expansion<MultiIndex, FloatDPApproximation> approximation_terms(bounds_source);
    Expansion<MultiIndex, FloatDPUpperInterval> upper_terms(bounds_source);

    FloatDPError bounds_error(0u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), relative_interface.sweep(bounds_terms, bounds_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), bounds_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), bounds_terms.front().index(), MultiIndex({1u}));

    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), relative_interface.sweep(approximation_terms));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), approximation_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), approximation_terms.front().index(), MultiIndex({1u}));

    FloatDPError upper_error(0u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), relative_interface.sweep(upper_terms, upper_error));
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), upper_terms.size(), 1u);
    ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(AlgebraTesting::context(), upper_terms.front().index(), MultiIndex({1u}));
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("sweeper");
    suite.add("legacy tests", [](Ariadne::Testing::Context& ctx) {
        AlgebraTesting::CaseScope capture(ctx);

        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_threshold_sweeper());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_sweeper_handle());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_sweeper_overloads());
        ARIADNE_TESTING_EXPECT_NO_THROW_CTX(AlgebraTesting::context(), test_relative_sweeper());
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
