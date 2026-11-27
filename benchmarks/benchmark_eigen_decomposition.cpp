/***************************************************************************
 *            benchmark_eigen_decomposition.cpp
 *
 *  Copyright  2026  Luca Geretti
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
 *  along with Ariadne Algebra. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file benchmark_eigen_decomposition.cpp
 *  @brief Deterministic end-to-end eigenvalue/eigenvector benchmarks.
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "algebra/eigen_decomposition.hpp"
#include "numeric/numeric.hpp"
#include "ariadne-benchmarking.hpp"

namespace {
using namespace Ariadne;
namespace Bench = Ariadne::Benchmarking;

/** @brief User-configurable measurement settings. */
struct Selection {
    unsigned iterations = 3u;
    unsigned samples = 7u;
    std::vector<unsigned> sizes{3u, 8u, 16u};
    std::vector<unsigned> bits{53u, 128u};
    bool csv = false;
};

/** @brief Parse distinct integral configuration values. */
std::vector<unsigned> parse_values(const String& input, unsigned minimum,
                                   const char* dimension) {
    const std::string text(input);
    if(text.empty() || text.front() == ',' || text.back() == ',')
        throw std::invalid_argument(std::string("invalid ") + dimension + " list");
    std::istringstream stream(text);
    std::string token;
    std::vector<unsigned> values;
    while(std::getline(stream, token, ',')) {
        const unsigned value = Bench::parse_unsigned<unsigned>(String(token));
        if(value < minimum ||
           std::find(values.begin(), values.end(), value) != values.end())
            throw std::invalid_argument(std::string("invalid or duplicate ") + dimension);
        values.push_back(value);
    }
    if(values.empty())
        throw std::invalid_argument("empty configuration list");
    return values;
}

/** @brief Extend Benchmarking's common command line. */
List<ArgumentParser> arguments(Selection& selection) {
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "i", "iterations", "Complete calls per sample, positive",
        [&selection](const String& input) -> VoidFunction {
            const unsigned value = Bench::parse_unsigned<unsigned>(input);
            if(value == 0u) throw std::invalid_argument("iterations must be positive");
            return [&selection, value] { selection.iterations = value; };
        }));
    parsers.append(make_valued_argument_parser(
        "n", "samples", "Odd samples, at least three",
        [&selection](const String& input) -> VoidFunction {
            const unsigned value = Bench::parse_unsigned<unsigned>(input);
            if(value < 3u || value % 2u == 0u)
                throw std::invalid_argument("samples must be odd and at least three");
            return [&selection, value] { selection.samples = value; };
        }));
    parsers.append(make_valued_argument_parser(
        "z", "sizes", "Distinct matrix orders between 3 and 64",
        [&selection](const String& input) -> VoidFunction {
            auto values = parse_values(input, 3u, "size");
            if(std::any_of(values.begin(), values.end(),
                           [](unsigned v) { return v > 64u; }))
                throw std::invalid_argument("matrix size must not exceed 64");
            return [&selection, values] { selection.sizes = values; };
        }));
    parsers.append(make_valued_argument_parser(
        "b", "bits", "Distinct MP bit precisions, at least 53",
        [&selection](const String& input) -> VoidFunction {
            auto values = parse_values(input, 53u, "bit precision");
            return [&selection, values] { selection.bits = values; };
        }));
    parsers.append(make_unvalued_argument_parser(
        "", "csv", "Print CSV instead of a formatted table",
        [&selection] { selection.csv = true; }));
    return parsers;
}

/** @brief Workloads select different currently implemented algorithmic paths. */
enum class Workload { Dense2, Triangular, Separated, Rotated };

/** @brief Stable workload labels for cross-commit comparisons. */
const char* label(Workload w) {
    switch(w) {
        case Workload::Dense2: return "dense_2x2";
        case Workload::Triangular: return "upper_triangular";
        case Workload::Separated: return "dense_separated";
        case Workload::Rotated: return "dense_rotated";
        default: throw std::invalid_argument("unknown workload");
    }
}

/** @brief Construct exact rational problems at the requested precision.
 *  @details Per-call construction is excluded from timed regions.
 */
template<class X, class PR>
Matrix<X> matrix_for(Workload w, unsigned n, PR precision) {
    const X zero(0u, precision);
    const X one(1u, precision);
    Matrix<X> A(n, n, zero);
    if(w == Workload::Dense2) {
        A[0u][0u] = X(2u, precision);
        A[0u][1u] = one;
        A[1u][0u] = one;
        A[1u][1u] = X(2u, precision);
        return A;
    }
    if(w == Workload::Triangular) {
        const X off = one / X(32u * n, precision);
        for(unsigned i = 0u; i != n; ++i) {
            A[i][i] = X(4u * i + 1u, precision);
            for(unsigned j = i + 1u; j != n; ++j)
                A[i][j] = off;
        }
    } else if(w == Workload::Separated) {
        const X off = one / X(32u * n, precision);
        for(unsigned i = 0u; i != n; ++i)
            for(unsigned j = 0u; j != n; ++j)
                A[i][j] = i == j ? X(4u * i + 1u, precision) : off;
    } else if(w == Workload::Rotated) {
        const X off = one / X(512u * n, precision);
        for(unsigned i = 0u; i != n; ++i) {
            A[i][i] = X(4u * (i / 2u) + 3u, precision);
            for(unsigned j = i + 1u; j != n; ++j) {
                A[i][j] = off;
                A[j][i] = off;
            }
        }
        for(unsigned i = 0u; i + 1u < n; i += 2u) {
            A[i][i + 1u] = one;
            A[i + 1u][i] = one;
        }
    } else {
        throw std::invalid_argument("unknown workload");
    }
    return A;
}

/** @brief Measurement and its semantic outcome. */
struct Datum {
    std::string backend;
    std::string workload;
    unsigned n;
    unsigned bits;
    EigenDecompositionStatus status;
    Bench::Statistics stats;
};

/** @brief Human-readable outcome. */
const char* status_name(EigenDecompositionStatus status) {
    switch(status) {
        case EigenDecompositionStatus::SUCCESS: return "SUCCESS";
        case EigenDecompositionStatus::INCONCLUSIVE: return "INCONCLUSIVE";
        case EigenDecompositionStatus::UNSUPPORTED_DIMENSION: return "UNSUPPORTED";
        default: return "INVALID";
    }
}

/** @brief Benchmark full eigendecomposition; verify status outside timing. */
template<class P, class PR>
void measure(const char* backend, unsigned bits, PR precision,
             Workload workload, unsigned size, const Selection& options,
             const Bench::Runner& runner, std::vector<Datum>& results) {
    using X = FloatType<P, PR>;
    using Status = EigenDecompositionStatus;
    const Matrix<X> input = matrix_for<X>(workload, size, precision);
    const EigenDecomposition<P, PR> solve(precision);
    const auto pilot = solve(input);
    const Status expected = pilot.status;

    // An inconclusive result in the rotated case is a genuine outcome
    // and must remain explicit in reports and A/B comparisons.
    if(expected == Status::UNSUPPORTED_DIMENSION ||
       (workload != Workload::Rotated && expected != Status::SUCCESS) ||
       (expected == Status::SUCCESS && !pilot.data)) {
        throw std::runtime_error(std::string("unexpected outcome for ") +
                                 backend + "/" + label(workload));
    }

    Status observed = expected;
    const std::string name = std::string(backend) + "/" + label(workload) +
                             "/n=" + std::to_string(size);
    const auto samples = runner.run(name, options.samples, [] {}, [&] {
            for(unsigned i = 0u; i != options.iterations; ++i) {
                const auto output = solve(input);
                observed = output.status;
                Bench::keep_alive(output);
            }
        }, [&] {
            if(observed != expected)
                throw std::runtime_error("inconsistent benchmark status: " + name);
        }, static_cast<double>(options.iterations));

    auto statistics = Bench::summarize(samples);
    if(!std::isfinite(statistics.median) || statistics.median <= 0.0)
        throw std::runtime_error("invalid benchmark duration: " + name);
    results.push_back({backend, label(workload), size, bits,
                       expected, std::move(statistics)});
}

/** @brief Run each workload for the selected paradigm and precision. */
template<class P, class PR>
void measure_backend(const char* backend, unsigned bits, PR precision,
                     const Selection& selection, const Bench::Runner& runner,
                     std::vector<Datum>& results) {
    measure<P>(backend, bits, precision, Workload::Dense2, 2u,
               selection, runner, results);
    for(unsigned size: selection.sizes) {
        measure<P>(backend, bits, precision, Workload::Triangular, size,
                   selection, runner, results);
        measure<P>(backend, bits, precision, Workload::Separated, size,
                   selection, runner, results);
        measure<P>(backend, bits, precision, Workload::Rotated, size,
                   selection, runner, results);
    }
}

/** @brief Render a tabular or CSV result matrix. */
void report(const std::vector<Datum>& rows, const Selection& options,
            const Bench::SummaryFields& fields) {
    if(options.csv) {
        LOGGING_PRINTLN("backend,workload,n,bits,status,median_ns,mad_ns,min_ns,max_ns,p95_ns")
        for(const auto& row: rows) {
            const auto& s = row.stats;
            std::ostringstream record;
            record << std::setprecision(17)
                   << row.backend << ',' << row.workload << ','
                   << row.n << ',' << row.bits << ',' << status_name(row.status) << ','
                   << s.median << ',' << s.mad << ',' << s.minimum << ','
                   << s.maximum << ',' << s.p95;
            LOGGING_PRINTLN(record.str())
        }
    } else {
        Bench::Report table("EigenDecomposition: nanoseconds per call", fields);
        table.set_dimension_columns({
            {"Backend"}, {"Workload"},
            {"n", 0u, Bench::TableAlignment::Right},
            {"Bits", 0u, Bench::TableAlignment::Right},
            {"Outcome"}
        });
        auto style = Bench::TableOptions::ruled();
        style.group_columns = {0u, 1u};
        table.set_table_options(style);
        for(const auto& row: rows) {
            table.add_row({row.backend, row.workload, std::to_string(row.n),
                           std::to_string(row.bits), status_name(row.status)},
                          row.stats);
        }
        LOGGING_PRINTLN("End-to-end eigenvalue + eigenvector decomposition; "
                        "matrix construction excluded.")
        LOGGING_PRINTLN("iterations=" << options.iterations
                                    << ", samples=" << options.samples
                                    << " (warm-ups via --warmups).")
        table.print();
    }
}

} // namespace

int main(int argc, const char* argv[]) {
    try {
        Selection selection;
        Bench::CommandLine command_line;
        command_line.add(CommandLineContribution("algebra-eigen-benchmark",
                                                  arguments(selection)));
        if(!command_line.acquire(argc, argv))
            return 0;

        const std::size_t configurations =
            (1u + 3u * selection.sizes.size()) *
            (2u + 2u * selection.bits.size());
        const Bench::Runner runner(command_line, configurations);
        std::vector<Datum> results;
        results.reserve(configurations);

        measure_backend<ApproximateTag>("dp/approx", 53u, dp,
                                        selection, runner, results);
        measure_backend<ValidatedTag>("dp/validated", 53u, dp,
                                      selection, runner, results);
        for(unsigned bits: selection.bits) {
            const MultiplePrecision precision(bits);
            measure_backend<ApproximateTag>("mp/approx", bits, precision,
                                            selection, runner, results);
            measure_backend<ValidatedTag>("mp/validated", bits, precision,
                                          selection, runner, results);
        }
        if(results.size() != configurations)
            throw std::runtime_error("missing benchmark configurations");
        report(results, selection, command_line.summary_fields());
        return 0;
    } catch(const std::exception& error) {
        LOGGING_PRINTLN("benchmark_eigen_decomposition: " << error.what())
        return 1;
    }
}
