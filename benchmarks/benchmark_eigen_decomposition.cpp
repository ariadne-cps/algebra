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
#include <type_traits>
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
    bool profile = false;
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
    parsers.append(make_unvalued_argument_parser(
        "", "profile", "Collect separate candidate and verification timings",
        [&selection] { selection.profile = true; }));
    return parsers;
}

/** @brief Workloads select different currently implemented algorithmic paths. */
enum class Workload { Dense2, Triangular, Separated, Rotated, Nonnormal4, NonnormalInterval4 };

/** @brief Stable workload labels for cross-commit comparisons. */
const char* label(Workload w) {
    switch(w) {
        case Workload::Dense2: return "dense_2x2";
        case Workload::Triangular: return "upper_triangular";
        case Workload::Separated: return "dense_separated";
        case Workload::Rotated: return "dense_rotated";
        case Workload::Nonnormal4: return "dense_nonnormal_4x4";
        case Workload::NonnormalInterval4: return "dense_nonnormal_interval_4x4";
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
    if(w == Workload::Nonnormal4 || w == Workload::NonnormalInterval4) {
        if(n != 4u)
            throw std::invalid_argument("nonnormal workload has fixed order four");
        // Exact characteristic roots {1, 3, 5, 7}; original Gershgorin
        // disks overlap. Match the independent functional regression.
        const Int coefficients[4][4] = {
            {1, 2, 0, 0},
            {0, 3, 0, 0},
            {-8, 10, 5, 0},
            {-6, 2, 0, 7}
        };
        for(unsigned i = 0u; i != 4u; ++i)
            for(unsigned j = 0u; j != 4u; ++j)
                A[i][j] = X(coefficients[i][j], precision);
        if(w == Workload::NonnormalInterval4) {
            if constexpr(std::is_same_v<typename X::Paradigm, ValidatedTag>) {
                A[0u][0u] += X(0, 1, precision) / X(1048576u, precision);
            } else {
                throw std::invalid_argument("interval workload requires ValidatedTag");
            }
        }
        return A;
    }
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
    EigenDecompositionStrategy strategy;
    Bench::Statistics stats;
    double candidate_median_ns = 0.0;
    double verification_median_ns = 0.0;
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

/** @brief Name the *actual* last attempted or successful algorithm. */
const char* strategy_name(EigenDecompositionStrategy strategy) {
    switch(strategy) {
        case EigenDecompositionStrategy::NONE: return "none";
        case EigenDecompositionStrategy::ANALYTIC: return "analytic";
        case EigenDecompositionStrategy::TRIANGULAR: return "triangular";
        case EigenDecompositionStrategy::JACOBI: return "jacobi";
        case EigenDecompositionStrategy::GERSHGORIN: return "gershgorin";
        case EigenDecompositionStrategy::SIMILARITY: return "similarity";
        case EigenDecompositionStrategy::KRAWCZYK: return "krawczyk";
        default: return "invalid";
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
    EigenDecompositionDiagnostics pilot_diagnostics;
    const auto pilot = options.profile
        ? solve(input, &pilot_diagnostics)
        : solve(input);
    const Status expected = pilot.status;

    // An inconclusive result in the rotated case is a genuine outcome.
    // Nonnormal regression cases require successful validated certification,
    // so a failure is never silently measured as a "speedup".
    if(expected == Status::UNSUPPORTED_DIMENSION ||
       (workload != Workload::Rotated && expected != Status::SUCCESS) ||
       (expected == Status::SUCCESS && !pilot.data)) {
        throw std::runtime_error(std::string("unexpected outcome for ") +
                                 backend + "/" + label(workload));
    }

    Status observed = expected;
    const std::string name = std::string(backend) + "/" + label(workload) +
                             "/n=" + std::to_string(size);
    Bench::Samples candidate_samples;
    Bench::Samples verification_samples;
    double round_candidate = 0.0;
    double round_verification = 0.0;
    unsigned completed_rounds = 0u;
    const auto samples = runner.run(name, options.samples, [&] {
            round_candidate = 0.0;
            round_verification = 0.0;
        }, [&] {
            for(unsigned i = 0u; i != options.iterations; ++i) {
                EigenDecompositionDiagnostics diagnostic;
                const auto output = options.profile
                    ? solve(input, &diagnostic)
                    : solve(input);
                observed = output.status;
                if(options.profile) {
                    if(diagnostic.strategy != pilot_diagnostics.strategy)
                        throw std::runtime_error("inconsistent benchmark strategy: " + name);
                    round_candidate += static_cast<double>(
                        diagnostic.candidate_time.count());
                    round_verification += static_cast<double>(
                        diagnostic.verification_time.count());
                }
                Bench::keep_alive(output);
            }
        }, [&] {
            if(observed != expected)
                throw std::runtime_error("inconsistent benchmark status: " + name);
            // The Runner validates warm-ups as well, but reports samples
            // only after warm-ups. Apply exactly the same filtering here.
            if(options.profile && completed_rounds >= runner.warmups()) {
                candidate_samples.add(round_candidate / options.iterations);
                verification_samples.add(round_verification / options.iterations);
            }
            ++completed_rounds;
        }, static_cast<double>(options.iterations));

    auto statistics = Bench::summarize(samples);
    if(!std::isfinite(statistics.median) || statistics.median <= 0.0)
        throw std::runtime_error("invalid benchmark duration: " + name);
    double candidate_median = 0.0;
    double verification_median = 0.0;
    if(options.profile) {
        if(candidate_samples.size() != options.samples ||
           verification_samples.size() != options.samples)
            throw std::runtime_error("incomplete profiling samples: " + name);
        candidate_median = Bench::summarize(candidate_samples).median;
        verification_median = Bench::summarize(verification_samples).median;
    }
    results.push_back({backend, label(workload), size, bits,
                       expected, pilot_diagnostics.strategy, std::move(statistics),
                       candidate_median, verification_median});
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
    // Compare the same non-normal point matrix in both paradigms.
    measure<P>(backend, bits, precision, Workload::Nonnormal4, 4u,
               selection, runner, results);
    if constexpr(std::is_same_v<P, ValidatedTag>) {
        measure<P>(backend, bits, precision, Workload::NonnormalInterval4, 4u,
                   selection, runner, results);
    }
}

/** @brief Format benchmark statistics as integral nanoseconds for display.
 *  @details Measurements remain full-precision doubles. Variance is in ns^2,
 *  and sample count remains an unscaled integer.
 */
std::string integer_statistic(const Bench::Statistics& result,
                              Bench::SummaryField field) {
    switch(field) {
        case Bench::SummaryField::Minimum:
            return Bench::Table::number(result.minimum, 0u);
        case Bench::SummaryField::Maximum:
            return Bench::Table::number(result.maximum, 0u);
        case Bench::SummaryField::Mean:
            return Bench::Table::number(result.mean, 0u);
        case Bench::SummaryField::MAD:
            return Bench::Table::number(result.mad, 0u);
        case Bench::SummaryField::Median:
            return Bench::Table::number(result.median, 0u);
        case Bench::SummaryField::P95:
            return Bench::Table::number(result.p95, 0u);
        case Bench::SummaryField::Stddev:
            return Bench::Table::number(result.stddev, 0u);
        case Bench::SummaryField::Variance:
            return Bench::Table::number(result.variance, 0u);
        case Bench::SummaryField::Count:
            return std::to_string(result.count);
        default:
            throw std::invalid_argument("unknown benchmark statistic");
    }
}

/** @brief Render a tabular or CSV result matrix. */
void report(const std::vector<Datum>& rows, const Selection& options,
            const Bench::SummaryFields& fields) {
    if(options.csv) {
        if(options.profile) {
            LOGGING_PRINTLN("backend,workload,n,bits,status,median_ns,mad_ns,min_ns,max_ns,p95_ns,strategy,candidate_median_ns,verification_median_ns")
        } else {
            LOGGING_PRINTLN("backend,workload,n,bits,status,median_ns,mad_ns,min_ns,max_ns,p95_ns")
        }
        for(const auto& row: rows) {
            const auto& s = row.stats;
            std::ostringstream record;
            record << std::fixed << std::setprecision(0)
                   << row.backend << ',' << row.workload << ','
                   << row.n << ',' << row.bits << ',' << status_name(row.status) << ','
                   << s.median << ',' << s.mad << ',' << s.minimum << ','
                   << s.maximum << ',' << s.p95;
            if(options.profile) {
                record << ',' << strategy_name(row.strategy) << ','
                       << row.candidate_median_ns << ','
                       << row.verification_median_ns;
            }
            LOGGING_PRINTLN(record.str())
        }
    } else {
        std::vector<Bench::TableColumn> columns = {
            {"Backend"}, {"Workload"},
            {"n", 0u, Bench::TableAlignment::Right},
            {"Bits", 0u, Bench::TableAlignment::Right},
            {"Outcome"}
        };
        if(options.profile) {
            columns.push_back({"Strategy"});
            columns.push_back({"Candidate ns", 0u, Bench::TableAlignment::Right});
            columns.push_back({"Verify ns", 0u, Bench::TableAlignment::Right});
        }
        for(Bench::SummaryField field: fields) {
            columns.push_back({std::string(Bench::summary_field_heading(field)),
                               0u, Bench::TableAlignment::Right});
        }
        auto style = Bench::TableOptions::ruled();
        style.group_columns = {0u, 1u};
        Bench::Table table(std::move(columns), style);
        for(const auto& row: rows) {
            std::vector<std::string> cells = {
                row.backend, row.workload, std::to_string(row.n),
                std::to_string(row.bits), status_name(row.status)
            };
            if(options.profile) {
                cells.push_back(strategy_name(row.strategy));
                cells.push_back(Bench::Table::number(row.candidate_median_ns, 0u));
                cells.push_back(Bench::Table::number(row.verification_median_ns, 0u));
            }
            for(Bench::SummaryField field: fields) {
                cells.push_back(integer_statistic(row.stats, field));
            }
            table.add_row(std::move(cells));
        }
        LOGGING_PRINTLN("End-to-end eigenvalue + eigenvector decomposition; "
                        "matrix construction excluded.")
        LOGGING_PRINTLN("iterations=" << options.iterations
                                    << ", samples=" << options.samples
                                    << " (warm-ups via --warmups).")
        if(options.profile) {
            LOGGING_PRINTLN("Profiling adds timer overhead; candidate and verification times "
                            "exclude dispatch/bookkeeping and need not sum to total.")
        }
        LOGGING_PRINTLN("EigenDecomposition: nanoseconds per call")
        LOGGING_PRINTLN(table.str())
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
                (2u + 2u * selection.bits.size()) +
            3u * (1u + selection.bits.size());
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
