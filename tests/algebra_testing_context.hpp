#ifndef ARIADNE_ALGEBRA_TESTING_CONTEXT_HPP
#define ARIADNE_ALGEBRA_TESTING_CONTEXT_HPP
#include "ariadne-testing.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
namespace AlgebraTesting {
inline thread_local Ariadne::Testing::Context* active = nullptr;
inline Ariadne::Testing::Context& context() {
    if (!active) throw std::logic_error("No active algebra testing context");
    return *active;
}
struct CaseScope {
    Ariadne::Testing::Context& ctx;
    Ariadne::Testing::Context* previous;
    std::ostringstream output, errors;
    std::streambuf* oldout;
    std::streambuf* olderr;
    explicit CaseScope(Ariadne::Testing::Context& c)
        : ctx(c), previous(active),
          oldout(std::cout.rdbuf(output.rdbuf())),
          olderr(std::cerr.rdbuf(errors.rdbuf())) { active = &c; }
    CaseScope(CaseScope const&) = delete;
    CaseScope& operator=(CaseScope const&) = delete;
    ~CaseScope() noexcept {
        std::cout.rdbuf(oldout);
        std::cerr.rdbuf(olderr);
        active = previous;
        try {
            if (!output.str().empty()) ctx.note("stdout: " + output.str());
            if (!errors.str().empty()) ctx.note("stderr: " + errors.str());
        } catch (...) {}
    }
};

// Compatibility assertions intentionally preserve the acceptance semantics of
// the pre-migration ARIADNE_TEST_EQUAL(S) and ARIADNE_TEST_ASSERT macros.
// "decide" is an approximate/domain decision, NOT a mathematical proof.
// Use the standard Ariadne::Testing evidence macros for new tests instead.
template<class T>
bool legacy_decision(T const& value) {
    if constexpr (requires(T const& v) { decide(v); })
        return static_cast<bool>(decide(value));
    else
        return static_cast<bool>(value);
}

template<class A, class B>
void decided_equal(Ariadne::Testing::Context& context, A const& actual,
                   B const& expected, char const* description,
                   char const* actual_name, char const* expected_name) {
    bool const decision = legacy_decision(actual == expected);
    if (!decision) {
        context.inspect(actual_name, actual);
        context.inspect(expected_name, expected);
    }
    context.check([&] { return decision; }, description);
}

template<class A, class B>
void decided_same(Ariadne::Testing::Context& context, A const& actual,
                  B const& expected, char const* description,
                  char const* actual_name, char const* expected_name) {
    bool const decision = legacy_decision(same(actual, expected));
    if (!decision) {
        context.inspect(actual_name, actual);
        context.inspect(expected_name, expected);
    }
    context.check([&] { return decision; }, description);
}

template<class T>
void legacy_assert(Ariadne::Testing::Context& context, T const& proposition,
                   char const* expression) {
    using Ariadne::Testing::Detail::evidence;
    using Ariadne::Testing::Evidence;
    auto const observed = evidence(proposition);
    if (observed == Evidence::Indeterminate)
        context.warning(std::string("Legacy assertion indeterminate: ") + expression);
    context.expect_evidence(
        [&] { return observed; },
        Ariadne::Testing::EvidenceExpectation::PossiblyTrue, expression);
}
} // namespace AlgebraTesting

// These names deliberately mark legacy decision-based checks. Avoid using them
// in new tests that require definite three-valued mathematical evidence.
#define ARIADNE_ALGEBRA_DECIDED_EQUAL_CTX(context, actual, expected) \
    do { \
        auto&& a_ariadne = (actual); \
        auto&& b_ariadne = (expected); \
        ::AlgebraTesting::decided_equal((context), a_ariadne, b_ariadne, \
            #actual " == " #expected, #actual, #expected); \
    } while (false)

#define ARIADNE_ALGEBRA_DECIDED_SAME_CTX(context, actual, expected) \
    do { \
        auto&& a_ariadne = (actual); \
        auto&& b_ariadne = (expected); \
        ::AlgebraTesting::decided_same((context), a_ariadne, b_ariadne, \
            #actual " same " #expected, #actual, #expected); \
    } while (false)

#define ARIADNE_ALGEBRA_DECIDED_REL_CTX(context, actual, op, expected) \
    do { \
        auto&& a_ariadne = (actual); \
        auto&& b_ariadne = (expected); \
        bool const result_ariadne = ::AlgebraTesting::legacy_decision(a_ariadne op b_ariadne); \
        if (!result_ariadne) { \
            (context).inspect(#actual, a_ariadne); \
            (context).inspect(#expected, b_ariadne); \
        } \
        (context).check([&] { return result_ariadne; }, \
                        #actual " " #op " " #expected); \
    } while (false)

#define ARIADNE_ALGEBRA_LEGACY_ASSERT_CTX(context, expression) \
    ::AlgebraTesting::legacy_assert((context), (expression), #expression)

#define ARIADNE_ALGEBRA_DECIDED_BINARY_CTX(context, predicate, left, right) \
    do { \
        auto&& a_ariadne = (left); \
        auto&& b_ariadne = (right); \
        bool const result_ariadne = ::AlgebraTesting::legacy_decision(predicate(a_ariadne, b_ariadne)); \
        (context).check([&] { return result_ariadne; }, #predicate "(" #left ", " #right ")"); \
    } while (false)

#define ARIADNE_ALGEBRA_DECIDED_NOT_BINARY_CTX(context, predicate, left, right) \
    do { \
        auto&& a_ariadne = (left); \
        auto&& b_ariadne = (right); \
        bool const result_ariadne = ::AlgebraTesting::legacy_decision(!predicate(a_ariadne, b_ariadne)); \
        (context).check([&] { return result_ariadne; }, "!" #predicate "(" #left ", " #right ")"); \
    } while (false)

#endif
