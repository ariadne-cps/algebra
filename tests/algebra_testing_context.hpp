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
}
#endif
