/***************************************************************************
 *            polynomial_submodule.cpp
 *
 *  Python bindings for polynomial representations in Algebra.
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

#include "algebra/chebyshev_polynomial.hpp"
#include "algebra/expansion.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/polynomial.hpp"
#include "algebra/vector.hpp"
#include "numeric/numeric.hpp"
#include "pybind11.hpp"

#include "algebra-utilities.hpp"

using namespace Ariadne;

namespace Ariadne {

template<>
struct PythonTemplateName<MultivariatePolynomial> {
    static std::string get() { return "MultivariatePolynomial"; }
};

template<class X>
struct PythonClassName<MultivariatePolynomial<X>> {
    static std::string get() { return python_template_class_name<X>("MultivariatePolynomial"); }
};

} // namespace Ariadne

namespace {

template<class X>
pybind11::class_<MultivariatePolynomial<X>>
    export_multivariate_polynomial(pybind11::module& module) {
    using P = MultivariatePolynomial<X>;
    using PR = typename X::PrecisionType;

    pybind11::class_<P> polynomial_class(module, python_class_name<P>().c_str());

    polynomial_class.def(pybind11::init<P>());
    polynomial_class.def(pybind11::init<SizeType, PR>());
    polynomial_class.def(pybind11::init<Expansion<MultiIndex, X>>());

    polynomial_class.def_static("constant", (P (*)(SizeType, const X&))&P::constant);
    polynomial_class.def_static("variable", (P (*)(SizeType, SizeType, PR))&P::variable);
    polynomial_class.def_static("coordinate", (P (*)(SizeType, SizeType, PR))&P::variable);
    polynomial_class.def_static("variables", [](SizeType as, PR pr) {
        const auto vars = P::variables(as, pr);
        pybind11::list result;
        for (SizeType i = 0; i != vars.size(); ++i) {
            result.append(pybind11::cast(vars[i]));
        }
        return result;
    });

    polynomial_class.def("argument_size", &P::argument_size);
    polynomial_class.def("degree", &P::degree);
    polynomial_class.def("number_of_terms", &P::number_of_terms);
    polynomial_class.def("value", &P::value);
    polynomial_class.def("expansion", (const Expansion<MultiIndex, X>& (P::*)() const) & P::expansion);
    polynomial_class.def("__getitem__", [](const P& p, const MultiIndex& a) { return X(p[a]); });
    polynomial_class.def("__setitem__", [](P& p, const MultiIndex& a, const X& x) { p[a] = x; });
    polynomial_class.def("__call__", [](const P& p, const Vector<X>& x) { return evaluate(p, x); });

    define_algebra(module, polynomial_class);
    polynomial_class.def("__str__", &__cstr__<P>);

    module.def("derivative", [](const P& p, SizeType i) { return derivative(p, i); });
    module.def("antiderivative", [](const P& p, SizeType i) { return antiderivative(p, i); });
    module.def("truncate", [](const P& p, DegreeType d) { return truncate(p, d); });
    module.def("evaluate", [](const P& p, const Vector<X>& x) { return evaluate(p, x); });

    export_vector<P>(module, python_template_class_name<X>("MultivariatePolynomialVector"));

    return polynomial_class;
}

template<class X>
pybind11::class_<UnivariateChebyshevPolynomial<X>>
    export_univariate_chebyshev_polynomial(pybind11::module& module) {
    using C = UnivariateChebyshevPolynomial<X>;
    using PR = typename X::PrecisionType;

    const std::string name = python_template_class_name<X>("UnivariateChebyshevPolynomial");
    pybind11::class_<C> chebyshev_class(module, name.c_str());

    chebyshev_class.def(pybind11::init<PR>());
    chebyshev_class.def_static("constant", (C (*)(const X&))&C::constant);
    chebyshev_class.def_static("coordinate", (C (*)(PR))&C::coordinate);
    chebyshev_class.def_static("basis", (C (*)(SizeType, PR))&C::basis);

    chebyshev_class.def("degree", &C::degree);
    chebyshev_class.def("sup_norm", &C::sup_norm);
    chebyshev_class.def("__getitem__", &C::operator[]);
    chebyshev_class.def("__call__", [](const C& p, const X& x) { return p(x); });

    define_algebra(module, chebyshev_class);
    chebyshev_class.def("__str__", &__cstr__<C>);

    module.def("evaluate", [](const C& p, const X& x) { return p(x); });
    module.def("sup_norm", [](const C& p) { return p.sup_norm(); });

    return chebyshev_class;
}

template<class X>
pybind11::class_<MultivariateChebyshevPolynomial<X>>
    export_multivariate_chebyshev_polynomial(pybind11::module& module) {
    using C = MultivariateChebyshevPolynomial<X>;
    using PR = typename X::PrecisionType;

    const std::string name = python_template_class_name<X>("MultivariateChebyshevPolynomial");
    pybind11::class_<C> chebyshev_class(module, name.c_str());

    chebyshev_class.def(pybind11::init<SizeType, PR>());
    chebyshev_class.def_static("constant", (C (*)(SizeType, const X&))&C::constant);
    chebyshev_class.def_static("coordinate", (C (*)(SizeType, SizeType, PR))&C::coordinate);
    chebyshev_class.def_static("basis", (C (*)(SizeType, SizeType, SizeType, PR))&C::basis);

    chebyshev_class.def("argument_size", &C::argument_size);
    chebyshev_class.def("sup_norm", &C::sup_norm);
    chebyshev_class.def("__call__", [](const C& p, const Vector<X>& x) { return p(x); });

    define_algebra(module, chebyshev_class);
    chebyshev_class.def("__str__", &__cstr__<C>);

    module.def("evaluate", [](const C& p, const Vector<X>& x) { return p(x); });
    module.def("sup_norm", [](const C& p) { return p.sup_norm(); });

    return chebyshev_class;
}

} // namespace

void polynomial_submodule(pybind11::module& module) {
    export_multivariate_polynomial<FloatDPApproximation>(module);
    export_multivariate_polynomial<FloatDPBounds>(module);
    export_multivariate_polynomial<FloatMPApproximation>(module);
    export_multivariate_polynomial<FloatMPBounds>(module);

    template_<MultivariatePolynomial> polynomial_template(module, "MultivariatePolynomial");
    polynomial_template.instantiate<FloatDPApproximation>();
    polynomial_template.instantiate<FloatDPBounds>();
    polynomial_template.instantiate<FloatMPApproximation>();
    polynomial_template.instantiate<FloatMPBounds>();

    export_univariate_chebyshev_polynomial<FloatDPApproximation>(module);
    export_univariate_chebyshev_polynomial<FloatMPApproximation>(module);
    export_multivariate_chebyshev_polynomial<FloatDPApproximation>(module);
    export_multivariate_chebyshev_polynomial<FloatMPApproximation>(module);
}
