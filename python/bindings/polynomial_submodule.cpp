/***************************************************************************
 *            polynomial_submodule.cpp
 *
 *  Python bindings for polynomial representations in Algebra.
 ****************************************************************************/

#include "pybind11.hpp"

#include "algebra-utilities.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/polynomial.hpp"
#include "algebra/chebyshev_polynomial.hpp"

using namespace Ariadne;

namespace Ariadne {

template<> struct PythonTemplateName<MultivariatePolynomial> {
    static std::string get() { return "MultivariatePolynomial"; }
};

template<class X> struct PythonClassName<MultivariatePolynomial<X>> {
    static std::string get() { return python_template_class_name<X>("MultivariatePolynomial"); }
};

} // namespace Ariadne

namespace {

template<class X>
pybind11::class_<MultivariatePolynomial<X>>
export_multivariate_polynomial(pybind11::module& module) {
    using P=MultivariatePolynomial<X>;
    using PR=typename X::PrecisionType;

    pybind11::class_<P> polynomial_class(module,python_class_name<P>().c_str());

    polynomial_class.def(pybind11::init<P>());
    polynomial_class.def(pybind11::init<SizeType,PR>());
    polynomial_class.def(pybind11::init<Expansion<MultiIndex,X>>());

    polynomial_class.def_static("constant",(P(*)(SizeType,X const&))&P::constant);
    polynomial_class.def_static("variable",(P(*)(SizeType,SizeType,PR))&P::variable);
    polynomial_class.def_static("coordinate",(P(*)(SizeType,SizeType,PR))&P::variable);
    polynomial_class.def_static("variables",[](SizeType as, PR pr) {
        auto const vars=P::variables(as,pr);
        pybind11::list result;
        for(SizeType i=0; i!=vars.size(); ++i) {
            result.append(pybind11::cast(vars[i]));
        }
        return result;
    });

    polynomial_class.def("argument_size",&P::argument_size);
    polynomial_class.def("degree",&P::degree);
    polynomial_class.def("number_of_terms",&P::number_of_terms);
    polynomial_class.def("value",&P::value);
    polynomial_class.def("expansion",(Expansion<MultiIndex,X>const&(P::*)()const)&P::expansion);
    polynomial_class.def("__getitem__",[](P const& p, MultiIndex const& a) { return X(p[a]); });
    polynomial_class.def("__setitem__",[](P& p, MultiIndex const& a, X const& x) { p[a]=x; });
    polynomial_class.def("__call__",[](P const& p, Vector<X> const& x) { return evaluate(p,x); });

    define_algebra(module,polynomial_class);
    polynomial_class.def("__str__",&__cstr__<P>);

    module.def("derivative",[](P const& p, SizeType i) { return derivative(p,i); });
    module.def("antiderivative",[](P const& p, SizeType i) { return antiderivative(p,i); });
    module.def("truncate",[](P const& p, DegreeType d) { return truncate(p,d); });
    module.def("evaluate",[](P const& p, Vector<X> const& x) { return evaluate(p,x); });

    export_vector<P>(module,python_template_class_name<X>("MultivariatePolynomialVector"));

    return polynomial_class;
}

template<class X>
pybind11::class_<UnivariateChebyshevPolynomial<X>>
export_univariate_chebyshev_polynomial(pybind11::module& module) {
    using C=UnivariateChebyshevPolynomial<X>;
    using PR=typename X::PrecisionType;

    const std::string name=python_template_class_name<X>("UnivariateChebyshevPolynomial");
    pybind11::class_<C> chebyshev_class(module,name.c_str());

    chebyshev_class.def(pybind11::init<PR>());
    chebyshev_class.def_static("constant",(C(*)(X const&))&C::constant);
    chebyshev_class.def_static("coordinate",(C(*)(PR))&C::coordinate);
    chebyshev_class.def_static("basis",(C(*)(SizeType,PR))&C::basis);

    chebyshev_class.def("degree",&C::degree);
    chebyshev_class.def("sup_norm",&C::sup_norm);
    chebyshev_class.def("__getitem__",&C::operator[]);
    chebyshev_class.def("__call__",[](C const& p, X const& x) { return p(x); });

    define_algebra(module,chebyshev_class);
    chebyshev_class.def("__str__",&__cstr__<C>);

    module.def("evaluate",[](C const& p, X const& x) { return p(x); });
    module.def("sup_norm",[](C const& p) { return p.sup_norm(); });

    return chebyshev_class;
}

template<class X>
pybind11::class_<MultivariateChebyshevPolynomial<X>>
export_multivariate_chebyshev_polynomial(pybind11::module& module) {
    using C=MultivariateChebyshevPolynomial<X>;
    using PR=typename X::PrecisionType;

    const std::string name=python_template_class_name<X>("MultivariateChebyshevPolynomial");
    pybind11::class_<C> chebyshev_class(module,name.c_str());

    chebyshev_class.def(pybind11::init<SizeType,PR>());
    chebyshev_class.def_static("constant",(C(*)(SizeType,X const&))&C::constant);
    chebyshev_class.def_static("coordinate",(C(*)(SizeType,SizeType,PR))&C::coordinate);
    chebyshev_class.def_static("basis",(C(*)(SizeType,SizeType,SizeType,PR))&C::basis);

    chebyshev_class.def("argument_size",&C::argument_size);
    chebyshev_class.def("sup_norm",&C::sup_norm);
    chebyshev_class.def("__call__",[](C const& p, Vector<X> const& x) { return p(x); });

    define_algebra(module,chebyshev_class);
    chebyshev_class.def("__str__",&__cstr__<C>);

    module.def("evaluate",[](C const& p, Vector<X> const& x) { return p(x); });
    module.def("sup_norm",[](C const& p) { return p.sup_norm(); });

    return chebyshev_class;
}

} // namespace

void polynomial_submodule(pybind11::module& module) {
    export_multivariate_polynomial<FloatDPApproximation>(module);
    export_multivariate_polynomial<FloatDPBounds>(module);
    export_multivariate_polynomial<FloatMPApproximation>(module);
    export_multivariate_polynomial<FloatMPBounds>(module);

    template_<MultivariatePolynomial> polynomial_template(module,"MultivariatePolynomial");
    polynomial_template.instantiate<FloatDPApproximation>();
    polynomial_template.instantiate<FloatDPBounds>();
    polynomial_template.instantiate<FloatMPApproximation>();
    polynomial_template.instantiate<FloatMPBounds>();

    export_univariate_chebyshev_polynomial<FloatDPApproximation>(module);
    export_univariate_chebyshev_polynomial<FloatMPApproximation>(module);
    export_multivariate_chebyshev_polynomial<FloatDPApproximation>(module);
    export_multivariate_chebyshev_polynomial<FloatMPApproximation>(module);
}
