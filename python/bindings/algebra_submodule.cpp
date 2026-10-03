/***************************************************************************
 *            algebra_submodule.cpp
 *
 *  Python bindings for core Algebra data structures.
 ****************************************************************************/

#include "pybind11.hpp"
#include <pybind11/stl.h>

#include "algebra-utilities.hpp"

#include "numeric/numeric.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/expansion.tpl.hpp"
#include "algebra/series.hpp"

using namespace Ariadne;

namespace {

MultiIndex multi_index_from_python(pybind11::sequence const& seq) {
    MultiIndex a(static_cast<SizeType>(seq.size()));
    for(SizeType i=0; i!=a.size(); ++i) {
        a.set(i,pybind11::cast<DegreeType>(seq[i]));
    }
    return a;
}

MultiIndex multi_index_from_object(pybind11::handle obj) {
    if(pybind11::isinstance<MultiIndex>(obj)) {
        return pybind11::cast<MultiIndex>(obj);
    }
    return multi_index_from_python(pybind11::reinterpret_borrow<pybind11::sequence>(obj));
}

template<class X>
Expansion<MultiIndex,X>
expansion_from_python(pybind11::dict const& terms, typename X::PrecisionType pr) {
    if(terms.empty()) {
        throw std::invalid_argument("An Expansion constructed from a dictionary needs at least one term.");
    }

    auto first=terms.begin();
    MultiIndex first_index=multi_index_from_object(first->first);
    Expansion<MultiIndex,X> result(first_index.size(),pr);

    for(auto item : terms) {
        MultiIndex index=multi_index_from_object(item.first);
        if(index.size()!=result.argument_size()) {
            throw std::invalid_argument("All Expansion indices must have the same argument size.");
        }
        result.append(index,from_python_object_or_literal<X>(item.second));
    }
    return result;
}

template<class X>
pybind11::class_<Expansion<MultiIndex,X>>
export_expansion(pybind11::module& module) {
    using E=Expansion<MultiIndex,X>;
    using PR=typename X::PrecisionType;

    const std::string name=python_template_class_name<X>("Expansion");
    pybind11::class_<E> expansion_class(module,name.c_str());

    expansion_class.def(pybind11::init<SizeType,PR>());
    expansion_class.def(pybind11::init(&expansion_from_python<X>));

    expansion_class.def("argument_size",&E::argument_size);
    expansion_class.def("number_of_terms",&E::number_of_terms);
    expansion_class.def("graded_sort",&E::graded_sort);
    expansion_class.def("__getitem__",[](E const& e, MultiIndex const& a) { return X(e[a]); });
    expansion_class.def("__setitem__",[](E& e, MultiIndex const& a, X const& x) { e.set(a,x); });
    expansion_class.def("__str__",&__cstr__<E>);
    expansion_class.def("__repr__",[](E const& e) {
        StringStream ss;
        ss << python_representation(e);
        return ss.str();
    });

    return expansion_class;
}

template<class X>
pybind11::class_<Series<X>>
export_series(pybind11::module& module) {
    using S=Series<X>;

    const std::string name=python_template_class_name<X>("Series");
    pybind11::class_<S> series_class(module,name.c_str());

    series_class.def_static("reciprocal",[](X const& c) { return S(Rec(),c); });
    series_class.def_static("sqrt",[](X const& c) { return S(Sqrt(),c); });
    series_class.def_static("exp",[](X const& c) { return S(Exp(),c); });
    series_class.def_static("log",[](X const& c) { return S(Log(),c); });
    series_class.def_static("sin",[](X const& c) { return S(Sin(),c); });
    series_class.def_static("cos",[](X const& c) { return S(Cos(),c); });
    series_class.def_static("tan",[](X const& c) { return S(Tan(),c); });
    series_class.def_static("tanh",[](X const& c) { return S(Tanh(),c); });
    series_class.def_static("atan",[](X const& c) { return S(Atan(),c); });

    series_class.def("__getitem__",[](S const& s, DegreeType d) { return X(s[d]); });
    series_class.def("coefficients",[](S const& s, DegreeType d) {
        pybind11::list result;
        for(DegreeType i=0; i<=d; ++i) {
            result.append(pybind11::cast(s[i]));
        }
        return result;
    });
    series_class.def("__str__",&__cstr__<S>);
    series_class.def("__repr__",&__cstr__<S>);

    return series_class;
}

} // namespace

void algebra_submodule(pybind11::module& module) {
    pybind11::class_<MultiIndex> multi_index_class(module,"MultiIndex");
    multi_index_class.def(pybind11::init<SizeType>());
    multi_index_class.def(pybind11::init(&multi_index_from_python));
    multi_index_class.def_static("unit",&MultiIndex::unit);
    multi_index_class.def("__len__",&MultiIndex::size);
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def("__str__",&__cstr__<MultiIndex>);
    multi_index_class.def("__repr__",[](MultiIndex const& a) {
        StringStream ss;
        ss << "MultiIndex((";
        for(SizeType i=0; i!=a.size(); ++i) {
            if(i!=0) { ss << ","; }
            ss << Int(a[i]);
        }
        if(a.size()==1u) { ss << ","; }
        ss << "))";
        return ss.str();
    });

    export_expansion<FloatDPApproximation>(module);
    export_expansion<FloatDPBounds>(module);
    export_expansion<FloatMPApproximation>(module);
    export_expansion<FloatMPBounds>(module);

    export_series<FloatDPApproximation>(module);
    export_series<FloatDPBounds>(module);
    export_series<FloatMPApproximation>(module);
    export_series<FloatMPBounds>(module);
}
