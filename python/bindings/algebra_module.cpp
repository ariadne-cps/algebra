#include "pybind11.hpp"

void foundation_submodule(pybind11::module& module);
void numeric_submodule(pybind11::module& module);
void interval_submodule(pybind11::module& module);
void algebra_submodule(pybind11::module& module);
void linear_algebra_submodule(pybind11::module& module);
void differentiation_submodule(pybind11::module& module);

PYBIND11_MODULE(pyariadne, module) {
    foundation_submodule(module);
    numeric_submodule(module);
    interval_submodule(module);
    algebra_submodule(module);
    linear_algebra_submodule(module);
    differentiation_submodule(module);
}
