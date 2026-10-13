#include "pybind11.hpp"

void algebra_submodule(pybind11::module& module);
void linear_algebra_submodule(pybind11::module& module);
void polynomial_submodule(pybind11::module& module);
void differentiation_submodule(pybind11::module& module);

void algebra_module(pybind11::module& module) {
    algebra_submodule(module);
    linear_algebra_submodule(module);
    polynomial_submodule(module);
    differentiation_submodule(module);
}
