#include "pybind11.hpp"

void interval_submodule(pybind11::module& module);
void linear_algebra_submodule(pybind11::module& module);
void differentiation_submodule(pybind11::module& module);

PYBIND11_MODULE(pyariadne, module) {
    interval_submodule(module);
    linear_algebra_submodule(module);
    differentiation_submodule(module);
}
