#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "LCP.hpp"
#include "Relay.hpp"
#include "relay_cst.h"
#include "SolverOptions.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(nonsmooth_formulations, m) {
  m.doc() = "Siconos nonsmooth formulations module";

  py::class_<siconos::nonsmooth_formulations::OneStepNSProblem,
             std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>(
      m, "OneStepNSProblem")
      .def("setSolverId", &siconos::nonsmooth_formulations::OneStepNSProblem::setSolverId,
           "Set the solver ID for the problem")
      .def("numericsSolverOptions",
           &siconos::nonsmooth_formulations::OneStepNSProblem::numericsSolverOptions,
           py::return_value_policy::reference_internal, "Access to solver options object");

  py::class_<siconos::nonsmooth_formulations::LinearOSNS,
             std::shared_ptr<siconos::nonsmooth_formulations::LinearOSNS>,
             siconos::nonsmooth_formulations::OneStepNSProblem>(m, "LinearOSNS");

  py::class_<siconos::nonsmooth_formulations::LCP,
             std::shared_ptr<siconos::nonsmooth_formulations::LCP>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "LCP")
      .def(py::init<int>(), py::arg("numericsSolverId") = 200);

  py::class_<siconos::nonsmooth_formulations::Relay,
             std::shared_ptr<siconos::nonsmooth_formulations::Relay>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "Relay")
      .def(py::init<int>(), py::arg("numericsSolverId") = 306);
}