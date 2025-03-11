#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "FrictionContact.hpp"
#include "GlobalFrictionContact.hpp"
#include "GlobalRollingFrictionContact.hpp"
#include "LCP.hpp"
#include "Relay.hpp"
#include "RollingFrictionContact.hpp"
#include "SolverOptions.h"
#include "relay_cst.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(nonsmooth_formulations, m) {
  m.doc() = "Siconos nonsmooth formulations module";

  py::module_ numerics = py::module_::import("siconos.numerics");
  py::object solver_ids = numerics.attr("solver_ids");

  py::class_<siconos::nonsmooth_formulations::OneStepNSProblem,
             std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>(
      m, "OneStepNSProblem")
      .def("setSolverId", &siconos::nonsmooth_formulations::OneStepNSProblem::setSolverId,
           "Set the solver ID for the problem")
      .def("numericsSolverOptions",
           &siconos::nonsmooth_formulations::OneStepNSProblem::numericsSolverOptions,
           py::return_value_policy::reference_internal, "Access to solver options object")
      .def("setMaxSize", &siconos::nonsmooth_formulations::OneStepNSProblem::setMaxSize,
           "Set the maximum size of the problem")
      .def("setNumericsVerboseMode",
           &siconos::nonsmooth_formulations::OneStepNSProblem::setNumericsVerboseMode,
           "Set the verbose mode for the numerics solver")
      .def("getSizeOutput", &siconos::nonsmooth_formulations::OneStepNSProblem::getSizeOutput,
           "Get the size of the output of the problem");

  py::class_<siconos::nonsmooth_formulations::LinearOSNS,
             std::shared_ptr<siconos::nonsmooth_formulations::LinearOSNS>,
             siconos::nonsmooth_formulations::OneStepNSProblem>(m, "LinearOSNS")
      .def("setMStorageType", &siconos::nonsmooth_formulations::LinearOSNS::setMStorageType)
      .def("setKeepLambdaAndYState",
           &siconos::nonsmooth_formulations::LinearOSNS::setKeepLambdaAndYState);

  py::class_<siconos::nonsmooth_formulations::LCP,
             std::shared_ptr<siconos::nonsmooth_formulations::LCP>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "LCP")
      .def(py::init<int>(),
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_LCP_LEMKE"));

  py::class_<siconos::nonsmooth_formulations::Relay,
             std::shared_ptr<siconos::nonsmooth_formulations::Relay>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "Relay")
      .def(py::init<int>(),
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_RELAY_AVI_CAOFERRIS"));

  py::class_<siconos::nonsmooth_formulations::FrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::FrictionContact>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "FrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb") = 3,
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_FRICTION_3D_NSGS"))
      .def(py::init<int, std::shared_ptr<SolverOptions>>(), py::arg("dimPb"),
           py::arg("options"));

  py::class_<siconos::nonsmooth_formulations::RollingFrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::RollingFrictionContact>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "RollingFrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb") = 3,
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_ROLLING_FRICTION_3D_NSGS"))
      .def(py::init<int, std::shared_ptr<SolverOptions>>(), py::arg("dimPb"),
           py::arg("options"));

  py::class_<siconos::nonsmooth_formulations::GlobalFrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::GlobalFrictionContact>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "GlobalFrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb"),
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_GLOBAL_FRICTION_3D_NSGS"));

  py::class_<siconos::nonsmooth_formulations::GlobalRollingFrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::GlobalRollingFrictionContact>,
             siconos::nonsmooth_formulations::GlobalFrictionContact>(
      m, "GlobalRollingFrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb"),
           py::arg("numericsSolverId") =
               solver_ids.attr("SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR"));
}