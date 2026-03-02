#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "FrictionContact.hpp"
// #include "FrictionContactProblem.h"
#include "GenericMechanical.hpp"
#include "GlobalFrictionContact.hpp"
#include "GlobalFrictionContactProblem.h"
#include "GlobalRollingFrictionContact.hpp"
#include "GlobalRollingFrictionContactProblem.h"
#include "LCP.hpp"
#include "LinearOSNS.hpp"
#include "MLCP.hpp"
#include "MLCPProjectOnConstraints.hpp"
#include "Relay.hpp"
#include "RollingFrictionContact.hpp"
#include "RollingFrictionContactProblem.h"
#include "SolverOptions.h"
// #include "Relay_options.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(nonsmooth_formulations, m) {
  m.doc() = "Siconos nonsmooth formulations module";

  py::module_ numerics = py::module_::import("siconos.numerics");
  py::object solver_ids = numerics.attr("solver_ids");

  py::enum_<siconos::nonsmooth_formulations::LinearOSNSAssemblyType>(m,
                                                                     "LinearOSNSAssemblyType")
      .value("REDUCED_BLOCK",
             siconos::nonsmooth_formulations::LinearOSNSAssemblyType::REDUCED_BLOCK,
             "Reduced block formulation")
      .value("GLOBAL", siconos::nonsmooth_formulations::LinearOSNSAssemblyType::GLOBAL,
             "Global formulation")
      .value("REDUCED_DIRECT",
             siconos::nonsmooth_formulations::LinearOSNSAssemblyType::REDUCED_DIRECT,
             "Reduced direct formulation")
      .value("GLOBAL_REDUCED",
             siconos::nonsmooth_formulations::LinearOSNSAssemblyType::GLOBAL_REDUCED,
             "Global reduced formulation")
      .export_values();

  py::class_<siconos::nonsmooth_formulations::OneStepNSProblem,
             std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>(
      m, "OneStepNSProblem")
      .def("setSolverId", &siconos::nonsmooth_formulations::OneStepNSProblem::setSolverId,
           "Set the solver ID for the problem")
      .def("numericsSolverOptions",
           &siconos::nonsmooth_formulations::OneStepNSProblem::numericsSolverOptions,
           py::return_value_policy::reference_internal, "Access to solver options object")
      .def("setMaxSize", &siconos::nonsmooth_formulations::OneStepNSProblem::setMaxSize,
           "Set the maximum size of the nonsmooth problem")
      .def("setNumericsVerboseMode",
           &siconos::nonsmooth_formulations::OneStepNSProblem::setNumericsVerboseMode,
           "Set the verbose mode for the numerics solver")
      .def("postCompute", &siconos::nonsmooth_formulations::OneStepNSProblem::postCompute,
           "post treatment for output of the solver")
      .def_property_readonly("getSizeOutput",
                             &siconos::nonsmooth_formulations::OneStepNSProblem::getSizeOutput,
                             "Get the size of the output of the problem");

  py::class_<siconos::nonsmooth_formulations::LinearOSNS,
             std::shared_ptr<siconos::nonsmooth_formulations::LinearOSNS>,
             siconos::nonsmooth_formulations::OneStepNSProblem>(m, "LinearOSNS")
      .def("setMStorageType", &siconos::nonsmooth_formulations::LinearOSNS::setMStorageType)
      .def("setKeepLambdaAndYState",
           &siconos::nonsmooth_formulations::LinearOSNS::setKeepLambdaAndYState)
      .def("setAssemblyType", &siconos::nonsmooth_formulations::LinearOSNS::setAssemblyType)
      .def("preCompute", &siconos::nonsmooth_formulations::LinearOSNS::preCompute,
           py::arg("time"), "build problem coefficients (if required)");

  py::class_<siconos::nonsmooth_formulations::LCP,
             std::shared_ptr<siconos::nonsmooth_formulations::LCP>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "LCP")
      .def(py::init<int>(), py::arg("numericsSolverId") = solver_ids.attr("SICONOS_LCP_LEMKE"))
      .def("solve", &siconos::nonsmooth_formulations::LCP::solve);

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
           py::arg("options"))
      .def("updateMu", &siconos::nonsmooth_formulations::FrictionContact::updateMu)
      .def("solve", &siconos::nonsmooth_formulations::FrictionContact::solve);

  py::class_<siconos::nonsmooth_formulations::RollingFrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::RollingFrictionContact>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "RollingFrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb") = 3,
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_ROLLING_FRICTION_3D_NSGS"))
      .def(py::init<int, std::shared_ptr<SolverOptions>>(), py::arg("dimPb"),
           py::arg("options"))
      .def("updateMu", &siconos::nonsmooth_formulations::RollingFrictionContact::updateMu)
      .def("solve", &siconos::nonsmooth_formulations::RollingFrictionContact::solve);

  py::class_<siconos::nonsmooth_formulations::GlobalFrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::GlobalFrictionContact>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "GlobalFrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb"),
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_GLOBAL_FRICTION_3D_NSGS"))
      .def(py::init<int, std::shared_ptr<SolverOptions>>(), py::arg("dimPb"),
           py::arg("options"))
      .def("updateMu", &siconos::nonsmooth_formulations::GlobalFrictionContact::updateMu)
      .def("solve", &siconos::nonsmooth_formulations::GlobalFrictionContact::solve);

  py::class_<siconos::nonsmooth_formulations::GenericMechanical,
             std::shared_ptr<siconos::nonsmooth_formulations::GenericMechanical>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "GenericMechanical")
      .def(py::init<int>(),
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_ONECONE_NSN"))
      .def(py::init<std::shared_ptr<SolverOptions>>(), py::arg("options"));

  py::class_<siconos::nonsmooth_formulations::MLCP,
             std::shared_ptr<siconos::nonsmooth_formulations::MLCP>,
             siconos::nonsmooth_formulations::LinearOSNS>(m, "MLCP")
      .def(py::init<int>(), py::arg("numericsSolverId") = solver_ids.attr("SICONOS_MLCP_ENUM"))
      .def(py::init<std::shared_ptr<SolverOptions>>(), py::arg("options"))
      .def("solve", &siconos::nonsmooth_formulations::MLCP::solve);

  py::class_<siconos::nonsmooth_formulations::MLCPProjectOnConstraints,
             std::shared_ptr<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>,
             siconos::nonsmooth_formulations::MLCP>(m, "MLCPProjectOnConstraints")
      .def(py::init<int, double>(),
           py::arg("numericsSolverId") = solver_ids.attr("SICONOS_MLCP_ENUM"),
           py::arg("alpha") = 1.0)
      .def(py::init<std::shared_ptr<SolverOptions>, double>(), py::arg("options"),
           py::arg("alpha") = 1.0);

  py::class_<siconos::nonsmooth_formulations::GlobalRollingFrictionContact,
             std::shared_ptr<siconos::nonsmooth_formulations::GlobalRollingFrictionContact>,
             siconos::nonsmooth_formulations::GlobalFrictionContact>(
      m, "GlobalRollingFrictionContact")
      .def(py::init<int, int>(), py::arg("dimPb") = 5,
           py::arg("numericsSolverId") =
               solver_ids.attr("SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR"))
      .def(py::init<int, std::shared_ptr<SolverOptions>>(), py::arg("dimPb"),
           py::arg("options"))
      .def("updateMu",
           &siconos::nonsmooth_formulations::GlobalRollingFrictionContact::updateMu)
      .def("solve", &siconos::nonsmooth_formulations::GlobalRollingFrictionContact::solve);
}
