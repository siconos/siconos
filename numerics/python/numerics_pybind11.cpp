/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "NM_types.h"
#include "NumericsVerbose.h"
#include "SolverOptions.h"

namespace py = pybind11;

// Forward declarations

// Friction contact problems wrappers
void wrap_friction_contact(py::module_ &m, py::module_ &params, py::module_ &solver_ids);
// LCP
void wrap_lcp(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

void wrap_mlcp(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

void wrap_relay(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

void wrap_generic_mechanical(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

py::array_t<int> get_iparam(SolverOptions &options) {
  return py::array_t<int>({options.iSize}, {sizeof(int)}, options.iparam, py::cast(&options));
}

py::array_t<double> get_dparam(SolverOptions &options) {
  return py::array_t<double>({options.dSize}, {sizeof(double)}, options.dparam,
                             py::cast(&options));
}

PYBIND11_MODULE(_numerics, m) {
  m.doc() = "Siconos numerics - Nonsmooth problems solvers toolbox";

  py::class_<SolverOptions, std::shared_ptr<SolverOptions>>(m, "SolverOptions")
      .def(py::init<>())
      .def_readwrite("solverId", &SolverOptions::solverId)
      .def_readonly("iSize", &SolverOptions::iSize)
      .def_readonly("dSize", &SolverOptions::dSize)
      .def_readwrite("isSet", &SolverOptions::isSet)
      .def_property_readonly("iparam", &get_iparam)
      .def_property_readonly("dparam", &get_dparam)
      .def_readwrite("filterOn", &SolverOptions::filterOn)
      .def("print", [](SolverOptions &options) { solver_options_print(&options); })
      .def("__repr__", [](const SolverOptions &options) {
        std::ostringstream oss;
        oss << "SolverOptions (ID: " << options.solverId << ")";
        solver_options_print(const_cast<SolverOptions *>(&options));
        return oss.str();
      });

  m.def("solver_options_create", &solver_options_create,
        py::return_value_policy::take_ownership, py::arg("solverId"));
  m.def("solver_options_get_internal_solver", &solver_options_get_internal_solver,
        py::return_value_policy::take_ownership, py::arg("options"), py::arg("id"));
  m.def("numerics_set_verbose", &numerics_set_verbose);
  m.def("solver_options_id_to_name", &solver_options_id_to_name);
  m.def("solver_options_update_internal", &solver_options_update_internal, py::arg("options"),
        py::arg("internal_solver_number"), py::arg("solver_id"));

  py::module_ params = m.def_submodule(
      "params", "Parameter names in numerics (storage types, param for solvers ...)");

  py::enum_<NumericsMatrix_types>(params, "NumericsMatrix_types",
                                  "Types of storage for NumericsMatrix")
      .value("NM_DENSE", NM_DENSE, "Dense format")
      .value("NM_SPARSE_BLOCK", NM_SPARSE_BLOCK, "Sparse block format")
      .value("NM_SPARSE", NM_SPARSE, "Compressed column format")
      .value("NM_UNKNOWN", NM_UNKNOWN, "Unset. Used in NM_null")
      .export_values();

  py::enum_<SICONOS_IPARAM>(params, "SICONOS_IPARAM_enum", "Some value for iparam index")
      .value("SICONOS_IPARAM_MAX_ITER", SICONOS_IPARAM::SICONOS_IPARAM_MAX_ITER,
             "Maximum iterations")
      .value("SICONOS_IPARAM_ITER_DONE", SICONOS_IPARAM::SICONOS_IPARAM_ITER_DONE,
             "Iterations done")
      .value("SICONOS_IPARAM_PREALLOC", SICONOS_IPARAM::SICONOS_IPARAM_PREALLOC,
             "Preallocate memory")
      .value("SICONOS_IPARAM_NSGS_SHUFFLE", SICONOS_IPARAM::SICONOS_IPARAM_NSGS_SHUFFLE,
             "Shuffle for NSGS")
      .value("SICONOS_IPARAM_ERROR_EVALUATION",
             SICONOS_IPARAM::SICONOS_IPARAM_ERROR_EVALUATION, "Error evaluation")
      .value("SICONOS_IPARAM_PATHSEARCH_STACKSIZE",
             SICONOS_IPARAM::SICONOS_IPARAM_PATHSEARCH_STACKSIZE, "Path search stack size")
      .export_values();

  py::enum_<SICONOS_DPARAM>(params, "SICONOS_DPARAM_enum", "Some values for double parameter index")
      .value("SICONOS_DPARAM_TOL", SICONOS_DPARAM::SICONOS_DPARAM_TOL, "Tolerance parameter")
      .value("SICONOS_DPARAM_RESIDU", SICONOS_DPARAM::SICONOS_DPARAM_RESIDU,
             "Residual parameter")
      .export_values();

  py::module_ solver_ids = m.def_submodule(
      "solver_ids", "List of ids for all registered solvers in siconos.numerics");

  wrap_friction_contact(m, params, solver_ids);

  wrap_lcp(m, params, solver_ids);

  wrap_mlcp(m, params, solver_ids);

  wrap_relay(m, params, solver_ids);

  wrap_generic_mechanical(m, params, solver_ids);
}
