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
#include "numerics_errors.h"

namespace py = pybind11;

/**
 * @brief Exception class for numerics solver errors
 */
class NumericsException : public std::exception {
 private:
  int error_code_;
  std::string message_;

 public:
  NumericsException(int code, const std::string &msg) : error_code_(code), message_(msg) {}

  const char *what() const noexcept override { return message_.c_str(); }

  int error_code() const { return error_code_; }
};

/**
 * @brief Check solver return code and throw Python exception on error
 *
 * This function should be called after solver operations to convert
 * C error codes into Python exceptions.
 *
 * @param info The error code returned by the solver
 * @param operation Name of the operation (for error message)
 */
inline void check_solver_error(int info, const char *operation) {
  if (info != 0) {
    std::string msg = std::string("Numerics error in ") + operation + ": " +
                      numerics_error_string((NumericsError)info);
    throw NumericsException(info, msg);
  }
}

// Forward declarations

// Friction contact problems wrappers
void wrap_friction_contact(py::module_ &m, py::module_ &params, py::module_ &solver_ids);
// LCP
void wrap_lcp(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

void wrap_mlcp(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

void wrap_relay(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

void wrap_generic_mechanical(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

// Plasticity
void wrap_plasticity(py::module_ &m, py::module_ &params, py::module_ &solver_ids);

py::array_t<int> get_iparam(SolverOptions &options) {
  return py::array_t<int>({options.iSize}, {sizeof(int)}, options.iparam, py::cast(&options));
}

py::array_t<double> get_dparam(SolverOptions &options) {
  return py::array_t<double>({options.dSize}, {sizeof(double)}, options.dparam,
                             py::cast(&options));
}

PYBIND11_MODULE(_numerics, m) {
  m.doc() = "Siconos numerics - Nonsmooth problems solvers toolbox";

  // Register custom exception translator
  py::register_exception<NumericsException>(m, "NumericsError");

  // Export error codes to Python
  py::enum_<::NumericsError>(m, "ErrorCode", "Numerics error codes")
      .value("OK", NUMERICS_OK, "Success")
      .value("NULL_POINTER", NUMERICS_ERR_NULL_POINTER, "Null pointer error")
      .value("INVALID_ARGUMENT", NUMERICS_ERR_INVALID_ARGUMENT, "Invalid argument")
      .value("INVALID_OPTION", NUMERICS_ERR_INVALID_OPTION, "Invalid option")
      .value("INVALID_SOLVER", NUMERICS_ERR_INVALID_SOLVER, "Invalid solver")
      .value("NOT_IMPLEMENTED", NUMERICS_ERR_NOT_IMPLEMENTED, "Not implemented")
      .value("MEMORY", NUMERICS_ERR_MEMORY, "Memory allocation failed")
      .value("FILE_IO", NUMERICS_ERR_FILE_IO, "File I/O error")
      .value("MAX_ITER", NUMERICS_ERR_MAX_ITER, "Maximum iterations reached")
      .value("DIVERGENCE", NUMERICS_ERR_DIVERGENCE, "Solver diverged")
      .value("STAGNATION", NUMERICS_ERR_STAGNATION, "Solver stagnated")
      .value("LINEAR_SOLVER", NUMERICS_ERR_LINEAR_SOLVER, "Linear solver failed")
      .value("LOCAL_SOLVER", NUMERICS_ERR_LOCAL_SOLVER, "Local solver failed")
      .value("PROJECTION", NUMERICS_ERR_PROJECTION, "Projection failed")
      .value("INFEASIBLE", NUMERICS_ERR_INFEASIBLE, "Problem infeasible")
      .value("UNBOUNDED", NUMERICS_ERR_UNBOUNDED, "Problem unbounded")
      .value("ILL_CONDITIONED", NUMERICS_ERR_ILL_CONDITIONED, "Ill-conditioned problem")
      .value("SINGULAR_MATRIX", NUMERICS_ERR_SINGULAR_MATRIX, "Singular matrix")
      .value("UNKNOWN", NUMERICS_ERR_UNKNOWN, "Unknown error")
      .export_values();

  // Function to get error string
  m.def(
      "error_string", [](int code) { return numerics_error_string((::NumericsError)code); },
      "Get error message for error code", py::arg("error_code"));

  py::class_<SolverOptions, std::shared_ptr<SolverOptions>>(m, "SolverOptions")
      .def(py::init([](int solver_id) {
             SolverOptions *opts = solver_options_create(solver_id);

             if (!opts) throw std::runtime_error("Error during Solver options creation");

             return std::shared_ptr<SolverOptions>(opts, [](SolverOptions *p) {
               if (p) solver_options_delete(p);
             });
           }),
           py::arg("solver_id"))

      .def("__enter__", [](std::shared_ptr<SolverOptions> self) { return self; })
      .def("__exit__",
           [](std::shared_ptr<SolverOptions>, py::object, py::object, py::object) {})
      .def_readwrite("solverId", &SolverOptions::solverId)
      .def_readonly("iSize", &SolverOptions::iSize)
      .def_readonly("dSize", &SolverOptions::dSize)
      .def_readwrite("isSet", &SolverOptions::isSet)
      .def_property_readonly("iparam", &get_iparam)
      .def_property_readonly("dparam", &get_dparam)
      .def_readwrite("filterOn", &SolverOptions::filterOn)
      .def("print", [](SolverOptions *self) { solver_options_print(self); })
      .def(
          "get_internal_solver",
          [](SolverOptions *self, size_t num) {
            return solver_options_get_internal_solver(self, num);
          },
          py::return_value_policy::reference_internal)
      .def("update_internal",
           [](SolverOptions *self, size_t num, int id) {
             return solver_options_update_internal(self, num, id);
           })
      .def_property_readonly(
          "name",
          [](SolverOptions *self) { return solver_options_id_to_name(self->solverId); })

      .def("__repr__", [](const SolverOptions *self) {
        std::ostringstream oss;
        oss << "SolverOptions (ID: " << self->solverId << ")";
        solver_options_print(const_cast<SolverOptions *>(self));
        return oss.str();
      });

  m.def("numerics_set_verbose", &numerics_set_verbose);
  m.def("solver_options_id_to_name", &solver_options_id_to_name);

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

  py::enum_<SICONOS_DPARAM>(params, "SICONOS_DPARAM_enum",
                            "Some values for double parameter index")
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

  wrap_plasticity(m, params, solver_ids);
}
