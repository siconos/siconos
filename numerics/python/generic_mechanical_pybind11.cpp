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

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include "GenericMechanical_cst.h"

namespace py = pybind11;

void wrap_generic_mechanical(py::module_ &m, py::module_ &params, py::module_ &solver_ids) {
  // GENERIC_MECHANICAL_SOLVER enum

  py::enum_<GENERIC_MECHANICAL_SOLVER>(solver_ids, "GENERIC_MECHANICAL_SOLVER",
                                       "Generic Mechanical solvers enum")
      .value("SICONOS_GENERIC_MECHANICAL_NSGS",
             GENERIC_MECHANICAL_SOLVER::SICONOS_GENERIC_MECHANICAL_NSGS,
             "Generic Mechanical NSGS solver")
      .export_values();

  py::enum_<GENERIC_MECHANICAL_IPARAM>(params, "GENERIC_MECHANICAL_IPARAM",
                                       "Generic Mechanical IPARAM enum")
      .value("SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED",
             GENERIC_MECHANICAL_IPARAM::SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED,
             "Reduced mode flag")
      .value("SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH",
             GENERIC_MECHANICAL_IPARAM::SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH,
             "Line search flag")
      .export_values();

  py::enum_<GENERIC_MECHANICAL_DPARAM>(params, "GENERIC_MECHANICAL_DPARAM",
                                       "Generic Mechanical DPARAM enum")
      .value("SICONOS_DPARAM_GMP_ERROR_LS",
             GENERIC_MECHANICAL_DPARAM::SICONOS_DPARAM_GMP_ERROR_LS,
             "Error threshold for line search")
      .value("SICONOS_DPARAM_GMP_COEFF_LS",
             GENERIC_MECHANICAL_DPARAM::SICONOS_DPARAM_GMP_COEFF_LS,
             "Coefficient for line search")
      .export_values();

  py::enum_<GENERIC_MECHANICAL_ISREDUCED>(params, "GENERIC_MECHANICAL_ISREDUCED",
                                          "Generic Mechanical ISREDUCED enum")
      .value("SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS,
             "GS on all blocks")
      .value("SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES,
             "Substituted equalities")
      .value("SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES,
             "Assemblated equalities")
      .value("SICONOS_GENERIC_MECHANICAL_MLCP_LIKE",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_MLCP_LIKE,
             "Solve like MLCP")
      .export_values();
}