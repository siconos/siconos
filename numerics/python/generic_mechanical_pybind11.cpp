// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include "GenericMechanical_cst.h"

namespace py = pybind11;

void wrap_generic_mechanical(py::module_& m, py::module_& params, py::module_& solver_ids) {
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
             GENERIC_MECHANICAL_IPARAM::SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_RELAXATION,
             "Line search flag")
      .export_values();

  py::enum_<GENERIC_MECHANICAL_DPARAM>(params, "GENERIC_MECHANICAL_DPARAM",
                                       "Generic Mechanical DPARAM enum")
      .value("SICONOS_DPARAM_GMP_ERROR_LS",
             GENERIC_MECHANICAL_DPARAM::SICONOS_DPARAM_GMP_RELAXATION_ERROR,
             "Error threshold for line search")
      .value("SICONOS_DPARAM_GMP_COEFF_LS",
             GENERIC_MECHANICAL_DPARAM::SICONOS_DPARAM_GMP_RELAXATION_COEFF,
             "Coefficient for line search")
      .export_values();

  py::enum_<GENERIC_MECHANICAL_ISREDUCED>(params, "GENERIC_MECHANICAL_ISREDUCED",
                                          "Generic Mechanical ISREDUCED enum")
      .value("SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS,
             "GS on all blocks")
      .value("SICONOS_GENERIC_MECHANICAL_SUBSTITUTE_EQUALITIES",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_SUBSTITUTE_EQUALITIES,
             "Substituted equalities")
      .value("SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES,
             "Assemblated equalities")
      .value("SICONOS_GENERIC_MECHANICAL_MLCP_LIKE",
             GENERIC_MECHANICAL_ISREDUCED::SICONOS_GENERIC_MECHANICAL_MLCP_LIKE,
             "Solve like MLCP")
      .export_values();
}
