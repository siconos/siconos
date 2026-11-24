/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include "Relay_options.h"

namespace py = pybind11;

void wrap_relay(py::module_& m, py::module_& params, py::module_& solver_ids) {
  // RELAY_SOLVER enum

  py::enum_<RELAY_SOLVER>(solver_ids, "RELAY_SOLVER", "Relay solvers enum")
      .value("SICONOS_RELAY_PGS", RELAY_SOLVER::SICONOS_RELAY_PGS, "PGS Relay solver")
      .value("SICONOS_RELAY_ENUM", RELAY_SOLVER::SICONOS_RELAY_ENUM, "Enum Relay solver")
      .value("SICONOS_RELAY_PATH", RELAY_SOLVER::SICONOS_RELAY_PATH, "Path Relay solver")
      .value("SICONOS_RELAY_LEMKE", RELAY_SOLVER::SICONOS_RELAY_LEMKE, "Lemke Relay solver")
      .value("SICONOS_RELAY_AVI_CAOFERRIS", RELAY_SOLVER::SICONOS_RELAY_AVI_CAOFERRIS,
             "AVI CaoFerris Relay solver")
      .value("SICONOS_RELAY_AVI_CAOFERRIS_TEST",
             RELAY_SOLVER::SICONOS_RELAY_AVI_CAOFERRIS_TEST,
             "Test for AVI CaoFerris Relay solver")
      .export_values();
}
