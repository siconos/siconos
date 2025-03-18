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

#include "mlcp_cst.h"

namespace py = pybind11;

void wrap_mlcp(py::module_ &m, py::module_ &params, py::module_ &solver_ids) {
  // MLCP solvers

  // Exposing the SICONOS_IPARAM_MLCP enum for integer parameters in MLCP solvers
  py::enum_<SICONOS_IPARAM_MLCP>(params, "SICONOS_IPARAM_MLCP",
                                 "Integer parameters for MLCP solvers")
      .value("SICONOS_IPARAM_MLCP_PGS_EXPLICIT",
             SICONOS_IPARAM_MLCP::SICONOS_IPARAM_MLCP_PGS_EXPLICIT, "Explicit PGS")
      .value("SICONOS_IPARAM_MLCP_PGS_SUM_ITER",
             SICONOS_IPARAM_MLCP::SICONOS_IPARAM_MLCP_PGS_SUM_ITER, "Sum Iterations for PGS")
      .value("SICONOS_IPARAM_MLCP_ENUM_USE_DGELS",
             SICONOS_IPARAM_MLCP::SICONOS_IPARAM_MLCP_ENUM_USE_DGELS,
             "Activate to use dgels rather than dgesv in MLCP driver")
      .value("SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS",
             SICONOS_IPARAM_MLCP::SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS,
             "Number of possible configurations")
      .value("SICONOS_IPARAM_MLCP_UPDATE_REQUIRED",
             SICONOS_IPARAM_MLCP::SICONOS_IPARAM_MLCP_UPDATE_REQUIRED,
             "True if the problem needs update")
      .export_values();

  // Exposing the SICONOS_DPARAM_MLCP enum for double precision parameters in MLCP solvers
  py::enum_<SICONOS_DPARAM_MLCP>(params, "SICONOS_DPARAM_MLCP",
                                 "Double precision parameters for MLCP solvers")
      .value("SICONOS_DPARAM_MLCP_PGS_SUM_ERRORS",
             SICONOS_DPARAM_MLCP::SICONOS_DPARAM_MLCP_PGS_SUM_ERRORS, "Sum of errors for PGS")
      .value("SICONOS_DPARAM_MLCP_RHO", SICONOS_DPARAM_MLCP::SICONOS_DPARAM_MLCP_RHO,
             "Rho parameter")
      .value("SICONOS_DPARAM_MLCP_OMEGA", SICONOS_DPARAM_MLCP::SICONOS_DPARAM_MLCP_OMEGA,
             "Omega parameter")
      .value("SICONOS_DPARAM_MLCP_SIGN_TOL_NEG",
             SICONOS_DPARAM_MLCP::SICONOS_DPARAM_MLCP_SIGN_TOL_NEG,
             "Negative tolerance for the direct solver")
      .value("SICONOS_DPARAM_MLCP_SIGN_TOL_POS",
             SICONOS_DPARAM_MLCP::SICONOS_DPARAM_MLCP_SIGN_TOL_POS,
             "Positive tolerance for the direct solver")
      .export_values();

  py::enum_<MLCP_SOLVER>(solver_ids, "MLCP_SOLVER",
                         "Mixed Linear Complementarity Problem (MLCP) solvers")
      .value("SICONOS_MLCP_PGS", MLCP_SOLVER::SICONOS_MLCP_PGS,
             "Projected Gauss-Seidel solver")
      .value("SICONOS_MLCP_RPGS", MLCP_SOLVER::SICONOS_MLCP_RPGS,
             "Regularized Projected Gauss-Seidel solver")
      .value("SICONOS_MLCP_PSOR", MLCP_SOLVER::SICONOS_MLCP_PSOR,
             "Preconditioned Successive Over-Relaxation solver")
      .value("SICONOS_MLCP_RPSOR", MLCP_SOLVER::SICONOS_MLCP_RPSOR,
             "Regularized Preconditioned Successive Over-Relaxation solver")
      .value("SICONOS_MLCP_PATH", MLCP_SOLVER::SICONOS_MLCP_PATH, "PATH solver")
      .value("SICONOS_MLCP_ENUM", MLCP_SOLVER::SICONOS_MLCP_ENUM, "Enum solver")
      .value("SICONOS_MLCP_SIMPLEX", MLCP_SOLVER::SICONOS_MLCP_SIMPLEX, "Simplex solver")
      .value("SICONOS_MLCP_DIRECT_ENUM", MLCP_SOLVER::SICONOS_MLCP_DIRECT_ENUM,
             "Direct Enum solver")
      .value("SICONOS_MLCP_PATH_ENUM", MLCP_SOLVER::SICONOS_MLCP_PATH_ENUM, "PATH Enum solver")
      .value("SICONOS_MLCP_DIRECT_SIMPLEX", MLCP_SOLVER::SICONOS_MLCP_DIRECT_SIMPLEX,
             "Direct Simplex solver")
      .value("SICONOS_MLCP_DIRECT_PATH", MLCP_SOLVER::SICONOS_MLCP_DIRECT_PATH,
             "Direct PATH solver")
      .value("SICONOS_MLCP_DIRECT_PATH_ENUM", MLCP_SOLVER::SICONOS_MLCP_DIRECT_PATH_ENUM,
             "Direct PATH Enum solver")
      .value("SICONOS_MLCP_FB", MLCP_SOLVER::SICONOS_MLCP_FB, "Fixed-Point solver")
      .value("SICONOS_MLCP_DIRECT_FB", MLCP_SOLVER::SICONOS_MLCP_DIRECT_FB,
             "Direct Fixed-Point solver")
      .value("SICONOS_MLCP_PGS_SBM", MLCP_SOLVER::SICONOS_MLCP_PGS_SBM,
             "Projected Gauss-Seidel with Successive Block Minimization solver")
      .value("SICONOS_MLCP_LCP_LEMKE", MLCP_SOLVER::SICONOS_MLCP_LCP_LEMKE,
             "Lemke's algorithm solver for LCP")
      .export_values();
}