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

#include "lcp_cst.h"

namespace py = pybind11;

void wrap_lcp(py::module_ &m, py::module_ &params, py::module_ &solver_ids) {
  // LCP_SOLVER enum
  py::enum_<LCP_SOLVER>(solver_ids, "LCP_SOLVER", "LCP solvers enum")
      .value("SICONOS_LCP_LEMKE", LCP_SOLVER::SICONOS_LCP_LEMKE, "Lemke LCP solver")
      .value("SICONOS_LCP_NSGS_SBM", LCP_SOLVER::SICONOS_LCP_NSGS_SBM, "NSGS SBM LCP solver")
      .value("SICONOS_LCP_PGS", LCP_SOLVER::SICONOS_LCP_PGS, "PGS LCP solver")
      .value("SICONOS_LCP_CPG", LCP_SOLVER::SICONOS_LCP_CPG, "CPG LCP solver")
      .value("SICONOS_LCP_LATIN", LCP_SOLVER::SICONOS_LCP_LATIN, "Latin LCP solver")
      .value("SICONOS_LCP_LATIN_W", LCP_SOLVER::SICONOS_LCP_LATIN_W, "Latin W LCP solver")
      .value("SICONOS_LCP_QP", LCP_SOLVER::SICONOS_LCP_QP, "QP LCP solver")
      .value("SICONOS_LCP_NSQP", LCP_SOLVER::SICONOS_LCP_NSQP, "NSQP LCP solver")
      .value("SICONOS_LCP_NEWTONMIN", LCP_SOLVER::SICONOS_LCP_NEWTONMIN,
             "Newtonmin LCP solver")
      .value("SICONOS_LCP_NEWTON_FB_FBLSA", LCP_SOLVER::SICONOS_LCP_NEWTON_FB_FBLSA,
             "Newton FB FBLSA LCP solver")
      .value("SICONOS_LCP_PSOR", LCP_SOLVER::SICONOS_LCP_PSOR, "PSOR LCP solver")
      .value("SICONOS_LCP_RPGS", LCP_SOLVER::SICONOS_LCP_RPGS, "RPGS LCP solver")
      .value("SICONOS_LCP_PATH", LCP_SOLVER::SICONOS_LCP_PATH, "Path LCP solver")
      .value("SICONOS_LCP_ENUM", LCP_SOLVER::SICONOS_LCP_ENUM, "Enum LCP solver")
      .value("SICONOS_LCP_AVI_CAOFERRIS", LCP_SOLVER::SICONOS_LCP_AVI_CAOFERRIS,
             "AVI CaoFerris LCP solver")
      .value("SICONOS_LCP_PIVOT", LCP_SOLVER::SICONOS_LCP_PIVOT, "Pivot LCP solver")
      .value("SICONOS_LCP_BARD", LCP_SOLVER::SICONOS_LCP_BARD, "Bard LCP solver")
      .value("SICONOS_LCP_MURTY", LCP_SOLVER::SICONOS_LCP_MURTY, "Murty LCP solver")
      .value("SICONOS_LCP_NEWTON_MIN_FBLSA", LCP_SOLVER::SICONOS_LCP_NEWTON_MIN_FBLSA,
             "Newton Min FBLSA LCP solver")
      .value("SICONOS_LCP_PATHSEARCH", LCP_SOLVER::SICONOS_LCP_PATHSEARCH,
             "Pathsearch LCP solver")
      .value("SICONOS_LCP_PIVOT_LUMOD", LCP_SOLVER::SICONOS_LCP_PIVOT_LUMOD,
             "Pivot Lumod LCP solver")
      .value("SICONOS_LCP_GAMS", LCP_SOLVER::SICONOS_LCP_GAMS, "GAMS LCP solver")
      .value("SICONOS_LCP_CONVEXQP_PG", LCP_SOLVER::SICONOS_LCP_CONVEXQP_PG,
             "ConvexQP PG LCP solver")
      .export_values();

  py::enum_<SICONOS_LCP_IPARAM>(params, "SICONOS_LCP_IPARAM_enum", "LCP IPARAM enum")
      .value("SICONOS_LCP_IPARAM_NSGS_ITERATIONS_SUM",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_NSGS_ITERATIONS_SUM,
             "Sum of local solver iterations")
      .value("SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE,
             "Type of pivoting methods")
      .value("SICONOS_LCP_IPARAM_SKIP_TRIVIAL",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_SKIP_TRIVIAL, "Skip trivial solution")
      .value("SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS,
             "Number of possible solutions")
      .value("SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM, "Current enum")
      .value("SICONOS_LCP_IPARAM_ENUM_SEED", SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_ENUM_SEED,
             "Seed for enum start")
      .value("SICONOS_LCP_IPARAM_ENUM_USE_DGELS",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_ENUM_USE_DGELS,
             "Use DGELS instead of DGESV")
      .value("SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS",
             SICONOS_LCP_IPARAM::SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS,
             "Activate multiple solutions search")
      .export_values();

  py::enum_<SICONOS_LCP_DPARAM>(params, "SICONOS_LCP_DPARAM_enum", "LCP DPARAM enum")
      .value("SICONOS_LCP_DPARAM_RHO", SICONOS_LCP_DPARAM::SICONOS_LCP_DPARAM_RHO,
             "Relaxation or regularization parameter")
      .value("SICONOS_LCP_DPARAM_NSGS_LOCAL_ERROR_SUM",
             SICONOS_LCP_DPARAM::SICONOS_LCP_DPARAM_NSGS_LOCAL_ERROR_SUM,
             "Sum of local error values")
      .value("SICONOS_LCP_DPARAM_LATIN_PARAMETER",
             SICONOS_LCP_DPARAM::SICONOS_LCP_DPARAM_LATIN_PARAMETER, "Latin parameter")
      .export_values();

  py::enum_<SICONOS_LCP_SKIP_TRIVIAL>(params, "SICONOS_LCP_SKIP_TRIVIAL_enum",
                                      "LCP skip trivial enum")
      .value("SICONOS_LCP_SKIP_TRIVIAL_NO",
             SICONOS_LCP_SKIP_TRIVIAL::SICONOS_LCP_SKIP_TRIVIAL_NO, "No trivial solution skip")
      .value("SICONOS_LCP_SKIP_TRIVIAL_YES",
             SICONOS_LCP_SKIP_TRIVIAL::SICONOS_LCP_SKIP_TRIVIAL_YES,
             "Yes trivial solution skip")
      .export_values();

  py::enum_<SICONOS_LCP_PIVOT_TYPE>(params, "SICONOS_LCP_PIVOT_TYPE_enum", "LCP pivot type enum")
      .value("SICONOS_LCP_PIVOT_BARD", SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_BARD,
             "Bard pivoting")
      .value("SICONOS_LCP_PIVOT_LEAST_INDEX",
             SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_LEAST_INDEX, "Least index pivoting")
      .value("SICONOS_LCP_PIVOT_LEMKE", SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_LEMKE,
             "Lemke pivoting")
      .value("SICONOS_LCP_PIVOT_PATHSEARCH",
             SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_PATHSEARCH, "Pathsearch pivoting")
      .export_values();
}
