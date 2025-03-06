// # Siconos is a program dedicated to modeling, simulation and control
// # of non smooth dynamical systems.
// #
// # Copyright 2024 INRIA.
// #
// # Licensed under the Apache License, Version 2.0 (the "License");
// # you may not use this file except in compliance with the License.
// # You may obtain a copy of the License at
// #
// # http://www.apache.org/licenses/LICENSE-2.0
// #
// # Unless required by applicable law or agreed to in writing, software
// # distributed under the License is distributed on an "AS IS" BASIS,
// # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// # See the License for the specific language governing permissions and
// # limitations under the License.
// #
// #

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Friction_cst.h"
#include "GenericMechanical_cst.h"
#include "NM_types.h"
#include "SolverOptions.h"
#include "lcp_cst.h"
#include "mlcp_cst.h"
#include "relay_cst.h"

namespace py = pybind11;

// Getter pour iparam qui retourne un tableau numpy lié directement au pointeur C++
py::array_t<int> get_iparam(SolverOptions &options) {
  return py::array_t<int>(
      {options.iSize},      // Dimension du tableau (taille de iparam)
      {sizeof(int)},        // Stride (décalage entre éléments)
      options.iparam,       // Pointeur vers les données C++ (iparam)
      py::cast(&options));  // Lier le cycle de vie de l'objet Python à l'objet C++
}

// Getter pour dparam qui retourne un tableau numpy lié directement au pointeur C++
py::array_t<double> get_dparam(SolverOptions &options) {
  return py::array_t<double>(
      {options.dSize},      // Dimension du tableau (taille de dparam)
      {sizeof(double)},     // Stride (décalage entre éléments)
      options.dparam,       // Pointeur vers les données C++ (dparam)
      py::cast(&options));  // Lier le cycle de vie de l'objet Python à l'objet C++
}

PYBIND11_MODULE(_numerics, m) {
  py::class_<SolverOptions, std::shared_ptr<SolverOptions>>(m, "SolverOptions")
      // Membres simples
      .def(py::init<>())  // Constructeur par défaut
      .def_readwrite("solverId", &SolverOptions::solverId)
      .def_readonly("iSize", &SolverOptions::iSize)
      .def_readonly("dSize", &SolverOptions::dSize)
      .def_readwrite("isSet", &SolverOptions::isSet)
      .def_property_readonly("iparam", &get_iparam)  // Exposer iparam avec accès direct
      .def_property_readonly("dparam", &get_dparam)  // Exposer dparam avec accès direct
      .def_readwrite("filterOn", &SolverOptions::filterOn)
      .def("print",
           [](SolverOptions &options) {
             solver_options_print(&options);  // Utiliser la fonction existante pour afficher
           })
      .def("__repr__", [](const SolverOptions &options) {
        std::ostringstream oss;
        oss << "SolverOptions (ID: " << options.solverId << ")";
        solver_options_print(const_cast<SolverOptions *>(&options));
        return oss.str();
      });

  m.def("solver_options_create", &solver_options_create,
        py::return_value_policy::take_ownership, py::arg("solverId"));

  py::module_ params = m.def_submodule(
      "params", "Parameter names in numerics (storage types, param for solvers ...)");

  py::enum_<NumericsMatrix_types>(params, "NumericsMatrix_types",
                                  "Types of storage for NumericsMatrix")
      .value("NM_DENSE", NM_DENSE, "Dense format")
      .value("NM_SPARSE_BLOCK", NM_SPARSE_BLOCK, "Sparse block format")
      .value("NM_SPARSE", NM_SPARSE, "Compressed column format")
      .value("NM_UNKNOWN", NM_UNKNOWN, "Unset. Used in NM_null")
      .export_values();

  py::enum_<SICONOS_IPARAM>(params, "SICONOS_IPARAM", "Some value for iparam index")
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

  py::enum_<SICONOS_DPARAM>(params, "SICONOS_DPARAM", "Some values for double parameter index")
      .value("SICONOS_DPARAM_TOL", SICONOS_DPARAM::SICONOS_DPARAM_TOL, "Tolerance parameter")
      .value("SICONOS_DPARAM_RESIDU", SICONOS_DPARAM::SICONOS_DPARAM_RESIDU,
             "Residual parameter")
      .export_values();

  py::module_ solver_ids = m.def_submodule(
      "solver_ids", "List of ids for all registered solvers in siconos.numerics");

  // Exposing the FRICTION_SOLVER enum using .value
  py::enum_<FRICTION_SOLVER>(solver_ids, "FRICTION_SOLVER",
                             "Friction solver types for 2D and 3D frictional contact problems")
      // 2D Frictional Contact solvers
      .value("SICONOS_FRICTION_2D_NSGS", FRICTION_SOLVER::SICONOS_FRICTION_2D_NSGS,
             "2D Non-smooth Gauss Seidel")
      .value("SICONOS_FRICTION_2D_CPG", FRICTION_SOLVER::SICONOS_FRICTION_2D_CPG,
             "2D CPG solver")
      .value("SICONOS_FRICTION_2D_LEMKE", FRICTION_SOLVER::SICONOS_FRICTION_2D_LEMKE,
             "2D Lemke solver")
      .value("SICONOS_FRICTION_2D_ENUM", FRICTION_SOLVER::SICONOS_FRICTION_2D_ENUM,
             "2D Friction solver enum")

      // 3D frictional contact solvers on local formulation
      .value("SICONOS_FRICTION_3D_NSGS", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGS,
             "3D Non-smooth Gauss Seidel")
      .value("SICONOS_FRICTION_3D_NSGSV", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGSV,
             "3D Non-smooth Gauss Seidel-velocity")
      .value("SICONOS_FRICTION_3D_PROX", FRICTION_SOLVER::SICONOS_FRICTION_3D_PROX,
             "3D Proximal solver")
      .value("SICONOS_FRICTION_3D_TFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_TFP,
             "3D Tresca solver")
      .value("SICONOS_FRICTION_3D_NSN_AC", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC,
             "3D Non-smooth Newton Alart-Curnier")
      .value("SICONOS_FRICTION_3D_DSFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_DSFP,
             "3D De Saxce fixed point")
      .value("SICONOS_FRICTION_3D_VI_FPP", FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP,
             "3D VI formulation, fixed point projection")
      .value("SICONOS_FRICTION_3D_VI_EG", FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_EG,
             "3D VI formulation, Extra-gradient")
      .value("SICONOS_FRICTION_3D_HP", FRICTION_SOLVER::SICONOS_FRICTION_3D_HP,
             "3D Hyperplane projection")
      .value("SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint,
             "3D Fischer Burmeister fixed point")
      .value("SICONOS_FRICTION_3D_FPP", FRICTION_SOLVER::SICONOS_FRICTION_3D_FPP,
             "3D Fixed point projection")
      .value("SICONOS_FRICTION_3D_EG", FRICTION_SOLVER::SICONOS_FRICTION_3D_EG,
             "3D Extra-gradient")
      .value("SICONOS_FRICTION_3D_NSN_FB", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_FB,
             "3D Non-smooth Newton Fischer Burmeister")
      .value("SICONOS_FRICTION_3D_GAMS_PATH", FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATH,
             "3D GAMS/Path (Ferris)")
      .value("SICONOS_FRICTION_3D_GAMS_PATHVI",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATHVI, "3D GAMS/Path VI formulation")

      // 3D Frictional Contact solvers for one contact (used mainly inside NSGS solvers)
      .value("SICONOS_FRICTION_3D_ONECONTACT_NSN",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN,
             "3D One contact Non-smooth Newton Alart-Curnier")
      .value("SICONOS_FRICTION_3D_ONECONTACT_NSN_GP",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN_GP,
             "3D One contact Non-smooth Newton Alart-Curnier, damped")
      .value("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone,
             "3D One contact Projection on cone")
      .value(
          "SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration",
          FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
          "3D One contact Projection on cone with local iteration")
      .value(
          "SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization",
          FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization,
          "3D One contact Projection on cone with regularization")
      .value(
          "SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization",
          FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization,
          "3D One contact Projection on cone with diagonalization")
      .value("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity,
             "3D One contact Projection on cone, velocity")
      .value("SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID,
             "3D Frictional contact hybrid solver for one contact")
      .value("SICONOS_FRICTION_3D_VI_FPP_Cylinder",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP_Cylinder,
             "3D VI formulation, Fixed Point Projection Cylinder")
      .value("SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER,
             "3D Convex QP Projection on Cylinder")
      .value("SICONOS_FRICTION_3D_ACLMFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_ACLMFP,
             "Alart-Curnier fixed point, local formulation")
      .value("SICONOS_FRICTION_3D_SOCLCP", FRICTION_SOLVER::SICONOS_FRICTION_3D_SOCLCP,
             "Second-order Cone LCP, local formulation")
      .value("SICONOS_FRICTION_3D_GAMS_LCP_PATH",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATH,
             "GAMS/PATH (Ferris) LCP, local formulation")
      .value("SICONOS_FRICTION_3D_GAMS_LCP_PATHVI",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATHVI,
             "VI formulation, GAMS/PATH (Ferris) LCP, local formulation")
      .value("SICONOS_FRICTION_3D_NSN_NM", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_NM,
             "Non-smooth Newton, natural map, local formulation")
      .value("SICONOS_FRICTION_3D_NSN_AC_TEST",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC_TEST, "NSN AC Test solver")
      .value("SICONOS_FRICTION_3D_PFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_PFP,
             "Panagiotopoulos, fixed point, local formulation")
      .value("SICONOS_FRICTION_3D_ADMM", FRICTION_SOLVER::SICONOS_FRICTION_3D_ADMM,
             "ADMM local formulation")

      // Fischer Burmeister/Path, Glocker formulation, one contact solver
      .value("SICONOS_FRICTION_3D_NCPGlockerFBPATH",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBPATH,
             "3D Fischer Burmeister, Glocker formulation")
      .value("SICONOS_FRICTION_3D_NCPGlockerFBNewton",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBNewton,
             "3D Fischer Burmeister, Glocker formulation, Newton")
      .value("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC,
             "3D One contact Quartic solver")
      .value("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU,
             "3D One contact Quartic solver, NU")
      .value("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder,
             "3D One contact Projection on cylinder")
      .value("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration",
             FRICTION_SOLVER::
                 SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration,
             "3D One contact Projection on cylinder with local iteration")

      /** 3D Frictional contact local solvers on global formulation */
      .value("SICONOS_GLOBAL_FRICTION_3D_NSGS_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS_WR,
             "Global 3D Frictional Contact solver NSGS (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR,
             "Global 3D Frictional Contact solver NSGSV (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_PROX_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_PROX_WR,
             "Global 3D Frictional Contact solver PROX (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_DSFP_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_DSFP_WR,
             "Global 3D Frictional Contact solver DSFP (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_TFP_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_TFP_WR,
             "Global 3D Frictional contact solver TFP (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSGS",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS,
             "Global 3D Frictional contact solver NSGS")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR,
             "Global 3D Non-smooth Newton Alart-Curnier solver (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSN_AC",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC,
             "Global 3D Non-smooth Newton Alart-Curnier solver")
      .value("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH,
             "Global 3D GAMS/Path (Ferris) solver")
      .value("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI,
             "Global 3D GAMS/Path VI solver")

      /** VI formulations */
      .value("SICONOS_GLOBAL_FRICTION_3D_VI_FPP",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_FPP,
             "Global 3D VI formulation, Fixed Point Projection")
      .value("SICONOS_GLOBAL_FRICTION_3D_VI_EG",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_EG,
             "Global 3D VI formulation, Extra-gradient")
      .value("SICONOS_GLOBAL_FRICTION_3D_ACLMFP",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ACLMFP,
             "Global 3D Alart-Curnier fixed point")
      .value("SICONOS_GLOBAL_FRICTION_3D_ADMM",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM, "Global 3D ADMM solver")
      .value("SICONOS_GLOBAL_FRICTION_3D_ADMM_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM_WR,
             "Global 3D ADMM solver (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_IPM", FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_IPM,
             "Global 3D Interior-point method solver")

      // Rolling Friction solvers (3D and 2D)
      .value("SICONOS_ROLLING_FRICTION_3D_NSGS",
             FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_NSGS,
             "3D Non-smooth Gauss Seidel, local formulation")
      .value("SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone",
             FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone,
             "3D Rolling friction one contact Projection on Cone")
      .value("SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration",
             FRICTION_SOLVER::
                 SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
             "3D Rolling friction one contact Projection on Cone with local iteration")
      .value("SICONOS_ROLLING_FRICTION_3D_ADMM",
             FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ADMM,
             "3D Rolling friction ADMM solver")

      // Rolling friction solvers for 2D problems
      .value("SICONOS_ROLLING_FRICTION_2D_NSGS",
             FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_NSGS,
             "2D Non-smooth Gauss Seidel, local formulation")
      .value("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone",
             FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone,
             "2D Rolling friction one contact Projection on Cone")
      .value("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration",
             FRICTION_SOLVER::
                 SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration,
             "2D Rolling friction one contact Projection on Cone with local iteration")

      // Global Rolling Friction solvers for 3D problems
      .value("SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR,
             "Global 3D Rolling friction solver NSGS (wrapped)")
      .value("SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM",
             FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM,
             "Global 3D Rolling friction solver IPM")
      .export_values();

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

  py::enum_<SICONOS_LCP_IPARAM>(params, "SICONOS_LCP_IPARAM", "LCP IPARAM enum")
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

  py::enum_<SICONOS_LCP_DPARAM>(params, "SICONOS_LCP_DPARAM", "LCP DPARAM enum")
      .value("SICONOS_LCP_DPARAM_RHO", SICONOS_LCP_DPARAM::SICONOS_LCP_DPARAM_RHO,
             "Relaxation or regularization parameter")
      .value("SICONOS_LCP_DPARAM_NSGS_LOCAL_ERROR_SUM",
             SICONOS_LCP_DPARAM::SICONOS_LCP_DPARAM_NSGS_LOCAL_ERROR_SUM,
             "Sum of local error values")
      .value("SICONOS_LCP_DPARAM_LATIN_PARAMETER",
             SICONOS_LCP_DPARAM::SICONOS_LCP_DPARAM_LATIN_PARAMETER, "Latin parameter")
      .export_values();

  py::enum_<SICONOS_LCP_SKIP_TRIVIAL>(params, "SICONOS_LCP_SKIP_TRIVIAL",
                                      "LCP skip trivial enum")
      .value("SICONOS_LCP_SKIP_TRIVIAL_NO",
             SICONOS_LCP_SKIP_TRIVIAL::SICONOS_LCP_SKIP_TRIVIAL_NO, "No trivial solution skip")
      .value("SICONOS_LCP_SKIP_TRIVIAL_YES",
             SICONOS_LCP_SKIP_TRIVIAL::SICONOS_LCP_SKIP_TRIVIAL_YES,
             "Yes trivial solution skip")
      .export_values();

  py::enum_<SICONOS_LCP_PIVOT_TYPE>(params, "SICONOS_LCP_PIVOT_TYPE", "LCP pivot type enum")
      .value("SICONOS_LCP_PIVOT_BARD", SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_BARD,
             "Bard pivoting")
      .value("SICONOS_LCP_PIVOT_LEAST_INDEX",
             SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_LEAST_INDEX, "Least index pivoting")
      .value("SICONOS_LCP_PIVOT_LEMKE", SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_LEMKE,
             "Lemke pivoting")
      .value("SICONOS_LCP_PIVOT_PATHSEARCH",
             SICONOS_LCP_PIVOT_TYPE::SICONOS_LCP_PIVOT_PATHSEARCH, "Pathsearch pivoting")
      .export_values();

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
