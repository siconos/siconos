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

// #include <cassert>

#include "FrictionContactProblem.h"
#include "Friction_cst.h"
#include "GlobalFrictionContactProblem.h"
#include "GlobalRollingFrictionContactProblem.h"
#include "NumericsMatrix_pybind11.h"
#include "RollingFrictionContactProblem.h"
namespace py = pybind11;

void wrap_friction_contact(py::module_ &m, py::module_ &params, py::module_ &solver_ids) {
  py::class_<FrictionContactProblem>(m, "FrictionContactProblem")
      .def(py::init([](int dimension, int numberOfContacts) {
             return frictionContactProblem_new_with_data(dimension, numberOfContacts, nullptr,
                                                         nullptr, nullptr);
           }),
           py::arg("dimension"), py::arg("numberOfContacts"))
      .def("__del__",
           [](FrictionContactProblem *self) {
             if (self) {
               frictionContactProblem_free(self);
             }
           })
      .def_readwrite("dimension", &FrictionContactProblem::dimension)
      .def_readwrite("numberOfContacts", &FrictionContactProblem::numberOfContacts)
      .def_property(
          "M",
          [](const FrictionContactProblem &self) {
            return get_matrix(self, &FrictionContactProblem::M);
          },
          [](FrictionContactProblem &self, py::array_t<double> array) {
            set_matrix(self, &FrictionContactProblem::M, array);
            assert(self.M->size0 == self.dimension * self.numberOfContacts);
            assert(self.M->size1 == self.dimension * self.numberOfContacts);
          })
      .def_property(
          "q",
          [](const FrictionContactProblem &self) {
            return get_array(self, self.q, self.numberOfContacts * self.dimension);
          },
          [](FrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.q, arr, self.numberOfContacts * self.dimension);
          },
          "Vector q (shared memory)")
      .def_property(
          "mu",
          [](const FrictionContactProblem &self) {
            return get_array(self, self.mu, self.numberOfContacts);
          },
          [](FrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.mu, arr, self.numberOfContacts);
          },
          "Vector of friction coeffs");

  py::class_<RollingFrictionContactProblem>(m, "RollingFrictionContactProblem")
      .def(py::init([](int dimension, int numberOfContacts) {
             return rollingFrictionContactProblem_new_with_data(
                 dimension, numberOfContacts, nullptr, nullptr, nullptr, nullptr);
           }),
           py::arg("dimension"), py::arg("numberOfContacts"))
      .def("__del__",
           [](RollingFrictionContactProblem *self) {
             if (self) {
               rollingFrictionContactProblem_free(self);
             }
           })
      .def_readwrite("dimension", &RollingFrictionContactProblem::dimension)
      .def_readwrite("numberOfContacts", &RollingFrictionContactProblem::numberOfContacts)
      .def_property(
          "M",
          [](const RollingFrictionContactProblem &self) {
            return get_matrix(self, &RollingFrictionContactProblem::M);
          },
          [](RollingFrictionContactProblem &self, py::array_t<double> array) {
            set_matrix(self, &RollingFrictionContactProblem::M, array);
            assert(self.M->size0 == self.dimension * self.numberOfContacts);
            assert(self.M->size1 == self.dimension * self.numberOfContacts);
          })
      .def_property(
          "q",
          [](const RollingFrictionContactProblem &self) {
            return get_array(self, self.q, self.numberOfContacts * self.dimension);
          },
          [](RollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.q, arr, self.numberOfContacts * self.dimension);
          },
          "Vector q (shared memory)")
      .def_property(
          "mu",
          [](const RollingFrictionContactProblem &self) {
            return get_array(self, self.mu, self.numberOfContacts);
          },
          [](RollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.mu, arr, self.numberOfContacts);
          },
          "Vector of friction coeffs")
      .def_property(
          "mu_r",
          [](const RollingFrictionContactProblem &self) {
            return get_array(self, self.mu_r, self.numberOfContacts);
          },
          [](RollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.mu_r, arr, self.numberOfContacts);
          },
          "Vector of friction coeffs");
  ;

  py::class_<GlobalFrictionContactProblem, std::shared_ptr<GlobalFrictionContactProblem>>(
      m, "GlobalFrictionContactProblem")
      .def(py::init([](int dimension, int numberOfContacts) {
        auto pb = globalFrictionContactProblem_new();
        pb->dimension = dimension;
        pb->numberOfContacts = numberOfContacts;
        return pb;
      }))
      .def("__del__",
           [](GlobalFrictionContactProblem *self) {
             if (self) {
               globalFrictionContact_free(self);
             }
           })

      .def_readwrite("dimension", &GlobalFrictionContactProblem::dimension)
      .def_readwrite("numberOfContacts", &GlobalFrictionContactProblem::numberOfContacts)
      .def_property(
          "M",
          [](const GlobalFrictionContactProblem &self) {
            return get_matrix(self, &GlobalFrictionContactProblem::M);
          },
          [](GlobalFrictionContactProblem &self, py::array_t<double> array) {
            set_matrix(self, &GlobalFrictionContactProblem::M, array);
            assert(self.M->size0 == self.M->size1);
          })
      .def_property(
          "H",
          [](const GlobalFrictionContactProblem &self) {
            return get_matrix(self, &GlobalFrictionContactProblem::H);
          },
          [](GlobalFrictionContactProblem &self, py::array_t<double> array) {
            set_matrix(self, &GlobalFrictionContactProblem::H, array);
            assert(self.H->size1 == self.dimension * self.numberOfContacts);
          })
      .def_property(
          "q",
          [](const GlobalFrictionContactProblem &self) {
            assert(self.M &&
                   "M must be properly set before any access to q");  // the only way to get q
                                                                      // size ...
            return get_array(self, self.q, self.M->size0);
          },
          [](GlobalFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.q, arr);  //, self.M.size0, might be unknown ...
          },
          "Vector q (shared memory)")
      //       .def_property_readonly(
      //           "norm_q", [](const GlobalFrictionContactProblem &self) { return self.norm_q;
      //           })
      .def_property(
          "b",
          [](const GlobalFrictionContactProblem &self) {
            return get_array(self, self.b, self.dimension * self.numberOfContacts);
          },
          [](GlobalFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.b, arr, self.dimension * self.numberOfContacts);
          })
      //       .def_property_readonly(
      //           "norm_b", [](const GlobalFrictionContactProblem &self) { return self.norm_b;
      //           })
      .def_property(
          "mu",
          [](const GlobalFrictionContactProblem &self) {
            return get_array(self, self.mu, self.numberOfContacts);
          },
          [](GlobalFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.mu, arr, self.numberOfContacts);
          },
          "Vector of friction coeffs");

  py::class_<GlobalRollingFrictionContactProblem>(m, "GlobalRollingFrictionContactProblem")
      .def(py::init([](int dimension, int numberOfContacts) {
        return globalRollingFrictionContactProblem_new_with_data(
            dimension, numberOfContacts, nullptr, nullptr, nullptr, nullptr);
      }))
      .def("__del__",
           [](GlobalRollingFrictionContactProblem *self) {
             if (self) {
               globalRollingFrictionContactProblem_free(self);
             }
           })
      .def_readwrite("dimension", &GlobalRollingFrictionContactProblem::dimension)
      .def_readwrite("numberOfContacts",
                     &GlobalRollingFrictionContactProblem::numberOfContacts)
      .def_property(
          "M",
          [](const GlobalRollingFrictionContactProblem &self) {
            return get_matrix(self, &GlobalRollingFrictionContactProblem::M);
          },
          [](GlobalRollingFrictionContactProblem &self, py::array_t<double> array) {
            set_matrix(self, &GlobalRollingFrictionContactProblem::M, array);
            assert(self.M->size0 == self.M->size1);
          })
      .def_property(
          "H",
          [](const GlobalRollingFrictionContactProblem &self) {
            return get_matrix(self, &GlobalRollingFrictionContactProblem::H);
          },
          [](GlobalRollingFrictionContactProblem &self, py::array_t<double> array) {
            set_matrix(self, &GlobalRollingFrictionContactProblem::H, array);
            assert(self.H->size1 == self.dimension * self.numberOfContacts);
          })
      .def_property(
          "q",
          [](const GlobalRollingFrictionContactProblem &self) {
            assert(self.M &&
                   "M must be properly set before any access to q");  // the only way to get q
                                                                      // size ...

            return get_array(self, self.q, self.M->size0);
          },
          [](GlobalRollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.q, arr);  //, self.M.size0, might be unknown ...
          },
          "Vector q (shared memory)")
      .def_property(
          "b",
          [](const GlobalRollingFrictionContactProblem &self) {
            return get_array(self, self.b, self.dimension * self.numberOfContacts);
          },
          [](GlobalRollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.b, arr, self.dimension * self.numberOfContacts);
          })
      .def_property(
          "mu",
          [](const GlobalRollingFrictionContactProblem &self) {
            return get_array(self, self.mu, self.numberOfContacts);
          },
          [](GlobalRollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.mu, arr, self.numberOfContacts);
          },
          "Vector of friction coeffs")
      .def_property(
          "mu_r",
          [](const GlobalRollingFrictionContactProblem &self) {
            return get_array(self, self.mu_r, self.numberOfContacts);
          },
          [](GlobalRollingFrictionContactProblem &self, py::array_t<double> arr) {
            set_array(self, self.mu_r, arr, self.numberOfContacts);
          },
          "Vector of friction coeffs");
  ;

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

  py::enum_<SICONOS_FRICTION_3D_IPARAM>(params, "SICONOS_FRICTION_3D_IPARAM")
      .value("SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY",
             SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY)
      .value("SICONOS_FRICTION_3D_IPARAM_RESCALING", SICONOS_FRICTION_3D_IPARAM_RESCALING)
      .value("SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE",
             SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE)
      .value("SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER",
             SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER)
      .value("SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION",
             SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION)
      .value("SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY",
             SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY)
      .value("SICONOS_FRICTION_3D_NUMBER_OF_CONTACTS", SICONOS_FRICTION_3D_NUMBER_OF_CONTACTS)
      .export_values();

  py::enum_<SICONOS_FRICTION_INTERNAL_ERROR_STRATEGY>(
      params, "SICONOS_FRICTION_INTERNAL_ERROR_STRATEGY")
      .value("SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE",
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE)
      .value("SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE",
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE)
      .value("SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT",
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_RESCALING_ENUM>(params, "SICONOS_FRICTION_3D_RESCALING_ENUM")
      .value("SICONOS_FRICTION_3D_RESCALING_NO", SICONOS_FRICTION_3D_RESCALING_NO)
      .value("SICONOS_FRICTION_3D_RESCALING_SCALAR", SICONOS_FRICTION_3D_RESCALING_SCALAR)
      .value("SICONOS_FRICTION_3D_RESCALING_BALANCING_M",
             SICONOS_FRICTION_3D_RESCALING_BALANCING_M)
      .value("SICONOS_FRICTION_3D_RESCALING_BALANCING_MH",
             SICONOS_FRICTION_3D_RESCALING_BALANCING_MH)
      .value("SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT",
             SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_RESCALING_CONE_ENUM>(params,
                                                     "SICONOS_FRICTION_3D_RESCALING_CONE_ENUM")
      .value("SICONOS_FRICTION_3D_RESCALING_CONE_NO", SICONOS_FRICTION_3D_RESCALING_CONE_NO)
      .value("SICONOS_FRICTION_3D_RESCALING_CONE_YES", SICONOS_FRICTION_3D_RESCALING_CONE_YES)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_DPARAM>(params, "SICONOS_FRICTION_3D_DPARAM")
      .value("SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO",
             SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_IPARAM>(params, "SICONOS_FRICTION_3D_NSGS_IPARAM")
      .value("SICONOS_FRICTION_3D_NSGS_RELAXATION", SICONOS_FRICTION_3D_NSGS_RELAXATION)
      .value("SICONOS_FRICTION_3D_NSGS_SHUFFLE", SICONOS_FRICTION_3D_NSGS_SHUFFLE)
      .value("SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED", SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED)
      .value("SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT",
             SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT)
      .value("SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION",
             SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_DPARAM>(params, "SICONOS_FRICTION_3D_NSGS_DPARAM")
      .value("SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE",
             SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM>(
      params, "SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM")
      .value("SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION",
             SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION>(
      params, "SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION")
      .value("SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE",
             SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE)
      .value("SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE",
             SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ENUM>(
      params, "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ENUM")
      .value("SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL",
             SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL)
      .value("SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT",
             SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
      .value("SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL",
             SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
      .value("SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE",
             SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_SHUFFLE_ENUM>(params,
                                                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_ENUM")
      .value("SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE", SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE)
      .value("SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE", SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE)
      .value("SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP",
             SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_RELAXATION_ENUM>(
      params, "SICONOS_FRICTION_3D_NSGS_RELAXATION_ENUM")
      .value("SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE",
             SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE)
      .value("SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE",
             SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_ENUM>(
      params, "SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_ENUM")
      .value("SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE",
             SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE)
      .value("SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE",
             SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_IPARAM>(params, "SICONOS_FRICTION_3D_NSN_IPARAM")
      .value("SICONOS_FRICTION_3D_NSN_RHO_STRATEGY", SICONOS_FRICTION_3D_NSN_RHO_STRATEGY)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION", SICONOS_FRICTION_3D_NSN_FORMULATION)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH", SICONOS_FRICTION_3D_NSN_LINESEARCH)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER",
             SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER)
      .value("SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER", SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP",
             SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER",
             SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER)
      .value("SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED",
             SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED)
      .value("SICONOS_FRICTION_3D_NSN_MPI_COM", SICONOS_FRICTION_3D_NSN_MPI_COM)
      .export_values();

  py::enum_<SICONOS_FC3D_NSN_LINEAR_SOLVER>(params, "SICONOS_FC3D_NSN_LINEAR_SOLVER")
      .value("SICONOS_FRICTION_3D_NSN_USE_CSLUSOL", SICONOS_FRICTION_3D_NSN_USE_CSLUSOL)
      .value("SICONOS_FRICTION_3D_NSN_USE_MUMPS", SICONOS_FRICTION_3D_NSN_USE_MUMPS)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_DPARAM>(params, "SICONOS_FRICTION_3D_NSN_DPARAM")
      .value("SICONOS_FRICTION_3D_NSN_RHO", SICONOS_FRICTION_3D_NSN_RHO)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_ENUM>(
      params, "SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_ENUM")
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_FORMULATION_ENUM>(
      params, "SICONOS_FRICTION_3D_NSN_FORMULATION_ENUM")
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD",
             SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD",
             SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED",
             SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED",
             SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_NULL",
             SICONOS_FRICTION_3D_NSN_FORMULATION_NULL)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_LINESEARCH_ENUM>(params,
                                                     "SICONOS_FRICTION_3D_NSN_LINESEARCH_ENUM")
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE",
             SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO",
             SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_NO", SICONOS_FRICTION_3D_NSN_LINESEARCH_NO)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_HYBRID_ENUM>(params, "SICONOS_FRICTION_3D_NSN_HYBRID_ENUM")
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_PROXIMAL_IPARAM>(params, "SICONOS_FRICTION_3D_PROXIMAL_IPARAM")
      .value("SICONOS_FRICTION_3D_FP_ERROR_STRATEGY", SICONOS_FRICTION_3D_FP_ERROR_STRATEGY)
      .value("SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE",
             SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE)
      .value("SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION",
             SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION)
      .value("SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY",
             SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_PROXIMAL_DPARAM>(params, "SICONOS_FRICTION_3D_PROXIMAL_DPARAM")
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA",
             SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA)
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA",
             SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA)
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU", SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU)
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION",
             SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_PROXIMAL>(params, "SICONOS_FRICTION_3D_PROXIMAL")
      .value("SICONOS_FRICTION_3D_PROXIMAL_PROX", SICONOS_FRICTION_3D_PROXIMAL_PROX)
      .value("SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION",
             SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_IPARAM_ENUM>(params,
                                                  "SICONOS_FRICTION_3D_ADMM_IPARAM_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY",
             SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO",
             SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION",
             SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO",
             SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S",
             SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H", SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_DPARAM_ENUM>(params,
                                                  "SICONOS_FRICTION_3D_ADMM_DPARAM_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_RHO", SICONOS_FRICTION_3D_ADMM_RHO)
      .value("SICONOS_FRICTION_3D_ADMM_RESTART_ETA", SICONOS_FRICTION_3D_ADMM_RESTART_ETA)
      .value("SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU",
             SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU)
      .value("SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI",
             SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_ACCELERATION_ENUM>(
      params, "SICONOS_FRICTION_3D_ADMM_ACCELERATION_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION",
             SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION)
      .value("SICONOS_FRICTION_3D_ADMM_ACCELERATION", SICONOS_FRICTION_3D_ADMM_ACCELERATION)
      .value("SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART",
             SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_SYMMETRY_ENUM>(params,
                                                    "SICONOS_FRICTION_3D_ADMM_SYMMETRY_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY",
             SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_SYMMETRIZE", SICONOS_FRICTION_3D_ADMM_SYMMETRIZE)
      .value("SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_STORAGE_ENUM>(params,
                                                   "SICONOS_FRICTION_3D_ADMM_STORAGE_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE", SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE)
      .value("SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_ENUM>(
      params, "SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO",
             SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO)
      .value("SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES",
             SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_UPDATE_S_ENUM>(params,
                                                    "SICONOS_FRICTION_3D_ADMM_UPDATE_S_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES", SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
      .value("SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO", SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_FULL_H_ENUM>(params,
                                                  "SICONOS_FRICTION_3D_ADMM_FULL_H_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_FULL_H_NO", SICONOS_FRICTION_3D_ADMM_FULL_H_NO)
      .value("SICONOS_FRICTION_3D_ADMM_FULL_H_YES", SICONOS_FRICTION_3D_ADMM_FULL_H_YES)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_ENUM>(
      params, "SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT",
             SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT)
      .value("SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING",
             SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING)
      .value("SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING",
             SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_ENUM>(
      params, "SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_ENUM")
      .value("SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN",
             SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN)
      .value("SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF",
             SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF)
      .value("SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES",
             SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM_ENUM>(params, "SICONOS_FRICTION_3D_IPM_IPARAM_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING",
             SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO",
             SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE",
             SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING",
             SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S",
             SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD",
             SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM", SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY",
             SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE",
             SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_DPARAM_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_DPARAM_ENUM>(params, "SICONOS_FRICTION_3D_IPM_DPARAM_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1",
             SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1)
      .value("SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2",
             SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2)
      .value("SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3",
             SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3)
      .value("SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1",
             SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1)
      .value("SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2",
             SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_STORAGE_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_STORAGE_ENUM>(params,
                                                  "SICONOS_FRICTION_3D_IPM_STORAGE_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_KEEP_STORAGE", SICONOS_FRICTION_3D_IPM_KEEP_STORAGE)
      .value("SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO",
             SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO)
      .value("SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES",
             SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP",
             SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      .value("SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F",
             SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QP2",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QP2)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QPH)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_NO",
             SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_NO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_YES",
             SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_YES)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_ENUM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_ENUM>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_ENUM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO",
             SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES",
             SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES)
      .export_values();
}
