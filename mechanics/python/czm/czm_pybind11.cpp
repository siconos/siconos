/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "BinaryCohesiveNSL.hpp"
#include "CohesiveZoneModelNIFNSL.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_czm, m) {
  // Import the modeling module to get access to NewtonImpactFrictionNSL base class
  py::module_ modeling_module = py::module_::import("siconos.modeling");

  m.doc() = "Siconos mechanics.czm module - Cohesive Zone Models";

  // Expose the ShapeType enum
  py::enum_<siconos::mechanics::czm::BinaryCohesiveNSL::ShapeType>(m, "ShapeType")
      .value("DOOR_SHAPE", siconos::mechanics::czm::BinaryCohesiveNSL::ShapeType::DOOR_SHAPE,
             "Step function shape (constant traction until failure)")
      .value("TRIANGLE_SHAPE", siconos::mechanics::czm::BinaryCohesiveNSL::ShapeType::TRIANGLE_SHAPE,
             "Linear degradation shape (linear softening)")
      .export_values();

  // Expose the InternalVariables enum
  py::enum_<siconos::mechanics::czm::BinaryCohesiveNSL::InternalVariables>(m, "InternalVariables")
      .value("BETA_SURFACE", siconos::mechanics::czm::BinaryCohesiveNSL::BETA_SURFACE,
             "Damage parameter (0=broken, 1=intact) and surface area")
      .value("R_COHESION", siconos::mechanics::czm::BinaryCohesiveNSL::R_COHESION,
             "Cohesion force vector")
      .value("DISPLACEMENT_JUMP", siconos::mechanics::czm::BinaryCohesiveNSL::DISPLACEMENT_JUMP,
             "Current displacement jump")
      .value("COHESIVE_POINT_1", siconos::mechanics::czm::BinaryCohesiveNSL::COHESIVE_POINT_1,
             "Contact point on body 1")
      .value("COHESIVE_POINT_2", siconos::mechanics::czm::BinaryCohesiveNSL::COHESIVE_POINT_2,
             "Contact point on body 2")
      .value("NORMAL", siconos::mechanics::czm::BinaryCohesiveNSL::NORMAL,
             "Normal vector")
      .value("TANGENT_1", siconos::mechanics::czm::BinaryCohesiveNSL::TANGENT_1,
             "First tangent vector")
      .value("TANGENT_2", siconos::mechanics::czm::BinaryCohesiveNSL::TANGENT_2,
             "Second tangent vector")
      .value("INITIAL_DISPLACEMENT_JUMP", siconos::mechanics::czm::BinaryCohesiveNSL::INITIAL_DISPLACEMENT_JUMP,
             "Initial displacement jump")
      .value("INITIAL_RELATIVE_COHESIVE_POINT_1", siconos::mechanics::czm::BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_1,
             "Initial relative contact point 1")
      .value("INITIAL_RELATIVE_COHESIVE_POINT_2", siconos::mechanics::czm::BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_2,
             "Initial relative contact point 2")
      .value("INITIAL_RELATIVE_NORMAL", siconos::mechanics::czm::BinaryCohesiveNSL::INITIAL_RELATIVE_NORMAL,
             "Initial relative normal")
      .value("INITIAL_RELATIVE_TANGENT_1", siconos::mechanics::czm::BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1,
             "Initial relative tangent 1")
      .value("INITIAL_RELATIVE_TANGENT_2", siconos::mechanics::czm::BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2,
             "Initial relative tangent 2")
      .value("INTERNAL_VARIABLE_LENGTH", siconos::mechanics::czm::BinaryCohesiveNSL::INTERNAL_VARIABLE_LENGTH,
             "Total number of internal variables")
      .export_values();


  // Expose the BinaryCohesiveNSL class
  py::class_<siconos::mechanics::czm::BinaryCohesiveNSL,
             siconos::modeling::CohesiveZoneModelNIFNSL,
             py::smart_holder>(m, "BinaryCohesiveNSL")
      // Constructor with just size
      .def(py::init<siconos::algebra::Index>(),
           py::arg("size"),
           "Construct a BinaryCohesiveNSL with just size")
      // Constructor with en, et, mu, sigma_c, delta_c, size
      .def(py::init<double, double, double, double, double, siconos::algebra::Index>(),
           py::arg("en"), py::arg("et"), py::arg("mu"),
           py::arg("sigma_c"), py::arg("delta_c"), py::arg("size"),
           "Construct a BinaryCohesiveNSL with restitution coefficients and cohesive parameters")
      // Constructor with en, et, mu, sigma_c, delta_c, size, shape_type
      .def(py::init<double, double, double, double, double,
                    siconos::algebra::Index,
                    siconos::mechanics::czm::BinaryCohesiveNSL::ShapeType>(),
           py::arg("en"), py::arg("et"), py::arg("mu"),
           py::arg("sigma_c"), py::arg("delta_c"), py::arg("size"),
           py::arg("shape_type"),
           "Construct a BinaryCohesiveNSL with shape type")
      // Getters and setters
      .def_property("sigma_c",
                    &siconos::mechanics::czm::BinaryCohesiveNSL::sigma_c,
                    &siconos::mechanics::czm::BinaryCohesiveNSL::setSigma_c,
                    "Cohesive resistance to traction (critical stress)")
      .def_property("delta_c",
                    &siconos::mechanics::czm::BinaryCohesiveNSL::delta_c,
                    &siconos::mechanics::czm::BinaryCohesiveNSL::setDelta_c,
                    "Critical displacement for failure")
      .def_property_readonly("shape_type",
                    &siconos::mechanics::czm::BinaryCohesiveNSL::shape_type,
			     "Shape type of the cohesive law");
      // // Methods
      // .def("initializeInternalVariables",
      //      &siconos::mechanics::czm::BinaryCohesiveNSL::initializeInternalVariables,
      //      py::arg("inter"),
      //      "Initialize internal variables for this cohesive law")
      // .def("updateInternalVariables",
      //      &siconos::mechanics::czm::BinaryCohesiveNSL::updateInternalVariables,
      //      py::arg("inter"),
      //      "Update internal variables after each time step")
      // .def("beta",
      //      &siconos::mechanics::czm::BinaryCohesiveNSL::beta,
      //      py::arg("inter"),
      //      "Get the damage parameter beta (1.0 = intact, 0.0 = broken)")
      // .def("displayInternalVariables",
      //      [](siconos::mechanics::czm::BinaryCohesiveNSL& self,
      //         siconos::algebra::blocks::SharedVector& internalVariables) {
      //        self.displayInternalVariables(internalVariables);
      //      },
      //      py::arg("internalVariables"),
      //      "Display internal variables for debugging");
}
