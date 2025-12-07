// # Siconos is a program dedicated to modeling, simulation and control
// # of non smooth dynamical systems.
// #
// # Copyright 2025 INRIA.
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

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "SiconosBulletCollisionManager.hpp"
#include "SiconosBulletOptions.hpp"
#include "SiconosContactor.hpp"
#include "SecondOrderDS.hpp"
namespace py = pybind11;

PYBIND11_MODULE(_bullet, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");
  py::module_ collision_module = py::module_::import("siconos.mechanics.collision");

  m.doc() = "Siconos mechanics.collision.bullet module";

  py::enum_<siconos::collision::bullet::SiconosBulletDimension>(m, "SiconosBulletDimension",
                                                                "Bullet dim")
      .value("ThreeD", siconos::collision::bullet::SiconosBulletDimension::ThreeD, "3D case")
      .value("TwoD", siconos::collision::bullet::SiconosBulletDimension::TwoD, "2D case")
      .export_values()
      .def("__str__",
           [](siconos::collision::bullet::SiconosBulletDimension d) {
             return d == siconos::collision::bullet::SiconosBulletDimension::ThreeD ? "ThreeD"
                                                                                    : "TwoD";
           })
      .def("__repr__", [](siconos::collision::bullet::SiconosBulletDimension d) {
        return d == siconos::collision::bullet::SiconosBulletDimension::ThreeD ? "\"ThreeD\""
                                                                               : "\"TwoD\"";
      });
  ;

  py::class_<siconos::collision::bullet::SiconosBulletOptions,
             std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions>>(
      m, "SiconosBulletOptions")
      .def(py::init<>())
      .def_readwrite("dimension", &siconos::collision::bullet::SiconosBulletOptions::dimension,
                     "Dimension setting")
      .def_readwrite(
          "contactBreakingThreshold",
          &siconos::collision::bullet::SiconosBulletOptions::contactBreakingThreshold,
          "Contact breaking threshold")
      .def_readwrite(
          "contactProcessingThreshold",
          &siconos::collision::bullet::SiconosBulletOptions::contactProcessingThreshold,
          "Contact processing threshold")
      .def_readwrite("worldScale",
                     &siconos::collision::bullet::SiconosBulletOptions::worldScale,
                     "World scale factor")
      .def_readwrite("useAxisSweep3",
                     &siconos::collision::bullet::SiconosBulletOptions::useAxisSweep3,
                     "Use Axis Sweep 3")
      .def_readwrite(
          "clearOverlappingPairCache",
          &siconos::collision::bullet::SiconosBulletOptions::clearOverlappingPairCache,
          "Clear overlapping pair cache")
      .def_readwrite("perturbationIterations",
                     &siconos::collision::bullet::SiconosBulletOptions::perturbationIterations,
                     "Number of perturbation iterations")
      .def_readwrite("minimumPointsPerturbationThreshold",
                     &siconos::collision::bullet::SiconosBulletOptions::
                         minimumPointsPerturbationThreshold,
                     "Minimum points perturbation threshold")
      .def_readwrite("enableSatConvex",
                     &siconos::collision::bullet::SiconosBulletOptions::enableSatConvex,
                     "Enable SAT convex")
      .def_readwrite(
          "enablePolyhedralContactClipping",
          &siconos::collision::bullet::SiconosBulletOptions::enablePolyhedralContactClipping,
          "Enable polyhedral contact clipping")
      .def_readwrite("Depth2D", &siconos::collision::bullet::SiconosBulletOptions::Depth2D,
                     "2D depth value")
      .def_readwrite(
          "extrapolationCoefficient",
          &siconos::collision::bullet::SiconosBulletOptions::extrapolationCoefficient,
          "Extrapolation coefficient");

  py::class_<siconos::collision::bullet::SiconosBulletStatistics>(m, "SiconosBulletStatistics")
      .def(py::init<>())  // Constructeur par défaut
      .def_readwrite(
          "new_interactions_created",
          &siconos::collision::bullet::SiconosBulletStatistics::new_interactions_created,
          "New interactions created")
      .def_readwrite("existing_interactions_processed",
                     &siconos::collision::bullet::SiconosBulletStatistics::
                         existing_interactions_processed,
                     "Existing interactions processed")
      .def_readwrite(
          "interaction_warnings",
          &siconos::collision::bullet::SiconosBulletStatistics::interaction_warnings,
          "Interaction warnings")
      .def_readwrite(
          "interaction_destroyed",
          &siconos::collision::bullet::SiconosBulletStatistics::interaction_destroyed,
          "Interaction destroyed");

  py::class_<siconos::collision::bullet::SiconosBulletCollisionManager,
             std::shared_ptr<siconos::collision::bullet::SiconosBulletCollisionManager>,
             siconos::collision::SiconosCollisionManager>(m, "SiconosBulletCollisionManager")
      .def(py::init<std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions>>(),
           py::arg("options"))
      .def(py::init<>())
      .def("addStaticBody",
           &siconos::collision::bullet::SiconosBulletCollisionManager::addStaticBody,
           py::arg("cs"), py::arg("position"), py::arg("number") = 0)
     .def("removeStaticBody",
           &siconos::collision::bullet::SiconosBulletCollisionManager::removeStaticBody,
	   py::arg("body") = 0)
      .def("removeBody",
           &siconos::collision::bullet::SiconosBulletCollisionManager::removeBody,
           py::arg("ds") = 0)
      .def("useEqualityConstraints",
           &siconos::collision::bullet::SiconosBulletCollisionManager::useEqualityConstraints)
      .def("statistics",
           &siconos::collision::bullet::SiconosBulletCollisionManager::statistics);
}
