/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "BodyShapeRecord.hpp"
#include "Circle.hpp"
#include "Contact2d3DR.hpp"
#include "Contact2dR.hpp"
#include "Contact5DR.hpp"
#include "ContactR.hpp"
#include "Disk.hpp"
#include "InteractionManager.hpp"
#include "Lagrangian2d2DR.hpp"
#include "Lagrangian2d3DR.hpp"
#include "NewtonEuler3DR.hpp"
#include "NewtonEuler5DR.hpp"
#include "RigidBody2dDS.hpp"
#include "RigidBodyDS.hpp"
#include "SiconosCollisionManager.hpp"
#include "SiconosShape.hpp"
#include "SpaceFilter.hpp"
#include "StaticBody.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_collision, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");

  m.doc() = "Siconos mechanics.collision module";

  py::module_ module_simulation = py::module_::import("siconos.simulation");

  m.def("example_function", []() { return "This is an example function in mechanics"; });

  py::module_ bullet_module =
      m.def_submodule("bullet", "submodule bullet for mechanics.collision module");

  py::class_<siconos::collision::SiconosShape, py::smart_holder>(m, "SiconosShape")
      .def("setInsideMargin", &siconos::collision::SiconosShape::setInsideMargin)
      .def("setOutsideMargin", &siconos::collision::SiconosShape::setOutsideMargin);

  py::class_<siconos::collision::SiconosSphere, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosSphere")
      .def(py::init<float>(), py::arg("radius"))
      .def("radius", &siconos::collision::SiconosSphere::radius)
      .def("setRadius", &siconos::collision::SiconosSphere::setRadius);

  py::class_<siconos::collision::SiconosBox, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosBox")
      .def(py::init<double, double, double>(), py::arg("width"), py::arg("height"),
           py::arg("depth"))
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>>(), py::arg("dimensions"));

  py::class_<siconos::collision::SiconosCylinder, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosCylinder")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCylinder::radius)
      .def("setRadius", &siconos::collision::SiconosCylinder::setRadius)
      .def("length", &siconos::collision::SiconosCylinder::length)
      .def("setLength", &siconos::collision::SiconosCylinder::setLength);

  py::class_<siconos::collision::SiconosCone, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosCone")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCone::radius)
      .def("setRadius", &siconos::collision::SiconosCone::setRadius)
      .def("length", &siconos::collision::SiconosCone::length)
      .def("setLength", &siconos::collision::SiconosCone::setLength);

  py::class_<siconos::collision::SiconosCapsule, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosCapsule")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCapsule::radius)
      .def("setRadius", &siconos::collision::SiconosCapsule::setRadius)
      .def("length", &siconos::collision::SiconosCapsule::length)
      .def("setLength", &siconos::collision::SiconosCapsule::setLength);

  py::class_<siconos::collision::SiconosPlane, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosPlane")
      .def(py::init());

  py::class_<siconos::collision::SiconosConvexHull, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosConvexHull")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::arg("vertices"));  // TODO : WARNING with vertices argument, it is shard ptr

  py::class_<siconos::collision::SiconosDisk, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosDisk")
      .def(py::init<float>(), py::arg("radius"))
      .def("radius", &siconos::collision::SiconosDisk::radius)
      .def("setRadius", &siconos::collision::SiconosDisk::setRadius);

  py::class_<siconos::collision::SiconosBox2d, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosBox2d")
      .def(py::init<double, double>(), py::arg("width"), py::arg("height"));

  py::class_<siconos::collision::SiconosConvexHull2d, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosConvexHull2d")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::arg("vertices"));  // TODO : WARNING with vertices argument, it is shard ptr

  py::class_<siconos::collision::SiconosContactor, py::smart_holder>(m, "SiconosContactor")
      .def(py::init<std::shared_ptr<siconos::collision::SiconosShape>>(), py::arg("shape"))
      .def(py::init<std::shared_ptr<siconos::collision::SiconosShape>,
                    const siconos::algebra::SiconosVector &, int>(),
           py::arg("shape"), py::arg("offset"), py::arg("collision_group") = 0);

  py::class_<siconos::collision::SiconosContactorSet, py::smart_holder>(m,
                                                                        "SiconosContactorSet")
      .def(py::init<>())
      .def("append",
           py::overload_cast<std::shared_ptr<siconos::collision::SiconosContactor>>(
               &siconos::collision::SiconosContactorSet::append),
           py::arg("contactor"));

  py::class_<siconos::collision::SiconosMesh, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosMesh")
      .def(py::init<std::vector<unsigned int>, Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(), py::arg("indices"), py::arg("vertices"));

  py::class_<siconos::collision::SiconosHeightMap, siconos::collision::SiconosShape,
             py::smart_holder>(m, "SiconosHeightMap")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>, double, double>(),
           py::keep_alive<1, 1>(), py::arg("height_data"), py::arg("length_x"),
           py::arg("length_y"));

  py::class_<siconos::collision::native::SpaceFilter, siconos::simulation::InteractionManager,
             py::smart_holder>(m, "SpaceFilter")
      .def(py::init<>())
      .def("setBBoxfactor", &siconos::collision::native::SpaceFilter::setBBoxfactor)
      .def("setCellsize", &siconos::collision::native::SpaceFilter::setCellsize)
      .def("insertLine", &siconos::collision::native::SpaceFilter::insertLine);

  py::class_<siconos::collision::native::bodies::CircularDS, siconos::modeling::LagrangianDS,
             py::smart_holder>(m, "CircularDS")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 4>(), py::keep_alive<1, 5>(), py::arg("radius"), py::arg("mass"),
           py::arg("q0"), py::arg("v0"))
      .def("getMassValue", &siconos::collision::native::bodies::CircularDS::getMassValue);

  py::class_<siconos::collision::native::bodies::Circle,
             siconos::collision::native::bodies::CircularDS, py::smart_holder>(m, "Circle")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 4>(), py::keep_alive<1, 5>(), py::arg("radius"), py::arg("mass"),
           py::arg("q0"), py::arg("v0"));

  py::class_<siconos::collision::native::bodies::Disk,
             siconos::collision::native::bodies::CircularDS, py::smart_holder>(m, "Disk")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 4>(), py::keep_alive<1, 5>(), py::arg("radius"), py::arg("mass"),
           py::arg("q0"), py::arg("v0"));

  py::class_<siconos::collision::SiconosCollisionManager,
             siconos::simulation::InteractionManager, py::smart_holder>(
      m, "SiconosCollisionManager")
      .def("insertNonSmoothLaw",
           &siconos::collision::SiconosCollisionManager::insertNonSmoothLaw);

  py::class_<siconos::collision::StaticBody, py::smart_holder>(m, "StaticBody")
      .def(py::init<>())
      .def_readonly("number", &siconos::collision::StaticBody::number);

  py::class_<siconos::collision::BodyShapeRecord, py::smart_holder>(m, "BodyShapeRecord")
      .def("display", &siconos::collision::BodyShapeRecord::display)
      .def("__repr__",
           [](const siconos::collision::BodyShapeRecord &a) {
             a.display();
             return "\n";
           })
      .def_readonly("staticBody", &siconos::collision::BodyShapeRecord::staticBody)
      .def_readonly("ds", &siconos::collision::BodyShapeRecord::ds);

  py::class_<siconos::collision::RigidBodyDS, siconos::modeling::NewtonEulerDS,
             py::smart_holder>(m, "RigidBodyDS")
      //       .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector7>,
      //                     Eigen::Ref<siconos::algebra::SiconosVector6>, double,
      //                     Eigen::Ref<siconos::algebra::SiconosMatrix33>>(),
      //            py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
      //            alive
      //                                     // as long as object is referenced
      //            py::keep_alive<1, 3>(), py::keep_alive<1, 5>(), py::arg("position"),
      //            py::arg("velocity"), py::arg("mass"), py::arg("inertia"))
      .def(py::init<const siconos::algebra::SiconosVector7 &,
                    const siconos::algebra::SiconosVector6 &, double,
                    const siconos::algebra::SiconosMatrix33 &>(),
           py::arg("q0"), py::arg("twist0"), py::arg("mass"), py::arg("inertia"))
      //     py::arg("copy_t"))

      .def("setUseContactorInertia", &siconos::collision::RigidBodyDS::setUseContactorInertia,
           py::arg("useContactorInertia"))
      .def("setContactors", &siconos::collision::RigidBodyDS::setContactors,
           py::arg("contactors"))
      .def("contactors", &siconos::collision::RigidBodyDS::contactors);

  py::class_<siconos::collision::RigidBody2dDS, siconos::modeling::LagrangianLinearTIDS,
             py::smart_holder>(m, "RigidBody2dDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>, double, double>(),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"),
           py::arg("mass"), py::arg("inertia"))
      .def("setUseContactorInertia",
           &siconos::collision::RigidBody2dDS::setUseContactorInertia)
      .def("setContactors", &siconos::collision::RigidBody2dDS::setContactors)
      .def("setAllowSelfCollide", &siconos::collision::RigidBody2dDS::setAllowSelfCollide)
      .def_property_readonly("useContactorInertia",
                             &siconos::collision::RigidBody2dDS::useContactorInertia)
      .def_property_readonly("allowSelfCollide",
                             &siconos::collision::RigidBody2dDS::allowSelfCollide)
      .def_property_readonly("scalarMass", &siconos::collision::RigidBody2dDS::scalarMass);

  py::class_<siconos::collision::ContactR, siconos::modeling::NewtonEuler3DR,
             py::smart_holder>(m, "ContactR")
      .def_readonly("bodyShapeRecordA", &siconos::collision::ContactR::bodyShapeRecordA)
      .def_readonly("bodyShapeRecordB", &siconos::collision::ContactR::bodyShapeRecordB);

  py::class_<siconos::collision::Contact5DR, siconos::modeling::NewtonEuler5DR,
             py::smart_holder>(m, "Contact5DR")
      .def_readonly("bodyShapeRecordA", &siconos::collision::Contact5DR::bodyShapeRecordA)
      .def_readonly("bodyShapeRecordB", &siconos::collision::Contact5DR::bodyShapeRecordB);

  py::class_<siconos::collision::Contact2d3DR, siconos::modeling::Lagrangian2d3DR,
             py::smart_holder>(m, "Contact2d3DR")
      .def_readonly("bodyShapeRecordA", &siconos::collision::Contact2d3DR::bodyShapeRecordA)
      .def_readonly("bodyShapeRecordB", &siconos::collision::Contact2d3DR::bodyShapeRecordB);

  py::class_<siconos::collision::Contact2dR, siconos::modeling::Lagrangian2d2DR,
             py::smart_holder>(m, "Contact2dR")
      .def_readonly("bodyShapeRecordA", &siconos::collision::Contact2dR::bodyShapeRecordA)
      .def_readonly("bodyShapeRecordB", &siconos::collision::Contact2dR::bodyShapeRecordB);
}
