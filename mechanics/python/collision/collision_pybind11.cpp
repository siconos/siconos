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

#include "Circle.hpp"
#include "Disk.hpp"
#include "InteractionManager.hpp"
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

  py::class_<siconos::collision::SiconosShape,
             std::shared_ptr<siconos::collision::SiconosShape>>(m, "SiconosShape")
      .def("setInsideMargin", &siconos::collision::SiconosShape::setInsideMargin)
      .def("setOutsideMargin", &siconos::collision::SiconosShape::setOutsideMargin);

  py::class_<siconos::collision::SiconosSphere,
             std::shared_ptr<siconos::collision::SiconosSphere>,
             siconos::collision::SiconosShape>(m, "SiconosSphere")
      .def(py::init<float>(), py::arg("radius"))
      .def("radius", &siconos::collision::SiconosSphere::radius)
      .def("setRadius", &siconos::collision::SiconosSphere::setRadius);

  py::class_<siconos::collision::SiconosBox, std::shared_ptr<siconos::collision::SiconosBox>,
             siconos::collision::SiconosShape>(m, "SiconosBox")
      .def(py::init<double, double, double>(), py::arg("width"), py::arg("height"),
           py::arg("depth"))
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>>(), py::arg("dimensions"));

  py::class_<siconos::collision::SiconosCylinder,
             std::shared_ptr<siconos::collision::SiconosCylinder>,
             siconos::collision::SiconosShape>(m, "SiconosCylinder")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCylinder::radius)
      .def("setRadius", &siconos::collision::SiconosCylinder::setRadius)
      .def("length", &siconos::collision::SiconosCylinder::length)
      .def("setLength", &siconos::collision::SiconosCylinder::setLength);

  py::class_<siconos::collision::SiconosCone, std::shared_ptr<siconos::collision::SiconosCone>,
             siconos::collision::SiconosShape>(m, "SiconosCone")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCone::radius)
      .def("setRadius", &siconos::collision::SiconosCone::setRadius)
      .def("length", &siconos::collision::SiconosCone::length)
      .def("setLength", &siconos::collision::SiconosCone::setLength);

  py::class_<siconos::collision::SiconosCapsule,
             std::shared_ptr<siconos::collision::SiconosCapsule>,
             siconos::collision::SiconosShape>(m, "SiconosCapsule")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCapsule::radius)
      .def("setRadius", &siconos::collision::SiconosCapsule::setRadius)
      .def("length", &siconos::collision::SiconosCapsule::length)
      .def("setLength", &siconos::collision::SiconosCapsule::setLength);

  py::class_<siconos::collision::SiconosPlane,
             std::shared_ptr<siconos::collision::SiconosPlane>,
             siconos::collision::SiconosShape>(m, "SiconosPlane")
      .def(py::init());

  py::class_<siconos::collision::SiconosConvexHull,
             std::shared_ptr<siconos::collision::SiconosConvexHull>,
             siconos::collision::SiconosShape>(m, "SiconosConvexHull")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::arg("vertices"));  // TODO : WARNING with vertices argument, it is shard ptr

  py::class_<siconos::collision::SiconosDisk, std::shared_ptr<siconos::collision::SiconosDisk>,
             siconos::collision::SiconosShape>(m, "SiconosDisk")
      .def(py::init<float>(), py::arg("radius"))
      .def("radius", &siconos::collision::SiconosDisk::radius)
      .def("setRadius", &siconos::collision::SiconosDisk::setRadius);

  py::class_<siconos::collision::SiconosBox2d,
             std::shared_ptr<siconos::collision::SiconosBox2d>,
             siconos::collision::SiconosShape>(m, "SiconosBox2d")
      .def(py::init<double, double>(), py::arg("width"), py::arg("height"));

  py::class_<siconos::collision::SiconosConvexHull2d,
             std::shared_ptr<siconos::collision::SiconosConvexHull2d>,
             siconos::collision::SiconosShape>(m, "SiconosConvexHull2d")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::arg("vertices"));  // TODO : WARNING with vertices argument, it is shard ptr

  py::class_<siconos::collision::SiconosContactor,
             std::shared_ptr<siconos::collision::SiconosContactor>>(m, "SiconosContactor")
      .def(py::init<std::shared_ptr<siconos::collision::SiconosShape>,
                    const siconos::algebra::SiconosVector &, int>(),
           py::arg("shape"), py::arg("offset"), py::arg("collision_group") = 0);

  py::class_<siconos::collision::SiconosContactorSet,
             std::shared_ptr<siconos::collision::SiconosContactorSet>>(m,
                                                                       "SiconosContactorSet")
      .def(py::init<>())
      .def("append",
           py::overload_cast<std::shared_ptr<siconos::collision::SiconosContactor>>(
               &siconos::collision::SiconosContactorSet::append),
           py::arg("contactor"));

  py::class_<siconos::collision::SiconosMesh, std::shared_ptr<siconos::collision::SiconosMesh>,
             siconos::collision::SiconosShape>(m, "SiconosMesh");

  py::class_<siconos::collision::SiconosHeightMap,
             std::shared_ptr<siconos::collision::SiconosHeightMap>,
             siconos::collision::SiconosShape>(m, "SiconosHeightMap");

  py::class_<siconos::collision::native::SpaceFilter,
             std::shared_ptr<siconos::collision::native::SpaceFilter>>(m, "SpaceFilter");

  py::class_<siconos::collision::native::bodies::Circle,
             std::shared_ptr<siconos::collision::native::bodies::Circle>>(m, "Circle")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::arg("radius"), py::arg("mass"), py::arg("position"), py::arg("velocity"));

  py::class_<siconos::collision::native::bodies::Disk,
             std::shared_ptr<siconos::collision::native::bodies::Disk>>(m, "Disk")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::arg("radius"), py::arg("mass"), py::arg("position"), py::arg("velocity"));

  py::class_<siconos::collision::SiconosCollisionManager,
             std::shared_ptr<siconos::collision::SiconosCollisionManager>,
             siconos::simulation::InteractionManager>(m, "SiconosCollisionManager")
      .def("insertNonSmoothLaw",
           &siconos::collision::SiconosCollisionManager::insertNonSmoothLaw);

  py::class_<siconos::collision::StaticBody, std::shared_ptr<siconos::collision::StaticBody>>(
      m, "StaticBody")
      .def(py::init<>());

  py::class_<siconos::collision::RigidBodyDS, std::shared_ptr<siconos::collision::RigidBodyDS>,
             siconos::modeling::NewtonEulerDS>(m, "RigidBodyDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>, double,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 5>(), py::arg("position"),
           py::arg("velocity"), py::arg("mass"), py::arg("inertia"))
      .def("setUseContactorInertia", &siconos::collision::RigidBodyDS::setUseContactorInertia,
           py::arg("useContactorInertia"))
      .def("setContactors", &siconos::collision::RigidBodyDS::setContactors,
           py::arg("contactors"))
      .def("contactors", &siconos::collision::RigidBodyDS::contactors);

  py::class_<siconos::collision::RigidBody2dDS,
             std::shared_ptr<siconos::collision::RigidBody2dDS>,
             siconos::modeling::LagrangianLinearTIDS>(m, "RigidBody2dDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("position"),
           py::arg("velocity"), py::arg("mass"));
}