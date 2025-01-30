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

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "Circle.hpp"
#include "Disk.hpp"
#include "InteractionManager.hpp"
#include "SiconosBulletCollisionManager.hpp"
#include "SiconosBulletOptions.hpp"
#include "SiconosCollisionManager.hpp"
#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"
#include "SpaceFilter.hpp"
#include "StaticBody.hpp"
#include "RigidBody2dDS.hpp"
#include "RigidBodyDS.hpp"


namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<siconos::collision::SiconosContactor>>);

PYBIND11_MODULE(pymechanics, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");

  m.doc() = "Siconos mechanics module";

  py::module_ module_simulation = py::module_::import("siconos.simulation");

  m.def("example_function", []() { return "This is an example function in pymechanics"; });

  py::module_ joints_module =
      m.def_submodule("joints", "submodule joints for pymechanics module");  // TODO : fill in

  joints_module.def("example_function",
                    []() { return "This is an example function in joints"; });

  py::module_ collision_module =
      m.def_submodule("collision", "submodule collision for pymechanics module");

  py::module_ bullet_module = collision_module.def_submodule(
      "bullet", "submodule bullet for pymechanics.collision module");

  py::class_<siconos::collision::SiconosShape,
             std::shared_ptr<siconos::collision::SiconosShape>>(collision_module,
                                                                "SiconosShape")
      .def("setInsideMargin", &siconos::collision::SiconosShape::setInsideMargin)
      .def("setOutsideMargin", &siconos::collision::SiconosShape::setOutsideMargin);

  py::class_<siconos::collision::SiconosSphere,
             std::shared_ptr<siconos::collision::SiconosSphere>,
             siconos::collision::SiconosShape>(collision_module, "SiconosSphere")
      .def(py::init<float>(), py::arg("radius"))
      .def("radius", &siconos::collision::SiconosSphere::radius)
      .def("setRadius", &siconos::collision::SiconosSphere::setRadius);

  py::class_<siconos::collision::SiconosBox, std::shared_ptr<siconos::collision::SiconosBox>,
             siconos::collision::SiconosShape>(collision_module, "SiconosBox")
      .def(py::init<double, double, double>(), py::arg("width"), py::arg("height"),
           py::arg("depth"))
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>>(), py::arg("dimensions"));

  py::class_<siconos::collision::SiconosCylinder,
             std::shared_ptr<siconos::collision::SiconosCylinder>,
             siconos::collision::SiconosShape>(collision_module, "SiconosCylinder")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCylinder::radius)
      .def("setRadius", &siconos::collision::SiconosCylinder::setRadius)
      .def("length", &siconos::collision::SiconosCylinder::length)
      .def("setLength", &siconos::collision::SiconosCylinder::setLength);

  py::class_<siconos::collision::SiconosCone, std::shared_ptr<siconos::collision::SiconosCone>,
             siconos::collision::SiconosShape>(collision_module, "SiconosCone")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCone::radius)
      .def("setRadius", &siconos::collision::SiconosCone::setRadius)
      .def("length", &siconos::collision::SiconosCone::length)
      .def("setLength", &siconos::collision::SiconosCone::setLength);

  py::class_<siconos::collision::SiconosCapsule,
             std::shared_ptr<siconos::collision::SiconosCapsule>,
             siconos::collision::SiconosShape>(collision_module, "SiconosCapsule")
      .def(py::init<float, float>(), py::arg("radius"), py::arg("length"))
      .def("radius", &siconos::collision::SiconosCapsule::radius)
      .def("setRadius", &siconos::collision::SiconosCapsule::setRadius)
      .def("length", &siconos::collision::SiconosCapsule::length)
      .def("setLength", &siconos::collision::SiconosCapsule::setLength);

  py::class_<siconos::collision::SiconosPlane,
             std::shared_ptr<siconos::collision::SiconosPlane>,
             siconos::collision::SiconosShape>(collision_module, "SiconosPlane");

  py::class_<siconos::collision::SiconosConvexHull,
             std::shared_ptr<siconos::collision::SiconosConvexHull>,
             siconos::collision::SiconosShape>(collision_module, "SiconosConvexHull")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::arg("vertices"));  // TODO : WARNING with vertices argument, it is shard ptr

  py::class_<siconos::collision::SiconosDisk, std::shared_ptr<siconos::collision::SiconosDisk>,
             siconos::collision::SiconosShape>(collision_module, "SiconosDisk")
      .def(py::init<float>(), py::arg("radius"))
      .def("radius", &siconos::collision::SiconosDisk::radius)
      .def("setRadius", &siconos::collision::SiconosDisk::setRadius);

  py::class_<siconos::collision::SiconosBox2d,
             std::shared_ptr<siconos::collision::SiconosBox2d>,
             siconos::collision::SiconosShape>(collision_module, "SiconosBox2d")
      .def(py::init<double, double>(), py::arg("width"), py::arg("height"));

  py::class_<siconos::collision::SiconosConvexHull2d,
             std::shared_ptr<siconos::collision::SiconosConvexHull2d>,
             siconos::collision::SiconosShape>(collision_module, "SiconosConvexHull2d")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::arg("vertices"));  // TODO : WARNING with vertices argument, it is shard ptr

  py::class_<siconos::collision::SiconosContactor,
             std::shared_ptr<siconos::collision::SiconosContactor>>(collision_module,
                                                                    "SiconosContactor")
      .def(py::init<std::shared_ptr<siconos::collision::SiconosShape>,
                    Eigen::Ref<siconos::algebra::SiconosVector>, int>(),
           py::arg("shape"), py::arg("offset"), py::arg("collision_group"));

  py::class_<siconos::collision::SiconosContactorSet,
             std::shared_ptr<siconos::collision::SiconosContactorSet>>(collision_module,
                                                                       "SiconosContactorSet")
      .def(py::init<>())
      .def("append",
           py::overload_cast<std::shared_ptr<siconos::collision::SiconosContactor>>(
               &siconos::collision::SiconosContactorSet::append),
           py::arg("contactor"));

  py::class_<siconos::collision::SiconosMesh, std::shared_ptr<siconos::collision::SiconosMesh>,
             siconos::collision::SiconosShape>(collision_module, "SiconosMesh");

  py::class_<siconos::collision::SiconosHeightMap,
             std::shared_ptr<siconos::collision::SiconosHeightMap>,
             siconos::collision::SiconosShape>(collision_module, "SiconosHeightMap");

  py::class_<siconos::collision::native::SpaceFilter,
             std::shared_ptr<siconos::collision::native::SpaceFilter>>(collision_module,
                                                                       "SpaceFilter");

  py::class_<siconos::collision::native::bodies::Circle,
             std::shared_ptr<siconos::collision::native::bodies::Circle>>(collision_module,
                                                                          "Circle")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::arg("radius"), py::arg("mass"), py::arg("position"), py::arg("velocity"));

  py::class_<siconos::collision::native::bodies::Disk,
             std::shared_ptr<siconos::collision::native::bodies::Disk>>(collision_module,
                                                                        "Disk")
      .def(py::init<double, double, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::arg("radius"), py::arg("mass"), py::arg("position"), py::arg("velocity"));

  py::class_<siconos::collision::bullet::SiconosBulletOptions,
             std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions>>(
      bullet_module, "SiconosBulletOptions")
      .def(py::init<>());

  py::class_<siconos::collision::SiconosCollisionManager,
             std::shared_ptr<siconos::collision::SiconosCollisionManager>,
             siconos::simulation::InteractionManager>(collision_module,
                                                      "SiconosCollisionManager")
      .def("insertNonSmoothLaw",
           &siconos::collision::SiconosCollisionManager::insertNonSmoothLaw);

  py::class_<siconos::collision::bullet::SiconosBulletCollisionManager,
             std::shared_ptr<siconos::collision::bullet::SiconosBulletCollisionManager>,
             siconos::collision::SiconosCollisionManager>(bullet_module,
                                                          "SiconosBulletCollisionManager")
      .def(py::init<std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions>>(),
           py::arg("options"))
      .def("addStaticBody",
           &siconos::collision::bullet::SiconosBulletCollisionManager::addStaticBody,
           py::arg("cs"), py::arg("position"), py::arg("number"));

  py::class_<siconos::collision::StaticBody, std::shared_ptr<siconos::collision::StaticBody>>(
      collision_module, "StaticBody")
      .def(py::init<>());

  py::class_<siconos::collision::RigidBodyDS, std::shared_ptr<siconos::collision::RigidBodyDS>,
             siconos::modeling::NewtonEulerDS>(m, "RigidBodyDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>, double,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("position"),
           py::arg("velocity"), py::arg("mass"), py::arg("inertia"))
      .def("setUseContactorInertia", &siconos::collision::RigidBodyDS::setUseContactorInertia,
           py::arg("useContactorInertia"))
      .def("setContactors", &siconos::collision::RigidBodyDS::setContactors,
           py::arg("contactors"));

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