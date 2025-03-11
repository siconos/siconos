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

namespace py = pybind11;

PYBIND11_MODULE(_bullet, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");
  py::module_ collision_module = py::module_::import("siconos.mechanics.collision");

  m.doc() = "Siconos mechanics.collision.bullet module";

  py::class_<siconos::collision::bullet::SiconosBulletOptions,
             std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions>>(
      m, "SiconosBulletOptions")
      .def(py::init<>());

  py::class_<siconos::collision::bullet::SiconosBulletCollisionManager,
             std::shared_ptr<siconos::collision::bullet::SiconosBulletCollisionManager>,
             siconos::collision::SiconosCollisionManager>(m, "SiconosBulletCollisionManager")
      .def(py::init<std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions>>(),
           py::arg("options"))
      .def("addStaticBody",
           &siconos::collision::bullet::SiconosBulletCollisionManager::addStaticBody,
           py::arg("cs"), py::arg("position"), py::arg("number"));
}