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

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "LagrangianLinearTIR.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace py = pybind11;

void init_relations(py::module &m) {
  // Relation base class
  py::class_<siconos::modeling::Relation,
             std::shared_ptr<siconos::modeling::Relation>>(m, "Relation");

  // LagrangianR
  py::class_<siconos::modeling::LagrangianR, siconos::modeling::Relation,
             std::shared_ptr<siconos::modeling::LagrangianR>>(m, "LagrangianR");

  // FirstOrderR
  py::class_<siconos::modeling::FirstOrderR, siconos::modeling::Relation,
             std::shared_ptr<siconos::modeling::FirstOrderR>>(m, "FirstOrderR");
  
  // FirstOrderLinearTIR
  py::class_<siconos::modeling::FirstOrderLinearTIR, siconos::modeling::FirstOrderR,
             std::shard_ptr<siconos::modeling::FirstOrderLinearTIR>>(m, "FirstOrderLinearTIR");

  // LagrangianLinearTIR
  py::class_<siconos::modeling::LagrangianLinearTIR,
             siconos::modeling::LagrangianR,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIR>>(
      m, "LagrangianLinearTIR", py::buffer_protocol())
      // LLTIR(C,e)
      .def(py::init<std::shared_ptr<siconos::algebra::SiconosMatrix>,
                    std::shared_ptr<siconos::algebra::SiconosVector>>(),
           py::arg("C").none(false), py::arg("e").none(false))

      .def("__repr__", [](const siconos::modeling::LagrangianLinearTIR &a) {
        a.display();
        return "\n";
      });
}
