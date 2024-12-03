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
#include "FirstOrderLinearTIR.hpp"
#include "LagrangianLinearTIR.hpp"
namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_relations(py::module_ &m) {
  py::class_<siconos::modeling::Relation, std::shared_ptr<siconos::modeling::Relation>>(
      m, "Relation");

  py::class_<siconos::modeling::LagrangianLinearTIR, siconos::modeling::Relation,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIR>>(m, "LagrangianLinearTIR")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix> &>())
      .def("display", &siconos::modeling::LagrangianLinearTIR::display)
      .def("__repr__", [](const siconos::modeling::LagrangianLinearTIR &a) {
        a.display();
        return "\n";
      });

  // FirstOrderR
  py::class_<siconos::modeling::FirstOrderR, siconos::modeling::Relation,
             std::shared_ptr<siconos::modeling::FirstOrderR>>(m, "FirstOrderR");

  // FirstOrderLinearTIR
  py::class_<siconos::modeling::FirstOrderLinearTIR, siconos::modeling::FirstOrderR,
             std::shared_ptr<siconos::modeling::FirstOrderLinearTIR>>(m,
                                                                      "FirstOrderLinearTIR");
}
