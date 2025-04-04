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

#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "LagrangianLinearTIR.hpp"
#include "NewtonEulerR.hpp"
namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_relations(py::module_ &m) {
  // ============================ Relation BASE CLASS ===================================
  py::class_<siconos::modeling::Relation, std::shared_ptr<siconos::modeling::Relation>>(
      m, "Relation");

  // ============================ LagrangianR CLASS ===================================
  py::class_<siconos::modeling::LagrangianR, std::shared_ptr<siconos::modeling::LagrangianR>,
             siconos::modeling::Relation>(m, "LagrangianR");

  // ============================ LagrangianLinearTIR CLASS ==============================
  py::class_<siconos::modeling::LagrangianLinearTIR,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIR>,
             siconos::modeling::LagrangianR>(m, "LagrangianLinearTIR")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(), py::keep_alive<1, 2>())
      .def("display", &siconos::modeling::LagrangianLinearTIR::display)
      .def("__repr__", [](const siconos::modeling::LagrangianLinearTIR &a) {
        a.display();
        return "\n";
      });

  // ============================ FirstOrderR CLASS ==============================
  py::class_<siconos::modeling::FirstOrderR, std::shared_ptr<siconos::modeling::FirstOrderR>,
             siconos::modeling::Relation>(m, "FirstOrderR");

  // ============================ FirstOrderLinearR CLASS ==============================
  py::class_<siconos::modeling::FirstOrderLinearR,
             std::shared_ptr<siconos::modeling::FirstOrderLinearR>,
             siconos::modeling::FirstOrderR>(m, "FirstOrderLinearR")
      .def(py::init<>())
      .def("setConstantB", &siconos::modeling::FirstOrderLinearR::setConstantB,
           py::keep_alive<1, 2>(), "To define a constant B operator")
      .def("setConstantC", &siconos::modeling::FirstOrderLinearR::setConstantC,
           py::keep_alive<1, 2>(), "To define a constant C operator")
      .def("__repr__", [](const siconos::modeling::FirstOrderLinearR &a) {
        a.display();
        return "\n";
      });

  // ============================ FirstOrderLinearTIR CLASS ==============================
  py::class_<siconos::modeling::FirstOrderLinearTIR,
             std::shared_ptr<siconos::modeling::FirstOrderLinearTIR>,
             siconos::modeling::FirstOrderR>(m, "FirstOrderLinearTIR")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
      .def("setConstantD", &siconos::modeling::FirstOrderLinearTIR::setConstantD,
           py::keep_alive<1, 2>(), "To define a constant D operator")
      .def("display", &siconos::modeling::FirstOrderLinearTIR::display)
      .def("__repr__", [](const siconos::modeling::FirstOrderLinearTIR &a) {
        a.display();
        return "\n";
      });

  // ============================ NewtonEulerR CLASS ==============================
  py::class_<siconos::modeling::NewtonEulerR,
             std::shared_ptr<siconos::modeling::NewtonEulerR>,
             siconos::modeling::Relation>(m, "NewtonEulerR")
      .def(py::init<>())
      // .def("display", &siconos::modeling::NewtonEulerR::display)
      // .def("__repr__", [](const siconos::modeling::NewtonEulerR &a) {
      //   a.display();
      //   return "\n";
      // })
      ;
}
