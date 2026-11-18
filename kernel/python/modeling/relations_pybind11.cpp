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
#include "Lagrangian2dNormalVelocityR.hpp"
#include "Lagrangian2dVelocityAngularVelocityR.hpp"
#include "Lagrangian2dVelocityR.hpp"
#include "LagrangianLinearTIR.hpp"
#include "LagrangianR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonEulerVelocityAngularVelocityR.hpp"
#include "NewtonEulerVelocityR.hpp"

namespace py = pybind11;

// PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_relations(py::module_& m) {
  // ============================ Relation BASE CLASS ===================================
  auto relpy = py::class_<siconos::modeling::Relation, py::smart_holder>(m, "Relation");

  /** std::shared_ptr<siconos::modeling::Relation>*/

  // ============================ LagrangianR CLASS ===================================
  auto lagrpy = py::class_<siconos::modeling::LagrangianR, siconos::modeling::Relation,
                           py::smart_holder>(m, "LagrangianR");

  // ============================ LagrangianLinearTIR CLASS ==============================
  py::class_<siconos::modeling::LagrangianLinearTIR, siconos::modeling::LagrangianR,
             py::smart_holder>(m, "LagrangianLinearTIR")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>>(), py::keep_alive<1, 2>())
      .def("display", &siconos::modeling::LagrangianLinearTIR::display)
      .def("__repr__", [](const siconos::modeling::LagrangianLinearTIR& a) {
        a.display();
        return "\n";
      });

  auto lscler_py =
      py::class_<siconos::modeling::LagrangianScleronomousR, siconos::modeling::LagrangianR,
                 py::smart_holder>(m, "LagrangianScleronomousR");

  py::class_<siconos::modeling::Lagrangian2dVelocityR,
             siconos::modeling::LagrangianScleronomousR, py::smart_holder>(
      m, "Lagrangian2dVelocityR")
      .def("distance", &siconos::modeling::Lagrangian2dVelocityR::distance,
           "distance between contact points 1 and 2, with sign according to normal ");

  py::class_<siconos::modeling::Lagrangian2dVelocityAngularVelocityR,
             siconos::modeling::LagrangianScleronomousR, py::smart_holder>(
      m, "Lagrangian2dVelocityAngularVelocityR")
      .def("distance", &siconos::modeling::Lagrangian2dVelocityAngularVelocityR::distance,
           "distance between contact points 1 and 2, with sign according to normal ");
  py::class_<siconos::modeling::Lagrangian2dNormalVelocityR,
             siconos::modeling::LagrangianScleronomousR, py::smart_holder>(
      m, "Lagrangian2dNormalVelocityR")
      .def("distance", &siconos::modeling::Lagrangian2dNormalVelocityR::distance,
           "distance between contact points 1 and 2, with sign according to normal ");

  // ============================ FirstOrderR CLASS ==============================
  auto forpy = py::class_<siconos::modeling::FirstOrderR, siconos::modeling::Relation,
                          py::smart_holder>(m, "FirstOrderR");

  // ============================ FirstOrderLinearR CLASS ==============================
  py::class_<siconos::modeling::FirstOrderLinearR, siconos::modeling::FirstOrderR,
             py::smart_holder>(m, "FirstOrderLinearR")
      .def(py::init<>())
      .def("setConstantBAlias", &siconos::modeling::FirstOrderLinearR::setConstantB,
           py::keep_alive<1, 2>(), "To define a constant B operator")
      .def("setConstantCAlias", &siconos::modeling::FirstOrderLinearR::setConstantC,
           py::keep_alive<1, 2>(), "To define a constant C operator")
      .def("__repr__", [](const siconos::modeling::FirstOrderLinearR& a) {
        a.display();
        return "\n";
      });

  // ============================ FirstOrderLinearTIR CLASS ==============================
  py::class_<siconos::modeling::FirstOrderLinearTIR, siconos::modeling::FirstOrderR,
             py::smart_holder>(m, "FirstOrderLinearTIR")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix>,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
      .def("setConstantD", &siconos::modeling::FirstOrderLinearTIR::setConstantD,
           py::keep_alive<1, 2>(), "To define a constant D operator")
      .def("display", &siconos::modeling::FirstOrderLinearTIR::display)
      .def("__repr__", [](const siconos::modeling::FirstOrderLinearTIR& a) {
        a.display();
        return "\n";
      });

  // ============================ NewtonEulerR CLASS ==============================
  py::class_<siconos::modeling::NewtonEulerR, siconos::modeling::Relation, py::smart_holder>(
      m, "NewtonEulerR")
      // std::shared_ptr<siconos::modeling::NewtonEulerR>,
      .def(py::init<>())
      .def("haseVector", &siconos::modeling::NewtonEulerR::haseVector)
      // .def("display", &siconos::modeling::NewtonEulerR::display)
      // .def("__repr__", [](const siconos::modeling::NewtonEulerR &a) {
      //   a.display();
      //   return "\n";
      // })
      ;

  py::class_<siconos::modeling::NewtonEulerNormalVelocityR, siconos::modeling::NewtonEulerR,
             py::smart_holder>(m, "NewtonEulerNormalVelocityR")
      .def("distance", &siconos::modeling::NewtonEulerNormalVelocityR::distance,
           "distance between contact points 1 and 2, with sign according to normal ");

  auto ne3dR = py::class_<siconos::modeling::NewtonEulerVelocityR,
                          siconos::modeling::NewtonEulerNormalVelocityR, py::smart_holder>(
      m, "NewtonEulerVelocityR");
  auto ne5dR = py::class_<siconos::modeling::NewtonEulerVelocityAngularVelocityR,
                          siconos::modeling::NewtonEulerNormalVelocityR, py::smart_holder>(
      m, "NewtonEulerVelocityAngularVelocityR");

  // Former class names, kept for one release (deprecated).
  m.attr("Lagrangian2d1DR") = m.attr("Lagrangian2dNormalVelocityR");
  m.attr("Lagrangian2d2DR") = m.attr("Lagrangian2dVelocityR");
  m.attr("Lagrangian2d3DR") = m.attr("Lagrangian2dVelocityAngularVelocityR");
  m.attr("NewtonEuler1DR") = m.attr("NewtonEulerNormalVelocityR");
  m.attr("NewtonEuler3DR") = m.attr("NewtonEulerVelocityR");
}
