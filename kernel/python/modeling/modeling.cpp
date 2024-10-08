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
// #include <pybind11/stl.h>  // Pour permettre la conversion entre std::vector et les objets
// Python comme les listes

#include <functional>
#include <memory>
#include <span>

#include "dynamical_systems_wrapper.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(modeling, m) {
  // Optional docstring
  m.doc() = "Siconos modeling library";

  py::class_<siconos::modeling::DynamicalSystem,
             std::shared_ptr<siconos::modeling::DynamicalSystem>>(m, "DynamicalSystem");

  py::class_<siconos::modeling::SecondOrderDS,
             std::shared_ptr<siconos::modeling::SecondOrderDS>,
             siconos::modeling::DynamicalSystem>(m, "SecondOrderDS")
      .def("p", &siconos::modeling::SecondOrderDS::p_read,
           py::return_value_policy::reference_internal);

  py::class_<siconos::modeling::LagrangianDS, std::shared_ptr<siconos::modeling::LagrangianDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"))

      .def("q", &siconos::modeling::LagrangianDS::q_python,
           py::return_value_policy::reference_internal)
      .def("velocity", &siconos::modeling::LagrangianDS::velocity_python,
           py::return_value_policy::reference_internal)

      .def("setConstantFext", &siconos::modeling::LagrangianDS::setConstantFext,
           py::keep_alive<1, 2>(), "To define a constant external forces vector")

      .def(
          "setComputeFextFunction",
          [](siconos::modeling::LagrangianDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")

      .def("computeFext", &siconos::modeling::LagrangianDS::computeFext,
           "compute external forces")

      .def("fext", &siconos::modeling::LagrangianDS::fext, "current values of external forces")

      .def(
          "setComputeFintFunction",
          [](siconos::modeling::LagrangianDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFintFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
                    const Eigen::Ref<const siconos::algebra::SiconosVector> &position,
                    double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(velocity, position, time,
                    result);  // Call python func with a memory view ...
                });
          },

          "How to compute internal forces")
      .def("computeFint", &siconos::modeling::LagrangianDS::computeFint,
           "compute internal forces")
      .def_property_readonly("fint", &siconos::modeling::LagrangianDS::fint,
                             "internal forces current values")
      //      .def_property_readonly("mass", &siconos::modeling::SecondOrderDS::mass);
      .def("mass", &siconos::modeling::LagrangianDS::mass, "mass matrix");

  py::class_<siconos::modeling::LagrangianLinearTIDS,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>,
             siconos::modeling::LagrangianDS>(m, "LagrangianLinearTIDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("q0"), py::arg("v0"),
           py::arg("M"));
  // .def("__repr__", [](const siconos::modeling::LagrangianLinearTIDS &a) {
  //     (a.display());
  //   return "\n";
  // });

  py::class_<siconos::modeling::NewtonImpactNSL,
             std::shared_ptr<siconos::modeling::NewtonImpactNSL>>(m, "NewtonImpactNSL")
      .def(py::init<double>());

  py::class_<siconos::modeling::Relation, std::shared_ptr<siconos::modeling::Relation>>(
      m, "Relation");
  //     .def(py::init<siconos::modeling::RelationSubType,
  //                   siconos::modeling::RelationSubType>());

  py::class_<siconos::modeling::LagrangianLinearTIR,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIR>,
             siconos::modeling::Relation>(m, "LagrangianLinearTIR")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosMatrix> &>())
      .def("display", &siconos::modeling::LagrangianLinearTIR::display)
      .def("__repr__", [](const siconos::modeling::LagrangianLinearTIR &a) {
        a.display();
        return "\n";
      });

  py::class_<siconos::modeling::Interaction, std::shared_ptr<siconos::modeling::Interaction>>(
      m, "Interaction")
      .def(py::init<std::shared_ptr<siconos::modeling::NewtonImpactNSL>,
                    std::shared_ptr<siconos::modeling::Relation>>())
      .def("lambda_python", &siconos::modeling::Interaction::lambda_python,
           py::return_value_policy::reference_internal);

  py::class_<siconos::modeling::NonSmoothDynamicalSystem,
             std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>>(
      m, "NonSmoothDynamicalSystem")
      .def(py::init<double, double>())
      .def("insertDynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::insertDynamicalSystem)

      .def("link", &siconos::modeling::NonSmoothDynamicalSystem::link,
           "link an interaction to two dynamical systems", py::arg("inter"), py::arg("ds1"),
           py::arg("ds2") = std::shared_ptr<siconos::modeling::DynamicalSystem>())

      .def("__repr__", [](const siconos::modeling::NonSmoothDynamicalSystem &a) {
        a.display();
        return "\n";
      });
}
