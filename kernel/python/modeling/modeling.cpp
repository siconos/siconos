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

#include <memory>

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
      .def("p", &siconos::modeling::SecondOrderDS::p_python,
           py::return_value_policy::reference_internal);

  py::class_<siconos::modeling::LagrangianDS, std::shared_ptr<siconos::modeling::LagrangianDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianDS");

  py::class_<siconos::modeling::LagrangianLinearTIDS,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianLinearTIDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector> &,
                    Eigen::Ref<siconos::algebra::SiconosVector> &,
                    Eigen::Ref<siconos::algebra::SiconosMatrix> &>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive as
                                    // long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("q0"), py::arg("v0"),
           py::arg("M"))
      .def("setConstantFExt",
           static_cast<void (siconos::modeling::LagrangianLinearTIDS::*)(
               Eigen::Ref<siconos::algebra::SiconosVector> &)>(
               &siconos::modeling::LagrangianLinearTIDS::setConstantFExt),
           py::keep_alive<1, 2>())
      .def("setComputeFExtFunction",
           static_cast<void (siconos::modeling::LagrangianLinearTIDS::*)(
               const std::string &pluginPath, const std::string &functionName)>(
               &siconos::modeling::LagrangianLinearTIDS::setComputeFExtFunction))

      .def("q", &siconos::modeling::LagrangianDS::q_python,
           py::return_value_policy::reference_internal)
      .def("velocity", &siconos::modeling::LagrangianDS::velocity_python,
           py::return_value_policy::reference_internal)
      .def_property_readonly("mass", &siconos::modeling::LagrangianDS::mass_python);
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
