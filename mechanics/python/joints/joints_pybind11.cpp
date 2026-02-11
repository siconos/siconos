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

#include "NewtonEulerDS.hpp"
#include "NewtonEulerJointR.hpp"
#include "SiconosJoints.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_joints, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");

  m.doc() = "Siconos mechanics.joints module";

  // Base, abstract, class
  py::classh<siconos::joints::NewtonEulerJointR, siconos::modeling::NewtonEulerR>(
      m, "NewtonEulerJointR")
      .def("setAbsolute", &siconos::joints::NewtonEulerJointR::setAbsolute,
           "To set the absolute reference frame for the joint")
      .def("setBasePositions", &siconos::joints::NewtonEulerJointR::setBasePositions,
           "To define the base position of the joint")
      .def("numberOfConstraints", &siconos::joints::NewtonEulerJointR::numberOfConstraints,
           "To get the number of constraints in the joint")
      .def("setPoint",
           py::overload_cast<size_t, const Eigen::Ref<siconos::algebra::SiconosVector3>&>(
               &siconos::joints::NewtonEulerJointR::setPoint),
           py::arg("index"), py::arg("point"), "To define a point in the joint")
      .def("setAxis", &siconos::joints::NewtonEulerJointR::setAxis,
           "To define an axis in the joint")
      .def("projectVectorDoF", &siconos::joints::NewtonEulerJointR::projectVectorDoF,
           py::arg("v"), py::arg("q0"), py::arg("q1") = std::nullopt, py::arg("axis") = 0,
           py::arg("absoluteRef") = true,
           "Project a vector on the DoF axis at given configuration.")
      .def("setAllowSelfCollide", &siconos::joints::NewtonEulerJointR::setAllowSelfCollide,
           "To set the allowSelfCollide flag")
      .def("normalDoF", &siconos::joints::NewtonEulerJointR::normalDoF, py::arg("q0"),
           py::arg("q1") = std::nullopt, py::arg("axis") = 0, py::arg("absoluteRef") = true,
           "Return the axis of rotation as a SiconosVector3 (q passed directly)")
      .def("computehDoF", &siconos::joints::NewtonEulerJointR::computehDoF, py::arg("q1"),
           py::arg("q2") = std::nullopt, py::arg("y"), py::arg("axis") = 0,
           "To compute the DoF of the joint")
      .def("computeJachqDoF", &siconos::joints::NewtonEulerJointR::computeJachqDoF,
           py::arg("inter"), py::arg("q1"), py::arg("q2") = std::nullopt, py::arg("jachq"),
           py::arg("axis") = 0)
      .def("numberOfConstraints", &siconos::joints::NewtonEulerJointR::numberOfConstraints,
           "To get the number of constraints in the joint")
      .def("computeh",
           py::overload_cast<
               const Eigen::Ref<const siconos::algebra::SiconosVector>&,
               const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>&,
               Eigen::Ref<siconos::algebra::SiconosVector>>(
               &siconos::joints::NewtonEulerJointR::computeh),
           py::arg("q1"), py::arg("q2") = std::nullopt, py::arg("y"),
           "to compute the output y = h(q) of the Relation");

  py::classh<siconos::joints::CylindricalJointR, siconos::joints::NewtonEulerJointR>(
      m, "CylindricalJointR")
      .def(py::init<>(), "Default constructor for CylindricalJointR");

  py::classh<siconos::joints::PrismaticJointR, siconos::joints::NewtonEulerJointR>(
      m, "PrismaticJointR")
      .def(py::init<>(), "Default constructor for PrismaticJointR")
      .def("numberOfConstraints", &siconos::joints::PrismaticJointR::numberOfConstraints,
           "To get the number of constraints in the joint");

  py::classh<siconos::joints::FixedJointR, siconos::joints::NewtonEulerJointR>(m,
                                                                               "FixedJointR")
      .def(py::init<>(), "Default constructor for FixedJointR");

  py::classh<siconos::joints::KneeJointR, siconos::joints::NewtonEulerJointR>(m, "KneeJointR")
      .def(py::init<>(), "Default constructor for KneeJointR");

  py::classh<siconos::joints::PivotJointR, siconos::joints::KneeJointR>(m, "PivotJointR")
      .def(py::init<>(), "Default constructor for PivotJointR");

  py::classh<siconos::joints::CouplerJointR, siconos::joints::NewtonEulerJointR>(
      m, "CouplerJointR")
      .def(py::init<std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int,
                    std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int, double,
                    std::optional<siconos::algebra::SiconosVector>, unsigned int,
                    std::optional<siconos::algebra::SiconosVector>, unsigned int>(),
           py::arg("joint1"), py::arg("dof1"), py::arg("joint2"), py::arg("dof2"),
           py::arg("ratio"), py::arg("ref1") = std::nullopt, py::arg("ref1_index") = 0,
           py::arg("ref2") = std::nullopt, py::arg("ref2_index") = 0)
      .def(py::init<std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int,
                    std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int, double,
                    std::shared_ptr<siconos::modeling::NewtonEulerDS>, unsigned int,
                    std::shared_ptr<siconos::modeling::NewtonEulerDS>, unsigned int>(),
           py::arg("joint1"), py::arg("dof1"), py::arg("joint2"), py::arg("dof2"),
           py::arg("ratio"), py::arg("refds1"), py::arg("ref1_index") = 0, py::arg("refds2"),
           py::arg("ref2_index") = 0);

  py::classh<siconos::joints::JointStopR, siconos::modeling::NewtonEulerR>(m, "JointStopR")
      .def(py::init<std::shared_ptr<siconos::joints::NewtonEulerJointR>, double, bool,
                    unsigned int>(),
           py::arg("joint"), py::arg("pos"), py::arg("dir"), py::arg("axis") = 0);

  py::classh<siconos::joints::JointFrictionR, siconos::modeling::NewtonEulerR>(
      m, "JointFrictionR")
      .def(py::init<std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int>(),
           py::arg("joint"), py::arg("axis"));
}
