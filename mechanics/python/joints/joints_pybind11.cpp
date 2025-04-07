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
#include "SiconosJoints.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_joints, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");

  m.doc() = "Siconos mechanics.joints module";

  // Base, abstract, class
  py::class_<siconos::joints::NewtonEulerJointR,
             std::shared_ptr<siconos::joints::NewtonEulerJointR>,
             siconos::modeling::NewtonEulerR>(m, "NewtonEulerJointR")
      .def("setAbsolute", &siconos::joints::NewtonEulerJointR::setAbsolute,
           "To set the absolute reference frame for the joint");

  py::class_<siconos::joints::CylindricalJointR,
             std::shared_ptr<siconos::joints::CylindricalJointR>,
             siconos::joints::NewtonEulerJointR>(m, "CylindricalJointR")
      .def(py::init<>())
      .def("setPoint", &siconos::joints::CylindricalJointR::setPoint,
           "To define a point in the joint")
      .def("setAxis", &siconos::joints::CylindricalJointR::setAxis,
           "To define an axis in the joint")
      .def("setBasePositions", &siconos::joints::CylindricalJointR::setBasePositions,
           "To define the base position of the joint")
      .def("numberOfConstraints", &siconos::joints::CylindricalJointR::numberOfConstraints,
           "To get the number of constraints in the joint");

  py::class_<siconos::joints::PrismaticJointR,
             std::shared_ptr<siconos::joints::PrismaticJointR>,
             siconos::joints::NewtonEulerJointR>(m, "PrismaticJointR")
      .def(py::init<const Eigen::Ref<siconos::algebra::SiconosVector3>&, bool,
                    std::shared_ptr<siconos::modeling::NewtonEulerDS>,
                    std::shared_ptr<siconos::modeling::NewtonEulerDS>>(),
           py::arg("axis"), py::arg("absoluteRef"), py::arg("ds1") = nullptr,
           py::arg("ds2") = nullptr);

  py::class_<siconos::joints::CouplerJointR, std::shared_ptr<siconos::joints::CouplerJointR>,
             siconos::joints::NewtonEulerJointR>(m, "CouplerJointR")
      .def(py::init<std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int,
                    std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int, double,
                    std::optional<siconos::algebra::SiconosVector>, unsigned int,
                    std::optional<siconos::algebra::SiconosVector>, unsigned int>(),
           py::arg("joint1"), py::arg("dof1"), py::arg("joint2"), py::arg("dof2"),
           py::arg("ratio"), py::arg("ref1") = std::nullopt, py::arg("ref1_index") = 0,
           py::arg("ref2") = std::nullopt, py::arg("ref2_index") = 0)
      .def("setBasePositions", &siconos::joints::CouplerJointR::setBasePositions,
           "To define the base position of the joint");
}
