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

#include "MechanicsIO.hpp"
#include "NonSmoothDynamicalSystem.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_io, m) {
  m.doc() = "Siconos io module";

  py::class_<siconos::io::MechanicsIO, std::shared_ptr<siconos::io::MechanicsIO>>(
      m, "MechanicsIO")
      .def(py::init<>())
      .def("contactInfo", &siconos::io::MechanicsIO::contactInfo)
      .def("positions", &siconos::io::MechanicsIO::positions)
      .def("velocities", &siconos::io::MechanicsIO::velocities)
      .def("contactPoints", &siconos::io::MechanicsIO::contactPoints)
      .def("contactContactWork", &siconos::io::MechanicsIO::contactContactWork);
}