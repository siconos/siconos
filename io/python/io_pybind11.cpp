// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "MechanicsIO.hpp"
#include "NonSmoothDynamicalSystem.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_io, m) {
  m.doc() = "Siconos io module";

  py::class_<siconos::io::MechanicsIO, py::smart_holder>(m, "MechanicsIO")
      .def(py::init<>())
      .def("contactInfo", &siconos::io::MechanicsIO::contactInfo)
      .def("positions", &siconos::io::MechanicsIO::positions, py::return_value_policy::move)
      .def("velocities", &siconos::io::MechanicsIO::velocities, py::return_value_policy::move)
      .def("contactPoints", &siconos::io::MechanicsIO::contactPoints,
           py::return_value_policy::move)
      .def("contactContactWork", &siconos::io::MechanicsIO::contactContactWork,
           py::return_value_policy::move, py::arg("nsds"), py::arg("index_Set") = 1,
           py::arg("omega") = 0.5, py::arg("tol") = 1e-8,
           "return the dissipation values  of all contact points as a matrix where"
           "each row corresponds to a contact, row[i] = id, normal contact work, "
           "tangent contact work, friction dissipation, contact status")
      .def("contactInternalVariables", &siconos::io::MechanicsIO::contactInternalVariables,
           py::return_value_policy::move, py::arg("nsds"), py::arg("index_Set") = 1,
           "return internal variables of all contact points as a matrix where"
           "each row corresponds to a contact");
}
