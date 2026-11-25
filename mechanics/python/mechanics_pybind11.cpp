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
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MODULE(_mechanics, m) {
  m.doc() = "Siconos mechanics module";

  py::module_ module_simulation = py::module_::import("siconos.simulation");

  m.def("example_function", []() { return "This is an example function in mechanics"; });
}
