/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "OccBody.hpp"
#include "OccSpaceFilter.hpp"
#include "OccTimeStepping.hpp"
#include "TimeDiscretisation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_occ, m) {
  m.doc() = "Siconos mechanics occ (OpenCascade) module";

  m.def("example_function", []() { return "This is an example function in mechanics.occ"; });

  py::module_ modeling_module = py::module_::import("siconos.modeling");  // For NewtonEulerDS

  py::module_ simulation_module =
      py::module_::import("siconos.simulation");  // For TimeStepping

  py::module_ mechanics_module = py::module_::import("siconos.mechanics");  // For SpaceFilter

  py::class_<siconos::mechanics::occ::OccSpaceFilter,
             std::shared_ptr<siconos::mechanics::occ::OccSpaceFilter>,
             siconos::collision::native::SpaceFilter>(m, "OccSpaceFilter")
      .def(py::init<>());

  py::class_<siconos::mechanics::occ::OccTimeStepping,
             std::shared_ptr<siconos::mechanics::occ::OccTimeStepping>,
             siconos::simulation::TimeStepping>(m, "OccTimeStepping")
      .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation>>(),
           py::arg("nsds"), py::arg("td"));

  py::class_<siconos::mechanics::occ::OccBody,
             std::shared_ptr<siconos::mechanics::occ::OccBody>,
             siconos::modeling::NewtonEulerDS>(m, "OccBody")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>, double,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 5>(), py::arg("q0"), py::arg("twist0"),
           py::arg("mass"), py::arg("inertia"));
}