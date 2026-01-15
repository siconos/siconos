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
#include <pybind11/stl.h>

#include "BoundaryCondition.hpp"
#include "FixedBC.hpp"
#include "HarmonicBC.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

using Indices = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
void wrap_boundaryconditions(py::module_ &m) {
  py::class_<siconos::modeling::BoundaryCondition,
             std::shared_ptr<siconos::modeling::BoundaryCondition>>(m, "BoundaryCondition")
      .def_property_readonly("size", &siconos::modeling::BoundaryCondition::size);

  py::class_<siconos::modeling::FixedBC, std::shared_ptr<siconos::modeling::FixedBC>,
             siconos::modeling::BoundaryCondition>(m, "FixedBC")
      .def(py::init<std::vector<int>>(), py::arg("Indices"));

  py::class_<siconos::modeling::HarmonicBC, std::shared_ptr<siconos::modeling::HarmonicBC>,
             siconos::modeling::BoundaryCondition>(m, "HarmonicBC")
      .def(py::init<std::vector<int>, double, double, double, double>(), py::arg("Indices"),
           py::arg("a"), py::arg("b"), py::arg("omega"), py::arg("phi"))
      .def(py::init<std::vector<int>, Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::keep_alive<1, 5>(),
           py::keep_alive<1, 6>(), py::arg("Indices"), py::arg("a"), py::arg("b"),
           py::arg("omega"), py::arg("phi"));
}
