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

#include <pybind11/detail/using_smart_holder.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "BoundaryCondition.hpp"
#include "HarmonicBC.hpp"
#include "SiconosVector.hpp"

namespace py = pybind11;

// PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

using Indices = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
void wrap_boundaryconditions(py::module_ &m) {
  py::class_<siconos::modeling::BoundaryCondition, py::smart_holder>(m, "BoundaryCondition")
      .def(py::init<std::vector<siconos::algebra::Index>>(), py::arg("Indices"))
      .def(py::init<std::vector<siconos::algebra::Index>,
                    const Eigen::Ref<const siconos::algebra::SiconosVector> &>(),
           py::arg("Indices"),
           py::arg("prescribedVelocities").noconvert())  // failing rather than copying

      .def_property_readonly("size", &siconos::modeling::BoundaryCondition::size)
      .def_property_readonly(
          "velocityIndices",
          [](const siconos::modeling::BoundaryCondition &self) {
            auto span = self.velocityIndices();
            return py::array_t<siconos::algebra::Index>(span.size(),     // shape
                                                        span.data(),     // data pointer
                                                        py::cast(&self)  // to keep self alive
            );
          },
          "indices where conditions apply")

      .def_property_readonly("prescribedVelocity",
                             &siconos::modeling::BoundaryCondition::prescribedVelocity,
                             "values applied on constrained indices")
      .def("computePrescribedVelocity",
           &siconos::modeling::BoundaryCondition::computePrescribedVelocity,
           "Compute velocity values on prescribed indices");

  py::class_<siconos::modeling::HarmonicBC, siconos::modeling::BoundaryCondition,
             py::smart_holder>(m, "HarmonicBC")
      .def(py::init<std::vector<int>, double, double, double, double>(), py::arg("Indices"),
           py::arg("a"), py::arg("b"), py::arg("omega"), py::arg("phi"))
      .def(
          py::init<std::vector<int>, const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                   const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                   const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                   const Eigen::Ref<const siconos::algebra::SiconosVector> &>(),
          //     py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
          //     py::keep_alive<1, 5>(), py::keep_alive<1, 6>(),
          py::arg("Indices"), py::arg("a").noconvert(), py::arg("b").noconvert(),
          py::arg("omega").noconvert(), py::arg("phi").noconvert());
}
