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
#include <pybind11/stl.h>  // To allow conversion between std::vector and Python objects

#include <functional>
#include <memory>
#include <span>

#include "ReferenceClasses.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(pb11_template, m) {
  // Optional docstring
  m.doc() = "Reference for developpers - How to write a pybind11 wrapper for Siconos";

  py::class_<siconos::internal::devel_model::ClassA,
             std::shared_ptr<siconos::internal::devel_model::ClassA>>(m, "ClassA")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector> &>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::arg("vectorName1"))

      .def("setConstantVectorName2",
           &siconos::internal::devel_model::ClassA::setConstantVectorName2,
           py::keep_alive<1, 2>())

      .def(
          "setComputeVectorName2Function",
          [](siconos::internal::devel_model::ClassA &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeVectorName2Function(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")

      .def("computeVectorName2", &siconos::internal::devel_model::ClassA::computeVectorName2,
           "compute external forces")

      .def("vectorName2", &siconos::internal::devel_model::ClassA::vectorName2_view,
           "current value of external forces")

      .def("setConstantVectorNameDirect",
           &siconos::internal::devel_model::ClassA::setConstantVectorNameDirect,
           py::keep_alive<1, 2>())

      .def(
          "setComputeVectorNameDirectFunction",
          [](siconos::internal::devel_model::ClassA &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeVectorNameDirectFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);
                });
          },
          "How to compute external forces")

      .def("computeVectorNameDirect",
           &siconos::internal::devel_model::ClassA::computeVectorNameDirect,
           "compute external forces")

      .def("vectorNameDirect", &siconos::internal::devel_model::ClassA::vectorNameDirect_view,
           "current value of external forces")

      .def("setConstantVectorNameSpan",
           &siconos::internal::devel_model::ClassA::setConstantVectorNameSpan,
           py::keep_alive<1, 2>())

      .def(
          "setComputeVectorNameSpanFunction",
          [](siconos::internal::devel_model::ClassA &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeVectorNameSpanFunction([f](double val, std::span<double> span) {
              // Build numpy pointer, sharing memory with std::span. No copy !
              auto array =
                  py::array_t<double>(span.size(), span.data(), py::cast(span.data()));
              f(val, array);  // Call python func with a memory view ...
            });
          },
          "How to compute external forces")

      .def("computeVectorNameSpan",
           &siconos::internal::devel_model::ClassA::computeVectorNameSpan,
           "compute external forces")

      .def("vectorNameSpan", &siconos::internal::devel_model::ClassA::vectorNameSpan_view,
           "current value of external forces")

      .def("vectorName3", &siconos::internal::devel_model::ClassA::vectorName3_python,
           py::return_value_policy::reference_internal)

      .def("setConstantMatrixName",
           &siconos::internal::devel_model::ClassA::setConstantMatrixName,
           py::keep_alive<1, 2>())

      .def(
          "setComputeMatrixNameFunction",
          [](siconos::internal::devel_model::ClassA &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeMatrixNameFunction(
                [f](Eigen::Ref<siconos::algebra::MapVectorType> pos, double val,
                    Eigen::Ref<siconos::algebra::MapType> result) {
                  f(val, pos, result);  // Call python func with a memory view ...
                });
          },
          "How to compute xxx matrix operator")

      .def("computeMatrixName", &siconos::internal::devel_model::ClassA::computeMatrixName,
           "compute mass matrix")

      .def_property_readonly("matrixName",
                             &siconos::internal::devel_model::ClassA::matrixName_view);
}
