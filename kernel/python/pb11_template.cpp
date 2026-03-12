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
           py::arg("vector1"))

      .def("setConstantVector2", &siconos::internal::devel_model::ClassA::setConstantVector2,
           py::keep_alive<1, 2>())

      .def(
          "setComputeVector2Function",
          [](siconos::internal::devel_model::ClassA &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeVector2Function(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")

      .def("computeVector2", &siconos::internal::devel_model::ClassA::computeVector2,
           "compute external forces")

      .def("vector2", &siconos::internal::devel_model::ClassA::vector2,
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

      .def("var", &siconos::internal::devel_model::ClassA::var_read,
           py::return_value_policy::reference_internal)

      .def("setConstantMatrix1", &siconos::internal::devel_model::ClassA::setConstantMatrix1,
           py::keep_alive<1, 2>())

      .def(
          "setComputeMatrix1Function",
          [](siconos::internal::devel_model::ClassA &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeMatrix1Function([f](Eigen::Ref<siconos::algebra::MapVectorType> pos,
                                               double val,
                                               Eigen::Ref<siconos::algebra::MapType> result) {
              f(val, pos, result);  // Call python func with a memory view ...
            });
          },
          "How to compute xxx matrix operator")

      .def("computeMatrix1", &siconos::internal::devel_model::ClassA::computeMatrix1,
           "compute mass matrix")

      .def_property_readonly("matrix1", &siconos::internal::devel_model::ClassA::matrix1);
}
