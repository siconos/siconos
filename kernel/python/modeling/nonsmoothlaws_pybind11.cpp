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

// #include <functional>
// #include <memory>
// #include <span>

#include "ComplementarityConditionNSL.hpp"
#include "EqualityConditionNSL.hpp"
#include "FremondImpactFrictionNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "RelayNSL.hpp"
// #include <pybind11/stl.h>  // Pour permettre la conversion entre std::vector et les objets
// Python comme les listes

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_nonsmoothlaws(py::module_& m) {
  py::classh<siconos::modeling::NonSmoothLaw>(m, "NonSmoothLaw")
      .def_property_readonly("size", &siconos::modeling::NonSmoothLaw::size);

  // nsl(e)
  py::classh<siconos::modeling::NewtonImpactNSL, siconos::modeling::NonSmoothLaw>(
      m, "NewtonImpactNSL")
      .def(py::init<double>())

      // nsl(size, e)
      .def(py::init<unsigned int, double>(), py::arg("size") = 1, py::arg("e") = 0.)
      // Access to restitution coefficient
      .def_property("e", &siconos::modeling::NewtonImpactNSL::e,
                    &siconos::modeling::NewtonImpactNSL::setE);

  py::classh<siconos::modeling::NewtonImpactFrictionNSL, siconos::modeling::NonSmoothLaw>(
      m, "NewtonImpactFrictionNSL")
      .def(py::init<unsigned int>(), py::arg("size") = 1)
      .def(py::init<double, double, double, unsigned int>(), py::arg("en") = 0.,
           py::arg("et") = 0., py::arg("mu") = 0., py::arg("size") = 1);

  py::classh<siconos::modeling::NewtonImpactRollingFrictionNSL,
             siconos::modeling::NonSmoothLaw>(m, "NewtonImpactRollingFrictionNSL")
      .def(py::init<unsigned int>(), py::arg("size") = 1)
      .def(py::init<double, double, double, double, unsigned int>(), py::arg("en") = 0.,
           py::arg("et") = 0., py::arg("mu") = 0., py::arg("muR") = 0., py::arg("size") = 1);

  py::classh<siconos::modeling::FremondImpactFrictionNSL, siconos::modeling::NonSmoothLaw>(
      m, "FremondImpactFrictionNSL")
      .def(py::init<unsigned int>(), py::arg("size") = 1)
      .def(py::init<double, double, double, unsigned int>(), py::arg("en") = 0.,
           py::arg("et") = 0., py::arg("mu") = 0., py::arg("size") = 1);

  py::classh<siconos::modeling::ComplementarityConditionNSL, siconos::modeling::NonSmoothLaw>(
      m, "ComplementarityConditionNSL")
      .def(py::init<unsigned int>());

  py::classh<siconos::modeling::RelayNSL, siconos::modeling::NonSmoothLaw>(m, "RelayNSL")
      .def(py::init<unsigned int, double, double>());

  py::classh<siconos::modeling::EqualityConditionNSL, siconos::modeling::NonSmoothLaw>(
      m, "EqualityConditionNSL")
      .def(py::init<unsigned int>());
}
