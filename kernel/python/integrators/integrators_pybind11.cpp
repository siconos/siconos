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
 */#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "EulerMoreauOSI.hpp"
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "MoreauJeanGOSI.hpp"
#include "MoreauJeanOSI.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(integrators, m) {
  // Optional docstring
  m.doc() = "Siconos one-step integrators";

  // OSI base class
  py::class_<siconos::integrators::OneStepIntegrator,
             std::shared_ptr<siconos::integrators::OneStepIntegrator>>(m, "OneStepIntegrator");

  // MoreauJean
  py::class_<siconos::integrators::MoreauJeanOSI,
             std::shared_ptr<siconos::integrators::MoreauJeanOSI>,
             siconos::integrators::OneStepIntegrator>(m, "MoreauJeanOSI")
      .def(py::init<double, double>(), py::arg("theta") = 0.5,
           py::arg("gamma") = std::numeric_limits<double>::quiet_NaN())
      .def("setConstraintActivationThreshold",
           &siconos::integrators::MoreauJeanOSI::setConstraintActivationThreshold)
      .def_property("theta",
                    &siconos::integrators::MoreauJeanOSI::theta,    // getter
                    &siconos::integrators::MoreauJeanOSI::setTheta  // setter
                    )
      .def_property("gamma",
                    &siconos::integrators::MoreauJeanOSI::gamma,    // getter
                    &siconos::integrators::MoreauJeanOSI::setGamma  // setter
                    )
      .def("__repr__", [](const siconos::integrators::MoreauJeanOSI &a) {
        a.display();
        return "\n";
      });

  py::class_<siconos::integrators::MoreauJeanGOSI,
             std::shared_ptr<siconos::integrators::MoreauJeanGOSI>,
             siconos::integrators::MoreauJeanOSI>(m, "MoreauJeanGOSI")
      .def(py::init<double, double>(), py::arg("theta") = 0.5,
           py::arg("gamma") = std::numeric_limits<double>::quiet_NaN());

  py::class_<siconos::integrators::EulerMoreauOSI,
             std::shared_ptr<siconos::integrators::EulerMoreauOSI>,
             siconos::integrators::OneStepIntegrator>(m, "EulerMoreauOSI")
      .def(py::init<double>(), py::arg("theta"))
      .def(py::init<double, double>(), py::arg("theta"), py::arg("gamma"))
      .def_property("theta",
                    &siconos::integrators::EulerMoreauOSI::theta,    // getter
                    &siconos::integrators::EulerMoreauOSI::setTheta  // setter
                    )
      .def_property("gamma",
                    &siconos::integrators::EulerMoreauOSI::gamma,    // getter
                    &siconos::integrators::EulerMoreauOSI::setGamma  // setter
                    )
      .def_property("useGamma", &siconos::integrators::EulerMoreauOSI::useGamma,
                    &siconos::integrators::EulerMoreauOSI::setUseGamma)
      .def_property("useGammaForRelation",
                    &siconos::integrators::EulerMoreauOSI::useGammaForRelation,
                    &siconos::integrators::EulerMoreauOSI::setUseGammaForRelation)

      .def("__repr__", [](const siconos::integrators::EulerMoreauOSI &a) {
        a.display();
        return "\n";
      });

  py::class_<siconos::integrators::MoreauJeanDirectProjectionOSI,
             std::shared_ptr<siconos::integrators::MoreauJeanDirectProjectionOSI>,
             siconos::integrators::MoreauJeanOSI>(m, "MoreauJeanDirectProjectionOSI")
      .def(py::init<double, double>(), py::arg("theta") = 0.5,
           py::arg("gamma") = std::numeric_limits<double>::quiet_NaN())
      .def(py::init<double>(), py::arg("theta"));
}
