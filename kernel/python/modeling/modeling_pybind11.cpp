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

#include "DynamicalSystem.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothLaw.hpp"
#include "Relation.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_dynamical_systems(py::module_ &m);
void wrap_nonsmoothlaws(py::module_ &m);
void wrap_relations(py::module_ &m);

std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactions(
    std::shared_ptr<siconos::graphs::InteractionsGraph> dsg) {
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> r =
      std::vector<std::shared_ptr<siconos::modeling::Interaction>>();
  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = dsg->vertices(); vi != viend; ++vi) {
    r.push_back(dsg->bundle(*vi));
  };
  return r;
};

PYBIND11_MODULE(modeling, m) {
     // py::module_ simulation_module = py::module_::import("siconos.simulation");

  // Optional docstring
  m.doc() = "Siconos modeling library";

  wrap_dynamical_systems(m);
  wrap_nonsmoothlaws(m);
  wrap_relations(m);

  // CLASSES with no Derived classes
  py::class_<siconos::modeling::Interaction, std::shared_ptr<siconos::modeling::Interaction>>(
      m, "Interaction")
      .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothLaw>,
                    std::shared_ptr<siconos::modeling::Relation>>())
      .def("lambda_python", &siconos::modeling::Interaction::lambda_python,
           py::return_value_policy::reference_internal)
      .def("y", &siconos::modeling::Interaction::y_python,
           py::return_value_policy::reference_internal)
      .def("relation", &siconos::modeling::Interaction::relation,
           py::return_value_policy::reference_internal);

  m.def("interactions", &interactions, py::arg("graph"),
        "Return a list of Interaction objects from an InteractionsGraph");

  py::class_<siconos::modeling::NonSmoothDynamicalSystem,
             std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>>(
      m, "NonSmoothDynamicalSystem")
      .def(py::init<double, double>())
      .def("insertDynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::insertDynamicalSystem)
      .def("getNumberOfDS",
           &siconos::modeling::NonSmoothDynamicalSystem::getNumberOfDS,
	   "get the number of DS contained in the NonSmoothDynamicalSystem")
      .def("dynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystem,
	   "get the Ds number nb")
      .def("removeDynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::removeDynamicalSystem,
	   py::arg("ds"),
	   "remove the Ds")

      .def("dynamicalSystemsVector",
           &siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystemsVector,
	   "get the array of all DSs contained in the NonSmoothDynamicalSystem")

      .def("link", &siconos::modeling::NonSmoothDynamicalSystem::link,
           "link an interaction to two dynamical systems", py::arg("inter"), py::arg("ds1"),
           py::arg("ds2") = std::shared_ptr<siconos::modeling::DynamicalSystem>())

      .def("setTitle", &siconos::modeling::NonSmoothDynamicalSystem::setTitle, "set DS title")

      .def("topology", &siconos::modeling::NonSmoothDynamicalSystem::topology,
           "display the topology of the system")

      .def("setName",
           py::overload_cast<std::shared_ptr<siconos::modeling::DynamicalSystem>,
                             const std::string &>(
               &siconos::modeling::NonSmoothDynamicalSystem::setName),
           "set DS name")
      .def("setName",
           py::overload_cast<std::shared_ptr<siconos::modeling::Interaction>,
                             const std::string &>(
               &siconos::modeling::NonSmoothDynamicalSystem::setName),
           "set Interaction name")
      .def("displayDynamicalSystems",
           &siconos::modeling::NonSmoothDynamicalSystem::displayDynamicalSystems,
           "Print all dynamical systems infos")
      .def("__repr__", [](const siconos::modeling::NonSmoothDynamicalSystem &a) {
        a.display();
        return "\n";
      });
}
