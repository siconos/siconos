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

#include <TopoDS_Shape.hxx>

#include "OccBody.hpp"
#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "OccContactShape.hpp"
#include "OccSpaceFilter.hpp"
#include "OccTimeStepping.hpp"
#include "OccUtils.hpp"
#include "TimeDiscretisation.hpp"

namespace py = pybind11;
void py_occ_move(siconos::mechanics::occ::OccContactShape& ocs,
                 const std::array<double, 7>& q) {
  // We do not want to expose ocs.data in python to avoid
  // confusion between pythonocc-core swig objects and ours
  // so we use this wrapper
  siconos::mechanics::occ::occ_move(ocs.data(), q);
}

PYBIND11_MODULE(_occ, m) {
  m.doc() = "Siconos mechanics occ (OpenCascade) module";

  m.def("example_function", []() { return "This is an example function in mechanics.occ"; });

  py::module_ modeling_module = py::module_::import("siconos.modeling");  // For NewtonEulerDS

  py::module_ simulation_module =
      py::module_::import("siconos.simulation");  // For TimeStepping

  py::module_ mechanics_module = py::module_::import("siconos.mechanics");  // For SpaceFilter

  auto pyocc_module = py::module::import("OCC.Core.TopoDS");

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
           py::arg("mass"), py::arg("inertia"))
      .def("addContactShape", &siconos::mechanics::occ::OccBody::addContactShape,
           "associate the body with a contact shape");

  m.def("occ_move", &py_occ_move, "Move a TopoDS_Shape using a translation and rotation array",
        py::arg("shape"), py::arg("q"));

  py::class_<siconos::mechanics::occ::OccContactShape,
             std::shared_ptr<siconos::mechanics::occ::OccContactShape>>(m, "OccContactShape")
      .def(py::init([](py::object shape_obj) {
             if (!py::hasattr(shape_obj, "this"))
               throw std::runtime_error(
                   "Expected a pythonocc-core object with a 'this' attribute");

             // Try to interact between pythonocc-core object (wrapped from swig)
             // and ours ...

             // .this is a SwigPyObject which has a C++ pointer as an int
             py::object this_attr = shape_obj.attr("this");

             // Call __long__ (Python 3) or __int__ (fallback) to get C++ address
             uintptr_t ptr_val = 0;
             if (py::hasattr(this_attr, "__int__"))
               ptr_val = this_attr.attr("__int__")().cast<uintptr_t>();
             else if (py::hasattr(this_attr, "__long__"))
               ptr_val = this_attr.attr("__long__")().cast<uintptr_t>();
             else
               throw std::runtime_error("Unable to extract C++ pointer from SWIG object");

             auto native = reinterpret_cast<TopoDS_Shape*>(ptr_val);
             if (!native) throw std::runtime_error("Null pointer received from SWIG object");

             return std::make_shared<siconos::mechanics::occ::OccContactShape>(*native);
           }),
           py::arg("shape"));

  py::class_<siconos::mechanics::occ::OccContactFace,
             std::shared_ptr<siconos::mechanics::occ::OccContactFace>,
             siconos::mechanics::occ::OccContactShape>(m, "OccContactFace")
      .def(py::init<const siconos::mechanics::occ::OccContactShape&, int>(), py::arg("shape"),
           py::arg("index"));

  py::class_<siconos::mechanics::occ::OccContactEdge,
             std::shared_ptr<siconos::mechanics::occ::OccContactEdge>,
             siconos::mechanics::occ::OccContactShape>(m, "OccContactEdge")
      .def(py::init<const siconos::mechanics::occ::OccContactShape&, int>(), py::arg("shape"),
           py::arg("index"));
}