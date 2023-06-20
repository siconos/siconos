#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "NewtonImpactNSL.hpp"

namespace py = pybind11;

void init_nonsmoothlaws(py::module &m)
{
  // Relation base class
  py::class_<siconos::modeling::NonSmoothLaw, std::shared_ptr<siconos::modeling::NonSmoothLaw>>(
      m, "NonSmoothLaw")
    .def_property_readonly("size", &siconos::modeling::NonSmoothLaw::size);


  // NewtonImpactNSL
  py::class_<siconos::modeling::NewtonImpactNSL, siconos::modeling::NonSmoothLaw,
             std::shared_ptr<siconos::modeling::NewtonImpactNSL>>(m, "NewtonImpactNSL")
      // nsl(size, e)
      .def(py::init<unsigned int, double>(), py::arg("size") = 1, py::arg("e") = 0.)
      // Access to restitution coefficient
      .def_property("e", &siconos::modeling::NewtonImpactNSL::e,
                    &siconos::modeling::NewtonImpactNSL::setE)
      // print
      .def("__repr__", [](const siconos::modeling::NewtonImpactNSL &a) {
        a.display();
        return "\n";
      });
}
