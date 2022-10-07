#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "LagrangianLinearTIR.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"

namespace py = pybind11;

void init_relations(py::module &m)
{
  // Relation base class
  py::class_<siconos::modeling::Relation, std::shared_ptr<siconos::modeling::Relation>>(
      m, "Relation");

  // LagrangianR
  py::class_<siconos::modeling::LagrangianR, siconos::modeling::Relation,
             std::shared_ptr<siconos::modeling::LagrangianR>>(m, "LagrangianR");

  // LagrangianLinearTIR
  py::class_<siconos::modeling::LagrangianLinearTIR, siconos::modeling::LagrangianR,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIR>>(m, "LagrangianLinearTIR",
                                                                      py::buffer_protocol())
      // LLTIR(C,e)
      .def(py::init<std::shared_ptr<siconos::algebra::SimpleMatrix>,
                    std::shared_ptr<siconos::algebra::SiconosVector>>(),
           py::arg("C").none(false), py::arg("e").none(false))

    
      .def("__repr__",
           [](const siconos::modeling::LagrangianLinearTIR &a) {
             a.display();
             return "\n";
           });
}
