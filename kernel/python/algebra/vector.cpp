#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SiconosVector.hpp"
namespace py = pybind11;

void init_siconos_vector(py::module &m)
{
  py::class_<siconos::algebra::SiconosVector,
             std::shared_ptr<siconos::algebra::SiconosVector>>(m, "SiconosVector",
                                                               py::buffer_protocol())
      .def_buffer([](siconos::algebra::SiconosVector &vec) -> py::buffer_info {
        return py::buffer_info(vec.getArray(),   /* Pointer to buffer */
                               sizeof(double), /* Size of one scalar */
                               py::format_descriptor<double>::format(), /* Python struct-style
                                                                           format descriptor */
                               1,               /* Number of dimensions */
                               {vec.size()},      /* Buffer dimensions */
                               {sizeof(double)} /* Strides (in bytes) for each index */
        );
      })
      .def(py::init<>())
      .def(py::init<unsigned int, siconos::algebra::UblasType>(), py::arg("size"),
           py::arg("storage_type")=siconos::algebra::UblasType::DENSE)
      .def("size", py::overload_cast<>(&siconos::algebra::SiconosVector::size, py::const_))
      .def("__repr__", [](const siconos::algebra::SiconosVector &a) {
        a.display();
        return "siconos.algebra.SiconosVector \n";
      });
      // .def("print",
      //      py::overload_cast<>(&siconos::algebra::SiconosVector::display, py::const_));
}



