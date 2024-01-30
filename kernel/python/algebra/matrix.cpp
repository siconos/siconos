#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SiconosMatrix.hpp"
#include "SimpleMatrix.hpp"
namespace py = pybind11;

// void init_siconos_matrix(py::module &m) {
//   py::class_<siconos::algebra::SiconosMatrix,
//              std::shared_ptr<siconos::algebra::SiconosMatrix>>(m,
//                                                                "SiconosMatrix");

//   py::class_<siconos::algebra::SimpleMatrix, siconos::algebra::SiconosMatrix,
//              std::shared_ptr<siconos::algebra::SimpleMatrix>>(
//       m, "SimpleMatrix", py::buffer_protocol())
//       .def_buffer([](siconos::algebra::SimpleMatrix &mat) -> py::buffer_info {
//         return py::buffer_info(
//             mat.data(),                              /* Pointer to buffer */
//             sizeof(double),                          /* Size of one scalar */
//             py::format_descriptor<double>::format(), /* Python struct-style
//                                                         format descriptor */
//             2,                                       /* Number of dimensions */
//             {mat.size(0), mat.size(1)},              /* Buffer dimensions */
//             {sizeof(double) * mat.size(1), sizeof(double)}
//             /* Strides (in bytes) for each index */
//         );
//       })
//       //.def(py::init<>())
//       .def(py::init<unsigned int, unsigned int>(), py::arg("row"),
//            py::arg("col"))
//       .def("size", &siconos::algebra::SimpleMatrix::size)
//       // .def("print",
//       // py::overload_cast<>(&siconos::algebra::SimpleMatrix::display,
//       // py::const_))
//       .def("__repr__", [](const siconos::algebra::SiconosMatrix &a) {
//         a.display();
//         return "siconos.algebra.SimpleMatrix \n";
//       });
// }
