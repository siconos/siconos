#include "matrix_wrapper.h"
#include "vector_wrapper.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <memory>


namespace py = pybind11;




PYBIND11_MODULE(algebra, m) {

    m.doc() = "Siconos types wrapper";

    pybind11::class_<MatrixWrapper, std::shared_ptr<MatrixWrapper>>(m, "MatrixWrapper")
        .def(pybind11::init<siconos::algebra::SiconosMatrix&>())
        .def("get_matrix", &MatrixWrapper::get_matrix, py::return_value_policy::reference);

    py::class_<MyVectorClass>(m, "MyVectorClass")
    .def(py::init<std::vector<float>&>())
    .def("print_v", &MyVectorClass::print_v);
    // .def_readwrite("contents", &MyVectorClass::contents);

    py::class_<MyEigenClass>(m, "MyEigenClass")
    .def(py::init<>())
    .def_readwrite("contents", &MyEigenClass::contents);

    pybind11::class_<VectorWrapper, std::shared_ptr<VectorWrapper>>(m, "VectorWrapper")
        .def(pybind11::init<siconos::algebra::SiconosVector&>())
        .def("get_vector", &VectorWrapper::get_vector, py::return_value_policy::reference_internal)
        .def("print_vector", &VectorWrapper::print_vector);

    pybind11::class_<MyVectorWrapper, std::shared_ptr<MyVectorWrapper>>(m, "MyVectorWrapper")
        .def(pybind11::init<MyVectorClass&>())
        .def("print_vector", &MyVectorWrapper::print_vector);

    pybind11::class_<MyEigenWrapper, std::shared_ptr<MyEigenWrapper>>(m, "MyEigenWrapper")
        .def(pybind11::init<MyEigenClass&>())
        .def("change_vector", &MyEigenWrapper::change_vector)
        .def("print_vector", &MyEigenWrapper::print_vector);
}