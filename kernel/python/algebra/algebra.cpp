#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <span>

#include "SiconosAlgebraTypes.hpp"

namespace py = pybind11;

// https://github.com/pybind/pybind11/issues/1042

void init_siconos_vector(py::module&);
void init_siconos_matrix(py::module&);

/**
 * \brief Returns span<T> from py:array_T<T>. Efficient as zero-copy.
 * \tparam T Type.
 * \param passthrough Numpy array.
 * \return Span<T> that with a clean and safe reference to contents of Numpy array.
 */
template <class T>
inline std::span<T> toSpan(const py::array_t<T>& passthrough)
{
  py::buffer_info passthroughBuf = passthrough.request();
  if (passthroughBuf.ndim != 1) {
    throw std::runtime_error("Error. Number of dimensions must be one");
  }
  size_t length = passthroughBuf.shape[0];
  T* passthroughPtr = static_cast<T*>(passthroughBuf.ptr);
  std::span<T> passthroughSpan(passthroughPtr, length);
  return passthroughSpan;
}


PYBIND11_MODULE(algebra, m)
{
  // Optional docstring
  m.doc() = "Siconos algebra library";

  // Export storage types enum
  py::enum_<siconos::algebra::UblasType>(m, "UblasType")
      .value("block", siconos::algebra::UblasType::BLOCK)
      .value("dense", siconos::algebra::UblasType::DENSE)
      .value("sparse", siconos::algebra::UblasType::SPARSE);

  init_siconos_vector(m);

  init_siconos_matrix(m);

  m.def("numpy_to_span", [](py::array q0) {
    py::print("Np conversion");
    std::vector<double> v = {};
    auto tmp = toSpan<double>(q0);
    // v.assign(tmp.begin(), tmp.end());
    py::print("vector size ... ", v.size());
    // m.cast<std::vector<double>>(q0);
    return true;  // m.cast<Eigen::Ref<Eigen::MatrixXd>>()(1, 0);
  });
  m.def("numpy_to_vec", [](py::array q0) {
    py::print("Np conversion");
    std::vector<double> v = {};
    auto tmp = toVec<double>(q0);
    // v.assign(tmp.begin(), tmp.end());
    py::print("vector size ... ", v.size());
    // m.cast<std::vector<double>>(q0);
    return true;  // m.cast<Eigen::Ref<Eigen::MatrixXd>>()(1, 0);
  });
}
