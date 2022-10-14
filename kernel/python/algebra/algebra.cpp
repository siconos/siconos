#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <span>
#include <vector>

#include "SiconosAlgebraTypes.hpp"

namespace py = pybind11;

void init_siconos_vector(py::module&);
void init_siconos_matrix(py::module&);

// https://github.com/pybind/pybind11/issues/1042
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

// helper function to avoid making a copy when returning a py::array_t
// author: https://github.com/YannickJadoul
// source: https://github.com/pybind/pybind11/issues/1042#issuecomment-642215028
template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray(Sequence&& seq)
{
  auto size = seq.size();
  auto data = seq.data();
  std::unique_ptr<Sequence> seq_ptr = std::make_unique<Sequence>(std::move(seq));
  auto capsule = py::capsule(seq_ptr.get(), [](void* p) {
    std::unique_ptr<Sequence>(reinterpret_cast<Sequence*>(p));
  });
  seq_ptr.release();
  return py::array(size, data, capsule);
}


template <typename T> std::vector<T> make_vector(std::size_t size) {
  std::vector<T> v(size, 0);
  std::iota(v.begin(), v.end(), 0);
  return v;
}

static py::array_t<double> vector_as_array_nocopy(std::size_t size)
{
  auto temp_vector = make_vector<double>(size);
  
  return as_pyarray(std::move(temp_vector));
}
static py::array_t<double> vector_as_array(std::size_t size)
{
  auto temp_vector = make_vector<double>(size);
  return as_pyarray(std::move(temp_vector));
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

  m.def("vector_as_array_nocopy", &vector_as_array_nocopy,
        "Returns a vector of 16-bit ints as a NumPy array without making a "
        "copy of the data");

  m.def("vector_as_array", &vector_as_array, "");

  
  // m.def("numpy_to_vec", [](py::array q0) {
  //   py::print("Np conversion");
  //   // std::vector<double> v = {};
  //   auto tmp = toSpan<double>(q0);
  //   using vector_adaptor =
  //       boost::numeric::ublas::vector<double, shallow_array_adaptor<double> >;

  //   // double a[size];
  //   vector_adaptor v(shallow_array_adaptor<double>(tmp.size(), tmp.data()));  //
  //   vec.data()));

  //   // v.assign(tmp.begin(), tmp.end());
  //   py::print("vector size ... ", v.size());
  //   // m.cast<std::vector<double>>(q0);
  //   return v;  // m.cast<Eigen::Ref<Eigen::MatrixXd>>()(1, 0);
  // });
}
