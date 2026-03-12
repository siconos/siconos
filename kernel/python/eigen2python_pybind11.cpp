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

#include "eigen2python_pybind11.hpp"

#include "SiconosMatrix.hpp"

namespace py = pybind11;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
using namespace pybind11::literals;  // to use "_a"

/**
 * @brief convert a scipy array (sparse) to a SiconosSparseMatrix (Eigen)
 *
 * @param csc the input array
 * @return the expected sparse matrix
 *
 * @note This function creates (and so allocate memory) for a SiconosSparseMatrix
 */
siconos::algebra::SiconosSparseMatrix siconos::pybind11_utils::csc_to_eigen(py::object csc) {
  auto data =
      csc.attr("data").cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();
  auto indices = csc.attr("indices")
                     .cast<py::array_t<siconos::algebra::SparseIndex,
                                       py::array::c_style | py::array::forcecast>>();
  auto indptr = csc.attr("indptr")
                    .cast<py::array_t<siconos::algebra::SparseIndex,
                                      py::array::c_style | py::array::forcecast>>();
  auto shape =
      csc.attr("shape")
          .cast<std::pair<siconos::algebra::SparseIndex, siconos::algebra::SparseIndex>>();

  siconos::algebra::SparseIndex rows = shape.first, cols = shape.second;
  siconos::algebra::SparseIndex nnz = static_cast<siconos::algebra::SparseIndex>(data.size());

  siconos::algebra::SiconosSparseMatrix M(rows, cols);
  M.reserve(nnz);
  auto data_ptr = data.mutable_data();
  auto idx_ptr = indices.mutable_data();
  auto indptr_ptr = indptr.mutable_data();

  for (siconos::algebra::SparseIndex j = 0; j < cols; ++j) {
    auto start = indptr_ptr[j];
    auto end = indptr_ptr[j + 1];
    for (auto k = start; k < end; ++k) {
      M.insert(idx_ptr[k], j) = data_ptr[k];
    }
  }
  M.makeCompressed();
  return M;  // RVO
}

/**
 * @brief convert a scipy array (sparse) to a Map onto a SiconosSparseMatrix (Eigen)
 *
 * @param csc the input array
 * @return the expected sparse matrix (as a shared pointer)
 *
 * @note No memory allocation!
 */
std::shared_ptr<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>
siconos::pybind11_utils::csc_to_eigen_map(py::object csc) {
  auto data =
      csc.attr("data").cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();
  auto indices = csc.attr("indices")
                     .cast<py::array_t<siconos::algebra::SparseIndex,
                                       py::array::c_style | py::array::forcecast>>();
  auto indptr = csc.attr("indptr")
                    .cast<py::array_t<siconos::algebra::SparseIndex,
                                      py::array::c_style | py::array::forcecast>>();
  // auto shape = csc.attr("shape").cast<std::pair<int, int>>();
  auto shape = csc.attr("shape").cast<std::pair<ssize_t, ssize_t>>();

  // make_shared with custom deleter? not needed: we want shared_ptr<Map>
  using MapType = Eigen::Map<siconos::algebra::SiconosSparseMatrix>;
  //   auto map_ptr =
  //       std::shared_ptr<MapType>(new Eigen::Map<siconos::algebra::SiconosSparseMatrix>(
  //           shape.first, shape.second, static_cast<int>(data.size()), indptr.mutable_data(),
  //           indices.mutable_data(), data.mutable_data()));
  auto map_ptr = std::make_shared<MapType>(shape.first, shape.second, indptr.size() - 1,
                                           indptr.mutable_data(), indices.mutable_data(),
                                           data.mutable_data());

  // Note: the numpy arrays (csc.data/indices/indptr) must stay alive.
  // The binding will attach 'csc' to the Python instance (see below).
  return map_ptr;
}

/**
 * @brief print info regarding a python array. Debug tool.
 *
 * @param arr the input array
 *
 * @note No memory allocation!
 */
void siconos::pybind11_utils::inspect_array(py::array arr) {
  py::buffer_info info = arr.request();

  std::cout << "ndim    = " << info.ndim << std::endl;
  std::cout << "shape   = [";
  for (size_t i = 0; i < info.shape.size(); ++i) {
    std::cout << info.shape[i];
    if (i < info.shape.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "strides = [";
  for (size_t i = 0; i < info.strides.size(); ++i) {
    std::cout << info.strides[i];
    if (i < info.strides.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "itemsize = " << info.itemsize << " bytes" << std::endl;
  std::cout << "format   = " << info.format << std::endl;
  std::cout << "ptr      = " << info.ptr << std::endl;
  std::cout << "values = [";
  if (info.format == py::format_descriptor<double>::format()) {
    double* ptr = static_cast<double*>(info.ptr);
    for (size_t i = 0; i < info.shape[0]; ++i) {
      std::cout << ptr[i] << " ";
    }
  } else {  // if (info.format == py::format_descriptor<int64_t>::format()) {
    int64_t* ptr = static_cast<int64_t*>(info.ptr);
    for (size_t i = 0; i < info.shape[0]; ++i) {
      std::cout << ptr[i] << " ";
    }
  }
  // else {
  // std::cout << "<unsupported type>";
  // }
  std::cout << "]" << std::endl;
}

//    py::object siconos::pybind11_utils::eigensparse_to_scipy(
//     const siconos::algebra::SiconosSparseMatrix& M, py::handle owner) {
//   py::array_t<double> data({M.nonZeros()}, {sizeof(double)}, M.valuePtr(),
//                            py::capsule(M.valuePtr()));
//   data.attr("flags").attr("writeable") = false;  // readonly
//   py::array_t<siconos::algebra::SparseIndex> indices(
//       {M.nonZeros()}, {sizeof(siconos::algebra::SparseIndex)}, M.innerIndexPtr(),
//       py::capsule(M.innerIndexPtr()));
//   indices.attr("flags").attr("writeable") = false;  // readonly
//   py::array_t<siconos::algebra::SparseIndex> indptr(
//       {M.outerSize() + 1}, {sizeof(siconos::algebra::SparseIndex)}, M.outerIndexPtr(),
//       py::capsule(M.outerIndexPtr()));
//   indptr.attr("flags").attr("writeable") = false;  // readonly
//   py::module scipy_sparse = py::module::import("scipy.sparse");
//   py::tuple shape = py::make_tuple(M.rows(), M.cols());
//   py::object csc = scipy_sparse.attr("csc_matrix")(py::make_tuple(data, indices, indptr),
//                                                    shape, "copy"_a = false);
//   // No copy ... assuming parameters have the proper types (complient with scipy)

//   if (!owner.is_none()) {
//     csc.attr("_cpp_owner") = owner;
//   }
//   return csc;
// }

/**
 * @brief Creates a sparse python array (scipy) from a SiconosSparseMatrix (eigen)
 *
 * @param M the input matrix
 * @param owner the object which calls this function
 *
 * @note No memory allocation!
 */
py::object siconos::pybind11_utils::make_readonly_csc_array(
    const siconos::algebra::SiconosSparseMatrix& M, py::handle owner) {
  auto data = siconos::pybind11_utils::make_readonly_array(M.valuePtr(), M.nonZeros(), owner);

  auto indices =
      siconos::pybind11_utils::make_readonly_array(M.innerIndexPtr(), M.nonZeros(), owner);
  auto indptr = siconos::pybind11_utils::make_readonly_array(M.outerIndexPtr(),
                                                             M.outerSize() + 1, owner);

  py::object sp = py::module_::import("scipy.sparse");
  py::tuple shape = py::make_tuple(M.rows(), M.cols());

  return sp.attr("csc_array")(py::make_tuple(data, indices, indptr), shape);
}
