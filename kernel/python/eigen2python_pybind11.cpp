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

// template <typename T>
py::object siconos::pybind11_utils::eigensparse_to_scipy(
    const siconos::algebra::SiconosSparseMatrix& M, py::handle owner) {
  py::array_t<double> data({M.nonZeros()}, {sizeof(double)}, M.valuePtr(),
                           py::capsule(M.valuePtr()));
  data.attr("flags").attr("writeable") = false;  // readonly
  py::array_t<siconos::algebra::SparseIndex> indices(
      {M.nonZeros()}, {sizeof(siconos::algebra::SparseIndex)}, M.innerIndexPtr(),
      py::capsule(M.innerIndexPtr()));
  indices.attr("flags").attr("writeable") = false;  // readonly
  py::array_t<siconos::algebra::SparseIndex> indptr(
      {M.outerSize() + 1}, {sizeof(siconos::algebra::SparseIndex)}, M.outerIndexPtr(),
      py::capsule(M.outerIndexPtr()));
  indptr.attr("flags").attr("writeable") = false;  // readonly
  py::module scipy_sparse = py::module::import("scipy.sparse");
  py::tuple shape = py::make_tuple(M.rows(), M.cols());
  py::object csc = scipy_sparse.attr("csc_matrix")(py::make_tuple(data, indices, indptr),
                                                   shape, "copy"_a = false);
  // No copy ... assuming parameters have the proper types (complient with scipy)

  if (!owner.is_none()) {
    csc.attr("_cpp_owner") = owner;
  }
  return csc;
}

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

// inline py::object make_readonly_csc_array_capsule(
//     const siconos::algebra::SiconosSparseMatrix& M) {  //, py::handle owner) {
//   // créer les arrays numpy readonly

//   auto capsule_data = py::capsule(M.valuePtr(), [](void* /*p*/) {});

//   py::array arr_data(py::buffer_info(const_cast<double*>(M.valuePtr()), sizeof(double),
//                                      py::format_descriptor<double>::format(), 1,
//                                      {M.nonZeros()}, {sizeof(double)}),
//                      capsule_data);

//   arr_data.attr("flags").attr("writeable") = false;
//   std::cout << "DATA OK \n";
//   auto capsule_indices = py::capsule(M.innerIndexPtr(), [](void* /*p*/) {});

//   //   py::array_t<int64_t> arr_indices(
//   //       py::buffer_info(const_cast<int*>(M.innerIndexPtr()), sizeof(int64_t),
//   //                       py::format_descriptor<int64_t>::format(), 1, {M.nonZeros()},
//   //                       {sizeof(int64_t)}),
//   //       capsule_indices);
//   //   py::array_t<int64_t> arr_indices({M.nonZeros()}, {sizeof(int64_t)},
//   //                                    const_cast<int*>(M.innerIndexPtr()),
//   //                                    py::capsule(M.innerIndexPtr()));
//   //  py::buffer_info(const_cast<int*>(M.innerIndexPtr()), sizeof(int64_t),
//   //                  py::format_descriptor<int64_t>::format(), 1, {M.nonZeros()},
//   //                  {sizeof(int64_t)}),
//   //  capsule_indices);

//   //  arr_indices.attr("flags").attr("writeable") = false;
//   std::cout << "INDICES OK \n";

//   auto capsule_indptr = py::capsule(M.outerIndexPtr(), [](void* /*p*/) {});

//   py::array_t<int64_t> arr_indptr(
//       py::buffer_info(const_cast<int*>(M.outerIndexPtr()), sizeof(int64_t),
//                       py::format_descriptor<int64_t>::format(), 1, {M.outerSize() + 1},
//                       {sizeof(int64_t)}),
//       capsule_indptr);

//   arr_indptr.attr("flags").attr("writeable") = false;
//   std::cout << "INDPTR OK \n";

//   std::cout << "sizeof(valuePtr) = " << sizeof(*M.valuePtr()) << "\n";
//   std::cout << "sizeof(innerIndexPtr) = " << sizeof(*M.innerIndexPtr()) << "\n";
//   std::cout << "sizeof(outerIndexPtr) = " << sizeof(*M.outerIndexPtr()) << "\n";

//   // importer scipy.sparse
//   py::object sp = py::module_::import("scipy.sparse");
//   py::tuple shape = py::make_tuple(M.rows(), M.cols());

//   return sp.attr("csc_array")(py::make_tuple(arr_data, arr_indptr, arr_indptr), shape);
// }

// inline py::object make_readonly_csc_array_map(
//     const Eigen::Map<siconos::algebra::SiconosSparseMatrix>& M, py::handle owner) {
//   auto data =
//       py::array(py::buffer_info(const_cast<double*>(M.valuePtr()), sizeof(double),
//                                 py::format_descriptor<double>::format(), 1,
//                                 {static_cast<ssize_t>(M.nonZeros())}, {sizeof(double)}),
//                 owner);
//   data.attr("flags").attr("writeable") = false;

//   auto indices = py::array(py::buffer_info(const_cast<int*>(M.innerIndexPtr()),
//   sizeof(int),
//                                            py::format_descriptor<int>::format(), 1,
//                                            std::vector<ssize_t>{sizeof(int)}),
//                            owner);
//   indices.attr("flags").attr("writeable") = false;

//   auto indptr = py::array(py::buffer_info(const_cast<int*>(M.outerIndexPtr()),
//   sizeof(int),
//                                           py::format_descriptor<int>::format(), 1,
//                                           std::vector<ssize_t>{sizeof(int)}),
//                           owner);
//   indptr.attr("flags").attr("writeable") = false;

//   py::object sp = py::module_::import("scipy.sparse");
//   py::tuple shape = py::make_tuple(M.rows(), M.cols());

//   return sp.attr("csc_array")(py::make_tuple(data, indices, indptr), shape);
// }
