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

#ifndef EIGENPYTHON_PB11
#define EIGENPYTHON_PB11

/*! \file eigen2python_pybind11.hpp
  Helpers (class and free functions) to deal with Eigen sparse or dense matrices and vector
  wrap to python
*/
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SiconosMatrix.hpp"

namespace siconos::pybind11_utils {

/**
 * @brief Convert a SciPy CSC (Compressed Sparse Column) matrix into an Eigen sparse matrix.
 *
 * This function reconstructs an `Eigen::SparseMatrix<double>` from a Python
 * `scipy.sparse.csc_matrix` or `scipy.sparse.csc_array` object. The conversion
 * is performed by directly reading the internal arrays `data`, `indices`,
 * and `indptr` of the SciPy object, and inserting the corresponding triplets
 * into a newly allocated Eigen sparse matrix.
 *
 * @param[in] csc
 *   A Python object representing a SciPy CSC sparse matrix or array
 *   (i.e. `scipy.sparse.csc_matrix` or `scipy.sparse.csc_array`).
 *
 * @return
 *   A newly constructed `siconos::algebra::SiconosSparseMatrix`
 *   (alias of `Eigen::SparseMatrix<double>`) containing the same data
 *   as the input CSC object.
 *
 * @details
 *   - The conversion is **deep**: data is copied from Python into a new C++ matrix.
 *   - The matrix is built column by column, using the CSC index structure.
 *   - The result is compressed before returning (`makeCompressed()`).
 *
 * @note
 *   This function assumes the CSC object stores data in standard
 *   SciPy-compatible arrays (`data`, `indices`, and `indptr`) and
 *   uses 0-based indexing.
 *
 * @warning
 *   - No type checking is performed beyond the array casting.
 *   - Modifying the returned Eigen matrix does *not* affect the original Python object.
 */
siconos::algebra::SiconosSparseMatrix csc_to_eigen(pybind11::object csc);

/**
 * @brief Create an Eigen::Map view of a SciPy CSC sparse matrix without copying data.
 *
 * This function creates a shared pointer to an `Eigen::Map` object referencing
 * the internal NumPy arrays of a SciPy CSC sparse matrix (`csc_matrix` or `csc_array`).
 * The resulting Eigen map allows read/write access directly into the Python-owned
 * memory, with **zero copy** overhead.
 *
 * @param[in] csc
 *   A Python object representing a SciPy CSC sparse matrix or array
 *   (i.e. `scipy.sparse.csc_matrix` or `scipy.sparse.csc_array`).
 *
 * @return
 *   A `std::shared_ptr<Eigen::Map<SiconosSparseMatrix>>` referencing
 *   the same data as the input CSC object.
 *
 * @details
 *   - The function extracts the arrays `data`, `indices`, and `indptr`
 *     using pybind11, then builds an `Eigen::Map` over them.
 *   - The returned `shared_ptr` keeps the map alive, but the actual
 *     memory remains owned by Python — it is crucial that the Python
 *     `csc` object remains alive while the map is in use.
 *   - This allows **bidirectional interoperability**: modifications
 *     in C++ are visible in Python, and vice versa.
 *
 * @warning
 *   - No memory copy occurs; lifetime safety depends on Python’s garbage collector.
 *   - If the Python CSC object is destroyed, the Eigen::Map becomes invalid.
 *   - Ensure that the Python object is attached as an owner (e.g. `_cpp_owner`)
 *     when exposing the mapped matrix via pybind11.
 */
std::shared_ptr<Eigen::Map<siconos::algebra::SiconosSparseMatrix>> csc_to_eigen_map(
    pybind11::object csc);

/**
 * @brief Print diagnostic information about a NumPy array.
 *
 *
 * @param[in] arr
 *   A `pybind11::array` object representing a NumPy array.
 *
 * @note
 *   - This function does not modify the input array.
 *   - Only the first dimension of the array is iterated for displaying values.
 *   - Printing large arrays may produce significant console output.
 *
 */
void inspect_array(pybind11::array arr);

pybind11::object eigensparse_to_scipy(const siconos::algebra::SiconosSparseMatrix& M,
                                      pybind11::handle owner = pybind11::none());

/**
 * @brief Create a non-writeable (read-only) NumPy array view on an existing C++ buffer.
 *
 * This function wraps a raw C++ pointer as a `pybind11::array` without copying the data.
 * The array is marked as read-only from the Python side, ensuring that Python code
 * cannot modify the underlying C++ memory.
 *
 * The lifetime of the buffer is tied to the provided @p owner Python object.
 * As long as the @p owner exists in Python, the NumPy array remains valid.
 *
 * @tparam T
 *   The element type of the array (e.g. `double`, `int`, etc.).
 *
 * @param[in] ptr
 *   Pointer to the first element of the C++ buffer to expose.
 *
 * @param[in] size
 *   Number of elements in the buffer (1D array length).
 *
 * @param[in] owner
 *   A Python handle (e.g. `py::cast(self)`) whose lifetime controls
 *   the validity of the underlying C++ buffer. When the owner is destroyed,
 *   the array becomes invalid.
 *
 * @return
 *   A `pybind11::array` object referencing the existing C++ data,
 *   marked as non-writeable (`flags.writeable = False`).
 *
 * @note
 *   - No copy of the data is performed.
 *   - The returned array must not outlive the @p owner object.
 */
template <typename T>
pybind11::array make_readonly_array(const T* ptr, ssize_t size, pybind11::handle owner) {
  pybind11::array arr(
      pybind11::buffer_info(const_cast<T*>(ptr),  // Remove const
                            sizeof(T), pybind11::format_descriptor<T>::format(),
                            1,           // 1D array
                            {size},      // shape
                            {sizeof(T)}  // strides
                            ),
      owner);
  arr.attr("flags").attr("writeable") = false;  // readonly
  // arr.setflags(false);
  // arr.attr("setflags")(false);
  return arr;
}

/**
 * @brief Construct a SciPy CSC sparse array from an Eigen::SparseMatrix without copying data.
 *
 * This function builds a read-only `scipy.sparse.csc_array` (or `csc_matrix`)
 * using the internal compressed column storage (CCS/CSC) data of an
 * `Eigen::SparseMatrix`. The resulting Python object directly views
 * the C++ memory, ensuring zero-copy interoperation.
 *
 * @param[in] M
 *   The Eigen sparse matrix (`Eigen::SparseMatrix<double>`) whose internal data
 *   pointers will be exposed.
 *
 * @param[in] owner
 *   A Python object controlling the lifetime of @p M.
 *
 * @return
 *   A `scipy.sparse.csc_array` Python object referencing the C++ matrix data.
 *
 * @note
 *   - No data copy is performed.
 *   - The resulting arrays are read-only from Python.
 *   - The function assumes column-major (`Eigen::ColMajor`) storage.
 *   - The returned object becomes invalid if the C++ matrix is destroyed
 *     before the Python owner.
 *
 * @warning
 *   The returned object should only be used for read-only access.
 *   Any modification attempt from Python will raise an exception.
 *
 */
pybind11::object make_readonly_csc_array(const siconos::algebra::SiconosSparseMatrix& M,
                                         pybind11::handle owner);

template <typename T>
pybind11::array make_strict_readonly_array(const T* ptr, ssize_t size) {
  auto capsule = pybind11::capsule(ptr, [](void* /*p*/) {});

  pybind11::array arr(
      pybind11::buffer_info(const_cast<T*>(ptr), sizeof(T),
                            pybind11::format_descriptor<T>::format(), 1, {size}, {sizeof(T)}),
      capsule);

  arr.attr("flags").attr("writeable") = false;
  return arr;
}

pybind11::object make_readonly_csc_array_capsule(
    const siconos::algebra::SiconosSparseMatrix& M);

}  // namespace siconos::pybind11_utils
#endif