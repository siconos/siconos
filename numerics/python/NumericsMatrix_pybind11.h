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
#ifndef NM_PB11_H
#define NM_PB11_H

/* NumericsMatrix_pybind11.h
    Some tools to deal with NumericsMatrix and double*
    has members of numerics structures when it comes
    to wrap things with pybind11
*/

// Example of use: see friction_contact_pybind11.cpp

// Works only with NM_DENSE.
// TODO: sparse matrices

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include "NumericsMatrix.h"

namespace py = pybind11;

static void inspect_array(py::array arr) {
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
}

template <typename T>
py::array_t<double> get_matrix(const T &self, NumericsMatrix *T::*matrix_attr) {
  NumericsMatrix *matrix = self.*matrix_attr;

  if (!matrix || !matrix->matrix0) {
    throw std::runtime_error("Input matrix has not been allocated.");
  }

  auto rows = matrix->size0;
  auto cols = matrix->size1;
  return py::array_t<double>({rows, cols},                             // Shape
                             {sizeof(double), sizeof(double) * rows},  // Strides
                             matrix->matrix0, py::capsule(matrix->matrix0));
}

template <typename T>
py::object get_matrix_sparse(const T &self, NumericsMatrix *T::*matrix_attr) {
  NumericsMatrix *matrix = self.*matrix_attr;
  if (!matrix || !matrix->matrix2 || !matrix->matrix2->csc) {
    throw std::runtime_error("input matrix is not allocated.");
  }

  CSparseMatrix *csc = matrix->matrix2->csc;
  // py::capsule data_capsule(csc->x, [](void *) {});
  // py::capsule indices_capsule(csc->i, [](void *) {});
  // py::capsule indptr_capsule(csc->p, [](void *) {});

  py::array_t<double> data({csc->nzmax}, {sizeof(double)}, csc->x, py::capsule(csc->x));
  // Note FP: capsule only for csc->x.
  // .cast for csc->i, csc->p
  // No copies.
  // This is the only way I found to ensure a correct behavior.
  // All other ways lead to strange results, with data corruption
  // after a sequence of print(M) in python ...
  //                                    .cast<py::array_t<int64_t>>();

  // py::array_t<int64_t> indptr = py::array_t<int64_t>({csc->n + 1}, {sizeof(int64_t)},
  // csc->p)
  //                                   .cast<py::array_t<int64_t>>();
  py::array_t<int64_t> indices({csc->nzmax}, {sizeof(int64_t)}, csc->i, py::capsule(csc->i));
  py::array_t<int64_t> indptr({csc->n + 1}, {sizeof(int64_t)}, csc->p, py::capsule(csc->p));

  // py::array_t<int64_t> indices = py::array_t<int64_t>({csc->nzmax}, {sizeof(int64_t)},
  // csc->i)
  //                                    .cast<py::array_t<int64_t>>();

  // py::array_t<int64_t> indptr = py::array_t<int64_t>({csc->n + 1}, {sizeof(int64_t)},
  // csc->p)
  //                                   .cast<py::array_t<int64_t>>();

  inspect_array(indices);
  inspect_array(data);
  // Build Python (scipy) csc (no copies)
  py::module scipy_sparse = py::module::import("scipy.sparse");
  py::object csc_matrix = scipy_sparse.attr("csc_matrix")(
      py::make_tuple(data, indices, indptr), py::make_tuple(csc->m, csc->n));
  // csc_matrix.attr("_keep_data") = data;
  // csc_matrix.attr("_keep_indices") = indices;
  // csc_matrix.attr("_keep_indptr") = indptr;

  return csc_matrix;
}
template <typename T>
void set_matrix(T &self, NumericsMatrix *T::*matrix_attr, py::array_t<double> array) {
  int size = self.dimension * self.numberOfContacts;

  if (array.ndim() != 2) {
    throw std::runtime_error("Input matrix must be a 2D matrix.");
  }
  NumericsMatrix *&matrix = self.*matrix_attr;
  if (!matrix) {
    matrix = new NumericsMatrix();
    matrix->storageType = NM_DENSE;
  }

  matrix->size0 = array.shape(0);
  matrix->size1 = array.shape(1);
  matrix->matrix0 = static_cast<double *>(array.mutable_data());
}

template <typename T>
void set_array(T &self, double *&array_ptr, py::array_t<double> arr, int expected_size = -1) {
  // expected size must be given when it's possible. If equal to -1, then no check
  auto buf = arr.request();
  if (expected_size > 0) {
    if (buf.ndim != 1 || buf.shape[0] != expected_size)
      throw std::runtime_error("Incorrect array size");
  }
  array_ptr = static_cast<double *>(buf.ptr);
}

template <typename T>
void set_matrix_sparse(T &self, NumericsMatrix *T::*matrix_attr, py::object pymat) {
  NumericsMatrix *&matrix = self.*matrix_attr;
  if (!matrix) {
    matrix = new NumericsMatrix();
    matrix->storageType = NM_SPARSE;
    matrix->matrix2 = new NumericsSparseMatrix();
  }

  NumericsSparseMatrix *sparseMat = matrix->matrix2;

  // Ensure CSC format (conversion). Copy if not? Check this
  // py::object csc_M = pymat.attr("tocsc")();
  py::object csc_M = pymat.attr("tocsc")();

  // Get numpy array. Warn FP: mind types (must be int64 in scipy because of cs.h types)
  py::array_t<double> data_array = pymat.attr("data").cast<py::array_t<double>>();
  py::array_t<int64_t> indptr_array = pymat.attr("indptr").cast<py::array_t<int64_t>>();
  py::array_t<int64_t> indices_array = pymat.attr("indices").cast<py::array_t<int64_t>>();
  auto shape = pymat.attr("shape").cast<std::pair<int, int>>();

  if (indices_array.itemsize() != sizeof(int64_t)) {
    throw std::runtime_error(
        "indices must be of type int64_t to be complient with SuiteSparse");
  }

  // Build CSparseMatrix
  sparseMat->csc = new CSparseMatrix();
  sparseMat->csc->nzmax = data_array.shape(0);
  sparseMat->csc->m = shape.first;
  sparseMat->csc->n = shape.second;

  inspect_array(indices_array);
  inspect_array(data_array);
  //  sparseMat->csc->p = static_cast<int64_t *>(indptr_array.mutable_data());
  //  sparseMat->csc->i = static_cast<int64_t *>(indices_array.mutable_data());
  // sparseMat->csc->p = static_cast<int64_t *>(indptr_array.request().ptr);
  // sparseMat->csc->i = static_cast<int64_t *>(indices_array.request().ptr);
  //  sparseMat->csc->x = static_cast<double *>(data_array.mutable_data());
  // sparseMat->csc->x = static_cast<double *>(data_array.request().ptr);
  sparseMat->csc->x =
      csc_M.attr("data")
          .cast<py::array_t<double, py::array::c_style | py::array::forcecast>>()
          .mutable_data();

  sparseMat->csc->i =
      csc_M.attr("indices")
          .cast<py::array_t<int64_t, py::array::c_style | py::array::forcecast>>()
          .mutable_data();

  sparseMat->csc->p =
      csc_M.attr("indptr")
          .cast<py::array_t<int64_t, py::array::c_style | py::array::forcecast>>()
          .mutable_data();

  sparseMat->csc->nz = -1;  // CSC format

  matrix->size0 = shape.first;
  matrix->size1 = shape.second;

  // Note FP: capsules and M.attr version below.
  // stuff finally not needed. Check this carefully
  //   // Capsules to ensure that data/memories are kept alive
  //   py::capsule data_capsule = py::reinterpret_steal<py::capsule>(
  //       PyCapsule_New(static_cast<void *>(data_array.mutable_data()), nullptr, nullptr));

  //   py::capsule indices_capsule = py::reinterpret_steal<py::capsule>(
  //       PyCapsule_New(static_cast<void *>(indices_array.mutable_data()), nullptr, nullptr));

  //   py::capsule indptr_capsule = py::reinterpret_steal<py::capsule>(
  //       PyCapsule_New(static_cast<void *>(indptr_array.mutable_data()), nullptr, nullptr));

  //   csc_M.attr("data").attr("setflags")(py::arg("write") = false);
  //   csc_M.attr("indices").attr("setflags")(py::arg("write") = false);
  //   csc_M.attr("indptr").attr("setflags")(py::arg("write") = false);
}

// Utility struct to keep buffers from numpy 'alive'
struct DenseMatrixMemoryOwner {
  py::array_t<double, py::array::f_style> data;
};

// Utility struct to keep buffers from scipy sparse 'alive'
struct SparseMatrixMemoryOwner {
  py::object scipy_ref;
  py::array data, indices, indptr;
};

/** Wrap the transfer from a numpy matrix to a XXpb matrix (e.g. M for FrictionContactProblem)

    \param self the struct representing the problem
    \param scipy_matrix python matrix
    \param mem a struct to keep memory alive in the wrapper
    \return the numerics matrix with some 'views' onto numpy memory (set matrix0 of
    NumericsMatrix)
 */
template <typename T>
NumericsMatrix *set_matrix_dense_new(T &self, py::object numpy_matrix,
                                     std::shared_ptr<DenseMatrixMemoryOwner> &mem) {
  if (!py::isinstance<py::array>(numpy_matrix))
    throw std::runtime_error("Expected a Numpy array");

  auto array = numpy_matrix.cast<py::array_t<double, py::array::f_style>>();

  if (array.ndim() != 2) throw std::runtime_error("Matrix must be 2D");

  mem = std::make_shared<DenseMatrixMemoryOwner>();
  mem->data = array;

  auto *nm = new NumericsMatrix;
  nm->size0 = array.shape(0);
  nm->size1 = array.shape(1);
  nm->matrix0 = static_cast<double *>(mem->data.mutable_data());
  nm->matrix2 = nullptr;
  nm->storageType = NM_DENSE;

  return nm;
}

/** Wrap the transfer from a scipy sparse to a XXpb matrix (e.g. M for FrictionContactProblem)

    \param self the struct representing the problem
    \param scipy_matrix python sparse matrix
    \param mem a struct to keep memory alive in the wrapper
    \return the numerics matrix with some 'views' onto scipy sparse memory (set matrix2 of
   NumericsMatrix)
 */
template <typename T>
NumericsMatrix *set_matrix_sparse_new(T &self, py::object scipy_matrix,
                                      std::shared_ptr<SparseMatrixMemoryOwner> &mem) {
  // Check input matrix validity
  if (!py::hasattr(scipy_matrix, "data") || !py::hasattr(scipy_matrix, "indices") ||
      !py::hasattr(scipy_matrix, "indptr"))
    throw std::runtime_error("Not a valid CSC scipy sparse matrix");

  // Link/views between scipy sparse input matrix and the 'memory owner'
  mem = std::make_shared<SparseMatrixMemoryOwner>();
  mem->scipy_ref = scipy_matrix;  // To avoid memory destruction
  mem->data = scipy_matrix.attr("data").cast<py::array_t<double, py::array::c_style>>();
  mem->indices = scipy_matrix.attr("indices").cast<py::array_t<int64_t, py::array::c_style>>();
  mem->indptr = scipy_matrix.attr("indptr").cast<py::array_t<int64_t, py::array::c_style>>();

  auto shape = scipy_matrix.attr("shape").cast<std::pair<ssize_t, ssize_t>>();

  // Build CSparse matrix with mem views
  auto *csc = new CSparseMatrix;
  csc->m = shape.first;
  csc->n = shape.second;
  csc->x = static_cast<double *>(mem->data.mutable_data());
  csc->i = static_cast<int64_t *>(mem->indices.mutable_data());
  csc->p = static_cast<int64_t *>(mem->indptr.mutable_data());
  csc->nzmax = mem->data.shape(0);
  csc->nz = -1;

  // Build NumericsMatrix with CSparse
  auto *nsm = new NumericsSparseMatrix;
  nsm->csc = csc;

  auto *nm = new NumericsMatrix;
  nm->size0 = csc->m;
  nm->size1 = csc->n;
  nm->matrix2 = nsm;
  nm->storageType = NM_SPARSE;

  return nm;
}

template <typename T>
py::array_t<double> get_array(const T &self, double *array_ptr, int size) {
  return py::array_t<double>({size},                            // Shape
                             {sizeof(double)},                  // Strides
                             array_ptr, py::capsule(array_ptr)  // Data pointer
  );
}

#endif