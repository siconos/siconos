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

namespace py = pybind11;

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
py::array_t<double> get_array(const T &self, double *array_ptr, int size) {
  return py::array_t<double>({size},                            // Shape
                             {sizeof(double)},                  // Strides
                             array_ptr, py::capsule(array_ptr)  // Data pointer
  );
}

#endif