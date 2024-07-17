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

/*! \file PlainObjectAddons.hpp

   Eigen PlainObjectBase class extensions.æ

   See https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html

*/

inline Index size(unsigned int index) const {
  if (index == 0) {
    return this->rows();
  } else if (index == 1) {
    return this->cols();
  } else {
    return 0;
  }
}

inline Index size() const { return this->rows() * this->cols(); }

void display() const { std::cout << *this << std::endl; }

inline void zero() { this->setZero(); }

inline Scalar norm2() { return this->norm(); }

// inline Scalar normInf() {
//     return this->cwiseAbs().rowwise().sum().maxCoeff();
// }

size_t nnz(double tol = 1e-14) {
  size_t nnz = 0;
  double* arr = this->data();
  for (size_t i = 0; i < this->size(0) * this->size(1); ++i) {
    if (fabs(arr[i]) > tol) {
      nnz++;
    }
  }
  return nnz;
}

inline Scalar vector_sum() {
  assert(this->cols() == 1);
  return this->sum();
}

/** copy the vector into an array
 *
 *  \param data the memory where to copy the data
 *  \return the number of element written (size of the vector)
 */
inline unsigned int copyData(double* data) {
  assert(this->cols() == 1);
  unsigned int size = this->size();
  if constexpr (std::is_same_v<Scalar, double>) {
    std::memcpy(data, this->data(),
                size * sizeof(Scalar));  // WARNING : Remember that
                                         // Eigen::Scalar has to be double
  } else {
    []<bool flag = false>() { static_assert(flag, "Eigen::Scalar has to be double"); }
    ();
  }
  return size;
}
