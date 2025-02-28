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

void display() const { std::cout << std::scientific << std::setprecision(6) << *this << "\n"; }

void displayT() const {
  std::cout << std::scientific << std::setprecision(6) << this->transpose() << "\n";
}

inline Scalar norm2() { return this->norm(); }

// inline Scalar normInf() {
//     return this->cwiseAbs().rowwise().sum().maxCoeff();
// }

size_t nnz(double tol = 1e-14) {
  size_t nnz = 0;
  double* arr = this->data();
  for (size_t i = 0; i < this->rows() * this->cols(); ++i) {
    if (fabs(arr[i]) > tol) {
      nnz++;
    }
  }
  return nnz;
}
