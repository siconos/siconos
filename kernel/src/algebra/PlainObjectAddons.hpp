// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file PlainObjectAddons.hpp

   Eigen PlainObjectBase class extensions.æ

   See https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html

*/

/** \return the number of non zero coefficients
 * \param tol tolerance value under which a coeff. is considered to be null
 */
size_t nonZeros(double tol = 1e-14) {
  // Warning: nonZeros to keep the same method name as in Eigen Sparse
  size_t nnz = 0;
  double* arr = this->data();
  for (size_t i = 0; i < this->rows() * this->cols(); ++i) {
    if (fabs(arr[i]) > tol) {
      nnz++;
    }
  }
  return nnz;
}
