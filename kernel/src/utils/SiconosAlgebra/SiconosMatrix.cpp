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

#include "SiconosMatrix.hpp"

#include "BlockMatrix.hpp"  // for BlockMatrix
// #include "NumericsToolsNamespace.h"  // for CSparseMatrix,  NSM_fix_csc
// #include "SiconosMatrixOp.hpp"       // for matrix operators declaration
#include "SiconosVector.hpp"  // for SiconosVector
#include "Tools.hpp"          // toString

bool siconos::algebra::fillTriplet(SiconosMatrix& m, CSparseMatrix* triplet, size_t row_off,
                                   size_t col_off, double tol) {
  assert(triplet);
  size_t nrow = m.rows();
  size_t ncol = m.cols();

  double* arr = m.data();
  for (size_t j = 0; j < ncol; ++j) {
    for (size_t i = 0; i < nrow; ++i) {
      // col-major

      CSparseMatrix_zentry(triplet, i + row_off, j + col_off, arr[i + j * nrow],
                           std::numeric_limits<double>::epsilon());
    }
  }
  return true;
}