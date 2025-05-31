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

/*! \file SiconosMatrixOp.hpp
  \brief Toolbox for operators acting on matrices
  For matrix-vector op, see SiconosMatrixVectorOp.hpp
*/
#ifndef SICOMAT_OPH
#define SICOMAT_OPH
#include <memory>

#include "BlockMatrix.hpp"
#include "SiconosMatrix.hpp"

namespace siconos::algebra {

/** computes C = A*B or C += AB if init = false.
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param[in,out] C a SiconosMatrix
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init = true);

}  // namespace siconos::algebra

#endif
