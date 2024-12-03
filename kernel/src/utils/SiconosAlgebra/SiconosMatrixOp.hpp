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

/** send data of the matrix to an ostream
 * \param os An output stream
 * \param bm a BlockMatrix
 * \return The same output stream
 */
// std::ostream& operator<<(std::ostream& os, const SiconosMatrix& sm);
// std::ostream& operator<<(std::ostream& os, const SiconosMatrix& sm);
std::ostream& operator<<(std::ostream& os, const BlockMatrix& sm);

// /** Compute the matrix exponential Exp = exp(A) for general matrices,
//   using scaling and Padé approximation. See expm.hpp.
//   \param A : input matrix
//   \param Exp : result = exp(A)
//   \param computeAndAdd : if true, result = result + exp(A)
// **/
// void expm(SiconosMatrix& A, SiconosMatrix& Exp, bool computeAndAdd = false);

/** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are
 *  given in dim and start.
 *  \param MIn a SPC::SiconosMatrix \param[in,out] MOut a std::shared_ptr<SiconosMatrix>
 *  \param dim an Index, dim[0], dim[1]: number of rows and columns of the sub-block
 *  \param start an Index, start[0], start[1]: position (row, column) of the first
 *  element of the sub-block in MIn start[2], start[3]: position (row, column) of the first
 *  element of the sub-block in MOut.
 */
void setBlock(const SiconosMatrix& MIn, std::shared_ptr<SiconosMatrix> MOut,
              const std::vector<std::size_t>& dim, const std::vector<std::size_t>& start);

// /** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are
//  *  given in dim and start.
//  *  \param MIn a SPC::SiconosMatrix \param[in,out] MOut a std::shared_ptr<SiconosMatrix>
//  *  \param dim an Index, dim[0], dim[1]: number of rows and columns of the sub-block
//  *  \param start an Index, start[0], start[1]: position (row, column) of the first
//  *  element of the sub-block in MIn start[2], start[3]: position (row, column) of the first
//  *  element of the sub-block in MOut.
//  */
// void setBlock(const SiconosMatrix& MIn, std::shared_ptr<MapType> MOut,
//               const std::vector<std::size_t>& dim, const std::vector<std::size_t>& start);

}  // namespace siconos::algebra

#endif
