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
/*! \file BlockCSRMatrix.hpp
Definition of a compressed row sparse block matrix of SiconosMatrix*
*/

#ifndef BLOCKCSRMATRIX_H
#define BLOCKCSRMATRIX_H

#include <boost/numeric/ublas/fwd.hpp>  // Boost forward declarations
#include <vector>

#include "NumericsFwd.h"             // for SparseBlockStructuredMatrix
#include "SiconosSerialization.hpp"  // for ACCEPT_SERIALIZATION
#include "SimulationGraphs.hpp"

namespace siconos::simulation {

/* with signed int typedef  boost::numeric::ublas::compressed_matrix<double*>
 * CompressedRowMat; */
/* cf
 * http://boost.2283326.n4.nabble.com/LU-decomposition-of-compressed-matrix-td3417929.html
 */

/**
   Definition of a compressed sparse row matrix of SiconosMatrix,
   used in siconos::nonsmooth_formulations::OneStepNSProblem to store the M matrix.

   This class defines a specific compressed row sparse storage for
   blocks matrices, each block being a SiconosMatrix*.

   It handles:

   - a SparseMat (boost-ublas) of SiconosMatrix*
   - a vector<SiconosMatrix*> which handles the non-null blocks

   - three vector<int> (std::vector<unsigned int>) to save non-null blocks position in
   row, columns and the list of the sizes of diagonal blocks.

   - two int, the number of blocks in a row and the number of non null blocks.

   Each block of the current object represents the connection between
   two coupled Interactions,

   (for example for Lagrangian
   systems, a single \f$ H W^{-1} H^t \f$ block or for first order
   systems \f$ hCW^{-1}B \f$ ...).

   This objects is built using an index set of std::shared_ptr<siconos::modeling::Interaction>,
   that represents the "active" constraints in the OSNS problem and a
   map<std::shared_ptr<siconos::modeling::Interaction> u1,
   <std::shared_ptr<siconos::modeling::Interaction> u2,
   std::shared_ptr<siconos::algebra::SiconosMatrix> block> >, block being the link between u1
   and u2. Only Interaction present in the index set are picked out in the map.

   A convert method is also implemented to create a
   SparseBlockStructuredMatrix which is Numerics-readable.

   As an example, consider the index set I={u1, u3, u5, u8} and the
   map where non null blocks are (ui,ui), (u1,u3), (u1,u8), (u3,u1),
   (u8,u1). Each block being a pointer to a 3x3 matrix.
   Then the resulting matrix has 4 X 4 blocks, with 8 non-null blocks and looks
   like:

   \f[

   M=\left\lbrace\begin{array}{cccc}
   b11 & b13 & 0 & b18 \\
   b31 & b22 & 0 & 0   \\
   0   & 0   & b33&0 \\
   b81 & 0   & 0 & b44
   \end{array}\right.
   \f]

   with nc = 4, nbNonNullBlocks = 8, RowPos = [0 0 0 1 1 2 3 3],
   RowCol = [0 1 3 0 1 2 0 3] and _diagsize0 = [3 6 9 12].

   We use std::vector (which may seems redundent with the double* of
   the numerics SparseBlockStructuredMatrix) because memory can be
   reserved during construction or initialized and then vectors are
   resized when the object is filled in. This avoid some call to
   malloc/free at each iteration.

*/
class BlockCSRMatrix {
 private:
  using CompressedRowMat = boost::numeric::ublas::compressed_matrix<
      double *, boost::numeric::ublas::basic_row_major<unsigned int>, 0,
      boost::numeric::ublas::unbounded_array<std::size_t>>;

  ACCEPT_SERIALIZATION(BlockCSRMatrix);

  /** Number of blocks rows (first dimension of the block matrix)*/
  siconos::graphs::InteractionsGraph::size_type _nr{0};

  /** Number of blocks columns (second dimension of the block matrix)*/
  siconos::graphs::InteractionsGraph::size_type _nc{0};

  /** Sparse-Block Boost Matrix. Each block is a SiconosMatrix**/
  std::shared_ptr<CompressedRowMat> _blockCSR{nullptr};

  /** Specific structure required when a (Numerics) solver block is used */
  std::shared_ptr<SparseBlockStructuredMatrix> _sparseBlockStructuredMatrix{nullptr};

  /** Vector used to save the sum of rows of diagonal blocks of M:
      _diagsize0[i] = _diagsize0[i-1] + ni, ni being the size of the
      diagonal block at row(block) i */
  std::shared_ptr<std::vector<unsigned int>> _diagsize0{nullptr};

  /** Vector used to save the sum of dim of diagonal blocks of M:
      _diagsize0[i] = _diagsize0[i-1] + ni, ni being the size of the
      diagonal block at row(block) i */
  std::shared_ptr<std::vector<unsigned int>> _diagsize1{nullptr};

  /** List of non null blocks positions (in row) */
  std::shared_ptr<std::vector<unsigned int>> rowPos{nullptr};

  /** List of non null blocks positions (in col) */
  std::shared_ptr<std::vector<unsigned int>> colPos{nullptr};

  // Rule of five
  BlockCSRMatrix() = delete;
  BlockCSRMatrix(const BlockCSRMatrix &) = delete;
  BlockCSRMatrix(BlockCSRMatrix &&) = delete;
  BlockCSRMatrix &operator=(const BlockCSRMatrix &) = delete;
  BlockCSRMatrix &operator=(BlockCSRMatrix &&) = delete;

 public:
  /** Constructor with dimension (number of blocks)
   *
   *  \param n number of blocks in a row/column (only square matrices allowed)
   */
  BlockCSRMatrix(siconos::graphs::InteractionsGraph::size_type n);

  /** Constructor from index set
   *
   *  \param indexSet the index set of the active constraints
   */
  BlockCSRMatrix(siconos::graphs::InteractionsGraph &indexSet);

  /** destructor
   */
  ~BlockCSRMatrix() noexcept = default;

  /** get size (in block-components)
   *
   *  \return unsigned int NumberOfBlocksInARow
   */
  inline auto numberOfBlocksInARow() const { return _nr; };

  /** get total number of non-null blocks
   *
   *  \return unsigned int
   */
  std::size_t getNbNonNullBlocks() const;

  /** \return the numerics-readable structure
   */
  inline std::shared_ptr<SparseBlockStructuredMatrix> getNumericsMatSparse() {
    return _sparseBlockStructuredMatrix;
  };

  /** get the ublas sparse mat
   *
   *  \return std::shared_ptr<CompressedRowMat>
   */
  inline std::shared_ptr<CompressedRowMat> getMSparse() { return _blockCSR; };

  /** get the dimension of the square-diagonal block number num
   *
   *  \param i block position
   *  \return unsigned int
   */
  std::vector<unsigned int>::value_type getSizeOfDiagonalBlock(int i) const {
    if (i == 0)
      return _diagsize0->at(0);
    else
      return (_diagsize0->at(i) - _diagsize0->at(i - 1));
  };

  /** get the index of blocks position (i=0 -> rows, i=1 -> columns)
   *
   *  \param i unsigned int, 0 for rows, 1 for columns
   *  \return std::shared_ptr<std::vector<unsigned int>>
   */
  inline std::shared_ptr<std::vector<unsigned int>> getPositionsIndex(bool i) {
    if (i)
      return rowPos;
    else
      return colPos;
  };

  /** fill the current class using an index set
   *
   *  \param indexSet set of the active constraints
   */
  void fill(siconos::graphs::InteractionsGraph &indexSet);

  /** fill the matrix with the Mass matrix
   *
   *  \warning only for NewtonEulerDS
   *
   *  \param indexSet of the active constraints
   */
  void fillW(siconos::graphs::InteractionsGraph &indexSet);

  /** fill the matrix with the H matrix
   *
   *  \warning only for NewtonEuler3DR
   *
   *  \param indexSet of the active constraints
   */
  void fillH(siconos::graphs::InteractionsGraph &indexSet);

  /** fill the numerics structure _sparseBlockStructuredMatrix using _blockCSR
   */
  void convert();

  /** display the current matrix
   */
  void display() const;
};
}  // namespace siconos::simulation
#endif
