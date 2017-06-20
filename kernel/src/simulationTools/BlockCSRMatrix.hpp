/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "SiconosNumerics.h"
#include "SimulationTypeDef.hpp"

/* with signed int typedef  boost::numeric::ublas::compressed_matrix<double*> CompressedRowMat; */
/* cf http://boost.2283326.n4.nabble.com/LU-decomposition-of-compressed-matrix-td3417929.html */
typedef ublas::compressed_matrix <
double*, ublas::basic_row_major<unsigned int>, 0, ublas::unbounded_array<std::size_t>
>  CompressedRowMat;


TYPEDEF_SPTR(CompressedRowMat)
TYPEDEF_SPTR(SparseBlockStructuredMatrix)

/** Definition of a compressed sparse row matrix of SiconosMatrix,
 * used in OneStepNSProblem to store the M matrix.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) 29/11/2007
 *
 * This class defines a specific compressed row sparse storage for
 * blocks matrices, each block being a SiconosMatrix*.
 *
 * It handles:
 * - a SparseMat (boost-ublas) of SiconosMatrix*
 * - a vector<SiconosMatrix*> which handles the non-null blocks

 * - three vector<int> (IndexInt) to save non-null blocks position in
     row, columns and the list of the sizes of diagonal blocks.

 * - two int, the number of blocks in a row and the number of non null blocks.
 *
 * Each block of the current object represents the connection between
 * two coupled Interactions, \n (for example for Lagrangian
 * systems, a single \f$ H W^{-1} H^t \f$ block or for first order
 * systems \f$ hCW^{-1}B \f$ ...) \n
 *
 * This objects is built using an index set of SP::Interaction,
 * that represents the "active" constraints in the OSNS problem and a
 * map<SP::Interaction u1, <SP::Interaction u2, SP::SiconosMatrix block> >, block being the link
 * between u1 and u2. Only Interaction present in the index set are picked out
 * in the map.
 *
 *  A convert method is also implemented to create a
 *  SparseBlockStructuredMatrix which is Numerics-readable.
 *
 * As an example, consider the index set I={u1, u3, u5, u8} and the
 * map where non null blocks are (ui,ui), (u1,u3), (u1,u8), (u3,u1),
 * (u8,u1).\n Each block being a pointer to a 3x3 matrix.\n Then the
 * resulting matrix has 4 X 4 blocks, with 8 non-null blocks and looks
 * like:
 *
 * \f{equation*}
   M=\left\lbrace\begin{array}{cccc}
    b11 & b13 & 0 & b18 \\
    b31 & b22 & 0 & 0 \\
    0   & 0   & b33&0 \\
    b81 & 0   & 0 & b44
    \end{array}\right.
    \f}
 *
 * with nc = 4, nbNonNullBlocks = 8, RowPos = [0 0 0 1 1 2 3 3],
 * RowCol = [0 1 3 0 1 2 0 3]\n and _diagsize0 = [3 6 9 12].
 *
 * We use stl::vector (which may seems redundent with the double* of
 * the numerics SparseBlockStructuredMatrix) because memory can be
 * reserved during construction or initialized and then vectors are
 * resized when the object is filled in. This avoid some call to
 * malloc/free at each iteration.
 *
 */
class BlockCSRMatrix
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BlockCSRMatrix);

  /** Number of blocks rows (first dimension of the block matrix)*/
  unsigned int _nr;

  /** Number of blocks columns (second dimension of the block matrix)*/
  unsigned int _nc;

  /** Sparse-Block Boost Matrix. Each block is a SiconosMatrix**/
  SP::CompressedRowMat _blockCSR;

  /** Specific structure required when a (Numerics) solver block is used */
  SP::SparseBlockStructuredMatrix _sparseBlockStructuredMatrix;

  /** Vector used to save the sum of rows of diagonal blocks of M:
      _diagsize0[i] = _diagsize0[i-1] + ni, ni being the size of the
      diagonal block at row(block) i */
  SP::IndexInt _diagsize0;

  /** Vector used to save the sum of dim of diagonal blocks of M:
      _diagsize0[i] = _diagsize0[i-1] + ni, ni being the size of the
      diagonal block at row(block) i */
  SP::IndexInt _diagsize1;

  /** List of non null blocks positions (in row) */
  SP::IndexInt rowPos;

  /** List of non null blocks positions (in col) */
  SP::IndexInt colPos;

  /** Private copy constructor => no copy nor pass by value */
  BlockCSRMatrix(const BlockCSRMatrix&);

  /** Private assignment -> forbidden 
   * \return  BlockCSRMatrix&
   */
  BlockCSRMatrix& operator=(const BlockCSRMatrix&);

public:

  /** Default constructor -> empty matrix
   */
  BlockCSRMatrix();

  /** Constructor with dimension (number of blocks)
      \param n number of blocks in a row/column (only square matrices allowed)
  */
  BlockCSRMatrix(unsigned int n);

  /** Constructor from index set
      \param indexSet the index set of the active constraints
  */
  BlockCSRMatrix(SP::InteractionsGraph indexSet);

  /** destructor
   */
  ~BlockCSRMatrix();

  /** get size (in block-components) 
   * \return unsigned int NumberOfBlocksInARow
   */
  inline unsigned int numberOfBlocksInARow() const
  {
    return _nr;
  };

  /** get total number of non-null blocks
   * \return unsigned int
   */
  unsigned int getNbNonNullBlocks() const;

  /** get the numerics-readable structure
   * \return SP::SparseBlockStructuredMatrix
   */
  inline SP::SparseBlockStructuredMatrix getNumericsMatSparse()
  {
    return _sparseBlockStructuredMatrix;
  };

  /** get the ublas sparse mat
   * \return SP::CompressedRowMat
   */
  inline SP::CompressedRowMat getMSparse()
  {
    return _blockCSR;
  };

  /** get the dimension of the square-diagonal block number num
   * \param i block position
   * \return unsigned int
  */
  unsigned int getSizeOfDiagonalBlock(int i) const
  {
    if (i == 0) return _diagsize0->at(0);
    else return (_diagsize0->at(i) - _diagsize0->at(i - 1));
  };

  /** get the index of blocks position (i=0 -> rows, i=1 -> columns)
   * \param i unsigned int, 0 for rows, 1 for columns
   * \return SP::IndexInt
   */
  inline SP::IndexInt getPositionsIndex(bool i)
  {
    if (i) return rowPos;
    else return colPos;
  };

  /** fill the current class using an index set
   *  \param indexSet set of the active constraints
   */
  void fill(SP::InteractionsGraph indexSet);


  /** fill the matrix with the Mass matrix 
   * \warning only for NewtonEulerDS
   * \param indexSet of the active constraints
   */
  void fillM(SP::InteractionsGraph indexSet);

  /** fill the matrix with the H matrix 
   * \warning only for NewtonEulerFrom3DLocalFrameR
   * \param indexSet of the active constraints
   */
  void fillH(SP::InteractionsGraph indexSet);

  /** fill the numerics structure _sparseBlockStructuredMatrix using _blockCSR */
  void convert();

  /** display the current matrix
   */
  void display() const;
};

#endif
