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

/*! \file BlockMatrix.hpp
  \brief Object to handle block-matrices.

*/

#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <memory>

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

/** "Block" matrix, ie container of matrices
 *
 * A BlockMatrix is a boost::ublas::compressed_matrix of std::shared_ptr<SiconosMatrix>.
 *
 * The blocks positions are given by two Index objects, tabRow and tabCol.
 *
 * If block 1 is n1xm1, block2 n2xm2, block3 n3xm3 ..., then:\n
 *  tabRow = [ n1 n1+n2 n1+n2+n3 ...] \n
 *  tabCol = [ m1 m1+m2 m1+m2+m3 ...] \n
 *
 */
class BlockMatrix {
 private:
  ACCEPT_SERIALIZATION(BlockMatrix);

  using BlocksMatrix = Eigen::Matrix<std::shared_ptr<SiconosMatrix>, Eigen::Dynamic,
                                     Eigen::Dynamic, Eigen::ColMajor>;

  /** A container of pointers to SiconosMatrix
   */
  std::shared_ptr<BlocksMatrix> _mat{nullptr};

  /** list of blocks dimension - tabRow[i] = tabRow[i-1] + ni, ni being the number of rows of
   * block i.
   */
  std::shared_ptr<std::vector<std::size_t>> _tabRow{nullptr};

  /** list of blocks dimension - tabCol[i] = tabCol[i-1] + ni, ni being the number of columns
   * of block i.
   */
  std::shared_ptr<std::vector<std::size_t>> _tabCol{nullptr};

  /** Number of rows (Warning: total number of scalar elements, not number of blocks) */
  unsigned int _dimRow{0};

  /** Number of columns (Warning: total number of scalar elements, not number of blocks) */
  unsigned int _dimCol{0};

  /** default constructor
   */
  BlockMatrix() = default;

 public:
  /** no-copy constructor
   *  \param m a SiconosMatrix
   */
  BlockMatrix(std::shared_ptr<SiconosMatrix> m);

  /** copy constructor
   *  \param m a MapMatrix
   */
  BlockMatrix(std::shared_ptr<MapType> m);

  /** Build from eigen view (shared-memory !) */
  BlockMatrix(Eigen::Ref<siconos::algebra::SiconosMatrix> input);

  // /** copy constructor
  //  *  \param m a SiconosMatrix
  //  */
  // BlockMatrix(const SiconosMatrix &m);

  /** copy constructor
   *  \param m a BlockMatrix
   */
  BlockMatrix(const BlockMatrix &m);

  /** constructor with a list of pointer to SiconosMatrix (!links with pointer, no copy!)
   *  \param m a vector of SiconosMatrix
   *  \param row number of blocks in a row
   *  \param col number of col in a row
   */
  BlockMatrix(const std::vector<std::shared_ptr<SiconosMatrix>> &m, unsigned int row,
              unsigned int col);

  /** contructor with a list of 4 pointer to SiconosMatrix (!links with pointer, no copy!)
   *  \param A block (0,0)
   *  \param B block (0,1)
   *  \param C block (1,0)
   *  \param D block (1,1)
   */
  BlockMatrix(std::shared_ptr<SiconosMatrix> A, std::shared_ptr<SiconosMatrix> B,
              std::shared_ptr<SiconosMatrix> C, std::shared_ptr<SiconosMatrix> D);

  /** destructor
   */
  ~BlockMatrix(void) noexcept;

  // inline bool checkSymmetry(double tol) const override { return false; };

  /** get the number of block (i=0, row, i=1 col)
   *  \param i unsigned int(i=0, row, i=1 col)
   *  \return an unsigned int
   */
  unsigned int numberOfBlocks(unsigned int i) const;

  /** return the address of the array of double values of the matrix
   *  \param row position for the required block ->useless for SiconosMatrix
   *  \param col position for the required block ->useless for SiconosMatrix
   *  \return double* : the pointer on the double array
   */
  double *getArray(unsigned int row = 0, unsigned int col = 0) const;

  /** sets all the values of the matrix to 0.0
   */
  void zero();

  /** Initialize the matrix with random values
   */
  void randomize();

  /** set an identity matrix
   */
  void setIdentity();

  /** \return the number of rows of the matrix */
  auto rows() const { return _dimRow; }

  /** \return the number of columns of the matrix */
  auto cols() const { return _dimCol; }

  /** compute the infinite norm of the Block matrix
   *  \return a double
   */
  double normInf() const;

  /** display data on standard output
   */
  void display() const;

  /** display data on standard output
   */
  void displayExpert(bool brief = true) const;

  friend std::ostream &operator<<(std::ostream &os, const BlockMatrix &bm);

  /** get or set the element matrix[i,j]
   *  \param i an unsigned int
   *  \param j an unsigned int
   *  \return the element matrix[i,j]
   */
  double &operator()(unsigned int i, unsigned int j);

  /** get or set the element matrix[i,j]
   *  \param i an unsigned int
   *  \param j an unsigned int
   *  \return the element matrix[i,j]
   */
  double operator()(unsigned int i, unsigned int j) const;

  /** return the element matrix[i,j]
   *  \param i an unsigned int
   *  \param j an unsigned int
   *  \return a double
   */
  double getValue(unsigned int i, unsigned int j) const;

  /** set the element matrix[i,j]
   *  \param i an unsigned int i
   *  \param j an unsigned int j
   *  \param value
   */
  void setValue(unsigned int i, unsigned int j, double value);

  /** get the vector tabRow
   *  \return a vector of int
   */
  inline std::vector<std::size_t> getTabRow() const { return *_tabRow; };

  /** get the vector tabCol
   *  \return a vector of int
   */
  inline std::vector<std::size_t> getTabCol() const { return *_tabCol; };

  /** get the vector tabRow
   *  \return a pointer to vector of int
   */
  inline const std::shared_ptr<std::vector<std::size_t>> tabRow() const { return _tabRow; };

  /** get the vector tabCol
   *  \return a pointer to vector of int
   */
  inline const std::shared_ptr<std::vector<std::size_t>> tabCol() const { return _tabCol; };

  /** get block at position row-col
   *  \param row unsigned int
   *  \param col unsigned int
   *  \return std::shared_ptr<SiconosMatrix> the requested block
   */
  std::shared_ptr<SiconosMatrix> block(unsigned int row = 0, unsigned int col = 0);

  /** get block at position row-col
   *  \param row unsigned int
   *  \param col unsigned int
   *  \return std::shared_ptr<SiconosMatrix> the requested block
   */
  std::shared_ptr<const SiconosMatrix> block(unsigned int row = 0, unsigned int col = 0) const;

  /** convert BlockMatrix to SiconosMatrix
   *  \return SiconosMatrix the converted matrix
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> toSiconosMatrix() const;

  /**
   *
   */
  void copyBlock(unsigned int i, unsigned int j,
                 std::shared_ptr<siconos::algebra::SiconosMatrix>);

  /** Set new block pointer
   *
   */
  void setBlock(unsigned int i, unsigned int j,
                std::shared_ptr<siconos::algebra::SiconosMatrix>);

  void updateNumericsMatrix() {
    THROW_EXCEPTION("BlockMatrix::updateNumericsMatrix(), not implemented fro BlockMatrix");
  };

  // friend class SiconosMatrix;
  friend SiconosMatrix &operator*=(SiconosMatrix &m, const double &s);
  friend SiconosMatrix &operator/=(SiconosMatrix &m, const double &s);

  /** number of non-zero in the matrix
   * \param tol the tolerance under which a number is considered zero
   */
  size_t nnz(double tol = 1.e-14);
};

}  // namespace siconos::algebra

#endif