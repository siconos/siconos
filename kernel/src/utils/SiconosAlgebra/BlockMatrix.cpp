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

#include "BlockMatrix.hpp"
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>  // subrange ...
#include "SiconosMatrixOp.hpp"  // For setBlock, isComparableto ...
#include "SiconosAlgebraTools.hpp"  // for isComparableTo
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

// =================================================
//                CONSTRUCTORS
// =================================================

siconos::algebra::BlockMatrix::BlockMatrix(const SiconosMatrix &m)
    : SiconosMatrix(UblasType::BLOCK), _dimRow(0), _dimCol(0)
{
  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  if (m.isBlock()) {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    unsigned int nbRows = m.numberOfBlocks(0);
    unsigned int nbCols = m.numberOfBlocks(1);
    _tabRow->reserve(nbRows);
    _tabCol->reserve(nbCols);

    // mat construction
    _mat = std::make_shared<BlocksMatrix>(nbRows, nbCols, nbRows * nbCols);

    unsigned int i, j;
    bool firstLoop = true;
    // We scan all the blocks of m ...
    for (auto it1 = mB._mat->begin1(); it1 != mB._mat->end1(); ++it1) {
      _dimRow += (*(it1.begin()))->size(0);
      _tabRow->push_back(_dimRow);
      for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
        i = it2.index1();
        j = it2.index2();
        if ((*it2)->isBlock())  // if the current matrix is a blockMatrix
          _mat->insert_element(i, j, std::make_shared<BlockMatrix>(**it2));
        else
          _mat->insert_element(i, j, std::make_shared<SimpleMatrix>(**it2));
        // _dimCol must be incremented only at first "column-loop"
        if (firstLoop) {
          _dimCol += (*it2)->size(1);
          _tabCol->push_back(_dimCol);
        }
      }
      firstLoop = false;
    }
  }
  else  // if m is a SimpleMatrix
  {
    _tabRow->reserve(1);
    _tabCol->reserve(1);
    // _mat construction
    _mat = std::make_shared<BlocksMatrix>(1, 1, 1);
    _mat->insert_element(0, 0, std::make_shared<SimpleMatrix>(m));

    _dimRow = m.size(0);
    _dimCol = m.size(1);
    _tabRow->push_back(_dimRow);
    _tabCol->push_back(_dimCol);
  }
}

siconos::algebra::BlockMatrix::BlockMatrix(const BlockMatrix &m)
    : SiconosMatrix(UblasType::BLOCK), _dimRow(0), _dimCol(0)
{
  unsigned int nbRows = m.numberOfBlocks(0);
  unsigned int nbCols = m.numberOfBlocks(1);
  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(nbRows);
  _tabCol->reserve(nbCols);

  // _mat construction
  _mat = std::make_shared<BlocksMatrix>(nbRows, nbCols, nbRows * nbCols);

  unsigned int i, j;
  // We scan all the blocks of m ...
  bool firstLoop = true;
  // We scan all the blocks of m ...
  for (auto it1 = m._mat->begin1(); it1 != m._mat->end1(); ++it1) {
    _dimRow += (*(it1.begin()))->size(0);
    _tabRow->push_back(_dimRow);
    for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
      i = it2.index1();
      j = it2.index2();
      if ((*it2)->isBlock())  // if the current _matrix is a blockMatrix
        _mat->insert_element(i, j, std::make_shared<BlockMatrix>(**it2));
      else
        _mat->insert_element(i, j, std::make_shared<SimpleMatrix>(**it2));

      // _dimCol must be incremented only at first "column-loop"
      if (firstLoop) {
        _dimCol += (*it2)->size(1);
        _tabCol->push_back(_dimCol);
      }
    }
    firstLoop = false;
  }
}

siconos::algebra::BlockMatrix::BlockMatrix(
    const std::vector<std::shared_ptr<SiconosMatrix>> &m, unsigned int row, unsigned int col)
    : SiconosMatrix(UblasType::BLOCK), _dimRow(0), _dimCol(0)
{
  if (m.size() != (row * col))
    THROW_EXCEPTION("number of blocks inconsistent with provided dimensions.");

  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(row);
  _tabCol->reserve(col);

  // _mat construction
  _mat = std::make_shared<BlocksMatrix>(row, col, row * col);

  unsigned int k = 0;
  bool firstRowLoop = true;
  bool firstColLoop = true;

  for (unsigned int i = 0; i < row; ++i) {
    for (unsigned int j = 0; j < col; ++j) {
      (*_mat)(i, j) = m[k++];

      // _dimCol must be incremented only at first "column-loop"
      if (firstColLoop) {
        _dimCol += m[k - 1]->size(1);
        _tabCol->push_back(_dimCol);
      }
      if (firstRowLoop) {
        _dimRow += m[k - 1]->size(0);
        _tabRow->push_back(_dimRow);
        firstRowLoop = false;
      }
    }
    firstColLoop = false;
    firstRowLoop = true;
  }
}

siconos::algebra::BlockMatrix::BlockMatrix(std::shared_ptr<SiconosMatrix> A,
                                           std::shared_ptr<SiconosMatrix> B,
                                           std::shared_ptr<SiconosMatrix> C,
                                           std::shared_ptr<SiconosMatrix> D)
    : SiconosMatrix(UblasType::BLOCK), _dimRow(0), _dimCol(0)
{
  if (A->size(0) != B->size(0) || C->size(0) != D->size(0) || A->size(1) != C->size(1) ||
      B->size(1) != D->size(1))
    THROW_EXCEPTION("inconsistent sizes between A, B, C or D SiconosMatrices.");

  // _mat = [ A B ]
  //       [ C D ]

  // _mat construction
  _mat = std::make_shared<BlocksMatrix>(2, 2, 4);

  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(2);
  _tabCol->reserve(2);

  (*_mat)(0, 0) = A;
  (*_mat)(0, 1) = B;
  (*_mat)(1, 0) = C;
  (*_mat)(1, 1) = D;
  _dimRow = A->size(0);
  _tabRow->push_back(_dimRow);
  _dimRow += C->size(0);
  _tabRow->push_back(_dimRow);
  _dimCol = A->size(1);
  _tabCol->push_back(_dimCol);
  _dimCol += B->size(1);
  _tabCol->push_back(_dimCol);
}

siconos::algebra::BlockMatrix::~BlockMatrix() noexcept
{
  _mat->clear();

  _tabRow->clear();
  _tabCol->clear();
}

// =================================================
//    get number of blocks
// =================================================

unsigned int siconos::algebra::BlockMatrix::numberOfBlocks(unsigned int dim) const
{
  if (dim == 0)
    return _tabRow->size();
  else
    return _tabCol->size();
}

// =================================================
//        get Ublas component (dense ...)
// =================================================

// return the boost dense _matrix of the block (i, j)
const siconos::algebra::DenseMat siconos::algebra::BlockMatrix::getDense(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);

  if (tmp->num() != UblasType::DENSE)
    THROW_EXCEPTION("the matrix at (row, col) is not a Dense matrix");

  return (tmp->getDense());
}

// return the boost triangular matrix of the block (i, j)
const siconos::algebra::TriangMat siconos::algebra::BlockMatrix::getTriang(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::TRIANGULAR) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Triangular matrix");
  }
  return (tmp->getTriang());
}

// return the boost symmetric matrix of the block (i, j)
const siconos::algebra::SymMat siconos::algebra::BlockMatrix::getSym(unsigned int row,
                                                                     unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::SYMMETRIC) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Symmmetric matrix");
  }
  return (tmp->getSym());
}

// return the boost sparse matrix of the block (i, j)
const siconos::algebra::SparseMat siconos::algebra::BlockMatrix::getSparse(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::SPARSE) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Sparse matrix");
  }
  return (tmp->getSparse());
}

// return the boost sparse matrix of the block (i, j)
const siconos::algebra::SparseCoordinateMat siconos::algebra::BlockMatrix::getSparseCoordinate(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::SPARSE_COORDINATE) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Sparse matrix");
  }
  return (tmp->getSparseCoordinate());
}
// return the boost banded matrix of the block (i, j)
const siconos::algebra::BandedMat siconos::algebra::BlockMatrix::getBanded(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::BANDED) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Banded matrix");
  }
  return (tmp->getBanded());
}

// return the boost zero matrix of the block (i, j)
const siconos::algebra::ZeroMat siconos::algebra::BlockMatrix::getZero(unsigned int row,
                                                                       unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::ZERO) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Zero matrix");
  }
  return (tmp->getZero());
}

// return the boost identity matrix of the block (i, j)
const siconos::algebra::IdentityMat siconos::algebra::BlockMatrix::getIdentity(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::IDENTITY) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Identity matrix");
  }
  return (tmp->getIdentity());
}

// The following functions return the corresponding pointers
siconos::algebra::DenseMat *siconos::algebra::BlockMatrix::dense(unsigned int row,
                                                                 unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::DENSE) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Dense matrix");
  }

  return (tmp->dense());
}

siconos::algebra::TriangMat *siconos::algebra::BlockMatrix::triang(unsigned int row,
                                                                   unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::TRIANGULAR) {
    THROW_EXCEPTION("The matrix at (row, col) is not a Triangular matrix");
  }
  return (tmp->triang());
}
siconos::algebra::SymMat *siconos::algebra::BlockMatrix::sym(unsigned int row,
                                                             unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::SYMMETRIC) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Symmmetric matrix");
  }
  return (tmp->sym());
}

siconos::algebra::SparseMat *siconos::algebra::BlockMatrix::sparse(unsigned int row,
                                                                   unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::SPARSE) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Sparse matrix");
  }
  return (tmp->sparse());
}
siconos::algebra::SparseCoordinateMat *siconos::algebra::BlockMatrix::sparseCoordinate(
    unsigned int row, unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::SPARSE_COORDINATE) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Sparse coordinate matrix");
  }
  return (tmp->sparseCoordinate());
}

siconos::algebra::BandedMat *siconos::algebra::BlockMatrix::banded(unsigned int row,
                                                                   unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::BANDED) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Banded matrix");
  }
  return (tmp->banded());
}

siconos::algebra::ZeroMat *siconos::algebra::BlockMatrix::zero_mat(unsigned int row,
                                                                   unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::ZERO) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Zero matrix");
  }
  return (tmp->zero_mat(row, col));
}

siconos::algebra::IdentityMat *siconos::algebra::BlockMatrix::identity(unsigned int row,
                                                                       unsigned int col) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(row, col);
  if (tmp->num() != UblasType::IDENTITY) {
    THROW_EXCEPTION("the matrix at (row, col) is not a Identity matrix");
  }
  return (tmp->identity());
}

double *siconos::algebra::BlockMatrix::getArray(unsigned int i, unsigned int j) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, j);
  return tmp->getArray();
}

// ===========================
//       fill matrix
// ===========================

void siconos::algebra::BlockMatrix::zero()
{
  for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
      (*it2)->zero();
    }
  }
}

void siconos::algebra::BlockMatrix::randomize()
{
  for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
      (*it2)->randomize();
    }
  }
}

void siconos::algebra::BlockMatrix::eye()
{
  for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
      if (it2.index1() == it2.index2())
        (*it2)->eye();
      else
        (*it2)->zero();
    }
  }
}

unsigned int siconos::algebra::BlockMatrix::size(unsigned int index) const
{
  if (index == 0)
    return _dimRow;
  else
    return _dimCol;
};

//=======================
// set matrix dimensions
//=======================

void siconos::algebra::BlockMatrix::resize(unsigned int, unsigned int, unsigned int,
                                           unsigned int, bool)
{
  THROW_EXCEPTION("forbidden for block matrices.");
}

//=======================
//       get norm
//=======================

double siconos::algebra::BlockMatrix::normInf() const
{
  double sum = 0, norm = 0;
  for (unsigned int i = 0; i < size(0); i++) {
    for (unsigned int j = 0; j < size(1); j++) {
      sum += (*this)(i, j);
    }
    if (fabs(sum) > norm) norm = fabs(sum);
    sum = 0;
  }
  return norm;
}

//=====================
// screen display
//=====================

void siconos::algebra::BlockMatrix::display(void) const
{
  std::cout << "==========> BlockMatrix (" << numberOfBlocks(0) << " X " << numberOfBlocks(1)
            << " blocks): \n";
  for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
      (*it2)->display();
    }
  }
  std::cout << "=============================================================================="
               "=============\n";
}
void siconos::algebra::BlockMatrix::displayExpert(bool brief) const
{
  std::cout << "==========> BlockMatrix (" << numberOfBlocks(0) << " X " << numberOfBlocks(1)
            << " blocks): \n";
  for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
      (*it2)->displayExpert(brief);
    }
  }
  std::cout << "=============================================================================="
               "=============\n";
}

//=====================
// convert to an ostream
//=====================

std::ostream &siconos::algebra::operator<<(std::ostream &os, const BlockMatrix &bm)
{
  os << "[" << bm.numberOfBlocks(0) << "," << bm.numberOfBlocks(1) << "](";
  for (auto it = bm._mat->begin1(); it != bm._mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
      if (it2 != it.begin()) os << ",";
      if (*it2)
        os << **it2;
      else
        os << "(nil)";
    }
  }
  os << ")";
  return os;
}

//=============================
// Elements access (get or set)
//=============================

double &siconos::algebra::BlockMatrix::operator()(unsigned int row, unsigned int col)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size()) nbRow++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size()) nbCol++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0) posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0) posCol -= (*_tabCol)[nbCol - 1];

  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(nbRow, nbCol);
  return (*tmp)(posRow, posCol);
}

double siconos::algebra::BlockMatrix::operator()(unsigned int row, unsigned int col) const
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size()) nbRow++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size()) nbCol++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0) posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0) posCol -= (*_tabCol)[nbCol - 1];

  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(nbRow, nbCol);
  return (*tmp)(posRow, posCol);
}

double siconos::algebra::BlockMatrix::getValue(unsigned int row, unsigned int col) const
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size()) nbRow++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size()) nbCol++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0) posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0) posCol -= (*_tabCol)[nbCol - 1];

  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(nbRow, nbCol);
  return (*tmp)(posRow, posCol);
}

void siconos::algebra::BlockMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size()) nbRow++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size()) nbCol++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0) posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0) posCol -= (*_tabCol)[nbCol - 1];

  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(nbRow, nbCol);
  (*tmp)(posRow, posCol) = value;
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void siconos::algebra::BlockMatrix::getRow(unsigned int r, SiconosVector &v) const
{
  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (r > _dimRow) THROW_EXCEPTION("row number is out of range");

  // Verification of the size of the result vector
  if (v.size() != _dimCol) THROW_EXCEPTION("inconsistent sizes");

  // Find the row-block number where "r" is
  while (r >= (*_tabRow)[numRow] && numRow < _tabRow->size()) numRow++;

  // Computation of the value of the index row into this block
  if (numRow != 0) posRow -= (*_tabRow)[numRow - 1];

  for (unsigned int j = 0; j < _tabCol->size(); j++) {
    start = stop;
    std::shared_ptr<SiconosMatrix> tmp = (*_mat)(numRow, j);
    stop += tmp->size(1);
    boost::numeric::ublas::subrange(*(v.dense()), start, stop) =
        boost::numeric::ublas::row(*(tmp->dense()), posRow);
  }
}

void siconos::algebra::BlockMatrix::getCol(unsigned int c, SiconosVector &v) const
{
  unsigned int numCol = 0, posCol = c, start = 0, stop = 0;

  if (c > _dimCol) THROW_EXCEPTION("column number is out of range");

  // Verification of the size of the result vector
  if (v.size() != _dimRow) THROW_EXCEPTION("inconsistent sizes");

  // Find the column-block number where "c" is
  while (c >= (*_tabCol)[numCol] && numCol < _tabCol->size()) numCol++;

  // Computation of the value of the index column into this block
  if (numCol != 0) posCol -= (*_tabCol)[numCol - 1];

  for (unsigned int i = 0; i < _tabRow->size(); i++) {
    start = stop;
    std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, numCol);
    stop += tmp->size(0);
    boost::numeric::ublas::subrange(*(v.dense()), start, stop) =
        boost::numeric::ublas::column(tmp->getDense(), posCol);
  }
}

void siconos::algebra::BlockMatrix::setRow(unsigned int r, const SiconosVector &v)
{
  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (v.size() != _dimCol) THROW_EXCEPTION("inconsistent sizes");

  while (r >= (*_tabRow)[numRow] && numRow < _tabRow->size()) numRow++;

  if (numRow != 0) posRow -= (*_tabRow)[numRow - 1];

  for (unsigned int j = 0; j < _tabCol->size(); j++) {
    start = stop;
    std::shared_ptr<SiconosMatrix> tmp = (*_mat)(numRow, j);
    stop += tmp->size(1);
    boost::numeric::ublas::row(*(tmp->dense()), posRow) =
        boost::numeric::ublas::subrange(*(v.dense()), start, stop);
  }
}

void siconos::algebra::BlockMatrix::setCol(unsigned int col, const SiconosVector &v)
{
  unsigned int numCol = 0, posCol = col, start = 0, stop = 0;

  if (v.size() != _dimRow) THROW_EXCEPTION("inconsistent sizes");

  while (col >= (*_tabCol)[numCol] && numCol < _tabCol->size()) numCol++;

  if (numCol != 0) posCol -= (*_tabCol)[numCol - 1];

  for (unsigned int i = 0; i < _tabRow->size(); i++) {
    start = stop;
    std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, numCol);
    stop += tmp->size(0);
    boost::numeric::ublas::column(*(tmp->dense()), posCol) =
        boost::numeric::ublas::subrange(*(v.dense()), start, stop);
  }
}

void siconos::algebra::BlockMatrix::addSimple(unsigned int &indRow, unsigned int &indCol,
                                              const SiconosMatrix &m)
{
  // Add a part of m (starting from (indRow,indCol) to the current matrix.
  // m must be a SimpleMatrix.

  // At the end of the present function, indRow (resp. indCol) is equal to indRow + the
  // corresponding dimension of the added sub-matrix.

  unsigned int row = m.size(0) - indRow;  // number of rows of the block to be added.
  unsigned int col = m.size(1) - indCol;  // number of columns of the block to be added.
  unsigned int initCol = indCol;

  if (row > _dimRow || col > _dimCol) THROW_EXCEPTION("invalid ranges");

  auto numM = m.num();

  // iterators through this
  unsigned int currentRow = 0, currentCol = 0;
  UblasType currentNum;
  for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
    for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
      if ((*it2)->isBlock())  // if the sub-block is also a BlockMatrix ...
        (std::static_pointer_cast<BlockMatrix>(*it2))->addSimple(indRow, indCol, m);

      else {
        currentCol = (*it2)->size(1);
        currentRow = (*it2)->size(0);
        currentNum = (*it2)->num();
        if (numM != currentNum) THROW_EXCEPTION("inconsistent types.");

        if (numM == UblasType::DENSE)
          noalias(*(*it2)->dense()) += boost::numeric::ublas::subrange(
              *m.dense(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::TRIANGULAR)
          noalias(*(*it2)->triang()) += boost::numeric::ublas::subrange(
              *m.triang(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::SYMMETRIC)
          noalias(*(*it2)->sym()) += boost::numeric::ublas::subrange(
              *m.sym(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::SPARSE)
          noalias(*(*it2)->sparse()) += boost::numeric::ublas::subrange(
              *m.sparse(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::BANDED)
          noalias(*(*it2)->banded()) += boost::numeric::ublas::subrange(
              *m.banded(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::ZERO) {
        }
        else
          THROW_EXCEPTION("inconsistent types.");
      }
      indCol += currentCol;
    }
    indRow += currentRow;
    indCol = initCol;
  }
}

void siconos::algebra::BlockMatrix::subSimple(unsigned int &indRow, unsigned int &indCol,
                                              const SiconosMatrix &m)
{
  // subtract a part of m (starting from (indRow,indCol) to the current matrix.
  // m must be a SimpleMatrix.

  // At the end of the present function, indRow (resp. indCol) is equal to indRow + the
  // corresponding dimension of the subtracted sub-matrix.

  unsigned int row = m.size(0) - indRow;  // number of rows of the block to be added.
  unsigned int col = m.size(1) - indCol;  // number of columns of the block to be added.
  unsigned int initCol = indCol;
  if (row > _dimRow || col > _dimCol) THROW_EXCEPTION("invalid ranges");

  auto numM = m.num();

  // iterators through this
  unsigned int currentRow = 0, currentCol = 0;
  UblasType currentNum;
  for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
    for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
      if ((*it2)->isBlock())  // if the sub-block is also a BlockMatrix ...
        (std::static_pointer_cast<BlockMatrix>(*it2))->subSimple(indRow, indCol, m);

      else {
        currentCol = (*it2)->size(1);
        currentRow = (*it2)->size(0);
        currentNum = (*it2)->num();
        if (numM != currentNum) THROW_EXCEPTION("inconsistent types.");

        if (numM == UblasType::DENSE)
          noalias(*(*it2)->dense()) -= boost::numeric::ublas::subrange(
              *m.dense(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::TRIANGULAR)
          noalias(*(*it2)->triang()) -= boost::numeric::ublas::subrange(
              *m.triang(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::SYMMETRIC)
          noalias(*(*it2)->sym()) -= boost::numeric::ublas::subrange(
              *m.sym(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::SPARSE)
          noalias(*(*it2)->sparse()) -= boost::numeric::ublas::subrange(
              *m.sparse(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::BANDED)
          noalias(*(*it2)->banded()) -= boost::numeric::ublas::subrange(
              *m.banded(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == UblasType::ZERO) {
        }
        else
          THROW_EXCEPTION("inconsistent types.");
      }
      indCol += currentCol;
    }
    indRow += currentRow;
    indCol = initCol;
  }
}

//===============
//  Assignment
//===============

siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator=(const SiconosMatrix &m)
{
  if (&m == this) return *this;  // auto-assignment.

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

  // Warning: we do not reallocate the blocks, but only copy the values. This means that
  // all blocks are already allocated and that dim of m and mat are to be consistent.
  // Thus, _tabRow and _tabCol remains unchanged.
  // If m and mat are not "block-consistent", we use the () operator for a component-wise copy.

  if (m.isBlock()) {
    if (siconos::algebra::isComparableTo(*this, m)) {
      const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
      // iterators through this
      // iterators through m
      auto itM1 = mB._mat->begin1();

      for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
        auto itM2 = itM1.begin();
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
          (**it2) = (**itM2);
          itM2++;  // increment column pos. in m.
        }
        itM1++;  // increment row pos. in m.
      }
    }
    else {
      for (unsigned int i = 0; i < _dimRow; ++i)
        for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) = m(i, j);
    }
  }
  else  // if m is a SimpleMatrix
  {
    unsigned int posRow = 0;
    unsigned int posCol = 0;
    std::vector<std::size_t> subDim(2);
    std::vector<std::size_t> subPos(4);

    for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
      for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
        // a sub-block of m is copied into this
        subDim[0] = (*it2)->size(0);
        subDim[1] = (*it2)->size(1);
        subPos[0] = posRow;
        subPos[1] = posCol;
        subPos[2] = 0;
        subPos[3] = 0;
        siconos::algebra::setBlock(m, *it2, subDim, subPos);
        posCol += subDim[1];
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }

  return *this;
}

siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator=(const BlockMatrix &m)
{
  if (&m == this) return *this;  // auto-assignment.

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

  // Warning: we do not reallocate the blocks, but only copy the values. This means that
  // all blocks are already allocated and that dim of m and mat are to be consistent.
  // Thus, _tabRow and _tabCol remains unchanged.
  // If m and mat are not "block-consistent", we use the () operator for a componet-wise copy.

  if (siconos::algebra::isComparableTo(*this, m)) {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    // iterators through this
    // iterators through m
    auto itM1 = mB._mat->begin1();

    for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
      auto itM2 = itM1.begin();
      for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
        (**it2) = (**itM2);
        itM2++;  // increment column pos. in m.
      }
      itM1++;  // increment row pos. in m.
    }
  }
  else {
    for (unsigned int i = 0; i < _dimRow; ++i)
      for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) = m(i, j);
  }
  return *this;
}

siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator=(const DenseMat &m)
{
  THROW_EXCEPTION("Not yet implemented.");
  return *this;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator+=(
    const SiconosMatrix &m)
{
  if (&m == this) {
    for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
      for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
        **it2 += **it2;
      }
    }
    return *this;
  }

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

  if (m.isBlock()) {
    if (siconos::algebra::isComparableTo(m, *this)) {
      const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
      // iterators through m
      auto itM1 = mB._mat->begin1();

      for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
        auto itM2 = itM1.begin();
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
          (**it2) += (**itM2);
          itM2++;  // increment column pos. in m.
        }
        itM1++;  // increment row pos. in m.
      }
    }
    else {
      for (unsigned int i = 0; i < _dimRow; ++i)
        for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) += m(i, j);
    }
  }
  else  // if m is a SimpleMatrix
  {
    unsigned int indRow = 0, indCol = 0;
    addSimple(indRow, indCol, m);  // a sub-block of m is added to each block of this.
  }
  return *this;
}

siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator-=(
    const SiconosMatrix &m)
{
  if (&m == this) {
    for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
      for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
        **it2 -= **it2;
      }
    }
    return *this;
  }

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

  if (m.isBlock()) {
    if (siconos::algebra::isComparableTo(m, *this)) {
      const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
      // iterators through m
      auto itM1 = mB._mat->begin1();

      for (auto it1 = _mat->begin1(); it1 != _mat->end1(); ++it1) {
        auto itM2 = itM1.begin();
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
          (**it2) -= (**itM2);
          itM2++;  // increment column pos. in m.
        }
        itM1++;  // increment row pos. in m.
      }
    }
    else {
      for (unsigned int i = 0; i < _dimRow; ++i)
        for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) -= m(i, j);
    }
  }
  else  // if m is a SimpleMatrix
  {
    unsigned int indRow = 0, indCol = 0;
    subSimple(indRow, indCol, m);  // a sub-block of m is subtracted to each block of this.
  }
  return *this;
}

void siconos::algebra::BlockMatrix::trans() { THROW_EXCEPTION("not yet implemented."); }

void siconos::algebra::BlockMatrix::trans(const SiconosMatrix &m)
{
  THROW_EXCEPTION("not yet implemented.");
}

void siconos::algebra::BlockMatrix::PLUFactorizationInPlace()
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}
void siconos::algebra::BlockMatrix::Factorize()
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}

void siconos::algebra::BlockMatrix::PLUInverseInPlace()
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}

void siconos::algebra::BlockMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}
void siconos::algebra::BlockMatrix::Solve(SiconosMatrix &B)
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}

void siconos::algebra::BlockMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}
void siconos::algebra::BlockMatrix::Solve(SiconosVector &B)
{
  THROW_EXCEPTION("not yet implemented for Block Matrices.");
}

std::shared_ptr<siconos::algebra::SiconosMatrix> siconos::algebra::BlockMatrix::block(
    unsigned int row, unsigned int col)
{
  return (*_mat)(row, col);
}

std::shared_ptr<const siconos::algebra::SiconosMatrix> siconos::algebra::BlockMatrix::block(
    unsigned int row, unsigned int col) const
{
  return std::shared_ptr<SiconosMatrix>((*_mat)(row, col));
}

size_t siconos::algebra::BlockMatrix::nnz(double tol)
{
  size_t nnz = 0;
  for (auto it = _mat->begin1(); it != _mat->end1(); ++it) {
    for (auto it2 = it.begin(); it2 != it.end(); ++it2) nnz += (**it2).nnz();
  }
  return nnz;
}
