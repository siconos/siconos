/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "SiconosMatrixOp.hpp"  // For setBlock, isComparableto ...
#include "SiconosAlgebraTools.hpp"  // for isComparableTo
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosException.hpp"

// =================================================
//                CONSTRUCTORS
// =================================================

siconos::algebra::BlockMatrix::BlockMatrix(const SiconosMatrix &m)
{
  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  
  _tabRow->reserve(1);
  _tabCol->reserve(1);
  // _mat construction
  _mat = std::make_shared<BlocksMatrix>(1, 1);
  _mat->setValue(0, 0, std::make_shared<SimpleMatrix>(m));

  _dimRow = m.size(0);
  _dimCol = m.size(1);
  _tabRow->push_back(_dimRow);
  _tabCol->push_back(_dimCol);
}

siconos::algebra::BlockMatrix::BlockMatrix(const BlockMatrix &m)
{
  unsigned int nbRows = m.numberOfBlocks(0);
  unsigned int nbCols = m.numberOfBlocks(1);
  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(nbRows);
  _tabCol->reserve(nbCols);

  // _mat construction
  _mat = std::make_shared<BlocksMatrix>();

  unsigned int i = 0, j = 0;
  // We scan all the blocks of m ...
  bool firstLoop = true;
  // We scan all the blocks of m ...
  for (auto row : m._mat->rowwise()) {
    _dimRow += row.size();
    _tabRow->push_back(_dimRow);
    for (auto col : row) {
        (*_mat)(i, j) = col;
      // _dimCol must be incremented only at first "column-loop"
      if (firstLoop) {
        _dimCol += col->size(1);
        _tabCol->push_back(_dimCol);
      }
      ++j;
    }
    ++i;
    firstLoop = false;
  }
}

siconos::algebra::BlockMatrix::BlockMatrix(
    const std::vector<std::shared_ptr<SiconosMatrix>> &m, unsigned int row, unsigned int col)
{
  if (m.size() != (row * col))
    THROW_EXCEPTION("number of blocks inconsistent with provided dimensions.");

  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(row);
  _tabCol->reserve(col);

  // _mat construction
  _mat = std::make_shared<BlocksMatrix>();

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
{
  if (A->size(0) != B->size(0) || C->size(0) != D->size(0) || A->size(1) != C->size(1) ||
      B->size(1) != D->size(1))
    THROW_EXCEPTION("inconsistent sizes between A, B, C or D SiconosMatrices.");

  // _mat = [ A B ]
  //       [ C D ]

  // _mat construction
  _mat = std::make_shared<BlocksMatrix>(2, 2);

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
  _mat->resize(0, 0);
  _mat = nullptr;

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


double *siconos::algebra::BlockMatrix::getArray(unsigned int i, unsigned int j) const
{
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, j);
  return tmp->data();
}

// ===========================
//       fill matrix
// ===========================

void siconos::algebra::BlockMatrix::zero()
{
  for (auto it : _mat->reshaped() ) {
    it->zero();
  }
}

void siconos::algebra::BlockMatrix::randomize()
{
  for (auto it : _mat->reshaped()) {
    // auto nrows = it->rows();
    // auto ncols = it->cols();
    it->setRandom();
  }
}

void siconos::algebra::BlockMatrix::eye()
{
  int i = 0, j = 0;
  for (auto row : _mat->rowwise()) {
    for (auto col : row) {
      if (i == j)
        col->setIdentity();
      else
        col->setZero();
      ++j;
    }
    ++i;
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
  for (auto row : _mat->rowwise()) {
    for (auto col : row) {
      col->display();
    }
  }
  std::cout << "=============================================================================="
               "=============\n";
}
void siconos::algebra::BlockMatrix::displayExpert(bool brief) const
{
  std::cout << "==========> BlockMatrix (" << numberOfBlocks(0) << " X " << numberOfBlocks(1)
            << " blocks): \n";
  for (auto row : _mat->rowwise()) {
    for (auto col : row) {
      // col->displayExpert(brief); // TODO
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
  for (auto row : bm._mat->rowwise()) {
    for (auto col : row) {
      if (col != *row.begin()) os << ",";
      if (col)
        os << *col;
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

// void siconos::algebra::BlockMatrix::getRow(unsigned int r, SiconosVector &v) const
// {
//   unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

//   if (r > _dimRow) THROW_EXCEPTION("row number is out of range");

//   // Verification of the size of the result vector
//   if (v.size() != _dimCol) THROW_EXCEPTION("inconsistent sizes");

//   // Find the row-block number where "r" is
//   while (r >= (*_tabRow)[numRow] && numRow < _tabRow->size()) numRow++;

//   // Computation of the value of the index row into this block
//   if (numRow != 0) posRow -= (*_tabRow)[numRow - 1];

//   for (unsigned int j = 0; j < _tabCol->size(); j++) {
//     start = stop;
//     std::shared_ptr<SiconosMatrix> tmp = (*_mat)(numRow, j);
//     stop += tmp->size(1);
//     boost::numeric::ublas::subrange(*(v.dense()), start, stop) =
//         boost::numeric::ublas::row(*(tmp->dense()), posRow);
//   }
// }

// void siconos::algebra::BlockMatrix::getCol(unsigned int c, SiconosVector &v) const
// {
//   unsigned int numCol = 0, posCol = c, start = 0, stop = 0;

//   if (c > _dimCol) THROW_EXCEPTION("column number is out of range");

//   // Verification of the size of the result vector
//   if (v.size() != _dimRow) THROW_EXCEPTION("inconsistent sizes");

//   // Find the column-block number where "c" is
//   while (c >= (*_tabCol)[numCol] && numCol < _tabCol->size()) numCol++;

//   // Computation of the value of the index column into this block
//   if (numCol != 0) posCol -= (*_tabCol)[numCol - 1];

//   for (unsigned int i = 0; i < _tabRow->size(); i++) {
//     start = stop;
//     std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, numCol);
//     stop += tmp->size(0);
//     boost::numeric::ublas::subrange(*(v.dense()), start, stop) =
//         boost::numeric::ublas::column(tmp->getDense(), posCol);
//   }
// }

// void siconos::algebra::BlockMatrix::setRow(unsigned int r, const SiconosVector &v)
// {
//   unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

//   if (v.size() != _dimCol) THROW_EXCEPTION("inconsistent sizes");

//   while (r >= (*_tabRow)[numRow] && numRow < _tabRow->size()) numRow++;

//   if (numRow != 0) posRow -= (*_tabRow)[numRow - 1];

//   for (unsigned int j = 0; j < _tabCol->size(); j++) {
//     start = stop;
//     std::shared_ptr<SiconosMatrix> tmp = (*_mat)(numRow, j);
//     stop += tmp->size(1);
//     boost::numeric::ublas::row(*(tmp->dense()), posRow) =
//         boost::numeric::ublas::subrange(*(v.dense()), start, stop);
//   }
// }

// void siconos::algebra::BlockMatrix::setCol(unsigned int col, const SiconosVector &v)
// {
//   unsigned int numCol = 0, posCol = col, start = 0, stop = 0;

//   if (v.size() != _dimRow) THROW_EXCEPTION("inconsistent sizes");

//   while (col >= (*_tabCol)[numCol] && numCol < _tabCol->size()) numCol++;

//   if (numCol != 0) posCol -= (*_tabCol)[numCol - 1];

//   for (unsigned int i = 0; i < _tabRow->size(); i++) {
//     start = stop;
//     std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, numCol);
//     stop += tmp->size(0);
//     boost::numeric::ublas::column(*(tmp->dense()), posCol) =
//         boost::numeric::ublas::subrange(*(v.dense()), start, stop);
//   }
// }

// void siconos::algebra::BlockMatrix::addSimple(unsigned int &indRow, unsigned int &indCol,
//                                               const SiconosMatrix &m)
// {
//   // Add a part of m (starting from (indRow,indCol) to the current matrix.
//   // m must be a SimpleMatrix.

//   // At the end of the present function, indRow (resp. indCol) is equal to indRow + the
//   // corresponding dimension of the added sub-matrix.

//   unsigned int row = m.size(0) - indRow;  // number of rows of the block to be added.
//   unsigned int col = m.size(1) - indCol;  // number of columns of the block to be added.
//   unsigned int initCol = indCol;

//   if (row > _dimRow || col > _dimCol) THROW_EXCEPTION("invalid ranges");

//   auto numM = m.num();

//   // iterators through this
//   unsigned int currentRow = 0, currentCol = 0;
//   UblasType currentNum;
//   for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//     for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//       if ((*it2)->isBlock())  // if the sub-block is also a BlockMatrix ...
//         (std::static_pointer_cast<BlockMatrix>(*it2))->addSimple(indRow, indCol, m);

//       else {
//         currentCol = (*it2)->size(1);
//         currentRow = (*it2)->size(0);
//         currentNum = (*it2)->num();
//         if (numM != currentNum) THROW_EXCEPTION("inconsistent types.");

//         if (numM == UblasType::DENSE)
//           noalias(*(*it2)->dense()) += boost::numeric::ublas::subrange(
//               *m.dense(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::TRIANGULAR)
//           noalias(*(*it2)->triang()) += boost::numeric::ublas::subrange(
//               *m.triang(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::SYMMETRIC)
//           noalias(*(*it2)->sym()) += boost::numeric::ublas::subrange(
//               *m.sym(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::SPARSE)
//           noalias(*(*it2)->sparse()) += boost::numeric::ublas::subrange(
//               *m.sparse(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::BANDED)
//           noalias(*(*it2)->banded()) += boost::numeric::ublas::subrange(
//               *m.banded(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::ZERO) {
//         }
//         else
//           THROW_EXCEPTION("inconsistent types.");
//       }
//       indCol += currentCol;
//     }
//     indRow += currentRow;
//     indCol = initCol;
//   }
// }

// void siconos::algebra::BlockMatrix::subSimple(unsigned int &indRow, unsigned int &indCol,
//                                               const SiconosMatrix &m)
// {
//   // subtract a part of m (starting from (indRow,indCol) to the current matrix.
//   // m must be a SimpleMatrix.

//   // At the end of the present function, indRow (resp. indCol) is equal to indRow + the
//   // corresponding dimension of the subtracted sub-matrix.

//   unsigned int row = m.size(0) - indRow;  // number of rows of the block to be added.
//   unsigned int col = m.size(1) - indCol;  // number of columns of the block to be added.
//   unsigned int initCol = indCol;
//   if (row > _dimRow || col > _dimCol) THROW_EXCEPTION("invalid ranges");

//   auto numM = m.num();

//   // iterators through this
//   unsigned int currentRow = 0, currentCol = 0;
//   UblasType currentNum;
//   for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//     for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//       if ((*it2)->isBlock())  // if the sub-block is also a BlockMatrix ...
//         (std::static_pointer_cast<BlockMatrix>(*it2))->subSimple(indRow, indCol, m);

//       else {
//         currentCol = (*it2)->size(1);
//         currentRow = (*it2)->size(0);
//         currentNum = (*it2)->num();
//         if (numM != currentNum) THROW_EXCEPTION("inconsistent types.");

//         if (numM == UblasType::DENSE)
//           noalias(*(*it2)->dense()) -= boost::numeric::ublas::subrange(
//               *m.dense(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::TRIANGULAR)
//           noalias(*(*it2)->triang()) -= boost::numeric::ublas::subrange(
//               *m.triang(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::SYMMETRIC)
//           noalias(*(*it2)->sym()) -= boost::numeric::ublas::subrange(
//               *m.sym(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::SPARSE)
//           noalias(*(*it2)->sparse()) -= boost::numeric::ublas::subrange(
//               *m.sparse(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::BANDED)
//           noalias(*(*it2)->banded()) -= boost::numeric::ublas::subrange(
//               *m.banded(), indRow, indRow + currentRow, indCol, indCol + currentCol);
//         else if (numM == UblasType::ZERO) {
//         }
//         else
//           THROW_EXCEPTION("inconsistent types.");
//       }
//       indCol += currentCol;
//     }
//     indRow += currentRow;
//     indCol = initCol;
//   }
// }

//===============
//  Assignment
//===============

// siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator=(const SiconosMatrix &m)
// {
//   if (&m == this) return *this;  // auto-assignment.

//   if (m.size(0) != _dimRow || m.size(1) != _dimCol)
//     THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

//   // Warning: we do not reallocate the blocks, but only copy the values. This means that
//   // all blocks are already allocated and that dim of m and mat are to be consistent.
//   // Thus, _tabRow and _tabCol remains unchanged.
//   // If m and mat are not "block-consistent", we use the () operator for a component-wise copy.

//   if (m.isBlock()) {
//     if (siconos::algebra::isComparableTo(*this, m)) {
//       const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
//       // iterators through this
//       // iterators through m
//       auto itM1 = mB._mat->begin();

//       for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//         auto itM2 = itM1.begin();
//         for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//           (**it2) = (**itM2);
//           itM2++;  // increment column pos. in m.
//         }
//         itM1++;  // increment row pos. in m.
//       }
//     }
//     else {
//       for (unsigned int i = 0; i < _dimRow; ++i)
//         for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) = m(i, j);
//     }
//   }
//   else  // if m is a SimpleMatrix
//   {
//     unsigned int posRow = 0;
//     unsigned int posCol = 0;
//     std::vector<std::size_t> subDim(2);
//     std::vector<std::size_t> subPos(4);

//     for (auto it = _mat->begin(); it != _mat->end(); ++it) {
//       for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
//         // a sub-block of m is copied into this
//         subDim[0] = (*it2)->size(0);
//         subDim[1] = (*it2)->size(1);
//         subPos[0] = posRow;
//         subPos[1] = posCol;
//         subPos[2] = 0;
//         subPos[3] = 0;
//         siconos::algebra::setBlock(m, *it2, subDim, subPos);
//         posCol += subDim[1];
//       }
//       posRow += (*it)->size(0);
//       posCol = 0;
//     }
//   }

//   return *this;
// }

// siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator=(const BlockMatrix &m)
// {
//   if (&m == this) return *this;  // auto-assignment.

//   if (m.size(0) != _dimRow || m.size(1) != _dimCol)
//     THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

//   // Warning: we do not reallocate the blocks, but only copy the values. This means that
//   // all blocks are already allocated and that dim of m and mat are to be consistent.
//   // Thus, _tabRow and _tabCol remains unchanged.
//   // If m and mat are not "block-consistent", we use the () operator for a componet-wise copy.

//   if (siconos::algebra::isComparableTo(*this, m)) {
//     const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
//     // iterators through this
//     // iterators through m
//     auto itM1 = mB._mat->begin();

//     for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//       auto itM2 = itM1.begin();
//       for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//         (**it2) = (**itM2);
//         itM2++;  // increment column pos. in m.
//       }
//       itM1++;  // increment row pos. in m.
//     }
//   }
//   else {
//     for (unsigned int i = 0; i < _dimRow; ++i)
//       for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) = m(i, j);
//   }
//   return *this;
// }


//=================================
// Op. and assignment (+=, -= ... )
//=================================

// siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator+=(
//     const SiconosMatrix &m)
// {
//   if (&m == this) {
//     for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//       for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//         **it2 += **it2;
//       }
//     }
//     return *this;
//   }

//   if (m.size(0) != _dimRow || m.size(1) != _dimCol)
//     THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

//   if (m.isBlock()) {
//     if (siconos::algebra::isComparableTo(m, *this)) {
//       const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
//       // iterators through m
//       auto itM1 = mB._mat->begin();

//       for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//         auto itM2 = itM1.begin();
//         for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//           (**it2) += (**itM2);
//           itM2++;  // increment column pos. in m.
//         }
//         itM1++;  // increment row pos. in m.
//       }
//     }
//     else {
//       for (unsigned int i = 0; i < _dimRow; ++i)
//         for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) += m(i, j);
//     }
//   }
//   else  // if m is a SimpleMatrix
//   {
//     unsigned int indRow = 0, indCol = 0;
//     addSimple(indRow, indCol, m);  // a sub-block of m is added to each block of this.
//   }
//   return *this;
// }

// siconos::algebra::BlockMatrix &siconos::algebra::BlockMatrix::operator-=(
//     const SiconosMatrix &m)
// {
//   if (&m == this) {
//     for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//       for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//         **it2 -= **it2;
//       }
//     }
//     return *this;
//   }

//   if (m.size(0) != _dimRow || m.size(1) != _dimCol)
//     THROW_EXCEPTION("Left and Right values have inconsistent sizes.");

//   if (m.isBlock()) {
//     if (siconos::algebra::isComparableTo(m, *this)) {
//       const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
//       // iterators through m
//       auto itM1 = mB._mat->begin();

//       for (auto it1 = _mat->begin(); it1 != _mat->end(); ++it1) {
//         auto itM2 = itM1.begin();
//         for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
//           (**it2) -= (**itM2);
//           itM2++;  // increment column pos. in m.
//         }
//         itM1++;  // increment row pos. in m.
//       }
//     }
//     else {
//       for (unsigned int i = 0; i < _dimRow; ++i)
//         for (unsigned int j = 0; j < _dimCol; ++j) (*this)(i, j) -= m(i, j);
//     }
//   }
//   else  // if m is a SimpleMatrix
//   {
//     unsigned int indRow = 0, indCol = 0;
//     subSimple(indRow, indCol, m);  // a sub-block of m is subtracted to each block of this.
//   }
//   return *this;
// }

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

std::shared_ptr<siconos::algebra::SiconosMatrix> siconos::algebra::BlockMatrix::toSiconosMatrix() const
{
  // get number of blocks in a row/col of m.
  auto m = std::make_shared<SiconosMatrix>(this->size(0), this->size(1));
  unsigned int posRow = 0;
  unsigned int posCol = 0;

  for (auto it : this->_mat->rowwise()) {
    for (auto it2 : it) {
      m->setBlock(posRow, posCol, *it2);
      posCol += it2->size(1);
    }
    posRow += it.size();
    posCol = 0;
  }
  return m;
}

size_t siconos::algebra::BlockMatrix::nnz(double tol)
{
  size_t nnz = 0;
  for (auto row : _mat->rowwise()) {
    for (auto col : row) nnz += (*col).nnz();
  }
  return nnz;
}
