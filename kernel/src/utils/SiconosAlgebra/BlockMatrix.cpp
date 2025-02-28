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

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // For setBlock, isComparableto ...
#include "SiconosVector.hpp"

// =================================================
//                CONSTRUCTORS
// =================================================

siconos::algebra::BlockMatrix::BlockMatrix(std::shared_ptr<SiconosMatrix> m) {
  _mat = std::make_shared<BlocksMatrix>(1, 1);

  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(1);
  _tabCol->reserve(1);

  (*_mat)(0, 0) = m;
  _dimRow = m->rows();
  _tabRow->push_back(_dimRow);
  _dimCol = m->cols();
  _tabCol->push_back(_dimCol);
}

siconos::algebra::BlockMatrix::BlockMatrix(std::shared_ptr<MapType> m) {
  _mat = std::make_shared<BlocksMatrix>(1, 1);

  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(1);
  _tabCol->reserve(1);

  m->display();
  (*_mat)(0, 0) = std::make_shared<siconos::algebra::SiconosMatrix>();
  *((*_mat)(0, 0)) = *m;  // TODOSAM : copy here, should not
  _dimRow = m->rows();
  _tabRow->push_back(_dimRow);
  _dimCol = m->cols();
  _tabCol->push_back(_dimCol);
}

siconos::algebra::BlockMatrix::BlockMatrix(Eigen::Ref<siconos::algebra::SiconosMatrix> input) {
  // Initialize _mat
  _mat = std::make_shared<BlocksMatrix>(1, 1);
  _tabRow = std::make_shared<std::vector<std::size_t>>();
  _tabCol = std::make_shared<std::vector<std::size_t>>();
  _tabRow->reserve(1);
  _tabCol->reserve(1);

  // Set block with memory shared (eigen map) with input matrix m
  (*_mat)(0, 0) = std::make_shared<siconos::algebra::SiconosMatrix>(
      Eigen::Map<SiconosMatrix>(input.data(), input.rows(), input.cols()));
  _dimRow = input.rows();
  _tabRow->push_back(_dimRow);
  _dimCol = input.cols();
  _tabCol->push_back(_dimCol);
}

// siconos::algebra::BlockMatrix::BlockMatrix(const SiconosMatrix &m) {
//   _tabRow = std::make_shared<std::vector<std::size_t>>();
//   _tabCol = std::make_shared<std::vector<std::size_t>>();

//   _tabRow->reserve(1);
//   _tabCol->reserve(1);
//   // _mat construction
//   _mat = std::make_shared<BlocksMatrix>(1, 1);
//   _mat->setValue(0, 0, std::make_shared<SiconosMatrix>(m));

//   _dimRow = m.rows();
//   _dimCol = m.cols();
//   _tabRow->push_back(_dimRow);
//   _tabCol->push_back(_dimCol);
// }

siconos::algebra::BlockMatrix::BlockMatrix(const BlockMatrix &m) {
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
        _dimCol += col->cols();
        _tabCol->push_back(_dimCol);
      }
      ++j;
    }
    ++i;
    firstLoop = false;
  }
}

siconos::algebra::BlockMatrix::BlockMatrix(
    const std::vector<std::shared_ptr<SiconosMatrix>> &m, unsigned int row, unsigned int col) {
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
        _dimCol += m[k - 1]->cols();
        _tabCol->push_back(_dimCol);
      }
      if (firstRowLoop) {
        _dimRow += m[k - 1]->rows();
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
                                           std::shared_ptr<SiconosMatrix> D) {
  if (A->rows() != B->rows() || C->rows() != D->rows() || A->cols() != C->cols() ||
      B->cols() != D->cols())
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
  _dimRow = A->rows();
  _tabRow->push_back(_dimRow);
  _dimRow += C->rows();
  _tabRow->push_back(_dimRow);
  _dimCol = A->cols();
  _tabCol->push_back(_dimCol);
  _dimCol += B->cols();
  _tabCol->push_back(_dimCol);
}

siconos::algebra::BlockMatrix::~BlockMatrix() noexcept {
  _mat->resize(0, 0);
  _mat = nullptr;

  _tabRow->clear();
  _tabCol->clear();
}

// =================================================
//    get number of blocks
// =================================================

unsigned int siconos::algebra::BlockMatrix::numberOfBlocks(unsigned int dim) const {
  if (dim == 0)
    return _tabRow->size();
  else
    return _tabCol->size();
}

// =================================================
//        get Ublas component (dense ...)
// =================================================

double *siconos::algebra::BlockMatrix::getArray(unsigned int i, unsigned int j) const {
  std::shared_ptr<SiconosMatrix> tmp = (*_mat)(i, j);
  return tmp->data();
}

// ===========================
//       fill matrix
// ===========================

void siconos::algebra::BlockMatrix::zero() {
  for (auto it : _mat->reshaped()) {
    it->setZero();
  }
}

void siconos::algebra::BlockMatrix::randomize() {
  for (auto it : _mat->reshaped()) {
    // auto nrows = it->rows();
    // auto ncols = it->cols();
    it->setRandom();
  }
}

void siconos::algebra::BlockMatrix::setIdentity() {
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

//=======================
//       get norm
//=======================

double siconos::algebra::BlockMatrix::normInf() const {
  double sum = 0, norm = 0;
  for (unsigned int i = 0; i < rows(); i++) {
    for (unsigned int j = 0; j < cols(); j++) {
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

void siconos::algebra::BlockMatrix::display(void) const {
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
void siconos::algebra::BlockMatrix::displayExpert(bool brief) const {
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

std::ostream &siconos::algebra::operator<<(std::ostream &os, const BlockMatrix &bm) {
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

double &siconos::algebra::BlockMatrix::operator()(unsigned int row, unsigned int col) {
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

double siconos::algebra::BlockMatrix::operator()(unsigned int row, unsigned int col) const {
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

double siconos::algebra::BlockMatrix::getValue(unsigned int row, unsigned int col) const {
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

void siconos::algebra::BlockMatrix::setValue(unsigned int row, unsigned int col,
                                             double value) {
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

std::shared_ptr<siconos::algebra::SiconosMatrix> siconos::algebra::BlockMatrix::block(
    unsigned int row, unsigned int col) {
  return (*_mat)(row, col);
}

std::shared_ptr<const siconos::algebra::SiconosMatrix> siconos::algebra::BlockMatrix::block(
    unsigned int row, unsigned int col) const {
  return std::shared_ptr<SiconosMatrix>((*_mat)(row, col));
}

void siconos::algebra::BlockMatrix::copyBlock(
    unsigned int i, unsigned int j, std::shared_ptr<siconos::algebra::SiconosMatrix> m) {
  (*_mat)(i, j)->resize(m->rows(), m->cols());
  *((*_mat)(i, j)) = *m;
}

void siconos::algebra::BlockMatrix::setBlock(
    unsigned int i, unsigned int j, std::shared_ptr<siconos::algebra::SiconosMatrix> m) {
  (*_mat)(i, j) = m;
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::algebra::BlockMatrix::toSiconosMatrix() const {
  // get number of blocks in a row/col of m.
  auto m = std::make_shared<SiconosMatrix>(this->rows(), this->cols());
  unsigned int posRow = 0;
  unsigned int posCol = 0;

  for (auto it : this->_mat->rowwise()) {
    for (auto it2 : it) {
      m->setBlock(posRow, posCol, *it2);
      posCol += it2->cols();
    }
    posRow += it.size();
    posCol = 0;
  }
  return m;
}

size_t siconos::algebra::BlockMatrix::nnz(double tol) {
  size_t nnz = 0;
  for (auto row : _mat->rowwise()) {
    for (auto col : row) nnz += (*col).nnz();
  }
  return nnz;
}
