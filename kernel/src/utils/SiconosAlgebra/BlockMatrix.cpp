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

#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrixException.hpp"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "SiconosAlgebra.hpp"
#include "Tools.hpp"

using  std::cout;
using std::endl;

// =================================================
//                CONSTRUCTORS
// =================================================

BlockMatrix::BlockMatrix(const SiconosMatrix &m): SiconosMatrix(0), _dimRow(0), _dimCol(0)
{
  _tabRow.reset(new Index());
  _tabCol.reset(new Index());
  if (m.isBlock())
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    unsigned int nbRows = m.numberOfBlocks(0);
    unsigned int nbCols = m.numberOfBlocks(1);
    _tabRow->reserve(nbRows);
    _tabCol->reserve(nbCols);

    // mat construction
    _mat.reset(new BlocksMat(nbRows, nbCols, nbRows * nbCols));

    unsigned int i, j;
    ConstBlocksIterator1 it1;
    ConstBlocksIterator2 it2;
    bool firstLoop = true;
    // We scan all the blocks of m ...
    for (it1 = mB._mat->begin1(); it1 != mB._mat->end1(); ++it1)
    {
      _dimRow += (*(it1.begin()))->size(0);
      _tabRow->push_back(_dimRow);
      for (it2 = it1.begin(); it2 != it1.end(); ++it2)
      {
        i = it2.index1();
        j = it2.index2();
        if ((*it2)->isBlock())  // if the current matrix is a blockMatrix
          _mat->insert_element(i, j, std11::shared_ptr<SiconosMatrix>(new BlockMatrix(**it2)));
        else
          _mat->insert_element(i, j, std11::shared_ptr<SiconosMatrix>(new SimpleMatrix(**it2)));
        // _dimCol must be incremented only at first "column-loop"
        if (firstLoop)
        {
          _dimCol += (*it2)->size(1);
          _tabCol->push_back(_dimCol);
        }
      }
      firstLoop = false;
    }
  }
  else // if m is a SimpleMatrix
  {
    _tabRow->reserve(1);
    _tabCol->reserve(1);
    // _mat construction
    _mat.reset(new BlocksMat(1, 1, 1));
    _mat->insert_element(0, 0, std11::shared_ptr<SiconosMatrix>(new SimpleMatrix(m)));

    _dimRow = m.size(0);
    _dimCol = m.size(1);
    _tabRow->push_back(_dimRow);
    _tabCol->push_back(_dimCol);
  }
}

BlockMatrix::BlockMatrix(const BlockMatrix &m): SiconosMatrix(0), _dimRow(0), _dimCol(0)
{
  unsigned int nbRows = m.numberOfBlocks(0);
  unsigned int nbCols = m.numberOfBlocks(1);
  _tabRow.reset(new Index());
  _tabCol.reset(new Index());
  _tabRow->reserve(nbRows);
  _tabCol->reserve(nbCols);

  // _mat construction
  _mat.reset(new BlocksMat(nbRows, nbCols, nbRows * nbCols));

  unsigned int i, j;
  // We scan all the blocks of m ...
  ConstBlocksIterator1 it1;
  ConstBlocksIterator2 it2;
  bool firstLoop = true;
  // We scan all the blocks of m ...
  for (it1 = m._mat->begin1(); it1 != m._mat->end1(); ++it1)
  {
    _dimRow += (*(it1.begin()))->size(0);
    _tabRow->push_back(_dimRow);
    for (it2 = it1.begin(); it2 != it1.end(); ++it2)
    {
      i = it2.index1();
      j = it2.index2();
      if ((*it2)->isBlock())  // if the current _matrix is a blockMatrix
        _mat->insert_element(i, j, std11::shared_ptr<SiconosMatrix>(new BlockMatrix(**it2)));
      else
        _mat->insert_element(i, j, std11::shared_ptr<SiconosMatrix>(new SimpleMatrix(**it2)));

      // _dimCol must be incremented only at first "column-loop"
      if (firstLoop)
      {
        _dimCol += (*it2)->size(1);
        _tabCol->push_back(_dimCol);
      }
    }
    firstLoop = false;
  }
}

BlockMatrix::BlockMatrix(const std::vector<SP::SiconosMatrix >& m, unsigned int row, unsigned int col):
  SiconosMatrix(0), _dimRow(0), _dimCol(0)
{
  if (m.size() != (row * col))
    SiconosMatrixException::selfThrow("BlockMatrix constructor from a vector<SiconosMatrix*>, number of blocks inconsistent with provided dimensions.");

  _tabRow.reset(new Index());
  _tabCol.reset(new Index());
  _tabRow->reserve(row);
  _tabCol->reserve(col);

  // _mat construction
  _mat.reset(new BlocksMat(row, col, row * col));

  unsigned int k = 0;
  bool firstRowLoop = true;
  bool firstColLoop = true;

  for (unsigned int i = 0; i < row; ++i)
  {
    for (unsigned int j = 0; j < col; ++j)
    {
      (*_mat)(i, j) = m[k++];

      // _dimCol must be incremented only at first "column-loop"
      if (firstColLoop)
      {
        _dimCol += m[k - 1]->size(1);
        _tabCol->push_back(_dimCol);
      }
      if (firstRowLoop)
      {
        _dimRow += m[k - 1]->size(0);
        _tabRow->push_back(_dimRow);
        firstRowLoop = false;
      }
    }
    firstColLoop = false;
    firstRowLoop = true;
  }
}

BlockMatrix::BlockMatrix(SP::SiconosMatrix A, SP::SiconosMatrix B, SP::SiconosMatrix C, SP::SiconosMatrix D):
  SiconosMatrix(0), _dimRow(0), _dimCol(0)
{
  if (A->size(0) != B->size(0) || C->size(0) != D->size(0) ||  A->size(1) != C->size(1) ||  B->size(1) != D->size(1))
    SiconosMatrixException::selfThrow("BlockMatrix constructor(A,B,C,D), inconsistent sizes between A, B, C or D SiconosMatrices.");

  // _mat = [ A B ]
  //       [ C D ]

  // _mat construction
  _mat.reset(new BlocksMat(2, 2, 4));

  _tabRow.reset(new Index());
  _tabCol.reset(new Index());
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

BlockMatrix::~BlockMatrix()
{

  _mat->clear();

  _tabRow->clear();
  _tabCol->clear();
}

// =================================================
//    get number of blocks
// =================================================

unsigned int BlockMatrix::numberOfBlocks(unsigned int dim) const
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
const DenseMat  BlockMatrix::getDense(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp;
  tmp = (*_mat)(row, col);

  if (tmp->num() != 1)
    SiconosMatrixException::selfThrow("DenseMat BlockMatrix::getDense(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Dense matrix");

  return (tmp->getDense());
}

// return the boost triangular matrix of the block (i, j)
const TriangMat BlockMatrix::getTriang(unsigned int row, unsigned int col) const
{
  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 2)
  {
    SiconosMatrixException::selfThrow("TriangMat BlockMatrix::getTriang(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Triangular matrix");
  }
  return (tmp->getTriang());
}

// return the boost symmetric matrix of the block (i, j)
const SymMat BlockMatrix::getSym(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 3)
  {
    SiconosMatrixException::selfThrow("SymMat BlockMatrix::getSym(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Symmmetric matrix");
  }
  return (tmp->getSym());
}

// return the boost sparse matrix of the block (i, j)
const SparseMat  BlockMatrix::getSparse(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 4)
  {
    SiconosMatrixException::selfThrow("SparseMat BlockMatrix::getSparse(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Sparse matrix");
  }
  return (tmp->getSparse());
}

// return the boost banded matrix of the block (i, j)
const BandedMat  BlockMatrix::getBanded(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 5)
  {
    SiconosMatrixException::selfThrow("BandedMat BlockMatrix::getBanded(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Banded matrix");
  }
  return (tmp->getBanded());
}

// return the boost zero matrix of the block (i, j)
const ZeroMat  BlockMatrix::getZero(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 5)
  {
    SiconosMatrixException::selfThrow("ZeroMat BlockMatrix::getZero(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Zero matrix");
  }
  return (tmp->getZero());
}

// return the boost identity matrix of the block (i, j)
const IdentityMat  BlockMatrix::getIdentity(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 5)
  {
    SiconosMatrixException::selfThrow("IdentityMat BlockMatrix::getIdentity(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Identity matrix");
  }
  return (tmp->getIdentity());
}

// The following functions return the corresponding pointers
DenseMat*  BlockMatrix::dense(unsigned int row, unsigned int col) const
{


  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 1)
  {
    SiconosMatrixException::selfThrow("DenseMat* BlockMatrix::dense(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Dense matrix");
  }

  return (tmp->dense());
}

TriangMat* BlockMatrix::triang(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 2)
  {
    SiconosMatrixException::selfThrow("TriangMat* BlockMatrix::triang(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Triangular matrix");
  }
  return (tmp->triang());
}
SymMat* BlockMatrix::sym(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 3)
  {
    SiconosMatrixException::selfThrow("SymMat* BlockMatrix::sym(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Symmmetric matrix");
  }
  return (tmp->sym());
}

SparseMat*  BlockMatrix::sparse(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 4)
  {
    SiconosMatrixException::selfThrow("SparseMat* BlockMatrix::sparse(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Sparse matrix");
  }
  return (tmp->sparse());
}

BandedMat*  BlockMatrix::banded(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 5)
  {
    SiconosMatrixException::selfThrow("BandedMat* BlockMatrix::banded(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Banded matrix");
  }
  return (tmp->banded());
}

ZeroMat*  BlockMatrix::zero_mat(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 6)
  {
    SiconosMatrixException::selfThrow("ZeroMat* BlockMatrix::zero(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Zero matrix");
  }
  return (tmp->zero_mat(row, col));
}

IdentityMat*  BlockMatrix::identity(unsigned int row, unsigned int col) const
{

  SP::SiconosMatrix tmp = (*_mat)(row, col);
  if (tmp->num() != 5)
  {
    SiconosMatrixException::selfThrow("IdentityMat* BlockMatrix::identity(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Identity matrix");
  }
  return (tmp->identity());
}

double* BlockMatrix::getArray(unsigned int i, unsigned int j) const
{
  SP::SiconosMatrix tmp = (*_mat)(i, j);
  return tmp->getArray();
}

// ===========================
//       fill matrix
// ===========================

void BlockMatrix::zero()
{
  BlocksMat::iterator1 it;
  BlocksMat::iterator2 it2;
  for (it = _mat->begin1(); it != _mat->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
    {
      (*it2)->zero();
    }
  }
}

void BlockMatrix::randomize()
{
  BlocksMat::iterator1 it;
  BlocksMat::iterator2 it2;
  for (it = _mat->begin1(); it != _mat->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
    {
      (*it2)->randomize();
    }
  }
}

void BlockMatrix::randomize_sym()
{
  SiconosMatrixException::selfThrow("BlockMatrix::randomize_sym : not yet implemented for block matrices.");
}

void BlockMatrix::eye()
{
  BlocksMat::iterator1 it;
  BlocksMat::iterator2 it2;

  for (it = _mat->begin1(); it != _mat->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
    {
      if (it2.index1() == it2.index2())
        (*it2)->eye();
      else
        (*it2)->zero();
    }
  }
}

unsigned int BlockMatrix::size(unsigned int index) const
{
  if (index == 0) return _dimRow;
  else return _dimCol;
};


//=======================
// set matrix dimensions
//=======================

void BlockMatrix::resize(unsigned int, unsigned int, unsigned int, unsigned int, bool)
{
  SiconosMatrixException::selfThrow("BlockMatrix::resize : forbidden for block matrices.");
}

//=======================
//       get norm
//=======================

double BlockMatrix::normInf()const
{
  double sum = 0, norm = 0;
  for (unsigned int i = 0; i < size(0); i++)
  {
    for (unsigned int j = 0; j < size(1); j++)
    {
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

void BlockMatrix::display(void)const
{
  std::cout << "==========> BlockMatrix (" << numberOfBlocks(0) << " X " << numberOfBlocks(1) << " blocks): " << std::endl;
  BlocksMat::iterator1 it;
  BlocksMat::iterator2 it2;
  for (it = _mat->begin1(); it != _mat->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
    {
      (*it2)->display();
    }
  }
  std::cout << "===========================================================================================" << std::endl;
}

//=====================
// convert to a string
//=====================

std::string BlockMatrix::toString() const
{
  return ::toString(*this);
}

//=====================
// convert to an ostream
//=====================

std::ostream& operator<<(std::ostream& os, const BlockMatrix& bm)
{
  BlocksMat::iterator1 it;
  BlocksMat::iterator2 it2;
  os << "[" << bm.numberOfBlocks(0) << "," << bm.numberOfBlocks(1) << "](";
  for (it = bm._mat->begin1(); it != bm._mat->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
    {
      if (it2 != it.begin()) os << ",";
      if (*it2) os << **it2; else os << "(nil)";
    }
  }
  os << ")";
  return os;
}

//=============================
// Elements access (get or set)
//=============================

double& BlockMatrix::operator()(unsigned int row, unsigned int col)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size())
    nbRow ++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0)
    posCol -= (*_tabCol)[nbCol - 1];


  SP::SiconosMatrix tmp = (*_mat)(nbRow, nbCol);
  return (*tmp)(posRow, posCol);
}

double BlockMatrix::operator()(unsigned int row, unsigned int col) const
{

  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size())
    nbRow ++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0)
    posCol -= (*_tabCol)[nbCol - 1];

  SP::SiconosMatrix tmp = (*_mat)(nbRow, nbCol);
  return (*tmp)(posRow, posCol);
}

double BlockMatrix::getValue(unsigned int row, unsigned int col) const
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size())
    nbRow ++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0)
    posCol -= (*_tabCol)[nbCol - 1];


  SP::SiconosMatrix tmp = (*_mat)(nbRow, nbCol);
  return (*tmp)(posRow, posCol);
}

void BlockMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= (*_tabRow)[nbRow] && nbRow < _tabRow->size())
    nbRow ++;

  while (col >= (*_tabCol)[nbCol] && nbCol < _tabCol->size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= (*_tabRow)[nbRow - 1];
  if (nbCol != 0)
    posCol -= (*_tabCol)[nbCol - 1];

  SP::SiconosMatrix tmp = (*_mat)(nbRow, nbCol);
  (*tmp)(posRow, posCol) = value;
}

//============================================
// Access (get or set) to blocks of elements
//============================================

// void BlockMatrix::getBlock (unsigned int row, unsigned int col, SiconosMatrix * m) const
// {
//   SiconosMatrixException::selfThrow("BlockMatrix::getBlock, not yet implemented or useless for BlockMatrices.");
// }

// void BlockMatrix::setBlock(unsigned int row, unsigned int col, const SiconosMatrix *m)
// {
//   // Set current matrix elements, starting from row row_min and column col_min, with the values of the matrix m.
//   // m may be a BlockMatrix.

//   if(m == this)
//     SiconosMatrixException::selfThrow("BlockMatrix::setBlock(pos,..., m): m = this.");

//   if(m->isBlock ())
//     SiconosMatrixException::selfThrow("BlockMatrix::setBlock of a block into an other block is forbidden.");

//   if(row > _dimRow || col > _dimCol )
//     SiconosMatrixException::selfThrow("BlockMatrix::setBlock(i,j,m), i or j is out of range.");

//   // Check dim
//   if( tmp->size(0)!=m->size(0) || tmp->size(1) != m->size(1) )
//     SiconosMatrixException::selfThrow("BlockMatrix::setBlock(x,y,m), block(x,y) of current matrix and m have inconsistent sizes.");
//   *((*_mat)(row,col)) = *m; // copy
// }

void BlockMatrix::getRow(unsigned int r, SiconosVector &v) const
{
  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (r > _dimRow)
    SiconosMatrixException::selfThrow("BlockMatrix:getRow : row number is out of range");

  // Verification of the size of the result vector
  if (v.size() != _dimCol)
    SiconosMatrixException::selfThrow("BlockMatrix:getRow : inconsistent sizes");

  // Find the row-block number where "r" is
  while (r >= (*_tabRow)[numRow] && numRow < _tabRow->size())
    numRow ++;

  // Computation of the value of the index row into this block
  if (numRow != 0)
    posRow -= (*_tabRow)[numRow - 1];

  for (unsigned int j = 0; j < _tabCol->size(); j++)
  {
    start = stop;
    SP::SiconosMatrix tmp = (*_mat)(numRow, j);
    stop += tmp->size(1);
    ublas::subrange(*(v.dense()), start, stop) = ublas::row(*(tmp->dense()), posRow);
  }
}

void BlockMatrix::getCol(unsigned int c, SiconosVector &v) const
{
  unsigned int numCol = 0, posCol = c, start = 0, stop = 0;

  if (c > _dimCol)
    SiconosMatrixException::selfThrow("BlockMatrix:getCol : column number is out of range");

  // Verification of the size of the result vector
  if (v.size() != _dimRow)
    SiconosMatrixException::selfThrow("BlockMatrix:getcol : inconsistent sizes");

  // Find the column-block number where "c" is
  while (c >= (*_tabCol)[numCol] && numCol < _tabCol->size())
    numCol ++;

  // Computation of the value of the index column into this block
  if (numCol != 0)
    posCol -= (*_tabCol)[numCol - 1];

  for (unsigned int i = 0; i < _tabRow->size(); i++)
  {
    start = stop;
    SP::SiconosMatrix tmp = (*_mat)(i, numCol);
    stop += tmp->size(0);
    ublas::subrange(*(v.dense()), start, stop) = ublas::column(tmp->getDense(), posCol);
  }
}

void BlockMatrix::setRow(unsigned int r, const SiconosVector &v)
{

  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (v.size() != _dimCol)
    SiconosMatrixException::selfThrow("BlockMatrix:setRow : inconsistent sizes");

  while (r >= (*_tabRow)[numRow] && numRow < _tabRow->size())
    numRow ++;

  if (numRow != 0)
    posRow -= (*_tabRow)[numRow - 1];

  for (unsigned int j = 0; j < _tabCol->size(); j++)
  {
    start = stop;
    SP::SiconosMatrix tmp = (*_mat)(numRow, j);
    stop += tmp->size(1);
    ublas::row(*(tmp->dense()), posRow) = ublas::subrange(*(v.dense()), start, stop);
  }
}

void BlockMatrix::setCol(unsigned int col, const SiconosVector &v)
{

  unsigned int numCol = 0, posCol = col, start = 0, stop = 0;

  if (v.size() != _dimRow)
    SiconosMatrixException::selfThrow("BlockMatrix:setCol : inconsistent sizes");

  while (col >= (*_tabCol)[numCol] && numCol < _tabCol->size())
    numCol ++;

  if (numCol != 0)
    posCol -= (*_tabCol)[numCol - 1];

  for (unsigned int i = 0; i < _tabRow->size(); i++)
  {
    start = stop;
    SP::SiconosMatrix tmp = (*_mat)(i, numCol);
    stop += tmp->size(0);
    ublas::column(*(tmp->dense()), posCol) = ublas::subrange(*(v.dense()), start, stop);
  }
}

void BlockMatrix::addSimple(unsigned int& indRow, unsigned int& indCol, const SiconosMatrix& m)
{
  // Add a part of m (starting from (indRow,indCol) to the current matrix.
  // m must be a SimpleMatrix.

  // At the end of the present function, indRow (resp. indCol) is equal to indRow + the corresponding dimension of the added sub-matrix.

  unsigned int row = m.size(0) - indRow; // number of rows of the block to be added.
  unsigned int col = m.size(1) - indCol; // number of columns of the block to be added.
  unsigned int initCol = indCol;

  if (row > _dimRow || col > _dimCol) SiconosMatrixException::selfThrow("BlockMatrix::addSimple : invalid ranges");

  unsigned int numM = m.num();

  // iterators through this
  BlocksMat::iterator1 it1;
  BlocksMat::iterator2 it2;
  unsigned int currentRow = 0 , currentCol = 0, currentNum;
  for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
  {
    for (it2 = it1.begin(); it2 != it1.end(); ++it2)
    {
      if ((*it2)->isBlock())  // if the sub-block is also a BlockMatrix ...
        (std11::static_pointer_cast<BlockMatrix>(*it2))->addSimple(indRow, indCol, m);

      else
      {
        currentCol = (*it2)->size(1);
        currentRow = (*it2)->size(0);
        currentNum = (*it2)->num();
        if (numM != currentNum) SiconosMatrixException::selfThrow("BlockMatrix::addSimple : inconsistent types.");

        if (numM == 1)
          noalias(*(*it2)->dense()) += ublas::subrange(*m.dense(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 2)
          noalias(*(*it2)->triang()) += ublas::subrange(*m.triang(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 3)
          noalias(*(*it2)->sym()) += ublas::subrange(*m.sym(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 4)
          noalias(*(*it2)->sparse()) += ublas::subrange(*m.sparse(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 5)
          noalias(*(*it2)->banded()) += ublas::subrange(*m.banded(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 6) {}
        else
          SiconosMatrixException::selfThrow("BlockMatrix::addSimple : inconsistent types.");
      }
      indCol += currentCol;
    }
    indRow += currentRow;
    indCol = initCol;
  }
}

void BlockMatrix::subSimple(unsigned int& indRow, unsigned int& indCol, const SiconosMatrix& m)
{
  // subtract a part of m (starting from (indRow,indCol) to the current matrix.
  // m must be a SimpleMatrix.

  // At the end of the present function, indRow (resp. indCol) is equal to indRow + the corresponding dimension of the subtracted sub-matrix.

  unsigned int row = m.size(0) - indRow; // number of rows of the block to be added.
  unsigned int col = m.size(1) - indCol; // number of columns of the block to be added.
  unsigned int initCol = indCol;
  if (row > _dimRow || col > _dimCol) SiconosMatrixException::selfThrow("BlockMatrix::addSimple : invalid ranges");

  unsigned int numM = m.num();

  // iterators through this
  BlocksMat::iterator1 it1;
  BlocksMat::iterator2 it2;
  unsigned int currentRow = 0, currentCol = 0, currentNum;
  for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
  {
    for (it2 = it1.begin(); it2 != it1.end(); ++it2)
    {
      if ((*it2)->isBlock())  // if the sub-block is also a BlockMatrix ...
        (std11::static_pointer_cast<BlockMatrix>(*it2))->subSimple(indRow, indCol, m);

      else
      {
        currentCol = (*it2)->size(1);
        currentRow = (*it2)->size(0);
        currentNum = (*it2)->num();
        if (numM != currentNum) SiconosMatrixException::selfThrow("BlockMatrix::addSimple : inconsistent types.");

        if (numM == 1)
          noalias(*(*it2)->dense()) -= ublas::subrange(*m.dense(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 2)
          noalias(*(*it2)->triang()) -= ublas::subrange(*m.triang(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 3)
          noalias(*(*it2)->sym()) -= ublas::subrange(*m.sym(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 4)
          noalias(*(*it2)->sparse()) -= ublas::subrange(*m.sparse(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 5)
          noalias(*(*it2)->banded()) -= ublas::subrange(*m.banded(), indRow, indRow + currentRow, indCol, indCol + currentCol);
        else if (numM == 6) {}
        else
          SiconosMatrixException::selfThrow("BlockMatrix::addSimple : inconsistent types.");
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

BlockMatrix& BlockMatrix::operator = (const SiconosMatrix &m)
{
  if (&m == this) return *this; // auto-assignment.

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    SiconosMatrixException::selfThrow("operator = (const SiconosMatrix&): Left and Right values have inconsistent sizes.");

  // Warning: we do not reallocate the blocks, but only copy the values. This means that
  // all blocks are already allocated and that dim of m and mat are to be consistent.
  // Thus, _tabRow and _tabCol remains unchanged.
  // If m and mat are not "block-consistent", we use the () operator for a component-wise copy.

  if (m.isBlock())
  {
    if (isComparableTo(*this, m))
    {
      const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
      // iterators through this
      BlocksIterator1 it1;
      BlocksIterator2 it2;
      // iterators through m
      ConstBlocksIterator1 itM1 = mB._mat->begin1();
      ConstBlocksIterator2 itM2;

      for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
      {
        itM2 = itM1.begin();
        for (it2 = it1.begin(); it2 != it1.end(); ++it2)
        {
          (**it2) = (**itM2);
          itM2++; // increment column pos. in m.
        }
        itM1++; // increment row pos. in m.
      }
    }
    else
    {
      for (unsigned int i = 0; i < _dimRow; ++i)
        for (unsigned int j = 0; j < _dimCol; ++j)
          (*this)(i, j) = m(i, j);
    }
  }
  else // if m is a SimpleMatrix
  {
    BlocksIterator1 it;
    BlocksIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;
    Index subDim(2);
    Index subPos(4);

    for (it = _mat->begin1(); it != _mat->end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        // a sub-block of m is copied into this
        subDim[0] = (*it2)->size(0);
        subDim[1] = (*it2)->size(1);
        subPos[0] = posRow;
        subPos[1] = posCol;
        subPos[2] = 0;
        subPos[3] = 0;
        setBlock(createSPtrConstSiconosMatrix(m), *it2, subDim, subPos);
        posCol += subDim[1];
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }

  return *this;
}

BlockMatrix& BlockMatrix::operator = (const BlockMatrix &m)
{
  if (&m == this) return *this; // auto-assignment.

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    SiconosMatrixException::selfThrow("operator = (const SiconosMatrix&): Left and Right values have inconsistent sizes.");

  // Warning: we do not reallocate the blocks, but only copy the values. This means that
  // all blocks are already allocated and that dim of m and mat are to be consistent.
  // Thus, _tabRow and _tabCol remains unchanged.
  // If m and mat are not "block-consistent", we use the () operator for a componet-wise copy.

  if (isComparableTo(*this, m))
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    // iterators through this
    BlocksIterator1 it1;
    BlocksIterator2 it2;
    // iterators through m
    ConstBlocksIterator1 itM1 = mB._mat->begin1();
    ConstBlocksIterator2 itM2;

    for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
    {
      itM2 = itM1.begin();
      for (it2 = it1.begin(); it2 != it1.end(); ++it2)
      {
        (**it2) = (**itM2);
        itM2++; // increment column pos. in m.
      }
      itM1++; // increment row pos. in m.
    }
  }
  else
  {
    for (unsigned int i = 0; i < _dimRow; ++i)
      for (unsigned int j = 0; j < _dimCol; ++j)
        (*this)(i, j) = m(i, j);
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator = (const DenseMat &m)
{
  SiconosMatrixException::selfThrow("BlockMatrix:operator = DenseMat - Not yet implemented.");
  return *this;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

BlockMatrix& BlockMatrix::operator += (const SiconosMatrix &m)
{
  if (&m == this)
  {
    BlocksMat::iterator1 it1;
    BlocksMat::iterator2 it2;
    for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
    {
      for (it2 = it1.begin(); it2 != it1.end(); ++it2)
      {
        **it2 += **it2;
      }
    }
    return *this;
  }

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    SiconosMatrixException::selfThrow("BlockMatrix::operator += Left and Right values have inconsistent sizes.");

  if (m.isBlock())
  {
    if (isComparableTo(m, *this))
    {
      const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
      // iterators through this
      BlocksMat::iterator1 it1;
      BlocksMat::iterator2 it2;
      // iterators through m
      BlocksMat::const_iterator1 itM1 = mB._mat->begin1();
      BlocksMat::const_iterator2 itM2;

      for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
      {
        itM2 = itM1.begin();
        for (it2 = it1.begin(); it2 != it1.end(); ++it2)
        {
          (**it2) += (**itM2);
          itM2++; // increment column pos. in m.
        }
        itM1++; // increment row pos. in m.
      }
    }
    else
    {
      for (unsigned int i = 0; i < _dimRow; ++i)
        for (unsigned int j = 0; j < _dimCol; ++j)
          (*this)(i, j) += m(i, j);
    }
  }
  else // if m is a SimpleMatrix
  {
    unsigned int indRow = 0, indCol = 0;
    addSimple(indRow, indCol, m); // a sub-block of m is added to each block of this.
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator -= (const SiconosMatrix &m)
{
  if (&m == this)
  {
    BlocksMat::iterator1 it1;
    BlocksMat::iterator2 it2;
    for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
    {
      for (it2 = it1.begin(); it2 != it1.end(); ++it2)
      {
        **it2 -= **it2;
      }
    }
    return *this;
  }

  if (m.size(0) != _dimRow || m.size(1) != _dimCol)
    SiconosMatrixException::selfThrow("BlockMatrix::operator += Left and Right values have inconsistent sizes.");

  if (m.isBlock())
  {
    if (isComparableTo(m, *this))
    {
      const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
      // iterators through this
      BlocksMat::iterator1 it1;
      BlocksMat::iterator2 it2;
      // iterators through m
      BlocksMat::const_iterator1 itM1 = mB._mat->begin1();
      BlocksMat::const_iterator2 itM2;

      for (it1 = _mat->begin1(); it1 != _mat->end1(); ++it1)
      {
        itM2 = itM1.begin();
        for (it2 = it1.begin(); it2 != it1.end(); ++it2)
        {
          (**it2) -= (**itM2);
          itM2++; // increment column pos. in m.
        }
        itM1++; // increment row pos. in m.
      }
    }
    else
    {
      for (unsigned int i = 0; i < _dimRow; ++i)
        for (unsigned int j = 0; j < _dimCol; ++j)
          (*this)(i, j) -= m(i, j);
    }
  }
  else // if m is a SimpleMatrix
  {
    unsigned int indRow = 0, indCol = 0;
    subSimple(indRow, indCol, m); // a sub-block of m is subtracted to each block of this.
  }
  return *this;
}

void BlockMatrix::trans()
{
  SiconosMatrixException::selfThrow("BlockMatrix::trans(): not yet implemented.");
}

void BlockMatrix::trans(const SiconosMatrix &m)
{
  SiconosMatrixException::selfThrow("BlockMatrix::trans(M): not yet implemented.");
}

void BlockMatrix::PLUFactorizationInPlace()
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUFactorizationInPlace: not yet implemented for Block Matrices.");
}

void BlockMatrix::PLUInverseInPlace()
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUInverseInPlace: not yet implemented for Block Matrices.");
}

void BlockMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUForwardBackwardInPlace: not yet implemented for Block Matrices.");
}

void BlockMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUForwardBackwardInPlace: not yet implemented for Block Matrices.");
}

SP::SiconosMatrix BlockMatrix::block(unsigned int row, unsigned int col)
{
  return (*_mat)(row, col);
}

SPC::SiconosMatrix BlockMatrix::block(unsigned int row, unsigned int col) const
{
  return std11::shared_ptr<SiconosMatrix>((*_mat)(row, col));
}

size_t BlockMatrix::nnz(double tol)
{
  size_t nnz = 0;
  BlocksMat::iterator1 it;
  BlocksMat::iterator2 it2;
  for (it = _mat->begin1(); it != _mat->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
      nnz += (**it2).nnz();
  }
  return nnz;
}

