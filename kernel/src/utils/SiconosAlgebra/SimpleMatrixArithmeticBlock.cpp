/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "SimpleMatrix.hpp"
#include "BlockMatrix.hpp"
#include "SiconosAlgebra.hpp"

using namespace Siconos;

void SimpleMatrix::addBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // add m to current matrix elements, starting from row row_min and column col_min, to the values of the matrix m.
  // m may be a BlockMatrix.

  if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    THROW_EXCEPTION("SimpleMatrix::addBlock(pos,..., m) forbidden for zero or identity matrix.");

  if(&m == this)
    THROW_EXCEPTION("SimpleMatrix::addBlock(pos,..., m): m = this.");

  if(row_min >= size(0))
    THROW_EXCEPTION("SimpleMatrix::addBlock(row,col): row is out of range");

  if(col_min >= size(1))
    THROW_EXCEPTION("SimpleMatrix::addBloc(row,col)k: col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if(row_max > size(0))
    THROW_EXCEPTION("SimpleMatrix::addBlock(row,col,m): m.row + row is out of range.");

  if(col_max > size(1))
    THROW_EXCEPTION("SimpleMatrix::addBlock(row,col,m): m.col + col is out of range.");

  Siconos::UBLAS_TYPE numM = m.num();

  if(numM == Siconos::BLOCK)  // if m is a block matrix ...
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for(it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for(it2 = it.begin(); it2 != it.end(); ++it2)
      {
        addBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if(numM == Siconos::ZERO)  // if m = 0
  {
    // nothing to do !
  }
  else // if m is a SimpleMatrix
  {
    if(_num == Siconos::DENSE)
    {
      switch(numM)
      {
      case Siconos::DENSE:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.dense());
        break;
      case Siconos::TRIANGULAR:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.triang());
        break;
      case Siconos::SYMMETRIC:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.sym());
        break;
      case Siconos::SPARSE:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.sparse());
        break;
      case Siconos::BANDED:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.banded());
        break;
      case Siconos::IDENTITY:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.identity());
        break;
      default:
        THROW_EXCEPTION("SimpleMatrix::addBlock(...,m): wrong matrix type for m.");
        break;
      }
    }
    else
      THROW_EXCEPTION("SimpleMatrix::addBlock(...): implemented only for dense matrices.");
    resetFactorizationFlags();
  }
}

void SimpleMatrix::subBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // sub m to current matrix elements, starting from row row_min and column col_min, to the values of the matrix m.
  // m may be a BlockMatrix.

  if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    THROW_EXCEPTION("SimpleMatrix::subBlock(pos,..., m) forbidden for zero or identity matrix.");

  if(&m == this)
    THROW_EXCEPTION("SimpleMatrix::subBlock(pos,..., m): m = this.");

  if(row_min >= size(0))
    THROW_EXCEPTION("SimpleMatrix::subBlock(row,col): row is out of range");

  if(col_min >= size(1))
    THROW_EXCEPTION("SimpleMatrix::subBlock(row,col): col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if(row_max > size(0))
    THROW_EXCEPTION("SimpleMatrix::subBlock(row,col,m): m.row + row is out of range.");

  if(col_max > size(1))
    THROW_EXCEPTION("SimpleMatrix::subBlock(row,col,m): m.col + col is out of range.");

  Siconos::UBLAS_TYPE numM = m.num();

  if(numM == 0)  // if m is a block matrix ...
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for(it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for(it2 = it.begin(); it2 != it.end(); ++it2)
      {
        subBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if(numM == Siconos::ZERO)  // if m = 0
  {
    // nothing to do !
  }
  else // if m is a SimpleMatrix
  {
    if(_num == Siconos::DENSE)
    {
      switch(numM)
      {
      case Siconos::DENSE:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.dense());
        break;
      case Siconos::TRIANGULAR:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.triang());
        break;
      case Siconos::SYMMETRIC:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.sym());
        break;
      case Siconos::SPARSE:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.sparse());
        break;
      case Siconos::BANDED:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.banded());
        break;
      case Siconos::IDENTITY:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.identity());
        break;
      default:
        THROW_EXCEPTION("SimpleMatrix::subBlock(...,m): wrong matrix type for m.");
        break;
      }
    }
    else
      THROW_EXCEPTION("SimpleMatrix::subBlock(...): implemented only for dense matrices.");
    resetFactorizationFlags();
  }
}


