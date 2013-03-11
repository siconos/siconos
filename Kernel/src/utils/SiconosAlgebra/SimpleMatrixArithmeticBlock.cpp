/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(pos,..., m) forbidden for zero or identity matrix.");

  if (&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(pos,..., m): m = this.");

  if (row_min >= dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(row,col): row is out of range");

  if (col_min >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBloc(row,col)k: col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if (row_max > dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(row,col,m): m.row + row is out of range.");

  if (col_max > dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(row,col,m): m.col + col is out of range.");

  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block matrix ...
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for (it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        addBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if (numM == 6) // if m = 0
  {
    // nothing to do !
  }
  else // if m is a SimpleMatrix
  {
    if (num == 1)
    {
      switch (numM)
      {
      case 1:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.dense());
        break;
      case 2:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.triang());
        break;
      case 3:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.sym());
        break;
      case 4:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.sparse());
        break;
      case 5:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.banded());
        break;
      case 7:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.identity());
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(...,m): wrong matrix type for m.");
        break;
      }
    }
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(...): implemeted only for dense matrices.");
    resetLU();
  }
}

void SimpleMatrix::subBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // sub m to current matrix elements, starting from row row_min and column col_min, to the values of the matrix m.
  // m may be a BlockMatrix.

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(pos,..., m) forbidden for zero or identity matrix.");

  if (&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(pos,..., m): m = this.");

  if (row_min >= dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col): row is out of range");

  if (col_min >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col): col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if (row_max > dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col,m): m.row + row is out of range.");

  if (col_max > dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col,m): m.col + col is out of range.");

  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block matrix ...
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for (it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        subBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if (numM == 6) // if m = 0
  {
    // nothing to do !
  }
  else // if m is a SimpleMatrix
  {
    if (num == 1)
    {
      switch (numM)
      {
      case 1:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.dense());
        break;
      case 2:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.triang());
        break;
      case 3:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.sym());
        break;
      case 4:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.sparse());
        break;
      case 5:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.banded());
        break;
      case 7:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.identity());
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(...,m): wrong matrix type for m.");
        break;
      }
    }
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(...): implemeted only for dense matrices.");
    resetLU();
  }
}


