/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"

using namespace Siconos;

//=============================
// Elements access (get or set)
//=============================

double SimpleMatrix::getValue(unsigned int row, unsigned int col) const
{
  if(row >= size(0) || col >= size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix:getValue(index): Index out of range");

  if(_num == Siconos::DENSE)
    return (*mat.Dense)(row, col);
  else if(_num == Siconos::TRIANGULAR)
    return (*mat.Triang)(row, col);
  else if(_num == Siconos::SYMMETRIC)
    return (*mat.Sym)(row, col);
  else if(_num == Siconos::SPARSE)
  {
    double * d = (*mat.Sparse).find_element(row, col);
    if(d)
      return *d;
    else
      return 0.0;
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    double * d = (*mat.SparseCoordinate).find_element(row, col);
    if(d)
      return *d;
    else
      return 0.0;
  }
  else if(_num == Siconos::BANDED)
    return (*mat.Banded)(row, col);
  else if(_num == Siconos::ZERO)
    return 0;
  else //if (_num == Siconos::IDENTITY)
    return(row == col);
}

void SimpleMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  if(row >= size(0) || col >= size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix:setValue: Index out of range");

  if(_num == Siconos::DENSE)
    (*mat.Dense)(row, col) = value;
  else if(_num == Siconos::TRIANGULAR)
    (*mat.Triang)(row, col) = value;
  else if(_num == Siconos::SYMMETRIC)
    (*mat.Sym)(row, col) = value ;
  else if(_num == Siconos::SPARSE)
  {
    double * d = (*mat.Sparse).find_element(row, col);
    if(d)
    {
      *d = value;
    }
    else
    {
      (*mat.Sparse).insert_element(row, col, value);
    }
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    // double * d = (*mat.Sparse).find_element(row, col);
    // if (d)
    // {
    //   *d = value;
    // }
    // else
    // {
    (*mat.SparseCoordinate).insert_element(row, col, value);
    // }
  }

  else if(_num == Siconos::BANDED)
    (*mat.Banded)(row, col) = value;
  else if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("SimpleMatrix:setValue: forbidden for Identity or Zero type matrices.");
  resetLU();

}

//============================================
// Access (get or set) to blocks of elements
//============================================

void SimpleMatrix::setBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // Set current matrix elements, starting from row row_min and column col_min, with the values of the matrix m.
  // m may be a BlockMatrix.

  // Exceptions ...
  if(&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(pos,..., m): m = this.");

  if(row_min >= size(0))
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col): row is out of range");

  if(col_min >= size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col): col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if(row_max > size(0))
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col,m): m.row + row is out of range.");

  if(col_max > size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col,m): m.col + col is out of range.");

  unsigned int numM = m.num();

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
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = col_min;
    }
  }
  else // if m is a SimpleMatrix
  {
    if(numM != _num)
      SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(i,j,m), inconsistent types.");

    if(_num == Siconos::DENSE)
      noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) = *(m.dense());
    else if(_num == Siconos::TRIANGULAR)
      noalias(ublas::subrange(*mat.Triang, row_min, row_max, col_min, col_max)) = *(m.triang());
    else if(_num == Siconos::SYMMETRIC)
      noalias(ublas::subrange(*mat.Sym, row_min, row_max, col_min, col_max)) = *(m.sym());
    else if(_num == Siconos::SPARSE)
      noalias(ublas::subrange(*mat.Sparse, row_min, row_max, col_min, col_max)) = *(m.sparse());
    else if(_num == Siconos::BANDED)
      noalias(ublas::subrange(*mat.Banded, row_min, row_max, col_min, col_max)) = *(m.banded());
    else // if(_num == Siconos::ZERO) or _num == Siconos::IDENTITY nothing to do
    {}
    resetLU();
  }
}

void setBlock(SPC::SiconosMatrix  MIn, SP::SiconosMatrix MOut, const Index& dim, const Index& start)
{
  // To copy a subBlock of MIn into a subBlock of MOut.
  // dim[0], dim[1]: number of rows and columns of the sub-block
  // start[0], start[1]: position (row, column) of the first element of the subBlock in MIn
  // start[2], start[3]: position (row, column) of the first element of the subBlock in MOut

  if(MIn == MOut)  // useless op => nothing to be done.
  {}// SiconosVectorException::selfThrow("");
  else
  {
    unsigned int numIn = MIn->num();
    unsigned int numOut = MOut->num();

    if(numOut == 6 || numOut == 7)  // if MOut = 0 or Identity => read-only
      SiconosMatrixException::selfThrow("matrices, setBlock(MIn, MOut...): MOut is read-only (zero or identity matrix?).");

    // Check dimension
    Index MDim(4); // dim. of matrices MIn and MOut.
    MDim[0] = MIn->size(0);
    MDim[1] = MIn->size(1);
    MDim[2] = MOut->size(0);
    MDim[3] = MOut->size(1);

    for(unsigned int i = 0; i < 4 ; ++i)
      if(start[i] >= MDim[i])
        SiconosMatrixException::selfThrow("matrices, setBlock(MIn, ...): sub-block indices are out of range.");

    // index position of the last element in subBlock ...
    Index end(4);
    end[0] = dim[0] + start[0];
    end[1] = dim[1] + start[1];
    end[2] = dim[0] + start[2];
    end[3] = dim[1] + start[3];

    for(unsigned int i = 0; i < 4 ; ++i)
      if(end[i] > MDim[i])
        SiconosMatrixException::selfThrow("matrices, setBlock(MIn, ...): sub-block indices are out of range.");

    // Elements from row/col start[i] to row/col (end[i]-1) will be copied.

    // If both matrices MIn and MOut are block, exception.
    if(numIn == 0 && numOut == 0)
      SiconosMatrixException::selfThrow("matrices, setBlock(MIn, MOut ...): not yet implemented for MIn and MOut both BlockMatrix. Try to use setBlock on the sub-matrices?");

    else if(numOut == 0)  // if MOut is a BlockMatrix.
    {

      // Steps:
      // A - Find the blocks of MOut that "own" indices start[2] and end[2] ie
      //     the first and last sub-block to be set in a block-column
      //         --> numbers blockStart0 and blockEnd0
      // B - Find the  Block of MOut that "owns" index start[3] and end[3] ie
      //     the first sub-block to be set in a block-row
      //         --> numbers blockStart1 and blockEnd1
      //
      //        => The blocks concerned in MOut, are those between (block) rows blockStart0 and blockEnd0
      //           and (block) columns blockStart1 and blockEnd1.
      //
      // C - Loop through the concerned blocks (name = currentBlock) of MOut and call setBlock(MIn, currentBlock, subSize, currentPos).
      //     subSize: dim. of the considered sub-block of currentBlock to be set
      //     currentPos: same as "start" vector but for currentBlock
      //

      // A - Block-Row position: we look for the block of MOut that include index start[2] and end[2].
      //
      unsigned int blockStart0 = 0;
      SPC::Index tab = MOut->tabRow();
      while(start[2] >= (*tab)[blockStart0] && blockStart0 < tab->size())
        blockStart0 ++;
      // Relative position in the block blockStart0 of the first element to be set.
      unsigned int posOut0 = start[2];
      if(blockStart0 != 0)
        posOut0 -= (*tab)[blockStart0 - 1];

      unsigned int blockEnd0 = blockStart0;
      while(end[2] > (*tab)[blockEnd0] && blockEnd0 < tab->size())
        blockEnd0 ++;

      // Size of the last sub-block in the column of block
      unsigned int lastBlockSize0 = end[2];
      if(blockEnd0 != 0)
        lastBlockSize0 -= (*tab)[blockEnd0 - 1];

      // B - Block-Col position: we look for the block of MOut that include index start[3] and end[3].
      unsigned int blockStart1 = 0;
      tab = MOut->tabCol();
      while(start[3] >= (*tab)[blockStart1] && blockStart1 < tab->size())
        blockStart1 ++;
      // Relative position in the block blockStart1 of the first element to be set.
      unsigned int posOut1 = start[3];
      if(blockStart1 != 0)
        posOut1 -= (*tab)[blockStart1 - 1];

      unsigned int blockEnd1 = blockStart1;
      while(end[3] > (*tab)[blockEnd1] && blockEnd1 < tab->size())
        blockEnd1 ++;

      // Size of the last sub-block in the row of block
      unsigned int lastBlockSize1 = end[3];
      if(blockEnd1 != 0)
        lastBlockSize1 -= (*tab)[blockEnd1 - 1];

      //C - Next, 3 steps for each row:
      // - set first sub-block in the row (number blockStart1)
      // - set all other blocks in the row except the last one
      // - set last block (number blockEnd1)
      // Same process for other rows ...

      // The current considered block
      SP::SiconosMatrix   currentBlock = MOut->block(blockStart0, blockStart1);

      // dim of the subBlock of currentBlock to be set.
      Index subSize(2);
      // indices of the first element of MIn (resp. currentBlock) to be read (resp. set)  (same as start for MIn and MOut).
      Index currentPos(4);

      // currentBlock position in MOut.
      unsigned int numRow = blockStart0;
      unsigned int numCol = blockStart1;

      // Init currentPos
      // row and col position for first element to be read in MIn,
      currentPos[0] = start[0];
      currentPos[1] = start[1];
      // row and col position for first element in sub-block of Mout (namely currentBlock).
      currentPos[2] = posOut0;
      currentPos[3] = posOut1;

      while(numRow != blockEnd0 + 1)
      {

        while(numCol != blockEnd1 + 1)
        {
          // Get the block of MOut from which a sub-block will be set ...
          currentBlock = MOut->block(numRow, numCol);

          // Set subSize[0], dim (rows) and subSize[1], dim (columns) of the sub-block.
          // subSize[0] is only required for the first block in the row, after it remains constant.
          subSize[1] = currentBlock->size(1);

          // Warning: test "a" must be done before test "b"
          if(numCol == blockEnd1)  // if last column of blocks -> test "a"
            subSize[1] = lastBlockSize1;

          if(numCol == blockStart1)  // -> test "b"
          {
            subSize[1] -= posOut1;
            subSize[0] = currentBlock->size(0);
            if(numRow == blockEnd0)  // if last row of blocks
              subSize[0] = lastBlockSize0;
            if(numRow == blockStart0)  // if first row of blocks
              subSize[0] -= posOut0;
          }

          // Set sub-block
          setBlock(MIn, currentBlock, subSize, currentPos);

          // Update currentPos:
          // col position for first element to be read in MIn,
          currentPos[1] += subSize[1] ;
          // col position for first element to be set in sub-block.
          currentPos[3] = 0;
          numCol++;
        }

        numCol = blockStart1;
        numRow++;

        // Update currentPos:
        // row position for first element to be read in MIn,
        currentPos[0] += subSize[0] ;
        // col position for first element to be read in MIn,
        currentPos[1] = start[1] ;
        // row position for first element to be set in sub-block.
        currentPos[2] = 0;
        // col position for first element to be set in sub-block.
        currentPos[3] = posOut1;

      }

    }
    else if(numIn == 0)  // If MIn is a BlockMatrix.
    {

      // Same process as for numOut == 0

      unsigned int blockStart0 = 0;
      SPC::Index tab = MIn->tabRow();
      while(start[0] >= (*tab)[blockStart0] && blockStart0 < tab->size())
        blockStart0 ++;
      // Relative position in the block blockStart0 of the first element to be set.
      unsigned int posOut0 = start[0];
      if(blockStart0 != 0)
        posOut0 -= (*tab)[blockStart0 - 1];

      unsigned int blockEnd0 = blockStart0;
      while(end[0] > (*tab)[blockEnd0] && blockEnd0 < tab->size())
        blockEnd0 ++;

      // Size of the last sub-block in the column of block
      unsigned int lastBlockSize0 = end[0];
      if(blockEnd0 != 0)
        lastBlockSize0 -= (*tab)[blockEnd0 - 1];

      // B - Block-Col position: we look for the block of MOut that include index start[3] and end[3].
      unsigned int blockStart1 = 0;
      tab = MIn->tabCol();
      while(start[1] >= (*tab)[blockStart1] && blockStart1 < tab->size())
        blockStart1 ++;
      // Relative position in the block blockStart1 of the first element to be set.
      unsigned int posOut1 = start[1];
      if(blockStart1 != 0)
        posOut1 -= (*tab)[blockStart1 - 1];

      unsigned int blockEnd1 = blockStart1;
      while(end[1] > (*tab)[blockEnd1] && blockEnd1 < tab->size())
        blockEnd1 ++;

      // Size of the last sub-block in the row of block
      unsigned int lastBlockSize1 = end[1];
      if(blockEnd1 != 0)
        lastBlockSize1 -= (*tab)[blockEnd1 - 1];

      //C - Next, 3 steps for each row:
      // - set first sub-block in the row (number blockStart1)
      // - set all other blocks in the row except the last one
      // - set last block (number blockEnd1)
      // Same process for other rows ...

      // The current considered block
      SPC::SiconosMatrix  currentBlock = MIn->block(blockStart0, blockStart1);

      // dim of the subBlock of currentBlock to be set.
      Index subSize(2);
      // indices of the first element of MIn (resp. currentBlock) to be read (resp. set)  (same as start for MIn and MOut).
      Index currentPos(4);

      // currentBlock position in MOut.
      unsigned int numRow = blockStart0;
      unsigned int numCol = blockStart1;

      // Init currentPos
      // row and col position for first element to be read in MIn,
      currentPos[0] = posOut0;
      currentPos[1] = posOut1;
      // row and col position for first element in sub-block of Mout (namely currentBlock).
      currentPos[2] = start[2];
      currentPos[3] = start[3];

      while(numRow != blockEnd0 + 1)
      {

        while(numCol != blockEnd1 + 1)
        {
          // Get the block of MOut from which a sub-block will be set ...
          currentBlock = MIn->block(numRow, numCol);

          // Set subSize[0], dim (rows) and subSize[1], dim (columns) of the sub-block.
          // subSize[0] is only required for the first block in the row, after it remains constant.
          subSize[1] = currentBlock->size(1);
          // Warning: test "a" must be done before test "b"
          if(numCol == blockEnd1)  // if last column of blocks -> test "a"
            subSize[1] = lastBlockSize1;

          if(numCol == blockStart1)  // -> test "b"
          {
            subSize[1] -= posOut1;
            subSize[0] = currentBlock->size(0);
            if(numRow == blockEnd0)  // if last row of blocks
              subSize[0] = lastBlockSize0;
            if(numRow == blockStart0)  // if first row of blocks
              subSize[0] -= posOut0;
          }

          // Set sub-block
          setBlock(currentBlock, MOut, subSize, currentPos);

          // Update currentPos:
          // col position for first element to be read in MIn,
          currentPos[1] = 0 ;
          // col position for first element to be set in sub-block.
          currentPos[3] += subSize[1];
          numCol++;
        }

        numCol = blockStart1;
        numRow++;

        // Update currentPos:
        // row position for first element to be read in MIn,
        currentPos[0] = 0;
        // col position for first element to be read in MIn,
        currentPos[1] = posOut1;
        // row position for first element to be set in sub-block.
        currentPos[2] += subSize[0] ;
        // col position for first element to be set in sub-block.
        currentPos[3] = start[3];

      }
      MOut->resetLU();

    }
    else // neither MIn nor MOut is a BlockMatrix.
    {
      switch(numIn)
      {
      case 1:
        if(numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->dense(), start[0], end[0], start[1], end[1]);
        break;

      case 2:
        if(numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->triang(), start[0], end[0], start[1], end[1]);
        break;

      case 3:
        if(numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sym(), start[0], end[0], start[1], end[1]);
        break;

      case 4:
        if(numOut == 1)
          noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sparse(), start[0], end[0], start[1], end[1]);
        else if(numOut == 4)
          noalias(ublas::subrange(*MOut->sparse(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sparse(), start[0], end[0], start[1], end[1]);
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case 5:
        if(numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->banded(), start[0], end[0], start[1], end[1]);
        break;

      case 6:
        if(numOut == 1)
          ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if(numOut == 2)
          ublas::subrange(*MOut->triang(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if(numOut == 4)
          ublas::subrange(*MOut->sparse(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if(numOut == 5)
          ublas::subrange(*MOut->banded(), start[2], end[2], start[3], end[3]) *= 0.0;
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case 7:
        if(numOut == 1)
          noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->identity(), start[0], end[0], start[1], end[1]);
        else if(numOut == 4)
          noalias(ublas::subrange(*MOut->sparse(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->identity(), start[0], end[0], start[1], end[1]);
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      default:
        SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;
      }
      MOut->resetLU();
    }
  }
}


void SimpleMatrix::getRow(unsigned int r, SiconosVector &vOut) const
{
  // Get row number r of current matrix and copy it into vOut.
  if(r >= size(0))
    SiconosMatrixException::selfThrow("getRow(row): row is out of range");

  if(vOut.size() != size(1))
    SiconosMatrixException::selfThrow("getRow(row,v): inconsistent sizes between this and v.");

  if(_num == Siconos::IDENTITY)  // identity matrix
  {
    vOut.zero();
    vOut(r) = 1.0;
  }
  else if(_num == Siconos::ZERO)  // Zero matrix
    vOut.zero();
  else
  {
    unsigned int numV = vOut.num();
    if(numV == 1)
    {
      if(_num == Siconos::DENSE)
      {
        noalias(*(vOut.dense())) = ublas::row(*mat.Dense, r);
      }
      else if(_num == Siconos::TRIANGULAR)
      {
        noalias(*(vOut.dense())) = ublas::row(*mat.Triang, r);
      }
      else if(_num == Siconos::SYMMETRIC)
      {
        noalias(*(vOut.dense())) = ublas::row(*mat.Sym, r);
      }
      else if(_num == Siconos::SPARSE)
      {
        noalias(*(vOut.dense())) = ublas::row(*mat.Sparse, r);
      }
      else //if(_num == Siconos::BANDED){
        noalias(*(vOut.dense())) = ublas::row(*mat.Banded, r);
    }
    else // if numV == 4
    {
      if(_num == Siconos::SPARSE)
      {
        noalias(*(vOut.sparse())) = ublas::row(*mat.Sparse, r);
      }
      else
        SiconosMatrixException::selfThrow("getRow(row,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void SimpleMatrix::setRow(unsigned int r, const SiconosVector& vIn)
{
  // Set row number r of current matrix with vIn.
  unsigned int numV = vIn.num();
  if(r >= size(0))
    SiconosMatrixException::selfThrow("setRow(row): row is out of range");

  if(vIn.size() != size(1))
    SiconosMatrixException::selfThrow("setRow(row,v): inconsistent sizes between this and v.");

  if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("setRow(row,v): current matrix is read-only (zero or identity).");

  {
    if(_num == Siconos::DENSE)
    {
      if(numV == 1)
      {
        noalias(ublas::row(*mat.Dense, r)) = *vIn.dense();
      }
      else if(numV == 4)
      {
        noalias(ublas::row(*mat.Dense, r)) = *vIn.sparse();
      }
    }
    else if(_num == Siconos::SPARSE && numV == 4)
      noalias(ublas::row(*mat.Sparse, r)) = *vIn.sparse();
    else
      SiconosMatrixException::selfThrow("setRow(row,v): inconsistent types between current matrix and v.");
  }

  resetLU();
}

void SimpleMatrix::getCol(unsigned int r, SiconosVector &vOut)const
{
  // Get column number r of current matrix and copy it into vOut.
  if(r >= size(1))
    SiconosMatrixException::selfThrow("getCol(col): col is out of range");

  if(vOut.size() != size(0))
    SiconosMatrixException::selfThrow("getCol(col,v): inconsistent sizes between this and v.");

  if(_num == Siconos::IDENTITY)  // identity matrix
  {
    vOut.zero();
    vOut(r) = 1.0;
  }
  else if(_num == Siconos::ZERO)  // Zero matrix
    vOut.zero();
  else
  {
    unsigned int numV = vOut.num();

    if(numV == 1)
    {

      if(_num == Siconos::DENSE)
      {
        noalias(*(vOut.dense())) = ublas::column(*mat.Dense, r);
      }
      else if(_num == Siconos::TRIANGULAR)
      {
        noalias(*(vOut.dense())) = ublas::column(*mat.Triang, r);
      }
      else if(_num == Siconos::SYMMETRIC)
      {
        noalias(*(vOut.dense())) = ublas::column(*mat.Sym, r);
      }
      else if(_num == Siconos::SPARSE)
      {
        noalias(*(vOut.dense())) = ublas::column(*mat.Sparse, r);
      }
      else //if(_num == Siconos:BANDED){
        noalias(*(vOut.dense())) = ublas::column(*mat.Banded, r);
    }
    else // if _numV == 4
    {
      if(_num == Siconos::SPARSE)
      {
        noalias(*(vOut.sparse())) = ublas::column(*mat.Sparse, r);
      }
      else
        SiconosMatrixException::selfThrow("getCol(col,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void SimpleMatrix::setCol(unsigned int r, const SiconosVector &vIn)
{
  // Set column number r of current matrix with vIn.
  unsigned int numV = vIn.num();
  if(r >= size(1))
    SiconosMatrixException::selfThrow("setCol(col): col is out of range");

  if(vIn.size() != size(0))
    SiconosMatrixException::selfThrow("setCol(col,v): inconsistent sizes between this and v.");

  if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("setCol(col,v): current matrix is read-only (zero or identity).");

  {
    if(_num == Siconos::DENSE)
    {
      if(numV == 1)
      {
        noalias(ublas::column(*mat.Dense, r)) = *vIn.dense();
      }
      else if(numV == 4)
      {
        noalias(ublas::column(*mat.Dense, r)) = *vIn.sparse();
      }
    }
    else if(_num == Siconos::SPARSE && numV == 4)
      noalias(ublas::column(*mat.Sparse, r)) = *vIn.sparse();
    else
      SiconosMatrixException::selfThrow("setCol(col,v): inconsistent types between current matrix and v.");
  }

  resetLU();
}

void SimpleMatrix::getSubRow(unsigned int r, unsigned int pos, SP::SiconosVector vOut) const
{
  // Get row number r of current matrix, starting from element at position pos, and copy it into vOut.
  if(r >= size(0))
    SiconosMatrixException::selfThrow("getSubRow(row,pos,v): row is out of range");

  if(vOut->size() > size(1) - pos)
    SiconosMatrixException::selfThrow("getSubRow(row,pos,v): inconsistent sizes between this and v.");

  if(_num == Siconos::IDENTITY)  // identity matrix
  {
    vOut->zero();
    if(r >= pos)
      (*vOut)(r - pos) = 1.0;
  }
  else if(_num == Siconos::ZERO)  // Zero matrix
    vOut->zero();
  else
  {
    unsigned int numV = vOut->num();
    unsigned int nbEl = vOut->size();

    if(numV == 1)
    {
      if(_num == Siconos::DENSE)
      {
        //      noalias(*(vOut->dense())) = ublas::row(ublas::subrange(*mat.Dense, r, r+1,pos, endPos),0);
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if(_num == Siconos::TRIANGULAR)
      {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<TriangMat >(*mat.Triang, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if(_num == Siconos::SYMMETRIC)
      {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SymMat >(*mat.Sym, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if(_num == Siconos::SPARSE)
      {
        // #ifdef BOOST_LIMITATION
        //         SiconosMatrixException("SimpleMatrix::getSubRow warning - ublas::matrix_vector_slice<SparseMat> does not exist for your boost distribution and your architecture.");
        // #else
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
        // #endif
      }
      else //if(_num == Siconos::BANDED){
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<BandedMat >(*mat.Banded, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
    }
    else // if numV == 4
    {
      if(_num == Siconos::SPARSE)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubRow warning - ublas::matrix_vector_slice<SparseMat> does not exist for your boost distribution and your architecture.");
#else
        noalias(*(vOut->sparse())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
#endif
      }
      else
        SiconosMatrixException::selfThrow("getSubRow(row,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }

}

void SimpleMatrix::setSubRow(unsigned int r, unsigned int pos, SP::SiconosVector vIn)
{
  // Set row number r, starting from element at position pos, of current matrix with vIn.
  unsigned int numV = vIn->num();
  if(r >= size(0))
    SiconosMatrixException::selfThrow("setSubRow(row): row is out of range");

  if(vIn->size() > size(1) - pos)
    SiconosMatrixException::selfThrow("setSubRow(row,v): inconsistent sizes between this and v.");

  if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("setSubRow(row,v): current matrix is read-only (zero or identity).");

  {
    unsigned int nbEl = vIn->size();
    if(_num == Siconos::DENSE)
    {
      if(numV == 1)
      {
        noalias(ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl))) = *vIn->dense();
      }
      else if(numV == 4)
      {
        ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl)) = *vIn->sparse();
      }
    }
    else if(_num == Siconos::SPARSE && numV == 4)
#ifdef BOOST_LIMITATION
      SiconosMatrixException("SimpleMatrix::setSubRow warning - ublas::matrix_vector_slice<SparseMat> does not exist for your boost distribution and your architecture.");
#else
      ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl)) = *vIn->sparse();
#endif
    else
      SiconosMatrixException::selfThrow("setSubRow(row,v): inconsistent types between current matrix and v.");
    resetLU();
  }

}

void SimpleMatrix::getSubCol(unsigned int r, unsigned int pos, SP::SiconosVector vOut) const
{
  // Get col _number r of current matrix, starting from element at position pos, and copy it into vOut.
  if(r >= size(1))
    SiconosMatrixException::selfThrow("getSubCol(col,pos,v): col is out of range");

  if(vOut->size() > size(0) - pos)
    SiconosMatrixException::selfThrow("getSubCol(col,pos,v): inconsistent sizes between this and v.");

  if(_num == Siconos::IDENTITY)  // identity matrix
  {
    vOut->zero();
    if(r >= pos)
      (*vOut)(r - pos) = 1.0;
  }
  else if(_num == Siconos::ZERO)  // Zero matrix
    vOut->zero();
  else
  {
    unsigned int numV = vOut->num();
    unsigned int nbEl = vOut->size();

    if(numV == 1)
    {
      if(_num == Siconos::DENSE)
      {
        //      noalias(*(vOut->dense())) = ublas::row(ublas::subrange(*mat.Dense, r, r+1,pos, endPos),0);
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if(_num == Siconos::TRIANGULAR)
      {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<TriangMat >(*mat.Triang, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if(_num == Siconos::SYMMETRIC)
      {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SymMat >(*mat.Sym, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if(_num == Siconos::SPARSE)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubCol warning - ublas::matrix_vector_slice<SparseMat> does not exist for your boost distribution and your architecture.");
#else
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
#endif
      }
      else //if(_num == Siconos::BANDED){
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<BandedMat >(*mat.Banded, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
    }
    else // if numV == 4
    {
      if(_num == Siconos::SPARSE)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubCol warning - ublas::matrix_vector_slice<SparseMat> does not exist for your boost distribution and your architecture.");
#else
        noalias(*(vOut->sparse())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
#endif
      }
      else
        SiconosMatrixException::selfThrow("getSubCol(col,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }

}

void SimpleMatrix::setSubCol(unsigned int r, unsigned int pos, SP::SiconosVector vIn)
{
  // Set column number r, starting from element at position pos, of current matrix with vIn.
  unsigned int numV = vIn->num();
  if(r >= size(1))
    SiconosMatrixException::selfThrow("setSubCol(col): col is out of range");

  if(vIn->size() > size(0) - pos)
    SiconosMatrixException::selfThrow("setSubCol(col,v): inconsistent sizes between this and v.");

  if(_num == Siconos::ZERO || _num == Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("setSubCol(col,v): current matrix is read-only (zero or identity).");

  {
    unsigned int nbEl = vIn->size();
    if(_num == Siconos::DENSE)
    {
      if(numV == 1)
      {
        noalias(ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl))) = *vIn->dense();
      }
      else if(numV == 4)
      {
        ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl)) = *vIn->sparse();
      }
    }
    else if(_num == Siconos::SPARSE && numV == 4)
#ifdef BOOST_LIMITATION
      SiconosMatrixException("SimpleMatrix::setSubCol warning - ublas::matrix_vector_slice<SparseMat> does not exist for your boost distribution and your architecture.");
#else
      ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl)) = *vIn->sparse();
#endif
    else
      SiconosMatrixException::selfThrow("setSubCol(row,v): inconsistent types between current matrix and v.");
    resetLU();
  }
}


