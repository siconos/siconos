/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
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

#include "SiconosMatrixSetBlock.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "SiconosMatrix.hpp"

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

    if(numOut == Siconos::ZERO || numOut == Siconos::IDENTITY)  // if MOut = 0 or Identity => read-only
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
      case Siconos::DENSE:
        if(numOut != Siconos::DENSE)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->dense(), start[0], end[0], start[1], end[1]);
        break;

      case Siconos::TRIANGULAR:
        if(numOut != Siconos::DENSE)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->triang(), start[0], end[0], start[1], end[1]);
        break;

      case Siconos::SYMMETRIC:
        if(numOut != Siconos::DENSE)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sym(), start[0], end[0], start[1], end[1]);
        break;

      case Siconos::SPARSE:
        if(numOut == Siconos::DENSE)
          noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sparse(), start[0], end[0], start[1], end[1]);
        else if(numOut == Siconos::SPARSE)
          noalias(ublas::subrange(*MOut->sparse(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sparse(), start[0], end[0], start[1], end[1]);
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case Siconos::BANDED:
        if(numOut != Siconos::DENSE)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->banded(), start[0], end[0], start[1], end[1]);
        break;

      case Siconos::ZERO:
        if(numOut == Siconos::DENSE)
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

      case Siconos::IDENTITY:
        if(numOut == Siconos::DENSE)
          noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->identity(), start[0], end[0], start[1], end[1]);
        else if(numOut == Siconos::SPARSE)
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
