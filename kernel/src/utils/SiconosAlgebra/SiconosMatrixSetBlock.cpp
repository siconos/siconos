/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "SiconosException.hpp"

void setBlock(SPC::SiconosMatrix  input_matrix, SP::SiconosMatrix output_matrix, const Index& dim, const Index& start)
{
  // To copy a subBlock of input_matrix into a subBlock of output_matrix.
  // dim[0], dim[1]: number of rows and columns of the sub-block
  // start[0], start[1]: position (row, column) of the first element of the subBlock in input_matrix
  // start[2], start[3]: position (row, column) of the first element of the subBlock in output_matrix

  if(input_matrix == output_matrix)  // useless op => nothing to be done.
    return;

  Siconos::UBLAS_TYPE numIn = input_matrix->num();
  Siconos::UBLAS_TYPE numOut = output_matrix->num();

  if(numOut == Siconos::ZERO || numOut == Siconos::IDENTITY)  // if output_matrix = 0 or Identity => read-only
    THROW_EXCEPTION("output_matrix is read-only (zero or identity matrix?).");

  // Check dimension
  Index MDim(4); // dim. of matrices input_matrix and output_matrix.
  MDim[0] = input_matrix->size(0);
  MDim[1] = input_matrix->size(1);
  MDim[2] = output_matrix->size(0);
  MDim[3] = output_matrix->size(1);

  for(unsigned int i = 0; i < 4 ; ++i)
    if(start[i] >= MDim[i])
      THROW_EXCEPTION("matrices, setBlock(input_matrix, ...): sub-block indices are out of range.");

  // index position of the last element in subBlock ...
  Index end(4);
  end[0] = dim[0] + start[0];
  end[1] = dim[1] + start[1];
  end[2] = dim[0] + start[2];
  end[3] = dim[1] + start[3];

  for(unsigned int i = 0; i < 4 ; ++i)
    if(end[i] > MDim[i])
      THROW_EXCEPTION("sub-block last indices are out of range.");

  // Elements from row/col start[i] to row/col (end[i]-1) will be copied.

  // If both matrices input_matrix and output_matrix are block, exception.
  if(numIn == Siconos::BLOCK && numOut == Siconos::BLOCK)
    THROW_EXCEPTION("not yet implemented for input_matrix and output_matrix both BlockMatrix. Try to use setBlock on the sub-matrices?");

  if(numOut == Siconos::BLOCK)  // if output_matrix is a BlockMatrix.
  {

    // Steps:
    // A - Find the blocks of output_matrix that "own" indices start[2] and end[2] ie
    //     the first and last sub-block to be set in a block-column
    //         --> numbers blockStart0 and blockEnd0
    // B - Find the  Block of output_matrix that "owns" index start[3] and end[3] ie
    //     the first sub-block to be set in a block-row
    //         --> numbers blockStart1 and blockEnd1
    //
    //        => The blocks concerned in output_matrix, are those between (block) rows blockStart0 and blockEnd0
    //           and (block) columns blockStart1 and blockEnd1.
    //
    // C - Loop through the concerned blocks (name = currentBlock) of output_matrix and call setBlock(input_matrix, currentBlock, subSize, currentPos).
    //     subSize: dim. of the considered sub-block of currentBlock to be set
    //     currentPos: same as "start" vector but for currentBlock
    //

    // A - Block-Row position: we look for the block of output_matrix that include index start[2] and end[2].
    //
    unsigned int blockStart0 = 0;
    SPC::Index tab = output_matrix->tabRow();
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

    // B - Block-Col position: we look for the block of output_matrix that include index start[3] and end[3].
    unsigned int blockStart1 = 0;
    tab = output_matrix->tabCol();
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
    SP::SiconosMatrix   currentBlock = output_matrix->block(blockStart0, blockStart1);

    // dim of the subBlock of currentBlock to be set.
    Index subSize(2);
    // indices of the first element of input_matrix (resp. currentBlock) to be read (resp. set)  (same as start for input_matrix and output_matrix).
    Index currentPos(4);

    // currentBlock position in output_matrix.
    unsigned int numRow = blockStart0;
    unsigned int numCol = blockStart1;

    // Init currentPos
    // row and col position for first element to be read in input_matrix,
    currentPos[0] = start[0];
    currentPos[1] = start[1];
    // row and col position for first element in sub-block of Mout (namely currentBlock).
    currentPos[2] = posOut0;
    currentPos[3] = posOut1;

    while(numRow != blockEnd0 + 1)
    {

      while(numCol != blockEnd1 + 1)
      {
        // Get the block of output_matrix from which a sub-block will be set ...
        currentBlock = output_matrix->block(numRow, numCol);

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
        setBlock(input_matrix, currentBlock, subSize, currentPos);

        // Update currentPos:
        // col position for first element to be read in input_matrix,
        currentPos[1] += subSize[1] ;
        // col position for first element to be set in sub-block.
        currentPos[3] = 0;
        numCol++;
      }

      numCol = blockStart1;
      numRow++;

      // Update currentPos:
      // row position for first element to be read in input_matrix,
      currentPos[0] += subSize[0] ;
      // col position for first element to be read in input_matrix,
      currentPos[1] = start[1] ;
      // row position for first element to be set in sub-block.
      currentPos[2] = 0;
      // col position for first element to be set in sub-block.
      currentPos[3] = posOut1;

    }

  }
  else if(numIn == Siconos::BLOCK)  // If input_matrix is a BlockMatrix.
  {

    // Same process as for numOut == 0

    unsigned int blockStart0 = 0;
    SPC::Index tab = input_matrix->tabRow();
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

    // B - Block-Col position: we look for the block of output_matrix that include index start[3] and end[3].
    unsigned int blockStart1 = 0;
    tab = input_matrix->tabCol();
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
    SPC::SiconosMatrix  currentBlock = input_matrix->block(blockStart0, blockStart1);

    // dim of the subBlock of currentBlock to be set.
    Index subSize(2);
    // indices of the first element of input_matrix (resp. currentBlock) to be read (resp. set)  (same as start for input_matrix and output_matrix).
    Index currentPos(4);

    // currentBlock position in output_matrix.
    unsigned int numRow = blockStart0;
    unsigned int numCol = blockStart1;

    // Init currentPos
    // row and col position for first element to be read in input_matrix,
    currentPos[0] = posOut0;
    currentPos[1] = posOut1;
    // row and col position for first element in sub-block of Mout (namely currentBlock).
    currentPos[2] = start[2];
    currentPos[3] = start[3];

    while(numRow != blockEnd0 + 1)
    {

      while(numCol != blockEnd1 + 1)
      {
        // Get the block of output_matrix from which a sub-block will be set ...
        currentBlock = input_matrix->block(numRow, numCol);

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
        setBlock(currentBlock, output_matrix, subSize, currentPos);

        // Update currentPos:
        // col position for first element to be read in input_matrix,
        currentPos[1] = 0 ;
        // col position for first element to be set in sub-block.
        currentPos[3] += subSize[1];
        numCol++;
      }

      numCol = blockStart1;
      numRow++;

      // Update currentPos:
      // row position for first element to be read in input_matrix,
      currentPos[0] = 0;
      // col position for first element to be read in input_matrix,
      currentPos[1] = posOut1;
      // row position for first element to be set in sub-block.
      currentPos[2] += subSize[0] ;
      // col position for first element to be set in sub-block.
      currentPos[3] = start[3];

    }
    output_matrix->resetFactorizationFlags();

  }
  else // neither input_matrix nor output_matrix is a BlockMatrix.
  {
    if(numOut == Siconos::DENSE)
    {
      ublas::matrix_range<DenseMat> out_range(*output_matrix->dense(),
                                              ublas::range(start[2],end[2]),
                                              ublas::range(start[3], end[3]));
      if(numIn == Siconos::DENSE)
      {
        ublas::matrix_range<DenseMat> in_range(*input_matrix->dense(),
                                               ublas::range(start[0],end[0]),
                                               ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::SYMMETRIC)
      {
        ublas::matrix_range<SymMat> in_range(*input_matrix->sym(),
                                             ublas::range(start[0],end[0]),
                                             ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::SPARSE)
      {
        ublas::matrix_range<SparseMat> in_range(*input_matrix->sparse(),
                                                ublas::range(start[0],end[0]),
                                                ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::IDENTITY)
      {
        ublas::matrix_range<IdentityMat> in_range(*input_matrix->identity(),
                                                  ublas::range(start[0],end[0]),
                                                  ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::ZERO)
        out_range *= 0.;
      else
        THROW_EXCEPTION("unconsistent types between input_matrix and output_matrix.");
    }
    else if(numOut == Siconos::SPARSE)
    {
      ublas::matrix_range<SparseMat> out_range(*output_matrix->sparse(),
                                               ublas::range(start[2],end[2]),
                                               ublas::range(start[3], end[3]));
      if(numIn == Siconos::DENSE)
      {
        ublas::matrix_range<DenseMat> in_range(*input_matrix->dense(),
                                               ublas::range(start[0],end[0]),
                                               ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::SYMMETRIC)
      {
        ublas::matrix_range<SymMat> in_range(*input_matrix->sym(),
                                             ublas::range(start[0],end[0]),
                                             ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::SPARSE)
      {
        ublas::matrix_range<SparseMat> in_range(*input_matrix->sparse(),
                                                ublas::range(start[0],end[0]),
                                                ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::IDENTITY)
      {
        ublas::matrix_range<IdentityMat> in_range(*input_matrix->identity(),
                                                  ublas::range(start[0],end[0]),
                                                  ublas::range(start[1], end[1]));
        noalias(out_range) = in_range;
      }
      else if(numIn == Siconos::ZERO)
        out_range *= 0.;
      else
        THROW_EXCEPTION("unconsistent types between input_matrix and output_matrix.");
    }

    output_matrix->resetFactorizationFlags();
  }
}
