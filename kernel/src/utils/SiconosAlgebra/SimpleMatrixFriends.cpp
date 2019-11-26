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

#include "SimpleMatrixFriends.hpp"
#include <boost/numeric/ublas/io.hpp>    
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include "SiconosAlgebra.hpp"
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"


const SimpleMatrix matrix_pow(const SimpleMatrix& m, unsigned int power)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, pow function: not yet implemented for BlockMatrix.");
  if ( m.size(0) != m.size(1))
    SiconosMatrixException::selfThrow("matrix_pow(SimpleMatrix), matrix is not square.");

  if (power > 0)
  {
    unsigned int num = m.num();
    if (num == 1)
    {
      DenseMat p = *m.dense();
      for (unsigned int i = 1; i < power; i++)
        p = prod(p, *m.dense());
      return p;
    }
    else if (num == 2)
    {
      TriangMat t = *m.triang();
      for (unsigned int i = 1; i < power; i++)
        t = prod(t, *m.triang());
      return t;
    }
    else if (num == 3)
    {
      SymMat s = *m.sym();
      for (unsigned int i = 1; i < power; i++)
        s = prod(s, *m.sym());
      return s;
    }
    else if (num == 4)
    {
      SparseMat sp = *m.sparse();
      for (unsigned int i = 1; i < power; i++)
        sp = prod(sp, *m.sparse());
      return sp;
    }
    else if (num == 5)
    {
      DenseMat b = *m.banded();
      for (unsigned int i = 1; i < power; i++)
        b = prod(b, *m.banded());
      return b;
    }
    else if (num == 6)
    {
      ZeroMat z(m.size(0), m.size(1));
      return z;
    }
    else // if (num==7)
    {
      IdentityMat I(m.size(0), m.size(1));;
      return I;
    }
  }
  else// if(power == 0)
  {
    IdentityMat I = ublas::identity_matrix<double>(m.size(0), m.size(1));
    return I;
  }
}

/*
The following code inverts the matrix input using LU-decomposition with backsubstitution of unit vectors. Reference: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling & Flannery.

you can solve Ax=b using three lines of ublas code:

permutation_matrix<> piv;
lu_factorize(A, piv);
lu_substitute(A, piv, x);

*/
#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

// REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/* Matrix inversion routine.
Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T, class U, class V>
bool InvertMatrix(const ublas::matrix<T, U, V>& input, ublas::matrix<T, U, V>& inverse)
{
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
// create a working copy of the input
  matrix<T, U, V> A(input);
// create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

// perform LU-factorization
  int res = lu_factorize(A,pm);
  if (res != 0) return false;

// create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<T>(A.size1()));

// backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

#endif //INVERT_MATRIX_HPP

void invertMatrix(const SimpleMatrix& input, SimpleMatrix& output)
{
  InvertMatrix(*input.dense(), *output.dense());
}


/* XXX Find out if we can use an elementwise ublas operation */
SP::SiconosVector compareMatrices(const SimpleMatrix& data, const SimpleMatrix& ref)
{
  SimpleMatrix diff(data.size(0), data.size(1));
  SP::SiconosVector res(new SiconosVector(data.size(1)));
  diff = data - ref;
  for (unsigned int i = 0; i < data.size(0); ++i)
  {
    for (unsigned int j = 0; j < data.size(1); ++j)
      diff(i, j) /= 1 + fabs(ref(i, j));
  }
  diff.normInfByColumn(res);
  return res;

}

void setBlock(SPC::SiconosMatrix  MIn, SP::SiconosMatrix MOut, const Index& dim, const Index& start)
{
  // To copy a subBlock of MIn into a subBlock of MOut.
  // dim[0], dim[1]: number of rows and columns of the sub-block
  // start[0], start[1]: position (row, column) of the first element of the subBlock in MIn
  // start[2], start[3]: position (row, column) of the first element of the subBlock in MOut

  if (MIn == MOut) // useless op => nothing to be done.
  {}// SiconosVectorException::selfThrow("");
  else
  {
    unsigned int numIn = MIn->num();
    unsigned int numOut = MOut->num();

    if (numOut == 6 || numOut == 7) // if MOut = 0 or Identity => read-only
      SiconosMatrixException::selfThrow("matrices, setBlock(MIn, MOut...): MOut is read-only (zero or identity matrix?).");

    // Check dimension
    Index MDim(4); // dim. of matrices MIn and MOut.
    MDim[0] = MIn->size(0);
    MDim[1] = MIn->size(1);
    MDim[2] = MOut->size(0);
    MDim[3] = MOut->size(1);

    for (unsigned int i = 0; i < 4 ; ++i)
      if (start[i] >= MDim[i])
        SiconosMatrixException::selfThrow("matrices, setBlock(MIn, ...): sub-block indices are out of range.");

    // index position of the last element in subBlock ...
    Index end(4);
    end[0] = dim[0] + start[0];
    end[1] = dim[1] + start[1];
    end[2] = dim[0] + start[2];
    end[3] = dim[1] + start[3];

    for (unsigned int i = 0; i < 4 ; ++i)
      if (end[i] > MDim[i])
        SiconosMatrixException::selfThrow("matrices, setBlock(MIn, ...): sub-block indices are out of range.");

    // Elements from row/col start[i] to row/col (end[i]-1) will be copied.

    // If both matrices MIn and MOut are block, exception.
    if (numIn == 0 && numOut == 0)
      SiconosMatrixException::selfThrow("matrices, setBlock(MIn, MOut ...): not yet implemented for MIn and MOut both BlockMatrix. Try to use setBlock on the sub-matrices?");

    else if (numOut == 0) // if MOut is a BlockMatrix.
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
      while (start[2] >= (*tab)[blockStart0] && blockStart0 < tab->size())
        blockStart0 ++;
      // Relative position in the block blockStart0 of the first element to be set.
      unsigned int posOut0 = start[2];
      if (blockStart0 != 0)
        posOut0 -= (*tab)[blockStart0 - 1];

      unsigned int blockEnd0 = blockStart0;
      while (end[2] > (*tab)[blockEnd0] && blockEnd0 < tab->size())
        blockEnd0 ++;

      // Size of the last sub-block in the column of block
      unsigned int lastBlockSize0 = end[2];
      if (blockEnd0 != 0)
        lastBlockSize0 -= (*tab)[blockEnd0 - 1];

      // B - Block-Col position: we look for the block of MOut that include index start[3] and end[3].
      unsigned int blockStart1 = 0;
      tab = MOut->tabCol();
      while (start[3] >= (*tab)[blockStart1] && blockStart1 < tab->size())
        blockStart1 ++;
      // Relative position in the block blockStart1 of the first element to be set.
      unsigned int posOut1 = start[3];
      if (blockStart1 != 0)
        posOut1 -= (*tab)[blockStart1 - 1];

      unsigned int blockEnd1 = blockStart1;
      while (end[3] > (*tab)[blockEnd1] && blockEnd1 < tab->size())
        blockEnd1 ++;

      // Size of the last sub-block in the row of block
      unsigned int lastBlockSize1 = end[3];
      if (blockEnd1 != 0)
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

      while (numRow != blockEnd0 + 1)
      {

        while (numCol != blockEnd1 + 1)
        {
          // Get the block of MOut from which a sub-block will be set ...
          currentBlock = MOut->block(numRow, numCol);

          // Set subSize[0], dim (rows) and subSize[1], dim (columns) of the sub-block.
          // subSize[0] is only required for the first block in the row, after it remains constant.
          subSize[1] = currentBlock->size(1);

          // Warning: test "a" must be done before test "b"
          if (numCol == blockEnd1) // if last column of blocks -> test "a"
            subSize[1] = lastBlockSize1;

          if (numCol == blockStart1) // -> test "b"
          {
            subSize[1] -= posOut1;
            subSize[0] = currentBlock->size(0);
            if (numRow == blockEnd0) // if last row of blocks
              subSize[0] = lastBlockSize0;
            if (numRow == blockStart0) // if first row of blocks
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
    else if (numIn == 0) // If MIn is a BlockMatrix.
    {

      // Same process as for numOut == 0

      unsigned int blockStart0 = 0;
      SPC::Index tab = MIn->tabRow();
      while (start[0] >= (*tab)[blockStart0] && blockStart0 < tab->size())
        blockStart0 ++;
      // Relative position in the block blockStart0 of the first element to be set.
      unsigned int posOut0 = start[0];
      if (blockStart0 != 0)
        posOut0 -= (*tab)[blockStart0 - 1];

      unsigned int blockEnd0 = blockStart0;
      while (end[0] > (*tab)[blockEnd0] && blockEnd0 < tab->size())
        blockEnd0 ++;

      // Size of the last sub-block in the column of block
      unsigned int lastBlockSize0 = end[0];
      if (blockEnd0 != 0)
        lastBlockSize0 -= (*tab)[blockEnd0 - 1];

      // B - Block-Col position: we look for the block of MOut that include index start[3] and end[3].
      unsigned int blockStart1 = 0;
      tab = MIn->tabCol();
      while (start[1] >= (*tab)[blockStart1] && blockStart1 < tab->size())
        blockStart1 ++;
      // Relative position in the block blockStart1 of the first element to be set.
      unsigned int posOut1 = start[1];
      if (blockStart1 != 0)
        posOut1 -= (*tab)[blockStart1 - 1];

      unsigned int blockEnd1 = blockStart1;
      while (end[1] > (*tab)[blockEnd1] && blockEnd1 < tab->size())
        blockEnd1 ++;

      // Size of the last sub-block in the row of block
      unsigned int lastBlockSize1 = end[1];
      if (blockEnd1 != 0)
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

      while (numRow != blockEnd0 + 1)
      {

        while (numCol != blockEnd1 + 1)
        {
          // Get the block of MOut from which a sub-block will be set ...
          currentBlock = MIn->block(numRow, numCol);

          // Set subSize[0], dim (rows) and subSize[1], dim (columns) of the sub-block.
          // subSize[0] is only required for the first block in the row, after it remains constant.
          subSize[1] = currentBlock->size(1);
          // Warning: test "a" must be done before test "b"
          if (numCol == blockEnd1) // if last column of blocks -> test "a"
            subSize[1] = lastBlockSize1;

          if (numCol == blockStart1) // -> test "b"
          {
            subSize[1] -= posOut1;
            subSize[0] = currentBlock->size(0);
            if (numRow == blockEnd0) // if last row of blocks
              subSize[0] = lastBlockSize0;
            if (numRow == blockStart0) // if first row of blocks
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
      switch (numIn)
      {
      case 1:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->dense(), start[0], end[0], start[1], end[1]);
        break;

      case 2:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->triang(), start[0], end[0], start[1], end[1]);
        break;

      case 3:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sym(), start[0], end[0], start[1], end[1]);
        break;

      case 4:
        if (numOut == 1)
          noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sparse(), start[0], end[0], start[1], end[1]);
        else if (numOut == 4)
          noalias(ublas::subrange(*MOut->sparse(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->sparse(), start[0], end[0], start[1], end[1]);
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case 5:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->banded(), start[0], end[0], start[1], end[1]);
        break;

      case 6:
        if (numOut == 1)
          ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if (numOut == 2)
          ublas::subrange(*MOut->triang(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if (numOut == 4)
          ublas::subrange(*MOut->sparse(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if (numOut == 5)
          ublas::subrange(*MOut->banded(), start[2], end[2], start[3], end[3]) *= 0.0;
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case 7:
        if (numOut == 1)
          noalias(ublas::subrange(*MOut->dense(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->identity(), start[0], end[0], start[1], end[1]);
        else if (numOut == 4)
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

void prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int startRow = 0;
  VectorOfVectors::const_iterator it;
  // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
  for (it = y.begin(); it != y.end(); ++it)
  {
    A.private_prod(startRow, x, **it, init); ///// VIEWED
    startRow += (*it)->size();
  }
}

void subprod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, const Index& coord, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Number of the subvector of x that handles element at position coord[4]
  std::size_t firstBlockNum = x.getNumVectorAtPos(coord[4]);
  // Number of the subvector of x that handles element at position coord[5]
  unsigned int lastBlockNum = x.getNumVectorAtPos(coord[5]);
  Index subCoord = coord;
  SPC::SiconosVector  tmp = x[firstBlockNum];
  std::size_t subSize =  tmp->size(); // Size of the sub-vector
  const SP::Index xTab = x.tabIndex();
  if (firstBlockNum != 0)
  {
    subCoord[4] -= (*xTab)[firstBlockNum - 1];
    subCoord[5] =  std::min(coord[5] - (*xTab)[firstBlockNum - 1], subSize);
  }
  else
    subCoord[5] =  std::min(coord[5], subSize);

  if (firstBlockNum == lastBlockNum)
  {
    subprod(A, *tmp, y, subCoord, init);
  }
  else
  {
    unsigned int xPos = 0 ; // Position in x of the current sub-vector of x
    bool firstLoop = true;
    subCoord[3] = coord[2] + subCoord[5] - subCoord[4];
    for (VectorOfVectors::const_iterator it = x.begin(); it != x.end(); ++it)
    {
      if ((*it)->num() == 0)
        SiconosMatrixException::selfThrow("subprod(A,x,y) error: not yet implemented for x block of blocks ...");
      if (xPos >= firstBlockNum && xPos <= lastBlockNum)
      {
        tmp = x[xPos];
        if (firstLoop)
        {
          subprod(A, *tmp, y, subCoord, init);
          firstLoop = false;
        }
        else
        {
          subCoord[2] += subCoord[5] - subCoord[4]; // !! old values for 4 and 5
          subSize = tmp->size();
          subCoord[4] = 0;
          subCoord[5] = std::min(coord[5] - (*xTab)[xPos - 1], subSize);
          subCoord[3] = subCoord[2] + subCoord[5] - subCoord[4];
          subprod(A, *tmp, y, subCoord, false);
        }
      }
      xPos++;
    }
  }
}


// void private_addprod(const SiconosMatrix& A, unsigned int startRow, unsigned int startCol, const BlockVector& x, SiconosVector& y)
// {
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   assert(!A.isBlock() && "private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

//   VectorOfVectors::const_iterator it;
//   unsigned int startColBis = startCol;
//   for (it = x.begin(); it != x.end(); ++it)
//   {
//     private_addprod(A, startRow, startColBis, **it, y);
//     startColBis += (*it)->size();
//   }

// }

// x block, y siconos
// void private_prod(const SiconosMatrix& A, unsigned int startRow, const BlockVector& x, SiconosVector& y, bool init)
// {
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

//   if (init) // y = subA * x , else y += subA * x
//     y.zero();
//   private_addprod(A, startRow, 0, x, y);
// }

// // x and y blocks
// void private_prod(SPC::SiconosMatrix A, const unsigned int startRow, SPC::BlockVector x, SP::BlockVector y, bool init)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   unsigned int row = startRow;
//   VectorOfVectors::const_iterator it;
//   for (it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(*A, row, *x, **it, init);
//     row += (*it)->size();
//   }
// }


// // x and y blocks
// void private_prod(SPC::SiconosMatrix A, const unsigned int startRow, SPC::SiconosVector x, SP::BlockVector y, bool init)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   unsigned int row = startRow;
//   VectorOfVectors::const_iterator it;
//   for (it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(*A, row, *x, **it, init);
//     row += (*it)->size();
//   }
// }


// void private_addprod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   VectorOfVectors::const_iterator it;
//   unsigned int startColBis = startCol;
//   for (it = x->begin(); it != x->end(); ++it)
//   {
//     private_addprod((*it), A, startRow, startColBis, y);
//     startColBis += (*it)->size();
//   }

// }

// /// OK
// void private_prod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::BlockVector  y, bool init)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   unsigned int col = startCol;
//   VectorOfVectors::const_iterator it;
//   for (it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(x, A, col, *it, init);
//     col += (*it)->size();
//   }
// }

// /// OK
// void private_prod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::SiconosVector  y, bool init)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY
//   if (init) // y = subA * x , else y += subA * x
//     y->zero();
//   private_addprod(x, A, startCol, 0 , y);

// }

// void private_prod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::BlockVector  y, bool init)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   unsigned int col = startCol;
//   VectorOfVectors::const_iterator it;
//   for (it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(x, A, col, *it, init);
//     col += (*it)->size();
//   }
// }

// void private_addprod(double a, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SPC::SiconosVector x, SP::SiconosVector y)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   if (A->isBlock())
//     SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

//   // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
//   // Then computation of y = subA*x + y.
//   unsigned int numA = A->num();
//   unsigned int numY = y->num();
//   unsigned int numX = x->num();
//   unsigned int sizeX = x->size();
//   unsigned int sizeY = y->size();

//   if (numX != numY)
//     SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

//   if (numY == 1 && numX == 1)
//   {

//     assert(y->dense() != x->dense());

//     if (numA == 1)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->dense(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else if (numA == 2)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->triang(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else if (numA == 3)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->sym(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else if (numA == 4)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else //if(numA==5)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->banded(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//   }
//   else // x and y sparse
//   {
//     if (numA == 4)
//       *y->sparse() += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->sparse());
//     else
//       SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
//   }

// }

// void private_prod(double a, SPC::SiconosMatrix A, unsigned int startRow, SPC::SiconosVector x, SP::SiconosVector  y, bool init)
// {
//   assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

//   if (init) // y = subA * x , else y += subA * x
//     y->zero();
//   private_addprod(a, A, startRow, 0, x, y);

// }

