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

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"

using namespace Siconos;



void subprod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, const Index& coord, bool init)
{
  // To compute subY = subA * subX in an "optimized" way (in comparison with y = prod(A,x) )
  // or subY += subA*subX if init = false.

  // coord is [r0A r1A c0A c1A r0x r1x r0y r1y]
  //
  // subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A and c1A-1.
  // The same for x and y with rix and riy.

  // Check dims
  unsigned int rowA = coord[1] - coord[0];
  unsigned int colA = coord[3] - coord[2];
  unsigned int dimX = coord[5] - coord[4];
  unsigned int dimY = coord[7] - coord[6];
  if (colA != dimX)
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: inconsistent sizes between A and x.");

  if (rowA != dimY)
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: inconsistent sizes between A and y.");

  if (dimX > x.size() || dimY > y.size() || rowA > A.size(0) || colA > A.size(1))
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: input index too large.");

  unsigned int numA = A.num();
  unsigned int numX = x.num();
  unsigned int numY = y.num();

  if (numA == 0)  // If A,x or y is Block
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: not yet implemented for A block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
    {
      if (numY == 1)
        ublas::subrange(*y.dense(), coord[6], coord[7]) *= 0.0;
      else //if(numY==4)
        ublas::subrange(*y.sparse(), coord[6], coord[7]) *= 0.0;
    }
    //else nothing
  }
  else if (numA == 7) // A = identity
  {
    if (!init)
      ublas::subrange(*y.dense(), coord[6], coord[7]) += ublas::subrange(*x.dense(), coord[4], coord[5]);
    else
    {
      // if x and y do not share memory (ie are different objects)
      if (&x != &y)
        noalias(ublas::subrange(*y.dense(), coord[6], coord[7])) = ublas::subrange(*x.dense(), coord[4], coord[5]);

      // else nothing
    }
  }

  else // A is not 0 or identity
  {
    {
      if (init)
      {
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[4], coord[5]));

            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));

            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[4], coord[5]));
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));

              if (numY == 1)
              {
                ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
                noalias(subY) = ublas::prod(subA, subX);
              }
              else
              {
                ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[6], coord[7]));
                noalias(subY) = ublas::prod(subA, subX);
              }
#endif
            }
            else //if(numA==5)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> and vector_range<SparseVect> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
          }
        }
      }
      else // += case
      {
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[4], coord[5]));

            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));

            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[4], coord[5]));
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              if (numY == 1)
              {
                ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
                noalias(subY) += ublas::prod(subA, subX);
              }
              else
              {
                ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[6], coord[7]));
                noalias(subY) += ublas::prod(subA, subX);
              }
#endif
            }
            else //if(numA==5)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
          }
        }
      }
    }
  }
}


