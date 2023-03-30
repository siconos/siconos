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

// #include <boost/numeric/bindings/std/vector.hpp>
// #include <boost/numeric/bindings/ublas/matrix.hpp>
// #include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include "BlockMatrix.hpp"
#include "BlockVector.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for subprod
#include "SimpleMatrix.hpp"

namespace ublas = boost::numeric::ublas;

void siconos::algebra::subprod(const SiconosMatrix &A, const SiconosVector &x,
                               SiconosVector &y, const std::vector<std::size_t> &coord,
                               bool init) {
  // To compute subY = subA * subX in an "optimized" way (in comparison with y = prod(A,x) )
  // or subY += subA*subX if init = false.

  // coord is [r0A r1A c0A c1A r0x r1x r0y r1y]
  //
  // subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A
  // and c1A-1. The same for x and y with rix and riy.

  // Check dims
  auto rowA = coord[1] - coord[0];
  auto colA = coord[3] - coord[2];
  auto dimX = coord[5] - coord[4];
  auto dimY = coord[7] - coord[6];
  if (colA != dimX) THROW_EXCEPTION("inconsistent sizes between A and x.");

  if (rowA != dimY) THROW_EXCEPTION("inconsistent sizes between A and y.");

  if (dimX > x.size() || dimY > y.size() || rowA > A.size(0) || colA > A.size(1))
    THROW_EXCEPTION("input index too large.");

  auto numA = A.num();
  auto numX = x.num();
  auto numY = y.num();

  if (numA == UblasType::BLOCK)  // If A,x or y is Block
    THROW_EXCEPTION("not yet implemented for A block matrices.");

  if (numA == UblasType::ZERO)  // A = 0
  {
    if (init) {
      if (numY == UblasType::DENSE)
        ublas::subrange(*y.dense(), coord[6], coord[7]) *= 0.0;
      else  // if(numY==UblasType::SPARSE)
        ublas::subrange(*y.sparse(), coord[6], coord[7]) *= 0.0;
    }
    // else nothing
  } else if (numA == UblasType::IDENTITY)  // A = identity
  {
    if (!init)
      ublas::subrange(*y.dense(), coord[6], coord[7]) +=
          ublas::subrange(*x.dense(), coord[4], coord[5]);
    else {
      // if x and y do not share memory (ie are different objects)
      if (&x != &y)
        noalias(ublas::subrange(*y.dense(), coord[6], coord[7])) =
            ublas::subrange(*x.dense(), coord[4], coord[5]);

      // else nothing
    }
  }

  else  // A is not 0 or identity
  {
    {
      if (init) {
        if (&x != &y)  // if no common memory between x and y.
        {
          if (numX == UblasType::DENSE) {
            ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[4], coord[5]));

            if (numY != UblasType::DENSE)
              THROW_EXCEPTION("y (output) must be a dense vector.");
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));

            if (numA == UblasType::DENSE) {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
          } else  // if(numX == UblasType::SPARSE)
          {
            ublas::vector_range<SparseVect> subX(*x.sparse(),
                                                 ublas::range(coord[4], coord[5]));
            if (numY != UblasType::DENSE && numA != UblasType::SPARSE)
              THROW_EXCEPTION("y (output) must be a dense vector.");

            if (numA == UblasType::DENSE) {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));

              if (numY == UblasType::DENSE) {
                ublas::vector_range<DenseVect> subY(*y.dense(),
                                                    ublas::range(coord[6], coord[7]));
                noalias(subY) = ublas::prod(subA, subX);
              } else {
                ublas::vector_range<SparseVect> subY(*y.sparse(),
                                                     ublas::range(coord[6], coord[7]));
                noalias(subY) = ublas::prod(subA, subX);
              }
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
          }
        } else  // if x and y are the same object => alias
        {
          if (numX == UblasType::DENSE) {
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[4], coord[5]));
            if (numA == UblasType::DENSE) {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> and vector_range<SparseVect> does not exist "
                  "for your boost distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
          } else  // if(numX == UblasType::SPARSE)
          {
            ublas::vector_range<SparseVect> subY(*y.sparse(),
                                                 ublas::range(coord[4], coord[5]));
            if (numA == UblasType::DENSE) {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
          }
        }
      } else  // += case
      {
        if (&x != &y)  // if no common memory between x and y.
        {
          if (numX == UblasType::DENSE) {
            ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[4], coord[5]));

            if (numY != UblasType::DENSE)
              THROW_EXCEPTION("y (output) must be a dense vector.");
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[6], coord[7]));

            if (numA == UblasType::DENSE) {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
          } else  // if(numX == UblasType::SPARSE)
          {
            ublas::vector_range<SparseVect> subX(*x.sparse(),
                                                 ublas::range(coord[4], coord[5]));
            if (numY != UblasType::DENSE && numA != UblasType::SPARSE)
              THROW_EXCEPTION("y (output) must be a dense vector.");

            if (numA == UblasType::DENSE) {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              if (numY == UblasType::DENSE) {
                ublas::vector_range<DenseVect> subY(*y.dense(),
                                                    ublas::range(coord[6], coord[7]));
                noalias(subY) += ublas::prod(subA, subX);
              } else {
                ublas::vector_range<SparseVect> subY(*y.sparse(),
                                                     ublas::range(coord[6], coord[7]));
                noalias(subY) += ublas::prod(subA, subX);
              }
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::vector_range<DenseVect> subY(*y.dense(),
                                                  ublas::range(coord[6], coord[7]));
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
          }
        } else  // if x and y are the same object => alias
        {
          if (numX == UblasType::DENSE) {
            ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[4], coord[5]));
            if (numA == UblasType::DENSE) {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
          } else  // if(numX == UblasType::SPARSE)
          {
            ublas::vector_range<SparseVect> subY(*y.sparse(),
                                                 ublas::range(coord[4], coord[5]));
            if (numA == UblasType::DENSE) {
              ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                                 ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            } else if (numA == UblasType::TRIANGULAR) {
              ublas::matrix_range<TriangMat> subA(*A.triang(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            } else if (numA == UblasType::SYMMETRIC) {
              ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                               ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            } else if (numA == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
              THROW_EXCEPTION(
                  "ublas::matrix_range<SparseMat> does not exist for your boost "
                  "distribution and your architecture.");
#else
              ublas::matrix_range<SparseMat> subA(*A.sparse(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
#endif
            } else  // if(numA==UblasType::BANDED)
            {
              ublas::matrix_range<BandedMat> subA(*A.banded(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
          }
        }
      }
    }
  }
}

void siconos::algebra::subprod(const SiconosMatrix &A, const BlockVector &x, SiconosVector &y,
                               const std::vector<std::size_t> &coord, bool init) {
  assert(!(A.isFactorized()) && "A is Factorized in prod !!");

  // Number of the subvector of x that handles element at position coord[4]
  auto firstBlockNum = x.getNumVectorAtPos(coord[4]);
  // Number of the subvector of x that handles element at position coord[5]
  auto lastBlockNum = x.getNumVectorAtPos(coord[5]);
  auto subCoord = coord;
  auto tmp = x.vector(firstBlockNum);
  auto subSize = static_cast<size_t>(tmp->size());  // Size of the sub-vector
  auto xTab = x.tabIndex();
  if (firstBlockNum != 0) {
    subCoord[4] -= (*xTab)[firstBlockNum - 1];
    subCoord[5] = std::min(coord[5] - (*xTab)[firstBlockNum - 1], subSize);
  } else
    subCoord[5] = std::min(coord[5], subSize);

  if (firstBlockNum == lastBlockNum) {
    subprod(A, *tmp, y, subCoord, init);
  } else {
    decltype(firstBlockNum) xPos = 0;  // Position in x of the current sub-vector of x
    bool firstLoop = true;
    subCoord[3] = coord[2] + subCoord[5] - subCoord[4];
    for (const auto &vec : x) {
      if (vec->num() == UblasType::BLOCK)
        THROW_EXCEPTION("not yet implemented for x block of blocks ...");
      if (xPos >= firstBlockNum && xPos <= lastBlockNum) {
        tmp = x.vector(xPos);
        if (firstLoop) {
          subprod(A, *tmp, y, subCoord, init);
          firstLoop = false;
        } else {
          subCoord[2] += subCoord[5] - subCoord[4];  // !! old values for 4 and 5
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
