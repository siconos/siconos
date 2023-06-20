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

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "BlockMatrix.hpp"
#include "SiconosConfig.h"
#include "SiconosException.hpp"
#include "SimpleMatrix.hpp"

double &siconos::algebra::SimpleMatrix::operator()(unsigned int row, unsigned int col)
{
  assert((row < size(0) && col < size(1)) && "SimpleMatrix:operator(): Index out of range");

  if (_num == UblasType::DENSE)
    return (*mat.Dense)(row, col);
  else if (_num == UblasType::TRIANGULAR)
    return (*mat.Triang)(row, col);
  else if (_num == UblasType::SYMMETRIC)
    return (*mat.Sym)(row, col);
  else if (_num == UblasType::SPARSE) {
    double *d = (*mat.Sparse).find_element(row, col);
    if (d == nullptr) THROW_EXCEPTION("Index out of range");
    double &ref = *d;
    return ref;
  }
  else if (_num == UblasType::SPARSE_COORDINATE) {
    double *d = (*mat.SparseCoordinate).find_element(row, col);
    if (d == nullptr) THROW_EXCEPTION("Index out of range");
    double &ref = *d;
    return ref;
  }
  else if (_num == UblasType::BANDED)
    return (*mat.Banded)(row, col);
  else if (_num == UblasType::ZERO)
    return const_cast<double &>((*mat.Zero)(row, col));
  else if (_num == UblasType::IDENTITY)
    return const_cast<double &>((*mat.Identity)(row, col));
  else {
    THROW_EXCEPTION("invalid type of matrix");
  }
}

double siconos::algebra::SimpleMatrix::operator()(unsigned int row, unsigned int col) const
{
  assert((row < size(0) && col < size(1)) && "SimpleMatrix:operator(): Index out of range");

  if (_num == UblasType::DENSE)
    return (*mat.Dense)(row, col);
  else if (_num == UblasType::TRIANGULAR)
    return (*mat.Triang)(row, col);
  else if (_num == UblasType::SYMMETRIC)
    return (*mat.Sym)(row, col);
  else if (_num == UblasType::SPARSE)
    return (*mat.Sparse)(row, col);
  else if (_num == UblasType::SPARSE_COORDINATE)
    return (*mat.SparseCoordinate)(row, col);
  else if (_num == UblasType::BANDED)
    return (*mat.Banded)(row, col);
  else if (_num == UblasType::ZERO)
    return 0.0;
  else if (_num == UblasType::IDENTITY)
    return (row == col);
  else {
    THROW_EXCEPTION("invalid type of matrix");
  }
}

//=============
// Assignment
//=============

siconos::algebra::SimpleMatrix &siconos::algebra::SimpleMatrix::operator=(
    const SiconosMatrix &m)
{
  if (&m == this) return *this;  // auto-assignment.

  auto numM = m.num();

  if (size(0) != m.size(0) || size(1) != m.size(1)) {
    resize(m.size(0), m.size(1));
  }

  if (numM == UblasType::ZERO) {
    zero();
    return *this;
  }

  if (numM == UblasType::IDENTITY) {
    eye();
    return *this;
  }

  if (numM == UblasType::BLOCK) {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    unsigned int posRow = 0;
    unsigned int posCol = 0;

    for (auto it = mB._mat->begin1(); it != mB._mat->end1(); ++it) {
      for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else {
    switch (_num) {
      case UblasType::DENSE:
        switch (numM) {
          case UblasType::DENSE:
            noalias(*(mat.Dense)) = *m.dense();
            break;
          case UblasType::TRIANGULAR:
            noalias(*(mat.Dense)) = *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.Dense)) = *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.Dense)) = *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.Dense)) = *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.Dense)) = *m.banded();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::TRIANGULAR:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.Triang)) = *m.triang();
            break;
          default:
            THROW_EXCEPTION("assignment of a bad type of matrix into a triangular one.");
            break;
        }
        break;
      case UblasType::SYMMETRIC:
        if (numM == UblasType::IDENTITY || numM == UblasType::SYMMETRIC)
          noalias(*(mat.Sym)) = *m.sym();
        else
          THROW_EXCEPTION("bad assignment of matrix (symmetric one = dense or ...)");
        break;
      case UblasType::SPARSE:
        switch (numM) {
          case UblasType::DENSE:
            noalias(*(mat.Sparse)) = *m.dense();
            break;
          case UblasType::TRIANGULAR:
            noalias(*(mat.Sparse)) = *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.Sparse)) = *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.Sparse)) = *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.Sparse)) = *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.Sparse)) = *m.banded();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::SPARSE_COORDINATE:
        switch (numM) {
          case UblasType::DENSE:
            noalias(*(mat.SparseCoordinate)) = *m.dense();
            break;
          case UblasType::TRIANGULAR:
            noalias(*(mat.SparseCoordinate)) = *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.SparseCoordinate)) = *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.SparseCoordinate)) = *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.SparseCoordinate)) = *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.SparseCoordinate)) = *m.banded();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::BANDED:
        switch (numM) {
          case UblasType::BANDED:
            noalias(*(mat.Banded)) = *m.banded();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      default:
        THROW_EXCEPTION("invalid type of matrix");
        break;
    }
    resetFactorizationFlags();
  }
  return *this;
}

siconos::algebra::SimpleMatrix &siconos::algebra::SimpleMatrix::operator=(
    const SimpleMatrix &m)
{
  if (&m == this) return *this;  // auto-assignment.

  auto numM = m.num();

  if (size(0) != m.size(0) || size(1) != m.size(1)) resize(m.size(0), m.size(1));

  if (numM == UblasType::ZERO) {
    zero();
    return *this;
  }
  else if (numM == UblasType::IDENTITY) {
    eye();
    return *this;
  }

  switch (_num) {
    case UblasType::DENSE:
      switch (numM) {
        case UblasType::DENSE:
          noalias(*(mat.Dense)) = *m.dense();
          break;
        case UblasType::TRIANGULAR:
          noalias(*(mat.Dense)) = *m.triang();
          break;
        case UblasType::SYMMETRIC:
          noalias(*(mat.Dense)) = *m.sym();
          break;
        case UblasType::SPARSE:
          noalias(*(mat.Dense)) = *m.sparse();
          break;
        case UblasType::SPARSE_COORDINATE:
          noalias(*(mat.Dense)) = *m.sparseCoordinate();
          break;
        case UblasType::BANDED:
          noalias(*(mat.Dense)) = *m.banded();
          break;
        default:
          THROW_EXCEPTION("invalid type of matrix");
          break;
      }
      break;
    case UblasType::TRIANGULAR:
      switch (numM) {
        case UblasType::TRIANGULAR:
          noalias(*(mat.Triang)) = *m.triang();
          break;
        default:
          THROW_EXCEPTION("assignment of a bad type of matrix into a triangular one.");
          break;
      }
      break;
    case UblasType::SYMMETRIC:
      if (numM == UblasType::SYMMETRIC)
        noalias(*(mat.Sym)) = *m.sym();
      else
        THROW_EXCEPTION("bad assignment of matrix (symmetric one = dense or ...)");
      break;
    case UblasType::SPARSE:
      switch (numM) {
        case UblasType::DENSE:
          noalias(*(mat.Sparse)) = *m.dense();
          break;
        case UblasType::TRIANGULAR:
          noalias(*(mat.Sparse)) = *m.triang();
          break;
        case UblasType::SYMMETRIC:
          noalias(*(mat.Sparse)) = *m.sym();
          break;
        case UblasType::SPARSE:
          noalias(*(mat.Sparse)) = *m.sparse();
          break;
        case UblasType::SPARSE_COORDINATE:
          noalias(*(mat.Sparse)) = *m.sparseCoordinate();
	  break;
        case UblasType::BANDED:
          noalias(*(mat.Sparse)) = *m.banded();
          break;
        default:
          THROW_EXCEPTION("invalid type of matrix");
          break;
      }
      break;
    case UblasType::SPARSE_COORDINATE:
      switch (numM) {
        case UblasType::DENSE:
          noalias(*(mat.SparseCoordinate)) = *m.dense();
          break;
        case UblasType::TRIANGULAR:
          noalias(*(mat.SparseCoordinate)) = *m.triang();
          break;
        case UblasType::SYMMETRIC:
          noalias(*(mat.SparseCoordinate)) = *m.sym();
          break;
        case UblasType::SPARSE:
          noalias(*(mat.SparseCoordinate)) = *m.sparse();
          break;
        case UblasType::SPARSE_COORDINATE:
          noalias(*(mat.SparseCoordinate)) = *m.sparseCoordinate();
          break;
        case UblasType::BANDED:
          noalias(*(mat.SparseCoordinate)) = *m.banded();
          break;
        default:
          THROW_EXCEPTION("invalid type of matrix");
          break;
      }
      break;

    case UblasType::BANDED:
      switch (numM) {
        case UblasType::BANDED:
          noalias(*(mat.Banded)) = *m.banded();
          break;
        default:
          THROW_EXCEPTION("invalid type of matrix");
          break;
      }
      break;
    default:
      THROW_EXCEPTION("invalid type of matrix");
      break;
  }
  resetFactorizationFlags();
  return *this;
}

siconos::algebra::SimpleMatrix &siconos::algebra::SimpleMatrix::operator=(const DenseMat &m)
{
  if (_num != UblasType::DENSE) THROW_EXCEPTION("the current matrix is not dense.");

  if (size(0) != m.size1() || size(1) != m.size2()) THROW_EXCEPTION("Inconsistent sizes.");

  noalias(*(mat.Dense)) = m;

  resetFactorizationFlags();
  return *this;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

siconos::algebra::SimpleMatrix &siconos::algebra::SimpleMatrix::operator+=(
    const SiconosMatrix &m)
{
  auto numM = m.num();
  if (numM == UblasType::ZERO)  // m = 0
    return *this;

  if (&m == this)  // auto-assignment
  {
    switch (_num) {
      case UblasType::DENSE:
        *mat.Dense += *mat.Dense;
        break;
      case UblasType::TRIANGULAR:
        *mat.Triang += *mat.Triang;
        break;
      case UblasType::SYMMETRIC:
        *mat.Sym += *mat.Sym;
        break;
      case UblasType::SPARSE:
        *mat.Sparse += *mat.Sparse;
        break;
      case UblasType::SPARSE_COORDINATE:
        *mat.SparseCoordinate += *mat.SparseCoordinate;
        break;
      case UblasType::BANDED:
        *mat.Banded += *mat.Banded;
        break;
      default:
        THROW_EXCEPTION("invalid type of matrix");
    }
    resetFactorizationFlags();
    return *this;
  }

  if (size(0) != m.size(0) || size(1) != m.size(1)) resize(m.size(0), m.size(1));

  if (numM == UblasType::BLOCK) {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    unsigned int posRow = 0;
    unsigned int posCol = 0;
    // We scan all the blocks of m ...
    for (auto it1 = mB._mat->begin1(); it1 != mB._mat->end1(); ++it1) {
      for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
        addBlock(posRow, posCol, **it2);  // Each block of m is added into this.
        posCol += (*it2)->size(1);
      }
      posRow += (*it1)->size(0);
      posCol = 0;
    }
  }
  else  // if m is a SimpleMatrix
  {
    switch (_num) {
      case UblasType::DENSE:
        switch (numM) {
          case UblasType::DENSE:
            noalias(*(mat.Dense)) += *m.dense();
            break;
          case UblasType::TRIANGULAR:
            noalias(*(mat.Dense)) += *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.Dense)) += *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.Dense)) += *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.Dense)) += *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.Dense)) += *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Dense)) += *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::TRIANGULAR:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.Triang)) += *m.triang();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Triang)) += *m.identity();
            break;
          default:
            THROW_EXCEPTION("Operation not allowed (add in place into a triangular matrix).");
            break;
        }
        break;
      case UblasType::SYMMETRIC:
        if (numM == UblasType::SYMMETRIC)
          noalias(*(mat.Sym)) += *m.sym();
        else if (numM == UblasType::IDENTITY)
          noalias(*(mat.Sym)) += *m.identity();
        else
          THROW_EXCEPTION("bad assignment of matrix (symmetric one = dense or ...)");
        break;
      case UblasType::SPARSE:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.Sparse)) += *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.Sparse)) += *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.Sparse)) += *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.Sparse)) += *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.Sparse)) += *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Sparse)) += *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::SPARSE_COORDINATE:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.SparseCoordinate)) += *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.SparseCoordinate)) += *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.SparseCoordinate)) += *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.SparseCoordinate)) += *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.SparseCoordinate)) += *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.SparseCoordinate)) += *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::BANDED:
        switch (numM) {
          case UblasType::BANDED:
            noalias(*(mat.Banded)) += *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Banded)) += *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      default:
        THROW_EXCEPTION("invalid type of matrix");
        break;
    }
    resetFactorizationFlags();
  }
  return *this;
}

siconos::algebra::SimpleMatrix &siconos::algebra::SimpleMatrix::operator-=(
    const SiconosMatrix &m)
{
  auto numM = m.num();
  if (numM == UblasType::ZERO)  // m = 0
    return *this;

  if (&m == this)  // auto-assignment
  {
    switch (_num) {
      case UblasType::DENSE:
        *mat.Dense -= *mat.Dense;
        break;
      case UblasType::TRIANGULAR:
        *mat.Triang -= *mat.Triang;
        break;
      case UblasType::SYMMETRIC:
        *mat.Sym -= *mat.Sym;
        break;
      case UblasType::SPARSE:
        *mat.Sparse -= *mat.Sparse;
        break;
      case UblasType::SPARSE_COORDINATE:
        *mat.SparseCoordinate -= *mat.SparseCoordinate;
        break;
      case UblasType::BANDED:
        *mat.Banded -= *mat.Banded;
        break;
      default:
        THROW_EXCEPTION("invalid type of matrix");
    }
    resetFactorizationFlags();
    return *this;
  }
  if (size(0) != m.size(0) || size(1) != m.size(1)) THROW_EXCEPTION("inconsistent sizes.");

  if (numM == UblasType::BLOCK)  // m is a BlockMatrix
  {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    unsigned int posRow = 0;
    unsigned int posCol = 0;
    // We scan all the blocks of m ...
    for (auto it1 = mB._mat->begin1(); it1 != mB._mat->end1(); ++it1) {
      for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
        subBlock(posRow, posCol, **it2);  // Each block of m is added into this.
        posCol += (*it2)->size(1);
      }
      posRow += (*it1)->size(0);
      posCol = 0;
    }
  }
  else  // if m is a SimpleMatrix
  {
    switch (_num) {
      case UblasType::DENSE:
        switch (numM) {
          case UblasType::DENSE:
            noalias(*(mat.Dense)) -= *m.dense();
            break;
          case UblasType::TRIANGULAR:
            noalias(*(mat.Dense)) -= *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.Dense)) -= *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.Dense)) -= *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.Dense)) -= *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.Dense)) -= *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Dense)) -= *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::TRIANGULAR:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.Triang)) -= *m.triang();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Triang)) -= *m.identity();
            break;
          default:
            THROW_EXCEPTION("Operation not allowed on triangular matrix.");
            break;
        }
        break;
      case UblasType::SYMMETRIC:
        if (numM == UblasType::SYMMETRIC)
          noalias(*(mat.Sym)) -= *m.sym();
        else if (numM == UblasType::IDENTITY)
          noalias(*(mat.Sym)) -= *m.identity();
        else
          THROW_EXCEPTION("bad assignment of matrix (symmetric one = dense or ...)");
        break;
      case UblasType::SPARSE:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.Sparse)) -= *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.Sparse)) -= *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.Sparse)) -= *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.Sparse)) -= *m.sparse();
            break;
          case UblasType::BANDED:
            noalias(*(mat.Sparse)) -= *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Sparse)) -= *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      case UblasType::SPARSE_COORDINATE:
        switch (numM) {
          case UblasType::TRIANGULAR:
            noalias(*(mat.SparseCoordinate)) -= *m.triang();
            break;
          case UblasType::SYMMETRIC:
            noalias(*(mat.SparseCoordinate)) -= *m.sym();
            break;
          case UblasType::SPARSE:
            noalias(*(mat.SparseCoordinate)) -= *m.sparse();
            break;
          case UblasType::SPARSE_COORDINATE:
            noalias(*(mat.SparseCoordinate)) -= *m.sparseCoordinate();
            break;
          case UblasType::BANDED:
            noalias(*(mat.SparseCoordinate)) -= *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.SparseCoordinate)) -= *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;

      case UblasType::BANDED:
        switch (numM) {
          case UblasType::BANDED:
            noalias(*(mat.Banded)) -= *m.banded();
            break;
          case UblasType::IDENTITY:
            noalias(*(mat.Banded)) -= *m.identity();
            break;
          default:
            THROW_EXCEPTION("invalid type of matrix");
            break;
        }
        break;
      default:
        THROW_EXCEPTION("invalid type of matrix");
        break;
    }
    resetFactorizationFlags();
  }
  return *this;
}
