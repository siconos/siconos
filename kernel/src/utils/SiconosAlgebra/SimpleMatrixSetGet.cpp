/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "BlockMatrix.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

namespace ublas = boost::numeric::ublas;

double siconos::algebra::SimpleMatrix::getValue(unsigned int row, unsigned int col) const
{
  if (row >= size(0) || col >= size(1)) THROW_EXCEPTION("Index out of range");

  if (_num == UblasType::DENSE)
    return (*mat.Dense)(row, col);
  else if (_num == UblasType::TRIANGULAR)
    return (*mat.Triang)(row, col);
  else if (_num == UblasType::SYMMETRIC)
    return (*mat.Sym)(row, col);
  else if (_num == UblasType::SPARSE) {
    double *d = (*mat.Sparse).find_element(row, col);
    if (d)
      return *d;
    else
      return 0.0;
  }
  else if (_num == UblasType::SPARSE_COORDINATE) {
    double *d = (*mat.SparseCoordinate).find_element(row, col);
    if (d)
      return *d;
    else
      return 0.0;
  }
  else if (_num == UblasType::BANDED)
    return (*mat.Banded)(row, col);
  else if (_num == UblasType::ZERO)
    return 0;
  else  // if (_num == UblasType::IDENTITY)
    return (row == col);
}

void siconos::algebra::SimpleMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  if (row >= size(0) || col >= size(1)) THROW_EXCEPTION("Index out of range");

  if (_num == UblasType::DENSE)
    (*mat.Dense)(row, col) = value;
  else if (_num == UblasType::TRIANGULAR)
    (*mat.Triang)(row, col) = value;
  else if (_num == UblasType::SYMMETRIC)
    (*mat.Sym)(row, col) = value;
  else if (_num == UblasType::SPARSE) {
    double *d = (*mat.Sparse).find_element(row, col);
    if (d) {
      *d = value;
    }
    else {
      (*mat.Sparse).insert_element(row, col, value);
    }
  }
  else if (_num == UblasType::SPARSE_COORDINATE) {
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

  else if (_num == UblasType::BANDED)
    (*mat.Banded)(row, col) = value;
  else if (_num == UblasType::ZERO || _num == UblasType::IDENTITY)
    THROW_EXCEPTION("orbidden for Identity or Zero type matrices.");
  resetFactorizationFlags();
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void siconos::algebra::SimpleMatrix::setBlock(unsigned int row_min, unsigned int col_min,
                                              const SiconosMatrix &m)
{
  // Set current matrix elements, starting from row row_min and column col_min, with the values
  // of the matrix m. m may be a BlockMatrix.

  if (&m == this) THROW_EXCEPTION("m = this.");

  if (row_min >= size(0)) THROW_EXCEPTION("row is out of range");

  if (col_min >= size(1)) THROW_EXCEPTION("col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if (row_max > size(0)) THROW_EXCEPTION("row is out of range.");

  if (col_max > size(1)) THROW_EXCEPTION("m.col + col is out of range.");

  auto numM = m.num();

  if (numM == UblasType::BLOCK)  // if m is a block matrix ...
  {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    auto posRow = row_min;
    auto posCol = col_min;

    for (auto it = mB._mat->begin1(); it != mB._mat->end1(); ++it) {
      for (auto &it2 : it) {
        setBlock(posRow, posCol, *it2);
        posCol += it2->size(1);
      }
      posRow += (*it)->size(0);
      posCol = col_min;
    }
  }
  else  // if m is a SimpleMatrix
  {
    if (numM != _num) THROW_EXCEPTION("Inconsistent matrix types.");

    if (_num == UblasType::DENSE)
      noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) = *(m.dense());
    else if (_num == UblasType::TRIANGULAR)
      noalias(ublas::subrange(*mat.Triang, row_min, row_max, col_min, col_max)) =
          *(m.triang());
    else if (_num == UblasType::SYMMETRIC)
      noalias(ublas::subrange(*mat.Sym, row_min, row_max, col_min, col_max)) = *(m.sym());
    else if (_num == UblasType::SPARSE)
      noalias(ublas::subrange(*mat.Sparse, row_min, row_max, col_min, col_max)) =
          *(m.sparse());
    else if (_num == UblasType::BANDED)
      noalias(ublas::subrange(*mat.Banded, row_min, row_max, col_min, col_max)) =
          *(m.banded());
    else  // if(_num == UblasType::ZERO) or _num == UblasType::IDENTITY nothing to do
    {
    }
    resetFactorizationFlags();
  }
}

void siconos::algebra::SimpleMatrix::getRow(unsigned int r, SiconosVector &vOut) const
{
  // Get row number r of current matrix and copy it into vOut.
  if (r >= size(0)) THROW_EXCEPTION("row is out of range");

  if (vOut.size() != size(1)) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::IDENTITY)  // identity matrix
  {
    vOut.zero();
    vOut(r) = 1.0;
  }
  else if (_num == UblasType::ZERO)  // Zero matrix
    vOut.zero();
  else {
    auto numV = vOut.num();
    if (numV == UblasType::DENSE) {
      if (_num == UblasType::DENSE) {
        noalias(*(vOut.dense())) = ublas::row(*mat.Dense, r);
      }
      else if (_num == UblasType::TRIANGULAR) {
        noalias(*(vOut.dense())) = ublas::row(*mat.Triang, r);
      }
      else if (_num == UblasType::SYMMETRIC) {
        noalias(*(vOut.dense())) = ublas::row(*mat.Sym, r);
      }
      else if (_num == UblasType::SPARSE) {
        noalias(*(vOut.dense())) = ublas::row(*mat.Sparse, r);
      }
      else  // if(_num == UblasType::BANDED){
        noalias(*(vOut.dense())) = ublas::row(*mat.Banded, r);
    }
    else  // if numV == UblasType::SPARSE
    {
      if (_num == UblasType::SPARSE) {
        noalias(*(vOut.sparse())) = ublas::row(*mat.Sparse, r);
      }
      else
        THROW_EXCEPTION("inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void siconos::algebra::SimpleMatrix::setRow(unsigned int r, const SiconosVector &vIn)
{
  // Set row number r of current matrix with vIn.
  auto numV = vIn.num();
  if (r >= size(0)) THROW_EXCEPTION("row is out of range");

  if (vIn.size() != size(1)) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::ZERO || _num == UblasType::IDENTITY)
    THROW_EXCEPTION("current matrix is read-only (zero or identity).");

  {
    if (_num == UblasType::DENSE) {
      if (numV == UblasType::DENSE) {
        noalias(ublas::row(*mat.Dense, r)) = *vIn.dense();
      }
      else if (numV == UblasType::SPARSE) {
        noalias(ublas::row(*mat.Dense, r)) = *vIn.sparse();
      }
    }
    else if (_num == UblasType::SPARSE && numV == UblasType::SPARSE)
      noalias(ublas::row(*mat.Sparse, r)) = *vIn.sparse();
    else
      THROW_EXCEPTION("inconsistent types between current matrix and v.");
  }

  resetFactorizationFlags();
}

void siconos::algebra::SimpleMatrix::getCol(unsigned int r, SiconosVector &vOut) const
{
  // Get column number r of current matrix and copy it into vOut.
  if (r >= size(1)) THROW_EXCEPTION("col is out of range");

  if (vOut.size() != size(0)) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::IDENTITY)  // identity matrix
  {
    vOut.zero();
    vOut(r) = 1.0;
  }
  else if (_num == UblasType::ZERO)  // Zero matrix
    vOut.zero();
  else {
    auto numV = vOut.num();

    if (numV == UblasType::DENSE) {
      if (_num == UblasType::DENSE) {
        noalias(*(vOut.dense())) = ublas::column(*mat.Dense, r);
      }
      else if (_num == UblasType::TRIANGULAR) {
        noalias(*(vOut.dense())) = ublas::column(*mat.Triang, r);
      }
      else if (_num == UblasType::SYMMETRIC) {
        noalias(*(vOut.dense())) = ublas::column(*mat.Sym, r);
      }
      else if (_num == UblasType::SPARSE) {
        noalias(*(vOut.dense())) = ublas::column(*mat.Sparse, r);
      }
      else  // if(_num == Siconos:BANDED){
        noalias(*(vOut.dense())) = ublas::column(*mat.Banded, r);
    }
    else  // if _numV == UblasType::SPARSE
    {
      if (_num == UblasType::SPARSE) {
        noalias(*(vOut.sparse())) = ublas::column(*mat.Sparse, r);
      }
      else
        THROW_EXCEPTION("inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void siconos::algebra::SimpleMatrix::setCol(unsigned int r, const SiconosVector &vIn)
{
  // Set column number r of current matrix with vIn.
  auto numV = vIn.num();
  if (r >= size(1)) THROW_EXCEPTION("col is out of range");

  if (vIn.size() != size(0)) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::ZERO || _num == UblasType::IDENTITY)
    THROW_EXCEPTION("current matrix is read-only (zero or identity).");

  {
    if (_num == UblasType::DENSE) {
      if (numV == UblasType::DENSE) {
        noalias(ublas::column(*mat.Dense, r)) = *vIn.dense();
      }
      else if (numV == UblasType::SPARSE) {
        noalias(ublas::column(*mat.Dense, r)) = *vIn.sparse();
      }
    }
    else if (_num == UblasType::SPARSE && numV == UblasType::SPARSE)
      noalias(ublas::column(*mat.Sparse, r)) = *vIn.sparse();
    else
      THROW_EXCEPTION("nconsistent types between current matrix and v.");
  }

  resetFactorizationFlags();
}

void siconos::algebra::SimpleMatrix::getSubRow(unsigned int r, unsigned int pos,
                                               std::shared_ptr<SiconosVector> vOut) const
{
  // Get row number r of current matrix, starting from element at position pos, and copy it
  // into vOut.
  if (r >= size(0)) THROW_EXCEPTION("row is out of range");

  if (vOut->size() > size(1) - pos) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::IDENTITY)  // identity matrix
  {
    vOut->zero();
    if (r >= pos) (*vOut)(r - pos) = 1.0;
  }
  else if (_num == UblasType::ZERO)  // Zero matrix
    vOut->zero();
  else {
    auto numV = vOut->num();
    auto nbEl = vOut->size();

    if (numV == UblasType::DENSE) {
      if (_num == UblasType::DENSE) {
        //      noalias(*(vOut->dense())) = ublas::row(ublas::subrange(*mat.Dense, r, r+1,pos,
        //      endPos),0);
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<DenseMat>(
            *mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if (_num == UblasType::TRIANGULAR) {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<TriangMat>(
            *mat.Triang, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if (_num == UblasType::SYMMETRIC) {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SymMat>(
            *mat.Sym, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if (_num == UblasType::SPARSE) {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SparseMat>(
            *mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else  // if(_num == UblasType::BANDED){
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<BandedMat>(
            *mat.Banded, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
    }
    else  // if numV == UblasType::SPARSE
    {
      if (_num == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
        THROW_EXCEPTION(
            "ublas::matrix_vector_slice<SparseMat> does not exist for your boost "
            "distribution and your architecture.");
#else
        noalias(*(vOut->sparse())) = ublas::matrix_vector_slice<SparseMat>(
            *mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
#endif
      }
      else
        THROW_EXCEPTION("inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void siconos::algebra::SimpleMatrix::setSubRow(unsigned int r, unsigned int pos,
                                               std::shared_ptr<SiconosVector> vIn)
{
  // Set row number r, starting from element at position pos, of current matrix with vIn.
  auto numV = vIn->num();
  if (r >= size(0)) THROW_EXCEPTION("row is out of range");

  if (vIn->size() > size(1) - pos) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::ZERO || _num == UblasType::IDENTITY)
    THROW_EXCEPTION("current matrix is read-only (zero or identity).");

  {
    auto nbEl = vIn->size();
    if (_num == UblasType::DENSE) {
      if (numV == UblasType::DENSE) {
        noalias(ublas::matrix_vector_slice<DenseMat>(
            *mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl))) = *vIn->dense();
      }
      else if (numV == UblasType::SPARSE) {
        ublas::matrix_vector_slice<DenseMat>(*mat.Dense, ublas::slice(r, 0, nbEl),
                                             ublas::slice(pos, 1, nbEl)) = *vIn->sparse();
      }
    }
    else if (_num == UblasType::SPARSE && numV == UblasType::SPARSE)
#ifdef BOOST_LIMITATION
      THROW_EXCEPTION(
          "ublas::matrix_vector_slice<SparseMat> does not exist for your boost "
          "distribution and your architecture.");
#else
      ublas::matrix_vector_slice<SparseMat>(*mat.Sparse, ublas::slice(r, 0, nbEl),
                                            ublas::slice(pos, 1, nbEl)) = *vIn->sparse();
#endif
    else
      THROW_EXCEPTION("inconsistent types between current matrix and v.");
    resetFactorizationFlags();
  }
}

void siconos::algebra::SimpleMatrix::getSubCol(unsigned int r, unsigned int pos,
                                               std::shared_ptr<SiconosVector> vOut) const
{
  // Get col _number r of current matrix, starting from element at position pos, and copy it
  // into vOut.
  if (r >= size(1)) THROW_EXCEPTION("col is out of range");

  if (vOut->size() > size(0) - pos) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::IDENTITY)  // identity matrix
  {
    vOut->zero();
    if (r >= pos) (*vOut)(r - pos) = 1.0;
  }
  else if (_num == UblasType::ZERO)  // Zero matrix
    vOut->zero();
  else {
    auto numV = vOut->num();
    auto nbEl = vOut->size();

    if (numV == UblasType::DENSE) {
      if (_num == UblasType::DENSE) {
        //      noalias(*(vOut->dense())) = ublas::row(ublas::subrange(*mat.Dense, r, r+1,pos,
        //      endPos),0);
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<DenseMat>(
            *mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if (_num == UblasType::TRIANGULAR) {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<TriangMat>(
            *mat.Triang, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if (_num == UblasType::SYMMETRIC) {
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SymMat>(
            *mat.Sym, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if (_num == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
        THROW_EXCEPTION(
            "ublas::matrix_vector_slice<SparseMat> does not exist for your boost "
            "distribution and your architecture.");
#else
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<SparseMat>(
            *mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
#endif
      }
      else  // if(_num == UblasType::BANDED){
        noalias(*(vOut->dense())) = ublas::matrix_vector_slice<BandedMat>(
            *mat.Banded, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
    }
    else  // if numV == UblasType::SPARSE
    {
      if (_num == UblasType::SPARSE) {
#ifdef BOOST_LIMITATION
        THROW_EXCEPTION(
            "ublas::matrix_vector_slice<SparseMat> does not exist for your boost "
            "distribution and your architecture.");
#else
        noalias(*(vOut->sparse())) = ublas::matrix_vector_slice<SparseMat>(
            *mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
#endif
      }
      else
        THROW_EXCEPTION("inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void siconos::algebra::SimpleMatrix::setSubCol(unsigned int r, unsigned int pos,
                                               std::shared_ptr<SiconosVector> vIn)
{
  // Set column number r, starting from element at position pos, of current matrix with vIn.
  auto numV = vIn->num();
  if (r >= size(1)) THROW_EXCEPTION("col is out of range");

  if (vIn->size() > size(0) - pos) THROW_EXCEPTION("inconsistent sizes between this and v.");

  if (_num == UblasType::ZERO || _num == UblasType::IDENTITY)
    THROW_EXCEPTION("current matrix is read-only (zero or identity).");

  {
    auto nbEl = vIn->size();
    if (_num == UblasType::DENSE) {
      if (numV == UblasType::DENSE) {
        noalias(ublas::matrix_vector_slice<DenseMat>(
            *mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl))) = *vIn->dense();
      }
      else if (numV == UblasType::SPARSE) {
        ublas::matrix_vector_slice<DenseMat>(*mat.Dense, ublas::slice(pos, 1, nbEl),
                                             ublas::slice(r, 0, nbEl)) = *vIn->sparse();
      }
    }
    else if (_num == UblasType::SPARSE && numV == UblasType::SPARSE)
#ifdef BOOST_LIMITATION
      THROW_EXCEPTION(
          "ublas::matrix_vector_slice<SparseMat> does not exist for your boost "
          "distribution and your architecture.");
#else
      ublas::matrix_vector_slice<SparseMat>(*mat.Sparse, ublas::slice(pos, 1, nbEl),
                                            ublas::slice(r, 0, nbEl)) = *vIn->sparse();
#endif
    else
      THROW_EXCEPTION("inconsistent types between current matrix and v.");
    resetFactorizationFlags();
  }
}
