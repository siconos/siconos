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

#include "SiconosMatrix.hpp"
#include <assert.h>                                   // for assert
#include <math.h>                                     // for fabs
#include <algorithm>                                  // for max, min, lower...
#include <boost/numeric/ublas/detail/config.hpp>      // for noalias, noalia...
#include <boost/numeric/ublas/detail/iterator.hpp>    // for bidirectional_i...
#include <boost/numeric/ublas/matrix_expression.hpp>  // for matrix_vector_b...
#include <boost/numeric/ublas/matrix_proxy.hpp>       // for matrix_range
#include <boost/numeric/ublas/matrix_sparse.hpp>      // for compressed_matr...
#include <boost/numeric/ublas/storage.hpp>            // for unbounded_array
#include <boost/numeric/ublas/vector.hpp>             // for vector
#include <memory>                                     // for __shared_ptr_ac...
#include <utility>                                    // for swap
#include <vector>                                     // for vector, operator==
#include "SiconosAlgebra.hpp"
#include "BlockMatrix.hpp"                            // for BlockMatrix
#include "CSparseMatrix.h"                            // for CSparseMatrix
#include "SiconosVector.hpp"                          // for SiconosVector


// Constructor with the type-number
SiconosMatrix::SiconosMatrix(unsigned int type): _num(type)
{}

const SP::Index SiconosMatrix::tabRow() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::tabRow() : not implemented for this type of matrix (Simple?) reserved to BlockMatrix.");
  // fake to avoid error on warning.
  return SP::Index();
}

const SP::Index SiconosMatrix::tabCol() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::tabCol() : not implemented for this type of matrix (Simple?) reserved to BlockMatrix.");
  // fake to avoid error on warning.
  return SP::Index();
}



//=====================
// matrices comparison
//=====================
bool isComparableTo(const  SiconosMatrix& m1, const  SiconosMatrix& m2)
{
  // return:
  // - true if one of the matrices is a Simple and if they have the same dimensions.
  // - true if both are block but with blocks which are facing each other of the same size.
  // - false in other cases

  if((!m1.isBlock() || !m2.isBlock()) && (m1.size(0) == m2.size(0)) && (m1.size(1) == m2.size(1)))
    return true;

  const SP::Index I1R = m1.tabRow();
  const SP::Index I2R = m2.tabRow();
  const SP::Index I1C = m1.tabCol();
  const SP::Index I2C = m2.tabCol();

  return ((*I1R == *I2R) && (*I1C == *I2C));
}

SiconosMatrix& operator *=(SiconosMatrix& m, const double& s)
{
  if(m._num == 0)  // BlockMatrix
  {
    BlockMatrix& mB = static_cast<BlockMatrix&>(m);
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    for(it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for(it2 = it.begin(); it2 != it.end(); ++it2)
        (**it2) *= s;
    }
  }
  else if(m._num == Siconos::DENSE)
    *m.dense() *= s;
  else if(m._num == Siconos::TRIANGULAR)
    *m.triang() *= s;
  else if(m._num == Siconos::SYMMETRIC)
    *m.sym() *= s;
  else if(m._num == Siconos::SPARSE)
    *m.sparse() *= s;
  else if(m._num == Siconos::SPARSE_COORDINATE)
    *m.sparseCoordinate() *= s;
  else if(m._num == Siconos::BANDED)
    *m.banded() *= s;
  else if(m._num == Siconos::ZERO) {}  // nothing!
  else //if(_num == 7)
    SiconosMatrixException::selfThrow(" SP::SiconosMatrix = (double) : invalid type of matrix");

  return m;
}

SiconosMatrix& operator /=(SiconosMatrix& m, const double& s)
{
  if(m._num == 0)  // BlockMatrix
  {
    BlockMatrix& mB = static_cast<BlockMatrix&>(m);
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    for(it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for(it2 = it.begin(); it2 != it.end(); ++it2)
        (**it2) /= s;
    }
  }
  else if(m._num == Siconos::DENSE)
    *m.dense() /= s;
  else if(m._num == Siconos::TRIANGULAR)
    *m.triang() /= s;
  else if(m._num == Siconos::SYMMETRIC)
    *m.sym() /= s;
  else if(m._num == Siconos::SPARSE)
    *m.sparse() /= s;
  else if(m._num == Siconos::SPARSE_COORDINATE)
    *m.sparseCoordinate() /= s;
  else if(m._num == Siconos::BANDED)
    *m.banded() /= s;
  else if(m._num == Siconos::ZERO) {}  // nothing!
  else //if(_num == 7)
    SiconosMatrixException::selfThrow(" SiconosMatrix *= (double) : invalid type of matrix");

  return m;
}

size_t SiconosMatrix::nnz(double tol)
{
  size_t nnz = 0;
  if(_num == Siconos::DENSE)  //dense
  {
    double* arr = getArray();
    for(size_t i = 0; i < size(0)*size(1); ++i)
    {
      if(fabs(arr[i]) > tol)
      {
        nnz++;
      }
    }
  }
  else if(_num == Siconos::SPARSE)
  {
    nnz = sparse()->nnz();
  }
  else
  {
    SiconosMatrixException::selfThrow("SiconosMatrix::nnz not implemented for the given matrix type");
  }

  return nnz;

}

bool SiconosMatrix::fillCSC(CSparseMatrix* csc, size_t row_off, size_t col_off, double tol)
{
  assert(csc);
  double* Mx = csc->x; // data
  CS_INT* Mi = csc->i; // row indx
  CS_INT* Mp = csc->p; // column pointers

  assert(Mp[col_off] >= 0);
  size_t nz = csc->p[col_off];

  size_t nrow = size(0);
  size_t ncol = size(1);

  CS_INT pval = Mp[col_off];

  if(_num == Siconos::DENSE)  //dense
  {
    double* arr = getArray();
    for(size_t j = 0, joff = col_off; j < ncol; ++j)
    {
      for(size_t i = 0; i < nrow; ++i)
      {
        // col-major
        double elt_val = arr[i + j*nrow];
        // std::cout << " a(i=" << i << ",j=" << j << ") = "<< elt_val << std::endl;
        if(fabs(elt_val) > tol)
        {
          Mx[pval] = elt_val;
          Mi[pval] = i + row_off;
          // std::cout << "Mx[" <<pval <<"] = " << Mx[pval]<<   std::endl;
          // std::cout << "Mp[" <<pval <<"] = " << Mi[pval]<<   std::endl;
          ++pval;
        }
      }
      // std::cout << "joff" << joff << std::endl;
      Mp[++joff] = pval;

    }
  }
  else if(_num == Siconos::SPARSE)
  {
    const Index& ptr = sparse()->index1_data();
    const Index& indx = sparse()->index2_data();
    const ublas::unbounded_array<double>& vals = sparse()->value_data();

    size_t nnz =  sparse()->nnz();

    assert(ptr.size() == ncol + 1);
    assert(indx.size() == nnz);
    assert(vals.size() == nnz);

    for(size_t i = 0; i < nnz; ++i)
    {
      Mx[pval] = vals[i];
      Mi[pval++] = row_off + indx[i];
    }
    for(size_t j = 1, joff = col_off + 1; j < ncol+1; ++j, ++joff)
    {
      Mp[joff] = nz + ptr[j];
    }
  }
  else
  {
    SiconosMatrixException::selfThrow("SiconosMatrix::fillCSC not implemented for the given matrix type");
  }

  return true;
}

bool SiconosMatrix::fillCSC(CSparseMatrix* csc, double tol)
{
  assert(csc);
  double* Mx = csc->x; // data
  CS_INT* Mi = csc->i; // row indx
  CS_INT* Mp = csc->p; // column pointers

  size_t nrow = size(0);
  size_t ncol = size(1);

  CS_INT pval = 0;

  if(_num == Siconos::DENSE)  //dense
  {
    double* arr = getArray();
    for(size_t j = 0, joff = 0; j < ncol; ++j)
    {
      for(size_t i = 0; i < nrow; ++i)
      {
        // col-major
        double elt_val = arr[i + j*nrow];
        // std::cout << " a(i=" << i << ",j=" << j << ") = "<< elt_val << std::endl;
        if(fabs(elt_val) > tol)
        {
          Mx[pval] = elt_val;
          Mi[pval] = i;
          // std::cout << "Mx[" <<pval <<"] = " << Mx[pval]<<   std::endl;
          // std::cout << "Mp[" <<pval <<"] = " << Mi[pval]<<   std::endl;
          ++pval;
        }
      }
      // std::cout << "joff" << joff << std::endl;
      Mp[++joff] = pval;

    }
  }
  else if(_num == Siconos::SPARSE)
  {
    const Index& ptr = sparse()->index1_data();
    const Index& indx = sparse()->index2_data();
    const ublas::unbounded_array<double>& vals = sparse()->value_data();

    size_t nnz =  sparse()->nnz();

    assert(ptr.size() == ncol + 1);
    assert(indx.size() >= nnz);
    assert(vals.size() >= nnz);

    for(size_t i = 0; i < nnz; ++i)
    {
      Mx[pval] = vals[i];
      Mi[pval++] = indx[i];
    }
    for(size_t j = 0; j < ncol+1; ++j)
    {
      Mp[j] = ptr[j];
    }
  }
  else
  {
    SiconosMatrixException::selfThrow("SiconosMatrix::fillCSC not implemented for the given matrix type");
  }

  return true;
}

bool SiconosMatrix::fillTriplet(CSparseMatrix* triplet, size_t row_off, size_t col_off, double tol)
{
  assert(triplet);
  size_t nrow = size(0);
  size_t ncol = size(1);

  if(_num == Siconos::DENSE)  //dense
  {
    double* arr = getArray();
    for(size_t j = 0; j < ncol; ++j)
    {
      for(size_t i = 0; i < nrow; ++i)
      {
        // col-major

        CSparseMatrix_zentry(triplet, i + row_off, j + col_off, arr[i + j*nrow]);
      }
    }
  }
  else
  {
    SiconosMatrixException::selfThrow("SiconosMatrix::fillCSC not implemented for the given matrix type");
  }

  return true;
}


std::ostream& operator<<(std::ostream& os, const SiconosMatrix& sm)
{
  os << sm.toString();
  return os;
}


void SiconosMatrix::private_prod(unsigned int startRow, const SiconosVector& x, SiconosVector& y, bool init) const
{
  assert(!(isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

  if(init)  // y = subA * x , else y += subA * x
    y.zero();
  private_addprod(startRow, 0, x, y);
}

/** Computation of y = subA.x

    where subA is a sub-matrix of A, subA[0,0] = A[startRow, startCol]
 */
void SiconosMatrix::private_addprod(unsigned startRow, unsigned int startCol, const SiconosVector& x, SiconosVector& y) const
{
  assert(!(isPLUFactorized()) && "A is PLUFactorized in prod !!");
  assert(!isBlock() && "private_addprod(start,x,y) error: not yet implemented for block matrix.");

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
  // Then computation of y = subA*x + y.
  unsigned int numA = num();
  unsigned int numY = y.num();
  unsigned int numX = x.num();
  unsigned int sizeX = x.size();
  unsigned int sizeY = y.size();

  assert(numX == numY && "private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

  if(numY == Siconos::DENSE && numX == Siconos::DENSE)
  {

    assert(y.dense() != x.dense());

    if(numA == Siconos::DENSE)
      noalias(*y.dense()) += prod(ublas::subrange(*dense(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else if(numA == Siconos::TRIANGULAR)
      noalias(*y.dense()) += prod(ublas::subrange(*triang(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else if(numA == Siconos::SYMMETRIC)
      noalias(*y.dense()) += prod(ublas::subrange(*sym(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else if(numA == Siconos::SPARSE)
      noalias(*y.dense()) += prod(ublas::subrange(*sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else //if(numA==Siconos::BANDED)
      noalias(*y.dense()) += prod(ublas::subrange(*banded(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
  }
  else // x and y sparse
  {
    if(numA == Siconos::SPARSE)
      *y.sparse() += prod(ublas::subrange(*sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.sparse());
    else
      SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
  }
}
