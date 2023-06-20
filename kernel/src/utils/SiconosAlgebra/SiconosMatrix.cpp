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

#include "SiconosMatrix.hpp"

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "BlockMatrix.hpp"           // for BlockMatrix
#include "NumericsToolsNamespace.h"  // for CSparseMatrix,  NSM_fix_csc
#include "SiconosMatrixOp.hpp"       // for matrix operators declaration
#include "SiconosVector.hpp"         // for SiconosVector
#include "Tools.hpp"                 // toString

namespace ublas = boost::numeric::ublas;

siconos::algebra::SiconosMatrix& siconos::algebra::operator*=(SiconosMatrix& m,
                                                              const double& s) {
  if (m._num == UblasType::BLOCK)  // BlockMatrix
  {
    BlockMatrix& mB = static_cast<BlockMatrix&>(m);
    for (auto it = mB._mat->begin1(); it != mB._mat->end1(); ++it) {
      for (auto it2 = it.begin(); it2 != it.end(); ++it2) (**it2) *= s;
    }
  } else if (m._num == UblasType::DENSE)
    *m.dense() *= s;
  else if (m._num == UblasType::TRIANGULAR)
    *m.triang() *= s;
  else if (m._num == UblasType::SYMMETRIC)
    *m.sym() *= s;
  else if (m._num == UblasType::SPARSE)
    *m.sparse() *= s;
  else if (m._num == UblasType::SPARSE_COORDINATE)
    *m.sparseCoordinate() *= s;
  else if (m._num == UblasType::BANDED)
    *m.banded() *= s;
  else if (m._num == UblasType::ZERO) {
  }     // nothing!
  else  // if(_num == 7)
    THROW_EXCEPTION("invalid type of matrix");

  return m;
}

siconos::algebra::SiconosMatrix& siconos::algebra::operator/=(SiconosMatrix& m,
                                                              const double& s) {
  if (m._num == UblasType::BLOCK)  // BlockMatrix
  {
    BlockMatrix& mB = static_cast<BlockMatrix&>(m);
    for (auto it = mB._mat->begin1(); it != mB._mat->end1(); ++it) {
      for (auto it2 = it.begin(); it2 != it.end(); ++it2) (**it2) /= s;
    }
  } else if (m._num == UblasType::DENSE)
    *m.dense() /= s;
  else if (m._num == UblasType::TRIANGULAR)
    *m.triang() /= s;
  else if (m._num == UblasType::SYMMETRIC)
    *m.sym() /= s;
  else if (m._num == UblasType::SPARSE)
    *m.sparse() /= s;
  else if (m._num == UblasType::SPARSE_COORDINATE)
    *m.sparseCoordinate() /= s;
  else if (m._num == UblasType::BANDED)
    *m.banded() /= s;
  else if (m._num == UblasType::ZERO) {
  }     // nothing!
  else  // if(_num == 7)
    THROW_EXCEPTION("invalid type of matrix");

  return m;
}

size_t siconos::algebra::SiconosMatrix::nnz(double tol) {
  size_t nnz = 0;
  if (_num == UblasType::DENSE)  // dense
  {
    double* arr = getArray();
    for (size_t i = 0; i < size(0) * size(1); ++i) {
      if (fabs(arr[i]) > tol) {
        nnz++;
      }
    }
  } else if (_num == UblasType::SPARSE) {
    nnz = sparse()->nnz();
  } else
    THROW_EXCEPTION("not implemented for the given matrix type");

  return nnz;
}

bool siconos::algebra::SiconosMatrix::fillCSC(CSparseMatrix* csc, size_t row_off,
                                              size_t col_off, double tol) {
  assert(csc);
  double* Mx = csc->x;  // data
  CS_INT* Mi = csc->i;  // row indx
  CS_INT* Mp = csc->p;  // column pointers

  assert(Mp[col_off] >= 0);
  size_t nz = csc->p[col_off];

  size_t nrow = size(0);
  size_t ncol = size(1);

  CS_INT pval = Mp[col_off];

  if (_num == UblasType::DENSE)  // dense
  {
    double* arr = getArray();
    for (size_t j = 0, joff = col_off; j < ncol; ++j) {
      for (size_t i = 0; i < nrow; ++i) {
        // col-major
        double elt_val = arr[i + j * nrow];
        // std::cout << " a(i=" << i << ",j=" << j << ") = "<< elt_val << std::endl;
        if (fabs(elt_val) > tol) {
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
  } else if (_num == UblasType::SPARSE) {
    const auto& ptr = sparse()->index1_data();
    const auto& indx = sparse()->index2_data();
    const ublas::unbounded_array<double>& vals = sparse()->value_data();

    size_t nnz = sparse()->nnz();

    assert(ptr.size() == ncol + 1);
    assert(indx.size() == nnz);
    assert(vals.size() == nnz);

    for (size_t i = 0; i < nnz; ++i) {
      Mx[pval] = vals[i];
      Mi[pval++] = row_off + indx[i];
    }
    for (size_t j = 1, joff = col_off + 1; j < ncol + 1; ++j, ++joff) {
      Mp[joff] = nz + ptr[j];
    }
  } else {
    THROW_EXCEPTION("not implemented for the given matrix type");
  }

  return true;
}

bool siconos::algebra::SiconosMatrix::fillCSC(CSparseMatrix* csc, double tol) {
  assert(csc);
  double* Mx = csc->x;  // data
  CS_INT* Mi = csc->i;  // row indx
  CS_INT* Mp = csc->p;  // column pointers

  size_t nrow = size(0);
  size_t ncol = size(1);

  CS_INT pval = 0;

  if (_num == UblasType::DENSE)  // dense
  {
    double* arr = getArray();
    for (size_t j = 0, joff = 0; j < ncol; ++j) {
      for (size_t i = 0; i < nrow; ++i) {
        // col-major
        double elt_val = arr[i + j * nrow];
        // std::cout << " a(i=" << i << ",j=" << j << ") = "<< elt_val << std::endl;
        if (fabs(elt_val) > tol) {
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
  } else if (_num == UblasType::SPARSE) {
    const auto& ptr = sparse()->index1_data();
    const auto& indx = sparse()->index2_data();
    const ublas::unbounded_array<double>& vals = sparse()->value_data();

    size_t nnz = sparse()->nnz();

    assert(ptr.size() == ncol + 1);
    assert(indx.size() >= nnz);
    assert(vals.size() >= nnz);

    for (size_t i = 0; i < nnz; ++i) {
      Mx[pval] = vals[i];
      Mi[pval++] = indx[i];
    }
    for (size_t j = 0; j < ncol + 1; ++j) {
      Mp[j] = ptr[j];
    }
  } else {
    THROW_EXCEPTION("not implemented for the given matrix type");
  }

  return true;
}

bool siconos::algebra::SiconosMatrix::fromCSC(CSparseMatrix* csc) {
  assert(csc);

  NSM_sort_csc(csc);

  double* Mx = csc->x;  // data
  CS_INT* Mi = csc->i;  // row indx
  CS_INT* Mp = csc->p;  // column pointers
  CS_INT n = csc->n;
  // CS_INT m = csc->m;

  // size_t nnz = csc->p[n];

  if (_num == UblasType::SPARSE) {
    sparse()->clear();
    CS_INT pval = 0;
    // push_back in order should be in constant time
    // http://www.guwi17.de/ublas/matrix_sparse_usage.html
    for (CS_INT col = 0; col < n; col++) {
      for (CS_INT p = Mp[col]; p < Mp[col + 1]; p++) {
        sparse()->push_back(Mi[pval], col, Mx[pval]);
        pval++;
      }
    }

    // not able to work directly on the contents of the sparse matrix.
    // sparse()->resize(m, nnz, false);
    // //sparse()->set_filled(nnz, n+1);
    // std::cout << sparse()->size1() << std::endl;
    // std::cout << sparse()->size2() << std::endl;

    // Index& ptr = sparse()->index1_data();
    // Index& indx = sparse()->index2_data();
    // ublas::unbounded_array<double>& vals = sparse()->value_data();

    // // ptr.resize(n+1);
    // // indx.resize(nnz);
    // // vals.resize(nnz);

    // std::cout << sparse()->filled1() << std::endl;
    // std::cout << sparse()->filled2() << std::endl;

    // // CS_INT pval=0;

    // // for(size_t i = 0; i < nnz; ++i)
    // // {
    // //   vals[i] = Mx[pval];
    // //   indx[i] = Mi[pval++];
    // // }
    // for(size_t i = 0; i < nnz; ++i)
    // {
    //   printf("vals[%i] = %e\t", i, vals[i]);
    //   printf("indx[%i] = %i\t", i, indx[i]);
    // }

    // // for(size_t j = 0 ; j < n+1; ++j)
    // // {
    // //   ptr[j] = Mp[j]  ;
    // // }
    // for(size_t j = 0 ; j < n+1; ++j) printf("ptr[%i] = %i\t", j, ptr[j]);
  } else if (_num == UblasType::DENSE) {
    CS_INT pval = 0;
    // push_back in order should be in constant time
    // http://www.guwi17.de/ublas/matrix_sparse_usage.html
    for (CS_INT col = 0; col < n; col++) {
      for (CS_INT p = Mp[col]; p < Mp[col + 1]; p++) {
        setValue(Mi[pval], col, Mx[pval]);
        pval++;
      }
    }
  } else {
    THROW_EXCEPTION("not implemented for the given matrix type");
  }
  return true;
}

bool siconos::algebra::SiconosMatrix::fillTriplet(CSparseMatrix* triplet, size_t row_off,
                                                  size_t col_off, double tol) {
  assert(triplet);
  size_t nrow = size(0);
  size_t ncol = size(1);

  if (_num == UblasType::DENSE)  // dense
  {
    double* arr = getArray();
    for (size_t j = 0; j < ncol; ++j) {
      for (size_t i = 0; i < nrow; ++i) {
        // col-major

        CSparseMatrix_zentry(triplet, i + row_off, j + col_off, arr[i + j * nrow],
                             std::numeric_limits<double>::epsilon());
      }
    }
  } else {
    THROW_EXCEPTION("not implemented for the given matrix type");
  }

  return true;
}

std::ostream& siconos::algebra::operator<<(std::ostream& os, const SiconosMatrix& sm) {
  os << siconos::tools::toString(sm);
  return os;
}

void siconos::algebra::SiconosMatrix::private_prod(unsigned int startRow,
                                                   const SiconosVector& x, SiconosVector& y,
                                                   bool init) const {
  assert(!(isFactorized()) && "A is Factorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of
  // index (row) startRow and startRow + sizeY

  if (init)  // y = subA * x , else y += subA * x
    y.zero();
  private_addprod(startRow, 0, x, y);
}

/** Computation of y = subA.x

    where subA is a sub-matrix of A, subA[0,0] = A[startRow, startCol]
 */
void siconos::algebra::SiconosMatrix::private_addprod(unsigned startRow, unsigned int startCol,
                                                      const SiconosVector& x,
                                                      SiconosVector& y) const {
  assert(!(isFactorized()) && "A is Factorized in prod !!");
  assert(!isBlock() &&
         "private_addprod(start,x,y) error: not yet implemented for block matrix.");

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and
  // between columns startCol and (startCol+sizeX). Then computation of y = subA*x + y.
  auto numA = num();
  auto numY = y.num();
  auto numX = x.num();
  auto sizeX = x.size();
  auto sizeY = y.size();

  assert(numX == numY &&
         "private_addprod(A,start,x,y) error: not yet implemented for x and y of different "
         "types.");

  if (numY == UblasType::DENSE && numX == UblasType::DENSE) {
    assert(y.dense() != x.dense());

    if (numA == UblasType::DENSE)
      noalias(*y.dense()) += prod(
          ublas::subrange(*dense(), startRow, startRow + sizeY, startCol, startCol + sizeX),
          *x.dense());
    else if (numA == UblasType::TRIANGULAR)
      noalias(*y.dense()) += prod(
          ublas::subrange(*triang(), startRow, startRow + sizeY, startCol, startCol + sizeX),
          *x.dense());
    else if (numA == UblasType::SYMMETRIC)
      noalias(*y.dense()) +=
          prod(ublas::subrange(*sym(), startRow, startRow + sizeY, startCol, startCol + sizeX),
               *x.dense());
    else if (numA == UblasType::SPARSE)
      noalias(*y.dense()) += prod(
          ublas::subrange(*sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX),
          *x.dense());
    else  // if(numA==UblasType::BANDED)
      noalias(*y.dense()) += prod(
          ublas::subrange(*banded(), startRow, startRow + sizeY, startCol, startCol + sizeX),
          *x.dense());
  } else  // x and y sparse
  {
    if (numA == UblasType::SPARSE)
      *y.sparse() += prod(
          ublas::subrange(*sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX),
          *x.sparse());
    else
      THROW_EXCEPTION("not yet implemented for x, y  sparse and A not sparse.");
  }
}
