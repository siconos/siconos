/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "SiconosAlgebra.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "BlockMatrix.hpp"
#include "SparseMatrix.h"

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

  if ((!m1.isBlock() || !m2.isBlock()) && (m1.size(0) == m2.size(0)) && (m1.size(1) == m2.size(1)))
    return true;

  const SP::Index I1R = m1.tabRow();
  const SP::Index I2R = m2.tabRow();
  const SP::Index I1C = m1.tabCol();
  const SP::Index I2C = m2.tabCol();

  return ((*I1R == *I2R) && (*I1C == *I2C));
}

SiconosMatrix& operator *=(SiconosMatrix& m, const double& s)
{
  if (m._num == 0) // BlockMatrix
  {
    BlockMatrix& mB = static_cast<BlockMatrix&>(m);
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    for (it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
        (**it2) *= s;
    }
  }
  else if (m._num == 1)
    *m.dense() *= s;
  else if (m._num == 2)
    *m.triang() *= s;
  else if (m._num == 3)
    *m.sym() *= s;
  else if (m._num == 4)
    *m.sparse() *= s;
  else if (m._num == 5)
    *m.banded() *= s;
  else if (m._num == 6) {} // nothing!
  else //if(_num == 7)
    SiconosMatrixException::selfThrow(" SP::SiconosMatrix = (double) : invalid type of matrix");

  return m;
}

SiconosMatrix& operator /=(SiconosMatrix& m, const double& s)
{
  if (m._num == 0) // BlockMatrix
  {
    BlockMatrix& mB = static_cast<BlockMatrix&>(m);
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    for (it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
        (**it2) /= s;
    }
  }
  else if (m._num == 1)
    *m.dense() /= s;
  else if (m._num == 2)
    *m.triang() /= s;
  else if (m._num == 3)
    *m.sym() /= s;
  else if (m._num == 4)
    *m.sparse() /= s;
  else if (m._num == 5)
    *m.banded() /= s;
  else if (m._num == 6) {} // nothing!
  else //if(_num == 7)
    SiconosMatrixException::selfThrow(" SiconosMatrix *= (double) : invalid type of matrix");

  return m;
}

size_t SiconosMatrix::nnz(double tol)
{
  size_t nnz = 0;
  if (_num == 1) //dense
  {
    double* arr = getArray();
    for (size_t i = 0; i < size(0)*size(1); ++i)
    {
      if (fabs(arr[i]) > tol) { nnz++; }
    }
  }
  else if (_num == 4)
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
  csi* Mi = csc->i; // row indx
  csi* Mp = csc->p; // column pointers

  assert(Mp[col_off] >= 0);
  size_t nz = csc->p[col_off];

  size_t nrow = size(0);
  size_t ncol = size(1);

  csi pval = Mp[col_off];

  if (_num == 1) //dense
  {
    double* arr = getArray();
    for (size_t j = 0, joff = col_off; j < ncol; ++j)
    {
      for (size_t i = 0; i < nrow; ++i)
      {
        // col-major
        double elt_val = arr[i + j*nrow];
        // std::cout << " a(i=" << i << ",j=" << j << ") = "<< elt_val << std::endl;
        if (fabs(elt_val) > tol)
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
  else if (_num == 4)
  {
    const Index& ptr = sparse()->index1_data();
    const Index& indx = sparse()->index2_data();
    const ublas::unbounded_array<double>& vals = sparse()->value_data();

    size_t nnz =  sparse()->nnz();

    assert(ptr.size() == ncol + 1);
    assert(indx.size() == nnz);
    assert(vals.size() == nnz);

    for (size_t i = 0; i < nnz; ++i)
    {
      Mx[pval] = vals[i];
      Mi[pval++] = row_off + indx[i];
    }
    for (size_t j = 1, joff = col_off + 1; j < ncol+1; ++j, ++joff)
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

bool SiconosMatrix::fillTriplet(CSparseMatrix* triplet, size_t row_off, size_t col_off, double tol)
{
  assert(triplet);
  size_t nrow = size(0);
  size_t ncol = size(1);

  if (_num == 1) //dense
  {
    double* arr = getArray();
    for (size_t j = 0; j < ncol; ++j)
    {
      for (size_t i = 0; i < nrow; ++i)
      {
        // col-major

        cs_zentry(triplet, i + row_off, j + col_off, arr[i + j*nrow] );
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
