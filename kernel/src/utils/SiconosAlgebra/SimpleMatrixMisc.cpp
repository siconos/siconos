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


#include "determinant.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"

using namespace Siconos;

//=======================
//       get norm
//=======================

double SimpleMatrix::normInf() const
{
  if (_num == 1)
    return norm_inf(*mat.Dense);
  else if (_num == 2)
    return norm_inf(*mat.Triang);
  else if (_num == 3)
    return norm_inf(*mat.Sym);
  else if (_num == 4)
    return norm_inf(*mat.Sparse);
  else if (_num == 5)
    return norm_inf(*mat.Banded);
  else if (_num == 6)
    return 0;
  else // if(_num==7)
    return 1;
}

void SimpleMatrix::normInfByColumn(SP::SiconosVector vIn) const
{
  if (_num == 1)
  {
    if (vIn->size() != size(1))
      RuntimeException::selfThrow("SimpleMatrix::normInfByColumn: the given vector does not have the right length");
    DenseVect tmpV = DenseVect(size(0));
    for (unsigned int i = 0; i < size(1); i++)
    {
       ublas::noalias(tmpV) = ublas::column(*mat.Dense, i);
       (*vIn)(i) = norm_inf(tmpV);
    }
  }
  else
    RuntimeException::selfThrow("SimpleMatrix::normInfByColumn: not implemented for data other than DenseMat");
}
//=======================
//       determinant
//=======================

double SimpleMatrix::det() const
{
  if (_num == 1)
    return determinant(*mat.Dense);
  else if (_num == 2)
    return determinant(*mat.Triang);
  else if (_num == 3)
    return determinant(*mat.Sym);
  else if (_num == 4)
    return determinant(*mat.Sparse);
  else if (_num == 5)
    return determinant(*mat.Banded);
  else if (_num == 6)
    return 0;
  else // if(_num==7)
    return 1;
}


void SimpleMatrix::trans()
{
  switch (_num)
  {
  case 1:
    *mat.Dense = ublas::trans(*mat.Dense);
    break;
  case 2:
    SiconosMatrixException::selfThrow("SimpleMatrix::trans() failed, the matrix is triangular matrix and can not be transposed in place.");
    break;
  case 3:
    break;
  case 4:
    *mat.Sparse = ublas::trans(*mat.Sparse);
    break;
  case 5:
    *mat.Banded = ublas::trans(*mat.Banded);
    break;
  case 6:
    break;
  case 7:
    break;
  }
  resetLU();
}

void SimpleMatrix::trans(const SiconosMatrix &m)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, not yet implemented for m being a BlockMatrix.");


  if (&m == this)
    trans();//SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, m = this, use this->trans().");
  else
  {
    unsigned int numM = m.num();
    switch (numM)
    {
    case 1:
      if (_num != 1)
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a dense matrix into another type.");
      noalias(*mat.Dense) = ublas::trans(*m.dense());
      break;
    case 2:
      if (_num != 1)
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a triangular matrix into a non-dense one.");
      noalias(*mat.Dense) = ublas::trans(*m.triang());
      break;
    case 3:
      *this = m;
      break;
    case 4:
      if (_num == 1)
        noalias(*mat.Dense) = ublas::trans(*m.sparse());
      else if (_num == 4)
        noalias(*mat.Sparse) = ublas::trans(*m.sparse());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a sparse matrix into a forbidden type (not dense nor sparse).");
      break;
    case 5:
      if (_num == 1)
        noalias(*mat.Dense) = ublas::trans(*m.banded());
      else if (_num == 5)
        noalias(*mat.Banded) = ublas::trans(*m.banded());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a banded matrix into a forbidden type (not dense nor banded).");
      break;
    case 6:
      *this = m;
      break;
    case 7:
      *this = m;
    }
    // unsigned int tmp = _dimRow;
    // _dimRow = _dimCol;
    // _dimCol = tmp;
    resetLU();
  }
}







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

