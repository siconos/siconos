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

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "BlockMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "SiconosAlgebra.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "determinant.hpp"

using namespace siconos;

//=======================
//       get norm
//=======================

double SimpleMatrix::normInf() const {
  if (_num == DENSE)
    return norm_inf(*mat.Dense);
  else if (_num == TRIANGULAR)
    return norm_inf(*mat.Triang);
  else if (_num == SYMMETRIC)
    return norm_inf(*mat.Sym);
  else if (_num == SPARSE)
    return norm_inf(*mat.Sparse);
  else if (_num == SPARSE_COORDINATE)
    return norm_inf(*mat.SparseCoordinate);
  else if (_num == BANDED)
    return norm_inf(*mat.Banded);
  else if (_num == ZERO)
    return 0;
  else if (_num == IDENTITY)
    return 1;

  THROW_EXCEPTION("Matrix type not supported");
  return std::numeric_limits<double>::infinity();
}

void SimpleMatrix::normInfByColumn(SP::SiconosVector vIn) const {
  if (_num == DENSE) {
    if (vIn->size() != size(1))
      THROW_EXCEPTION("the given vector does not have the right length");
    DenseVect tmpV = DenseVect(size(0));
    for (unsigned int i = 0; i < size(1); i++) {
      ublas::noalias(tmpV) = ublas::column(*mat.Dense, i);
      (*vIn)(i) = norm_inf(tmpV);
    }
  } else
    THROW_EXCEPTION("not implemented for data other than DenseMat");
}
//=======================
//       determinant
//=======================

double SimpleMatrix::det() const {
  if (_num == DENSE)
    return siconos::externals::ublas::determinant(*mat.Dense);
  else if (_num == TRIANGULAR)
    return siconos::externals::ublas::determinant(*mat.Triang);
  else if (_num == SYMMETRIC)
    return siconos::externals::ublas::determinant(*mat.Sym);
  else if (_num == SPARSE)
    return siconos::externals::ublas::determinant(*mat.Sparse);
  else if (_num == SPARSE_COORDINATE)
    return siconos::externals::ublas::determinant(*mat.Sparse);
  else if (_num == BANDED)
    return siconos::externals::ublas::determinant(*mat.Banded);
  else if (_num == ZERO)
    return 0;
  else if (_num == IDENTITY)
    return 1;
  THROW_EXCEPTION("Matrix type not supported");
  return std::numeric_limits<double>::infinity();
}

void SimpleMatrix::trans() {
  switch (_num) {
    case DENSE:
      *mat.Dense = ublas::trans(*mat.Dense);
      break;
    case TRIANGULAR:
      THROW_EXCEPTION(
          "failed, the matrix is triangular matrix and can not be transposed in place.");
      break;
    case SYMMETRIC:
      break;
    case SPARSE:
      *mat.Sparse = ublas::trans(*mat.Sparse);
      break;
    case SPARSE_COORDINATE:
      *mat.Sparse = ublas::trans(*mat.Sparse);
      break;
    case BANDED:
      *mat.Banded = ublas::trans(*mat.Banded);
      break;
    case siconos::ZERO:
      break;
    case siconos::IDENTITY:
      break;
    default:
      THROW_EXCEPTION("Matrix type not supported");
  }
  resetFactorizationFlags();
}

void SimpleMatrix::trans(const SiconosMatrix& m) {
  if (m.isBlock()) THROW_EXCEPTION("not yet implemented for m being a BlockMatrix.");

  if (&m == this)
    trans();
  else {
    siconos::UBLAS_TYPE numM = m.num();
    switch (numM) {
      case DENSE:
        if (_num != DENSE)
          THROW_EXCEPTION("try to transpose a dense matrix into another type.");
        noalias(*mat.Dense) = ublas::trans(*m.dense());
        break;
      case TRIANGULAR:
        if (_num != DENSE)
          THROW_EXCEPTION("try to transpose a triangular matrix into a non-dense one.");
        noalias(*mat.Dense) = ublas::trans(*m.triang());
        break;
      case SYMMETRIC:
        *this = m;
        break;
      case SPARSE:
        if (_num == DENSE)
          noalias(*mat.Dense) = ublas::trans(*m.sparse());
        else if (_num == SPARSE)
          noalias(*mat.Sparse) = ublas::trans(*m.sparse());
        else if (_num == SPARSE_COORDINATE)
          noalias(*mat.SparseCoordinate) = ublas::trans(*m.sparse());
        else
          THROW_EXCEPTION(
              "try to transpose a sparse matrix into a forbidden type (not dense nor "
              "sparse).");
        break;
      case SPARSE_COORDINATE:
        if (_num == DENSE)
          noalias(*mat.Dense) = ublas::trans(*m.sparseCoordinate());
        else if (_num == SPARSE)
          noalias(*mat.Sparse) = ublas::trans(*m.sparseCoordinate());
        else if (_num == SPARSE_COORDINATE)
          noalias(*mat.SparseCoordinate) = ublas::trans(*m.sparseCoordinate());
        else
          THROW_EXCEPTION(
              "try to transpose a sparse coordinate matrix into a forbidden type (not dense "
              "nor sparse coordinate).");
        break;
      case BANDED:
        if (_num == DENSE)
          noalias(*mat.Dense) = ublas::trans(*m.banded());
        else if (_num == BANDED)
          noalias(*mat.Banded) = ublas::trans(*m.banded());
        else
          THROW_EXCEPTION(
              "try to transpose a banded matrix into a forbidden type (not dense nor "
              "banded).");
        break;
      case ZERO:
        *this = m;
        break;
      case IDENTITY:
        *this = m;
        break;
      default:
        THROW_EXCEPTION("");
    }
    resetFactorizationFlags();
  }
}

/*
The following code inverts the matrix input using LU-decomposition with backsubstitution of
unit vectors. Reference: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling &
Flannery.

you can solve Ax=b using three lines of ublas code:

permutation_matrix<> piv;
lu_factorize(A, piv);
lu_substitute(A, piv, x);

*/
#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

// REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

/* Matrix inversion routine.
Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template <class T, class U, class V>
bool InvertMatrix(const ublas::matrix<T, U, V>& input, ublas::matrix<T, U, V>& inverse) {
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<T, U, V> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0) return false;

  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<T>(A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

#endif  // INVERT_MATRIX_HPP

// Note FP: never used. Comment before removal ?
// void invertMatrix(const SimpleMatrix& input, SimpleMatrix& output)
// {
//   InvertMatrix(*input.dense(), *output.dense());
// }
