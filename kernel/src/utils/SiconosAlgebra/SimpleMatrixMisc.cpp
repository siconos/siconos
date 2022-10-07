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
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "BlockMatrix.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "determinant.hpp"

namespace ublas = boost::numeric::ublas;

//=======================
//       get norm
//=======================

double siconos::algebra::SimpleMatrix::normInf() const
{
  if (_num == UblasType::DENSE)
    return norm_inf(*mat.Dense);
  else if (_num == UblasType::TRIANGULAR)
    return norm_inf(*mat.Triang);
  else if (_num == UblasType::SYMMETRIC)
    return norm_inf(*mat.Sym);
  else if (_num == UblasType::SPARSE)
    return norm_inf(*mat.Sparse);
  else if (_num == UblasType::SPARSE_COORDINATE)
    return norm_inf(*mat.SparseCoordinate);
  else if (_num == UblasType::BANDED)
    return norm_inf(*mat.Banded);
  else if (_num == UblasType::ZERO)
    return 0;
  else if (_num == UblasType::IDENTITY)
    return 1;

  THROW_EXCEPTION("Matrix type not supported");
  return std::numeric_limits<double>::infinity();
}

void siconos::algebra::SimpleMatrix::normInfByColumn(std::shared_ptr<SiconosVector> vIn) const
{
  if (_num == UblasType::DENSE) {
    if (vIn->size() != size(1))
      THROW_EXCEPTION("the given vector does not have the right length");
    DenseVect tmpV = DenseVect(size(0));
    for (unsigned int i = 0; i < size(1); i++) {
      ublas::noalias(tmpV) = ublas::column(*mat.Dense, i);
      (*vIn)(i) = norm_inf(tmpV);
    }
  }
  else
    THROW_EXCEPTION("not implemented for data other than DenseMat");
}
//=======================
//       siconos::externals::ublas::determinant
//=======================

double siconos::algebra::SimpleMatrix::det() const
{
  if (_num == UblasType::DENSE)
    return siconos::externals::ublas::determinant(*mat.Dense);
  else if (_num == UblasType::TRIANGULAR)
    return siconos::externals::ublas::determinant(*mat.Triang);
  else if (_num == UblasType::SYMMETRIC)
    return siconos::externals::ublas::determinant(*mat.Sym);
  else if (_num == UblasType::SPARSE)
    return siconos::externals::ublas::determinant(*mat.Sparse);
  else if (_num == UblasType::SPARSE_COORDINATE)
    return siconos::externals::ublas::determinant(*mat.Sparse);
  else if (_num == UblasType::BANDED)
    return siconos::externals::ublas::determinant(*mat.Banded);
  else if (_num == UblasType::ZERO)
    return 0;
  else if (_num == UblasType::IDENTITY)
    return 1;
  THROW_EXCEPTION("Matrix type not supported");
  return std::numeric_limits<double>::infinity();
}

void siconos::algebra::SimpleMatrix::trans()
{
  switch (_num) {
    case UblasType::DENSE:
      *mat.Dense = ublas::trans(*mat.Dense);
      break;
    case UblasType::TRIANGULAR:
      THROW_EXCEPTION(
          "failed, the matrix is triangular matrix and can not be transposed in place.");
      break;
    case UblasType::SYMMETRIC:
      break;
    case UblasType::SPARSE:
      *mat.Sparse = ublas::trans(*mat.Sparse);
      break;
    case UblasType::SPARSE_COORDINATE:
      *mat.Sparse = ublas::trans(*mat.Sparse);
      break;
    case UblasType::BANDED:
      *mat.Banded = ublas::trans(*mat.Banded);
      break;
    case UblasType::ZERO:
      break;
    case UblasType::IDENTITY:
      break;
    default:
      THROW_EXCEPTION("Matrix type not supported");
  }
  resetFactorizationFlags();
}

void siconos::algebra::SimpleMatrix::trans(const SiconosMatrix &m)
{
  if (m.isBlock()) THROW_EXCEPTION("not yet implemented for m being a BlockMatrix.");

  if (&m == this)
    trans();
  else {
    auto numM = m.num();
    switch (numM) {
      case UblasType::DENSE:
        if (_num != UblasType::DENSE)
          THROW_EXCEPTION("try to transpose a dense matrix into another type.");
        noalias(*mat.Dense) = ublas::trans(*m.dense());
        break;
      case UblasType::TRIANGULAR:
        if (_num != UblasType::DENSE)
          THROW_EXCEPTION("try to transpose a triangular matrix into a non-dense one.");
        noalias(*mat.Dense) = ublas::trans(*m.triang());
        break;
      case UblasType::SYMMETRIC:
        *this = m;
        break;
      case UblasType::SPARSE:
        if (_num == UblasType::DENSE)
          noalias(*mat.Dense) = ublas::trans(*m.sparse());
        else if (_num == UblasType::SPARSE)
          noalias(*mat.Sparse) = ublas::trans(*m.sparse());
        else if (_num == UblasType::SPARSE_COORDINATE)
          noalias(*mat.SparseCoordinate) = ublas::trans(*m.sparse());
        else
          THROW_EXCEPTION(
              "try to transpose a sparse matrix into a forbidden type (not dense nor "
              "sparse).");
        break;
      case UblasType::SPARSE_COORDINATE:
        if (_num == UblasType::DENSE)
          noalias(*mat.Dense) = ublas::trans(*m.sparseCoordinate());
        else if (_num == UblasType::SPARSE)
          noalias(*mat.Sparse) = ublas::trans(*m.sparseCoordinate());
        else if (_num == UblasType::SPARSE_COORDINATE)
          noalias(*mat.SparseCoordinate) = ublas::trans(*m.sparseCoordinate());
        else
          THROW_EXCEPTION(
              "try to transpose a sparse coordinate matrix into a forbidden type "
              "(not dense nor sparse coordinate).");
        break;
      case UblasType::BANDED:
        if (_num == UblasType::DENSE)
          noalias(*mat.Dense) = ublas::trans(*m.banded());
        else if (_num == UblasType::BANDED)
          noalias(*mat.Banded) = ublas::trans(*m.banded());
        else
          THROW_EXCEPTION(
              "try to transpose a banded matrix into a forbidden type (not dense nor "
              "banded).");
        break;
      case UblasType::ZERO:
        *this = m;
        break;
      case UblasType::IDENTITY:
        *this = m;
        break;
      default:
        THROW_EXCEPTION("");
    }
    resetFactorizationFlags();
  }
}
