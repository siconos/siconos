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

#include "SiconosConfig.h"

#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

//#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace lapack = boost::numeric::bindings::lapack;


#include "SiconosVector.hpp"
#include "cholesky.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"

using namespace Siconos;

void SimpleMatrix::PLUFactorizationInPlace()
{
  if (_isPLUFactorized)
  {
    std::cout << "SimpleMatrix::PLUFactorizationInPlace warning: this matrix is already PLUFactorized. " << std::endl;
    return;
  }
  if (_num == 1)
  {
    if (!_ipiv)
      _ipiv.reset(new VInt(size(0)));
    else
      _ipiv->resize(size(0));
    int info = lapack::getrf(*mat.Dense, *_ipiv);
    if (info != 0)
    {
      _isPLUFactorized = false;
      SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace failed: the matrix is singular.");
    }
    else _isPLUFactorized = true;
  }
  else
  {
    int info = cholesky_decompose(*sparse());
    // \warning: VA 24/11/2010: work only for symmetric matrices. Should be replaced by efficient implementatation (e.g. mumps )
    if (info != 0)
    {
      display();
      _isPLUFactorized = false;
      std::cout << "Problem in Cholesky Decomposition for the row number" << info   << std::endl;
      SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace failed. ");
    }
    else _isPLUFactorized = true;
  }

}

void SimpleMatrix::PLUInverseInPlace()
{
  if (!_isPLUFactorized)
    PLUFactorizationInPlace();
  if (_num != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices.");

#if defined(HAVE_ATLAS) && defined(OUTSIDE_FRAMEWORK_BLAS)
  int info = lapack::getri(*mat.Dense, *_ipiv);   // solve from factorization

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace failed, the matrix is singular.");

  _isPLUInversed = true;
#else
  SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace not implemented with lapack.");
#endif
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(B) failed at solving Ax = B. Not yet implemented for a BlockMatrix B.");
  int info = 0;

  if (_num == 1)
  {
    if (!_isPLUFactorized) // call gesv => LU-factorize+solve
    {
      // solve system:
      if (!_ipiv)
        _ipiv.reset(new VInt(size(0)));
      else
        _ipiv->resize(size(0));
      info = lapack::gesv(*mat.Dense, *_ipiv, *(B.dense()));
      _isPLUFactorized = true;

      /*
        ublas::vector<double> S(std::max(size(0),size(1)));
        ublas::matrix<double, ublas::column_major> U(size(0),size(1));
        ublas::matrix<double, ublas::column_major> VT(size(0),size(1));

        int ierr = lapack::gesdd(*mat.Dense, S, U, VT);
        printf("info = %d, ierr = %d, emax = %f, emin = %f , cond = %f\n",info,ierr,S(0),S(2),S(0)/S(2));
      */
      // B now contains solution:
    }
    else // call getrs: only solve using previous lu-factorization
      if (B.num() == 1)
        info = lapack::getrs(*mat.Dense, *_ipiv, *(B.dense()));
      else
        SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices in RHS.");
  }
  else
  {
    if (!_isPLUFactorized) // call first PLUFactorizationInPlace
    {
      PLUFactorizationInPlace();
    }
    // and then solve
    if (B.num() == 1)
    {
      std::cout << "Version Sparse-Dense" << std::endl;
      inplace_solve(*sparse(), *(B.dense()), ublas::lower_tag());
      inplace_solve(ublas::trans(*sparse()), *(B.dense()), ublas::upper_tag());
      std::cout << "End Version Sparse-Dense" << std::endl;
    }
    else if (B.num() == 4)
    {
      std::cout << "Start Version Sparse-SParse" << std::endl;
      inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
      inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());
      std::cout << "End Version Sparse-Sparse" << std::endl;
    }
    else
      SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense ans sparse matrices in RHS.");
    info = 0 ;
  }
  //  SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices.");



  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(V) failed. Not yet implemented for V being a BlockVector.");


  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.dense()); // Conversion of vector to matrix. Temporary solution.
  int info;

  if (_num == 1)
  {
    if (!_isPLUFactorized) // call gesv => LU-factorize+solve
    {
      // solve system:
      if (!_ipiv)
        _ipiv.reset(new VInt(size(0)));
      else
        _ipiv->resize(size(0));

      info = lapack::gesv(*mat.Dense, *_ipiv, tmpB);
      _isPLUFactorized = true;

      /*
        ublas::matrix<double> COPY(*mat.Dense);
        ublas::vector<double> S(std::max(size(0),size(1)));
        ublas::matrix<double, ublas::column_major> U(size(0),size(1));
        ublas::matrix<double, ublas::column_major> VT(size(0),size(1));

        int ierr = lapack::gesdd(COPY, S, U, VT);
        printf("info = %d, ierr = %d, emax = %f, emin = %f , cond = %f\n",info,ierr,S(0),S(2),S(0)/S(2));
      */
      // B now contains solution:
    }
    else // call getrs: only solve using previous lu-factorization
      info = lapack::getrs(*mat.Dense, *_ipiv, tmpB);
  }
  else
  {
    if (!_isPLUFactorized) // call first PLUFactorizationInPlace
    {
      PLUFactorizationInPlace();
    }
    // and then solve
    inplace_solve(*sparse(), tmpB, ublas::lower_tag());
    inplace_solve(ublas::trans(*sparse()), tmpB, ublas::upper_tag());
    info = 0;
  }
  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
  else
  {
    noalias(*(B.dense())) = ublas::column(tmpB, 0);
  }
}

void SimpleMatrix::resetLU()
{
  if (_ipiv) _ipiv->clear();
  _isPLUFactorized = false;
  _isPLUInversed = false;
}

void SimpleMatrix::resetQR()
{
  _isQRFactorized = false;

}
// const SimpleMatrix operator * (const SimpleMatrix & A, const SimpleMatrix& B )
// {
//   return (DenseMat)prod(*A.dense() , *B.dense());
//   //  return A;
// }

void SimpleMatrix::SolveByLeastSquares(SiconosMatrix &B)
{
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::SolveByLeastSquares(Siconos Matrix &B) failed. Not yet implemented for M being a BlockMatrix.");
  int info = 0;
#ifdef USE_OPTIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, *(B.dense()), lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, *(B.dense()), lapack::minimal_workspace());
#endif
  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::SolveByLeastSquares failed.");
}




void SimpleMatrix::SolveByLeastSquares(SiconosVector &B)
{
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::SolveByLeastSquares(SiconosVector &B) failed. Not yet implemented for V being a BlockVector.");

  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.dense()); // Conversion of vector to matrix. Temporary solution.
  int info = 0;

#ifdef USE_OPTIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, tmpB, lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, tmpB, lapack::minimal_workspace());
#endif
  if (info != 0)
  {
    std::cout << "info = " << info << std::endl;
    SiconosMatrixException::selfThrow("SimpleMatrix::SolveByLeastSquares failed.");
  }
  else
  {
    noalias(*(B.dense())) = ublas::column(tmpB, 0);
  }

}

/*
void polePlacement(const SiconosMatrix& A, const SiconosVector& B, SiconosVector& P, bool transpose)
{
  unsigned int n = A.size(0);
  DenseMat AA(n, n);
  DenseMat Q(n, n);
  DenseVect tau(n);
  DenseVect BB(n);
  noalias(AA) = (*A.dense());
  lapack::gehrd(1, n, AA, tau);
  lapack::orghr(n, 1, n, Q, tau);
  noalias(BB) = prod(Q, *B.dense());
}
*/
