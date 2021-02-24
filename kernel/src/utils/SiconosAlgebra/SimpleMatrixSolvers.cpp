/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"

//#define DEBUG_MESSAGES
#include "siconos_debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#include <cs.h>
#endif

using namespace Siconos;

void SimpleMatrix::PLUFactorizationInPlace()
{
  if (_isPLUFactorized)
  {
    std::cout << "SimpleMatrix::PLUFactorizationInPlace warning: this matrix is already PLUFactorized. " << std::endl;
    return;
  }
  if (_num == Siconos::DENSE)
  {
    if (!_ipiv)
      _ipiv.reset(new VInt(size(0)));
    else
      _ipiv->resize(size(0));
    int info = lapack::getrf(*mat.Dense, *_ipiv);
    if (info != 0)
    {
      _isPLUFactorized = false;
      _isPLUFactorizedInPlace = true;
      THROW_EXCEPTION("SimpleMatrix::PLUFactorizationInPlace failed: the matrix is singular.");
    }
    else
    {
      _isPLUFactorized = true;
      _isPLUFactorizedInPlace = true;
    }
  }
  else
  {
    int info = cholesky_decompose(*sparse());
    // \warning: VA 24/11/2010: work only for symmetric matrices. Should be replaced by efficient implementatation (e.g. mumps )
    if (info != 0)
    {
      display();
      _isPLUFactorized = false;
      _isPLUFactorizedInPlace = false;
      std::cout << "Problem in Cholesky Decomposition for the row number" << info   << std::endl;
      THROW_EXCEPTION("SimpleMatrix::PLUFactorizationInPlace failed. ");
    }
    else
    {
      _isPLUFactorized = true;
      _isPLUFactorizedInPlace = false;
    }
  }

}



void SimpleMatrix::PLUInverseInPlace()
{
  if(!_isPLUFactorized)
    PLUFactorizationInPlace();
  if(_num != Siconos::DENSE)
    THROW_EXCEPTION(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices.");

#if defined(HAS_LAPACK_dgetri)
  int info = lapack::getri(*mat.Dense, *_ipiv);   // solve from factorization

  if(info != 0)
    THROW_EXCEPTION("SimpleMatrix::PLUInverseInPlace failed, the matrix is singular.");

  _isPLUInversed = true;
#else
  THROW_EXCEPTION("SimpleMatrix::PLUInverseInPlace not implemented with lapack.");
#endif
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  if(B.isBlock())
    THROW_EXCEPTION("SimpleMatrix PLUForwardBackwardInPlace(B) failed at solving Ax = B. Not yet implemented for a BlockMatrix B.");
  int info = 0;

  if(_num == Siconos::DENSE)
  {
    if(!_isPLUFactorized)  // call gesv => LU-factorize+solve
    {
      // solve system:
      if(!_ipiv)
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
      if(B.num() == DENSE)
        info = lapack::getrs(*mat.Dense, *_ipiv, *(B.dense()));
      else
        THROW_EXCEPTION(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices in RHS.");

  }
  else
  {
    if(!_isPLUFactorized)  // call first PLUFactorizationInPlace
    {
      PLUFactorizationInPlace();
    }
    // and then solve
    if(B.num() == Siconos::DENSE)
    {
      inplace_solve(*sparse(), *(B.dense()), ublas::lower_tag());
      inplace_solve(ublas::trans(*sparse()), *(B.dense()), ublas::upper_tag());
    }
    else if(B.num() == Siconos::SPARSE)
    {
      inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
      inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());
    }
    else
      THROW_EXCEPTION(" SimpleMatrix::PLUInverseInPlace: only implemented for dense ans sparse matrices in RHS.");
    info = 0 ;
  }
  //  THROW_EXCEPTION(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices.");



  if(info != 0)
    THROW_EXCEPTION("SimpleMatrix::PLUForwardBackwardInPlace failed.");
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.dense()); // Conversion of vector to matrix. Temporary solution.
  int info;

  if(_num == Siconos::DENSE)
  {
    if(!_isPLUFactorized)  // call gesv => LU-factorize+solve
    {
      // solve system:
      if(!_ipiv)
        _ipiv.reset(new VInt(size(0)));
      else
        _ipiv->resize(size(0));

      info = lapack::gesv(*mat.Dense, *_ipiv, tmpB);
      _isPLUFactorized = true;

      /*
        ublas::matrix<double> COPY(*mat.Dense);
        ublas::vector<double> S(std::max(size(0),size(1)));
        ublas::matrix<double, ublas::column_major> U(size(0),size(1));
`        ublas::matrix<double, ublas::column_major> VT(size(0),size(1));

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
    if(!_isPLUFactorized)  // call first PLUFactorizationInPlace
    {
      PLUFactorizationInPlace();
    }
    // and then solve
    inplace_solve(*sparse(), tmpB, ublas::lower_tag());
    inplace_solve(ublas::trans(*sparse()), tmpB, ublas::upper_tag());
    info = 0;
  }
  if(info != 0)
    THROW_EXCEPTION("SimpleMatrix::PLUForwardBackwardInPlace failed.");

  noalias(*(B.dense())) = ublas::column(tmpB, 0);

}

void SimpleMatrix::resetLU()
{
  if(_ipiv) _ipiv->clear();
  _isPLUFactorized = false;
  _isPLUFactorizedInPlace = false;
  _isPLUInversed = false;
}

void SimpleMatrix::resetCholesky()
{
  _isCholeskyFactorized = false;
  _isCholeskyFactorizedInPlace = false;
}

void SimpleMatrix::resetQR()
{
  _isQRFactorized = false;

}

void SimpleMatrix::resetFactorizationFlags()
{
  resetLU();
  resetCholesky();
  resetQR();
}


// const SimpleMatrix operator * (const SimpleMatrix & A, const SimpleMatrix& B )
// {
//   return (DenseMat)prod(*A.dense() , *B.dense());
//   //  return A;
// }

void SimpleMatrix::SolveByLeastSquares(SiconosMatrix &B)
{
  if(B.isBlock())
    THROW_EXCEPTION("SimpleMatrix::SolveByLeastSquares(Siconos Matrix &B) failed. Not yet implemented for M being a BlockMatrix.");
  int info = 0;
#ifdef USE_OPTIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, *(B.dense()), lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, *(B.dense()), lapack::minimal_workspace());
#endif
  if(info != 0)
    THROW_EXCEPTION("SimpleMatrix::SolveByLeastSquares failed.");
}




void SimpleMatrix::SolveByLeastSquares(SiconosVector &B)
{
  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.dense()); // Conversion of vector to matrix. Temporary solution.
  int info = 0;

#ifdef USE_OPTIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, tmpB, lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, tmpB, lapack::minimal_workspace());
#endif
  if(info != 0)
  {
    std::cout << "info = " << info << std::endl;
    THROW_EXCEPTION("SimpleMatrix::SolveByLeastSquares failed.");
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


void SimpleMatrix::Factorize()
{
  DEBUG_BEGIN("void SimpleMatrix::Factorize()\n");
  if(isFactorized())
  {
    std::cout << "SimpleMatrix::PLUFactorize warning: this matrix is already Factorized. " << std::endl;
    return;
  }

  /* set the numericsMatrix */

  updateNumericsMatrix();
  NumericsMatrix * NM = numericsMatrix();

  /* Factorization calling the right method in Numerics */
  int info =1;
  if (isSymmetric())
  {
    if (isPositiveDefinite()) // Cholesky Factorization
    {
      //std::cout << "Cholesky Factorize"<< std::endl;
      info  = NM_Cholesky_factorize(NM);

      if(info != 0)
      {
        _isCholeskyFactorized = false;
        THROW_EXCEPTION("SimpleMatrix::Factorize failed (Cholesky)");
      }
      else
      {
        _isCholeskyFactorized = true;
      }
    }
    else  //  LDLT Factorization
    {
      THROW_EXCEPTION("SimpleMatrix::Factorize failed: LDL^T not yet implemented.");
    }
  }
  else //  LU Factorization  by default
  {
    info  = NM_LU_factorize(NM);

    if(info != 0)
    {
      _isPLUFactorized = false;
      if(_num == DENSE)
        _isPLUFactorizedInPlace = true;
      THROW_EXCEPTION("SimpleMatrix::PLUFactorize failed: the matrix is singular.");
    }
    else
    {
      _isPLUFactorized = true;
       if(_num == DENSE)
         _isPLUFactorizedInPlace = true;
    }
  }
  DEBUG_END("void SimpleMatrix::Factorize()\n");
}

//#define SPARSE_RHS_COPY_TO_DENSE 1

void SimpleMatrix::Solve(SiconosMatrix &B)
{
  if(B.isBlock())
    THROW_EXCEPTION("SimpleMatrix Solve(B) failed at solving Ax = B. Not yet implemented for a BlockMatrix B.");

  int info = 1;
  if(!isFactorized())
  {
    Factorize();
  }
  // and then solve
  NumericsMatrix * NM = _numericsMatrix.get();


#ifdef SPARSE_RHS_COPY_TO_DENSE

  double * b;
  SP::SimpleMatrix Bdense;

  if(B.num() == DENSE)
  {
    b = B.getArray();
  }
  else if(B.num() == SPARSE)
  {
    // First way.
    // We copy to dense since our sparse solver is not able to take
    // into account for sparse r.h.s, yet.
    Bdense.reset (new SimpleMatrix(size(0),size(1)));
    * Bdense= B ;                                                // copy to dense
    b = &(*Bdense->getArray());

    // Second way
    // use inplace_solve of ublas (see above with SolveInPlace)
    // For that, we need to fill our factorization given by NM_LU_factorize
    // into a ublas sparse matrix
    // inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
    // inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());

  }
  else
    THROW_EXCEPTION(" SimpleMatrix::Solve: only implemented for dense and sparse matrices in RHS.");

#else
  B.updateNumericsMatrix();
  NumericsMatrix * NM_B = B.numericsMatrix();
  //NM_display(NM_B);
#endif



  if (isSymmetric())
  {
    if (isPositiveDefinite()) // Cholesky Solving
    {
      //std::cout << "Cholesky Solve"<< std::endl;
#ifdef SPARSE_RHS_COPY_TO_DENSE
      info  = NM_Cholesky_solve(NM, b, B.size(1));
#else
      info  = NM_Cholesky_solve_matrix_rhs(NM, NM_B);
#endif
      if(info != 0)
      {
        THROW_EXCEPTION("SimpleMatrix::Solve failed (Cholesky)");
      }
    }
    else  //  LDLT Factorization
    {
      THROW_EXCEPTION("SimpleMatrix::Solve failed: LDL^T not yet implemented.");
    }
  }
  else //  LU Factorization  by default
  {
#ifdef SPARSE_RHS_COPY_TO_DENSE
    info  = NM_LU_solve(NM, b, B.size(1));
#else
    info  = NM_LU_solve_matrix_rhs(NM, NM_B);
#endif
    if(info != 0)
    {
      THROW_EXCEPTION("SimpleMatrix::PLUFactorize failed: the matrix is singular.");
    }
  }

  if(B.num() == SPARSE)
  {
#ifdef SPARSE_RHS_COPY_TO_DENSE
    B = *Bdense ; // we copy back to sparse.
    //B.displayExpert();
#else
    /* we need to fill back again */
    //B.fromCSC(NM_csc(NM_B));
#endif
  }


  if(info != 0)
    THROW_EXCEPTION("SimpleMatrix::Solve failed.");
}

void SimpleMatrix::Solve(SiconosVector &B)
{
  DEBUG_BEGIN("SimpleMatrix::PLUSolve(SiconosVector &B)\n");

  if(!isFactorized())
  {
    Factorize();
  }

  // and then solve
  int info =1;
  NumericsMatrix * NM;
  double * b;
  SP::SiconosVector Bdense;


  if(B.num() == DENSE)
  {
    NM = _numericsMatrix.get();
    b = B.getArray();
  }
  else if(B.num() == SPARSE)
  {
    // First way. We copy to dense since our sparse solver is not able to take into account sparse r.h.s
    Bdense.reset (new SiconosVector(size(0)));
    * Bdense= B ;
    b = &(*Bdense->getArray());
    NM = _numericsMatrix.get();

    // Second way use inplace_solve of ublas
    // For that, we need to fill our factorization given by NM_LU_factorize into a ublas sparse matrix
    //inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
    //inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());
  }
  else
    THROW_EXCEPTION(" SimpleMatrix::Solve: only implemented for dense and sparse matrices in RHS.");


  if (isSymmetric())
  {
    if (isPositiveDefinite()) // Cholesky Factorization
    {
      //std::cout << "Cholesky Solve"<< std::endl;
      info  = NM_Cholesky_solve(NM, b, 1);

      if(info != 0)
      {
        THROW_EXCEPTION("SimpleMatrix::Solve failed (Cholesky)");
      }
    }
    else  //  LDLT Factorization
    {
      THROW_EXCEPTION("SimpleMatrix::Solve failed: LDL^T not yet implemented.");
    }
  }
  else //  LU Factorization  by default
  {
    info  = NM_LU_solve(NM, b, 1);

    if(info != 0)
    {
      THROW_EXCEPTION("SimpleMatrix::Solve failed (LU)");
    }
  }

  if(B.num() == SPARSE)
  {
    B = *Bdense ;                                                // we copy back to sparse.
  }





  if(info != 0)
    THROW_EXCEPTION("SimpleMatrix::Solve failed.");
  // else
  // {
  //   noalias(*(B.dense())) = ublas::column(tmpB, 0);
  // }
  DEBUG_END("SimpleMatrix::Solve(SiconosVector &B)\n");
}
