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
#include "debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

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
      _isPLUFactorizedInPlace = true;
      SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace failed: the matrix is singular.");
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
      SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace failed. ");
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
  if(_num != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices.");

#if defined(HAS_LAPACK_dgetri)
  int info = lapack::getri(*mat.Dense, *_ipiv);   // solve from factorization

  if(info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace failed, the matrix is singular.");

  _isPLUInversed = true;
#else
  SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace not implemented with lapack.");
#endif
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  if(B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(B) failed at solving Ax = B. Not yet implemented for a BlockMatrix B.");
  int info = 0;

  if(_num == 1)
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
        SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices in RHS.");

  }
  else
  {
    if(!_isPLUFactorized)  // call first PLUFactorizationInPlace
    {
      PLUFactorizationInPlace();
    }
    // and then solve
    if(B.num() == 1)
    {
      inplace_solve(*sparse(), *(B.dense()), ublas::lower_tag());
      inplace_solve(ublas::trans(*sparse()), *(B.dense()), ublas::upper_tag());
    }
    else if(B.num() == 4)
    {
      inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
      inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());
    }
    else
      SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense ans sparse matrices in RHS.");
    info = 0 ;
  }
  //  SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace: only implemented for dense matrices.");



  if(info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  if(B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(V) failed. Not yet implemented for V being a BlockVector.");


  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.dense()); // Conversion of vector to matrix. Temporary solution.
  int info;

  if(_num == 1)
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
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
  else
  {
    noalias(*(B.dense())) = ublas::column(tmpB, 0);
  }
}

void SimpleMatrix::resetLU()
{
  if(_ipiv) _ipiv->clear();
  _isPLUFactorized = false;
  _isPLUFactorizedInPlace = false;
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
  if(B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::SolveByLeastSquares(Siconos Matrix &B) failed. Not yet implemented for M being a BlockMatrix.");
  int info = 0;
#ifdef USE_OPTIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, *(B.dense()), lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
  info += lapack::gels(*mat.Dense, *(B.dense()), lapack::minimal_workspace());
#endif
  if(info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::SolveByLeastSquares failed.");
}




void SimpleMatrix::SolveByLeastSquares(SiconosVector &B)
{
  if(B.isBlock())
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
  if(info != 0)
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


void SimpleMatrix::PLUFactorize()
{
  DEBUG_BEGIN("void SimpleMatrix::PLUFactorize()\n");
  if(_isPLUFactorized)
  {
    std::cout << "SimpleMatrix::PLUFactorize warning: this matrix is already PLUFactorized. " << std::endl;
    return;
  }
  if(_num == DENSE)
  {
    _numericsMatrix.reset(NM_new(),NM_clear);
    NumericsMatrix * NM = _numericsMatrix.get();
    double * data = (double*)(getArray());
    DEBUG_EXPR(NV_display(data,size(0)*size(1)););
    NM_fill(NM, NM_DENSE, size(0), size(1), data ); // Pointer link
    DEBUG_EXPR(NM_display(NM));

    int info  = NM_LU_factorize(NM);
    DEBUG_EXPR(NM_display(NM));

    if(info != 0)
    {
      _isPLUFactorized = false;
      _isPLUFactorizedInPlace = true;
      SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorize failed: the matrix is singular.");
    }
    else
    {
      _isPLUFactorized = true;
      _isPLUFactorizedInPlace = true;
    }
  }
  else
  {
    // We build a sparse matrix and we call numerics for the LU factorization of a sparse matrix.
    _numericsMatrix.reset(NM_create(NM_SPARSE, size(0), size(1)),NM_clear);
    NumericsMatrix * NM = _numericsMatrix.get();
    _numericsMatrix->matrix2->origin = NSM_CSC;
    NM_csc_alloc(NM, nnz());

    fillCSC(numericsSparseMatrix(NM)->csc, std::numeric_limits<double>::epsilon());

    DEBUG_EXPR(cs_print(numericsSparseMatrix(NM)->csc, 0););

    NM_LU_factorize(NM);

    NSM_linear_solver_params* p = NSM_linearSolverParams(NM);
    CSparseMatrix_factors* cs_lu_A =  (CSparseMatrix_factors*)p->linear_solver_data;

    int info = 1;
    if (cs_lu_A)
    {
      info = 0;

      DEBUG_EXPR(CSparseMatrix * L = cs_lu_A->N->L; cs_print(L,0););
      DEBUG_EXPR(CSparseMatrix * U = cs_lu_A->N->U; cs_print(U,0););
      DEBUG_EXPR(
        CS_INT * pinv = cs_lu_A->N->pinv;
        for (int i =0; i < size(0) ; i++)
          printf("pinv[%i] = %i\n", i, pinv[i])
        );
    }

    if(info != 0)
    {
      _isPLUFactorized = false;
      _isPLUFactorizedInPlace = false;
      std::cout << "Problem in LU Factorization" << std::endl;
      SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorize failed. ");
    }
    else
    {
      _isPLUFactorized = true;
      _isPLUFactorizedInPlace = false;
    }
  }
  DEBUG_END("void SimpleMatrix::PLUFactorize()\n");
}
void SimpleMatrix::PLUSolve(SiconosMatrix &B)
{
  if(B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUSolve(B) failed at solving Ax = B. Not yet implemented for a BlockMatrix B.");
  int info = 0;
  if(!_isPLUFactorized)  // call gesv => LU-factorize+solve
  {
    PLUFactorize();
  }
  // and then solve
  if(B.num() == DENSE)
  {
    NumericsMatrix * NM = _numericsMatrix.get();
    double * b = B.getArray();
    info = NM_LU_solve(NM, b, B.size(1));
  }
  else if(B.num() == SPARSE)
  {
    // First way.
    // We copy to dense since our sparse solver is not able to take
    // into account for sparse r.h.s, yet.
    SP::SimpleMatrix Bdense (new SimpleMatrix(size(0),size(1)));
    * Bdense= B ;                                                // copy to dense
    double * b = &(*Bdense->getArray());
    NumericsMatrix * NM = _numericsMatrix.get();
    info = NM_LU_solve(NM, b, B.size(1));       // we solve with dense R.H.S.
    DEBUG_EXPR(NV_display(b,size(0)*size(1)););
    DEBUG_EXPR(Bdense->display(););
    B = *Bdense ;                                                // we copy back to sparse.

    // Second way
    // use inplace_solve of ublas (see above with PLUSolveInPlace)
    // For that, we need to fill our factorization given by NM_LU_factorize
    // into a ublas sparse matrix
    // inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
    // inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());
  }
  else
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUSolve: only implemented for dense and sparse matrices in RHS.");

  if(info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUSolve failed.");
}
void SimpleMatrix::PLUSolve(SiconosVector &B)
{
  DEBUG_BEGIN("SimpleMatrix::PLUSolve(SiconosVector &B)\n");
  if(B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUSolve(V) failed. Not yet implemented for V being a BlockVector.");

  if(!_isPLUFactorized)  // call gesv => LU-factorize+solve
  {
    PLUFactorize();
  }


  // and then solve
  int info =1;
  if(B.num() == 0)
  {
    NumericsMatrix * NM = _numericsMatrix.get();
    double * b = B.getArray();
    DEBUG_EXPR(NV_display(b, B.size()););
    DEBUG_EXPR(B.display(););
    info = NM_LU_solve(NM, b, 1);
    DEBUG_EXPR(NV_display(b, B.size()););
    DEBUG_EXPR(B.display(););
  }
  else if(B.num() == 1)
  {
    // First way. We copy to dense since our sparse solver is not able to take into account sparse r.h.s
    SP::SiconosVector Bdense (new SiconosVector(size(0)));
    * Bdense= B ;
    double * b = &(*Bdense->getArray());
    NumericsMatrix * NM = _numericsMatrix.get();
    info = NM_LU_solve(NM, b, 1);
    B = *Bdense ;

    // Second way use inplace_solve of ublas
    // For that, we need to fill our factorization given by NM_LU_factorize into a ublas sparse matrix
    //inplace_solve(*sparse(), *(B.sparse()), ublas::lower_tag());
    //inplace_solve(ublas::trans(*sparse()), *(B.sparse()), ublas::upper_tag());
  }
  else
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUSolve: only implemented for dense and sparse matrices in RHS.");

  if(info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUSolve failed.");
  // else
  // {
  //   noalias(*(B.dense())) = ublas::column(tmpB, 0);
  // }
  DEBUG_END("SimpleMatrix::PLUSolve(SiconosVector &B)\n");
}
