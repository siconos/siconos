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



#include "SiconosConfig.h"
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/bindings/blas/level3.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "SiconosAlgebraProd.hpp" // for axpy_prod and prod
// Note Franck : sounds useless. It seems it's defined in bindings
// (to be checked, especially on windows)

// #define BIND_FORTRAN_LOWERCASE_UNDERSCORE

// needed for blas3
#include <assert.h>

namespace siconosBindings = boost::numeric::bindings;

// for ublas::axpy_prod, ...
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

// require for matrix stuff like value_type
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"
#include "SiconosAlgebraProd.hpp" // for prod

using namespace siconos;



//======================
// Product of matrices
//======================
// Note FP: this function is never used. We keep it for the record. Remove it later ?


// const SimpleMatrix prod(const SiconosMatrix &A, const SiconosMatrix& B)
// {
//   // To compute C = A * B
//   assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");


//   if((A.size(1) != B.size(0)))
//     THROW_EXCEPTION("Matrix function C=prod(A,B): inconsistent sizes");

//   siconos::UBLAS_TYPE numA = A.num();
//   siconos::UBLAS_TYPE numB = B.num();

//   // == TODO: implement block product ==
//   if(numA == 0 || numB == 0)
//     THROW_EXCEPTION("Matrix product ( C=prod(A,B) ): not yet implemented for BlockMatrix objects.");

//   if(numA == siconos::IDENTITY || numB == siconos::ZERO)  // A = identity or B = 0
//     return SimpleMatrix(B);

//   else if(numB == siconos::IDENTITY || numA == siconos::ZERO)  // B = identity or A = 0
//     return SimpleMatrix(A);

//   else // neither A or B is equal to identity or zero.
//   {
//     if(numB == siconos::DENSE)
//     {
//       if(numA == siconos::DENSE)
//       {
//         DenseMat p(A.size(0), B.size(1));
//         siconosBindings::blas::gemm(1.0, *A.dense(), *B.dense(), 1.0, p);
//         //      return (DenseMat)(prod(*A.dense(),*B.dense()));
//         return p;
//       }
//       else if(numA == siconos::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.dense()));
//       else if(numA == siconos::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.dense()));
//       else if(numA == siconos::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.dense()));
//       else if(numA == siconos::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.dense()));
//       else// if(numA==siconos::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.dense()));
//     }
//     else if(numB == siconos::TRIANGULAR)
//     {
//       if(numA == siconos::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.triang()));
//       else if(numA == siconos::TRIANGULAR)
//         return (TriangMat)(prod(*A.triang(), *B.triang()));
//       else if(numA == siconos::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.triang()));
//       else if(numA == siconos::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.triang()));
//       else if(numA == siconos::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.triang()));
//       else //if(numA==siconos::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.triang()));
//     }
//     else if(numB == siconos::SYMMETRIC)
//     {
//       if(numA == siconos::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.sym()));
//       else if(numA == siconos::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.sym()));
//       else if(numA == siconos::SYMMETRIC)
//         return (SymMat)(prod(*A.sym(), *B.sym()));
//       else if(numA == siconos::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.sym()));
//       else if(numA == siconos::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.sym()));
//       else // if (numA == siconos::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.sym()));
//     }
//     else if(numB == siconos::SPARSE)
//     {
//       if(numA == siconos::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.sparse()));
//       else if(numA == siconos::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.sparse()));
//       else if(numA == siconos::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.sparse()));
//       else if(numA == siconos::SPARSE)
//         return (SparseMat)(prod(*A.sparse(), *B.sparse()));
//       else if(numA == siconos::SPARSE_COORDINATE)
//         return (SparseMat)(prod(*A.sparseCoordinate(), *B.sparse()));
//       else //if(numA==siconos::BANDED){
//         return (DenseMat)(prod(*A.banded(), *B.sparse()));
//     }
//     else if(numB == siconos::SPARSE_COORDINATE)
//     {
//       if(numA == siconos::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.sparseCoordinate()));
//       else if(numA == siconos::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.sparseCoordinate()));
//       else if(numA == siconos::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.sparseCoordinate()));
//       else if(numA == siconos::SPARSE)
//         return (SparseMat)(prod(*A.sparse(), *B.sparseCoordinate()));
//       else if(numA == siconos::SPARSE_COORDINATE)
//         return (SparseMat)(prod(*A.sparseCoordinate(), *B.sparseCoordinate()));
//       else //if(numA==siconos::BANDED){
//         return (DenseMat)(prod(*A.banded(), *B.sparseCoordinate()));
//     }
//     else //if(numB==siconos::BANDED)
//     {
//       if(numA == siconos::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.banded()));
//       else if(numA == siconos::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.banded()));
//       else if(numA == siconos::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.banded()));
//       else if(numA == siconos::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.banded()));
//       else if(numA == siconos::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.banded()));
//       else //if(numA==siconos::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.banded()));
//     }
//   }
// }
/**

indexStart : indexStart[0] is the first raw, indexStart[1] is the first col
dim : dim[0] number of raw, dim[1] number of col
*/
// void zeroBlock(const SiconosMatrix& A, index indexStart, index dim){
//   ;
// }
// void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, int indexACol, bool init){
//   // To compute C[indexAcol::] = A * B

//   siconos::UBLAS_TYPE numA = A.num();
//   siconos::UBLAS_TYPE numB = B.num();
//   siconos::UBLAS_TYPE numC = C.num();
//   if (numA == 0 || numB == 0 || numC == 0)
//     THROW_EXCEPTION("Matrix function prod(A,B,C,index): inconsistent sizes");
//   // === if C is zero or identity => read-only ===
//   if (numC == siconos::Zero || numC == siconos::IDENTITY)
//     THROW_EXCEPTION("Matrix product ( prod(A,B,C,index) ): wrong type for resulting matrix C (read-only: zero or identity).");


//   if (numA == siconos::IDENTITY || numC == siconos::Zero) // A = identity or 0
//     THROW_EXCEPTION("Matrix function prod(A,B,C,index): numA == siconos::IDENTITY || numC == siconos::Zero not yet implemented");

//   int rawB = B.size(0);
//   int colB = B.size(1);

// }


void axpy_prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init)
{

  // To compute C = A * B (init = true) or C += A * B (init = false) using ublas axpy_prod.
  // High speedup for sparse matrices.
  // Warning FP: ublas::axpy_prod(A, B, C, init) with init = True is equivalent
  // to C = A*B with C.clear BEFORE product. So C==A or B must be forbidden.
  // See http://www.boost.org/doc/libs/1_63_0/libs/numeric/ublas/doc/products.html
  //

  if((A.size(1) != B.size(0)))
    THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): inconsistent sizes");

  if(A.size(0) != C.size(0) || B.size(1) != C.size(1))
    THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): inconsistent sizes");

  if(&A == &C || &B == &C)
    THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): C must be different from A and B.");

  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorized in place in prod !!");
  assert(!(B.isPLUFactorizedInPlace()) && "B is PLUFactorized in place in prod !!");
  if(!C.isBlock())
    C.resetFactorizationFlags();
  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numB = B.num();
  siconos::UBLAS_TYPE numC = C.num();
  // == TODO: implement block product ==
  if(numA == siconos::BLOCK || numB == siconos::BLOCK)
    THROW_EXCEPTION("Matrix product ( prod(A,B,C) ): not yet implemented for BlockMatrix objects.");

  // === if C is zero or identity => read-only ===
  if(numC == siconos::ZERO || numC == siconos::IDENTITY)
    THROW_EXCEPTION("Matrix product ( prod(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");


  if(numA == siconos::IDENTITY)  // A = identity ...
  {
    if(!init)
      C += B;
    else
      C = B; // if C and B are two different objects.
  }

  else if(numB == siconos::IDENTITY)  // B = identity
  {
    if(!init)
      C += A;
    else
      C = A; // if C and A are two different objects.
  }


  else if(numA == siconos::ZERO || numB == siconos::ZERO)  // if A or B = 0
  {
    if(init) C.zero();  // else nothing
  }
  else if(numC == siconos::BLOCK)  // if C is Block - Temp. solution
  {
    SimpleMatrix tmp(C);
    axpy_prod(A, B, tmp, init);
    C = tmp;
  }
  else // neither A or B is equal to identity or zero.
  {
    switch(numC)
    {
    case siconos::DENSE:
      if(numB == siconos::DENSE)
      {
        if(numA == siconos::DENSE)
          ublas::axpy_prod(*A.dense(), *B.dense(), *C.dense(), init);
        else if(numA == siconos::TRIANGULAR)
          ublas::axpy_prod(*A.triang(), *B.dense(), *C.dense(), init);
        else if(numA == siconos::SYMMETRIC)
          ublas::axpy_prod(*A.sym(), *B.dense(), *C.dense(), init);
        else if(numA == siconos::SPARSE)
          ublas::axpy_prod(*A.sparse(), *B.dense(), *C.dense(), init);
        else// if(numA==siconos::BANDED)
          ublas::axpy_prod(*A.banded(), *B.dense(), *C.dense(), init);
      }
      else if(numB == siconos::TRIANGULAR)
      {
        if(numA == siconos::DENSE)
          ublas::axpy_prod(*A.dense(), *B.triang(), *C.dense(), init);
        else if(numA == siconos::TRIANGULAR)
          ublas::axpy_prod(*A.triang(), *B.triang(), *C.dense(), init);
        else if(numA == siconos::SYMMETRIC)
          ublas::axpy_prod(*A.sym(), *B.triang(), *C.dense(), init);
        else if(numA == siconos::SPARSE)
          ublas::axpy_prod(*A.sparse(), *B.triang(), *C.dense(), init);
        else //if(numA==siconos::BANDED)
          ublas::axpy_prod(*A.banded(), *B.triang(), *C.dense(), init);
      }
      else if(numB == siconos::SYMMETRIC)
      {
        if(numA == siconos::DENSE)
          ublas::axpy_prod(*A.dense(), *B.sym(), *C.dense(), init);
        else if(numA == siconos::TRIANGULAR)
          ublas::axpy_prod(*A.triang(), *B.sym(), *C.dense(), init);
        else if(numA == siconos::SYMMETRIC)
          ublas::axpy_prod(*A.sym(), *B.sym(), *C.dense(), init);
        else if(numA == siconos::SPARSE)
          ublas::axpy_prod(*A.sparse(), *B.sym(), *C.dense(), init);
        else // if (numA == siconos::BANDED)
          ublas::axpy_prod(*A.banded(), *B.sym(), *C.dense(), init);
      }
      else if(numB == siconos::SPARSE)
      {
        if(numA == siconos::DENSE)
          ublas::axpy_prod(*A.dense(), *B.sparse(), *C.dense(), init);
        else if(numA == siconos::TRIANGULAR)
          ublas::axpy_prod(*A.triang(), *B.sparse(), *C.dense(), init);
        else if(numA == siconos::SYMMETRIC)
          ublas::axpy_prod(*A.sym(), *B.sparse(), *C.dense(), init);
        else if(numA == siconos::SPARSE)
          ublas::axpy_prod(*A.sparse(), *B.sparse(), *C.dense(), init);
        else //if(numA==siconos::BANDED){
          ublas::axpy_prod(*A.banded(), *B.sparse(), *C.dense(), init);
      }
      else //if(numB==siconos::BANDED)
      {
        if(numA == siconos::DENSE)
          ublas::axpy_prod(*A.dense(), *B.banded(), *C.dense(), init);
        else if(numA == siconos::TRIANGULAR)
          ublas::axpy_prod(*A.triang(), *B.banded(), *C.dense(), init);
        else if(numA == siconos::SYMMETRIC)
          ublas::axpy_prod(*A.sym(), *B.banded(), *C.dense(), init);
        else if(numA == siconos::SPARSE)
          ublas::axpy_prod(*A.sparse(), *B.banded(), *C.dense(), init);
        else //if(numA==siconos::BANDED)
          ublas::axpy_prod(*A.banded(), *B.banded(), *C.dense(), init);
      }
      break;
    case siconos::TRIANGULAR:
      // if(numA!= siconos::TRIANGULAR || numB != siconos::TRIANGULAR)
      THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      //ublas::axpy_prod(*A.triang(), *B.triang(),*C.triang(), init);
      break;
    case siconos::SYMMETRIC:
      //        if(numA!= siconos::SYMMETRIC || numB != siconos::SYMMETRIC)
      THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      //ublas::axpy_prod(*A.sym(), *B.sym(),*C.sym(),init);
      break;
    case siconos::SPARSE:
      if(numA != siconos::SPARSE || numB != siconos::SPARSE)
        THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      ublas::sparse_prod(*A.sparse(), *B.sparse(), *C.sparse(), init);
      break;
    default:
      THROW_EXCEPTION("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
    }
    if(!C.isBlock())
      C.resetFactorizationFlags();
  }
}


// Note FP: this function is never used. We keep it for the record. Remove it later ?
// void gemmtranspose(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
// {
//   if(A.isBlock() || B.isBlock() || C.isBlock())
//     THROW_EXCEPTION("gemm(...) not yet implemented for block matrices.");
//   siconos::UBLAS_TYPE numA = A.num();
//   siconos::UBLAS_TYPE numB = B.num();
//   siconos::UBLAS_TYPE numC = C.num();
//   if(numA != siconos::DENSE || numB != siconos::DENSE || numC != siconos::DENSE)
//     THROW_EXCEPTION("gemm(...) failed: reserved to dense matrices.");


//   assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");


//   siconosBindings::blas::gemm(a, siconosBindings::trans(*A.dense()), siconosBindings::trans(*B.dense()), b, *C.dense());


//   C.resetFactorizationFlags();
// }

// Note FP: this function is never used. We keep it for the record. Remove it later ?
// void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
// {
//   siconos::UBLAS_TYPE numA = A.num();
//   siconos::UBLAS_TYPE numB = B.num();
//   siconos::UBLAS_TYPE numC = C.num();
//   assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");
//   C.resetFactorizationFlags();

//   // At the time, only dense output allowed
//   DenseMat * tmpC = nullptr;
//   if(numA == 0 || numB == 0 || numC == 0)
//     THROW_EXCEPTION("gemm(...) not yet implemented for block matrices.");

//   if(numA == siconos::DENSE && numB == siconos::DENSE && numC == siconos::DENSE)
//     siconosBindings::blas::gemm(a, *A.dense(), *B.dense(), b, *C.dense());
//   else if(numA == siconos::DENSE && numB == siconos::DENSE && numC != siconos::DENSE)
//   {
//     // Copy C into tmpC ...
//     tmpC = new DenseMat(*C.dense());
//     siconosBindings::blas::gemm(a, *A.dense(), *B.dense(), b, *tmpC);
//     std::cout << *tmpC << std::endl;
//     noalias(*C.dense()) = *tmpC;
//     delete tmpC;
//   }
//   else
//     THROW_EXCEPTION("gemm(...) not yet implemented for these kinds of matrices.");
//   C.resetFactorizationFlags();
// }

