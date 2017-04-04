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
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/bindings/blas/level3.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
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

using namespace Siconos;



//======================
// Product of matrices
//======================

const SimpleMatrix prod(const SiconosMatrix &A, const SiconosMatrix& B)
{
  // To compute C = A * B
  assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function C=prod(A,B): inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();

  // == TODO: implement block product ==
  if (numA == 0 || numB == 0)
    SiconosMatrixException::selfThrow("Matrix product ( C=prod(A,B) ): not yet implemented for BlockMatrix objects.");

  if (numA == 7 || numB == 6) // A = identity or B = 0
    return SimpleMatrix(B);

  else if (numB == 7 || numA == 6) // B = identity or A = 0
    return SimpleMatrix(A);

  else // neither A or B is equal to identity or zero.
  {
    if (numB == 1)
    {
      if (numA == 1)
      {
        DenseMat p(A.size(0), B.size(1));
        siconosBindings::blas::gemm(1.0, *A.dense(), *B.dense(), 1.0, p);
        //      return (DenseMat)(prod(*A.dense(),*B.dense()));
        return p;
      }
      else if (numA == 2)
        return (DenseMat)(prod(*A.triang(), *B.dense()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.sym(), *B.dense()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.sparse(), *B.dense()));
      else// if(numA==5)
        return (DenseMat)(prod(*A.banded(), *B.dense()));
    }
    else if (numB == 2)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.dense(), *B.triang()));
      else if (numA == 2)
        return (TriangMat)(prod(*A.triang(), *B.triang()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.sym(), *B.triang()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.sparse(), *B.triang()));
      else //if(numA==5)
        return (DenseMat)(prod(*A.banded(), *B.triang()));
    }
    else if (numB == 3)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.dense(), *B.sym()));
      else if (numA == 2)
        return (DenseMat)(prod(*A.triang(), *B.sym()));
      else if (numA == 3)
        return (SymMat)(prod(*A.sym(), *B.sym()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.sparse(), *B.sym()));
      else // if (numA == 5)
        return (DenseMat)(prod(*A.banded(), *B.sym()));
    }
    else if (numB == 4)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.dense(), *B.sparse()));
      else if (numA == 2)
        return (DenseMat)(prod(*A.triang(), *B.sparse()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.sym(), *B.sparse()));
      else if (numA == 4)
        return (SparseMat)(prod(*A.sparse(), *B.sparse()));
      else //if(numA==5){
        return (DenseMat)(prod(*A.banded(), *B.sparse()));
    }
    else //if(numB==5)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.dense(), *B.banded()));
      else if (numA == 2)
        return (DenseMat)(prod(*A.triang(), *B.banded()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.sym(), *B.banded()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.sparse(), *B.banded()));
      else //if(numA==5)
        return (DenseMat)(prod(*A.banded(), *B.banded()));
    }
  }
}
/**

indexStart : indexStart[0] is the first raw, indexStart[1] is the first col
dim : dim[0] number of raw, dim[1] number of col
*/
// void zeroBlock(const SiconosMatrix& A, index indexStart, index dim){
//   ;
// }
// void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, int indexACol, bool init){
//   // To compute C[indexAcol::] = A * B

//   unsigned int numA = A.num();
//   unsigned int numB = B.num();
//   unsigned int numC = C.num();
//   if (numA == 0 || numB == 0 || numC == 0)
//     SiconosMatrixException::selfThrow("Matrix function prod(A,B,C,index): inconsistent sizes");
//   // === if C is zero or identity => read-only ===
//   if (numC == 6 || numC == 7)
//     SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C,index) ): wrong type for resulting matrix C (read-only: zero or identity).");


//   if (numA == 7 || numC == 6) // A = identity or 0
//     SiconosMatrixException::selfThrow("Matrix function prod(A,B,C,index): numA == 7 || numC == 6 not yet implemented");

//   int rawB = B.size(0);
//   int colB = B.size(1);

// }
void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init)
{
  // To compute C = A * B
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!" );
  assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!" );
  if(!C.isBlock())
    C.resetLU();

  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): inconsistent sizes");

  if (A.size(0) != C.size(0) || B.size(1) != C.size(1))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();
  unsigned int numC = C.num();

  // == TODO: implement block product ==
  if (numA == 0 || numB == 0)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): not yet implemented for BlockMatrix objects.");

  // === if C is zero or identity => read-only ===
  if (numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");


  if (numA == 7) // A = identity ...
  {
    if (init)
    {
      if (&C != &B) C = B; // if C and B are two different objects.
      // else nothing
    }
    else
      C += B;
  }

  else if (numB == 7) // B = identity
  {
    if (init)
    {
      if (&C != &A) C = A; // if C and A are two different objects.
      // else nothing
    }
    else
      C += A;
  }

  else if (numA == 6 || numB == 6) // if A or B = 0
  {
    if (init)
      C.zero();
    //else nothing
  }
  else if (numC == 0) // if C is Block - Temp. solution
  {
    SimpleMatrix tmp(C);
    prod(A, B, tmp, init);
    C = tmp;
  }
  else // neither A or B is equal to identity or zero.
  {
    if (init)
    {
      if (&C == &A) // if common memory between A and C
      {
        switch (numA)
        {
        case 1:
          if (numB == 1)
          {
            *C.dense()  = prod(*A.dense(), *B.dense());
            //siconosBindings::blas::gemm(1.0, *A.dense(), *B.dense(), 0.0, *C.dense());
          }
          else if (numB == 2)
            *C.dense()  = prod(*A.dense(), *B.triang());
          else if (numB == 3)
            *C.dense()  = prod(*A.dense(), *B.sym());
          else if (numB == 4)
            *C.dense()  = prod(*A.dense(), *B.sparse());
          else //if(numB==5)
            *C.dense() = prod(*A.dense(), *B.banded());
          break;
        case 2:
          if (numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.triang() = prod(*A.triang(), *B.triang());
          break;
        case 3:
          if (numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sym() = prod(*A.sym(), *B.sym());
          break;
        case 4:
          if (numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sparse() = prod(*A.sparse(), *B.sparse());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else if (&C == &B)
      {
        switch (numB)
        {
        case 1:
          if (numA == 1)
            *C.dense() = prod(*A.dense(), *B.dense());
          else if (numA == 2)
            *C.dense()  = prod(*A.triang(), *B.dense());
          else if (numA == 3)
            *C.dense()  = prod(*A.sym(), *B.dense());
          else if (numA == 4)
            *C.dense()  = prod(*A.sparse(), *B.dense());
          else //if(numB==5)
            *C.dense() = prod(*A.banded(), *B.dense());
          break;
        case 2:
          if (numA != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.triang() = prod(*A.triang(), *B.triang());
          break;
        case 3:
          if (numA != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sym() = prod(*A.sym(), *B.sym());
          break;
        case 4:
          if (numA != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sparse() = prod(*A.sparse(), *B.sparse());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else // if no alias between C and A or B.
      {
        switch (numC)
        {
        case 1:
          if (numB == 1)
          {
            if (numA == 1)
              noalias(*C.dense()) = prod(*A.dense(), *B.dense());
            else if (numA == 2)
              noalias(*C.dense()) = prod(*A.triang(), *B.dense());
            else if (numA == 3)
              noalias(*C.dense())  = prod(*A.sym(), *B.dense());
            else if (numA == 4)
              noalias(*C.dense()) = prod(*A.sparse(), *B.dense());
            else// if(numA==5)
              noalias(*C.dense())  = prod(*A.banded(), *B.dense());
          }
          else if (numB == 2)
          {
            if (numA == 1)
              noalias(*C.dense())  = prod(*A.dense(), *B.triang());
            else if (numA == 2)
              noalias(*C.dense())  = prod(*A.triang(), *B.triang());
            else if (numA == 3)
              noalias(*C.dense())  = prod(*A.sym(), *B.triang());
            else if (numA == 4)
              noalias(*C.dense())  = prod(*A.sparse(), *B.triang());
            else //if(numA==5)
              noalias(*C.dense())  = prod(*A.banded(), *B.triang());
          }
          else if (numB == 3)
          {
            if (numA == 1)
              noalias(*C.dense())  = prod(*A.dense(), *B.sym());
            else if (numA == 2)
              noalias(*C.dense())  = prod(*A.triang(), *B.sym());
            else if (numA == 3)
              noalias(*C.dense())  = prod(*A.sym(), *B.sym());
            else if (numA == 4)
              noalias(*C.dense())  = prod(*A.sparse(), *B.sym());
            else // if (numA == 5)
              noalias(*C.dense())  = prod(*A.banded(), *B.sym());
          }
          else if (numB == 4)
          {
            if (numA == 1)
              noalias(*C.dense()) = prod(*A.dense(), *B.sparse());
            else if (numA == 2)
              noalias(*C.dense()) = prod(*A.triang(), *B.sparse());
            else if (numA == 3)
              noalias(*C.dense()) = prod(*A.sym(), *B.sparse());
            else if (numA == 4)
              noalias(*C.dense()) = prod(*A.sparse(), *B.sparse());
            else //if(numA==5){
              noalias(*C.dense()) = prod(*A.banded(), *B.sparse());
          }
          else //if(numB==5)
          {
            if (numA == 1)
              noalias(*C.dense()) = prod(*A.dense(), *B.banded());
            else if (numA == 2)
              noalias(*C.dense()) = prod(*A.triang(), *B.banded());
            else if (numA == 3)
              noalias(*C.dense()) = prod(*A.sym(), *B.banded());
            else if (numA == 4)
              noalias(*C.dense()) = prod(*A.sparse(), *B.banded());
            else //if(numA==5)
              noalias(*C.dense()) = prod(*A.banded(), *B.banded());
          }
          break;
        case 2:
          if (numA != 2 || numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.triang()) = prod(*A.triang(), *B.triang());
          break;
        case 3:
          if (numA != 3 || numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.sym()) = prod(*A.sym(), *B.sym());
          break;
        case 4:
          if (numA != 4 || numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.sparse()) = prod(*A.sparse(), *B.sparse());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
    }
    else // += case
    {
      if (&C == &A) // if common memory between A and C
      {
        switch (numA)
        {
        case 1:
          if (numB == 1)
            *C.dense() += prod(*A.dense(), *B.dense());
          else if (numB == 2)
            *C.dense()  += prod(*A.dense(), *B.triang());
          else if (numB == 3)
            *C.dense()  += prod(*A.dense(), *B.sym());
          else if (numB == 4)
            *C.dense()  += prod(*A.dense(), *B.sparse());
          else //if(numB==5)
            *C.dense() += prod(*A.dense(), *B.banded());
          break;
        case 2:
          if (numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.triang() += prod(*A.triang(), *B.triang());
          break;
        case 3:
          if (numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sym() += prod(*A.sym(), *B.sym());
          break;
        case 4:
          if (numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sparse() += prod(*A.sparse(), *B.sparse());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else if (&C == &B)
      {
        switch (numB)
        {
        case 1:
          if (numA == 1)
            *C.dense() += prod(*A.dense(), *B.dense());
          else if (numA == 2)
            *C.dense()  += prod(*A.triang(), *B.dense());
          else if (numA == 3)
            *C.dense()  += prod(*A.sym(), *B.dense());
          else if (numA == 4)
            *C.dense()  += prod(*A.sparse(), *B.dense());
          else //if(numB==5)
            *C.dense() += prod(*A.banded(), *B.dense());
          break;
        case 2:
          if (numA != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.triang() += prod(*A.triang(), *B.triang());
          break;
        case 3:
          if (numA != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sym() += prod(*A.sym(), *B.sym());
          break;
        case 4:
          if (numA != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.sparse() += prod(*A.sparse(), *B.sparse());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else // if no alias between C and A or B.
      {
        switch (numC)
        {
        case 1:
          if (numB == 1)
          {
            if (numA == 1)
              noalias(*C.dense()) += prod(*A.dense(), *B.dense());
            else if (numA == 2)
              noalias(*C.dense()) += prod(*A.triang(), *B.dense());
            else if (numA == 3)
              noalias(*C.dense())  += prod(*A.sym(), *B.dense());
            else if (numA == 4)
              noalias(*C.dense()) += prod(*A.sparse(), *B.dense());
            else// if(numA==5)
              noalias(*C.dense())  += prod(*A.banded(), *B.dense());
          }
          else if (numB == 2)
          {
            if (numA == 1)
              noalias(*C.dense())  += prod(*A.dense(), *B.triang());
            else if (numA == 2)
              noalias(*C.dense())  += prod(*A.triang(), *B.triang());
            else if (numA == 3)
              noalias(*C.dense())  += prod(*A.sym(), *B.triang());
            else if (numA == 4)
              noalias(*C.dense())  += prod(*A.sparse(), *B.triang());
            else //if(numA==5)
              noalias(*C.dense())  += prod(*A.banded(), *B.triang());
          }
          else if (numB == 3)
          {
            if (numA == 1)
              noalias(*C.dense())  += prod(*A.dense(), *B.sym());
            else if (numA == 2)
              noalias(*C.dense())  += prod(*A.triang(), *B.sym());
            else if (numA == 3)
              noalias(*C.dense())  += prod(*A.sym(), *B.sym());
            else if (numA == 4)
              noalias(*C.dense())  += prod(*A.sparse(), *B.sym());
            else // if (numA == 5)
              noalias(*C.dense())  += prod(*A.banded(), *B.sym());
          }
          else if (numB == 4)
          {
            if (numA == 1)
              noalias(*C.dense()) += prod(*A.dense(), *B.sparse());
            else if (numA == 2)
              noalias(*C.dense()) += prod(*A.triang(), *B.sparse());
            else if (numA == 3)
              noalias(*C.dense()) += prod(*A.sym(), *B.sparse());
            else if (numA == 4)
              noalias(*C.dense()) += prod(*A.sparse(), *B.sparse());
            else //if(numA==5){
              noalias(*C.dense()) += prod(*A.banded(), *B.sparse());
          }
          else //if(numB==5)
          {
            if (numA == 1)
              noalias(*C.dense()) += prod(*A.dense(), *B.banded());
            else if (numA == 2)
              noalias(*C.dense()) += prod(*A.triang(), *B.banded());
            else if (numA == 3)
              noalias(*C.dense()) += prod(*A.sym(), *B.banded());
            else if (numA == 4)
              noalias(*C.dense()) += prod(*A.sparse(), *B.banded());
            else //if(numA==5)
              noalias(*C.dense()) += prod(*A.banded(), *B.banded());
          }
          break;
        case 2:
          if (numA != 2 || numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.triang()) += prod(*A.triang(), *B.triang());
          break;
        case 3:
          if (numA != 3 || numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.sym()) += prod(*A.sym(), *B.sym());
          break;
        case 4:
          if (numA != 4 || numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.sparse()) += prod(*A.sparse(), *B.sparse());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
    }
  if(!C.isBlock())
    C.resetLU();
  }
}


void axpy_prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init)
{

  // To compute C = A * B (init = true) or C += A * B (init = false) using ublas axpy_prod.
  // High speedup for sparse matrices.
  // Warning FP: ublas::axpy_prod(A, B, C, init) with init = True is equivalent
  // to C = A*B with C.clear BEFORE product. So C==A or B must be forbidden.
  // See http://www.boost.org/doc/libs/1_63_0/libs/numeric/ublas/doc/products.html
  //

  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): inconsistent sizes");

  if (A.size(0) != C.size(0) || B.size(1) != C.size(1))
    SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): inconsistent sizes");

  if(&A == &C || &B == &C)
    SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): C must be different from A and B.");

  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!" );
  assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!" );
  if(!C.isBlock())
    C.resetLU();
  unsigned int numA = A.num();
  unsigned int numB = B.num();
  unsigned int numC = C.num();
  // == TODO: implement block product ==
  if (numA == 0 || numB == 0)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): not yet implemented for BlockMatrix objects.");

  // === if C is zero or identity => read-only ===
  if (numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");


  if (numA == 7) // A = identity ...
  {
    if (!init)
      C += B;
    else
      C = B; // if C and B are two different objects.
  }

  else if (numB == 7) // B = identity
  {
    if (!init)
      C += A;
    else
      C = A; // if C and A are two different objects.
  }


  else if (numA == 6 || numB == 6) // if A or B = 0
  {
    if (init) C.zero(); // else nothing
  }
  else if (numC == 0) // if C is Block - Temp. solution
  {
    SimpleMatrix tmp(C);
    axpy_prod(A, B, tmp, init);
    C = tmp;
  }
  else // neither A or B is equal to identity or zero.
  {
    switch (numC)
    {
    case 1:
      if (numB == 1)
      {
        if (numA == 1)
          ublas::axpy_prod(*A.dense(), *B.dense(), *C.dense(), init);
        else if (numA == 2)
          ublas::axpy_prod(*A.triang(), *B.dense(), *C.dense(), init);
        else if (numA == 3)
          ublas::axpy_prod(*A.sym(), *B.dense(), *C.dense(), init);
        else if (numA == 4)
          ublas::axpy_prod(*A.sparse(), *B.dense(), *C.dense(), init);
        else// if(numA==5)
          ublas::axpy_prod(*A.banded(), *B.dense(), *C.dense(), init);
      }
      else if (numB == 2)
      {
        if (numA == 1)
          ublas::axpy_prod(*A.dense(), *B.triang(), *C.dense(), init);
        else if (numA == 2)
          ublas::axpy_prod(*A.triang(), *B.triang(), *C.dense(), init);
        else if (numA == 3)
          ublas::axpy_prod(*A.sym(), *B.triang(), *C.dense(), init);
        else if (numA == 4)
          ublas::axpy_prod(*A.sparse(), *B.triang(), *C.dense(), init);
        else //if(numA==5)
          ublas::axpy_prod(*A.banded(), *B.triang(), *C.dense(), init);
      }
      else if (numB == 3)
      {
        if (numA == 1)
          ublas::axpy_prod(*A.dense(), *B.sym(), *C.dense(), init);
        else if (numA == 2)
          ublas::axpy_prod(*A.triang(), *B.sym(), *C.dense(), init);
        else if (numA == 3)
          ublas::axpy_prod(*A.sym(), *B.sym(), *C.dense(), init);
        else if (numA == 4)
          ublas::axpy_prod(*A.sparse(), *B.sym(), *C.dense(), init);
        else // if (numA == 5)
          ublas::axpy_prod(*A.banded(), *B.sym(), *C.dense(), init);
      }
      else if (numB == 4)
      {
        if (numA == 1)
          ublas::axpy_prod(*A.dense(), *B.sparse(), *C.dense(), init);
        else if (numA == 2)
          ublas::axpy_prod(*A.triang(), *B.sparse(), *C.dense(), init);
        else if (numA == 3)
          ublas::axpy_prod(*A.sym(), *B.sparse(), *C.dense(), init);
        else if (numA == 4)
          ublas::axpy_prod(*A.sparse(), *B.sparse(), *C.dense(), init);
        else //if(numA==5){
          ublas::axpy_prod(*A.banded(), *B.sparse(), *C.dense(), init);
      }
      else //if(numB==5)
      {
        if (numA == 1)
          ublas::axpy_prod(*A.dense(), *B.banded(), *C.dense(), init);
        else if (numA == 2)
          ublas::axpy_prod(*A.triang(), *B.banded(), *C.dense(), init);
        else if (numA == 3)
          ublas::axpy_prod(*A.sym(), *B.banded(), *C.dense(), init);
        else if (numA == 4)
          ublas::axpy_prod(*A.sparse(), *B.banded(), *C.dense(), init);
        else //if(numA==5)
          ublas::axpy_prod(*A.banded(), *B.banded(), *C.dense(), init);
      }
      break;
    case 2:
      // if(numA!= 2 || numB != 2)
      SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      //ublas::axpy_prod(*A.triang(), *B.triang(),*C.triang(), init);
      break;
    case 3:
      //        if(numA!= 3 || numB != 3)
      SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      //ublas::axpy_prod(*A.sym(), *B.sym(),*C.sym(),init);
      break;
    case 4:
      if (numA != 4 || numB != 4)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      ublas::sparse_prod(*A.sparse(), *B.sparse(), *C.sparse(), init);
      break;
    default:
      SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
    }
    if(!C.isBlock())
      C.resetLU();
  }
}

void gemmtranspose(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");
  unsigned int numA = A.num();
  unsigned int numB = B.num();
  unsigned int numC = C.num();
  if (numA != 1 || numB != 1 || numC != 1)
    SiconosMatrixException::selfThrow("gemm(...) failed: reserved to dense matrices.");

  assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");


  siconosBindings::blas::gemm(a, siconosBindings::trans(*A.dense()), siconosBindings::trans(*B.dense()), b, *C.dense());



  C.resetLU();
}

void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  unsigned int numA = A.num();
  unsigned int numB = B.num();
  unsigned int numC = C.num();
  assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");
  C.resetLU();

  // At the time, only dense output allowed 
  DenseMat * tmpC = NULL;
  if (numA == 0 || numB == 0 || numC == 0)
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");

  if (numA == 1 && numB == 1 && numC == 1)
    siconosBindings::blas::gemm(a, *A.dense(), *B.dense(), b, *C.dense());
  else if (numA == 1 && numB == 1 && numC != 1)
  {
    // Copy C into tmpC ... 
    tmpC = new DenseMat(*C.dense());    
    siconosBindings::blas::gemm(a, *A.dense(), *B.dense(), b, *tmpC);
    std::cout << *tmpC << std::endl;
    noalias(*C.dense()) = *tmpC;
    delete tmpC;
  }
  else
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for these kinds of matrices.");
  C.resetLU();
}

void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool init)
{
  // To compute B = a * A (init = true) or B += a*A (init = false).
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!" );
  if(!B.isBlock())
    B.resetLU();

  if (&A == &B)
  {
    if (init) B *= a;
    else B *= (1.0 + a);
  }
  else
  {
    unsigned int numA = A.num();
    unsigned int numB = B.num();

    if (numB == 6 || numB == 7) // B = 0 or identity.
      SiconosMatrixException::selfThrow("scal(a,A,B) : forbidden for B being a zero or identity matrix.");

    if (numA == 6)
    {
      if (init) B.zero(); // else nothing
    }
    else if (numA == 7)
    {
      if (init)
      {
        B.eye();
        B *= a;
      }
      else
      {
        // Assuming B is square ...
        for (unsigned int i = 0; i < B.size(0); ++i)
          B(i, i) += a;
      }
    }
    else
    {
      if (numA == numB) // if A and B are of the same type ...
      {
        switch (numA)
        {

        case 0: // A and B are block
          if (isComparableTo(A, B))
          {
            const BlockMatrix& Aref = static_cast<const BlockMatrix&>(A);
            BlockMatrix& Bref = static_cast<BlockMatrix&>(B);
            BlocksIterator1 itB1;
            BlocksIterator2 itB2;
            ConstBlocksIterator1 itA1 = Aref._mat->begin1();
            ConstBlocksIterator2 itA2;
            for (itB1 = Bref._mat->begin1(); itB1 != Bref._mat->end1(); ++itB1)
            {
              itA2 = itA1.begin();
              for (itB2 = itB1.begin(); itB2 != itB1.end(); ++itB2)
              {
                scal(a, **itA2++, **itB2, init);
              }
              itA1++;
            }
          }
          else // if A and B are not "block-consistent"
          {
            if (init)
            {
              for (unsigned int i = 0; i < A.size(0); ++i)
                for (unsigned int j = 0; j < A.size(1); ++j)
                  B(i, j) = a * A(i, j);
            }
            else
            {
              for (unsigned int i = 0; i < A.size(0); ++i)
                for (unsigned int j = 0; j < A.size(1); ++j)
                  B(i, j) += a * A(i, j);
            }
          }
          break;

        case 1: // if both are dense
          if (init)
            noalias(*B.dense()) = a ** A.dense();
          else
            noalias(*B.dense()) += a ** A.dense();
          break;
        case 2:
          if (init)
            noalias(*B.triang()) = a ** A.triang();
          else
            noalias(*B.triang()) += a ** A.triang();
          break;
        case 3:
          if (init)
            noalias(*B.sym()) = a ** A.sym();
          else
            noalias(*B.sym()) += a ** A.sym();
          break;
        case 4:
          if (init)
            noalias(*B.sparse()) = a ** A.sparse();
          else
            noalias(*B.sparse()) += a ** A.sparse();
          break;
        case 5:
          if (init)
            noalias(*B.banded()) = a ** A.banded();
          else
            noalias(*B.banded()) += a ** A.banded();
          break;
        }
      }
      else // if A and B are of different types.
      {
        if (numA == 0 || numB == 0) // if A or B is block
        {
          if (init)
          {
            B = A;
            B *= a;
          }
          else
          {
            SimpleMatrix tmp(A);
            tmp *= a;
            B += tmp; // bof bof ...
          }
        }
        else
        {
          if (numB != 1)
            SiconosMatrixException::selfThrow("scal(a,A,B) failed. A and B types do not fit together.");

          if (init)
          {
            switch (numB)
            {
            case 1:
              noalias(*B.dense()) = a ** A.dense();
              break;
            case 2:
              noalias(*B.dense()) = a ** A.triang();
              break;
            case 3:
              noalias(*B.dense()) = a ** A.sym();
              break;
            case 4:
              noalias(*B.dense()) = a ** A.sparse();
              break;
            case 5:
              noalias(*B.dense()) = a ** A.banded();
              break;
            }
          }
          else

          {
            switch (numB)
            {
            case 1:
              noalias(*B.dense()) += a ** A.dense();
              break;
            case 2:
              noalias(*B.dense()) += a ** A.triang();
              break;
            case 3:
              noalias(*B.dense()) += a ** A.sym();
              break;
            case 4:
              noalias(*B.dense()) += a ** A.sparse();
              break;
            case 5:
              noalias(*B.dense()) += a ** A.banded();
              break;
            }
          }
        }
      }
    }
  }
}


