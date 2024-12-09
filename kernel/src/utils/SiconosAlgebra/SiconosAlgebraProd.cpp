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
#include "SiconosAlgebraProd.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "SiconosMatrix.hpp"
#include "SimpleMatrix.hpp"
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosException.hpp"

void prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y, bool init)
{
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

  unsigned int startRow = 0;
  VectorOfVectors::const_iterator it;
  // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
  //       // private_prod takes into account the fact that x and y[i] may be block vectors.
  for(it = y.begin(); it != y.end(); ++it)
  {
    A.private_prod(startRow, x, **it, init);
    startRow += (*it)->size();
  }

}

void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init)
{

  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");


  if(init)
    y.zero();
  unsigned int startRow = 0;
  unsigned int startCol = 0;
  // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]), with subA a submatrix of A,
  // starting from position startRow in rows and startCol in columns.
  // private_prod takes also into account the fact that each block of x can also be a block.
  VectorOfVectors::const_iterator it;
  for(it = x.begin(); it != x.end(); ++it)
  {
    A.private_addprod(startRow, startCol, **it, y);
    startCol += (*it)->size();
  }
}

void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += A*x if init = false.
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

  if(A.size(1) != x.size())
    THROW_EXCEPTION("inconsistent sizes between A and x.")

  if(A.size(0) != y.size())
    THROW_EXCEPTION("inconsistent sizes between A and y.");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numX = x.num();
  siconos::UBLAS_TYPE numY = y.num();

  if(numA == siconos::BLOCK)
    THROW_EXCEPTION("not yet implemented for block matrices.");

  if(numA == siconos::ZERO)
  {
    if(init)
      y.zero();
    //else nothing
  }

  else if(numA == siconos::IDENTITY)
  {
    if(!init)
      y += x;
    else
    {
      if(&x != &y) y = x ;  // if x and y do not share memory (ie are different objects)
      // else nothing
    }
  }

  else // A is not 0 or identity
  {

    // === First case: y is not a block vector ===
    if(init)
    {
      if(&x != &y)  // if no common memory between x and y.
      {
        if(numX == siconos::DENSE)
        {
          if(numY != siconos::DENSE)
            THROW_EXCEPTION("y (output) must be a dense vector.");

          assert(y.dense() != x.dense());

          if(numA == siconos::DENSE)
            noalias(*y.dense()) = ublas::prod(*A.dense(), *x.dense());
          else if(numA == siconos::TRIANGULAR)
            noalias(*y.dense()) = ublas::prod(*A.triang(), *x.dense());
          else if(numA == siconos::SYMMETRIC)
            noalias(*y.dense()) = ublas::prod(*A.sym(), *x.dense());
          else if(numA == siconos::SPARSE)
            noalias(*y.dense()) = ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==siconos::BANDED)
            noalias(*y.dense()) = ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == siconos::SPARSE)
        {
          if(numY != siconos::DENSE && numA != siconos::SPARSE)
            THROW_EXCEPTION(" y (output) must be a dense vector or A a sparse matrix.");

          if(numA == siconos::DENSE)
            noalias(*y.dense()) = ublas::prod(*A.dense(), *x.sparse());
          else if(numA == siconos::TRIANGULAR)
            noalias(*y.dense()) = ublas::prod(*A.triang(), *x.sparse());
          else if(numA == siconos::SYMMETRIC)
            noalias(*y.dense()) = ublas::prod(*A.sym(), *x.sparse());
          else if(numA == siconos::SPARSE)
          {
            if(numY == siconos::DENSE)
              noalias(*y.dense()) = ublas::prod(*A.sparse(), *x.sparse());
            else
              noalias(*y.sparse()) = ublas::prod(*A.sparse(), *x.sparse());
          }
          else //if(numA==siconos::BANDED)
            noalias(*y.dense()) = ublas::prod(*A.banded(), *x.sparse());
        }
      }
      else // if x and y are the same object => alias
      {
        if(numX == siconos::DENSE)
        {
          if(numA == siconos::DENSE)
            *y.dense() = ublas::prod(*A.dense(), *x.dense());
          else if(numA == siconos::TRIANGULAR)
            *y.dense() = ublas::prod(*A.triang(), *x.dense());
          else if(numA == siconos::SYMMETRIC)
            *y.dense() = ublas::prod(*A.sym(), *x.dense());
          else if(numA == siconos::SPARSE)
            *y.dense() = ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==siconos::BANDED)
            *y.dense() = ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == siconos::SPARSE)
        {
          if(numA == siconos::DENSE)
            *y.sparse() = ublas::prod(*A.dense(), *x.sparse());
          else if(numA == siconos::TRIANGULAR)
            *y.sparse() = ublas::prod(*A.triang(), *x.sparse());
          else if(numA == siconos::SYMMETRIC)
            *y.sparse() = ublas::prod(*A.sym(), *x.sparse());
          else if(numA == siconos::SPARSE)
            *y.sparse() = ublas::prod(*A.sparse(), *x.sparse());
          else //if(numA==siconos::BANDED)
            *y.sparse() = ublas::prod(*A.banded(), *x.sparse());
        }
      }
    }
    else // += case
    {
      if(&x != &y)  // if no common memory between x and y.
      {
        if(numX == siconos::DENSE)
        {
          if(numY != siconos::DENSE)
            THROW_EXCEPTION("y (output) must be a dense vector.");

          if(numA == siconos::DENSE)
            noalias(*y.dense()) += ublas::prod(*A.dense(), *x.dense());
          else if(numA == siconos::TRIANGULAR)
            noalias(*y.dense()) += ublas::prod(*A.triang(), *x.dense());
          else if(numA == siconos::SYMMETRIC)
            noalias(*y.dense()) += ublas::prod(*A.sym(), *x.dense());
          else if(numA == siconos::SPARSE)
            noalias(*y.dense()) += ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==siconos::BANDED)
            noalias(*y.dense()) += ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == siconos::SPARSE)
        {
          if(numY != siconos::DENSE && numA != siconos::SPARSE)
            THROW_EXCEPTION("y (output) must be a dense vector or A a sparse matrix.");

          if(numA == siconos::DENSE)
            noalias(*y.dense()) += ublas::prod(*A.dense(), *x.sparse());
          else if(numA == siconos::TRIANGULAR)
            noalias(*y.dense()) += ublas::prod(*A.triang(), *x.sparse());
          else if(numA == siconos::SYMMETRIC)
            noalias(*y.dense()) += ublas::prod(*A.sym(), *x.sparse());
          else if(numA == siconos::SPARSE)
          {
            if(numY == siconos::DENSE)
              noalias(*y.dense()) += ublas::prod(*A.sparse(), *x.sparse());
            else
              noalias(*y.sparse()) += ublas::prod(*A.sparse(), *x.sparse());
          }
          else //if(numA==siconos::BANDED)
            noalias(*y.dense()) += ublas::prod(*A.banded(), *x.sparse());
        }
      }
      else // if x and y are the same object => alias
      {
        if(numX == siconos::DENSE)
        {
          if(numA == siconos::DENSE)
            *y.dense() += ublas::prod(*A.dense(), *x.dense());
          else if(numA == siconos::TRIANGULAR)
            *y.dense() += ublas::prod(*A.triang(), *x.dense());
          else if(numA == siconos::SYMMETRIC)
            *y.dense() += ublas::prod(*A.sym(), *x.dense());
          else if(numA == siconos::SPARSE)
            *y.dense() += ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==siconos::BANDED)
            *y.dense() += ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == siconos::SPARSE)
        {
          if(numA == siconos::DENSE)
            *y.sparse() += ublas::prod(*A.dense(), *x.sparse());
          else if(numA == siconos::TRIANGULAR)
            *y.sparse() += ublas::prod(*A.triang(), *x.sparse());
          else if(numA == siconos::SYMMETRIC)
            *y.sparse() += ublas::prod(*A.sym(), *x.sparse());
          else if(numA == siconos::SPARSE)
            *y.sparse() += ublas::prod(*A.sparse(), *x.sparse());
          else //if(numA==siconos::BANDED)
            *y.sparse() += ublas::prod(*A.banded(), *x.sparse());
        }
      }
    }
  }
}

void prod(const SiconosVector& x, const SiconosMatrix& A, SiconosVector& y, bool init)
{
  // To compute y = trans(A) * x in an "optimized" way, if init = true
  // (or y = trans(A) * x + y if init = false
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

  if(A.size(0) != x.size())
    THROW_EXCEPTION("inconsistent sizes between A and x.");

  if(A.size(1) != y.size())
    THROW_EXCEPTION("inconsistent sizes between A and y.");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numX = x.num();
  siconos::UBLAS_TYPE numY = y.num();

  if(numA == siconos::BLOCK)
    THROW_EXCEPTION("not yet implemented for block matrices.");

  if(numA == siconos::ZERO)  // A = 0
  {
    if(init)
      y.zero();
    // else nothing
  }

  else if(numA == siconos::IDENTITY)  // A = identity
  {
    if(!init)
      y += x;
    else
    {
      if(&x != &y) y = x ;  // if x and y do not share memory (ie are different objects)
      // else nothing
    }
  }

  else // A is not 0 or identity
  {
    {
      if(init)
      {

        if(&x != &y)  // if no common memory between x and y.
        {
          if(numX == siconos::DENSE)
          {
            if(numY != siconos::DENSE)
              THROW_EXCEPTION("y (output) must be a dense vector.");

            if(numA == siconos::DENSE)
              noalias(*y.dense()) = ublas::prod(trans(*A.dense()), *x.dense());
            else if(numA == siconos::TRIANGULAR)
              noalias(*y.dense()) = ublas::prod(trans(*A.triang()), *x.dense());
            else if(numA == siconos::SYMMETRIC)
              noalias(*y.dense()) = ublas::prod(trans(*A.sym()), *x.dense());
            else if(numA == siconos::SPARSE)
              noalias(*y.dense()) = ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==siconos::BANDED)
              noalias(*y.dense()) = ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == siconos::SPARSE)
          {
            if(numY != siconos::DENSE && numA != siconos::SPARSE)
              THROW_EXCEPTION("y (output) must be a dense vector or A a sparse matrix.");
            if(numA == siconos::DENSE)
              noalias(*y.dense()) = ublas::prod(trans(*A.dense()), *x.sparse());
            else if(numA == siconos::TRIANGULAR)
              noalias(*y.dense()) = ublas::prod(trans(*A.triang()), *x.sparse());
            else if(numA == siconos::SYMMETRIC)
              noalias(*y.dense()) = ublas::prod(trans(*A.sym()), *x.sparse());
            else if(numA == siconos::SPARSE)
            {
              if(numY == siconos::DENSE)
                noalias(*y.dense()) = ublas::prod(trans(*A.sparse()), *x.sparse());
              else
                noalias(*y.sparse()) = ublas::prod(trans(*A.sparse()), *x.sparse());
            }
            else //if(numA==siconos::BANDED)
              noalias(*y.dense()) = ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
        else // if x and y are the same object => alias
        {
          if(numX == siconos::DENSE)
          {
            if(numA == siconos::DENSE)
              *y.dense() = ublas::prod(trans(*A.dense()), *x.dense());
            else if(numA == siconos::TRIANGULAR)
              *y.dense() = ublas::prod(trans(*A.triang()), *x.dense());
            else if(numA == siconos::SYMMETRIC)
              *y.dense() = ublas::prod(trans(*A.sym()), *x.dense());
            else if(numA == siconos::SPARSE)
              *y.dense() = ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==siconos::BANDED)
              *y.dense() = ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == siconos::SPARSE)
          {
            if(numA == siconos::DENSE)
              *y.sparse() = ublas::prod(trans(*A.dense()), *x.sparse());
            else if(numA == siconos::TRIANGULAR)
              *y.sparse() = ublas::prod(trans(*A.triang()), *x.sparse());
            else if(numA == siconos::SYMMETRIC)
              *y.sparse() = ublas::prod(trans(*A.sym()), *x.sparse());
            else if(numA == siconos::SPARSE)
              *y.sparse() = ublas::prod(trans(*A.sparse()), *x.sparse());
            else //if(numA==siconos::BANDED)
              *y.sparse() = ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
      }
      else // += case
      {

        if(&x != &y)  // if no common memory between x and y.
        {
          if(numX == siconos::DENSE)
          {
            if(numY != siconos::DENSE)
              THROW_EXCEPTION("y (output) must be a dense vector.");

            if(numA == siconos::DENSE)
              noalias(*y.dense()) += ublas::prod(trans(*A.dense()), *x.dense());
            else if(numA == siconos::TRIANGULAR)
              noalias(*y.dense()) += ublas::prod(trans(*A.triang()), *x.dense());
            else if(numA == siconos::SYMMETRIC)
              noalias(*y.dense()) += ublas::prod(trans(*A.sym()), *x.dense());
            else if(numA == siconos::SPARSE)
              noalias(*y.dense()) += ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==siconos::BANDED)
              noalias(*y.dense()) += ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == siconos::SPARSE)
          {
            if(numY != siconos::DENSE && numA != siconos::SPARSE)
              THROW_EXCEPTION("y (output) must be a dense vector or A a sparse matrix.");

            if(numA == siconos::DENSE)
              noalias(*y.dense()) += ublas::prod(trans(*A.dense()), *x.sparse());
            else if(numA == siconos::TRIANGULAR)
              noalias(*y.dense()) += ublas::prod(trans(*A.triang()), *x.sparse());
            else if(numA == siconos::SYMMETRIC)
              noalias(*y.dense()) += ublas::prod(trans(*A.sym()), *x.sparse());
            else if(numA == siconos::SPARSE)
            {
              if(numY == siconos::DENSE)
                noalias(*y.dense()) += ublas::prod(trans(*A.sparse()), *x.sparse());
              else
                noalias(*y.sparse()) += ublas::prod(trans(*A.sparse()), *x.sparse());
            }
            else //if(numA==siconos::BANDED)
              noalias(*y.dense()) += ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
        else // if x and y are the same object => alias
        {
          if(numX == siconos::DENSE)
          {
            if(numA == siconos::DENSE)
              *y.dense() += ublas::prod(trans(*A.dense()), *x.dense());
            else if(numA == siconos::TRIANGULAR)
              *y.dense() += ublas::prod(trans(*A.triang()), *x.dense());
            else if(numA == siconos::SYMMETRIC)
              *y.dense() += ublas::prod(trans(*A.sym()), *x.dense());
            else if(numA == siconos::SPARSE)
              *y.dense() += ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==siconos::BANDED)
              *y.dense() += ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == siconos::SPARSE)
          {
            if(numA == siconos::DENSE)
              *y.sparse() += ublas::prod(trans(*A.dense()), *x.sparse());
            else if(numA == siconos::TRIANGULAR)
              *y.sparse() += ublas::prod(trans(*A.triang()), *x.sparse());
            else if(numA == siconos::SYMMETRIC)
              *y.sparse() += ublas::prod(trans(*A.sym()), *x.sparse());
            else if(numA == siconos::SPARSE)
              *y.sparse() += ublas::prod(trans(*A.sparse()), *x.sparse());
            else //if(numA==siconos::BANDED)
              *y.sparse() += ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
      }
    }
  }
}

void prod(const SiconosVector& x, const SiconosMatrix& A, BlockVector& y, bool init)
{
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

  if(A.size(0) != x.size())
    THROW_EXCEPTION("inconsistent sizes between A and x.");

  if(A.size(1) != y.size())
    THROW_EXCEPTION("inconsistent sizes between A and y.");

  unsigned int pos = 0;
  VectorOfVectors::const_iterator it;
  // For Each subvector of y, y[i], computes y[i] = transpose(subA) x, subA being a submatrix of A corresponding to y[i] position.
  // private_prod takes into account the fact that x and y[i] may be block vectors.
  for(it = y.begin(); it != y.end(); ++it)
  {
    taxpy(createSPtrConstSiconosVector(x), createSPtrConstSiconosMatrix(A), pos, 0, *it, init);
    pos += (*it)->size();
  }
}


// ========== Products matrix - vector

SiconosVector prod(const SiconosMatrix& A, const SiconosVector& x)
{
  // To compute y = A * x
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

  if(A.size(1) != x.size())
    THROW_EXCEPTION("inconsistent sizes between A and x.");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numX = x.num();

  if(numA == siconos::BLOCK)  // if A is block ...
    THROW_EXCEPTION("Not implemented for block matrices.");

  if(numA == siconos::ZERO)  // A = 0
    return (DenseVect)(ublas::zero_vector<double>(x.size()));

  else if(numA == siconos::IDENTITY)  // A = Identity
    return x;

  else
  {
    if(numX == siconos::DENSE)
    {
      if(numA == siconos::DENSE)
        return (DenseVect)(prod(*A.dense(), *x.dense()));
      else if(numA == siconos::TRIANGULAR)
        return (DenseVect)(prod(*A.triang(), *x.dense()));
      else if(numA == siconos::SYMMETRIC)
        return (DenseVect)(prod(*A.sym(), *x.dense()));
      else if(numA == siconos::SPARSE)
        return (DenseVect)(prod(*A.sparse(), *x.dense()));
      else // if(numA==siconos::BANDED)
        return (DenseVect)(prod(*A.banded(), *x.dense()));
    }
    else //if(numX == siconos::SPARSE)
    {
      if(numA == siconos::DENSE)
        return (DenseVect)(prod(*A.dense(), *x.sparse()));
      else if(numA == siconos::TRIANGULAR)
        return (DenseVect)(prod(*A.triang(), *x.sparse()));
      else if(numA == siconos::SYMMETRIC)
        return (DenseVect)(prod(*A.sym(), *x.sparse()));
      else if(numA == siconos::SPARSE)
        return (DenseVect)(prod(*A.sparse(), *x.sparse()));
      else // if(numA==siconos::BANDED)
        return (DenseVect)(prod(*A.banded(), *x.sparse()));
    }
  }
}

const SimpleMatrix  prod(const SiconosMatrix& A, const SiconosMatrix& B)
{
  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numB = B.num();

  if (numA == numB)
  {
    SimpleMatrix  C(A.size(0),B.size(1), numA);
    prod(A, B, C);
    return C;
  }
  else
  {
    SimpleMatrix  C(A.size(0),B.size(1));
    prod(A, B, C);
    return C;
  }
}
void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init)
{
  // To compute C = A * B
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
  assert(!(B.isPLUFactorizedInPlace()) && "B is PLUFactorizedInPlace in prod !!");
  if(!C.isBlock())
    C.resetFactorizationFlags();

  if((A.size(1) != B.size(0)))
    THROW_EXCEPTION("inconsistent sizes between A and B");

  if(A.size(0) != C.size(0) || B.size(1) != C.size(1))
    THROW_EXCEPTION("inconsistent sizes between A and C or B and C.");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numB = B.num();
  siconos::UBLAS_TYPE numC = C.num();

  // == TODO: implement block product ==
  if(numA == siconos::BLOCK || numB == siconos::BLOCK)
    THROW_EXCEPTION("not yet implemented for block matrices");

  // === if C is zero or identity => read-only ===
  if(numC == siconos::ZERO || numC == siconos::IDENTITY)
    THROW_EXCEPTION("wrong type for resulting matrix C (read-only: zero or identity).");


  if(numA == siconos::IDENTITY)  // A = identity ...
  {
    if(init)
    {
      if(&C != &B) C = B;  // if C and B are two different objects.
      // else nothing
    }
    else
      C += B;
  }

  else if(numB == siconos::IDENTITY)  // B = identity
  {
    if(init)
    {
      if(&C != &A) C = A;  // if C and A are two different objects.
      // else nothing
    }
    else
      C += A;
  }

  else if(numA == siconos::ZERO || numB == siconos::ZERO)  // if A or B = 0
  {
    if(init)
      C.zero();
    //else nothing
  }
  else if(numC == siconos::BLOCK)  // if C is Block - Temp. solution
  {
    SimpleMatrix tmp(C);
    prod(A, B, tmp, init);
    C = tmp;
  }
  else // neither A or B is equal to identity or zero.
  {
    if(init)
    {
      if(&C == &A)  // if common memory between A and C
      {
        switch(numA)
        {
        case siconos::DENSE:
          if(numB == siconos::DENSE)
          {
            *C.dense()  = prod(*A.dense(), *B.dense());
            //siconosBindings::blas::gemm(1.0, *A.dense(), *B.dense(), 0.0, *C.dense());
          }
          else if(numB == siconos::TRIANGULAR)
            *C.dense()  = prod(*A.dense(), *B.triang());
          else if(numB == siconos::SYMMETRIC)
            *C.dense()  = prod(*A.dense(), *B.sym());
          else if(numB == siconos::SPARSE)
            *C.dense()  = prod(*A.dense(), *B.sparse());
          else //if(numB==siconos::BANDED)
            *C.dense() = prod(*A.dense(), *B.banded());
          break;
        case siconos::TRIANGULAR:
          if(numB != siconos::TRIANGULAR)
            THROW_EXCEPTION("wrong type for B (according to A type).");
          *C.triang() = prod(*A.triang(), *B.triang());
          break;
        case siconos::SYMMETRIC:
          if(numB != siconos::SYMMETRIC)
            THROW_EXCEPTION("wrong type for B (according to A type).");
          *C.sym() = prod(*A.sym(), *B.sym());
          break;
        case siconos::SPARSE:
          if(numB != siconos::SPARSE)
            THROW_EXCEPTION("wrong type for B (according to A type).");
          *C.sparse() = prod(*A.sparse(), *B.sparse());
          break;
        default:
          THROW_EXCEPTION("not implemented for A type.");
        }
      }
      else if(&C == &B)
      {
        switch(numB)
        {
        case siconos::DENSE:
          if(numA == siconos::DENSE)
            *C.dense() = prod(*A.dense(), *B.dense());
          else if(numA == siconos::TRIANGULAR)
            *C.dense()  = prod(*A.triang(), *B.dense());
          else if(numA == siconos::SYMMETRIC)
            *C.dense()  = prod(*A.sym(), *B.dense());
          else if(numA == siconos::SPARSE)
            *C.dense()  = prod(*A.sparse(), *B.dense());
          else //if(numB==siconos::BANDED)
            *C.dense() = prod(*A.banded(), *B.dense());
          break;
        case siconos::TRIANGULAR:
          if(numA != siconos::TRIANGULAR)
            THROW_EXCEPTION("wrong type for A (according to B type).");
          *C.triang() = prod(*A.triang(), *B.triang());
          break;
        case siconos::SYMMETRIC:
          if(numA != siconos::SYMMETRIC)
            THROW_EXCEPTION("wrong type for A (according to B type).");
          *C.sym() = prod(*A.sym(), *B.sym());
          break;
        case siconos::SPARSE:
          if(numA != siconos::SPARSE)
            THROW_EXCEPTION("wrong type for A (according to B type).");
          *C.sparse() = prod(*A.sparse(), *B.sparse());
          break;
        default:
          THROW_EXCEPTION("not implemented for B type.");
        }
      }
      else // if no alias between C and A or B.
      {
        switch(numC)
        {
        case siconos::DENSE:
          if(numB == siconos::DENSE)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense()) = prod(*A.dense(), *B.dense());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) = prod(*A.triang(), *B.dense());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense())  = prod(*A.sym(), *B.dense());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) = prod(*A.sparse(), *B.dense());
            else// if(numA==siconos::BANDED)
              noalias(*C.dense())  = prod(*A.banded(), *B.dense());
          }
          else if(numB == siconos::TRIANGULAR)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense())  = prod(*A.dense(), *B.triang());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense())  = prod(*A.triang(), *B.triang());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense())  = prod(*A.sym(), *B.triang());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense())  = prod(*A.sparse(), *B.triang());
            else //if(numA==siconos::BANDED)
              noalias(*C.dense())  = prod(*A.banded(), *B.triang());
          }
          else if(numB == siconos::SYMMETRIC)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense())  = prod(*A.dense(), *B.sym());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense())  = prod(*A.triang(), *B.sym());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense())  = prod(*A.sym(), *B.sym());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense())  = prod(*A.sparse(), *B.sym());
            else // if (numA == siconos::BANDED)
              noalias(*C.dense())  = prod(*A.banded(), *B.sym());
          }
          else if(numB == siconos::SPARSE)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense()) = prod(*A.dense(), *B.sparse());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) = prod(*A.triang(), *B.sparse());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense()) = prod(*A.sym(), *B.sparse());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) = prod(*A.sparse(), *B.sparse());
            else //if(numA==siconos::BANDED){
              noalias(*C.dense()) = prod(*A.banded(), *B.sparse());
          }
          else //if(numB==siconos::BANDED)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense()) = prod(*A.dense(), *B.banded());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) = prod(*A.triang(), *B.banded());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense()) = prod(*A.sym(), *B.banded());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) = prod(*A.sparse(), *B.banded());
            else //if(numA==siconos::BANDED)
              noalias(*C.dense()) = prod(*A.banded(), *B.banded());
          }
          break;
        case siconos::TRIANGULAR:
          if(numA != siconos::TRIANGULAR || numB != siconos::TRIANGULAR)
            THROW_EXCEPTION("wrong type for A or B (according to C type).");
          noalias(*C.triang()) = prod(*A.triang(), *B.triang());
          break;
        case siconos::SYMMETRIC:
          if(numA != siconos::SYMMETRIC || numB != siconos::SYMMETRIC)
            THROW_EXCEPTION("wrong type for A or B (according to C type).");
          noalias(*C.sym()) = prod(*A.sym(), *B.sym());
          break;
        case siconos::SPARSE:
          if(numA != siconos::SPARSE || numB != siconos::SPARSE)
            THROW_EXCEPTION("wrong type for A or B (according to C type).");
          noalias(*C.sparse()) = prod(*A.sparse(), *B.sparse());
          break;
        default:
          THROW_EXCEPTION("not implemented for C type.");
        }
      }
    }
    else // += case
    {
      if(&C == &A)  // if common memory between A and C
      {
        switch(numA)
        {
        case siconos::DENSE:
          if(numB == siconos::DENSE)
            *C.dense() += prod(*A.dense(), *B.dense());
          else if(numB == siconos::TRIANGULAR)
            *C.dense()  += prod(*A.dense(), *B.triang());
          else if(numB == siconos::SYMMETRIC)
            *C.dense()  += prod(*A.dense(), *B.sym());
          else if(numB == siconos::SPARSE)
            *C.dense()  += prod(*A.dense(), *B.sparse());
          else //if(numB==siconos::BANDED)
            *C.dense() += prod(*A.dense(), *B.banded());
          break;
        case siconos::TRIANGULAR:
          if(numB != siconos::TRIANGULAR)
            THROW_EXCEPTION("wrong type for B (according to A type).");
          *C.triang() += prod(*A.triang(), *B.triang());
          break;
        case siconos::SYMMETRIC:
          if(numB != siconos::SYMMETRIC)
            THROW_EXCEPTION("wrong type for B (according to A type).");
          *C.sym() += prod(*A.sym(), *B.sym());
          break;
        case siconos::SPARSE:
          if(numB != siconos::SPARSE)
            THROW_EXCEPTION("wrong type for B (according to A type).");
          *C.sparse() += prod(*A.sparse(), *B.sparse());
          break;
        default:
          THROW_EXCEPTION("not implemented for A type.");
        }
      }
      else if(&C == &B)
      {
        switch(numB)
        {
        case siconos::DENSE:
          if(numA == siconos::DENSE)
            *C.dense() += prod(*A.dense(), *B.dense());
          else if(numA == siconos::TRIANGULAR)
            *C.dense()  += prod(*A.triang(), *B.dense());
          else if(numA == siconos::SYMMETRIC)
            *C.dense()  += prod(*A.sym(), *B.dense());
          else if(numA == siconos::SPARSE)
            *C.dense()  += prod(*A.sparse(), *B.dense());
          else //if(numB==siconos::BANDED)
            *C.dense() += prod(*A.banded(), *B.dense());
          break;
        case siconos::TRIANGULAR:
          if(numA != siconos::TRIANGULAR)
            THROW_EXCEPTION("wrong type for A (according to B type).");
          *C.triang() += prod(*A.triang(), *B.triang());
          break;
        case siconos::SYMMETRIC:
          if(numA != siconos::SYMMETRIC)
            THROW_EXCEPTION("wrong type for A (according to B type).");
          *C.sym() += prod(*A.sym(), *B.sym());
          break;
        case siconos::SPARSE:
          if(numA != siconos::SPARSE)
            THROW_EXCEPTION("wrong type for A (according to B type).");
          *C.sparse() += prod(*A.sparse(), *B.sparse());
          break;
        default:
            THROW_EXCEPTION("not yet implemented for A type.");
        }
      }
      else // if no alias between C and A or B.
      {
        switch(numC)
        {
        case siconos::DENSE:
          if(numB == siconos::DENSE)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense()) += prod(*A.dense(), *B.dense());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) += prod(*A.triang(), *B.dense());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense())  += prod(*A.sym(), *B.dense());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) += prod(*A.sparse(), *B.dense());
            else// if(numA==siconos::BANDED)
              noalias(*C.dense())  += prod(*A.banded(), *B.dense());
          }
          else if(numB == siconos::TRIANGULAR)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense())  += prod(*A.dense(), *B.triang());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense())  += prod(*A.triang(), *B.triang());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense())  += prod(*A.sym(), *B.triang());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense())  += prod(*A.sparse(), *B.triang());
            else //if(numA==siconos::BANDED)
              noalias(*C.dense())  += prod(*A.banded(), *B.triang());
          }
          else if(numB == siconos::SYMMETRIC)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense())  += prod(*A.dense(), *B.sym());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense())  += prod(*A.triang(), *B.sym());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense())  += prod(*A.sym(), *B.sym());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense())  += prod(*A.sparse(), *B.sym());
            else // if (numA == BANDED)
              noalias(*C.dense())  += prod(*A.banded(), *B.sym());
          }
          else if(numB == siconos::SPARSE)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense()) += prod(*A.dense(), *B.sparse());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) += prod(*A.triang(), *B.sparse());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense()) += prod(*A.sym(), *B.sparse());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) += prod(*A.sparse(), *B.sparse());
            else //if(numA==siconos::BANDED){
              noalias(*C.dense()) += prod(*A.banded(), *B.sparse());
          }
          else //if(numB==siconos::BANDED)
          {
            if(numA == siconos::DENSE)
              noalias(*C.dense()) += prod(*A.dense(), *B.banded());
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) += prod(*A.triang(), *B.banded());
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense()) += prod(*A.sym(), *B.banded());
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) += prod(*A.sparse(), *B.banded());
            else //if(numA==siconos::BANDED)
              noalias(*C.dense()) += prod(*A.banded(), *B.banded());
          }
          break;
        case siconos::TRIANGULAR:
          if(numA != siconos::TRIANGULAR || numB != siconos::TRIANGULAR)
            THROW_EXCEPTION("wrong type for A or B (according to C type).");
          noalias(*C.triang()) += prod(*A.triang(), *B.triang());
          break;
        case siconos::SYMMETRIC:
          if(numA != siconos::SYMMETRIC || numB != siconos::SYMMETRIC)
            THROW_EXCEPTION("wrong type for A or B (according to C type).");
          noalias(*C.sym()) += prod(*A.sym(), *B.sym());
          break;
        case siconos::SPARSE:
          if(numA != siconos::SPARSE || numB != siconos::SPARSE)
            THROW_EXCEPTION("wrong type for A or B (according to C type).");
          noalias(*C.sparse()) += prod(*A.sparse(), *B.sparse());
          break;
        default:
            THROW_EXCEPTION("not implemented for C type).");
        }
      }
    }
    if(!C.isBlock())
      C.resetFactorizationFlags();
  }
}

void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = a*A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += a*A*x if init = false.
  assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

  if(A.size(1) != x.size())
    THROW_EXCEPTION("inconsistent sizes between A and x.");

  if(A.size(0) != y.size())
    THROW_EXCEPTION("inconsistent sizes between A and y.");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numX = x.num();
  siconos::UBLAS_TYPE numY = y.num();

  if(numA == siconos::BLOCK)  // If A is Block
    THROW_EXCEPTION("not yet implemented for block matrices.");

  if(numA == siconos::ZERO)  // A = 0
  {
    if(init)
      y.zero();
    //else nothing
  }

  else if(numA == siconos::IDENTITY)  // A = identity
  {
    scal(a, x, y, init);
  }

  else // A is not 0 or identity
  {

    // === First case: y is not a block vector ===
    {
      {
        if(init)
        {
          if(&x != &y)  // if no common memory between x and y.
          {
            if(numX == siconos::DENSE)
            {
              if(numY != siconos::DENSE)
                THROW_EXCEPTION("y (output) must be a dense vector.");

              if(numA == siconos::DENSE)
                noalias(*y.dense()) = a * ublas::prod(*A.dense(), *x.dense());
              else if(numA == siconos::TRIANGULAR)
                noalias(*y.dense()) = a * ublas::prod(*A.triang(), *x.dense());
              else if(numA == siconos::SYMMETRIC)
                noalias(*y.dense()) = a * ublas::prod(*A.sym(), *x.dense());
              else if(numA == siconos::SPARSE)
                noalias(*y.dense()) = a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==siconos::BANDED)
                noalias(*y.dense()) = a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == siconos::SPARSE)
            {
              if(numY != siconos::DENSE && numA != siconos::SPARSE)
               THROW_EXCEPTION("y (output) must be a dense vector or A a sparse matrix.");

              if(numA == siconos::DENSE)
                noalias(*y.dense()) = a * ublas::prod(*A.dense(), *x.sparse());
              else if(numA == siconos::TRIANGULAR)
                noalias(*y.dense()) = a * ublas::prod(*A.triang(), *x.sparse());
              else if(numA == siconos::SYMMETRIC)
                noalias(*y.dense()) = a * ublas::prod(*A.sym(), *x.sparse());
              else if(numA == siconos::SPARSE)
              {
                if(numY == siconos::DENSE)
                  noalias(*y.dense()) = a * ublas::prod(*A.sparse(), *x.sparse());
                else
                  noalias(*y.sparse()) = a * ublas::prod(*A.sparse(), *x.sparse());
              }
              else //if(numA==siconos::BANDED)
                noalias(*y.dense()) = a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
          else // if x and y are the same object => alias
          {
            if(numX == siconos::DENSE)
            {
              if(numA == siconos::DENSE)
                *y.dense() = a * ublas::prod(*A.dense(), *x.dense());
              else if(numA == siconos::TRIANGULAR)
                *y.dense() = a * ublas::prod(*A.triang(), *x.dense());
              else if(numA == siconos::SYMMETRIC)
                *y.dense() = a * ublas::prod(*A.sym(), *x.dense());
              else if(numA == siconos::SPARSE)
                *y.dense() = a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==siconos::BANDED)
                *y.dense() = a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == siconos::SPARSE)
            {
              if(numA == siconos::DENSE)
                *y.sparse() = a * ublas::prod(*A.dense(), *x.sparse());
              else if(numA == siconos::TRIANGULAR)
                *y.sparse() = a * ublas::prod(*A.triang(), *x.sparse());
              else if(numA == siconos::SYMMETRIC)
                *y.sparse() = a * ublas::prod(*A.sym(), *x.sparse());
              else if(numA == siconos::SPARSE)
                *y.sparse() = a * ublas::prod(*A.sparse(), *x.sparse());
              else //if(numA==siconos::BANDED)
                *y.sparse() = a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
        }
        else // += case
        {
          if(&x != &y)  // if no common memory between x and y.
          {
            if(numX == siconos::DENSE)
            {
              if(numY != siconos::DENSE)
               THROW_EXCEPTION("y (output) must be a dense vector.");

              if(numA == siconos::DENSE)
                noalias(*y.dense()) += a * ublas::prod(*A.dense(), *x.dense());
              else if(numA == siconos::TRIANGULAR)
                noalias(*y.dense()) += a * ublas::prod(*A.triang(), *x.dense());
              else if(numA == siconos::SYMMETRIC)
                noalias(*y.dense()) += a * ublas::prod(*A.sym(), *x.dense());
              else if(numA == siconos::SPARSE)
                noalias(*y.dense()) += a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==siconos::BANDED)
                noalias(*y.dense()) += a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == siconos::SPARSE)
            {
              if(numY != siconos::DENSE && numA != siconos::SPARSE)
                THROW_EXCEPTION("y (output) must be a dense vector or A a sparse matrix.");

              if(numA == siconos::DENSE)
                noalias(*y.dense()) += a * ublas::prod(*A.dense(), *x.sparse());
              else if(numA == siconos::TRIANGULAR)
                noalias(*y.dense()) += a * ublas::prod(*A.triang(), *x.sparse());
              else if(numA == siconos::SYMMETRIC)
                noalias(*y.dense()) += a * ublas::prod(*A.sym(), *x.sparse());
              else if(numA == siconos::SPARSE)
              {
                if(numY == siconos::DENSE)
                  noalias(*y.dense()) += a * ublas::prod(*A.sparse(), *x.sparse());
                else
                  noalias(*y.sparse()) += a * ublas::prod(*A.sparse(), *x.sparse());
              }
              else //if(numA==siconos::BANDED)
                noalias(*y.dense()) += a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
          else // if x and y are the same object => alias
          {
            if(numX == siconos::DENSE)
            {
              if(numA == siconos::DENSE)
                *y.dense() += a * ublas::prod(*A.dense(), *x.dense());
              else if(numA == siconos::TRIANGULAR)
                *y.dense() += a * ublas::prod(*A.triang(), *x.dense());
              else if(numA == siconos::SYMMETRIC)
                *y.dense() += a * ublas::prod(*A.sym(), *x.dense());
              else if(numA == siconos::SPARSE)
                *y.dense() += a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==siconos::BANDED)
                *y.dense() += a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == siconos::SPARSE)
            {
              if(numA == siconos::DENSE)
                *y.sparse() += a * ublas::prod(*A.dense(), *x.sparse());
              else if(numA == siconos::TRIANGULAR)
                *y.sparse() += a * ublas::prod(*A.triang(), *x.sparse());
              else if(numA == siconos::SYMMETRIC)
                *y.sparse() += a * ublas::prod(*A.sym(), *x.sparse());
              else if(numA == siconos::SPARSE)
                *y.sparse() += a * ublas::prod(*A.sparse(), *x.sparse());
              else //if(numA==siconos::BANDED)
                *y.sparse() += a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
        }
      }
    }
  }
}

void taxpy(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y, bool init)
{
  assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY
  if(init)  // y = subA * x , else y += subA * x
    y->zero();

  if(A->isBlock())
    THROW_EXCEPTION("not yet implemented for block matrix.");

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
  // Then computation of y = subA*x + y.
  siconos::UBLAS_TYPE numA = A->num();
  siconos::UBLAS_TYPE numY = y->num();
  siconos::UBLAS_TYPE numX = x->num();
  unsigned int sizeX = x->size();
  unsigned int sizeY = y->size();

  if(numX != numY)
    THROW_EXCEPTION("not yet implemented for x and y of different types.");

  if(numY == siconos::DENSE && numX == siconos::DENSE)
  {

    assert(y->dense() != x->dense());

    if(numA == siconos::DENSE)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->dense()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if(numA == siconos::TRIANGULAR)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->triang()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if(numA == siconos::SYMMETRIC)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->sym()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if(numA == siconos::SPARSE)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->sparse()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else //if(numA==siconos::BANDED)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->banded()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
  }
  else // x and y sparse
  {
    if(numA == siconos::SPARSE)
      *y->sparse() += prod(ublas::subrange(trans(*A->sparse()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->sparse());
    else
      THROW_EXCEPTION("not yet implemented for x, y  sparse and A not sparse.");
  }
}

