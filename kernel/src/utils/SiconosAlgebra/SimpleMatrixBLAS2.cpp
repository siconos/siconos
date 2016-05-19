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
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

// Note Franck : sounds useless. It seems it's defined in bindings
// (to be checked, especially on windows)

// #define BIND_FORTRAN_LOWERCASE_UNDERSCORE
namespace siconosBindings = boost::numeric::bindings;

// for ublas::axpy_prod, ...
#include <boost/numeric/ublas/operation.hpp>

// for matrix stuff like value_type
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"
#include "SiconosAlgebra.hpp"

using namespace Siconos;

// ========== Products matrix - vector //

// Computation of y = A*x
//
// Two specific functions are used to handle all the cases where x or y are blocks.
// All of their blocks can also be blocks. Then we use:
// - private_prod to "slice" A when y is block, ie according to its rows.
// - private_addprod to "slice" A when x is block, and to sum over the columns of blocks to compute y = sum subA x[i].

// The following function is private and used inside prod(...) public functions.
// It is required to deal with block vectors of blocks ( of blocks ...).
// It computes res = subA*x +res, subA being a submatrix of A (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
// If x is a block vector, it call the present function for all blocks.
const SiconosVector prod(const SiconosMatrix& A, const SiconosVector& x)
{
  // To compute y = A * x
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: inconsistent sizes.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();

  if (numA == 0) // if A is block ...
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: not yet implemented for block matrix.");

  if (numA == 6) // A = 0
    return (DenseVect)(ublas::zero_vector<double>(x.size()));

  else if (numA == 7) // A = Identity
    return x;

  else
  {
    if (numX == 1)
    {
      if (numA == 1)
        return (DenseVect)(prod(*A.dense(), *x.dense()));
      else if (numA == 2)
        return (DenseVect)(prod(*A.triang(), *x.dense()));
      else if (numA == 3)
        return (DenseVect)(prod(*A.sym(), *x.dense()));
      else if (numA == 4)
        return (DenseVect)(prod(*A.sparse(), *x.dense()));
      else // if(numA==5)
        return (DenseVect)(prod(*A.banded(), *x.dense()));
    }
    else //if(numX == 4)
    {
      if (numA == 1)
        return (DenseVect)(prod(*A.dense(), *x.sparse()));
      else if (numA == 2)
        return (DenseVect)(prod(*A.triang(), *x.sparse()));
      else if (numA == 3)
        return (DenseVect)(prod(*A.sym(), *x.sparse()));
      else if (numA == 4)
        return (DenseVect)(prod(*A.sparse(), *x.sparse()));
      else // if(numA==5)
        return (DenseVect)(prod(*A.banded(), *x.sparse()));
    }
  }
}


void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = a*A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += a*A*x if init = false.
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("prod(A,x,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
      y.zero();
    //else nothing
  }

  else if (numA == 7) // A = identity
  {
    scal(a, x, y, init);
  }

  else // A is not 0 or identity
  {

    // === First case: y is not a block vector ===
    {
      {
        if (init)
        {
          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.dense()) = a * ublas::prod(*A.dense(), *x.dense());
              else if (numA == 2)
                noalias(*y.dense()) = a * ublas::prod(*A.triang(), *x.dense());
              else if (numA == 3)
                noalias(*y.dense()) = a * ublas::prod(*A.sym(), *x.dense());
              else if (numA == 4)
                noalias(*y.dense()) = a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==5)
                noalias(*y.dense()) = a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.dense()) = a * ublas::prod(*A.dense(), *x.sparse());
              else if (numA == 2)
                noalias(*y.dense()) = a * ublas::prod(*A.triang(), *x.sparse());
              else if (numA == 3)
                noalias(*y.dense()) = a * ublas::prod(*A.sym(), *x.sparse());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.dense()) = a * ublas::prod(*A.sparse(), *x.sparse());
                else
                  noalias(*y.sparse()) = a * ublas::prod(*A.sparse(), *x.sparse());
              }
              else //if(numA==5)
                noalias(*y.dense()) = a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.dense() = a * ublas::prod(*A.dense(), *x.dense());
              else if (numA == 2)
                *y.dense() = a * ublas::prod(*A.triang(), *x.dense());
              else if (numA == 3)
                *y.dense() = a * ublas::prod(*A.sym(), *x.dense());
              else if (numA == 4)
                *y.dense() = a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==5)
                *y.dense() = a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.sparse() = a * ublas::prod(*A.dense(), *x.sparse());
              else if (numA == 2)
                *y.sparse() = a * ublas::prod(*A.triang(), *x.sparse());
              else if (numA == 3)
                *y.sparse() = a * ublas::prod(*A.sym(), *x.sparse());
              else if (numA == 4)
                *y.sparse() = a * ublas::prod(*A.sparse(), *x.sparse());
              else //if(numA==5)
                *y.sparse() = a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
        }
        else // += case
        {
          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.dense()) += a * ublas::prod(*A.dense(), *x.dense());
              else if (numA == 2)
                noalias(*y.dense()) += a * ublas::prod(*A.triang(), *x.dense());
              else if (numA == 3)
                noalias(*y.dense()) += a * ublas::prod(*A.sym(), *x.dense());
              else if (numA == 4)
                noalias(*y.dense()) += a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==5)
                noalias(*y.dense()) += a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.dense()) += a * ublas::prod(*A.dense(), *x.sparse());
              else if (numA == 2)
                noalias(*y.dense()) += a * ublas::prod(*A.triang(), *x.sparse());
              else if (numA == 3)
                noalias(*y.dense()) += a * ublas::prod(*A.sym(), *x.sparse());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.dense()) += a * ublas::prod(*A.sparse(), *x.sparse());
                else
                  noalias(*y.sparse()) += a * ublas::prod(*A.sparse(), *x.sparse());
              }
              else //if(numA==5)
                noalias(*y.dense()) += a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.dense() += a * ublas::prod(*A.dense(), *x.dense());
              else if (numA == 2)
                *y.dense() += a * ublas::prod(*A.triang(), *x.dense());
              else if (numA == 3)
                *y.dense() += a * ublas::prod(*A.sym(), *x.dense());
              else if (numA == 4)
                *y.dense() += a * ublas::prod(*A.sparse(), *x.dense());
              else //if(numA==5)
                *y.dense() += a * ublas::prod(*A.banded(), *x.dense());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.sparse() += a * ublas::prod(*A.dense(), *x.sparse());
              else if (numA == 2)
                *y.sparse() += a * ublas::prod(*A.triang(), *x.sparse());
              else if (numA == 3)
                *y.sparse() += a * ublas::prod(*A.sym(), *x.sparse());
              else if (numA == 4)
                *y.sparse() += a * ublas::prod(*A.sparse(), *x.sparse());
              else //if(numA==5)
                *y.sparse() += a * ublas::prod(*A.banded(), *x.sparse());
            }
          }
        }
      }
    }
  }
}

void prod(const SiconosVector& x, const SiconosMatrix& A, SiconosVector& y, bool init)
{
  // To compute y = trans(A) * x in an "optimized" way, if init = true
  // (or y = trans(A) * x + y if init = false
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A.size(0) != x.size())
    SiconosMatrixException::selfThrow("prod(x,A,y) error: inconsistent sizes between A and x.");

  if (A.size(1) != y.size())
    SiconosMatrixException::selfThrow("prod(x,A,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("prod(x,A,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
      y.zero();
    // else nothing
  }

  else if (numA == 7) // A = identity
  {
    if (!init)
      y += x;
    else
    {
      if (&x != &y) y = x ; // if x and y do not share memory (ie are different objects)
      // else nothing
    }
  }

  else // A is not 0 or identity
  {
    {
      if (init)
      {

        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              noalias(*y.dense()) = ublas::prod(trans(*A.dense()), *x.dense());
            else if (numA == 2)
              noalias(*y.dense()) = ublas::prod(trans(*A.triang()), *x.dense());
            else if (numA == 3)
              noalias(*y.dense()) = ublas::prod(trans(*A.sym()), *x.dense());
            else if (numA == 4)
              noalias(*y.dense()) = ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==5)
              noalias(*y.dense()) = ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == 4)
          {
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              noalias(*y.dense()) = ublas::prod(trans(*A.dense()), *x.sparse());
            else if (numA == 2)
              noalias(*y.dense()) = ublas::prod(trans(*A.triang()), *x.sparse());
            else if (numA == 3)
              noalias(*y.dense()) = ublas::prod(trans(*A.sym()), *x.sparse());
            else if (numA == 4)
            {
              if (numY == 1)
                noalias(*y.dense()) = ublas::prod(trans(*A.sparse()), *x.sparse());
              else
                noalias(*y.sparse()) = ublas::prod(trans(*A.sparse()), *x.sparse());
            }
            else //if(numA==5)
              noalias(*y.dense()) = ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            if (numA == 1)
              *y.dense() = ublas::prod(trans(*A.dense()), *x.dense());
            else if (numA == 2)
              *y.dense() = ublas::prod(trans(*A.triang()), *x.dense());
            else if (numA == 3)
              *y.dense() = ublas::prod(trans(*A.sym()), *x.dense());
            else if (numA == 4)
              *y.dense() = ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==5)
              *y.dense() = ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == 4)
          {
            if (numA == 1)
              *y.sparse() = ublas::prod(trans(*A.dense()), *x.sparse());
            else if (numA == 2)
              *y.sparse() = ublas::prod(trans(*A.triang()), *x.sparse());
            else if (numA == 3)
              *y.sparse() = ublas::prod(trans(*A.sym()), *x.sparse());
            else if (numA == 4)
              *y.sparse() = ublas::prod(trans(*A.sparse()), *x.sparse());
            else //if(numA==5)
              *y.sparse() = ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
      }
      else // += case
      {

        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              noalias(*y.dense()) += ublas::prod(trans(*A.dense()), *x.dense());
            else if (numA == 2)
              noalias(*y.dense()) += ublas::prod(trans(*A.triang()), *x.dense());
            else if (numA == 3)
              noalias(*y.dense()) += ublas::prod(trans(*A.sym()), *x.dense());
            else if (numA == 4)
              noalias(*y.dense()) += ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==5)
              noalias(*y.dense()) += ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == 4)
          {
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              noalias(*y.dense()) += ublas::prod(trans(*A.dense()), *x.sparse());
            else if (numA == 2)
              noalias(*y.dense()) += ublas::prod(trans(*A.triang()), *x.sparse());
            else if (numA == 3)
              noalias(*y.dense()) += ublas::prod(trans(*A.sym()), *x.sparse());
            else if (numA == 4)
            {
              if (numY == 1)
                noalias(*y.dense()) += ublas::prod(trans(*A.sparse()), *x.sparse());
              else
                noalias(*y.sparse()) += ublas::prod(trans(*A.sparse()), *x.sparse());
            }
            else //if(numA==5)
              noalias(*y.dense()) += ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            if (numA == 1)
              *y.dense() += ublas::prod(trans(*A.dense()), *x.dense());
            else if (numA == 2)
              *y.dense() += ublas::prod(trans(*A.triang()), *x.dense());
            else if (numA == 3)
              *y.dense() += ublas::prod(trans(*A.sym()), *x.dense());
            else if (numA == 4)
              *y.dense() += ublas::prod(trans(*A.sparse()), *x.dense());
            else //if(numA==5)
              *y.dense() += ublas::prod(trans(*A.banded()), *x.dense());
          }
          else //if(numX == 4)
          {
            if (numA == 1)
              *y.sparse() += ublas::prod(trans(*A.dense()), *x.sparse());
            else if (numA == 2)
              *y.sparse() += ublas::prod(trans(*A.triang()), *x.sparse());
            else if (numA == 3)
              *y.sparse() += ublas::prod(trans(*A.sym()), *x.sparse());
            else if (numA == 4)
              *y.sparse() += ublas::prod(trans(*A.sparse()), *x.sparse());
            else //if(numA==5)
              *y.sparse() += ublas::prod(trans(*A.banded()), *x.sparse());
          }
        }
      }
    }
  }
}

void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += A*x if init = false.
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("prod(A,x,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
      y.zero();
    //else nothing
  }

  else if (numA == 7) // A = identity
  {
    if (!init)
      y += x;
    else
    {
      if (&x != &y) y = x ; // if x and y do not share memory (ie are different objects)
      // else nothing
    }
  }

  else // A is not 0 or identity
  {

    // === First case: y is not a block vector ===
    if (init)
    {
      if (&x != &y) // if no common memory between x and y.
      {
        if (numX == 1)
        {
          if (numY != 1)
            SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

          assert(y.dense() != x.dense());

          if (numA == 1)
            noalias(*y.dense()) = ublas::prod(*A.dense(), *x.dense());
          else if (numA == 2)
            noalias(*y.dense()) = ublas::prod(*A.triang(), *x.dense());
          else if (numA == 3)
            noalias(*y.dense()) = ublas::prod(*A.sym(), *x.dense());
          else if (numA == 4)
            noalias(*y.dense()) = ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==5)
            noalias(*y.dense()) = ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == 4)
        {
          if (numY != 1 && numA != 4)
            SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

          if (numA == 1)
            noalias(*y.dense()) = ublas::prod(*A.dense(), *x.sparse());
          else if (numA == 2)
            noalias(*y.dense()) = ublas::prod(*A.triang(), *x.sparse());
          else if (numA == 3)
            noalias(*y.dense()) = ublas::prod(*A.sym(), *x.sparse());
          else if (numA == 4)
          {
            if (numY == 1)
              noalias(*y.dense()) = ublas::prod(*A.sparse(), *x.sparse());
            else
              noalias(*y.sparse()) = ublas::prod(*A.sparse(), *x.sparse());
          }
          else //if(numA==5)
            noalias(*y.dense()) = ublas::prod(*A.banded(), *x.sparse());
        }
      }
      else // if x and y are the same object => alias
      {
        if (numX == 1)
        {
          if (numA == 1)
            *y.dense() = ublas::prod(*A.dense(), *x.dense());
          else if (numA == 2)
            *y.dense() = ublas::prod(*A.triang(), *x.dense());
          else if (numA == 3)
            *y.dense() = ublas::prod(*A.sym(), *x.dense());
          else if (numA == 4)
            *y.dense() = ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==5)
            *y.dense() = ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == 4)
        {
          if (numA == 1)
            *y.sparse() = ublas::prod(*A.dense(), *x.sparse());
          else if (numA == 2)
            *y.sparse() = ublas::prod(*A.triang(), *x.sparse());
          else if (numA == 3)
            *y.sparse() = ublas::prod(*A.sym(), *x.sparse());
          else if (numA == 4)
            *y.sparse() = ublas::prod(*A.sparse(), *x.sparse());
          else //if(numA==5)
            *y.sparse() = ublas::prod(*A.banded(), *x.sparse());
        }
      }
    }
    else // += case
    {
      if (&x != &y) // if no common memory between x and y.
      {
        if (numX == 1)
        {
          if (numY != 1)
            SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

          if (numA == 1)
            noalias(*y.dense()) += ublas::prod(*A.dense(), *x.dense());
          else if (numA == 2)
            noalias(*y.dense()) += ublas::prod(*A.triang(), *x.dense());
          else if (numA == 3)
            noalias(*y.dense()) += ublas::prod(*A.sym(), *x.dense());
          else if (numA == 4)
            noalias(*y.dense()) += ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==5)
            noalias(*y.dense()) += ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == 4)
        {
          if (numY != 1 && numA != 4)
            SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

          if (numA == 1)
            noalias(*y.dense()) += ublas::prod(*A.dense(), *x.sparse());
          else if (numA == 2)
            noalias(*y.dense()) += ublas::prod(*A.triang(), *x.sparse());
          else if (numA == 3)
            noalias(*y.dense()) += ublas::prod(*A.sym(), *x.sparse());
          else if (numA == 4)
          {
            if (numY == 1)
              noalias(*y.dense()) += ublas::prod(*A.sparse(), *x.sparse());
            else
              noalias(*y.sparse()) += ublas::prod(*A.sparse(), *x.sparse());
          }
          else //if(numA==5)
            noalias(*y.dense()) += ublas::prod(*A.banded(), *x.sparse());
        }
      }
      else // if x and y are the same object => alias
      {
        if (numX == 1)
        {
          if (numA == 1)
            *y.dense() += ublas::prod(*A.dense(), *x.dense());
          else if (numA == 2)
            *y.dense() += ublas::prod(*A.triang(), *x.dense());
          else if (numA == 3)
            *y.dense() += ublas::prod(*A.sym(), *x.dense());
          else if (numA == 4)
            *y.dense() += ublas::prod(*A.sparse(), *x.dense());
          else //if(numA==5)
            *y.dense() += ublas::prod(*A.banded(), *x.dense());
        }
        else //if(numX == 4)
        {
          if (numA == 1)
            *y.sparse() += ublas::prod(*A.dense(), *x.sparse());
          else if (numA == 2)
            *y.sparse() += ublas::prod(*A.triang(), *x.sparse());
          else if (numA == 3)
            *y.sparse() += ublas::prod(*A.sym(), *x.sparse());
          else if (numA == 4)
            *y.sparse() += ublas::prod(*A.sparse(), *x.sparse());
          else //if(numA==5)
            *y.sparse() += ublas::prod(*A.banded(), *x.sparse());
        }
      }
    }
  }
}


void axpy_prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = A * x ( init = true) or y += A * x (init = false) using ublas::axpy_prod
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("axpy_prod(A,x,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init) y.zero(); // else nothing ...
  }

  else if (numA == 7) // A = identity
  {
    if (!init) y += x;
    else
    {
      if (&x != &y)
        y = x ; // if x and y do not share memory (ie are different objects)
    }
    // else nothing
  }

  else // A is not 0 or identity
  {
    {
      {
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              ublas::axpy_prod(*A.dense(), *x.dense(), *y.dense(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.triang(), *x.dense(), *y.dense(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.sym(), *x.dense(), *y.dense(), init);
            else if (numA == 4)
              ublas::axpy_prod(*A.sparse(), *x.dense(), *y.dense(), init);
            else //if(numA==5)
              ublas::axpy_prod(*A.banded(), *x.dense(), *y.dense(), init);
          }
          else //if(numX == 4)
          {
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("axpy_prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              ublas::axpy_prod(*A.dense(), *x.sparse(), *y.dense(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.triang(), *x.sparse(), *y.dense(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.sym(), *x.sparse(), *y.dense(), init);
            else if (numA == 4)
            {
              if (numY == 1)
                ublas::axpy_prod(*A.sparse(), *x.sparse(), *y.dense(), init);
              else
                ublas::axpy_prod(*A.sparse(), *x.sparse(), *y.sparse(), init);
            }
            else //if(numA==5)
              ublas::axpy_prod(*A.banded(), *x.sparse(), *y.dense(), init);
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            if (numA == 1)
              ublas::axpy_prod(*A.dense(), *x.dense(), *x.dense(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.triang(), *x.dense(), *x.dense(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.sym(), *x.dense(), *x.dense(), init);
            else if (numA == 4)
              ublas::axpy_prod(*A.sparse(), *x.dense(), *x.dense(), init);
            else //if(numA==5)
              ublas::axpy_prod(*A.banded(), *x.dense(), *x.dense(), init);
          }
          else //if(numX == 4)
          {
            if (numA == 1)
              ublas::axpy_prod(*A.dense(), *x.sparse(), *x.sparse(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.triang(), *x.sparse(), *x.sparse(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.sym(), *x.sparse(), *x.sparse(), init);
            else if (numA == 4)
              ublas::axpy_prod(*A.sparse(), *x.sparse(), *x.sparse(), init);
            else //if(numA==5)
              ublas::axpy_prod(*A.banded(), *x.sparse(), *x.sparse(), init);
          }
        }
      }
    }
  }
}

void gemvtranspose(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  siconosBindings::blas::gemv(a, siconosBindings::trans(*A.dense()), *x.dense(), b, *y.dense());
}

void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  siconosBindings::blas::gemv(a, *A.dense(), *x.dense(), b, *y.dense());
}
