/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "SimpleMatrixFriends.hpp"
#include "SiconosMatrix.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "SiconosAlgebraTools.hpp" // for isComparableTo

using Siconos::Algebra::isComparableTo;

void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool init)
{
  // To compute B = a * A (init = true) or B += a*A (init = false).
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");
  if(!B.isBlock())
    B.resetFactorizationFlags();

  if(&A == &B)
  {
    if(init) B *= a;
    else B *= (1.0 + a);
  }
  else
  {
    unsigned int numA = A.num();
    unsigned int numB = B.num();

    if(numB == Siconos::ZERO || numB == Siconos::IDENTITY)  // B = 0 or identity.
      SiconosMatrixException::selfThrow("scal(a,A,B) : forbidden for B being a zero or identity matrix.");

    if(numA == Siconos::ZERO)
    {
      if(init) B.zero();  // else nothing
    }
    else if(numA == Siconos::IDENTITY)
    {
      if(init)
      {
        B.eye();
        B *= a;
      }
      else
      {
        // Assuming B is square ...
        for(unsigned int i = 0; i < B.size(0); ++i)
          B(i, i) += a;
      }
    }
    else
    {
      if(numA == numB)  // if A and B are of the same type ...
      {
        switch(numA)
        {

        case 0: // A and B are block
          if(isComparableTo(A, B))
          {
            const BlockMatrix& Aref = static_cast<const BlockMatrix&>(A);
            BlockMatrix& Bref = static_cast<BlockMatrix&>(B);
            BlocksIterator1 itB1;
            BlocksIterator2 itB2;
            ConstBlocksIterator1 itA1 = Aref._mat->begin1();
            ConstBlocksIterator2 itA2;
            for(itB1 = Bref._mat->begin1(); itB1 != Bref._mat->end1(); ++itB1)
            {
              itA2 = itA1.begin();
              for(itB2 = itB1.begin(); itB2 != itB1.end(); ++itB2)
              {
                scal(a, **itA2++, **itB2, init);
              }
              itA1++;
            }
          }
          else // if A and B are not "block-consistent"
          {
            if(init)
            {
              for(unsigned int i = 0; i < A.size(0); ++i)
                for(unsigned int j = 0; j < A.size(1); ++j)
                  B(i, j) = a * A(i, j);
            }
            else
            {
              for(unsigned int i = 0; i < A.size(0); ++i)
                for(unsigned int j = 0; j < A.size(1); ++j)
                  B(i, j) += a * A(i, j);
            }
          }
          break;

        case Siconos::DENSE: // if both are dense
          if(init)
            noalias(*B.dense()) = a ** A.dense();
          else
            noalias(*B.dense()) += a ** A.dense();
          break;
        case Siconos::TRIANGULAR:
          if(init)
            noalias(*B.triang()) = a ** A.triang();
          else
            noalias(*B.triang()) += a ** A.triang();
          break;
        case Siconos::SYMMETRIC:
          if(init)
            noalias(*B.sym()) = a ** A.sym();
          else
            noalias(*B.sym()) += a ** A.sym();
          break;
        case Siconos::SPARSE:
          if(init)
            noalias(*B.sparse()) = a ** A.sparse();
          else
            noalias(*B.sparse()) += a ** A.sparse();
          break;
        case Siconos::BANDED:
          if(init)
            noalias(*B.banded()) = a ** A.banded();
          else
            noalias(*B.banded()) += a ** A.banded();
          break;
        }
      }
      else // if A and B are of different types.
      {
        if(numA == 0 || numB == 0)  // if A or B is block
        {
          if(init)
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
          if(numB != Siconos::DENSE)
            SiconosMatrixException::selfThrow("scal(a,A,B) failed. A and B types do not fit together.");

          if(init)
          {
            switch(numB)
            {
            case Siconos::DENSE:
              noalias(*B.dense()) = a ** A.dense();
              break;
            case Siconos::TRIANGULAR:
              noalias(*B.dense()) = a ** A.triang();
              break;
            case Siconos::SYMMETRIC:
              noalias(*B.dense()) = a ** A.sym();
              break;
            case Siconos::SPARSE:
              noalias(*B.dense()) = a ** A.sparse();
              break;
            case Siconos::BANDED:
              noalias(*B.dense()) = a ** A.banded();
              break;
            }
          }
          else

          {
            switch(numB)
            {
            case Siconos::DENSE:
              noalias(*B.dense()) += a ** A.dense();
              break;
            case Siconos::TRIANGULAR:
              noalias(*B.dense()) += a ** A.triang();
              break;
            case Siconos::SYMMETRIC:
              noalias(*B.dense()) += a ** A.sym();
              break;
            case Siconos::SPARSE:
              noalias(*B.dense()) += a ** A.sparse();
              break;
            case Siconos::BANDED:
              noalias(*B.dense()) += a ** A.banded();
              break;
            }
          }
        }
      }
    }
  }
}
