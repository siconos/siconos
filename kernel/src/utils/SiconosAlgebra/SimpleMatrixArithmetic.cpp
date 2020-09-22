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

#include "BlockMatrixIterators.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "SimpleMatrix.hpp"
#include "SiconosAlgebra.hpp"

using namespace Siconos;


void add(const SiconosMatrix & A, const SiconosMatrix& B, SiconosMatrix& C)
{
  // To compute C = A + B in an "optimized" way (in comparison with operator +)

  if((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");
  if((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();
  unsigned int numC = C.num();

  // === if C is zero or identity => read-only ===
  if(numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");

  // === common memory between A, B, C ===
  if(&A == &C)  // A and C have common memory
  {
    C += B;
  }
  else if(&B == &C)   // B and C have common memory
  {
    C += A;
  }
  else // No common memory between C and A or B.
  {
    if(numA == 6)  // A = 0
      C = B ;
    else if(numB == 6)  // B = 0
      C = A;
    else // A and B different from 0
    {
      if(numC == 0)  // if C is Block
      {
        if(numA != 0)  // A simple, whatever is B
        {
          C = A;
          C += B;
        }
        else  // A Block
        {
          C = B;
          C += A;
        }
      }
      else // if C is a SimpleMatrix
      {
        if(numA == numB && numA != 0)  // A and B are of the same type and NOT block
        {
          if(numC == numA)
          {
            if(numA == 1)
              noalias(*C.dense()) = *A.dense() + *B.dense();
            else if(numA == 2)
              noalias(*C.triang()) = *A.triang() + *B.triang();
            else if(numA == 3)
              noalias(*C.sym()) = *A.sym() + *B.sym();
            else if(numA == 4)
              noalias(*C.sparse()) = *A.sparse() + *B.sparse();
            else //if(numA==5)
              noalias(*C.banded()) = *A.banded() + *B.banded();
          }
          else // C and A of different types.
          {
            if(numC != 1)
              SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
            // Only dense matrices are allowed for output.

            if(numA == 1)
              noalias(*C.dense()) = *A.dense() + *B.dense();
            else if(numA == 2)
              noalias(*C.dense()) = *A.triang() + *B.triang();
            else if(numA == 3)
              noalias(*C.dense()) = *A.sym() + *B.sym();
            else if(numA == 4)
              noalias(*C.dense()) = *A.sparse() + *B.sparse();
            else //if(numA==5)
              noalias(*C.dense()) = *A.banded() + *B.banded();
          }
          C.resetFactorizationFlags();
        }
        else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
        {
          if(numC != 1)
            SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if(numA == 1)
            switch(numB)
            {
            case 2:
              noalias(*C.dense()) = *A.dense() + *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.dense() + *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.dense() + *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.dense() + *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.dense() + *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 2)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.triang() + *B.dense();
              break;
            case 3:
              noalias(*C.dense()) = *A.triang() + *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.triang() + *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.triang() + *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.triang() + *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 3)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.sym() + *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.sym() + *B.triang();
              break;
            case 4:
              noalias(*C.dense()) = *A.sym() + *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.sym() + *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.sym() + *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 4)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.sparse() + *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.sparse() + *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.sparse() + *B.sym();
              break;
            case 5:
              noalias(*C.dense()) = *A.sparse() + *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.sparse() + *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 5)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.banded() + *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.banded() + *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.banded() + *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.banded() + *B.sparse();
              break;
            case 7:
              noalias(*C.dense()) = *A.banded() + *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 7)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.identity() + *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.identity() + *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.identity() + *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.identity() + *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.identity() + *B.banded();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          C.resetFactorizationFlags();
        }
        else // A and/or B is Block
        {
          if(numA != 0)  // A Simple, whatever is B
          {
            C = A;
            C += B;
          }
          else // A Block
          {
            C = B;
            C += A;
          }
        }
      }
    }
  }

}

void sub(const SiconosMatrix & A, const SiconosMatrix& B, SiconosMatrix& C)
{
  // To compute C = A - B in an "optimized" way (in comparison with operator +)

  if((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");
  if((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();
  unsigned int numC = C.num();

  // === if C is zero or identity => read-only ===
  if(numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");

  // === common memory between A, B, C ===
  if(&A == &C)  // A and C have common memory
  {
    C -= B;
  }
  else if(&B == &C)   // B and C have common memory
  {
    if(numB == 0 || numA == 0)  // if A or B(C) is Block
    {
      C *= -1.0;
      C += A;
    }
    else
    {
      if(numC == 0)  // if C is Block
      {
        C = A;
        C -= B;
      }
      else // if C is a SimpleMatrix
      {
        if(numA == numB && numA != 0)  // A and B are of the same type and NOT block
        {
          if(numA == 1)
            *C.dense() = *A.dense() - *B.dense();
          else if(numA == 2)
            *C.triang() = *A.triang() - *B.triang();
          else if(numA == 3)
            *C.sym() = *A.sym() - *B.sym();
          else if(numA == 4)
            *C.sparse() = *A.sparse() - *B.sparse();
          else //if(numA==5)
            *C.banded() = *A.banded() - *B.banded();
        }
        else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
        {
          if(numC != 1)   // => numB == 1
            SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if(numA == 1)
            *C.dense() = *A.dense() - *B.dense();
          else if(numA == 2)
            *C.dense() = *A.triang() - *B.dense();
          else if(numA == 3)
            *C.dense() = *A.sym() - *B.dense();
          else if(numA == 4)
            *C.dense() = *A.sparse() - *B.dense();
          else if(numA == 5)
            *C.dense() = *A.banded() - *B.dense();
          else if(numA == 6)
            *C.dense() = *A.zero_mat() - *B.dense();
          else //if(numA==7)
            *C.dense() = *A.identity() - *B.dense();
        }
        else // A and/or B is Block
        {
          C = A;
          C -= B;
        }
        C.resetFactorizationFlags();
      }
    }
  }
  else // No common memory between C and A or B.
  {
    if(numB == 6)  // B = 0
      C = A;
    else // B different from 0
    {
      if(numC == 0)  // if C is Block
      {
        C = A;
        C -= B;
      }
      else // if C is a SimpleMatrix
      {
        if(numA == numB && numA != 0)  // A and B are of the same type and NOT block
        {
          if(numC == numA)
          {
            if(numA == 1)
              noalias(*C.dense()) = *A.dense() - *B.dense();
            else if(numA == 2)
              noalias(*C.triang()) = *A.triang() - *B.triang();
            else if(numA == 3)
              noalias(*C.sym()) = *A.sym() - *B.sym();
            else if(numA == 4)
              noalias(*C.sparse()) = *A.sparse() - *B.sparse();
            else //if(numA==5)
              noalias(*C.banded()) = *A.banded() - *B.banded();
          }
          else // C and A of different types.
          {
            if(numC != 1)
              SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
            // Only dense matrices are allowed for output.

            if(numA == 1)
              noalias(*C.dense()) = *A.dense() - *B.dense();
            else if(numA == 2)
              noalias(*C.dense()) = *A.triang() - *B.triang();
            else if(numA == 3)
              noalias(*C.dense()) = *A.sym() - *B.sym();
            else if(numA == 4)
              noalias(*C.dense()) = *A.sparse() - *B.sparse();
            else //if(numA==5)
              noalias(*C.dense()) = *A.banded() - *B.banded();
          }
          C.resetFactorizationFlags();
        }
        else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
        {
          if(numC != 1)
            SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if(numA == 1)
            switch(numB)
            {
            case 2:
              noalias(*C.dense()) = *A.dense() - *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.dense() - *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.dense() - *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.dense() - *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.dense() - *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 2)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.triang() - *B.dense();
              break;
            case 3:
              noalias(*C.dense()) = *A.triang() - *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.triang() - *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.triang() - *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.triang() - *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 3)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.sym() - *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.sym() - *B.triang();
              break;
            case 4:
              noalias(*C.dense()) = *A.sym() - *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.sym() - *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.sym() - *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 4)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.sparse() - *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.sparse() - *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.sparse() - *B.sym();
              break;
            case 5:
              noalias(*C.dense()) = *A.sparse() - *B.banded();
              break;
            case 7:
              noalias(*C.dense()) = *A.sparse() - *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 5)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.banded() - *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.banded() - *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.banded() - *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.banded() - *B.sparse();
              break;
            case 7:
              noalias(*C.dense()) = *A.banded() - *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 6)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.zero_mat() - *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.zero_mat() - *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.zero_mat() - *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.zero_mat() - *B.sparse();
              break;
            case 7:
              noalias(*C.dense()) = *A.zero_mat() - *B.identity();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == 7)
            switch(numB)
            {
            case 1:
              noalias(*C.dense()) = *A.identity() - *B.dense();
              break;
            case 2:
              noalias(*C.dense()) = *A.identity() - *B.triang();
              break;
            case 3:
              noalias(*C.dense()) = *A.identity() - *B.sym();
              break;
            case 4:
              noalias(*C.dense()) = *A.identity() - *B.sparse();
              break;
            case 5:
              noalias(*C.dense()) = *A.identity() - *B.banded();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          C.resetFactorizationFlags();
        }
        else // A and/or B is Block
        {
          C = A;
          C -= B;
        }
      }
    }
  }
}


