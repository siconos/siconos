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

#include "BlockMatrixIterators.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "SimpleMatrix.hpp"
#include "SiconosAlgebra.hpp"

using namespace siconos;


void add(const SiconosMatrix & A, const SiconosMatrix& B, SiconosMatrix& C)
{
  // To compute C = A + B in an "optimized" way (in comparison with operator +)

  if((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    THROW_EXCEPTION("Matrix addition: inconsistent sizes");
  if((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    THROW_EXCEPTION("Matrix addition: inconsistent sizes");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numB = B.num();
  siconos::UBLAS_TYPE numC = C.num();

  // === if C is zero or identity => read-only ===
  if(numC == siconos::ZERO || numC == siconos::IDENTITY)
    THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");

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
    if(numA == siconos::ZERO)  // A = 0
      C = B ;
    else if(numB == siconos::ZERO)  // B = 0
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
            if(numA == siconos::DENSE)
              noalias(*C.dense()) = *A.dense() + *B.dense();
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.triang()) = *A.triang() + *B.triang();
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.sym()) = *A.sym() + *B.sym();
            else if(numA == siconos::SPARSE)
              noalias(*C.sparse()) = *A.sparse() + *B.sparse();
            else //if(numA==siconos::BANDED)
              noalias(*C.banded()) = *A.banded() + *B.banded();
          }
          else // C and A of different types.
          {
            if(numC != siconos::DENSE)
              THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
            // Only dense matrices are allowed for output.

            if(numA == siconos::DENSE)
              noalias(*C.dense()) = *A.dense() + *B.dense();
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) = *A.triang() + *B.triang();
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense()) = *A.sym() + *B.sym();
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) = *A.sparse() + *B.sparse();
            else //if(numA==siconos::BANDED)
              noalias(*C.dense()) = *A.banded() + *B.banded();
          }
          C.resetFactorizationFlags();
        }
        else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
        {
          if(numC != siconos::DENSE)
            THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if(numA == siconos::DENSE)
            switch(numB)
            {
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.dense() + *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.dense() + *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.dense() + *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.dense() + *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.dense() + *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::TRIANGULAR)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.triang() + *B.dense();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.triang() + *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.triang() + *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.triang() + *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.triang() + *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::SYMMETRIC)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.sym() + *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.sym() + *B.triang();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.sym() + *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.sym() + *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.sym() + *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::SPARSE)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.sparse() + *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.sparse() + *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.sparse() + *B.sym();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.sparse() + *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.sparse() + *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::BANDED)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.banded() + *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.banded() + *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.banded() + *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.banded() + *B.sparse();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.banded() + *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::IDENTITY)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.identity() + *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.identity() + *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.identity() + *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.identity() + *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.identity() + *B.banded();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else
            THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
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
    THROW_EXCEPTION("Matrix addition: inconsistent sizes");
  if((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    THROW_EXCEPTION("Matrix addition: inconsistent sizes");

  siconos::UBLAS_TYPE numA = A.num();
  siconos::UBLAS_TYPE numB = B.num();
  siconos::UBLAS_TYPE numC = C.num();

  // === if C is zero or identity => read-only ===
  if(numC == siconos::ZERO || numC == siconos::IDENTITY)
    THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");

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
          if(numA == siconos::DENSE)
            *C.dense() = *A.dense() - *B.dense();
          else if(numA == siconos::TRIANGULAR)
            *C.triang() = *A.triang() - *B.triang();
          else if(numA == siconos::SYMMETRIC)
            *C.sym() = *A.sym() - *B.sym();
          else if(numA == siconos::SPARSE)
            *C.sparse() = *A.sparse() - *B.sparse();
          else //if(numA==siconos::BANDED)
            *C.banded() = *A.banded() - *B.banded();
        }
        else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
        {
          if(numC != siconos::DENSE)   // => numB == siconos::DENSE
            THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if(numA == siconos::DENSE)
            *C.dense() = *A.dense() - *B.dense();
          else if(numA == siconos::TRIANGULAR)
            *C.dense() = *A.triang() - *B.dense();
          else if(numA == siconos::SYMMETRIC)
            *C.dense() = *A.sym() - *B.dense();
          else if(numA == siconos::SPARSE)
            *C.dense() = *A.sparse() - *B.dense();
          else if(numA == siconos::BANDED)
            *C.dense() = *A.banded() - *B.dense();
          else if(numA == siconos::ZERO)
            *C.dense() = *A.zero_mat() - *B.dense();
          else //if(numA==siconos::IDENTITY)
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
    if(numB == siconos::ZERO)  // B = 0
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
            if(numA == siconos::DENSE)
              noalias(*C.dense()) = *A.dense() - *B.dense();
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.triang()) = *A.triang() - *B.triang();
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.sym()) = *A.sym() - *B.sym();
            else if(numA == siconos::SPARSE)
              noalias(*C.sparse()) = *A.sparse() - *B.sparse();
            else //if(numA==siconos::BANDED)
              noalias(*C.banded()) = *A.banded() - *B.banded();
          }
          else // C and A of different types.
          {
            if(numC != siconos::DENSE)
              THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
            // Only dense matrices are allowed for output.

            if(numA == siconos::DENSE)
              noalias(*C.dense()) = *A.dense() - *B.dense();
            else if(numA == siconos::TRIANGULAR)
              noalias(*C.dense()) = *A.triang() - *B.triang();
            else if(numA == siconos::SYMMETRIC)
              noalias(*C.dense()) = *A.sym() - *B.sym();
            else if(numA == siconos::SPARSE)
              noalias(*C.dense()) = *A.sparse() - *B.sparse();
            else //if(numA==siconos::BANDED)
              noalias(*C.dense()) = *A.banded() - *B.banded();
          }
          C.resetFactorizationFlags();
        }
        else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
        {
          if(numC != siconos::DENSE)
            THROW_EXCEPTION("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if(numA == siconos::DENSE)
            switch(numB)
            {
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.dense() - *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.dense() - *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.dense() - *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.dense() - *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.dense() - *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::TRIANGULAR)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.triang() - *B.dense();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.triang() - *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.triang() - *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.triang() - *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.triang() - *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::SYMMETRIC)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.sym() - *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.sym() - *B.triang();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.sym() - *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.sym() - *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.sym() - *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::SPARSE)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.sparse() - *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.sparse() - *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.sparse() - *B.sym();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.sparse() - *B.banded();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.sparse() - *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::BANDED)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.banded() - *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.banded() - *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.banded() - *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.banded() - *B.sparse();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.banded() - *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::ZERO)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.zero_mat() - *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.zero_mat() - *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.zero_mat() - *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.zero_mat() - *B.sparse();
              break;
            case siconos::IDENTITY:
              noalias(*C.dense()) = *A.zero_mat() - *B.identity();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if(numA == siconos::IDENTITY)
            switch(numB)
            {
            case siconos::DENSE:
              noalias(*C.dense()) = *A.identity() - *B.dense();
              break;
            case siconos::TRIANGULAR:
              noalias(*C.dense()) = *A.identity() - *B.triang();
              break;
            case siconos::SYMMETRIC:
              noalias(*C.dense()) = *A.identity() - *B.sym();
              break;
            case siconos::SPARSE:
              noalias(*C.dense()) = *A.identity() - *B.sparse();
              break;
            case siconos::BANDED:
              noalias(*C.dense()) = *A.identity() - *B.banded();
              break;
            default:
              THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
            }
          else
            THROW_EXCEPTION("Matrix function add(A,B,C): invalid type of matrix");
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


