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

#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"

using namespace Siconos;
// #define DEBUG_MESSAGES
#include "debug.h"

SP::SimpleMatrix operator * (const SP::SimpleMatrix A, const SP::SimpleMatrix B)
{
  SP::SimpleMatrix aux(new SimpleMatrix((DenseMat)prod(*(*A).dense() , *(*B).dense())));
  return aux;
}

const SimpleMatrix operator * (const SiconosMatrix & A, double a)
{
  // To compute B = a * A

  unsigned int numA = A.num();

  if (numA == 6) // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if (numA == 7)
  {
    return (DenseMat)(a**A.identity());
  }
  else if (numA == 0) // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp *= a;
    return tmp;
  }
  else if (numA == 1) // dense)
    return (DenseMat)(a** A.dense());
  else if (numA == 2)
    return (TriangMat)(a ** A.triang());
  else if (numA == 3)
    return (SymMat)(a ** A.sym());
  else if (numA == 4)
    return (SparseMat)(a ** A.sparse());
  else //if(numA==5)
    return (BandedMat)(a ** A.banded());
}

SimpleMatrix operator * (double a, const SiconosMatrix & A)
{
  // To compute B = a * A

  unsigned int numA = A.num();

  if (numA == 6) // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if (numA == 7)
  {
    return (DenseMat)(a**A.identity());
  }
  else if (numA == 0) // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp *= a;
    return tmp;
  }
  else if (numA == 1) // dense)
    return (DenseMat)(a** A.dense());
  else if (numA == 2)
    return (TriangMat)(a ** A.triang());
  else if (numA == 3)
    return (SymMat)(a ** A.sym());
  else if (numA == 4)
    return (SparseMat)(a ** A.sparse());
  else //if(numA==5)
    return (BandedMat)(a ** A.banded());
}

const SimpleMatrix operator / (const SiconosMatrix & A, double a)
{
  // To compute B = A/a

  if (a == 0.0)
    SiconosMatrixException::selfThrow(" Matrix, operator / , division by zero.");

  unsigned int numA = A.num();

  if (numA == 6) // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if (numA == 7)
  {
    return (DenseMat)(*A.identity() / a);
  }
  else if (numA == 0) // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp /= a;
    return tmp;
  }
  else if (numA == 1) // dense)
    return (DenseMat)(*A.dense() / a);
  else if (numA == 2)
    return (TriangMat)(*A.triang() / a);
  else if (numA == 3)
    return (SymMat)(*A.sym() / a);
  else if (numA == 4)
    return (SparseMat)(*A.sparse() / a);
  else //if(numA==5)
    return (BandedMat)(*A.banded() / a);
}

// const SimpleMatrix operator + (const  SimpleMatrix& A, const  SimpleMatrix& B){
//   return (DenseMat)(*A.dense() + *B.dense());
// }
SimpleMatrix operator + (const  SimpleMatrix& A, const  SimpleMatrix& B)
{

  return (DenseMat)(*A.dense() + *B.dense());
}

void operator +=(SP::SiconosMatrix A, SP::SimpleMatrix B)
{
  *A += *B;
}


SP::SimpleMatrix operator +(const SP::SimpleMatrix A, const SP::SimpleMatrix B)
{
  return SP::SimpleMatrix(new SimpleMatrix(*A + *B));
}



const SimpleMatrix operator + (const  SiconosMatrix& A, const  SiconosMatrix& B)
{
  // To compute C = A + B

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix operator + : inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();

  // == A or B equal to null ==
  if (numA == 6) // A = 0
  {
    if (numB == 6) // B = 0
      return SimpleMatrix(A.size(0), A.size(1));
    else
      return SimpleMatrix(B);
  }

  if (numB == 6)
    return SimpleMatrix(A);

  // == A and B different from 0 ==

  if (numA == numB && numA != 0) // all matrices are of the same type and NOT block
  {
    if (numA == 1)
      return (DenseMat)(*A.dense() + *B.dense());
    else if (numA == 2)
      return (TriangMat)(*A.triang() + *B.triang());
    else if (numA == 3)
      return (SymMat)(*A.sym() + *B.sym());
    else if (numA == 4)
    {
      SparseMat tmp(*A.sparse());
      tmp += *B.sparse();
      return tmp;
      // return (SparseMat)(*A.sparse() + *B.sparse());
    }
    else //if(numA==5)
    {
      BandedMat tmp(*A.banded());
      tmp += *B.banded();
      return tmp;
    }
  }
  else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
  {
    if (numA == 1)
    {
      if (numB == 2)
        return (DenseMat)(*A.dense() + *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.dense() + *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.dense() + *B.sparse());
      else if (numB == 5)
        return (DenseMat)(*A.dense() + *B.banded());
      else // if(numB ==7)
        return (DenseMat)(*A.dense() + *B.identity());
    }
    else if (numA == 2)
    {
      if (numB == 1)
        return (DenseMat)(*A.triang() + *B.dense());
      else if (numB == 3)
        return (DenseMat)(*A.triang() + *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.triang() + *B.sparse());
      else if (numB == 5)
        return (DenseMat)(*A.triang() + *B.banded());
      else // if(numB ==7:
        return (DenseMat)(*A.triang() + *B.identity());
    }
    else if (numA == 3)
    {
      if (numB == 1)
        return (DenseMat)(*A.sym() + *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.sym() + *B.triang());
      else if (numB == 4)
        return (DenseMat)(*A.sym() + *B.sparse());
      else if (numB == 5)
        return (DenseMat)(*A.sym() + *B.banded());
      else // if(numB ==7)
        return (DenseMat)(*A.sym() + *B.identity());
    }
    else if (numA == 4)
    {
      if (numB == 1)
        return (DenseMat)(*A.sparse() + *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.sparse() + *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.sparse() + *B.sym());
      else if (numB == 5)
        return (DenseMat)(*A.sparse() + *B.banded());
      else // if(numB ==7)
        return (DenseMat)(*A.sparse() + *B.identity());
    }

    else if (numA == 5)
    {
      if (numB == 1)
        return (DenseMat)(*A.banded() + *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.banded() + *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.banded() + *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.banded() + *B.sparse());
      else //if(numB ==7)
        return (DenseMat)(*A.banded() + *B.identity());
    }

    else //if(numA==7)
    {
      if (numB == 1)
        return (DenseMat)(*A.identity() + *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.identity() + *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.identity() + *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.identity() + *B.sparse());
      else //if(numB ==5)
        return (DenseMat)(*A.identity() + *B.banded());
    }
  }
  else if (numB != 0) // B Simple, whatever is A
  {
    SimpleMatrix tmp(B);
    tmp += A;
    return tmp;
  }
  else // B Block, A simple or block
  {
    SimpleMatrix tmp(A);
    tmp += B;
    return tmp;
  }
}

const SimpleMatrix operator - (const  SiconosMatrix& A, const  SiconosMatrix& B)
{
  // To compute C = A - B

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix operator -  : inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();


  // == B equal to null ==
  if (numB == 6)
    return SimpleMatrix(A);

  // == B different from 0 ==

  if (numA == numB && numA != 0) // all matrices are of the same type and NOT block
  {
    if (numA == 1)
      return (DenseMat)(*A.dense() - *B.dense());
    else if (numA == 2)
      return (TriangMat)(*A.triang() - *B.triang());
    else if (numA == 3)
      return (SymMat)(*A.sym() - *B.sym());
    else if (numA == 4)
    {
      SparseMat tmp(*A.sparse());
      tmp -= *B.sparse();
      return tmp;
      //return (SparseMat)(*A.sparse() - *B.sparse());
    }
    else //if(numA==5)
    {
      BandedMat tmp(*A.banded());
      tmp -= *B.banded();
      return tmp;
      //return (BandedMat)(*A.banded() - *B.banded());
    }
  }
  else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
  {
    if (numA == 1)
    {
      if (numB == 2)
        return (DenseMat)(*A.dense() - *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.dense() - *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.dense() - *B.sparse());
      else if (numB == 5)
        return (DenseMat)(*A.dense() - *B.banded());
      else // if(numB ==7)
        return (DenseMat)(*A.dense() - *B.identity());
    }
    else if (numA == 2)
    {
      if (numB == 1)
        return (DenseMat)(*A.triang() - *B.dense());
      else if (numB == 3)
        return (DenseMat)(*A.triang() - *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.triang() - *B.sparse());
      else if (numB == 5)
        return (DenseMat)(*A.triang() - *B.banded());
      else // if(numB ==7:
        return (DenseMat)(*A.triang() - *B.identity());
    }
    else if (numA == 3)
    {
      if (numB == 1)
        return (DenseMat)(*A.sym() - *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.sym() - *B.triang());
      else if (numB == 4)
        return (DenseMat)(*A.sym() - *B.sparse());
      else if (numB == 5)
        return (DenseMat)(*A.sym() - *B.banded());
      else // if(numB ==7)
        return (DenseMat)(*A.sym() - *B.identity());
    }
    else if (numA == 4)
    {
      if (numB == 1)
        return (DenseMat)(*A.sparse() - *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.sparse() - *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.sparse() - *B.sym());
      else if (numB == 5)
        return (DenseMat)(*A.sparse() - *B.banded());
      else // if(numB ==7)
        return (DenseMat)(*A.sparse() - *B.identity());
    }

    else if (numA == 5)
    {
      if (numB == 1)
        return (DenseMat)(*A.banded() - *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.banded() - *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.banded() - *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.banded() - *B.sparse());
      else //if(numB ==7)
        return (DenseMat)(*A.banded() - *B.identity());
    }

    else if (numA == 6)
    {
      if (numB == 1)
        return (DenseMat)(*A.zero_mat() - *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.zero_mat() - *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.zero_mat() - *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.zero_mat() - *B.sparse());
      else //if(numB ==7)
        return (DenseMat)(*A.zero_mat() - *B.identity());
    }
    else //if(numA==7)
    {
      if (numB == 1)
        return (DenseMat)(*A.identity() - *B.dense());
      else if (numB == 2)
        return (DenseMat)(*A.identity() - *B.triang());
      else if (numB == 3)
        return (DenseMat)(*A.identity() - *B.sym());
      else if (numB == 4)
        return (DenseMat)(*A.identity() - *B.sparse());
      else //if(numB ==5)
        return (DenseMat)(*A.identity() - *B.banded());
    }
  }
  else // A and/or B are/is Block
  {
    SimpleMatrix tmp(A);
    tmp -= B;
    return tmp;
  }


}

//========================
// Matrices comparison
//========================

bool operator == (const SiconosMatrix &m, const SiconosMatrix &x)
{
  //  if( ! isComparableTo( m, x))
  //    return false;
  // Warning : two block matrices may be "equal" but have blocks of different sizes.
  double norm = (m - x).normInf();
  DEBUG_PRINTF("norm = %12.8e \n", norm );
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n", std::numeric_limits<double>::epsilon() );
  DEBUG_EXPR(std::cout << std::boolalpha << (norm <= std::numeric_limits<double>::epsilon()) <<std::endl;);
  double atol = 1e-14;
  double rtol = std::numeric_limits<double>::epsilon();
  return (norm <= atol + rtol * x.normInf()) ;
}

bool operator != (const SiconosMatrix &m, const SiconosMatrix &x)
{
  double norm = (m - x).normInf();
  return (norm > std::numeric_limits<double>::epsilon());
}


