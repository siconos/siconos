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

#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"

#include "SiconosAlgebra.hpp"
#include "SimpleMatrixFriends.hpp" // for operators
using namespace Siconos;
//#define DEBUG_MESSAGES
#include "debug.h"

SP::SimpleMatrix operator * (const SP::SimpleMatrix A, const SP::SimpleMatrix B)
{
  SP::SimpleMatrix aux(new SimpleMatrix((DenseMat)prod(*(*A).dense(), *(*B).dense())));
  return aux;
}

const SimpleMatrix operator * (const SiconosMatrix & A, double a)
{
  // To compute B = a * A

  unsigned int numA = A.num();

  if(numA == ZERO)  // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if(numA == IDENTITY)
  {
    return (DenseMat)(a**A.identity());
  }
  else if(numA == 0)  // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp *= a;
    return tmp;
  }
  else if(numA == DENSE)  // dense)
    return (DenseMat)(a** A.dense());
  else if(numA == TRIANGULAR)
    return (TriangMat)(a ** A.triang());
  else if(numA == SYMMETRIC)
    return (SymMat)(a ** A.sym());
  else if(numA == SPARSE)
    return (SparseMat)(a ** A.sparse());
  else if(numA == BANDED)
    return (BandedMat)(a ** A.banded());
  else
  {
    SiconosMatrixException::selfThrow("SimpleMatrix::op* (const SimpleMatrix): invalid type of matrix");
    return 0;
  }
}

SimpleMatrix operator * (double a, const SiconosMatrix & A)
{
  // To compute B = a * A

  unsigned int numA = A.num();

  if(numA == ZERO)  // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if(numA == IDENTITY)
  {
    return (DenseMat)(a**A.identity());
  }
  else if(numA == 0)  // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp *= a;
    return tmp;
  }
  else if(numA == DENSE)  // dense)
    return (DenseMat)(a** A.dense());
  else if(numA == TRIANGULAR)
    return (TriangMat)(a ** A.triang());
  else if(numA == SYMMETRIC)
    return (SymMat)(a ** A.sym());
  else if(numA == SPARSE)
    return (SparseMat)(a ** A.sparse());
  else if(numA == SPARSE_COORDINATE)
    return (SparseCoordinateMat)(a ** A.sparseCoordinate());
  else if(numA == BANDED)
    return (BandedMat)(a ** A.banded());
  else
  {
    SiconosMatrixException::selfThrow("SimpleMatrix::op* (const SimpleMatrix): invalid type of matrix");
    return 0;
  }
}

const SimpleMatrix operator / (const SiconosMatrix & A, double a)
{
  // To compute B = A/a

  if(a == 0.0)
    SiconosMatrixException::selfThrow(" Matrix, operator / , division by zero.");

  unsigned int numA = A.num();

  if(numA == ZERO)  // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if(numA == IDENTITY)
  {
    return (DenseMat)(*A.identity() / a);
  }
  else if(numA == 0)  // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp /= a;
    return tmp;
  }
  else if(numA == DENSE)  // dense)
    return (DenseMat)(*A.dense() / a);
  else if(numA == TRIANGULAR)
    return (TriangMat)(*A.triang() / a);
  else if(numA == SYMMETRIC)
    return (SymMat)(*A.sym() / a);
  else if(numA == SPARSE)
    return (SparseMat)(*A.sparse() / a);
  else if(numA == BANDED)
    return (BandedMat)(*A.banded() / a);
  else
  {
    SiconosMatrixException::selfThrow("SimpleMatrix::op / (const SimpleMatrix): invalid type of matrix");
    return 0;
  }
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

  if((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix operator + : inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();

  // == A or B equal to null ==
  if(numA == ZERO)  // A = 0
  {
    if(numB == ZERO)  // B = 0
      return SimpleMatrix(A.size(0), A.size(1));
    else
      return SimpleMatrix(B);
  }

  if(numB == ZERO)
    return SimpleMatrix(A);

  // == A and B different from 0 ==

  if(numA == numB && numA != 0)  // all matrices are of the same type and NOT block
  {
    if(numA == DENSE)
      return (DenseMat)(*A.dense() + *B.dense());
    else if(numA == TRIANGULAR)
      return (TriangMat)(*A.triang() + *B.triang());
    else if(numA == SYMMETRIC)
      return (SymMat)(*A.sym() + *B.sym());
    else if(numA == SPARSE)
    {
      SparseMat tmp(*A.sparse());
      tmp += *B.sparse();
      return tmp;
      // return (SparseMat)(*A.sparse() + *B.sparse());
    }
    else if(numA == SPARSE_COORDINATE)
    {
      SparseMat tmp(*A.sparseCoordinate());
      tmp += *B.sparseCoordinate();
      return tmp;
    }
    else if(numA == BANDED)
    {
      BandedMat tmp(*A.banded());
      tmp += *B.banded();
      return tmp;
    }
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
  }
  else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
  {
    if(numA == DENSE)
    {
      if(numB == TRIANGULAR)
        return (DenseMat)(*A.dense() + *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.dense() + *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.dense() + *B.sparse());
      else if(numB == SPARSE_COORDINATE)
        return (DenseMat)(*A.dense() + *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.dense() + *B.banded());
      else if(numB == IDENTITY)
        return (DenseMat)(*A.dense() + *B.identity());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
    }
    else if(numA == TRIANGULAR)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.triang() + *B.dense());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.triang() + *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.triang() + *B.sparse());
      else if(numB == SPARSE_COORDINATE)
        return (DenseMat)(*A.triang() + *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.triang() + *B.banded());
      else if(numB == IDENTITY)
        return (DenseMat)(*A.triang() + *B.identity());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
    }
    else if(numA == SYMMETRIC)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.sym() + *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.sym() + *B.triang());
      else if(numB == SPARSE)
        return (DenseMat)(*A.sym() + *B.sparse());
      else if(numB == SPARSE_COORDINATE)
        return (DenseMat)(*A.sym() + *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.sym() + *B.banded());
      else if(numB == IDENTITY)
        return (DenseMat)(*A.sym() + *B.identity());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
    }
    else if(numA == SPARSE)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.sparse() + *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.sparse() + *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.sparse() + *B.sym());
      else if(numB == BANDED)
        return (DenseMat)(*A.sparse() + *B.banded());
      else if(numB ==IDENTITY)
        return (DenseMat)(*A.sparse() + *B.identity());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
    }

    else if(numA == BANDED)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.banded() + *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.banded() + *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.banded() + *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.banded() + *B.sparse());
      else if(numB == SPARSE_COORDINATE)
        return (DenseMat)(*A.banded() + *B.sparseCoordinate());
      else if(numB ==IDENTITY)
        return (DenseMat)(*A.banded() + *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
      }
    }

    else if(numA == IDENTITY)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.identity() + *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.identity() + *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.identity() + *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.identity() + *B.sparse());
      else if(numB == SPARSE_COORDINATE)
        return (DenseMat)(*A.identity() + *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.identity() + *B.banded());
    }
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
  }
  else if(numB != 0)  // B Simple, whatever is A
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
  SiconosMatrixException::selfThrow("SimpleMatrix::op + (const SimpleMatrix): invalid type of matrix");
  return 0;
}

const SimpleMatrix operator - (const  SiconosMatrix& A, const  SiconosMatrix& B)
{
  // To compute C = A - B

  if((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix operator -  : inconsistent sizes");

  unsigned int numA = A.num();
  unsigned int numB = B.num();


  // == B equal to null ==
  if(numB == ZERO)
    return SimpleMatrix(A);

  // == B different from 0 ==

  if(numA == numB && numA != 0)  // all matrices are of the same type and NOT block
  {
    if(numA == DENSE)
      return (DenseMat)(*A.dense() - *B.dense());
    else if(numA == TRIANGULAR)
      return (TriangMat)(*A.triang() - *B.triang());
    else if(numA == SYMMETRIC)
      return (SymMat)(*A.sym() - *B.sym());
    else if(numA == SPARSE)
    {
      SparseMat tmp(*A.sparse());
      tmp -= *B.sparse();
      return tmp;
      //return (SparseMat)(*A.sparse() - *B.sparse());
    }
    else if(numA==SPARSE_COORDINATE)
    {
      SparseCoordinateMat tmp(*A.sparseCoordinate());
      tmp -= *B.sparseCoordinate();
      return tmp;
    }
    else if(numA == BANDED)
    {
      BandedMat tmp(*A.banded());
      tmp -= *B.banded();
      return tmp;
      //return (BandedMat)(*A.banded() - *B.banded());
    }
    else
    {
      SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
      return 0;
    }
  }
  else if(numA != 0 && numB != 0 && numA != numB)  // A and B of different types and none is block
  {
    if(numA == DENSE)
    {
      if(numB == TRIANGULAR)
        return (DenseMat)(*A.dense() - *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.dense() - *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.dense() - *B.sparse());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.dense() - *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.dense() - *B.banded());
      else if(numB == IDENTITY)
        return (DenseMat)(*A.dense() - *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }
    else if(numA == TRIANGULAR)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.triang() - *B.dense());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.triang() - *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.triang() - *B.sparse());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.triang() - *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.triang() - *B.banded());
      else  if(numB == IDENTITY)
        return (DenseMat)(*A.triang() - *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }
    else if(numA == SYMMETRIC)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.sym() - *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.sym() - *B.triang());
      else if(numB == SPARSE)
        return (DenseMat)(*A.sym() - *B.sparse());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.sym() - *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.sym() - *B.banded());
      else  if(numB == IDENTITY)
        return (DenseMat)(*A.sym() - *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }
    else if(numA == SPARSE)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.sparse() - *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.sparse() - *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.sparse() - *B.sym());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.sparse() - *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.sparse() - *B.banded());
      else  if(numB == IDENTITY)
        return (DenseMat)(*A.sparse() - *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }

    else if(numA == BANDED)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.banded() - *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.banded() - *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.banded() - *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.banded() - *B.sparse());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.banded() - *B.sparseCoordinate());
      else if(numB == IDENTITY)
        return (DenseMat)(*A.banded() - *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }

    else if(numA == ZERO)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.zero_mat() - *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.zero_mat() - *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.zero_mat() - *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.zero_mat() - *B.sparse());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.zero_mat() - *B.sparseCoordinate());
      else if(numB == IDENTITY)
        return (DenseMat)(*A.zero_mat() - *B.identity());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }
    else if(numA == IDENTITY)
    {
      if(numB == DENSE)
        return (DenseMat)(*A.identity() - *B.dense());
      else if(numB == TRIANGULAR)
        return (DenseMat)(*A.identity() - *B.triang());
      else if(numB == SYMMETRIC)
        return (DenseMat)(*A.identity() - *B.sym());
      else if(numB == SPARSE)
        return (DenseMat)(*A.identity() - *B.sparse());
      else if(numB ==SPARSE_COORDINATE)
        return (DenseMat)(*A.identity() - *B.sparseCoordinate());
      else if(numB == BANDED)
        return (DenseMat)(*A.identity() - *B.banded());
      else
      {
        SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
        return 0;
      }
    }
    else
    {
      SiconosMatrixException::selfThrow("SimpleMatrix::op - (const SimpleMatrix): invalid type of matrix");
      return 0;
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
  DEBUG_EXPR((m - x).display());
  DEBUG_PRINTF("norm = %12.8e \n", norm);
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n", std::numeric_limits<double>::epsilon());
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
