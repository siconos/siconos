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

// #include "BlockMatrix.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixOp.hpp"  // For matrix operators declaration
#include "SimpleMatrix.hpp"
// #define DEBUG_MESSAGES
// #include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "siconos_debug.h"

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::algebra::operator*(
    const std::shared_ptr<SimpleMatrix> A, const std::shared_ptr<SimpleMatrix> B) {
  auto aux = std::make_shared<SimpleMatrix>((DenseMat)prod(*(*A).dense(), *(*B).dense()));
  return aux;
}

const siconos::algebra::SimpleMatrix siconos::algebra::operator*(const SiconosMatrix &A,
                                                                 double a) {
  // To compute B = a * A

  auto numA = A.num();

  if (numA == UblasType::ZERO)  // if A = 0
  {
    // DenseMat p(zero_matrix(A.size(0),A.size(1)));
    // return p;
    return A;
  } else if (numA == UblasType::IDENTITY) {
    return (DenseMat)(a * *A.identity());
  } else if (numA == UblasType::BLOCK)  // A block
  {
    SimpleMatrix tmp(A);  // ... copy ...
    tmp *= a;
    return tmp;
  } else if (numA == UblasType::DENSE)  // dense)
    return (DenseMat)(a * *A.dense());
  else if (numA == UblasType::TRIANGULAR)
    return (TriangMat)(a * *A.triang());
  else if (numA == UblasType::SYMMETRIC)
    return (SymMat)(a * *A.sym());
  else if (numA == UblasType::SPARSE)
    return (SparseMat)(a * *A.sparse());
  else if (numA == UblasType::BANDED)
    return (BandedMat)(a * *A.banded());
  else {
    THROW_EXCEPTION("invalid type of matrix");
  }
}

siconos::algebra::SimpleMatrix siconos::algebra::operator*(double a, const SiconosMatrix &A) {
  // To compute B = a * A

  auto numA = A.num();

  if (numA == UblasType::ZERO)  // if A = 0
  {
    // DenseMat p(zero_matrix(A.size(0),A.size(1)));
    // return p;
    return A;
  } else if (numA == UblasType::IDENTITY) {
    return (DenseMat)(a * *A.identity());
  } else if (numA == UblasType::BLOCK)  // A block
  {
    SimpleMatrix tmp(A);  // ... copy ...
    tmp *= a;
    return tmp;
  } else if (numA == UblasType::DENSE)  // dense)
    return (DenseMat)(a * *A.dense());
  else if (numA == UblasType::TRIANGULAR)
    return (TriangMat)(a * *A.triang());
  else if (numA == UblasType::SYMMETRIC)
    return (SymMat)(a * *A.sym());
  else if (numA == UblasType::SPARSE)
    return (SparseMat)(a * *A.sparse());
  else if (numA == UblasType::SPARSE_COORDINATE)
    return (SparseCoordinateMat)(a * *A.sparseCoordinate());
  else if (numA == UblasType::BANDED)
    return (BandedMat)(a * *A.banded());
  else {
    THROW_EXCEPTION("invalid type of matrix");
  }
}

const siconos::algebra::SimpleMatrix siconos::algebra::operator/(const SiconosMatrix &A,
                                                                 double a) {
  // To compute B = A/a

  if (a == 0.0) THROW_EXCEPTION("division by zero.");

  auto numA = A.num();

  if (numA == UblasType::ZERO)  // if A = 0
  {
    // DenseMat p(zero_matrix(A.size(0),A.size(1)));
    // return p;
    return A;
  } else if (numA == UblasType::IDENTITY) {
    return (DenseMat)(*A.identity() / a);
  } else if (numA == UblasType::BLOCK)  // A block
  {
    SimpleMatrix tmp(A);  // ... copy ...
    tmp /= a;
    return tmp;
  } else if (numA == UblasType::DENSE)  // dense)
    return (DenseMat)(*A.dense() / a);
  else if (numA == UblasType::TRIANGULAR)
    return (TriangMat)(*A.triang() / a);
  else if (numA == UblasType::SYMMETRIC)
    return (SymMat)(*A.sym() / a);
  else if (numA == UblasType::SPARSE)
    return (SparseMat)(*A.sparse() / a);
  else if (numA == UblasType::BANDED)
    return (BandedMat)(*A.banded() / a);
  else {
    THROW_EXCEPTION("invalid type of matrix");
  }
}

// const SimpleMatrix operator + (const  SimpleMatrix& A, const  SimpleMatrix& B){
//   return (DenseMat)(*A.dense() + *B.dense());
// }
// siconos::algebra::SimpleMatrix siconos::algebra::operator+(const SimpleMatrix &A,
//                                                            const SimpleMatrix &B)
// {
//   return (DenseMat)(*A.dense() + *B.dense());
// }

void siconos::algebra::operator+=(std::shared_ptr<SiconosMatrix> A,
                                  std::shared_ptr<SimpleMatrix> B) {
  *A += *B;
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::algebra::operator+(
    const std::shared_ptr<SimpleMatrix> A, const std::shared_ptr<SimpleMatrix> B) {
  return std::make_shared<SimpleMatrix>(*A + *B);
}

const siconos::algebra::SimpleMatrix siconos::algebra::operator+(const SiconosMatrix &A,
                                                                 const SiconosMatrix &B) {
  // To compute C = A + B

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    THROW_EXCEPTION("inconsistent sizes");

  auto numA = A.num();
  auto numB = B.num();

  // == A or B equal to null ==
  if (numA == UblasType::ZERO)  // A = 0
  {
    if (numB == UblasType::ZERO)  // B = 0
      return SimpleMatrix(A.size(0), A.size(1));
    else
      return SimpleMatrix(B);
  }

  if (numB == UblasType::ZERO) return SimpleMatrix(A);

  // == A and B different from 0 ==

  if (numA == numB &&
      numA != UblasType::BLOCK)  // all matrices are of the same type and NOT block
  {
    if (numA == UblasType::DENSE)
      return (DenseMat)(*A.dense() + *B.dense());
    else if (numA == UblasType::TRIANGULAR)
      return (TriangMat)(*A.triang() + *B.triang());
    else if (numA == UblasType::SYMMETRIC)
      return (SymMat)(*A.sym() + *B.sym());
    else if (numA == UblasType::SPARSE) {
      SparseMat tmp(*A.sparse());
      tmp += *B.sparse();
      return tmp;
      // return (SparseMat)(*A.sparse() + *B.sparse());
    } else if (numA == UblasType::SPARSE_COORDINATE) {
      SparseMat tmp(*A.sparseCoordinate());
      tmp += *B.sparseCoordinate();
      return tmp;
    } else if (numA == UblasType::BANDED) {
      BandedMat tmp(*A.banded());
      tmp += *B.banded();
      return tmp;
    } else
      THROW_EXCEPTION("invalid type of matrix");
  } else if (numA != UblasType::BLOCK && numB != UblasType::BLOCK &&
             numA != numB)  // A and B of different types and none is block
  {
    if (numA == UblasType::DENSE) {
      if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.dense() + *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.dense() + *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.dense() + *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.dense() + *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.dense() + *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.dense() + *B.identity());
      else
        THROW_EXCEPTION("invalid type of matrix");
    } else if (numA == UblasType::TRIANGULAR) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.triang() + *B.dense());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.triang() + *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.triang() + *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.triang() + *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.triang() + *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.triang() + *B.identity());
      else
        THROW_EXCEPTION("invalid type of matrix");
    } else if (numA == UblasType::SYMMETRIC) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.sym() + *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.sym() + *B.triang());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.sym() + *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.sym() + *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.sym() + *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.sym() + *B.identity());
      else
        THROW_EXCEPTION("invalid type of matrix");
    } else if (numA == UblasType::SPARSE) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.sparse() + *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.sparse() + *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.sparse() + *B.sym());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.sparse() + *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.sparse() + *B.identity());
      else
        THROW_EXCEPTION("invalid type of matrix");
    }

    else if (numA == UblasType::BANDED) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.banded() + *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.banded() + *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.banded() + *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.banded() + *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.banded() + *B.sparseCoordinate());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.banded() + *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    }

    else if (numA == UblasType::IDENTITY) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.identity() + *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.identity() + *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.identity() + *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.identity() + *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.identity() + *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.identity() + *B.banded());
    } else
      THROW_EXCEPTION("invalid type of matrix");
  } else if (numB != UblasType::BLOCK)  // B Simple, whatever is A
  {
    SimpleMatrix tmp(B);
    tmp += A;
    return tmp;
  } else  // B Block, A simple or block
  {
    SimpleMatrix tmp(A);
    tmp += B;
    return tmp;
  }
  THROW_EXCEPTION("invalid type of matrix");
}

const siconos::algebra::SimpleMatrix siconos::algebra::operator-(const SiconosMatrix &A,
                                                                 const SiconosMatrix &B) {
  // To compute C = A - B

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    THROW_EXCEPTION("inconsistent sizes");

  auto numA = A.num();
  auto numB = B.num();

  // == B equal to null ==
  if (numB == UblasType::ZERO) return SimpleMatrix(A);

  // == B different from 0 ==

  if (numA == numB &&
      numA != UblasType::BLOCK)  // all matrices are of the same type and NOT block
  {
    if (numA == UblasType::DENSE)
      return (DenseMat)(*A.dense() - *B.dense());
    else if (numA == UblasType::TRIANGULAR)
      return (TriangMat)(*A.triang() - *B.triang());
    else if (numA == UblasType::SYMMETRIC)
      return (SymMat)(*A.sym() - *B.sym());
    else if (numA == UblasType::SPARSE) {
      SparseMat tmp(*A.sparse());
      tmp -= *B.sparse();
      return tmp;
      // return (SparseMat)(*A.sparse() - *B.sparse());
    } else if (numA == UblasType::SPARSE_COORDINATE) {
      SparseCoordinateMat tmp(*A.sparseCoordinate());
      tmp -= *B.sparseCoordinate();
      return tmp;
    } else if (numA == UblasType::BANDED) {
      BandedMat tmp(*A.banded());
      tmp -= *B.banded();
      return tmp;
      // return (BandedMat)(*A.banded() - *B.banded());
    } else {
      THROW_EXCEPTION("invalid type of matrix");
    }
  } else if (numA != UblasType::BLOCK && numB != UblasType::BLOCK &&
             numA != numB)  // A and B of different types and none is block
  {
    if (numA == UblasType::DENSE) {
      if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.dense() - *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.dense() - *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.dense() - *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.dense() - *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.dense() - *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.dense() - *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    } else if (numA == UblasType::TRIANGULAR) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.triang() - *B.dense());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.triang() - *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.triang() - *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.triang() - *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.triang() - *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.triang() - *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    } else if (numA == UblasType::SYMMETRIC) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.sym() - *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.sym() - *B.triang());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.sym() - *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.sym() - *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.sym() - *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.sym() - *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    } else if (numA == UblasType::SPARSE) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.sparse() - *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.sparse() - *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.sparse() - *B.sym());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.sparse() - *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.sparse() - *B.banded());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.sparse() - *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    }

    else if (numA == UblasType::BANDED) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.banded() - *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.banded() - *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.banded() - *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.banded() - *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.banded() - *B.sparseCoordinate());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.banded() - *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    }

    else if (numA == UblasType::ZERO) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.zero_mat() - *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.zero_mat() - *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.zero_mat() - *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.zero_mat() - *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.zero_mat() - *B.sparseCoordinate());
      else if (numB == UblasType::IDENTITY)
        return (DenseMat)(*A.zero_mat() - *B.identity());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    } else if (numA == UblasType::IDENTITY) {
      if (numB == UblasType::DENSE)
        return (DenseMat)(*A.identity() - *B.dense());
      else if (numB == UblasType::TRIANGULAR)
        return (DenseMat)(*A.identity() - *B.triang());
      else if (numB == UblasType::SYMMETRIC)
        return (DenseMat)(*A.identity() - *B.sym());
      else if (numB == UblasType::SPARSE)
        return (DenseMat)(*A.identity() - *B.sparse());
      else if (numB == UblasType::SPARSE_COORDINATE)
        return (DenseMat)(*A.identity() - *B.sparseCoordinate());
      else if (numB == UblasType::BANDED)
        return (DenseMat)(*A.identity() - *B.banded());
      else {
        THROW_EXCEPTION("invalid type of matrix");
      }
    } else {
      THROW_EXCEPTION("invalid type of matrix");
    }
  } else  // A and/or B are/is Block
  {
    SimpleMatrix tmp(A);
    tmp -= B;
    return tmp;
  }
}

//========================
// Matrices comparison
//========================

bool siconos::algebra::operator==(const SiconosMatrix &m, const SiconosMatrix &x) {
  //  if( ! isComparableTo( m, x))
  //    return false;
  // Warning : two block matrices may be "equal" but have blocks of different sizes.
  double norm = (m - x).normInf();
  DEBUG_EXPR((m - x).display());
  DEBUG_PRINTF("norm = %12.8e \n", norm);
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n",
               std::numeric_limits<double>::epsilon());
  DEBUG_EXPR(std::cout << std::boolalpha << (norm <= std::numeric_limits<double>::epsilon())
                       << std::endl;);
  double atol = 1e-14;
  double rtol = std::numeric_limits<double>::epsilon();
  return (norm <= atol + rtol * x.normInf());
}

bool siconos::algebra::operator!=(const SiconosMatrix &m, const SiconosMatrix &x) {
  double norm = (m - x).normInf();
  return (norm > std::numeric_limits<double>::epsilon());
}
