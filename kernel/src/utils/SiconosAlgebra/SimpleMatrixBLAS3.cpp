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

#include <assert.h>
// for ublas::axpy_prod, ...
// #include <boost/numeric/ublas/operation.hpp>
// #include <boost/numeric/ublas/operation_sparse.hpp>
// require for matrix stuff like value_type
// #include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // for matrix operators declaration
#include "SiconosMatrix.hpp"

// namespace ublas = boost::numeric::ublas;
// namespace bindings_blas = boost::numeric::bindings::blas;

//======================
// Product of matrices
//======================
// Note FP: this function is never used. We keep it for the record. Remove it
// later ?

// const SiconosMatrix prod(const SiconosMatrix &A, const SiconosMatrix& B)
// {
//   // To compute C = A * B
//   assert(!(B.isPLUFactorized()) && "B is PLUFactorized in prod !!");
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   if((A.size(1) != B.size(0)))
//     THROW_EXCEPTION("Matrix function C=prod(A,B): inconsistent sizes");

//   auto numA = A.num();
//   auto numB = B.num();

//   // == TODO: implement block product ==
//   if(numA == 0 || numB == 0)
//     THROW_EXCEPTION("Matrix product ( C=prod(A,B) ): not yet implemented for
//     BlockMatrix objects.");

//   if(numA == UblasType::IDENTITY || numB == UblasType::ZERO)  // A = identity
//   or B = 0
//     return SiconosMatrix(B);

//   else if(numB == UblasType::IDENTITY || numA == UblasType::ZERO)  // B =
//   identity or A = 0
//     return SiconosMatrix(A);

//   else // neither A or B is equal to identity or zero.
//   {
//     if(numB == UblasType::DENSE)
//     {
//       if(numA == UblasType::DENSE)
//       {
//         DenseMat p(A.size(0), B.size(1));
//         bindings_blas::blas::gemm(1.0, *A.dense(), *B.dense(), 1.0, p);
//         //      return (DenseMat)(prod(*A.dense(),*B.dense()));
//         return p;
//       }
//       else if(numA == UblasType::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.dense()));
//       else if(numA == UblasType::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.dense()));
//       else if(numA == UblasType::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.dense()));
//       else if(numA == UblasType::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.dense()));
//       else// if(numA==UblasType::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.dense()));
//     }
//     else if(numB == UblasType::TRIANGULAR)
//     {
//       if(numA == UblasType::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.triang()));
//       else if(numA == UblasType::TRIANGULAR)
//         return (TriangMat)(prod(*A.triang(), *B.triang()));
//       else if(numA == UblasType::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.triang()));
//       else if(numA == UblasType::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.triang()));
//       else if(numA == UblasType::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.triang()));
//       else //if(numA==UblasType::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.triang()));
//     }
//     else if(numB == UblasType::SYMMETRIC)
//     {
//       if(numA == UblasType::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.sym()));
//       else if(numA == UblasType::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.sym()));
//       else if(numA == UblasType::SYMMETRIC)
//         return (SymMat)(prod(*A.sym(), *B.sym()));
//       else if(numA == UblasType::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.sym()));
//       else if(numA == UblasType::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.sym()));
//       else // if (numA == UblasType::BANDED)
//         return (DenseMat)(prod(*A.banded(), *B.sym()));
//     }
//     else if(numB == UblasType::SPARSE)
//     {
//       if(numA == UblasType::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.sparse()));
//       else if(numA == UblasType::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.sparse()));
//       else if(numA == UblasType::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.sparse()));
//       else if(numA == UblasType::SPARSE)
//         return (SparseMat)(prod(*A.sparse(), *B.sparse()));
//       else if(numA == UblasType::SPARSE_COORDINATE)
//         return (SparseMat)(prod(*A.sparseCoordinate(), *B.sparse()));
//       else //if(numA==UblasType::BANDED){
//         return (DenseMat)(prod(*A.banded(), *B.sparse()));
//     }
//     else if(numB == UblasType::SPARSE_COORDINATE)
//     {
//       if(numA == UblasType::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.sparseCoordinate()));
//       else if(numA == UblasType::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.sparseCoordinate()));
//       else if(numA == UblasType::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.sparseCoordinate()));
//       else if(numA == UblasType::SPARSE)
//         return (SparseMat)(prod(*A.sparse(), *B.sparseCoordinate()));
//       else if(numA == UblasType::SPARSE_COORDINATE)
//         return (SparseMat)(prod(*A.sparseCoordinate(),
//         *B.sparseCoordinate()));
//       else //if(numA==UblasType::BANDED){
//         return (DenseMat)(prod(*A.banded(), *B.sparseCoordinate()));
//     }
//     else //if(numB==UblasType::BANDED)
//     {
//       if(numA == UblasType::DENSE)
//         return (DenseMat)(prod(*A.dense(), *B.banded()));
//       else if(numA == UblasType::TRIANGULAR)
//         return (DenseMat)(prod(*A.triang(), *B.banded()));
//       else if(numA == UblasType::SYMMETRIC)
//         return (DenseMat)(prod(*A.sym(), *B.banded()));
//       else if(numA == UblasType::SPARSE)
//         return (DenseMat)(prod(*A.sparse(), *B.banded()));
//       else if(numA == UblasType::SPARSE_COORDINATE)
//         return (DenseMat)(prod(*A.sparseCoordinate(), *B.banded()));
//       else //if(numA==UblasType::BANDED)
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
// void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C,
// int indexACol, bool init){
//   // To compute C[indexAcol::] = A * B

//   auto numA = A.num();
//   auto numB = B.num();
//   auto numC = C.num();
//   if (numA == 0 || numB == 0 || numC == 0)
//     THROW_EXCEPTION("Matrix function prod(A,B,C,index): inconsistent sizes");
//   // === if C is zero or identity => read-only ===
//   if (numC == UblasType::Zero || numC == UblasType::IDENTITY)
//     THROW_EXCEPTION("Matrix product ( prod(A,B,C,index) ): wrong type for
//     resulting matrix C (read-only: zero or identity).");

//   if (numA == UblasType::IDENTITY || numC == UblasType::Zero) // A = identity
//   or 0
//     THROW_EXCEPTION("Matrix function prod(A,B,C,index): numA ==
//     UblasType::IDENTITY || numC
//     == UblasType::Zero not yet implemented");

//   int rawB = B.size(0);
//   int colB = B.size(1);

// }

void siconos::algebra::axpy_prod(const SiconosMatrix &A, const SiconosMatrix &B,
                                 SiconosMatrix &C, bool init) {
  // To compute C = A * B (init = true) or C += A * B (init = false) using ublas
  // axpy_prod. High speedup for sparse matrices. Warning FP:
  // ublas::axpy_prod(A, B, C, init) with init = True is equivalent to C = A*B
  // with C.clear BEFORE product. So C==A or B must be forbidden. See
  // http://www.boost.org/doc/libs/1_63_0/libs/numeric/ublas/doc/products.html
  //

  if (init == true) {
    C.setZero();
  }
  C += A * B;
}
