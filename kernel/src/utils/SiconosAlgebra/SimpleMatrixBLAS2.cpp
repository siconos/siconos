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

#include "SiconosConfig.h"
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace siconosBindings = boost::numeric::bindings;

// for ublas::axpy_prod, ...
#include <boost/numeric/ublas/operation.hpp>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"
#include "SiconosAlgebra.hpp"

using namespace Siconos;

// void axpy_prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
// {
//   // To compute y = A * x ( init = true) or y += A * x (init = false) using ublas::axpy_prod
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   if(A.size(1) != x.size())
//     SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

//   if(A.size(0) != y.size())
//     SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

//   unsigned int numA = A.num();
//   unsigned int numX = x.num();
//   unsigned int numY = y.num();

//   if(numA == 0)  // If A is Block
//     SiconosMatrixException::selfThrow("axpy_prod(A,x,y) error: not yet implemented for block matrices.");

//   if(numA == ZERO)  // A = 0
//   {
//     if(init) y.zero();  // else nothing ...
//   }

//   else if(numA == IDENTITY)  // A = identity
//   {
//     if(!init) y += x;
//     else
//     {
//       if(&x != &y)
//         y = x ; // if x and y do not share memory (ie are different objects)
//     }
//     // else nothing
//   }

//   else // A is not 0 or identity
//   {
//     {
//       {
//         if(&x != &y)  // if no common memory between x and y.
//         {
//           if(numX == DENSE)
//           {
//             if(numY != DENSE)
//               SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

//             if(numA == DENSE)
//               ublas::axpy_prod(*A.dense(), *x.dense(), *y.dense(), init);
//             else if(numA == TRIANGULAR)
//               ublas::axpy_prod(*A.triang(), *x.dense(), *y.dense(), init);
//             else if(numA == SYMMETRIC)
//               ublas::axpy_prod(*A.sym(), *x.dense(), *y.dense(), init);
//             else if(numA == SPARSE)
//               ublas::axpy_prod(*A.sparse(), *x.dense(), *y.dense(), init);
//             else //if(numA==BANDED)
//               ublas::axpy_prod(*A.banded(), *x.dense(), *y.dense(), init);
//           }
//           else //if(numX == SPARSE)
//           {
//             if(numY != DENSE && numA != SPARSE)
//               SiconosMatrixException::selfThrow("axpy_prod(A,x,y) error: y (output) must be a dense vector.");

//             if(numA == DENSE)
//               ublas::axpy_prod(*A.dense(), *x.sparse(), *y.dense(), init);
//             else if(numA == TRIANGULAR)
//               ublas::axpy_prod(*A.triang(), *x.sparse(), *y.dense(), init);
//             else if(numA == SYMMETRIC)
//               ublas::axpy_prod(*A.sym(), *x.sparse(), *y.dense(), init);
//             else if(numA == SPARSE)
//             {
//               if(numY == DENSE)
//                 ublas::axpy_prod(*A.sparse(), *x.sparse(), *y.dense(), init);
//               else
//                 ublas::axpy_prod(*A.sparse(), *x.sparse(), *y.sparse(), init);
//             }
//             else //if(numA==BANDED)
//               ublas::axpy_prod(*A.banded(), *x.sparse(), *y.dense(), init);
//           }
//         }
//         else // if x and y are the same object => alias
//         {
//           if(numX == DENSE)
//           {
//             if(numA == DENSE)
//               ublas::axpy_prod(*A.dense(), *x.dense(), *x.dense(), init);
//             else if(numA == TRIANGULAR)
//               ublas::axpy_prod(*A.triang(), *x.dense(), *x.dense(), init);
//             else if(numA == SYMMETRIC)
//               ublas::axpy_prod(*A.sym(), *x.dense(), *x.dense(), init);
//             else if(numA == SPARSE)
//               ublas::axpy_prod(*A.sparse(), *x.dense(), *x.dense(), init);
//             else //if(numA==BANDED)
//               ublas::axpy_prod(*A.banded(), *x.dense(), *x.dense(), init);
//           }
//           else //if(numX == SPARSE)
//           {
//             if(numA == DENSE)
//               ublas::axpy_prod(*A.dense(), *x.sparse(), *x.sparse(), init);
//             else if(numA == TRIANGULAR)
//               ublas::axpy_prod(*A.triang(), *x.sparse(), *x.sparse(), init);
//             else if(numA == SYMMETRIC)
//               ublas::axpy_prod(*A.sym(), *x.sparse(), *x.sparse(), init);
//             else if(numA == SPARSE)
//               ublas::axpy_prod(*A.sparse(), *x.sparse(), *x.sparse(), init);
//             else //if(numA==BANDED)
//               ublas::axpy_prod(*A.banded(), *x.sparse(), *x.sparse(), init);
//           }
//         }
//       }
//     }
//   }
// }

// void gemvtranspose(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
// {
//   if(A.isBlock())
//     SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   unsigned int numA = A.num();
//   unsigned int numX = x.num();
//   unsigned int numY = y.num();
//   if(numA != DENSE || numX != DENSE || numY != DENSE)
//     SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

//   siconosBindings::blas::gemv(a, siconosBindings::trans(*A.dense()), *x.dense(), b, *y.dense());
// }

// void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
// {
//   if(A.isBlock())
//     SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
//   assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

//   unsigned int numA = A.num();
//   unsigned int numX = x.num();
//   unsigned int numY = y.num();
//   if(numA != DENSE || numX != DENSE || numY != DENSE)
//     SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

//   siconosBindings::blas::gemv(a, *A.dense(), *x.dense(), b, *y.dense());
// }
