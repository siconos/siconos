/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "BlockMatrix.hpp"
#include "SiconosAlgebraTools.hpp"  // for isComparableTo
#include "SiconosAlgebraTypes.hpp"  // UblasType
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // For isComparableto
#include "SimpleMatrix.hpp"

void siconos::algebra::scal(double a, const SiconosMatrix &A, SiconosMatrix &B, bool init) {
  // To compute B = a * A (init = true) or B += a*A (init = false).
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");
  if (!B.isBlock()) B.resetFactorizationFlags();

  if (&A == &B) {
    if (init)
      B *= a;
    else
      B *= (1.0 + a);
  } else {
    auto numA = A.num();
    auto numB = B.num();

    if (numB == UblasType::ZERO || numB == UblasType::IDENTITY)  // B = 0 or identity.
      THROW_EXCEPTION("forbidden for B being a zero or identity matrix.");

    if (numA == UblasType::ZERO) {
      if (init) B.zero();  // else nothing
    } else if (numA == UblasType::IDENTITY) {
      if (init) {
        B.eye();
        B *= a;
      } else {
        // Assuming B is square ...
        for (unsigned int i = 0; i < B.size(0); ++i) B(i, i) += a;
      }
    } else {
      if (numA == numB)  // if A and B are of the same type ...
      {
        switch (numA) {
          case UblasType::BLOCK:  // A and B are block
            if (isComparableTo(A, B)) {
              const BlockMatrix &Aref = static_cast<const BlockMatrix &>(A);
              BlockMatrix &Bref = static_cast<BlockMatrix &>(B);
              auto itA1 = Aref._mat->begin1();
              for (auto itB1 = Bref._mat->begin1(); itB1 != Bref._mat->end1(); ++itB1) {
                auto itA2 = itA1.begin();
                for (auto itB2 = itB1.begin(); itB2 != itB1.end(); ++itB2) {
                  scal(a, **itA2++, **itB2, init);
                }
                itA1++;
              }
            } else  // if A and B are not "block-consistent"
            {
              if (init) {
                for (unsigned int i = 0; i < A.size(0); ++i)
                  for (unsigned int j = 0; j < A.size(1); ++j) B(i, j) = a * A(i, j);
              } else {
                for (unsigned int i = 0; i < A.size(0); ++i)
                  for (unsigned int j = 0; j < A.size(1); ++j) B(i, j) += a * A(i, j);
              }
            }
            break;

          case UblasType::DENSE:  // if both are dense
            if (init)
              noalias(*B.dense()) = a * *A.dense();
            else
              noalias(*B.dense()) += a * *A.dense();
            break;
          case UblasType::TRIANGULAR:
            if (init)
              noalias(*B.triang()) = a * *A.triang();
            else
              noalias(*B.triang()) += a * *A.triang();
            break;
          case UblasType::SYMMETRIC:
            if (init)
              noalias(*B.sym()) = a * *A.sym();
            else
              noalias(*B.sym()) += a * *A.sym();
            break;
          case UblasType::SPARSE:
            if (init)
              noalias(*B.sparse()) = a * *A.sparse();
            else
              noalias(*B.sparse()) += a * *A.sparse();
            break;
          case UblasType::BANDED:
            if (init)
              noalias(*B.banded()) = a * *A.banded();
            else
              noalias(*B.banded()) += a * *A.banded();
            break;
          default:
            THROW_EXCEPTION("Not implemented for A/B type.");
        }
      } else  // if A and B are of different types.
      {
        if (numA == UblasType::BLOCK || numB == UblasType::BLOCK)  // if A or B is block
        {
          if (init) {
            B = A;
            B *= a;
          } else {
            SimpleMatrix tmp(A);
            tmp *= a;
            B += tmp;  // bof bof ...
          }
        } else {
          if (numB != UblasType::DENSE)
            THROW_EXCEPTION("Inconsistent types between A and B (must be dense?)");

          if (init) {
            switch (numA) {
              case UblasType::DENSE:
                noalias(*B.dense()) = a * *A.dense();
                break;
              case UblasType::TRIANGULAR:
                noalias(*B.dense()) = a * *A.triang();
                break;
              case UblasType::SYMMETRIC:
                noalias(*B.dense()) = a * *A.sym();
                break;
              case UblasType::SPARSE:
                noalias(*B.dense()) = a * *A.sparse();
                break;
              case UblasType::BANDED:
                noalias(*B.dense()) = a * *A.banded();
                break;
              default:
                THROW_EXCEPTION("Not implemented for A type.");
            }
          } else

          {
            switch (numA) {
              case UblasType::DENSE:
                noalias(*B.dense()) += a * *A.dense();
                break;
              case UblasType::TRIANGULAR:
                noalias(*B.dense()) += a * *A.triang();
                break;
              case UblasType::SYMMETRIC:
                noalias(*B.dense()) += a * *A.sym();
                break;
              case UblasType::SPARSE:
                noalias(*B.dense()) += a * *A.sparse();
                break;
              case UblasType::BANDED:
                noalias(*B.dense()) += a * *A.banded();
                break;
              default:
                THROW_EXCEPTION("Not implemented for A type.");
            }
          }
        }
      }
    }
  }
}
