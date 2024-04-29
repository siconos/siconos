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

#include "SiconosException.hpp"
#include "SiconosMatrixOp.hpp"  // for matrix op. declarations
#include "SiconosMatrix.hpp"

void siconos::algebra::add(const SiconosMatrix &A, const SiconosMatrix &B,
                           SiconosMatrix &C) {
  // To compute C = A + B in an "optimized" way (in comparison with operator +)

  // === if C is zero or identity => read-only ===
  // if (numC == UblasType::ZERO || numC == UblasType::IDENTITY)
  //   THROW_EXCEPTION(
  //       "Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C "
  //       "(read-only: zero or identity).");

  // === common memory between A, B, C ===
  if (&A == &C)  // A and C have common memory
  {
    C += B;
  } else if (&B == &C)  // B and C have common memory
  {
    C += A;
  } else  // No common memory between C and A or B.
  {
    C = A + B;
  }
}

void siconos::algebra::sub(const SiconosMatrix &A, const SiconosMatrix &B,
                           SiconosMatrix &C) {
  // To compute C = A - B in an "optimized" way (in comparison with operator +)

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    THROW_EXCEPTION("Matrix addition: inconsistent sizes");
  if ((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    THROW_EXCEPTION("Matrix addition: inconsistent sizes");

  // === if C is zero or identity => read-only ===
  // if (numC == UblasType::ZERO || numC == UblasType::IDENTITY)
  //   THROW_EXCEPTION(
  //       "Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C "
  //       "(read-only: zero or identity).");

  // === common memory between A, B, C ===
  if (&A == &C)  // A and C have common memory
  {
    C -= B;
  } else if (&B == &C)  // B and C have common memory
  {
    C *= -1.0;
    C += A;
  } else  // No common memory between C and A or B.
  {
    C = A - B;
  }
}
