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

#include "SiconosVector.hpp"

#include "SiconosException.hpp"

void siconos::algebra::concatenateVectors(siconos::algebra::SiconosVector& target,
                                          const siconos::algebra::SiconosVector& a,
                                          const siconos::algebra::SiconosVector& b) {
  target.resize(a.size() + b.size());
  target.head(a.size()) = a;
  target.tail(b.size()) = b;
}

std::shared_ptr<siconos::algebra::SiconosVector> siconos::algebra::concatenateVectors(
    const siconos::algebra::SiconosVector& a, const siconos::algebra::SiconosVector& b) {
  std::shared_ptr<siconos::algebra::SiconosVector> tmp{nullptr};
  tmp = std::make_shared<siconos::algebra::SiconosVector>(a.size() + b.size());
  tmp->head(a.size()) = a;
  tmp->tail(b.size()) = b;
  return tmp;
}

bool siconos::algebra::orthoBaseFromVector(SiconosVector3& A, SiconosVector3& A1,
                                           SiconosVector3& A2) {
  double normA = A.norm();
  if (normA == 0.0) {
    // If A is null, we assign Nan to outputs  and return an error code
    A1 = SiconosVector3(std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN());
    A2 = SiconosVector3(std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN());
    return false;
  }

  // Normalize A
  A.normalize();

  double sign = std::copysign(1.0, A.z());
  const double a = -1.0 / (sign + A.z());
  const double b = A.x() * A.y() * a;

  // Build the orthonormal basis using a and b
  A1 << 1.0 + sign * A.x() * A.x() * a, sign * b, -sign * A.x();
  A2 << b, sign + A.y() * A.y() * a, -A.y();

  // and normalize them
  // A1.normalize();
  // A2.normalize();

  // Check norms ...
  assert(std::fabs(A1.norm() - 1.0) < 1e-14);
  assert(std::fabs(A.dot(A1)) < 1e-14);
  assert(std::fabs(A2.norm() - 1.0) < 1e-14);
  assert(std::fabs(A.dot(A2)) < 1e-14);
  assert(std::fabs(A1.dot(A2)) < 1e-14);

  return true;
}
