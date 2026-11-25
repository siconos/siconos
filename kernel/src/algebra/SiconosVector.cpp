/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
