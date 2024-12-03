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
#include "SiconosVectorOp.hpp"  // For setBlock decl. ...

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


void siconos::algebra::setBlock(const SiconosVector& vIn, std::shared_ptr<SiconosVector> vOut,
                                unsigned int sizeB, unsigned int startIn,
                                unsigned int startOut) {
  unsigned int endOut = startOut + sizeB;
  assert(vOut->size() >= endOut && "The output vector is too small");

  vOut->segment(startOut, sizeB).noalias() = vIn.segment(startIn, sizeB);
}
