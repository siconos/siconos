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

void siconos::algebra::concatenateVectors(
    siconos::algebra::SiconosVector& target,
    const siconos::algebra::SiconosVector& a,
    const siconos::algebra::SiconosVector& b) {
  target.resize(a.size() + b.size());
  target.head(a.size()) = a;
  target.tail(b.size()) = b;
}

void siconos::algebra::scal(double a, const SiconosVector& x, SiconosVector& y,
                            bool init) {
  // To compute y = a *x (init = true) or y += a*x (init = false)

  if (&x == &y) {
    if (init)
      y *= a;
    else {
      y *= (1.0 + a);
    }
  } else {
    unsigned int sizeX = x.size();
    unsigned int sizeY = y.size();

    if (sizeX != sizeY) THROW_EXCEPTION("sizes are not consistent.");

    if (init)
      y.noalias() = a * x;  // WARNING : add noalias
    else
      y.noalias() += a * x;  // WARNING : add noalias
  }
}

void siconos::algebra::getMin(const SiconosVector& V, double& minvalue,
                              unsigned int& idmin) {
  minvalue = V(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it) {
    if (V(it) < minvalue) {
      minvalue = V(it);
      idmin = it;
    }
  }
}

void siconos::algebra::getMax(const SiconosVector& V, double& maxvalue,
                              unsigned int& idmax) {
  maxvalue = V(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it) {
    if (V(it) > maxvalue) {
      maxvalue = V(it);
      idmax = it;
    };
  };
}

void siconos::algebra::sub(const SiconosVector& x, const SiconosVector& y,
                           SiconosVector& z) {
  // Computes z = x - y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size())
    THROW_EXCEPTION("inconsistent sizes");

  if (&z == &x)  // x and z are the same object.
  {
    z -= y;
  } else if (&z == &y)  // y and z are the same object
  {
    z = x - y;
  } else  // No common memory between x or y and z
  {
    z.noalias() = x - y;
  }
}

void siconos::algebra::setBlock(const SiconosVector& vIn,
                                std::shared_ptr<SiconosVector> vOut,
                                unsigned int sizeB, unsigned int startIn,
                                unsigned int startOut) {
  unsigned int endOut = startOut + sizeB;
  assert(vOut->size() >= endOut && "The output vector is too small");

  vOut->segment(startOut, sizeB).noalias() = vIn.segment(startIn, sizeB);
}

void siconos::algebra::cross_product(const SiconosVector& V1,
                                     const SiconosVector& V2,
                                     SiconosVector& VOUT) {
  if (V1.size() != 3 || V2.size() != 3 || VOUT.size() != 3)
    THROW_EXCEPTION("allowed only with dim 3.");

  double aux =
      V1.getValue(1) * V2.getValue(2) - V1.getValue(2) * V2.getValue(1);
  VOUT.setValue(0, aux);

  aux = V1.getValue(2) * V2.getValue(0) - V1.getValue(0) * V2.getValue(2);
  VOUT.setValue(1, aux);

  aux = V1.getValue(0) * V2.getValue(1) - V1.getValue(1) * V2.getValue(0);
  VOUT.setValue(2, aux);
}

void siconos::algebra::abs_wise(const SiconosVector& V, SiconosVector& Vabs) {
  for (unsigned int it = 0; it < V.size(); ++it) {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

double siconos::algebra::inner_prod(const SiconosVector& x,
                                    const SiconosVector& m) {
  if (x.size() != m.size()) THROW_EXCEPTION("inconsistent sizes");
  return x.dot(m);
}

void siconos::algebra::subscal(double a, const SiconosVector& x,
                               SiconosVector& y,
                               const std::vector<std::size_t>& coord,
                               bool init) {
  // To compute sub_y = a *sub_x (init = true) or sub_y += a*sub_x (init =
  // false) Coord  = [r0x r1x r0y r1y]; subX is the sub-vector of x, for row
  // numbers between r0x and r1x-1. The same for y with riy.

  // Check dimensions
  auto dimX = coord[1] - coord[0];
  auto dimY = coord[3] - coord[2];
  if (dimY != dimX)
    THROW_EXCEPTION("inconsistent sizes between (sub)x and (sub)y.");
  if (dimY > y.size() || dimX > x.size())
    THROW_EXCEPTION("input index too large.");

  if (&x == &y)  // if x and y are the same object
  {
    if (coord[0] == coord[2]) {
      if (init)
        y.segment(coord[2], dimY) *= a;
      else
        y.segment(coord[2], dimY) *= (1.0 + a);
    } else {
      if (init)
        y.segment(coord[2], dimY) = a * x.segment(coord[0], dimX);
      else
        y.segment(coord[2], dimY) += a * x.segment(coord[0], dimX);
    }
  } else {
    if (init)
      y.segment(coord[2], dimY).noalias() = a * x.segment(coord[0], dimX);
    else
      y.segment(coord[2], dimY).noalias() += a * x.segment(coord[0], dimX);
  }
}
