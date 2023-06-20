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

#include "HarmonicBC.hpp"

#include <math.h>  // for cos

#include "BoundaryCondition.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

siconos::modeling::HarmonicBC::HarmonicBC(
    Indices&& newVelocityIndices, std::shared_ptr<siconos::algebra::SiconosVector> a,
    std::shared_ptr<siconos::algebra::SiconosVector> b,
    std::shared_ptr<siconos::algebra::SiconosVector> omega,
    std::shared_ptr<siconos::algebra::SiconosVector> phi)
    : BoundaryCondition(std::move(newVelocityIndices)),
      _aV(a),
      _bV(b),
      _omegaV(omega),
      _phiV(phi) {
  DEBUG_BEGIN(
      "HarmonicBC::Harmonic((std::shared_ptr<std::vector<unsigned int>> newVelocityIndices,\
               std::shared_ptr<siconos::algebra::SiconosVector> a, std::shared_ptr<siconos::algebra::SiconosVector> b,                \
               std::shared_ptr<siconos::algebra::SiconosVector> omega, std::shared_ptr<siconos::algebra::SiconosVector> phi)\n");

  auto newsize = _velocityIndices.size();
  if (newsize != a->size() || newsize != b->size() || newsize != omega->size() ||
      newsize != phi->size())
    THROW_EXCEPTION(
        "HarmonicBC::HarmonicBC indices and vectors of data \
           (a,b,omega,phi) must be of the same size ");

  DEBUG_END(
      "HarmonicBC::Harmonic((std::shared_ptr<std::vector<unsigned int>> newVelocityIndices,\
           std::shared_ptr<siconos::algebra::SiconosVector> a, std::shared_ptr<siconos::algebra::SiconosVector> b,\
           std::shared_ptr<siconos::algebra::SiconosVector> omega, std::shared_ptr<siconos::algebra::SiconosVector> phi)\n");
};

void siconos::modeling::HarmonicBC::computePrescribedVelocity(double time) {
  DEBUG_BEGIN("HarmonicBC::computePrescribedVelocity(double time)\n");
  if (!_prescribedVelocity)
    _prescribedVelocity =
        std::make_shared<siconos::algebra::SiconosVector>(_velocityIndices.size());
  if (!_aV) {
    for (unsigned int k = 0; k < _velocityIndices.size(); k++) {
      _prescribedVelocity->setValue(k, _a + _b * cos(_omega * time + _phi));
      DEBUG_PRINTF("_prescribedVelocity[%i] at time  %e = %e \n", k, time,
                   _prescribedVelocity->getValue(k));
    }
  } else {
    for (unsigned int k = 0; k < _velocityIndices.size(); k++) {
      _prescribedVelocity->setValue(
          k, (*_aV)(k) + (*_bV)(k)*cos((*_omegaV)(k)*time + (*_phiV)(k)));
      DEBUG_PRINTF("_prescribedVelocity[%i] at time  %e = %e \n", k, time,
                   _prescribedVelocity->getValue(k));
    }
  }

  DEBUG_END("HarmonicBC::computePrescribedVelocity(double time)\n");
}
