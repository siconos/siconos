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

#include "BoundaryCondition.hpp"

#include <iostream>

#include "PluggedObject.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"

siconos::modeling::BoundaryCondition::BoundaryCondition(
    const Indices& newVelocityIndices)
    : _velocityIndices{newVelocityIndices},
      _prescribedVelocity{std::make_shared<siconos::algebra::SiconosVector>(
          _velocityIndices.size(), 0.)},
      _prescribedVelocityOld{std::make_shared<siconos::algebra::SiconosVector>(
          _velocityIndices.size(), 0.)} {}

siconos::modeling::BoundaryCondition::BoundaryCondition(
    Indices&& newVelocityIndices)
    : _velocityIndices{std::move(newVelocityIndices)},
      _prescribedVelocity{std::make_shared<siconos::algebra::SiconosVector>(
          _velocityIndices.size(), 0.)},
      _prescribedVelocityOld{std::make_shared<siconos::algebra::SiconosVector>(
          _velocityIndices.size(), 0.)} {}

siconos::modeling::BoundaryCondition::BoundaryCondition(
    Indices&& newVelocityIndices,
    std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues)
    : _velocityIndices{std::move(newVelocityIndices)},
      _prescribedVelocity(newVelocityValues),
      _prescribedVelocityOld{std::make_shared<siconos::algebra::SiconosVector>(
          *newVelocityValues)} {
  if (newVelocityIndices.size() != newVelocityValues->size())
    THROW_EXCEPTION(
        "siconos::modeling::BoundaryCondition::BoundaryCondition  "
        "constructor. "
        "velocityIndices and prescribedVelocity must have the same size");
}

siconos::modeling::BoundaryCondition::BoundaryCondition(
    const Indices& newVelocityIndices,
    std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues)
    : _velocityIndices{newVelocityIndices},
      _prescribedVelocity(newVelocityValues),
      _prescribedVelocityOld{std::make_shared<siconos::algebra::SiconosVector>(
          *newVelocityValues)} {
  if (newVelocityIndices.size() != newVelocityValues->size())
    THROW_EXCEPTION(
        "siconos::modeling::BoundaryCondition::BoundaryCondition  "
        "constructor. "
        "velocityIndices and prescribedVelocity must have the same size");
}

void siconos::modeling::BoundaryCondition::setComputePrescribedVelocityFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginPrescribedVelocity =
      std::make_shared<siconos::plugins::PluggedObject>();

  _pluginPrescribedVelocity->setComputeFunction(pluginPath, functionName);

  if (!_prescribedVelocity)
    _prescribedVelocity = std::make_shared<siconos::algebra::SiconosVector>(
        _velocityIndices.size());
}

void siconos::modeling::BoundaryCondition::computePrescribedVelocity(
    double time) {
  if (_pluginPrescribedVelocity && _pluginPrescribedVelocity->fPtr)
    ((FPtrPrescribedVelocity)_pluginPrescribedVelocity->fPtr)(
        time, _velocityIndices.size(), &(*_prescribedVelocity)(0));
}

void siconos::modeling::BoundaryCondition::display() {
  std::cout << "=====  BoundaryCondition display ===== " << std::endl;
  std::cout << "- indices : " << std::endl;
  for (auto i : _velocityIndices) {
    std::cout << i << " ";
  }
  std::cout << std::endl;

  std::cout << "- velocities : " << std::endl;
  if (_prescribedVelocity)
    _prescribedVelocity->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "=========================================================== "
            << std::endl;
}

void siconos::modeling::BoundaryCondition::appendIndex(unsigned int ind)
{
  if (find(_velocityIndices.begin(), _velocityIndices.end(), ind) == _velocityIndices.end())
    {
      _velocityIndices.push_back(ind);
      _prescribedVelocity->resize(size(), 0.);
      _prescribedVelocityOld->resize(size(), 0.);
    }
}
