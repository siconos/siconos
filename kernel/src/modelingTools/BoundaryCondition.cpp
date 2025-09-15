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

#include "BoundaryCondition.hpp"



// BoundaryCondition::BoundaryCondition()
// _velocityIndices(nullptr)
// {
// }


BoundaryCondition::BoundaryCondition(SP::UnsignedIntVector newVelocityIndices, SP::SiconosVector newVelocityValues): _velocityIndices(newVelocityIndices),  _prescribedVelocity(newVelocityValues)
{

  if(newVelocityIndices->size() != newVelocityValues->size())
    THROW_EXCEPTION("BoundaryCondition::BoundaryCondition  constructor. velocityIndices and prescribedVelocity must have the same size");
  _prescribedVelocityOld.reset(new SiconosVector(*newVelocityValues));
  _pluginPrescribedVelocity.reset(new PluggedObject());
}

BoundaryCondition::BoundaryCondition(SP::UnsignedIntVector newVelocityIndices): _velocityIndices(newVelocityIndices)
{
  _prescribedVelocityOld.reset(new SiconosVector(newVelocityIndices->size()));
  _pluginPrescribedVelocity.reset(new PluggedObject());
}


BoundaryCondition::~BoundaryCondition()
{
}

void BoundaryCondition::computePrescribedVelocity(double time)
{
  if(_pluginPrescribedVelocity->fPtr)
    ((FPtrPrescribedVelocity)_pluginPrescribedVelocity->fPtr)(time, _velocityIndices->size(), &(*_prescribedVelocity)(0));
}

void BoundaryCondition::display()
{
  std::cout << "=====  BoundaryCondition display ===== " <<std::endl;
  std::cout << "- indices : " <<std::endl;
  if(_velocityIndices)
  {
    for (unsigned int i : *_velocityIndices)
    {
      std::cout << i << " " ;
    }
    std::cout << std::endl;
  }
  else std::cout << "-> nullptr" <<std::endl;
  std::cout << "- velocities : " <<std::endl;
  if(_prescribedVelocity)
    _prescribedVelocity->display();
  else
    std::cout << "-> nullptr" <<std::endl;
  std::cout << "=========================================================== " <<std::endl;
}
