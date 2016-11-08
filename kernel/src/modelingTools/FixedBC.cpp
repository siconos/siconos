/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "FixedBC.hpp"
#include "BoundaryCondition.hpp"


#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"

FixedBC::FixedBC(SP::UnsignedIntVector newVelocityIndices) :
  BoundaryCondition (newVelocityIndices) 
{
};

FixedBC::~FixedBC()
{
}
void FixedBC::computePrescribedVelocity(double time)
{
  DEBUG_BEGIN("FixedBC::computePrescribedVelocity(double time)\n");
  if (!_prescribedVelocity) _prescribedVelocity.reset(new SiconosVector((unsigned int)_velocityIndices->size()));
  for (unsigned int k = 0 ; k < _velocityIndices->size(); k++)
  {
    _prescribedVelocity->setValue(k,0.0);
    DEBUG_PRINTF("_prescribedVelocity[%i] at time  %e = %e \n",k, time,_prescribedVelocity->getValue(k));
  }
  DEBUG_END("FixedBC::computePrescribedVelocity(double time)\n");
}
