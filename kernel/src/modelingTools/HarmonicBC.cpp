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
#include "BoundaryCondition.hpp"


// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

HarmonicBC::HarmonicBC(SP::UnsignedIntVector newVelocityIndices,
                       double a, double b,
                       double omega, double phi) :
  BoundaryCondition(newVelocityIndices),_a(a),_b(b),_omega(omega),_phi(phi)
{
};

HarmonicBC::HarmonicBC(SP::UnsignedIntVector newVelocityIndices,
                       SP::SiconosVector a, SP::SiconosVector b,
                       SP::SiconosVector omega, SP::SiconosVector phi):
  BoundaryCondition(newVelocityIndices),_aV(a),_bV(b),_omegaV(omega),_phiV(phi)
{
  DEBUG_BEGIN("HarmonicBC::Harmonic((SP::UnsignedIntVector newVelocityIndices,\
               SP::SiconosVector a, SP::SiconosVector b,                \
               SP::SiconosVector omega, SP::SiconosVector phi)\n");

  if(newVelocityIndices->size() != a->size()     ||
      newVelocityIndices->size() != b->size()     ||
      newVelocityIndices->size() != omega->size() ||
      newVelocityIndices->size() != phi->size())
    THROW_EXCEPTION("HarmonicBC::HarmonicBC indices and vectors of data \
           (a,b,omega,phi) must be of the same size ");


  DEBUG_END("HarmonicBC::Harmonic((SP::UnsignedIntVector newVelocityIndices,\
           SP::SiconosVector a, SP::SiconosVector b,\
           SP::SiconosVector omega, SP::SiconosVector phi)\n");


};


HarmonicBC::~HarmonicBC()
{
}
void HarmonicBC::computePrescribedVelocity(double time)
{
  DEBUG_BEGIN("HarmonicBC::computePrescribedVelocity(double time)\n");
  if(!_prescribedVelocity) _prescribedVelocity.reset(new SiconosVector((unsigned int)_velocityIndices->size()));
  if(!_aV)
  {
    for(unsigned int k = 0 ; k < _velocityIndices->size(); k++)
    {
      _prescribedVelocity->setValue(k,_a+_b*cos(_omega*time+_phi));
      DEBUG_PRINTF("_prescribedVelocity[%i] at time  %e = %e \n",k, time,_prescribedVelocity->getValue(k));
    }
  }
  else
  {
    for(unsigned int k = 0 ; k < _velocityIndices->size(); k++)
    {
      _prescribedVelocity->setValue(k,(*_aV)(k)+(*_bV)(k)*cos((*_omegaV)(k)*time+(*_phiV)(k)));
      DEBUG_PRINTF("_prescribedVelocity[%i] at time  %e = %e \n",k, time,_prescribedVelocity->getValue(k));
    }
  }

  DEBUG_END("HarmonicBC::computePrescribedVelocity(double time)\n");
}
