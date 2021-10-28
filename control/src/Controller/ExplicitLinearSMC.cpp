/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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


#include "ActuatorFactory.hpp"
#include "SiconosVector.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "ControlSensor.hpp"
#include "SiconosAlgebraProd.hpp"
#include "ExplicitLinearSMC.hpp"


ExplicitLinearSMC::ExplicitLinearSMC(SP::ControlSensor sensor): CommonSMC(EXPLICIT_LINEAR_SMC, sensor)
{
}

ExplicitLinearSMC::ExplicitLinearSMC(SP::ControlSensor sensor, SP::SimpleMatrix B): CommonSMC(EXPLICIT_LINEAR_SMC, sensor, B)
{
}


ExplicitLinearSMC::~ExplicitLinearSMC()
{
  _sigma.reset();
}

void ExplicitLinearSMC::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation &s)
{
  CommonSMC::initialize(nsds,s);

  _sigma.reset(new SiconosVector(_u->size()));
}

void ExplicitLinearSMC::actuate()
{
  if(!_noUeq)
  {
    computeUeq();
  }

  prod(*_Csurface, _sensor->y(), *_sigma);

  unsigned int sDim = _u->size();

  if(_D)  // we are using a saturation
  {
    for(unsigned int i = 0; i < sDim; i++)
    {
      if((*_sigma)(i) > (*_D)(i, i))
        (*_us)(i) = -_alpha;
      else if((*_sigma)(i) < -(*_D)(i, i))
        (*_us)(i) = _alpha;
      else
      {
        if((*_D)(i, i) != 0)
          (*_us)(i) = -(*_sigma)(i) / (*_D)(i, i);
        else
          (*_us)(i) = 0;
      }
    }
  }
  else
  {
    for(unsigned int i = 0; i < sDim; i++)
    {
      if((*_sigma)(i) > 0)
        (*_us)(i) = -_alpha;
      else if((*_sigma)(i) < 0)
        (*_us)(i) = _alpha;
      else
        (*_us)(i) = 0;
    }
  }

  *_lambda = *_us;
  *_u = *_us;
  *_u += *_ueq;
}

AUTO_REGISTER_ACTUATOR(EXPLICIT_LINEAR_SMC, ExplicitLinearSMC)
