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

#include "FirstOrderLinearTIDS.hpp"
#include "TimeStepping.hpp"
#include "Relay.hpp"
#include "EventsManager.hpp"

#include "LinearSMC.hpp"

#include "SiconosVector.hpp"
#include "ControlSensor.hpp"
#include "ZeroOrderHoldOSI.hpp"
#include "TimeDiscretisation.hpp"
#include "ActuatorFactory.hpp"

LinearSMC::LinearSMC(SP::ControlSensor sensor, unsigned int type):
  CommonSMC(type, sensor)
{
}

LinearSMC::LinearSMC(SP::ControlSensor sensor, SP::SimpleMatrix B, SP::SimpleMatrix D, unsigned int type):
  CommonSMC(type, sensor, B, D)
{
}

LinearSMC::~LinearSMC()
{
}

void LinearSMC::actuate()
{

  if (!_noUeq)
  {
    computeUeq();
    prod(*_B, *_ueq, *(_DS_SMC->b()));
  }

  *(_DS_SMC->x()) = _sensor->y();
  *(_DS_SMC->xMemory()->getSiconosVector(0)) = _sensor->y();

  Type::Siconos dsType = Type::value(*_DS_SMC);
  if (dsType == Type::FirstOrderNonLinearDS)
  {
    _DS_SMC->computef(_simulationSMC->startingTime(), _DS_SMC->x());
    _DS_SMC->swapInMemory();
//    _DS_SMC->computef(_simulationSMC->startingTime());
//    *_DS_SMC->fold() = *_DS_SMC->f();
  }

  _simulationSMC->computeOneStep();
  //  if (_indx > 0)
  {
    _simulationSMC->nextStep();
  }


  // discontinous part
  *_us = *_lambda;
  *_u = *_us;
  *_u += *_ueq;
  _indx++;

}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC, LinearSMC)
