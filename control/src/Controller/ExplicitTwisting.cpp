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

#include <iostream>

#include <FirstOrderLinearTIDS.hpp>
#include <TimeStepping.hpp>
#include <SiconosVector.hpp>
#include <SimpleMatrix.hpp>

#include "ExplicitTwisting.hpp"

#include "ControlSensor.hpp"
#include "ActuatorFactory.hpp"

ExplicitTwisting::ExplicitTwisting(SP::ControlSensor sensor, double gain, double beta):
  CommonSMC(EXPLICIT_TWISTING, sensor)
{
  _u.reset(new SiconosVector(2));
  if(beta <= 0.0 || beta >= 1.0)
  {
    std::cout << "ExplicitTwisting constructor: beta is not in (0, 1)" << std::endl;
  }

  _B.reset(new SimpleMatrix(2, 2));
  (*_B)(1, 0) = gain;
  (*_B)(1, 1) = gain*beta;

}

ExplicitTwisting::ExplicitTwisting(SP::ControlSensor sensor):
  CommonSMC(EXPLICIT_TWISTING, sensor)
{
  _u.reset(new SiconosVector(2));
}

void ExplicitTwisting::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation& s)
{
  // \TODO(Xhub) this is quite unnecessary
  CommonSMC::initialize(nsds, s);
}


ExplicitTwisting::~ExplicitTwisting()
{
}

void ExplicitTwisting::actuate()
{

  const SiconosVector& sigma = _sensor->y();

  // discontinous part
  _u->setValue(0, copysign(1., -sigma(0)));
  _u->setValue(1, copysign(1., -sigma(1)));
  *_us = *_u;
  _indx++;

}

AUTO_REGISTER_ACTUATOR(TWISTING, ExplicitTwisting)
