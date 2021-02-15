/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "PID.hpp"

#include "FirstOrderLinearTIDS.hpp"
#include "SiconosVector.hpp"
#include "ActuatorFactory.hpp"
#include "TimeDiscretisation.hpp"
#include "ControlSensor.hpp"
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "EventsManager.hpp"

#include <iostream>
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"


PID::PID(SP::ControlSensor sensor, SP::SimpleMatrix B): Actuator(PID_, sensor, B), _ref(0), _curDeltaT(0)
{
  _u.reset(new SiconosVector(1, 0));
}

PID::~PID()
{
}

void PID::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation& s)
{
  Actuator::initialize(nsds,s);

  _curDeltaT = s.currentTimeStep();

  // initialize _err
  _err.reset(new boost::circular_buffer<double> (3));
  for(unsigned int i = 0; i < 3; ++i)
    (*_err).push_front(0.0);
}

void PID::actuate()
{
  DEBUG_BEGIN("void PID::actuate()\n");
  /** \todo We have to distinguish two cases : linear or nonlinear
   *  support the nonlinear case
   */

  // Compute the new error

  (*_err).push_front(_ref - _sensor->y()(0));
  DEBUG_PRINTF("_curDeltaT = %g\n", _curDeltaT);
  DEBUG_PRINTF("_ref = %g\n", _ref);
  DEBUG_EXPR(_sensor->y().display(););
  DEBUG_EXPR(_u->display(););
  DEBUG_PRINTF("added term  = %g\n", ((*_K)(0) + (*_K)(2) / _curDeltaT + (*_K)(1) * _curDeltaT) * (*_err)[0] +
               (-(*_K)(0) - 2 * (*_K)(2) / _curDeltaT) * (*_err)[1] + (*_K)(2) / _curDeltaT * (*_err)[2]);
  // compute the new control and update it
  (*_u)(0) += ((*_K)(0) + (*_K)(2) / _curDeltaT + (*_K)(1) * _curDeltaT) * (*_err)[0] +
              (-(*_K)(0) - 2 * (*_K)(2) / _curDeltaT) * (*_err)[1] + (*_K)(2) / _curDeltaT * (*_err)[2];
  DEBUG_EXPR(_u->display(););
  DEBUG_END("void PID::actuate()\n");
}

void PID::setK(SP::SiconosVector K)
{
  // check dimensions ...
  if(K->size() != 3)
  {
    THROW_EXCEPTION("PID::setK - the size of K should be 3");
  }
  else
  {
    _K = K;
  }
}

void PID::setTimeDiscretisation(const TimeDiscretisation& td)
{
}

void PID::display() const
{
  Actuator::display();
  std::cout << "current error vector: ";
  std::cout << (*_err)[0] << " "  << (*_err)[1] << " " << (*_err)[2] << std::endl;
}
AUTO_REGISTER_ACTUATOR(PID_, PID)
