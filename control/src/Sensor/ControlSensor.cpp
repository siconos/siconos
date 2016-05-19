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

/*! \file ControlSensor.cpp
 * A generic control sensor
*/

#include "ControlSensor.hpp"
#include <cmath>
#include "TimeDiscretisation.hpp"
#include "SiconosVector.hpp"

void ControlSensor::initialize(const Model& m)
{
  Sensor::initialize(m);
//  if (_delay > 0)
//  {
//    if (_timeDiscretisation->getTDCase() != 2)
//    {
//       RuntimeException::selfThrow("ControlSensor::initialize the timediscretization should be of type 2");
//    }
//    else
//    {
//       double h = _timeDiscretisation->currentTimeStep();
//       double shift = _delay;
//       if (_delay >= h)
//       {
//         unsigned int size = ceil(_delay/h);
//         shift = fmod(_delay, h);
//         _bufferY.resize(size);
//       }
//      _timeDiscretisation->setT0(_timeDiscretisation->currentTime() + h - shift);
//    }
//  }
//  else if (_delay < 0)
//    RuntimeException::selfThrow("ControlSensor::initialize the delay value should be >= 0");
}

unsigned int ControlSensor::getYDim() const
{
  return _storedY->size();
}
