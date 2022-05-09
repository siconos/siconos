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

#include "Observer.hpp"
#include "ObserverEvent.hpp"
#include "TimeDiscretisation.hpp"
#include "Simulation.hpp"
#include "EventsManager.hpp"
#include <iostream>


Observer::Observer(): _type(0), _id("none")
{
}

Observer::Observer(unsigned int type, SP::ControlSensor sensor, const SiconosVector& xHat0, const std::string& newId):
  _type(type), _sensor(sensor), _id(newId)
{
  _xHat.reset(new SiconosVector(xHat0));
}

Observer::Observer(unsigned int type, SP::ControlSensor sensor, const SiconosVector& xHat0, SP::DynamicalSystem ds, const std::string& newId):
  _type(type), _DS(ds), _sensor(sensor), _id(newId)
{
  _xHat.reset(new SiconosVector(xHat0));
}

Observer::~Observer()
{
}

void Observer::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation& s)
{
  // Get the dimension of the output
  // XXX What if there is more than one sensor ...
  if(!_sensor)
  {
    THROW_EXCEPTION("Observer::initialize - the no ControlSensor was given");
  }
}

void Observer::display() const
{
  std::cout << "=====> Observer of type " << _type << ", named " << _id ;
  std::cout << std::endl;
}

void Observer::setTimeDiscretisation(const TimeDiscretisation& td)
{
  _td.reset(new TimeDiscretisation(td));
}
