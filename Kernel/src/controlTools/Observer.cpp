/* Siconos-Kernel, Copyright INRIA 2005-2012.
* Siconos is a program dedicated to modeling, simulation and control
* of non smooth dynamical systems.
* Siconos is a free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
* Siconos is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Siconos; if not, write to the Free Software
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*
* Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#include "Observer.hpp"
#include "ObserverEvent.hpp"
#include "Model.hpp"
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

void Observer::initialize(const Model& m)
{
  // Get the dimension of the output
  // XXX What if there is more than one sensor ...
  if (!_sensor)
  {
    RuntimeException::selfThrow("Observer::initialize - the no ControlSensor was given");
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
