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

#include "Sensor.hpp"
#include "SensorEvent.hpp"
#include "DynamicalSystem.hpp"
#include "TimeDiscretisation.hpp"
#include "EventFactory.hpp"
#include "Simulation.hpp"
#include <iostream>


Sensor::Sensor(): _type(0), _id("none")
{}

Sensor::Sensor(unsigned int type, SP::TimeDiscretisation t, SP::DynamicalSystem ds): _type(type), _id("none"), _DS(ds), _timeDiscretisation(t)
{
  _nDim = _DS->getN();
}

Sensor::~Sensor()
{}

void Sensor::initialize(const Model& m)
{
  // == Create an event linked to the present Sensor. ==
  Event& ev = m.simulation()->eventsManager()->insertEvent(SENSOR_EVENT, _timeDiscretisation->currentTime());
  static_cast<SensorEvent&>(ev).setSensorPtr(shared_from_this());
}

void Sensor::display() const
{
  std::cout << "=====> Sensor of type " << _type << ", named " << _id ;
  if (_DS)
    std::cout << " and linked to the DynamicalSystem number " << _DS->number() << "." <<std::endl;
  else
    std::cout << " and not linked to a DynamicalSystem." <<std::endl;
  std::cout << "======" <<std::endl ;
  std::cout <<std::endl;
}
