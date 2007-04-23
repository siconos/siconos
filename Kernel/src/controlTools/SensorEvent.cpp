/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "SensorEvent.h"
#include "Sensor.h"
#include "EventFactory.h"
using namespace std;
using namespace EventFactory;

// Default constructor
SensorEvent::SensorEvent(): Event(0.0, "SensorEvent")
{}

SensorEvent::SensorEvent(double time, const std::string& name): Event(time, name)
{}

SensorEvent::~SensorEvent()
{}

void SensorEvent::process(Simulation*)
{
  sensor->capture();
}

AUTO_REGISTER_EVENT("SensorEvent", SensorEvent);
