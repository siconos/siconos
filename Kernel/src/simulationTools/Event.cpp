/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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
#include "Event.h"
using namespace std;

// Default constructor (protected)
Event::Event(): timeOfEvent(DEFAULT_EVENT_TIME), type(DEFAULT_EVENT_TYPE)
{}

// copy constructor
Event::Event(const Event& newEvent): timeOfEvent(newEvent.getTimeOfEvent()), type(newEvent.getType() + "(copy)")
{}

Event::Event(const unsigned long int& time, const string& newType): timeOfEvent(time), type(newType)
{}

// Event(EventXML*, const std::string& ):timeOfEvent(0), type("undefined from xml")
//{}

Event::~Event()
{}

// const double Event::getTimeOfEvent() const
// {//todo
//   return 0.0;
// }

// void Event::setTimeOfEvent(const double& newVal)
// {//todo
// }

void Event::display() const
{
  cout << "===== Event data display =====" << endl;
  cout << " - Type: " << type << endl;
  cout << " - time (unsigned int format): " << timeOfEvent << endl;
  cout << "===== End of Event display =====" << endl;
}

//void Event::saveEventToXML()
//{
//  RuntimeException::selfThrow("saveEventToXML: not yet implemented");
//}
