/* Siconos version 1.0, Copyright INRIA 2005.
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

#include "EventDriven.h"
using namespace std;

// --- Default constructor ---
EventDriven::EventDriven(): Strategy()
{
  strategyType = EVENTDRIVEN_STRATEGY;
}

// --- XML constructor ---
EventDriven::EventDriven(StrategyXML* strxml, Model *newModel): Strategy(strxml, newModel)
{
  strategyType = EVENTDRIVEN_STRATEGY;
}

// --- Destructor ---
EventDriven::~EventDriven()
{}

void EventDriven::createStrategy(StrategyXML * newStrategyXML, Model * newModel)
{
  strategyxml = NULL;
  strategyType = EVENTDRIVEN_STRATEGY;
  model = newModel;
}


EventDriven* EventDriven::convert(Strategy *str)
{
  cout << "EventDriven::convert (Strategy *str)" << endl;
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}

