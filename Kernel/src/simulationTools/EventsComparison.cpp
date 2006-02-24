/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#include "EventsComparison.h"
using namespace std;

bool EventsComparison::operator()(Event* e1, Event* e2) const
{
  //  bool result = 0;
  unsigned long int t1 = e1->getTimeOfEvent();
  unsigned long int t2 = e2->getTimeOfEvent();
  // if time are different, sort according to time ...
  if (t1 == t2)
  {
    //      bool comp = 0;
    // Todo: Set type name in order to be consistent with string operator ">" ??
    string type1 = e1->getType();
    string type2 = e2->getType();
    return type1 > type2; // This means TimeDiscretisationEvent is set as anterior to NonSmoothEvent
  }
  else
    return mode == normal ? t1 < t2 : t2 < t1;
}
