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
#include "InteractionLink.h"

using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

InteractionLink::InteractionLink(Interaction* itOrig, Interaction* itLinked, const std::vector<DynamicalSystem*>& dsList):
  originInteraction(itOrig), linkedInteraction(itLinked), commonDS(dsList)
{}

InteractionLink::~InteractionLink()
{
  originInteraction = NULL;
  linkedInteraction = NULL;
  commonDS.resize(1, NULL);
}

void InteractionLink::display() const
{
  cout << " ===== interactionLink display ===== " << endl;
  cout << " Origin interaction number: " << endl;
  if (originInteraction != NULL)
    cout << originInteraction->getNumber();
  else
    cout << "-> NULL" << endl;
  cout << " linked interaction number: " << endl;
  if (linkedInteraction != NULL)
    cout << linkedInteraction->getNumber();
  else
    cout << "-> NULL" << endl;

  cout << " common dynamical systems: " << endl;
  for (unsigned int i = 0; i < commonDS.size(); i++)
    commonDS[i]->display();
  cout << " =================================== " << endl;
}

// default (private) constructor
InteractionLink::InteractionLink(): originInteraction(NULL), linkedInteraction(NULL)
{}

