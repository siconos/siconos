/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "UnitaryRelation.h"

using namespace std;

// --- CONSTRUCTORS ---

// Default (private) constructor
UnitaryRelation::UnitaryRelation(): mainInteraction(NULL), relativePosition(0)
{}

// Copy constructor
UnitaryRelation::UnitaryRelation(const UnitaryRelation& newI): mainInteraction(NULL), relativePosition(0)
{
}

// --- DESTRUCTOR ---
UnitaryRelation::~UnitaryRelation()
{
  mainInteraction = NULL;
}

SimpleVector* UnitaryRelation::getYPtr(const unsigned int& i) const
{
  //unsigned int nsLawSize = mainInteraction->getNonSmoothLawPtr()->getNsLawSize();

  // get the block that corresponds to the current relation in vector y of the main Interaction
  // Find a way to avoid copy? Modify BlockVector?

  return mainInteraction->getYPtr(1);
}

SimpleVector* UnitaryRelation::getLambdaPtr(const unsigned int& i) const
{
  //unsigned int nsLawSize = mainInteraction->getNonSmoothLawPtr()->getNsLawSize();

  // get the block that corresponds to the current relation in vector lambda of the main Interaction
  return mainInteraction->getLambdaPtr(1);
}
