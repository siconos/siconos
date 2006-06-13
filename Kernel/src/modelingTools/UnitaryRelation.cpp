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
#include "Interaction.h"

using namespace std;

// --- CONSTRUCTORS ---

// Default (private) constructor
UnitaryRelation::UnitaryRelation(): mainInteraction(NULL), relativePosition(0)
{}

// Copy constructor
UnitaryRelation::UnitaryRelation(const UnitaryRelation& newUR): mainInteraction(NULL), relativePosition(0)
{
  // Copy of Unitary relation is not allowed. Since they are identified/sorted in UnitaryRelationsSet thanks to their address (as pointers)
  // a copy could results in two differents objects that will correspond to the same relation.
  RuntimeException::selfThrow("UnitaryRelation copy constructor call: forbidden operation.");
}

// Data constructor
UnitaryRelation::UnitaryRelation(Interaction* inter, const unsigned int& pos): mainInteraction(inter), relativePosition(pos)
{}

// --- DESTRUCTOR ---
UnitaryRelation::~UnitaryRelation()
{
  mainInteraction = NULL;
}

const double UnitaryRelation::getY(const unsigned int& i) const
{
  // get the single value used to build indexSets
  // Warning: the relativePosition depends on NsLawSize and/or type.
  // This means that at the time, for the block of y that corresponds to the present relation, the first scalar value is used.
  // For example, for friction, normal part is in first position, followed by the tangential parts.
  return (*(mainInteraction->getYPtr(i)))(relativePosition);
}

const double UnitaryRelation::getLambda(const unsigned int& i) const
{
  // get the single value used to build indexSets
  return (*(mainInteraction->getLambdaPtr(1)))(relativePosition);
}
