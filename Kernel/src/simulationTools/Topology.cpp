/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "Topology.h"
#include "NonSmoothLaw.h"
#include "NonSmoothDynamicalSystem.h"
#include "Interaction.h"
#include "UnitaryRelation.h"

using namespace std;

const bool Topology::addInteractionInIndexSet(SP::Interaction inter)
{
  // Private function
  //
  // Creates UnitaryRelations corresponding to inter and add them into indexSet0

  // First, we get the number of relations in the interaction.
  // This corresponds to inter->getNumberOfRelations but since Interaction has not
  // been initialized yet, this value is not set and we need to get interaction size and nsLaw size.
  unsigned int nsLawSize = inter->getNonSmoothLawPtr()->getNsLawSize();
  unsigned int m = inter->getSizeOfY() / nsLawSize;
  unsigned int pos = 0; // relative position of the relation in the y vector of the Interaction
  CheckInsertUnitaryRelation checkUR;
  bool res = true; // output value. False if insertion of one of the relations fails.
  for (unsigned int i = 0; i < m; ++i)
  {
    // each UnitaryRelation is of size "nsLawSize", at position pos and of number i.
    checkUR = indexSet0->insert(SP::UnitaryRelation(new UnitaryRelation(inter, pos, i)));
    pos += nsLawSize;
    if (checkUR.second == false) res = false;
  }

  numberOfConstraints += m * nsLawSize;
  return res;
}

// Compute relative degrees map
void Topology::computeRelativeDegrees()
{
  // for each Unitary Relation relative degree vector depends on NonSmooth Law and relations

  relativeDegrees.clear();
  string nslawType;

  // loop through indexSet0
  UnitaryRelationsIterator it;
  for (it = indexSet0->begin(); it != indexSet0->end(); it++)
  {
    nslawType = (*it)->getNonSmoothLawType();
    if (nslawType == COMPLEMENTARITYCONDITIONNSLAW || nslawType == MIXEDCOMPLEMENTARITYCONDITIONNSLAW)
      relativeDegrees[*it] = 0;

    else if (nslawType == NEWTONIMPACTNSLAW)
    {
      relativeDegrees[*it] = 2;
      isTopologyTimeInvariant = false;
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {
      relativeDegrees[*it] = 2;
      isTopologyTimeInvariant = false;
    }
    else
      RuntimeException::selfThrow("Topology::computeRelativeDegree(...), not yet implemented for non smooth law of type" + nslawType);
  }
}

// --- CONSTRUCTORS/DESTRUCTOR ---

// default
Topology::Topology(): isTopologyUpToDate(false), isTopologyTimeInvariant(true), numberOfConstraints(0)
{}

// with InteractionSet
Topology::Topology(SP::InteractionsSet newInteractionSet): isTopologyUpToDate(false), isTopologyTimeInvariant(true), numberOfConstraints(0)
{
  allInteractions = newInteractionSet;
  indexSet0.reset(new UnitaryRelationsSet());
}

// destructor
Topology::~Topology()
{
}

const bool Topology::hasInteraction(SP::Interaction inter) const
{
  return allInteractions->isIn(inter);
}

const unsigned int Topology::getMaxRelativeDegree()
{
  if (relativeDegrees.empty())
    RuntimeException::selfThrow("Topology::getMaxRelativeDegree, non-existent value, since the relative degrees map is empty.");

  ConstIteratorForRelativeDegrees it = max_element(relativeDegrees.begin(), relativeDegrees.end());
  return(it->second);
}

const unsigned int Topology::getMinRelativeDegree()
{
  if (relativeDegrees.empty())
    RuntimeException::selfThrow("Topology::getMinRelativeDegree, non-existent value, since the relative degrees map is empty.");

  ConstIteratorForRelativeDegrees it = min_element(relativeDegrees.begin(), relativeDegrees.end());
  return(it->second);
}

void Topology::initialize()
{

  assert(allInteractions && "Topology : allInteractions is NULL");

  // -- Creates Unitary Relations and put them in indexSet0 ---
  // loop through interactions list (from NSDS)
  indexSet0->clear();
  InteractionsIterator it;
  for (it = allInteractions->begin()  ; it != allInteractions->end() ; ++it)
    addInteractionInIndexSet(*it);

  //-- Fills RelativeDegreesMaps in --
  computeRelativeDegrees();

  isTopologyUpToDate = true;
}

