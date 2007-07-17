/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

const bool Topology::addInteractionInIndexSet(Interaction * inter)
{
  // Creates UnitaryRelations corresponding to inter and add them into indexSet0

  // get number of relations in the interaction. This corresponds to inter->getNumberOfRelations but since Interaction has not
  // been initialized yet, this value is not set and we need to get interaction size and nsLaw size.
  unsigned int m = inter->getSizeOfY() / inter->getNonSmoothLawPtr()->getNsLawSize() ;
  unsigned int nsLawSize = inter->getNonSmoothLawPtr()->getNsLawSize();
  unsigned int pos = 0; // relative position of the relation in the y vector of the Interaction
  CheckInsertUnitaryRelation checkUR;
  bool res = true; // output value. False if insertion of one of the relations fails.
  for (unsigned int i = 0; i < m; ++i)
  {
    // each UnitaryRelation is of size "nsLawSize", at position pos and of number i.
    checkUR = indexSet0.insert(new UnitaryRelation(inter, pos, i));
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
  UnitaryRelationIterator it;
  for (it = indexSet0.begin(); it != indexSet0.end(); it++)
  {
    nslawType = (*it)->getNonSmoothLawType();
    if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
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
Topology::Topology(): isTopologyUpToDate(false), isTopologyTimeInvariant(true), nsds(NULL), numberOfConstraints(0)
{}

// with NonSmoothDynamicalSystem
Topology::Topology(NonSmoothDynamicalSystem* newNsds): isTopologyUpToDate(false), isTopologyTimeInvariant(true), nsds(newNsds), numberOfConstraints(0)
{}

// destructor
Topology::~Topology()
{
  // Clears all Unitary Relations of IndexSets[0]
  UnitaryRelationIterator it;
  for (it = indexSet0.begin(); it != indexSet0.end(); ++it)
    delete *it;

  nsds = NULL;
}

const bool Topology::hasInteraction(Interaction* inter) const
{
  return allInteractions.isInteractionIn(inter);
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
    RuntimeException::selfThrow("Topology::getMaxRelativeDegree, non-existent value, since the relative degrees map is empty.");

  ConstIteratorForRelativeDegrees it = min_element(relativeDegrees.begin(), relativeDegrees.end());
  return(it->second);
}

void Topology::initialize()
{
  if (nsds == NULL)
    RuntimeException::selfThrow("Topology::initialize, the topology is not linked to a non smooth dynamical system.");

  //-- Get all interactions --
  allInteractions = nsds->getInteractions() ;

  // -- Creates Unitary Relations and put them in indexSet0 ---
  // loop through interactions list (from NSDS)
  indexSet0.clear();
  InteractionsIterator it;
  for (it = allInteractions.begin()  ; it != allInteractions.end() ; ++it)
    addInteractionInIndexSet(*it);

  //-- Fills RelativeDegreesMaps in --
  computeRelativeDegrees();

  isTopologyUpToDate = true;
}

