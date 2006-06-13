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
#include "Topology.h"

using namespace std;

const bool Topology::addInteractionInIndexSet(Interaction * inter)
{
  unsigned int m = inter->getInteractionSize() ;
  unsigned int nsLawSize = inter->getNonSmoothLawPtr()->getNsLawSize();
  unsigned int pos = 0; // relative position of the relation in the y vector of the Interaction
  CheckInsertUnitaryRelation checkUR;
  bool res = true; // output value. False if insertion of one of the relations fails.
  for (unsigned int i = 0; i < m; ++i)
  {
    checkUR = indexSets[0].insert(new UnitaryRelation(inter, pos));
    pos = pos + nsLawSize;
    if (checkUR.second == false) res = false;
  }
  return res;
}

// --- CONSTRUCTORS/DESTRUCTOR ---

// default
Topology::Topology():
  effectiveSizeOutput(0), isTopologyUpToDate(false), isTopologyTimeInvariant(true), nsds(NULL)
{}

// with NonSmoothDynamicalSystem
Topology::Topology(NonSmoothDynamicalSystem* newNsds):
  effectiveSizeOutput(0), isTopologyUpToDate(false), isTopologyTimeInvariant(true), nsds(newNsds)
{
  // Creates indexSets vector. Its size depends on the relative degree of the system
  // For this first draft, we set size == 3

  indexSets.resize(3);

  InteractionsSet interSet = nsds->getInteractions() ;

  // -- Creates Unitary Relations and put them in indexSets[0] ---
  // loop through interactions list (from NSDS)
  InteractionsIterator it;
  for (it = interSet.begin()  ; it != interSet.end(); ++it)
    addInteractionInIndexSet(*it);

  updateTopology();
}

// destructor
Topology::~Topology()
{
  // Clears all Unitary Relations of IndexSets[0]
  UnitaryRelationIterator it;
  for (it = indexSets[0].begin(); it != indexSets[0].end(); ++it)
    delete *it;


  map< Interaction*, std::vector<InteractionLink*> >::iterator mapIt;
  // for all elements in the map ...
  for (mapIt = linkedInteractionMap.begin(); mapIt != linkedInteractionMap.end(); mapIt++)
  {
    // get vector of interactionLink
    vector<InteractionLink*> linkVector = (*mapIt).second;
    vector<InteractionLink*>::iterator linkVectorIt;
    // delete each interactionLink
    for (linkVectorIt = linkVector.begin(); linkVectorIt != linkVector.end(); linkVectorIt++)
    {
      delete(*linkVectorIt);
      *linkVectorIt = NULL ;
    }
    linkVector.clear();
  }
  linkedInteractionMap.clear();
  nsds = NULL;
}


Interaction* Topology::getInteractionPtrNumber(const int& nb) const
{
  if (! allInteractions.isInteractionIn(nb)) // if Interaction number nb is not in the set ...
    RuntimeException::selfThrow("Topology::getInteractionOnNumber(nb), Interaction number nb is not in the set.");

  return allInteractions.getInteraction(nb);
}

void Topology::setInteractions(const InteractionsSet& newVect)
{
  // clear old set
  allInteractions.clear();
  // copy the new one
  allInteractions = newVect;
}

const bool Topology::hasInteraction(Interaction* inter) const
{
  return allInteractions.isInteractionIn(inter);
}

void Topology::updateTopology()
{
  //-- Get all the interactions --
  InteractionsSet listInteractions = nsds->getInteractions();

  //-- Fill RelativeDegreesMaps in --
  computeRelativeDegreesMap();

  //-- Fill indexMinMap and indexMaxMap --
  computeIndexMinMap();
  // initialization of indexMax and effectiveIndexes
  indexMaxMap = indexMinMap; // by default, indexMax = indexMin

  // Note: all 'effective' values depend on time and should be updated during computation ('computeEffectiveOutput'
  // in OneStepNSProblem, this as soon as one relative degree differs from 0 or 1.
  // Thus, all default values given in the current function are those corresponding to r=0 or 1.


  // -- Compute sizeOutput  ( ie sum of all interactions sizes) --
  InteractionsIterator it;
  // loop through interactions list
  unsigned int size;
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
  {
    size = (*it)->getInteractionSize();

    // initialization of indexMax and effectiveIndexes

    // by default, all relations are effective
    effectiveIndexesMap[*it].resize(size);
    for (unsigned int i = 0; i < size; i++)
      effectiveIndexesMap[*it][i] = i;
  }

  // fill linkedInteractionMap and interactionPositionMap in:
  computeLinkedInteractionMap();
  //  computeInteractionPositionMap();
  computeInteractionEffectivePositionMap();

  if (isTopologyTimeInvariant)
    computeEffectiveSizeOutput();

  // else this value will be computed in OneStepNSProblem
  // after indexMap computation.

  isTopologyUpToDate = true;

}

void Topology::computeEffectiveSizeOutput()
{
  effectiveSizeOutput = 0;

  // Get all the interactions
  InteractionsSet listInteractions = nsds->getInteractions();

  InteractionsIterator it;
  // loop over interactions list
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
  {
    // loop over ouput vector
    unsigned int j;
    for (j = 0; j < (*it)->getInteractionSize(); j++)
      effectiveSizeOutput += (indexMaxMap[*it])[j] - (indexMinMap[*it])[j] + 1 ;
  }
}

unsigned int Topology::computeEffectiveSizeOutput(Interaction * inter)
{
  unsigned int sizeOutput = 0;
  unsigned int j;
  // loop over ouput vector
  for (j = 0; j < inter->getInteractionSize(); j++)
    sizeOutput += indexMaxMap[inter][j] - indexMinMap[inter][j] + 1 ;
  return sizeOutput;
}

// Linked interactions map computing:
void Topology::computeLinkedInteractionMap()
{
  DSSet dsOrig, dsLinked;

  // Get all the interactions
  InteractionsSet listInteractions = nsds->getInteractions();

  linkedInteractionMap.clear();

  InteractionsIterator itOrig;

  // -- loop over all the interactions of the NonSmoothDynamicalSystem --
  for (itOrig = listInteractions.begin(); itOrig != listInteractions.end(); itOrig++)
  {
    dsOrig = (*itOrig)->getDynamicalSystems();

    InteractionsIterator itLinked;
    // -- check all other interactions --
    for (itLinked = listInteractions.begin(); itLinked != listInteractions.end() && itLinked != itOrig; itLinked++)
    {
      // list of ds of the second interaction
      dsLinked = (*itLinked)->getDynamicalSystems();

      DSSet commonDS = intersection(dsOrig, dsLinked); // list of ds common to both interactions

      // built linkedInteractionMap
      if (commonDS.size() != 0)
        linkedInteractionMap[(*itOrig) ].push_back(new InteractionLink((*itOrig), (*itLinked), commonDS));
    }
  }
}

// interactionEffectivePosition map computing:
// in the 'effective' matrix, the size of each bloc (ie for a specific interaction) is the sum of indexMax(j)-indexMin(j)+1.
// j the relation number  => this value depends on time is should be updated during computation
void Topology::computeInteractionEffectivePositionMap()
{
  interactionEffectivePositionMap.clear();

  // Get all the interactions
  InteractionsSet listInteractions = nsds->getInteractions();

  InteractionsIterator it;
  unsigned int currentEffectivePosition = 0;

  // loop through interactions list
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
  {
    interactionEffectivePositionMap[*it] = currentEffectivePosition;
    currentEffectivePosition += computeEffectiveSizeOutput(*it);
  }
}

// Compute relative degrees map
void Topology::computeRelativeDegreesMap()
{
  // for each interaction relative degree vector depends on NonSmooth Law and relations

  relativeDegreesMap.clear();

  // Get all the interactions
  InteractionsSet listInteractions = nsds->getInteractions();

  // loop through interactions list
  InteractionsIterator it;
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
    relativeDegreesMap[(*it) ] = computeRelativeDegrees(*it) ;
}

// Compute relative degree of a single interaction
vector<unsigned int> Topology::computeRelativeDegrees(Interaction * inter)
{
  // get NonSmooth law, relation and their types
  NonSmoothLaw * nslaw = inter ->getNonSmoothLawPtr();
  string nslawType = nslaw -> getType();
  Relation * relation  = inter ->getRelationPtr();
  string relationType = relation->getType();

  vector<unsigned int> relativeDegree;
  unsigned int sizeInter = inter->getInteractionSize();

  // loop over various non smooth law types
  // \todo check isTimeInvariant properly
  if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
  {
    relativeDegree.resize(sizeInter, 0);
    // \todo compute properly
  }
  else if (nslawType == NEWTONIMPACTNSLAW)
  {
    relativeDegree.resize(sizeInter, 2);
    isTopologyTimeInvariant = false;
  }
  else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
  {
    relativeDegree.resize(sizeInter, 2);
    isTopologyTimeInvariant = false;
  }
  else
    RuntimeException::selfThrow("Topology::computeRelativeDegrees, not yet implemented for non smooth law of type" + nslawType);

  return relativeDegree;
}

// Compute indexMin map
void Topology::computeIndexMinMap()
{
  // for each interaction indexMin depends on NonSmooth Law

  indexMinMap.clear();

  // Get all the interactions
  InteractionsSet listInteractions = nsds->getInteractions();

  // loop through interactions list
  InteractionsIterator it;
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
    indexMinMap[(*it) ] = computeIndexMin(*it);

}

// Compute relative degree of a single interaction
vector<unsigned int> Topology::computeIndexMin(Interaction * inter)
{
  // get NonSmooth law, relation and their types
  NonSmoothLaw * nslaw = inter ->getNonSmoothLawPtr();
  string nslawType = nslaw -> getType();
  Relation * relation  = inter ->getRelationPtr();
  string relationType = relation->getType();

  vector<unsigned int> indexMin;
  unsigned int sizeInter = inter->getInteractionSize();

  // loop over various non smooth law types
  if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
    indexMin.resize(sizeInter, 0);

  else if (nslawType == NEWTONIMPACTNSLAW)
    indexMin.resize(sizeInter, 1);

  else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    indexMin.resize(sizeInter, 1);

  else
    RuntimeException::selfThrow("Topology::computeIndexMin, not yet implemented for non smooth law of type" + nslawType);

  return indexMin;
}

void Topology::updateIndexSets()
{
  for (unsigned int i = 1; i < indexSets.size() ; ++i)
    updateIndexSet(i);
}


void Topology::updateIndexSet(const unsigned int& i)
{
  if (i > indexSets.size())
    RuntimeException::selfThrow("Topology::updateIndexSet(i), indexSets[i] does not exist.");

  if (i == 0) // IndexSets[0] must not be updated by this function.
    RuntimeException::selfThrow("Topology::updateIndexSet(i=0), indexSets[0] can not be updated.");

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and update the indexSet[i]
  UnitaryRelationIterator it, itForFind;

  for (it = indexSets[i - 1].begin(); it != indexSets[i - 1].end(); ++it)
  {
    // check if current Unitary Relation (ie *it) is in indexSets[i]
    // if not itForFind = indexSets.end()
    itForFind = indexSets[i].find(*it);

    // if y[i-1] <=0, then the unitary relation is added in indexSets[i] (if it was not already there)
    // else if y[i-1] > 0 and if the unitary relation was in the set, it is removed.
    if ((*it)->getY(i - 1) <= 0 && itForFind == indexSets[i].end())
      indexSets[i].insert(*it);
    else if ((*it)->getY(i - 1) > 0 && itForFind != indexSets[i].end())
      indexSets[i].erase(*it);
  }
}

void Topology::updateIndexSetsWithDoubleCondition()
{

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and update the indexSet[i]
  UnitaryRelationIterator it, itForFind;

  for (it = indexSets[2].begin(); it != indexSets[2].end(); ++it)
  {
    double gamma = (*it)->getY(2);
    double F     = (*it)->getLambda(2);

    if (gamma > 0 && F < TOLERANCE)
      indexSets[2].erase(*it);
    else if (gamma < TOLERANCE && F < TOLERANCE) // undetermined case
      RuntimeException::selfThrow("Topology::updateIndexSetsWithDoubleCondition(), undetermined case.");
  }
}

