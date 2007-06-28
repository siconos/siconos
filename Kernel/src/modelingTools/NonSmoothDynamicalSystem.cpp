/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "NonSmoothDynamicalSystem.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "DynamicalSystemXML.h"
#include "Topology.h"
#include "Interaction.h"
#include "LagrangianLinearTIDS.h"
#include "FirstOrderLinearTIDS.h"

using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default (private) constructor (isBvp is optional, default = false)
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const bool& isBvp):
  BVP(isBvp), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{}

// xml constuctor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(NonSmoothDynamicalSystemXML* newNsdsxml):
  BVP(false), topology(NULL), nsdsxml(newNsdsxml), isTopologyAllocatedIn(false)
{
  if (nsdsxml == NULL)
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: xml constructor, xml file=NULL");

  // === DS Vector fill-in: we sweep the list of DSXML and for each of them, add a new DynamicalSystem in the set allDS. ===
  SetOfDSXML dsList = nsdsxml->getDynamicalSystemsXML();
  SetOfDSXMLIt it;
  CheckInsertDS checkDS;
  string type;

  for (it = dsList.begin(); it != dsList.end(); ++it)
  {
    type = (*it)->getType();
    if (type  == LAGRANGIAN_NON_LINEARDS_TAG)  // LagrangianDS
      checkDS = allDS.insert(new LagrangianDS(*it , this));
    else if (type == LAGRANGIAN_TIDS_TAG)  // Lagrangian Linear Time Invariant
      checkDS = allDS.insert(new LagrangianLinearTIDS(*it, this));
    else if (type == LINEAR_DS_TAG)  // Linear DS
      checkDS = allDS.insert(new FirstOrderLinearDS(*it, this));
    else if (type == LINEAR_TIDS_TAG)  // Linear DS
      checkDS = allDS.insert(new FirstOrderLinearTIDS(*it, this));
    else if (type == NON_LINEAR_DS_TAG)  // Non linear DS
      checkDS = allDS.insert(new FirstOrderNonLinearDS(*it, this));
    else RuntimeException::selfThrow("NonSmoothDynamicalSystem::xml constructor, wrong Dynamical System type" + type);
    // checkDS.first is an iterator that points to the DS inserted into the set.
    isDSAllocatedIn[*(checkDS.first)] = true ;
  }

  // ===  The same process is applied for Interactions ===
  SetOfInteractionsXML  interactionsList = nsdsxml->getInteractionsXML();
  SetOfInteractionsXMLIt it2;
  CheckInsertInteraction checkInter;
  for (it2 = interactionsList.begin(); it2 != interactionsList.end(); ++it2)
  {
    checkInter = allInteractions.insert(new Interaction(*it2, this));
    // checkInter.first is an iterator that points to the Interaction inserted into the set.
    isInteractionAllocatedIn[*(checkInter.first)] = true;
  }

  // === Checks that sets are not empty ===
  if (allDS.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(xml, ...): the set of DS is empty.");

  if (allInteractions.isEmpty())  // Note: empty Interactions set is allowed.
    cout << "Warning: NonSmoothDynamicalSystem:: constructor(xml, ...): the set of Interactions is empty." << endl;

  // === Builds topology ===
  topology = new Topology(this);
  isTopologyAllocatedIn = true;
}

// Constructor with one DS and one Interaction (optional, default = NULL).
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(DynamicalSystem* newDS, Interaction* newInteraction, const bool& isBVP):
  BVP(isBVP), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{

  // === Checks that sets are not empty ===
  if (newDS == NULL)
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(DynamicalSystem* ds...): ds == NULL.");

  // Note that Interaction == NULL is possible and has sense.

  allDS.insert(newDS);
  isDSAllocatedIn[newDS] = false;
  if (newInteraction != NULL)
  {
    allInteractions.insert(newInteraction);
    isInteractionAllocatedIn[newInteraction] = false;
  }

  // === build topology ===
  topology = new Topology(this);
  isTopologyAllocatedIn = true;
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(DynamicalSystemsSet& listOfDS, InteractionsSet& listOfInteractions, const bool& isBVP):
  BVP(isBVP), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{

  // === Checks that sets are not empty ===
  if (listOfDS.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(DynamicalSystemsSet, ...): the set of DS is empty.");

  if (listOfInteractions.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(...,InteractionsSet, ...): the set of Interactions is empty.");

  // === "copy" listOfDS/listOfInteractions in allDS/allInteractions ===
  // Warning: this a false copy, since = operator of those objects creates links between pointers of the sets.
  allDS = listOfDS;
  allInteractions = listOfInteractions;

  // === initialize isXXAllocatedIn objects ===
  DSIterator itDS;
  InteractionsIterator itInter;
  for (itDS = allDS.begin(); itDS != allDS.end(); ++itDS)
    isDSAllocatedIn[*itDS] = false;
  for (itInter = allInteractions.begin(); itInter != allInteractions.end(); ++itInter)
    isInteractionAllocatedIn[*itInter] = false;

  // === build topology ===
  topology = new Topology(this);
  isTopologyAllocatedIn = true;
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(DynamicalSystemsSet& listOfDS, const bool& isBVP):
  BVP(isBVP), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{

  // === Checks that sets are not empty ===
  if (listOfDS.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(DynamicalSystemsSet, ...): the set of DS is empty.");

  // === "copy" listOfDS/listOfInteractions in allDS/allInteractions ===
  // Warning: this a false copy, since = operator of those objects creates links between pointers of the sets.
  allDS = listOfDS;

  // === initialize isXXAllocatedIn objects ===
  DSIterator itDS;
  InteractionsIterator itInter;
  for (itDS = allDS.begin(); itDS != allDS.end(); ++itDS)
    isDSAllocatedIn[*itDS] = false;

  // === build topology ===
  topology = new Topology(this);
  isTopologyAllocatedIn = true;
}

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  // == delete DS ==
  DSIterator it;
  for (it = allDS.begin(); it != allDS.end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }

  allDS.clear();
  isDSAllocatedIn.clear();

  // == delete Interactions ==
  InteractionsIterator it2;
  for (it2 = allInteractions.begin(); it2 != allInteractions.end(); ++it2)
  {
    if (isInteractionAllocatedIn[*it2]) delete *it2;
  }

  allInteractions.clear();
  isInteractionAllocatedIn.clear();

  if (isTopologyAllocatedIn)
  {
    delete topology;
    topology = NULL;
  }
}

// === DynamicalSystems management ===

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtr(const int& nb) const
{
  // Mind that DS are sorted in a growing order according to their id number.
  DSIterator it = allDS.begin();
  for (int i = 0; i < nb; ++i)
    it++;

  if (it == allDS.end())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getDynamicalSystemPtr(nb) : nb is out of range");

  return *it;
}

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtrNumber(const int& nb) const
{
  if (! allDS.isDynamicalSystemIn(nb)) // if ds number nb is not in the set ...
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::getDynamicalSystemOnNumber(nb), DS number nb is not in the set.");

  return allDS.getDynamicalSystemPtr(nb);
}

void NonSmoothDynamicalSystem::setDynamicalSystems(const DynamicalSystemsSet& newVect)
{
  // clear old set
  DSIterator it;
  for (it = allDS.begin(); it != allDS.end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }

  allDS.clear();
  isDSAllocatedIn.clear();

  // copy the new one
  allDS = newVect;
  for (it = allDS.begin(); it != allDS.end(); ++it)
    isDSAllocatedIn[*it] = false;

  topology->setUpToDate(false);
}

const bool NonSmoothDynamicalSystem::hasDynamicalSystemNumber(const int& nb) const
{
  return allDS.isDynamicalSystemIn(nb);
}

const bool NonSmoothDynamicalSystem::hasDynamicalSystem(DynamicalSystem* ds) const
{
  return allDS.isDynamicalSystemIn(ds);
}

// === Interactions management ===

Interaction* NonSmoothDynamicalSystem::getInteractionPtr(const int& nb) const
{
  // Mind that Interactions are sorted in a growing order according to their id number.
  InteractionsIterator it = allInteractions.begin();
  for (int i = 0; i < nb; ++i)
    it++;

  if (it == allInteractions.end())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getInteractionPtr(nb) : nb is out of range");

  return *it;
}

Interaction* NonSmoothDynamicalSystem::getInteractionPtrNumber(const int& nb) const
{
  if (! allInteractions.isInteractionIn(nb)) // if Interaction number nb is not in the set ...
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::getInteractionOnNumber(nb), Interaction number nb is not in the set.");

  return allInteractions.getInteraction(nb);
}

void NonSmoothDynamicalSystem::setInteractions(const InteractionsSet& newVect)
{
  // clear old set
  InteractionsIterator it;
  for (it = allInteractions.begin(); it != allInteractions.end(); ++it)
  {
    if (isInteractionAllocatedIn[*it]) delete *it;
  }

  allInteractions.clear();
  isInteractionAllocatedIn.clear();

  // copy the new one
  allInteractions = newVect;
  for (it = allInteractions.begin(); it != allInteractions.end(); ++it)
    isInteractionAllocatedIn[*it] = false;

  topology->setUpToDate(false);
}

const bool NonSmoothDynamicalSystem::hasInteractionNumber(const int& nb) const
{
  return allInteractions.isInteractionIn(nb);
}

const bool NonSmoothDynamicalSystem::hasInteraction(Interaction* inter) const
{
  return allInteractions.isInteractionIn(inter);
}

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
  if (nsdsxml != NULL)
  {
    nsdsxml->setBVP(BVP);// no need to change the value of BVP, it mustn't change anyway

    DSIterator it;
    for (it = allDS.begin(); it != allDS.end(); ++it)
    {
      if ((*it)->getType() == LNLDS)
        (static_cast<LagrangianDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == LLTIDS)
        (static_cast<LagrangianLinearTIDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == FOLDS)
        (static_cast<FirstOrderLinearDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == FOLTIDS)
        (static_cast<FirstOrderLinearDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == FONLDS)
        (*it)->saveDSToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of DynamicalSystem");
    }

    InteractionsIterator it2;
    for (it2 = allInteractions.begin(); it2 != allInteractions.end(); ++it2)
      (*it2)->saveInteractionToXML();
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NonSmoothDynamicalSystemXML object doesn't exists");
}

void NonSmoothDynamicalSystem::display() const
{
  cout << " ===== Non Smooth Dynamical System display ===== " << endl;
  cout << "---> isBVP = " << BVP << endl;
  allDS.display();
  allInteractions.display();
  cout << "===================================================" << endl;
}

void NonSmoothDynamicalSystem::addDynamicalSystemPtr(DynamicalSystem *ds)
{
  ds->setNonSmoothDynamicalSystemPtr(this);
  allDS.insert(ds);
  isDSAllocatedIn[ds] = false;
  topology->setUpToDate(false);
}

void NonSmoothDynamicalSystem::addInteractionPtr(Interaction *inter)
{
  inter->setNonSmoothDynamicalSystemPtr(this);
  allInteractions.insert(inter);
  isInteractionAllocatedIn[inter] = false;
  // the topology should be updated
  topology->setUpToDate(false);
}

double NonSmoothDynamicalSystem::nsdsConvergenceIndicator()
{
  // calculate the max value of all DS convergence indicators
  double convergenceIndicator = (* allDS.begin()) -> dsConvergenceIndicator();
  double dsIndic ;
  DSIterator iter;
  for (iter = allDS.begin(); iter != allDS.end(); ++iter)
  {
    dsIndic = (*iter)->dsConvergenceIndicator();
    if (dsIndic > convergenceIndicator) convergenceIndicator = dsIndic;
  }
  return(convergenceIndicator);
}

