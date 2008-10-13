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
#include "NonSmoothDynamicalSystem.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "DynamicalSystemXML.h"
#include "Topology.h"
#include "Interaction.h"
#include "LagrangianLinearTIDS.h"
#include "FirstOrderLinearTIDS.h"

using namespace std;
using namespace DS;
using namespace RELATION;

// --- CONSTRUCTORS/DESTRUCTOR ---

// xml constuctor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(SP::NonSmoothDynamicalSystemXML newNsdsxml):
  BVP(false), nsdsxml(newNsdsxml)
{
  assert(nsdsxml && "NonSmoothDynamicalSystem:: xml constructor, xml file=NULL");

  // === DS Vector fill-in: we sweep the list of DSXML and for each of them, add a new DynamicalSystem in the set allDS. ===
  SetOfDSXML dsList = nsdsxml->getDynamicalSystemsXML();
  SetOfDSXML::iterator it;
  CheckInsertDS checkDS;
  DS::TYPES type;
  allDS.reset(new DynamicalSystemsSet());
  for (it = dsList.begin(); it != dsList.end(); ++it)
  {
    type = (*it)->getType();
    if (type  == LNLDS)  // LagrangianDS
      checkDS = allDS->insert(SP::LagrangianDS(new LagrangianDS(*it)));
    else if (type == LLTIDS)  // Lagrangian Linear Time Invariant
      checkDS = allDS->insert(SP::LagrangianLinearTIDS(new LagrangianLinearTIDS(*it)));
    else if (type == FOLDS)  // Linear DS
      checkDS = allDS->insert(SP::FirstOrderLinearDS(new FirstOrderLinearDS(*it)));
    else if (type == FOLTIDS)  // Linear Time Invariant DS
      checkDS = allDS->insert(SP::FirstOrderLinearTIDS(new FirstOrderLinearTIDS(*it)));
    else if (type == FONLDS)  // Non linear DS
      checkDS = allDS->insert(SP::FirstOrderNonLinearDS(new FirstOrderNonLinearDS(*it)));
    else RuntimeException::selfThrow("NonSmoothDynamicalSystem::xml constructor, wrong Dynamical System type" + type);
    // checkDS.first is an iterator that points to the DS inserted into the set.
  }

  // ===  The same process is applied for Interactions ===
  SetOfInteractionsXML  interactionsList = nsdsxml->getInteractionsXML();
  SetOfInteractionsXMLIt it2;
  CheckInsertInteraction checkInter;
  allInteractions.reset(new InteractionsSet());
  for (it2 = interactionsList.begin(); it2 != interactionsList.end(); ++it2)
  {
    checkInter = allInteractions->insert(SP::Interaction(new Interaction(*it2, allDS)));
    // checkInter.first is an iterator that points to the Interaction inserted into the set.
  }

  // === Checks that sets are not empty ===
  if (allDS->isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(xml, ...): the set of DS is empty.");

  if (allInteractions->isEmpty())  // Note: empty Interactions set is allowed.
    cout << "Warning: NonSmoothDynamicalSystem:: constructor(xml, ...): the set of Interactions is empty." << endl;

  // === Builds topology ===
  topology.reset(new Topology(allInteractions));
}

// Constructor with one DS and one Interaction (optional)
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(SP::DynamicalSystem newDS, SP::Interaction newInteraction, const bool& isBVP):
  BVP(isBVP)
{
  // === Checks that sets are not empty ===
  if (!newDS)
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(SP::DynamicalSystem ds...): ds == NULL.");

  // Note that Interaction == NULL is possible and has sense.

  allDS.reset(new DynamicalSystemsSet());
  allInteractions.reset(new InteractionsSet());
  allDS->insert(newDS);
  if (newInteraction)
  {
    allInteractions->insert(newInteraction);
  }

  // === build topology ===
  topology.reset(new Topology(allInteractions));
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(DynamicalSystemsSet& listOfDS, InteractionsSet& listOfInteractions, const bool& isBVP):
  BVP(isBVP)
{

  // === Checks that sets are not empty ===
  if (listOfDS.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(DynamicalSystemsSet, ...): the set of DS is empty.");

  if (listOfInteractions.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(...,InteractionsSet, ...): the set of Interactions is empty.");

  // === "copy" listOfDS/listOfInteractions in allDS/allInteractions ===
  // Warning: DS/Interactions are not copied but pointers are inserted into the corresponding set.
  allDS.reset(new DynamicalSystemsSet());
  allInteractions.reset(new InteractionsSet());
  InteractionsIterator itInter;
  DSIterator itDS;
  for (itDS = listOfDS.begin(); itDS != listOfDS.end(); ++itDS)
  {
    allDS->insert(*itDS);
  }

  for (itInter = listOfInteractions.begin(); itInter != listOfInteractions.end(); ++itInter)
  {
    allInteractions->insert(*itInter);
  }

  // === build topology ===
  topology.reset(new Topology(allInteractions));
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(DynamicalSystemsSet& listOfDS, const bool& isBVP):
  BVP(isBVP)
{

  // === Checks that sets are not empty ===
  if (listOfDS.isEmpty())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: constructor(DynamicalSystemsSet, ...): the set of DS is empty.");

  // === "copy" listOfDS/listOfInteractions in allDS/allInteractions ===
  // Warning: DS/Interactions are not copied but pointers are inserted into the corresponding set.
  allDS.reset(new DynamicalSystemsSet());
  allInteractions.reset(new InteractionsSet());
  DSIterator itDS;
  for (itDS = listOfDS.begin(); itDS != listOfDS.end(); ++itDS)
  {
    allDS->insert(*itDS);
  }

  // === build topology ===
  topology.reset(new Topology(allInteractions));
}

// === DynamicalSystems management ===

SP::DynamicalSystem NonSmoothDynamicalSystem::getDynamicalSystemPtrNumber(int nb) const
{
  // if ds number nb is not in the set ...
  assert(allDS->isIn(nb) && "NonSmoothDynamicalSystem::getDynamicalSystemOnNumber(nb), DS number nb is not in the set.");
  return allDS->getPtr(nb);
}

void NonSmoothDynamicalSystem::setDynamicalSystems(const DynamicalSystemsSet& newVect)
{

  // clear old set
  allDS->clear();

  // copy the new one
  for (DSIterator it = newVect.begin(); it != newVect.end(); ++it)
  {
    allDS->insert(*it);
  }

  topology->setUpToDate(false);
}

const bool NonSmoothDynamicalSystem::hasDynamicalSystemNumber(const int& nb) const
{
  return allDS->isIn(nb);
}

const bool NonSmoothDynamicalSystem::hasDynamicalSystem(SP::DynamicalSystem ds) const
{
  return allDS->isIn(ds);
}

// === Interactions management ===

SP::Interaction NonSmoothDynamicalSystem::getInteractionPtr(const int& nb) const
{
  // Mind that Interactions are sorted in a growing order according to their id number.
  InteractionsIterator it = allInteractions->begin();
  for (int i = 0; i < nb; ++i)
    it++;

  if (it == allInteractions->end())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getInteractionPtr(nb) : nb is out of range");

  return *it;
}

SP::Interaction NonSmoothDynamicalSystem::getInteractionPtrNumber(const int& nb) const
{
  if (! allInteractions->isIn(nb)) // if Interaction number nb is not in the set ...
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::getInteractionOnNumber(nb), Interaction number nb is not in the set.");

  return allInteractions->getPtr(nb);
}

void NonSmoothDynamicalSystem::setInteractions(const InteractionsSet& newVect)
{
  // clear old set
  allInteractions->clear();

  // copy the new one
  for (InteractionsIterator it = newVect.begin(); it != newVect.end(); ++it)
  {
    allInteractions->insert(*it);
  }

  topology->setUpToDate(false);
}

const bool NonSmoothDynamicalSystem::hasInteractionNumber(const int& nb) const
{
  return allInteractions->isIn(nb);
}

const bool NonSmoothDynamicalSystem::hasInteraction(SP::Interaction inter) const
{
  return allInteractions->isIn(inter);
}

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
  if (nsdsxml)
  {
    nsdsxml->setBVP(BVP);// no need to change the value of BVP, it mustn't change anyway

    DSIterator it;
    for (it = allDS->begin(); it != allDS->end(); ++it)
    {
      if ((*it)->getType() == LNLDS)
        (boost::static_pointer_cast<LagrangianDS>((*it)))->saveDSToXML();
      else if ((*it)->getType() == LLTIDS)
        (boost::static_pointer_cast<LagrangianLinearTIDS>((*it)))->saveDSToXML();
      else if ((*it)->getType() == FOLDS)
        (boost::static_pointer_cast<FirstOrderLinearDS>((*it)))->saveDSToXML();
      else if ((*it)->getType() == FOLTIDS)
        (boost::static_pointer_cast<FirstOrderLinearDS>((*it)))->saveDSToXML();
      else if ((*it)->getType() == FONLDS)
        (*it)->saveDSToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of DynamicalSystem");
    }

    InteractionsIterator it2;
    for (it2 = allInteractions->begin(); it2 != allInteractions->end(); ++it2)
      (*it2)->saveInteractionToXML();
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NonSmoothDynamicalSystemXML object doesn't exists");
}

void NonSmoothDynamicalSystem::display() const
{
  cout << " ===== Non Smooth Dynamical System display ===== " << endl;
  cout << "---> isBVP = " << BVP << endl;
  allDS->display();
  allInteractions->display();
  cout << "===================================================" << endl;
}

void NonSmoothDynamicalSystem::addDynamicalSystemPtr(SP::DynamicalSystem ds)
{
  allDS->insert(ds);
  topology->setUpToDate(false);
}

void NonSmoothDynamicalSystem::addInteractionPtr(SP::Interaction inter)
{
  allInteractions->insert(inter);
  // the topology should be updated
  topology->setUpToDate(false);
}

double NonSmoothDynamicalSystem::nsdsConvergenceIndicator()
{
  // calculate the max value of all DS convergence indicators
  double convergenceIndicator = (* allDS->begin()) -> dsConvergenceIndicator();
  double dsIndic ;
  DSIterator iter;
  for (iter = allDS->begin(); iter != allDS->end(); ++iter)
  {
    dsIndic = (*iter)->dsConvergenceIndicator();
    if (dsIndic > convergenceIndicator) convergenceIndicator = dsIndic;
  }
  return(convergenceIndicator);
}

