/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#include "OneStepIntegrator.h"
using namespace std;

// RuntimeCmp object used to sort DynamicalSystems in the ds set
// !!! \todo Find a way to avoid global variable ... !!!
RuntimeCmp<DynamicalSystem> compareDS(&DynamicalSystem::getNumberForSorting);

//-- Default constructor --
OneStepIntegrator::OneStepIntegrator(const string& id, Strategy* newS): integratorType(id), dsList(compareDS), sizeMem(1), strategyLink(newS), integratorXml(NULL)
{
  strategyLink->addOneStepIntegrator(this);
}

// --- Xml constructor ---
OneStepIntegrator::OneStepIntegrator(const string& id, OneStepIntegratorXML* osixml, Strategy* newS): integratorType(id), dsList(compareDS), sizeMem(1), strategyLink(newS), integratorXml(osixml)
{
  if (integratorXml == NULL)
    RuntimeException::selfThrow("OneStepIntegrator::xml constructor - OneStepIntegratorXML object == NULL");

  if (strategyLink == NULL)
    RuntimeException::selfThrow("OneStepIntegrator::xml constructor - StrategyLink == NULL");

  strategyLink->addOneStepIntegrator(this);

  // get a link to the NonSmoothDynamicalSystem.
  NonSmoothDynamicalSystem * nsds = strategyLink->getModelPtr()->getNonSmoothDynamicalSystemPtr();
  // load DS list if present
  if (osixml->hasDSList())
  {
    if (osixml->hasAllDS()) // if flag all=true is present -> get all ds from the nsds
    {
      // In nsds DS are saved in a vector.
      // Future version: saved them in a set? And then just call:
      // dsList = nsds->getDynamicalSystems();
      vector<DynamicalSystem*> tmpDS;
      vector<DynamicalSystem*>::iterator it;
      for (it = tmpDS.begin(); it != tmpDS.end(); ++it)
        dsList.insert(*it);
    }
    else
    {
      // get list of ds numbers implicate in the OSI
      vector<int> dsNumbers;
      osixml->getDSNumbers(dsNumbers);
      // get corresponding ds and insert them into the set.
      vector<int>::iterator it;
      for (it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
        dsList.insert(nsds->getDynamicalSystemPtrNumber(*it));
    }
  }

  // load interactions list if present
  if (osixml->hasInteractionsList())
  {
    if (osixml->hasAllInteractions()) // if flag all=true is present -> get all interactions from the nsds
    {
      // In nsds interactions are saved in a vector.
      // Future version: saved them in a set? And then just call:
      // interactionsList = nsds -> getInteractions();
      vector<Interaction*> tmpI;
      vector<Interaction*>::iterator it;
      for (it = tmpI.begin(); it != tmpI.end(); ++it)
        interactionsList.insert(*it);
    }
    else
    {
      // get list of interactions numbers implicate in the OSI
      vector<int> interactionsNumbers;
      osixml->getInteractionsNumbers(interactionsNumbers);
      // get corresponding interactions and insert them into the set.
      vector<int>::iterator it;
      for (it = interactionsNumbers.begin(); it != interactionsNumbers.end(); ++it)
        interactionsList.insert(nsds->getInteractionPtrNumber(*it));
    }
  }
}

// --- Constructor from a minimum set of data ---
OneStepIntegrator::OneStepIntegrator(const string& id, const dsSet& listOfDs, Strategy* newS):
  integratorType(id), dsList(compareDS), sizeMem(1), strategyLink(newS), integratorXml(NULL)
{
  strategyLink->addOneStepIntegrator(this);
}

// --- Destructor ---
OneStepIntegrator::~OneStepIntegrator()
{
  dsList.clear();
  interactionsList.clear();
  integratorXml = NULL;
}

void OneStepIntegrator::setDynamicalSystemsList(const dsSet& newSet)
{
  // Warning: pointers links between ds of newSet and dsList.
  dsIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
    dsList.insert(*it);
}

void OneStepIntegrator::setInteractionsList(const interactionSet& newSet)
{
  // Warning: pointers links between ds of newSet and dsList.
  interactionIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
    interactionsList.insert(*it);
}

void OneStepIntegrator::initialize()
{
  double t0 = strategyLink->getTimeDiscretisationPtr()->getT0();
  dsIterator it;
  for (it = dsList.begin(); it != dsList.end(); ++it)
    (*it)->initialize(t0, sizeMem);
}

void OneStepIntegrator::nextStep()
{
  dsIterator it;
  for (it = dsList.begin(); it != dsList.end(); ++it)
  {
    (*it)->swapInMemory();
    (*it)->getRPtr()->zero();
  }
}

void OneStepIntegrator::display() const
{
  cout << "==== OneStepIntegrator display =====" << endl;
  cout << "| integratorType : " << integratorType << endl;
  cout << "| sizeMem: " << sizeMem << endl;
  cout << "====================================" << endl;
}

void OneStepIntegrator::saveIntegratorToXML()
{
  //   if(integratorXml != 0)
  //     {
  //       vector<int> dsConcerned;
  //       dsConcerned.push_back(ds->getNumber());
  //       integratorXml->setDSConcerned( &dsConcerned );
  //     }
  //   else
  //RuntimeException::selfThrow("OneStepIntegrator::saveIntegratorToXML - OneStepIntegratorXML object = NULL");
  RuntimeException::selfThrow("OneStepIntegrator::saveIntegratorToXML - Not yet implemented.");
}

