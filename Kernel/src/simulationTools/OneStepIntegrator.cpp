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

#include "OneStepIntegrator.h"
#include "OneStepIntegratorXML.h"
#include "Simulation.h"
#include "Model.h"
#include "DynamicalSystem.h"
#include "NonSmoothDynamicalSystem.h"

using namespace std;

OneStepIntegrator::OneStepIntegrator(const string& id, SimulationSPtr newS):
  integratorType(id),
  sizeMem(1), simulationLink(newS)
{
  assert(simulationLink &&
         "OneStepIntegrator:: constructor(Id,simulation) - simulation == NULL");

  OSIDynamicalSystems.reset(new DynamicalSystemsSet());
  OSIInteractions.reset(new InteractionsSet());
}

// --- Xml constructor ---
OneStepIntegrator::OneStepIntegrator(const string& id, OneStepIntegratorXMLSPtr osixml, SimulationSPtr newS):
  integratorType(id), sizeMem(1), simulationLink(newS), integratorXml(osixml)
{
  if (!integratorXml)
    RuntimeException::selfThrow("OneStepIntegrator::xml constructor - OneStepIntegratorXML object == NULL");

  if (!simulationLink)
    RuntimeException::selfThrow("OneStepIntegrator::xml constructor - SimulationLink == NULL");

  simulationLink->addOneStepIntegratorPtr(shared_from_this());

  OSIDynamicalSystems.reset(new DynamicalSystemsSet());
  OSIInteractions.reset(new InteractionsSet());

  // get a link to the NonSmoothDynamicalSystem.
  SP::NonSmoothDynamicalSystem nsds = simulationLink->getModelPtr()->getNonSmoothDynamicalSystemPtr();
  // load DS list if present
  if (osixml->hasDSList())
  {
    if (osixml->hasAllDS()) // if flag all=true is present -> get all ds from the nsds
    {
      DSIterator it;
      for (it = nsds->dynamicalSystemsBegin(); it != nsds->dynamicalSystemsEnd(); ++it)
        OSIDynamicalSystems->insert(*it);
    }
    else
    {
      // get list of ds numbers implicate in the OSI
      vector<int> dsNumbers;
      osixml->getDSNumbers(dsNumbers);
      // get corresponding ds and insert them into the set.
      vector<int>::iterator it;
      for (it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
        OSIDynamicalSystems->insert(nsds->getDynamicalSystemPtrNumber(*it));
    }
  }

  // load interactions list if present
  if (osixml->hasInteractionsList())
  {
    if (osixml->hasAllInteractions()) // if flag all=true is present -> get all interactions from the nsds
    {
      InteractionsIterator it;
      for (it = nsds->interactionsBegin(); it != nsds->interactionsEnd(); ++it)
        OSIInteractions->insert(*it);
    }
    else
    {
      // get list of interactions numbers implicate in the OSI
      vector<int> interactionsNumbers;
      osixml->getInteractionsNumbers(interactionsNumbers);
      // get corresponding interactions and insert them into the set.
      vector<int>::iterator it;
      for (it = interactionsNumbers.begin(); it != interactionsNumbers.end(); ++it)
        OSIInteractions->insert(nsds->getInteractionPtrNumber(*it));
    }
  }
}

// --- Constructors from a minimum set of data ---
OneStepIntegrator::OneStepIntegrator(const string& id, const DynamicalSystemsSet& listOfDs, SimulationSPtr newS):
  integratorType(id), sizeMem(1), simulationLink(newS)
{
  if (!simulationLink)
    RuntimeException::selfThrow("OneStepIntegrator:: constructor(Id,listDS,simulation) - simulation == NULL");

  OSIDynamicalSystems.reset(new DynamicalSystemsSet());
  OSIInteractions.reset(new InteractionsSet());
  setDynamicalSystems(listOfDs);

}

// --- Destructor ---
OneStepIntegrator::~OneStepIntegrator()
{
}


void OneStepIntegrator::setDynamicalSystems(const DynamicalSystemsSet& newSet)
{
  // Warning: pointers links between ds of newSet and OSIDynamicalSystems.
  DSIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
    OSIDynamicalSystems->insert(*it);
}

void OneStepIntegrator::setInteractions(const InteractionsSet& newSet)
{
  // Warning: pointers links between ds of newSet and OSIDynamicalSystems.
  InteractionsIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
    OSIInteractions->insert(*it);
}

SiconosVectorSPtr OneStepIntegrator::getWorkX(SP::DynamicalSystem ds)
{
  if (workX.find(ds) == workX.end())
    RuntimeException::selfThrow("OneStepIntegrator::getWorkX(ds): this vector does not exists for ds.");
  return workX[ds];
}

void OneStepIntegrator::initialize()
{

  simulationLink->addOneStepIntegratorPtr(shared_from_this());

  double t0 = simulationLink->getModelPtr()->getT0();
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // Initialization of the dynamical systems belonging to this OSI
    (*it)->initialize(simulationLink->getType(), t0, sizeMem);

    // Register this DS and the OSI into OSIMap of the Simulation
    simulationLink->addInOSIMap(*it, shared_from_this());
  }
}

void OneStepIntegrator::saveInMemory()
{
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    (*it)->swapInMemory();
}

double OneStepIntegrator::computeResidu()
{
  RuntimeException::selfThrow("OneStepIntegrator::computeResidu not implemented for integrator of type " + integratorType);
  return 0.0;
}

void OneStepIntegrator::computeFreeState()
{
  RuntimeException::selfThrow("OneStepIntegrator::computeFreeState not implemented for integrator of type " + integratorType);
}

void OneStepIntegrator::resetNonSmoothPart()
{
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    (*it)->resetNonSmoothPart();
}

void OneStepIntegrator::display()
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

