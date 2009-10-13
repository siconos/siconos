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

OneStepIntegrator::OneStepIntegrator(const OSI::TYPES& id):
  integratorType(id), sizeMem(1)
{
  OSIDynamicalSystems.reset(new DynamicalSystemsSet());
  OSIInteractions.reset(new InteractionsSet());
}

// --- Xml constructor ---
OneStepIntegrator::OneStepIntegrator(const OSI::TYPES& id, SP::OneStepIntegratorXML osixml,
                                     SP::DynamicalSystemsSet dsList, SP::InteractionsSet interactionsList):
  integratorType(id), sizeMem(1), integratorXml(osixml)
{
  if (!integratorXml)
    RuntimeException::selfThrow("OneStepIntegrator::xml constructor - OneStepIntegratorXML object == NULL");

  OSIDynamicalSystems.reset(new DynamicalSystemsSet());
  OSIInteractions.reset(new InteractionsSet());

  // load DS list if present
  if (osixml->hasDSList())
  {
    assert(dsList && "OneStepIntegrator xml constructor: empty ds list from NSDS.");
    if (osixml->hasAllDS()) // if flag all=true is present -> get all ds from the nsds
      OSIDynamicalSystems->insert(dsList->begin(), dsList->end());

    else
    {
      // get list of ds numbers implicate in the OSI
      vector<int> dsNumbers;
      osixml->getDSNumbers(dsNumbers);
      // get corresponding ds and insert them into the set.
      for (vector<int>::iterator it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
        OSIDynamicalSystems->insert(dsList->getPtr(*it));
    }
  }

  // load interactions list if present
  if (osixml->hasInteractionsList())
  {
    assert(interactionsList && "OneStepIntegrator xml constructor: empty interaction list from NSDS.");
    if (osixml->hasAllInteractions()) // if flag all=true is present -> get all interactions
      OSIInteractions->insert(interactionsList->begin(), interactionsList->end());

    else
    {
      // get list of interactions numbers implicate in the OSI
      vector<int> interactionsNumbers;
      osixml->getInteractionsNumbers(interactionsNumbers);
      // get corresponding interactions and insert them into the set.
      for (vector<int>::iterator it = interactionsNumbers.begin(); it != interactionsNumbers.end(); ++it)
        OSIInteractions->insert(interactionsList->getPtr(*it));
    }
  }
}

// --- Constructors from a minimum set of data ---
OneStepIntegrator::OneStepIntegrator(const OSI::TYPES& id, const DynamicalSystemsSet& listOfDs):
  integratorType(id), sizeMem(1)
{
  OSIDynamicalSystems.reset(new DynamicalSystemsSet());
  OSIInteractions.reset(new InteractionsSet());
  setDynamicalSystems(listOfDs);
}

void OneStepIntegrator::setDynamicalSystems(const DynamicalSystemsSet& newSet)
{
  OSIDynamicalSystems->insert(newSet.begin(), newSet.end());
}

void OneStepIntegrator::setInteractions(const InteractionsSet& newSet)
{
  OSIInteractions->insert(newSet.begin(), newSet.end());
}

// SP::SiconosVector OneStepIntegrator::getWorkX(SP::DynamicalSystem ds)
// {
//   assert(workX.find(ds)!=workX.end()&&"OneStepIntegrator::getWorkX(ds): this vector does not exists for ds.");
//   return workX[ds];
// }

void OneStepIntegrator::initialize(SP::Simulation sim)
{
  // Connection to the simulation owner of this OSI
  assert(sim && "OneStepIntegrator::initialize(sim) error: sim is null.");
  simulationLink = sim;

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
  for_each(OSIDynamicalSystems->begin(), OSIDynamicalSystems->end(), boost::bind(&DynamicalSystem::swapInMemory, _1));
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
  for_each(OSIDynamicalSystems->begin(), OSIDynamicalSystems->end(), boost::bind(&DynamicalSystem::resetNonSmoothPart, _1));
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

