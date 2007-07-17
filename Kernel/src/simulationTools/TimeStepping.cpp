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

#include "TimeStepping.h"
#include "SimulationXML.h"
#include "OneStepNSProblemXML.h"
#include "Topology.h"
#include "LCP.h"
#include "FrictionContact2D.h"
#include "FrictionContact3D.h"
#include "Model.h"
#include "TimeDiscretisation.h"
#include "NonSmoothDynamicalSystem.h"
#include "UnitaryRelation.h"
#include "OneStepIntegrator.h"
#include "Interaction.h"
#include "EventsManager.h"

using namespace std;

// --- Default constructor ---
TimeStepping::TimeStepping(): Simulation("TimeStepping")
{}

TimeStepping::TimeStepping(TimeDiscretisation * td): Simulation(td, "TimeStepping")
{}

// --- XML constructor ---
TimeStepping::TimeStepping(SimulationXML* strxml, Model *newModel): Simulation(strxml, newModel, "TimeStepping")
{
  // === One Step NS Problem ===
  // For time stepping, only one non smooth problem is built.
  if (simulationxml->hasOneStepNSProblemXML())
  {
    // OneStepNSProblem - Memory allocation/construction
    string type = simulationxml->getOneStepNSProblemXMLPtr()->getNSProblemType();
    if (type == LCP_TAG)  // LCP
    {
      allNSProblems["timeStepping"] = new LCP(simulationxml->getOneStepNSProblemXMLPtr(), this);
      isNSProblemAllocatedIn[ allNSProblems["timeStepping"] ] = true;
    }
    else if (type == FrictionContact2D_TAG) // Friction 2D
    {
      allNSProblems["timeStepping"] = new FrictionContact2D(simulationxml->getOneStepNSProblemXMLPtr(), this);
      isNSProblemAllocatedIn[ allNSProblems["timeStepping"] ] = true;
    }
    else if (type == FrictionContact3D_TAG) // Friction 3D
    {
      allNSProblems["timeStepping"] = new FrictionContact3D(simulationxml->getOneStepNSProblemXMLPtr(), this);
      isNSProblemAllocatedIn[ allNSProblems["timeStepping"] ] = true;
    }
    else RuntimeException::selfThrow("TimeStepping::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");

    allNSProblems["timeStepping"]->setId("timeStepping");

    // Add QP and Relay cases when these classes will be fully implemented.
  }
}

// --- Destructor ---
TimeStepping::~TimeStepping()
{}

void TimeStepping::updateIndexSet(unsigned int i)
{
  if (i > indexSets.size())
    RuntimeException::selfThrow("TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");

  if (i == 0) // IndexSets[0] must not be updated by this function.
    RuntimeException::selfThrow("TimeStepping::updateIndexSet(i=0), indexSets[0] can not be updated.");

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and update the indexSet[i]
  UnitaryRelationIterator it, itForFind;

  double y;

  if (i == 1) // special case for Moreau time-stepping
  {
    double yp;
    double yDot;
    for (it = (indexSets[0])->begin(); it != (indexSets[0])->end(); ++it)
    {
      double h = getTimeStep();
      // checks if current Unitary Relation (ie *it) is already in indexSets[1]
      // (if not itForFind will be equal to indexSets.end())
      itForFind = (indexSets[1])->find(*it);
      y = (*it)->getYRef(0);
      yDot = (*it)->getYRef(1);
      yp = y + 0.5 * h * yDot;

      // if yp <=0, then the unitary relation is added in indexSets[1] (if it was not already there)
      // else if yp > 0 and if the unitary relation was in the set, it is removed.
      if (yp <= 0 && itForFind == (indexSets[1])->end())
        (indexSets[1])->insert(*it);

      else if (yp > 0 && itForFind != (indexSets[1])->end())
        (indexSets[1])->erase(*it);
    }
  }
  else
  {
    for (it = indexSets[i - 1]->begin(); it != indexSets[i - 1]->end(); ++it)
    {
      // check if current Unitary Relation (ie *it) is in indexSets[i]
      // (if not itForFind will be equal to indexSets.end())
      itForFind = indexSets[i]->find(*it);

      // Get y[i-1] double value
      y = (*it)->getYRef(i - 1);

      // if y[i-1] <=0, then the unitary relation is added in indexSets[i] (if it was not already there)
      // else if y[i-1] > 0 and if the unitary relation was in the set, it is removed.
      if (y <= 0 && itForFind == indexSets[i]->end())
        indexSets[i]->insert(*it);
      else if (y > 0 && itForFind != indexSets[i]->end())
        indexSets[i]->erase(*it);
    }
  }
}

void TimeStepping::addOneStepNSProblemPtr(OneStepNSProblem* osns)
{
  // A the time, a time stepping simulation can only have one non smooth problem.
  if (!allNSProblems.empty())
    RuntimeException::selfThrow("TimeStepping, addOneStepNSProblemPtr - A non smooth problem already exist. You can not have more than one.");

  string name = "timeStepping"; // osns->getId();
  osns->setId(name);
  allNSProblems[name] = osns;
  isNSProblemAllocatedIn[osns] = false;
}

void TimeStepping::initOSNS()
{
  if (!allNSProblems.empty()) // ie if some Interactions have been declared and a Non smooth problem built.
  {
    if (allNSProblems.size() > 1)
      RuntimeException::selfThrow("TimeStepping::initialize, at the time, a time stepping simulation can not have more than one non smooth problem.");

    // At the time, we consider that for all systems, levelMin is equal to the minimum value of the relative degree - 1
    // except for degree 0 case where we keep 0.

    levelMin = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getMinRelativeDegree();

    if (levelMin != 0)
      levelMin--;

    // === update all index sets ===
    updateIndexSets();

    // initialization of  OneStepNonSmoothProblem
    OSNSIterator itOsns;
    for (itOsns = allNSProblems.begin(); itOsns != allNSProblems.end(); ++itOsns)
    {
      (itOsns->second)->setLevels(levelMin, levelMax);
      (itOsns->second)->initialize();
    }
  }
}

void TimeStepping::initLevelMax()
{
  levelMax = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getMaxRelativeDegree();
  // Interactions initialization (here, since level depends on the type of simulation)
  // level corresponds to the number of Y and Lambda derivatives computed.

  if (levelMax != 0) // level max is equal to relative degree-1. But for relative degree 0 case, we keep 0 value for levelMax
    levelMax--;
}

void TimeStepping::nextStep()
{
  eventsManager->processEvents();
}


void TimeStepping::update(unsigned int levelInput)
{
  // compute input (lambda -> r)
  updateInput(levelInput);

  // compute state for each dynamical system

  OSIIterator itOSI;
  for (itOSI = allOSI.begin(); itOSI != allOSI.end() ; ++itOSI)
    (*itOSI)->updateState(levelInput);

  if (!allNSProblems.empty() && levelMin > 0)
  {
    updateOutput(0, levelMax);
  }
}

void TimeStepping::computeFreeState()
{
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->computeFreeState();
}

// compute simulation between current and next event.
// Initial DS/interaction state is given by memory vectors
// and final state is the one saved in DS/Interaction at the end of this function
void TimeStepping::computeOneStep()
{
  advanceToEvent();
}

void TimeStepping::advanceToEvent()
{
  // solve ...

  computeFreeState();
  if (!allNSProblems.empty())
  {
    updateIndexSets();
    computeOneStepNSProblem("timeStepping");
  }
  // update
  update(levelMin);
}

void TimeStepping::newtonSolve(double criterion, unsigned int maxStep)
{
  bool isNewtonConverge = false;
  unsigned int nbNewtonStep = 0; // number of Newton iterations

  while ((!isNewtonConverge) && (nbNewtonStep <= maxStep))
  {
    nbNewtonStep++;
    advanceToEvent();
    // Process all events simultaneous to nextEvent.
    //  eventsManager->process();
    isNewtonConverge = newtonCheckConvergence(criterion);
  }

  // Process NextEvent (Save OSI (DS) and OSNS (Interactions) states into Memory vectors ...)
  eventsManager->processEvents();

  if (!isNewtonConverge)
    cout << "Newton process stopped: reach max step number" << endl ;
}

bool TimeStepping::newtonCheckConvergence(double criterion)
{
  bool checkConvergence = false;
  // get the nsds indicator of convergence
  double nsdsConverge = model-> getNonSmoothDynamicalSystemPtr()->nsdsConvergenceIndicator();
  if (nsdsConverge < criterion) checkConvergence = true ;
  return(checkConvergence);
}

void TimeStepping::run(const std::string& opt, double criterion, unsigned int maxIter)
{
  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of " << simulationType << " simulation - This may take a while ... ====" << endl;
  while (eventsManager->hasNextEvent())
  {
    if (opt == "linear")
    {
      advanceToEvent();
      eventsManager->processEvents();
    }
    else if (opt == "Newton")
      newtonSolve(criterion, maxIter);
    else
      RuntimeException::selfThrow("TimeStepping::run(opt) failed. Unknow simulation option: " + opt);
    count++;
  }
  cout << "===== End of " << simulationType << "simulation. " << count << " events have been processed. ==== " << endl;
}

TimeStepping* TimeStepping::convert(Simulation *str)
{
  TimeStepping* ts = dynamic_cast<TimeStepping*>(str);
  return ts;
}

