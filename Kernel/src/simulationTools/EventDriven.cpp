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

#include "EventDriven.h"
#include "Lsodar.h"

using namespace std;

// --- Default constructor ---
EventDriven::EventDriven(Model* newModel): Simulation(newModel, "EventDriven")
{}

// --- XML constructor ---
EventDriven::EventDriven(SimulationXML* strxml, Model *newModel): Simulation(strxml, newModel, "EventDriven")
{
  // === One Step NS Problem ===
  // We read data in the xml output (mainly Interactions concerned and solver) and assign them to
  // both one step ns problems ("acceleration" and "impact").
  // At the time, only LCP is available for Event Driven.

  if (simulationxml->hasOneStepNSProblemXML())
  {
    // OneStepNSProblem - Memory allocation/construction
    string type = simulationxml->getOneStepNSProblemXMLPtr()->getNSProblemType();
    if (type == LCP_TAG)  // LCP
    {
      allNSProblems["acceleration"] = new LCP(simulationxml->getOneStepNSProblemXMLPtr(), this);
      isNSProblemAllocatedIn[ allNSProblems["acceleration"] ] = true;
      allNSProblems["impact"] = new LCP(simulationxml->getOneStepNSProblemXMLPtr(), this);
      isNSProblemAllocatedIn[ allNSProblems["impact"] ] = true;
    }
    else
      RuntimeException::selfThrow("EventDriven::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");

    allNSProblems["acceleration"]->setId("acceleration");
    allNSProblems["impact"]->setId("impact");
  }
}

// --- Destructor ---
EventDriven::~EventDriven()
{
  if (eventsManager != NULL) delete eventsManager;
  eventsManager = NULL;
}

void EventDriven::setEventsManagerPtr(EventsManager*)
{
  // TODO IF NECESSARY?
}

void EventDriven::updateIndexSet(const unsigned int i)
{
  if (i > indexSets.size())
    RuntimeException::selfThrow("Topology::updateIndexSet(i), indexSets[i] does not exist.");

  if (i == 0) // IndexSets[0] must not be updated by this function.
    RuntimeException::selfThrow("Topology::updateIndexSet(i=0), indexSets[0] can not be updated.");

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and update the indexSet[i]
  UnitaryRelationIterator it, itForFind;

  double y;
  for (it = indexSets[i - 1].begin(); it != indexSets[i - 1].end(); ++it)
  {
    // check if current Unitary Relation (ie *it) is in indexSets[i]
    // (if not itForFind will be equal to indexSets.end())
    itForFind = indexSets[i].find(*it);

    // Get y[i-1] double value
    y = (*it)->getYRef(i - 1);

    // if y[i-1] <=0, then the unitary relation is added in indexSets[i] (if it was not already there)
    // else if y[i-1] > 0 and if the unitary relation was in the set, it is removed.
    if (y <= 0 && itForFind == indexSets[i].end())
      indexSets[i].insert(*it);
    else if (y > 0 && itForFind != indexSets[i].end())
      indexSets[i].erase(*it);
  }
}

void EventDriven::updateIndexSetsWithDoubleCondition()
{

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and update the indexSet[i]
  UnitaryRelationIterator it, itForFind;

  for (it = indexSets[2].begin(); it != indexSets[2].end(); ++it)
  {
    double gamma = (*it)->getYRef(2);
    double F     = (*it)->getLambdaRef(2);

    if (gamma > 0 && F < TOLERANCE)
      indexSets[2].erase(*it);
    else if (gamma < TOLERANCE && F < TOLERANCE) // undetermined case
      RuntimeException::selfThrow("Topology::updateIndexSetsWithDoubleCondition(), undetermined case.");
  }
}

void EventDriven::initialize()
{
  Simulation::initialize();

  // === Events manager creation and initialization ===
  eventsManager = new EventsManager(DEFAULT_TICK, this); //
  eventsManager->initialize();

  // === OneStepNSProblem initialization. ===
  // First check that there are 2 osns: one "impact" and one "acceleration"
  if (allNSProblems.size() != 2)
    RuntimeException::selfThrow("EventDriven::initialize, an EventDriven simulation must have two non smooth problem. Here, there are " + allNSProblems.size());

  if (allNSProblems.find("impact") == allNSProblems.end()) // ie if the impact problem does not exist
    RuntimeException::selfThrow("EventDriven::initialize, an EventDriven simulation must have an 'impact' non smooth problem.");
  if (allNSProblems.find("acceleration") == allNSProblems.end()) // ie if the impact problem does not exist
    RuntimeException::selfThrow("EventDriven::initialize, an EventDriven simulation must have an 'acceleration' non smooth problem.");

  // WARNING: only for Lagrangian systems - To be reviewed for other ones.
  allNSProblems["impact"]->setLevels(1, 1);
  allNSProblems["impact"]->initialize();
  allNSProblems["acceleration"]->setLevels(2, 2);
  allNSProblems["acceleration"]->initialize();

  // === update all index sets ===
  updateIndexSets();
}

void EventDriven::computeF(OneStepIntegrator* osi, integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{
  // Check osi type: only lsodar is allowed.
  if (osi->getType() != "lsodar")
    RuntimeException::selfThrow("EventDriven::computeF(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  Lsodar * lsodar = static_cast<Lsodar*>(osi);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar->fillXWork(sizeOfX, x);

  double t = *time;

  // update the DS of the OSI.
  lsodar->updateState(t);

  // solve a LCP "acceleration"

  allNSProblems["acceleration"]->compute(t);
  allNSProblems["acceleration"]->updateOutput();
  allNSProblems["acceleration"]->updateInput();

  // Compute the right-hand side ( xdot = f + Tu in DS) for all the ds

  lsodar->computeRhs(t);

  DynamicalSystemsSet dsOfTheOsi = lsodar->getDynamicalSystems();

  DSIterator it;
  unsigned int i = 0;
  for (it = dsOfTheOsi.begin(); it != dsOfTheOsi.end(); ++it)
  {
    SiconosVector * xtmp2 = (*it)->getRhsPtr(); // Pointer link !
    for (unsigned int j = 0 ; j < (*it)->getDim() ; ++j)
      xdot[i++] = (*xtmp2)(j);
  }
}

void EventDriven::computeJacobianF(OneStepIntegrator* osi, integer *sizeOfX, doublereal *time, doublereal *x,  doublereal *jacob)
{

  if (osi->getType() != "lsodar")
    RuntimeException::selfThrow("EventDriven::computeF(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  Lsodar * lsodar = static_cast<Lsodar*>(osi);

  //   // Remark A: according to DLSODAR doc, each call to jacobian is preceded by a call to f with the same
  //   // arguments NEQ, T, and Y.  Thus to gain some efficiency, intermediate quantities shared by both calculations may be
  //   // saved in class members?
  //   cout <<"in jaco f: " <<  endl;

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  // fillXWork(x); // -> copy // Maybe this step is not necessary? because of remark A above

  // Compute the jacobian of the vector field according to x for the current ds
  double t = *time;
  lsodar->computeJacobianRhs(t);

  //   // Save jacobianX values from dynamical system into current jacob (in-out parameter)
  DynamicalSystemsSet dsOfTheOsi = lsodar->getDynamicalSystems();

  DSIterator it;
  unsigned int i = 0;
  for (it = dsOfTheOsi.begin(); it != dsOfTheOsi.end(); ++it)
  {
    SiconosMatrix * jacotmp = (*it)->getJacobianXFPtr(); // Pointer link !
    for (unsigned int j = 0 ; j < (*it)->getDim() ; ++j)
    {
      for (unsigned k = 0 ; k < (*it)->getDim() ; ++k)
        jacob[i++] = (*jacotmp)(k, j);
    }
  }
}

void EventDriven::computeG(OneStepIntegrator* osi)
{
  RuntimeException::selfThrow("EventDriven::computeG(osi, ...), not yet implemented");
}

// Run the whole simulation
void EventDriven::run()
{

  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of Event Driven simulation - This may take a while ... ====" << endl;
  while (eventsManager->hasNextEvent())
  {
    // Integrate system between "current" and "next" event of events manager.
    advanceToEvent();
    // update events
    eventsManager->processEvents();
    update();
    count++;
  }
  cout << "===== End of Event Driven simulation. " << count << " events have been processed. ==== " << endl;
}

void EventDriven::computeOneStep()
{
  // Integrate system between "current" and "next" event of events manager.
  advanceToEvent();
  // update events
  eventsManager->processEvents();
  update();
}

void EventDriven::update()
{
  // compute input (lambda -> r)
  OSNSIterator itOsns;
  for (itOsns = allNSProblems.begin(); itOsns != allNSProblems.end(); ++itOsns)
  {
    (itOsns->second)->updateInput();
    (itOsns->second)->updateOutput();
    (itOsns->second)->nextStep();
  }

  updateIndexSets();

  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->nextStep();
}

void EventDriven::advanceToEvent()
{

  // WARNING: this is supposed to work for only one OSI, including all the DS.
  // To be reviewed for multiple OSI case (if it has sense?).

  // ---> Step 1: integrate the smooth dynamics from current event to next event;
  // Current event = last accessed event.
  // Next event = next time step or first root of the 'g' function found by integrator (Lsodar)
  double tinit = eventsManager->getCurrentTime();
  double tend =  eventsManager->getNextTime();
  double tout = tend;
  bool isNewEventOccur = false;  // set to true if a new event occur during integration
  // call integrate method for each OSI, between tinit and tend.
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
  {
    bool iout = false;
    (*it)->integrate(tinit, tend, tout, iout); // integrate must return a flag telling if tend has been reached or not.
    // If not, tout is the real reached time.
    if (!iout)
    {
      isNewEventOccur = true;
      // Add an event into the events manager list
      bool isScheduleOk = eventsManager->scheduleEvent("NonSmoothEvent", tout);
      if (!isScheduleOk) cout << " EventDriven advanceToEvent warning: try to add an already existing event" << endl;
    }
  }

  // ---> Step 2: update Index sets according to temporary values obtained at previous step.

  updateIndexSets();  // This requires that y[i] values have been well computed and saved in Interactions.

  // ---> Step 3: solve impact LCP if IndexSet[1]\IndexSet[2] is not empty.

  if (!(indexSets[1] - indexSets[2]).isEmpty())
  {
    // solve the LCP-impact => y[1],lambda[1]
    computeOneStepNSProblem("impact"); // solveLCPImpact();
    // update indexSets
    updateIndexSet(1);
    updateIndexSet(2);

    // check that IndexSet[1]-IndexSet[2] is now empty
    if (!(indexSets[1] - indexSets[2]).isEmpty())
      RuntimeException::selfThrow("EventDriven advanceToEvent, error after impact-LCP solving.");
  }

  if (!((indexSets[2]).isEmpty()))
  {
    // solve LCP-acceleration
    computeOneStepNSProblem("acceleration"); //solveLCPAcceleration();
    // for all index in IndexSets[2], update the index set according to y[2] and/or lambda[2] sign.
    updateIndexSetsWithDoubleCondition();
  }

  model->setCurrentT(tout);

}

EventDriven* EventDriven::convert(Simulation *str)
{
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}
