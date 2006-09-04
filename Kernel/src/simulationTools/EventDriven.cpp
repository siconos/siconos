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
#include "LagrangianLinearTIDS.h"
#include "Lsodar.h"

using namespace std;

// --- Default constructor ---
EventDriven::EventDriven(Model* newModel): Simulation(newModel, "EventDriven"), istate(1), tinit(0), tend(0)
{}

// --- XML constructor ---
EventDriven::EventDriven(SimulationXML* strxml, Model *newModel): Simulation(strxml, newModel, "EventDriven"), istate(1), tinit(0), tend(0)
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

  double borneInf;
  if (i == 1)
    borneInf = -1000;
  else
    borneInf = -TOLERANCE;

  double y;
  for (it = indexSets[i - 1].begin(); it != indexSets[i - 1].end(); ++it)
  {
    // check if current Unitary Relation (ie *it) is in indexSets[i]
    // (if not itForFind will be equal to indexSets.end())
    itForFind = indexSets[i].find(*it);

    // Get y[i-1] double value
    y = (*it)->getYRef(i - 1);

    //       // if y[i-1] <=0, then the unitary relation is added in indexSets[i] (if it was not already there)
    //       // else if y[i-1] > 0 and if the unitary relation was in the set, it is removed.
    //       if( y <= 0 && itForFind==indexSets[i].end())
    //  indexSets[i].insert(*it);
    //       else if( y > 0 && itForFind!=indexSets[i].end())
    //  indexSets[i].erase(*it);

    if ((borneInf < y && y < TOLERANCE) &&  itForFind == indexSets[i].end())
      indexSets[i].insert(*it);
    else if ((-TOLERANCE > y || y > TOLERANCE) &&  itForFind != indexSets[i].end())
      indexSets[i].erase(*it);

    // Note if *it is already in the set, insert(*it) is supposed to have no effect. Same for erase when (*it) is not in the set.
    // => to be reviewed in UnitaryRelationsSet

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
      RuntimeException::selfThrow("EventDriven::updateIndexSetsWithDoubleCondition(), undetermined case.");
  }
}

void EventDriven::initialize()
{
  Simulation::initialize();

  InteractionsSet allInteractions = model->getNonSmoothDynamicalSystemPtr()->getInteractions();

  if (!allInteractions.isEmpty()) // ie if some Interactions have been declared
  {
    levelMax = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getMaxRelativeDegree();
    // Interactions initialization (here, since level depends on the type of simulation)
    // level corresponds to the number of Y and Lambda derivatives computed.

    if (levelMax == 0)
      levelMax++;

    InteractionsIterator it;
    for (it = allInteractions.begin(); it != allInteractions.end(); ++it)
      (*it)->initializeMemory(levelMax + 1);

    indexSets.resize(levelMax + 1);
    indexSets[0] = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getIndexSet0();
  }

  // initialization of the OneStepIntegrators
  OSIIterator itOsi;
  for (itOsi = allOSI.begin(); itOsi != allOSI.end(); ++itOsi)
    (*itOsi)->initialize();

  // === Events manager creation and initialization ===
  eventsManager = new EventsManager(DEFAULT_TICK, this); //
  eventsManager->initialize();
  tinit = eventsManager->getCurrentTime();
  tend =  eventsManager->getNextTime();

  if (!allNSProblems.empty()) // ie if some Interactions have been declared and a Non smooth problem built.
  {
    // === OneStepNSProblem initialization. ===
    // First check that there are 2 osns: one "impact" and one "acceleration"
    if (allNSProblems.size() != 2)
      RuntimeException::selfThrow(" EventDriven::initialize, \n an EventDriven simulation must have two non smooth problem.\n Here, there are " + allNSProblems.size());

    if (allNSProblems.find("impact") == allNSProblems.end()) // ie if the impact problem does not exist
      RuntimeException::selfThrow("EventDriven::initialize, an EventDriven simulation must have an 'impact' non smooth problem.");
    if (allNSProblems.find("acceleration") == allNSProblems.end()) // ie if the impact problem does not exist
      RuntimeException::selfThrow("EventDriven::initialize, an EventDriven simulation must have an 'acceleration' non smooth problem.");

    // At the time, we consider that for all systems, levelMin is equal to the minimum value of the relative degree
    levelMin = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getMinRelativeDegree();
    if (levelMin == 0)
      levelMin++;

    updateInput(levelMin);
    updateOutput(0, levelMax);

    // WARNING: only for Lagrangian systems - To be reviewed for other ones.
    allNSProblems["impact"]->setLevels(levelMin - 1, levelMax - 1);
    allNSProblems["impact"]->initialize();
    allNSProblems["acceleration"]->setLevels(levelMin, levelMax);
    allNSProblems["acceleration"]->initialize();

    // === update all index sets ===
    updateIndexSets();
  }
}

void EventDriven::computeF(OneStepIntegrator* osi, integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{
  // Check osi type: only lsodar is allowed.
  if (osi->getType() != "Lsodar")
    RuntimeException::selfThrow("EventDriven::computeF(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  Lsodar * lsodar = static_cast<Lsodar*>(osi);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar->fillXWork(sizeOfX, x);


  double t = *time;
  model->setCurrentT(t);

  // update the DS of the OSI.
  lsodar->updateState(2); // update based on the last saved values for the DS state, ie the ones computed by lsodar (x above)

  // solve a LCP "acceleration"
  if (!allNSProblems.empty())
  {
    if (!(allNSProblems["acceleration"]->getInteractions()).isEmpty())
    {
      allNSProblems["acceleration"]->compute(t);
      updateInput(2); // Necessary to compute DS state below
    }
    // Compute the right-hand side ( xdot = f + Tu + r in DS) for all the ds, with the new value of input.
    lsodar->computeRhs(t);
  }

  // Get the required value, ie xdot for output.
  DynamicalSystemsSet dsOfTheOsi = lsodar->getDynamicalSystems();
  SiconosVector * xtmp2; // The Right-Hand side
  DSIterator it;
  unsigned int i = 0;
  for (it = dsOfTheOsi.begin(); it != dsOfTheOsi.end(); ++it)
  {
    xtmp2 = (*it)->getRhsPtr(); // Pointer link !
    for (unsigned int j = 0 ; j < (*it)->getN() ; ++j) // Warning: getN, not getDim !!!!
      xdot[i++] = (*xtmp2)(j);
  }
  xtmp2 = NULL;
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
    for (unsigned int j = 0 ; j < (*it)->getN() ; ++j)
    {
      for (unsigned k = 0 ; k < (*it)->getDim() ; ++k)
        jacob[i++] = (*jacotmp)(k, j);
    }
  }
}

void EventDriven::computeG(OneStepIntegrator* osi, integer * sizeOfX, doublereal* time, doublereal* x, integer * ng, doublereal * gOut)
{
  UnitaryRelationIterator itUR;
  unsigned int nsLawSize, k = 0 ;
  SiconosVector * y = NULL, * lambda = NULL;

  Lsodar * lsodar = static_cast<Lsodar*>(osi);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar->fillXWork(sizeOfX, x); // That may be not necessary? Check if computeF is called for each computeG.

  // IN: - lambda[2] obtained during LCP call in computeF()
  //     - y[0]: need to be updated.

  if (!allNSProblems.empty())
    updateOutput(0, 0);

  // If UR is in IndexSet2, g = lambda[2], else g = y[0]
  for (itUR = indexSets[0].begin(); itUR != indexSets[0].end(); ++itUR)
  {
    nsLawSize = (*itUR)->getNonSmoothLawSize();
    if (indexSets[2].find(*itUR) != indexSets[2].end()) // ie if the current Unitary Relation is in indexSets[2]
    {
      lambda = (*itUR)->getLambdaPtr(2);
      for (unsigned int i = 0; i < nsLawSize; i++)
        gOut[k++] = (*lambda)(i);
      //cout << " +++++++++++++++ IN G type lambda" << gOut[k-1] << endl;
    }
    else
    {
      y = (*itUR)->getYPtr(0);
      for (unsigned int i = 0; i < nsLawSize; i++)
        gOut[k++] = (*y)(i);
      //cout << " +++++++++++++++ IN G type y" << gOut[k-1] << endl;
    }
  }
}

void EventDriven::updateImpactState()
{
  OSIIterator itOSI;
  OSNSIterator itOsns;
  // Compute input = R(lambda[1])
  updateInput(1);

  // Compute post-impact velocity
  for (itOSI = allOSI.begin(); itOSI != allOSI.end() ; ++itOSI)
    (*itOSI)->updateState(1);
}

// Run the whole simulation
void EventDriven::run()
{
  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of Event Driven simulation - This may take a while ... ====" << endl;
  while (eventsManager->hasNextEvent())
  {
    computeOneStep();
    count++;
  }
  cout << "===== End of Event Driven simulation. " << count << " events have been processed. ==== " << endl;
}

void EventDriven::computeOneStep()
{
  advanceToEvent();
  eventsManager->processEvents();
}

void EventDriven::update(const unsigned int levelInput)
{
  if (!allNSProblems.empty())
  {
    // compute input
    updateInput(levelInput);

    OSIIterator itOSI;
    for (itOSI = allOSI.begin(); itOSI != allOSI.end() ; ++itOSI)
      (*itOSI)->updateState(levelInput);

    updateOutput(levelInput, levelInput);

    updateIndexSet(levelInput);
  }
}

void EventDriven::nextStep()
{

  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->nextStep();

  OSNSIterator itOsns;
  for (itOsns = allNSProblems.begin(); itOsns != allNSProblems.end(); ++itOsns)
    (itOsns->second)->nextStep();
}

void EventDriven::advanceToEvent()
{

  // WARNING: this is supposed to work for only one OSI, including all the DS.
  // To be reviewed for multiple OSI case (if it has sense?).

  // ---> Step 1: integrate the smooth dynamics from current event to next event;
  // Current event = last accessed event.
  // Next event = next time step or first root of the 'g' function found by integrator (Lsodar)

  // if istate == 1 => first call. It this case we suppose that tinit and tend have been initialized before.

  if (istate == 2) // ie no root found at previous step
  {
    tinit = eventsManager->getCurrentTime();
    tend =  eventsManager->getNextTime();
  }
  else if (istate == 3) // ie a root has been found at previous step => no changes in tint/tend values (done by integrate)
  {
    //istate = 1;
    tinit = eventsManager->getCurrentTime();
    tend =  eventsManager->getNextTime();
  }

  double tout = tend;
  bool isNewEventOccur = false;  // set to true if a new event occur during integration

  // call integrate method for each OSI, between tinit and tend.
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
  {
    //      DynamicalSystemsSet ds = (*it)->getDynamicalSystems();
    //      (*ds.begin())->getXPtr()->display();


    (*it)->integrate(tinit, tend, tout, istate); // integrate must return a flag telling if tend has been reached or not.

    if (istate == 3)
    {
      isNewEventOccur = true;
      // Add an event into the events manager list
      bool isScheduleOk = eventsManager->scheduleEvent("NonSmoothEvent", tout);
      if (!isScheduleOk) cout << " EventDriven advanceToEvent warning: try to add an already existing event" << endl;
    }
  }

  model->setCurrentT(tout);
}

EventDriven* EventDriven::convert(Simulation *str)
{
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}
