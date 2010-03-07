/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "EventDriven.hpp"
#include "SimulationXML.hpp"
#include "OneStepNSProblemXML.hpp"
#include "SimulationXML.hpp"
#include "UnitaryRelation.hpp"
#include "Lsodar.hpp"
#include "LCP.hpp"
#include "Model.hpp"
#include "Interaction.hpp"
#include "EventsManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "DynamicalSystem.hpp"

using namespace std;

// --- XML constructor ---
EventDriven::EventDriven(SP::SimulationXML strxml, double t0, double T,
                         SP::DynamicalSystemsSet dsList,
                         SP::InteractionsSet interList):
  Simulation(strxml, t0, T, dsList, interList), istate(1)
{
  // === One Step NS Problem === We read data in the xml output
  // (mainly Interactions concerned and solver) and assign them to
  // both one step ns problems ("acceleration" and "impact").  At the
  // time, only LCP is available for Event Driven.

  if (_simulationxml->hasOneStepNSProblemXML()) // ie if OSNSList is not empty
  {
    SetOfOSNSPBXML OSNSList = _simulationxml->getOneStepNSProblemsXML();
    SP::OneStepNSProblemXML osnsXML;
    string type;
    // For EventDriven, two OSNSPb are required, "acceleration" and
    // "impact"

    if (OSNSList.size() != 2)
      RuntimeException::selfThrow("EventDriven::xml constructor - Wrong number of OSNS problems: 2 are required.");
    string id = "impact";
    for (SetOfOSNSPBXMLIt it = OSNSList.begin(); it != OSNSList.end(); ++it)
    {
      osnsXML = *it;
      type = osnsXML->getNSProblemType();
      if (type == LCP_TAG)  // LCP
      {
        (*_allNSProblems)[id].reset(new LCP(osnsXML));
      }
      else
        RuntimeException::selfThrow("EventDriven::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");
      id = "acceleration";
    }

    (*_allNSProblems)["acceleration"]->setId("acceleration");
    (*_allNSProblems)["impact"]->setId("impact");
  }
}

void EventDriven::updateIndexSet(unsigned int i)
{
  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  assert(i < topo->indexSetsSize() &&
         "EventDriven::updateIndexSet(i), indexSets[i] does not exist.");

  // IndexSets[0] must not be updated by this function.
  assert(i > 0  &&
         "EventDriven::updateIndexSet(i=0), indexSets[0] can not be updated.");

  assert(topo->indexSet(i));
  assert(topo->indexSet(i - 1));

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]

  SP::UnitaryRelationsGraph indexSeti = topo->indexSet(i);
  SP::UnitaryRelationsGraph indexSetip = topo->indexSet(i - 1);

  double y;

  if (i == 1) // Action is different when first set is treated
  {
    // We get jroot array, output of root information Warning! work
    // only if there is only one osi in the simulation.
    if (_allOSI->size() > 1)
      RuntimeException::selfThrow
      ("EventDriven::updateIndexSet(i), not yet implemented for several OneStepIntegrators in the same Simulation process.");
    SP::Lsodar lsodar =
      boost::static_pointer_cast<Lsodar>(*(_allOSI->begin()));
    SA::integer jroot = lsodar->getJroot();
    unsigned int nsLawSize; // size of each UR, which corresponds to
    // the related nsLaw size
    unsigned int absolutePosition = 0; // global position of the UR
    // in the vector of
    // constraints, ie in jroot.
    UnitaryRelationsGraph::VIterator uip, uipend, vpnext;
    bool out = true;
    boost::tie(uip, uipend) = indexSetip->vertices();
    for (vpnext = uip;
         uip != uipend; uip = vpnext)
    {
      ++vpnext;

      SP::UnitaryRelation urp = indexSetip->bundle(*uip);
      nsLawSize = urp->getNonSmoothLawSize();

      // 1 - If the UR is not yet in the set

      if (!indexSeti->is_vertex(urp))
      {
        // Each UR may handle several constraints and thus
        // corresponds to a set of position in jroot
        for (unsigned int pos = absolutePosition;
             pos < (absolutePosition + nsLawSize); ++pos)
        {
          if (jroot[pos] == 1)
          {
            // vertex and edges insertions
            indexSeti->copy_vertex(urp, *indexSetip);
            break; // if one, at least, of the nsLawSize
            // constraints is active, the UR is
            // inserted into the set.

          }
        }
      }
      else // if the UR was already in the set
      {
        out = true;
        for (unsigned int pos = absolutePosition;
             pos < (absolutePosition + nsLawSize); ++pos)
        {
          if (jroot[pos] == 1)
          {
            out = false;
            break; // if one, at least, of the nsLawSize
            // constraints is active, the UR remains
            // in the set.
          }
        }
        if (out)
        {
          indexSeti->remove_vertex(urp);
          urp->lambda(i)->zero();
        }

      }
      absolutePosition += nsLawSize; // step to next UR ...
    }
  }
  else
  {
    UnitaryRelationsGraph::VIterator uip, uipend, vpnext;
    boost::tie(uip, uipend) = indexSetip->vertices();
    for (vpnext = uip; uip != uipend; uip = vpnext)
    {
      ++vpnext;

      // check if current Unitary Relation (ie *it) is in
      // indexSets[i]
      SP::UnitaryRelation ur = indexSetip->bundle(*uip);
      // Get y[i-1] double value
      y = ur->getYRef(i - 1);


      // if y[i-1] <=0, then the unitary relation is added in
      // indexSets[i] (if it was not already there) else if y[i-1]
      // > 0 and if the unitary relation was in the set, it is
      // removed.
      if (!indexSeti->is_vertex(ur)) // If the UR is not yet in
        // the set ...
      {
        if (fabs(y) <= _tolerance)
          // vertex and edges insertion
          indexSeti->copy_vertex(ur, *indexSetip);
      }
      else  // if the UR is already in the set
      {
        if (fabs(y) > _tolerance)
        {
          indexSeti->remove_vertex(ur);
          ur->lambda(i)->zero();
        }
      }
    }
  }
}

void EventDriven::updateIndexSetsWithDoubleCondition()
{

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  SP::UnitaryRelationsGraph indexSet2 = topo->indexSet(2);

  UnitaryRelationsGraph::VIterator ui, uiend, vnext;
  boost::tie(ui, uiend) = indexSet2->vertices();

  for (vnext = ui; ui != uiend; ui = vnext)
  {
    ++vnext;

    SP::UnitaryRelation ur = indexSet2->bundle(*ui);
    double gamma = ur->getYRef(2);
    double F     = ur->getLambdaRef(2);

    if (gamma > 0 && fabs(F) < _tolerance)
      indexSet2->remove_vertex(ur);
    else if (fabs(gamma) < _tolerance && fabs(F) < _tolerance) // undetermined case
      RuntimeException::selfThrow
      ("EventDriven::updateIndexSetsWithDoubleCondition(), undetermined case.");
  }
}

void EventDriven::initOSNS()
{

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]
  UnitaryRelationsGraph::VIterator ui, uiend;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  SP::UnitaryRelationsGraph indexSet0 = topo->indexSet(0);

  // For each Unitary relation in I0 ...
  for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    indexSet0->bundle(*ui)->initialize("EventDriven");

  if (!_allNSProblems->empty()) // ie if some Interactions have been
    // declared and a Non smooth problem
    // built.
  {
    // === OneStepNSProblem initialization. === First check that
    // there are 2 osns: one "impact" and one "acceleration"
    if (_allNSProblems->size() != 2)
      RuntimeException::selfThrow
      (" EventDriven::initialize, \n an EventDriven simulation must have two non smooth problem.\n Here, there are "
       + _allNSProblems->size());

    if (_allNSProblems->find("impact") ==
        _allNSProblems->end())  // ie if the impact problem does not
      // exist
      RuntimeException::selfThrow
      ("EventDriven::initialize, an EventDriven simulation must have an 'impact' non smooth problem.");
    if (_allNSProblems->find("acceleration") ==
        _allNSProblems->end()) // ie if the acceleration-level problem
      // does not exist
      RuntimeException::selfThrow
      ("EventDriven::initialize, an EventDriven simulation must have an 'acceleration' non smooth problem.");

    // At the time, we consider that for all systems, levelMin is
    // equal to the minimum value of the relative degree
    _levelMin = model()->nonSmoothDynamicalSystem()
                ->topology()->minRelativeDegree();
    if (_levelMin == 0)
      _levelMin++;

    // === update all index sets ===
    updateIndexSets();

    // WARNING: only for Lagrangian systems - To be reviewed for
    // other ones.
    (*_allNSProblems)["impact"]->setLevels(_levelMin - 1, _levelMax - 1);
    (*_allNSProblems)["impact"]->initialize(shared_from_this());
    (*_allNSProblems)["acceleration"]->setLevels(_levelMin, _levelMax);
    (*_allNSProblems)["acceleration"]->initialize(shared_from_this());
  }
}

void EventDriven::initLevelMax()
{
  _levelMax = model()->nonSmoothDynamicalSystem()->
              topology()->maxRelativeDegree();
  // Interactions initialization (here, since level depends on the
  // type of simulation) level corresponds to the number of Y and
  // Lambda derivatives computed.
  if (_levelMax == 0)
    _levelMax++;
}

void EventDriven::computef(SP::OneStepIntegrator osi, integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{

  // computeF is supposed to fill xdot in, using the definition of the
  // dynamical systems belonging to the osi

  // Check osi type: only lsodar is allowed.
  if (osi->getType() != OSI::LSODAR)
    RuntimeException::selfThrow("EventDriven::computef(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  SP::Lsodar lsodar = boost::static_pointer_cast<Lsodar>(osi);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar->fillXWork(sizeOfX, x);

  double t = *time;
  model()->setCurrentTime(t);

  // solve a LCP at "acceleration" level if required
  if (!_allNSProblems->empty())
  {
    if (!((*_allNSProblems)["acceleration"]->interactions())->isEmpty())
    {
      (*_allNSProblems)["acceleration"]->compute(t);
      updateInput(2); // Necessary to compute DS state below
    }
    // Compute the right-hand side ( xdot = f + r in DS) for all the
    //ds, with the new value of input.  lsodar->computeRhs(t);
  }
  // update the DS of the OSI.
  lsodar->computeRhs(t);
  //  lsodar->updateState(2); // update based on the last saved values
  //  for the DS state, ie the ones computed by lsodar (x above)

  // Update Index sets?

  // Get the required value, ie xdot for output.
  SP::SiconosVector xtmp2; // The Right-Hand side
  DSIterator it;
  unsigned int i = 0;
  for (it = lsodar->dynamicalSystemsBegin(); it != lsodar->dynamicalSystemsEnd(); ++it)
  {
    xtmp2 = (*it)->rhs(); // Pointer link !
    for (unsigned int j = 0 ; j < (*it)->getN() ; ++j) // Warning: getN, not getDim !!!!
      xdot[i++] = (*xtmp2)(j);
  }
}

void EventDriven::computeJacobianfx(SP::OneStepIntegrator osi,
                                    integer *sizeOfX,
                                    doublereal *time,
                                    doublereal *x,
                                    doublereal *jacob)
{
  if (osi->getType() != OSI::LSODAR)
    RuntimeException::selfThrow("EventDriven::computef(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  SP::Lsodar lsodar = boost::static_pointer_cast<Lsodar>(osi);

  // Remark A: according to DLSODAR doc, each call to jacobian is
  // preceded by a call to f with the same arguments NEQ, T, and Y.
  // Thus to gain some efficiency, intermediate quantities shared by
  // both calculations may be saved in class members?  fill in xWork
  // vector (ie all the x of the ds of this osi) with x fillXWork(x);
  // -> copy
  // Maybe this step is not necessary?  because of
  // remark A above

  // Compute the jacobian of the vector field according to x for the
  // current ds
  double t = *time;
  model()->setCurrentTime(t);
  lsodar->computeJacobianRhs(t);

  // Save jacobianX values from dynamical system into current jacob
  // (in-out parameter)

  DSIterator it;
  unsigned int i = 0;
  for (it = lsodar->dynamicalSystemsBegin(); it != lsodar->dynamicalSystemsEnd(); ++it)
  {
    SP::SiconosMatrix jacotmp = (*it)->jacobianRhsx(); // Pointer link !
    for (unsigned int j = 0 ; j < (*it)->getN() ; ++j)
    {
      for (unsigned k = 0 ; k < (*it)->getDim() ; ++k)
        jacob[i++] = (*jacotmp)(k, j);
    }
  }
}

void EventDriven::computeg(SP::OneStepIntegrator osi,
                           integer * sizeOfX, doublereal* time,
                           doublereal* x, integer * ng,
                           doublereal * gOut)
{

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());
  UnitaryRelationsGraph::VIterator ui, uiend;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();
  SP::UnitaryRelationsGraph indexSet0 = topo->indexSet(0);
  SP::UnitaryRelationsGraph indexSet2 = topo->indexSet(2);

  unsigned int nsLawSize, k = 0 ;

  SP::SiconosVector y, lambda;

  SP::Lsodar lsodar = boost::static_pointer_cast<Lsodar>(osi);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar->fillXWork(sizeOfX, x); // That may be not necessary? Check if
  // computeF is called for each
  // computeg.

  double t = *time;
  model()->setCurrentTime(t);

  // IN: - lambda[2] obtained during LCP call in computef()
  //     - y[0]: need to be updated.

  if (!_allNSProblems->empty())
    updateOutput(0, 0);

  // If UR is in IndexSet2, g = lambda[2], else g = y[0]
  for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet0->bundle(*ui);
    nsLawSize = ur->getNonSmoothLawSize();
    if (indexSet2->is_vertex(ur)) // ie if the current Unitary
      // Relation is in indexSets[2]
    {
      // Get lambda at acc. level, solution of LCP acc, called
      // during computef().
      lambda = ur->lambda(2);
      for (unsigned int i = 0; i < nsLawSize; i++)
        gOut[k++] = (*lambda)(i);
    }
    else
    {
      y = ur->y(0);
      for (unsigned int i = 0; i < nsLawSize; i++)
        gOut[k++] = (*y)(i);
    }
  }
}

void EventDriven::updateImpactState()
{
  OSIIterator itOSI;
  // Compute input = R(lambda[1])
  updateInput(1);

  // Compute post-impact velocity
  for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    (*itOSI)->updateState(1);
}

void EventDriven::update(unsigned int levelInput)
{
  if (!_allNSProblems->empty())
  {
    // compute input (lambda -> r)
    updateInput(levelInput);

    // Update dynamical systems states
    OSIIterator itOSI;
    for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
      (*itOSI)->updateState(levelInput);

    // Update output (y)
    updateOutput(levelInput, levelInput);
  }
  // Warning: index sets are not updated in this function !!
}

void EventDriven::advanceToEvent()
{
  // WARNING: this is supposed to work for only one OSI, including all
  // the DS.  To be reviewed for multiple OSI case (if it has sense?).

  // ---> Step 1: integrate the smooth dynamics from current event to
  // next event; Starting event = last accessed event.  Next event =
  // next time step or first root of the 'g' function found by
  // integrator (Lsodar)

  // if istate == 1 => first call. It this case we suppose that _tinit
  // and _tend have been initialized before.

  if (istate == 2 || istate == 3)
  {
    _tinit = _eventsManager->startingTime();
    _tend =  _eventsManager->nextTime();
  }

  _tout = _tend;
  bool isNewEventOccur = false;  // set to true if a new event occur
  // during integration

  // call integrate method for each OSI, between _tinit and _tend.
  OSIIterator it;
  for (it = _allOSI->begin(); it != _allOSI->end(); ++it)
  {
    (*it)->integrate(_tinit, _tend, _tout, istate); // integrate must
    // return a flag
    // (istate) telling if
    // _tend has been
    // reached or not.

    if (_printStat)
    {
      statOut << " =================> Results after advanceToEvent <================= " << endl;
      statOut << " Starting time: " << _tinit << endl;

    }
    if (istate == 3) // ie if _tout is not equal to _tend: one or more
      // roots have been found.
    {
      isNewEventOccur = true;
      // Add an event into the events manager list
      _eventsManager->scheduleNonSmoothEvent(_tout);
      if (_printStat)
        statOut << " -----------> New non-smooth event at time " << _tout << endl;
    }
    if (_printStat)
    {
      SP::Lsodar lsodar = boost::static_pointer_cast<Lsodar>(*it);
      statOut << "Results at time " << _tout << ":" << endl;
      SA::integer iwork = lsodar->getIwork();
      statOut << "Number of steps: " << iwork[10] << ", number of f evaluations: " << iwork[11] << ", number of jacobianF eval.: " << iwork[12] << "." << endl;

    }
  }


  // Set model time to _tout
  model()->setCurrentTime(_tout);
  // Update all the index sets ...
  updateOutput(1, 1);
  updateIndexSets();
}

EventDriven* EventDriven::convert(Simulation *str)
{
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}
