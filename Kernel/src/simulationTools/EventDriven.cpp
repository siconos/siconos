/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "Lsodar.hpp"
#include "LCP.hpp"
#include "Model.hpp"
#include "Interaction.hpp"
#include "EventsManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "DynamicalSystem.hpp"
#include "LagrangianDS.hpp"
#include "EventFactory.hpp"
#include "BlockMatrix.hpp"

#define DEBUG_MESSAGES 1

#include <debug.h>

using namespace std;
using namespace RELATION;

// --- XML constructor ---
EventDriven::EventDriven(SP::SimulationXML strxml, double t0, double T,
                         SP::DynamicalSystemsSet dsList,
                         SP::InteractionsSet interList):
  Simulation(strxml, t0, T, dsList, interList), _istate(1)
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
    int id = SICONOS_OSNSP_ED_IMPACT;
    (*_allNSProblems).resize(SICONOS_OSNSP_ED_NUMBER);
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
      id = SICONOS_OSNSP_ED_ACCELERATION;
    }


  }
}
/** defaut constructor
 *  \param a pointer to a timeDiscretisation (linked to the model that owns this simulation)
 */
EventDriven::EventDriven(SP::TimeDiscretisation td): Simulation(td), _istate(1)
{
  (*_allNSProblems).resize(SICONOS_OSNSP_ED_NUMBER);
}

EventDriven::EventDriven(SP::TimeDiscretisation td, int nb): Simulation(td), _istate(1)
{
  (*_allNSProblems).resize(nb);
}

double EventDriven::TOL_ED = DEFAULT_TOL_ED;

void EventDriven::updateIndexSet(unsigned int i)
{
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  assert(i < topo->indexSetsSize() &&
         "EventDriven::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0  &&
         "EventDriven::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);
  assert(indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  DEBUG_PRINTF("update indexSets start : indexSet0 size : %d\n", (int)indexSet0->size());
  DEBUG_PRINTF("update IndexSets start : indexSet1 size : %d\n", (int)indexSet1->size());
  DEBUG_PRINTF("update IndexSets start : indexSet2 size : %d\n", (int)indexSet2->size());

  InteractionsGraph::VIterator uibegin, uipend, uip;
  cpp11ns::tie(uibegin, uipend) = indexSet0->vertices();
  // loop over all vextice of the indexSet[i-1]
  for (uip = uibegin; uip != uipend; ++uip)
  {
    SP::Interaction inter = indexSet0->bundle(*uip);
    if (i == 1) // IndexSet[1]
    {
      // if indexSet[1]=>getYRef(0): output y
      // if indexSet[2]=>getYRef(1): output ydot
      double y = inter->getYRef(0); // output to define the IndexSets at this Interaction
      /*
         if (i == 1)
         {
         cout << "Id of Interaction: " << inter->number() << endl;
         cout << "Output of level 0 at this Interaction: " << y << endl;
         cout << endl;
         }
         */
      if (y < -TOL_ED) // y[0] < 0
      {
        RuntimeException::selfThrow("EventDriven::updateIndexSet, output of level 0 must be positive!!! ");
      }
      // 1 - If the Interaction is not yet in the set
      if (!indexSet1->is_vertex(inter)) // Interaction is not yet in the indexSet[i]
      {
        if (fabs(y) <= TOL_ED)
        {
          // vertex and edges insertions
          indexSet1->copy_vertex(inter, *indexSet0);
        }
      }
      else // if the Interaction was already in the set
      {
        if (fabs(y) > TOL_ED)
        {
          indexSet1->remove_vertex(inter); // remove the Interaction from IndexSet[1]
          inter->lambda(1)->zero(); // reset the lambda[1] to zero
        }
      }
    }
    else if (i == 2) // IndexSet[2]
    {
      if (indexSet1->is_vertex(inter)) // Interaction is in the indexSet[1]
      {
        double y = inter->getYRef(1); // output of level 1 at this Interaction
        if (!indexSet2->is_vertex(inter)) // Interaction is not yet in the indexSet[2]
        {
          if (fabs(y) <= TOL_ED)
          {
            // vertex and edges insertions
            indexSet2->copy_vertex(inter, *indexSet0);
          }
        }
        else // if the Interaction was already in the set
        {
          if (fabs(y) > TOL_ED)
          {
            indexSet2->remove_vertex(inter); // remove the Interaction from IndexSet[1]
            inter->lambda(2)->zero(); // reset the lambda[i] to zero
          }
        }
      }
      else // Interaction is not in the indexSet[1]
      {
        if (indexSet2->is_vertex(inter)) // Interaction is in the indexSet[2]
        {
          indexSet2->remove_vertex(inter); // remove the Interaction from IndexSet[2]
          inter->lambda(2)->zero(); // reset the lambda[i] to zero
        }
      }
    }
    else
    {
      RuntimeException::selfThrow("EventDriven::updateIndexSet, IndexSet[i > 2] doesn't exist");
    }
  }

  DEBUG_PRINTF("update indexSets end : indexSet0 size : %d\n", (int)indexSet0->size());
  DEBUG_PRINTF("update IndexSets end : indexSet1 size : %d\n", (int)indexSet1->size());
  DEBUG_PRINTF("update IndexSets end : indexSet2 size : %d\n", (int)indexSet2->size());
}

void EventDriven::updateIndexSetsWithDoubleCondition()
{

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  // for all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  SP::InteractionsGraph indexSet2 = topo->indexSet(2);

  InteractionsGraph::VIterator ui, uiend, vnext;
  cpp11ns::tie(ui, uiend) = indexSet2->vertices();

  for (vnext = ui; ui != uiend; ui = vnext)
  {
    ++vnext;

    SP::Interaction inter = indexSet2->bundle(*ui);
    double gamma = inter->getYRef(2);
    double F     = inter->getLambdaRef(2);
    if (fabs(F) < TOL_ED)
      indexSet2->remove_vertex(inter);
    else if ((gamma < -TOL_ED) || (F < -TOL_ED))
      RuntimeException::selfThrow("EventDriven::updateIndexSetsWithDoubleCondition(), output[2] and lambda[2] for Interactionof indexSet[2] must be higher or equal to zero.");
    else if (((fabs(gamma) > TOL_ED) && (fabs(F) > TOL_ED)))
      RuntimeException::selfThrow("EventDriven::updateIndexSetsWithDoubleCondition(), something is wrong for the LCP resolution.");
  }
}

void EventDriven::initializeInteraction(SP::Interaction inter)
{

  RELATION::TYPES pbType = inter->getRelationType();
  if (pbType == Lagrangian)
  {
    //    inter->setDataXFromVelocity();
  }
  else
    RuntimeException::selfThrow("EventDriven::initializeInteractions(SP::interaction inter) - not implemented for Relation of type " + pbType);

}




void EventDriven::initOSNS()
{

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());
  // === DS Rhs initialization ===

  for (OSIIterator itosi = _allOSI->begin();
       itosi != _allOSI->end(); ++itosi)
  {
    for (DSIterator itds = (*itosi)->dynamicalSystems()->begin();
         itds != (*itosi)->dynamicalSystems()->end();
         ++itds)
    {
      (*itds)->initRhs(model()->t0());
    }
  }

  // for all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]
  InteractionsGraph::VIterator ui, uiend;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  SP::InteractionsGraph indexSet0 = topo->indexSet(0);

  // For each Interaction in I0 ...
  for (cpp11ns::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    // indexSet0->bundle(*ui)->initialize("EventDriven");
    initializeInteraction(indexSet0->bundle(*ui));
  }
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

    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_IMPACT])) // ie if the impact problem does not
      // exist
      RuntimeException::selfThrow
      ("EventDriven::initialize, an EventDriven simulation must have an 'impact' non smooth problem.");
    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION])) // ie if the acceleration-level problem
      // does not exist
      RuntimeException::selfThrow
      ("EventDriven::initialize, an EventDriven simulation must have an 'acceleration' non smooth problem.");



    // === update all index sets ===
    updateIndexSets();

    // WARNING: only for Lagrangian systems - To be reviewed for
    // other ones.
    (*_allNSProblems)[SICONOS_OSNSP_ED_IMPACT]->setLevels(1, 1);
    (*_allNSProblems)[SICONOS_OSNSP_ED_IMPACT]->initialize(shared_from_this());
    (*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->setLevels(2, 2);
    (*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->initialize(shared_from_this());


    //====================================== added by Son Nguyen ===================================
    // Detect NonSmoothEvent at the beginning of the simulation
    SP::InteractionsGraph indexSet1 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    assert(indexSet1);
    if (indexSet1->size() >  0) // There is one non-smooth event to be added
    {
      //       EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
      //       _ENonSmooth = regEvent.instantiate(model()->currentTime(),2);
      //       _allEvents.insert(_ENonSmooth);
      //      _hasNS = true;
      _eventsManager->scheduleNonSmoothEvent(model()->currentTime());
    };
    //===============================================================================================

  }
}


// void EventDriven::initLevelMin()
// {
//   // At the time, we consider that for all systems, levelMin is
//   // equal to the minimum value of the relative degree
//   _levelMin = model()->nonSmoothDynamicalSystem()
//               ->topology()->minRelativeDegree();
//   if(_levelMin==0)
//     _levelMin++;
// }



// void EventDriven::initLevelMax()
// {
//   _levelMax = model()->nonSmoothDynamicalSystem()->topology()->maxRelativeDegree();
//   // Interactions initialization (here, since level depends on the
//   // type of simulation) level corresponds to the number of Y and
//   // Lambda derivatives computed.
//   if(_levelMax==0)
//     _levelMax++;
// }

void EventDriven::computef(SP::OneStepIntegrator osi, integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{
  //std::cout << "EventDriven::computef -------------------------> start" <<std::endl;

  // computeF is supposed to fill xdot in, using the definition of the
  // dynamical systems belonging to the osi

  // Check osi type: only lsodar is allowed.
  if (osi->getType() != OSI::LSODAR)
    RuntimeException::selfThrow("EventDriven::computef(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  Lsodar& lsodar = static_cast<Lsodar&>(*osi);
  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar.fillXWork(sizeOfX, x);

  double t = *time;
  model()->setCurrentTime(t);
  // solve a LCP at "acceleration" level if required
  if (!_allNSProblems->empty())
  {
    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->interactions())->isEmpty())
    {
      // Update the state of the DS
      (*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->compute(t);
      updateInput(2); // Necessary to compute DS state below
    }
    // Compute the right-hand side ( xdot = f + r in DS) for all the
    //ds, with the new value of input.  lsodar->computeRhs(t);
  }
  // update the DS of the OSI.
  lsodar.computeRhs(t);
  //  for the DS state, ie the ones computed by lsodar (x above)
  // Update Index sets? No !!

  // Get the required value, ie xdot for output.
  DSIterator it;
  unsigned int i = 0;
  for (it = lsodar.dynamicalSystemsBegin(); it != lsodar.dynamicalSystemsEnd(); ++it)
  {
    if (Type::value(**it) == Type::LagrangianDS ||
        Type::value(**it) == Type::LagrangianLinearTIDS)
    {
      LagrangianDS& LDS = *cpp11ns::static_pointer_cast<LagrangianDS>(*it);
      SiconosVector& qDotTmp = *LDS.velocity();
      SiconosVector& qDotDotTmp = *LDS.acceleration();
      for (unsigned int j = 0 ; j < (*it)->getDim() ; ++j)
        xdot[i++] = qDotTmp(j);
      for (unsigned int j = 0 ; j < (*it)->getDim() ; ++j)
        xdot[i++] = qDotDotTmp(j);
    }
    else
    {
      SiconosVector& xtmp2 = *(*it)->rhs(); // Pointer link !
      for (unsigned int j = 0 ; j < (*it)->getN() ; ++j) // Warning: getN, not getDim !!!!
        xdot[i++] = xtmp2(j);
    }
  }

  //std::cout << "EventDriven::computef -------------------------> stop" <<std::endl;

}

void EventDriven::computeJacobianfx(SP::OneStepIntegrator osi,
                                    integer *sizeOfX,
                                    doublereal *time,
                                    doublereal *x,
                                    doublereal *jacob)
{
  if (osi->getType() != OSI::LSODAR)
    RuntimeException::selfThrow("EventDriven::computeJacobianfx(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  SP::Lsodar lsodar = cpp11ns::static_pointer_cast<Lsodar>(osi);

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
    if (Type::value(**it) == Type::LagrangianDS ||
        Type::value(**it) == Type::LagrangianLinearTIDS)
    {
      LagrangianDS& lds = *cpp11ns::static_pointer_cast<LagrangianDS>(*it);
      BlockMatrix& jacotmp = *lds.jacobianRhsx();
      for (unsigned int j = 0; j < (*it)->getN(); ++j)
      {
        for (unsigned int k = 0; k < (*it)->getDim(); ++k)
          jacob[i++] = jacotmp(k, j);
      }
    }
    else
    {
      SiconosMatrix& jacotmp = *(*it)->jacobianRhsx(); // Pointer link !
      for (unsigned int j = 0; j < (*it)->getN(); ++j)
      {
        for (unsigned int k = 0; k < (*it)->getDim(); ++k)
          jacob[i++] = jacotmp(k, j);
      }
    }
  }
}

void EventDriven::computeg(SP::OneStepIntegrator osi,
                           integer * sizeOfX, doublereal* time,
                           doublereal* x, integer * ng,
                           doublereal * gOut)
{
  //std::cout << "EventDriven::computeg start" <<std::endl;
  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());
  InteractionsGraph::VIterator ui, uiend;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);
  unsigned int nsLawSize, k = 0 ;
  SP::SiconosVector y, ydot, lambda;
  SP::Lsodar lsodar = cpp11ns::static_pointer_cast<Lsodar>(osi);

  // Solve LCP at acceleration level to calculate the lambda[2] at Interaction of indexSet[2]
  lsodar->fillXWork(sizeOfX, x);
  //
  double t = *time;
  model()->setCurrentTime(t);
  if (!_allNSProblems->empty())
  {
    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->interactions())->isEmpty())
    {
      (*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->compute(t);
    }
  };
  /*
     double * xdottmp = (double *)malloc(*sizeOfX*sizeof(double));
     computef(osi, sizeOfX,time,x,xdottmp);
     free(xdottmp);
     */

  // Update the output from level 0 to level 1
  updateOutput(0);
  updateOutput(1);
  //
  for (cpp11ns::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui);
    nsLawSize = inter->getNonSmoothLawSize();
    y = inter->y(0);   // output y at this Interaction
    ydot = inter->y(1); // output of level 1 at this Interaction
    lambda = inter->lambda(2); // input of level 2 at this Interaction
    if (!(indexSet2->is_vertex(inter))) // if Interaction is not in the indexSet[2]
    {
      for (unsigned int i = 0; i < nsLawSize; ++i)
      {
        if ((*y)(i) >= TOL_ED) // y[0] > 0
        {
          gOut[k] = (*y)(i);
        }
        else // y[0] = 0
        {
          if ((*ydot)(i) >= 0) // if y[1] >= 0;
          {
            gOut[k] = TOL_ED; // g = TOL_ED temporarily
          }
          else  // if y[1] < 0
          {
            gOut[k] = (*y)(i); // g = y[0]
          }
        }
        k++;
      }
    }
    else // If Interaction is in the indexSet[2]
    {
      for (unsigned int i = 0; i < nsLawSize; ++i)
      {
        gOut[k] = (*lambda)(i); // g = lambda[2]
        k++;
      }
    }
  }
  //std::cout << "EventDriven::computeg stop" <<std::endl;
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
  assert(levelInput == 1);
  updateImpactState();
  // Update output (y)
  updateOutput(levelInput);
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

  // if _istate == 1 => first call. It this case we suppose that _tinit
  // and _tend have been initialized before
  if (_istate == 2 || _istate == 3)
  {
    _tinit = _eventsManager->startingTime();
    _tend =  _eventsManager->nextTime();
  }


  _tout = _tend;
  bool isNewEventOccur = false;  // set to true if a new event occur during integration
  // call integrate method for each OSI, between _tinit and _tend.
  OSIIterator it;
  for (it = _allOSI->begin(); it != _allOSI->end(); ++it)
  {

    (*it)->resetNonSmoothPart();

    //====================================================================================
    //    std::cout << " Start of Lsodar integration" << std::endl;
    (*it)->integrate(_tinit, _tend, _tout, _istate); // integrate must
    //    std::cout << " End of Lsodar integration" << std::endl;
    // return a flag (_istate) telling if _tend has been  reached or not.
    //====================================================================================

    if (_printStat)
    {
      statOut << " =================> Results after advanceToEvent <================= " << endl;
      statOut << " Starting time: " << _tinit << endl;
      statOut << " _istate " << _istate << endl;
    }
    if (_istate == 3) // ie if _tout is not equal to _tend: one or more roots have been found.
    {
      isNewEventOccur = true;
      // Add an event into the events manager list
      _eventsManager->scheduleNonSmoothEvent(_tout);
      if (_printStat)
        statOut << " -----------> New non-smooth event at time " << _tout << endl;
    }
    if (_printStat)
    {
      SP::Lsodar lsodar = cpp11ns::static_pointer_cast<Lsodar>(*it);
      statOut << "Results at time " << _tout << ":" << endl;
      SA::integer iwork = lsodar->getIwork();
      SA::doublereal Rwork = lsodar->getRwork();
      statOut << "Number of steps: " << iwork[10] << ", number of f evaluations: " << iwork[11] << ", number of jacobianF eval.: " << iwork[12] << "." << endl;
    }
  }
  // Set model time to _tout
  model()->setCurrentTime(_tout);
  //update output[0], output[1]
  updateOutput(0);
  updateOutput(1);
  // Update all the index sets ...
  updateIndexSets();
  //update lambda[2], input[2] and indexSet[2] with double consitions for the case there is no new event added during time integration, otherwise, this
  // update is done when the new event is processed
  if (!isNewEventOccur)
  {
    if (!_allNSProblems->empty())
    {
      // Solve LCP at acceleration level
      if (!((*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->interactions())->isEmpty())
      {
        (*_allNSProblems)[SICONOS_OSNSP_ED_ACCELERATION]->compute(_tout);
        updateInput(2); //
        //  updateInput(1); // this is done to reset the nonsmoothinput at the level of impact
      }
      // update indexSet[2] with double condition
      updateIndexSetsWithDoubleCondition();
    }
  }
}
//
EventDriven* EventDriven::convert(Simulation *str)
{
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}
