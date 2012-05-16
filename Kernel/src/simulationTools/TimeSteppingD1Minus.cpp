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

#include "TimeSteppingD1Minus.hpp"
#include "D1MinusLinear.hpp"
#include "TimeDiscretisation.hpp"
#include "Topology.hpp"
//#include "Interaction.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "debug.h"

using namespace std;
using namespace RELATION;

void TimeSteppingD1Minus::initOSNS()
{
  // initialize OSNS for InteractionsGraph from Topology
  assert(model()->nonSmoothDynamicalSystem()->topology()->isUpToDate());
  SP::Topology topo =  model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);

  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    initializeInteraction(indexSet0->bundle(*ui));
  }

  // there is at least one OSNP
  if (!_allNSProblems->empty())
  {
    if (_allNSProblems->size() != 2)
      RuntimeException::selfThrow("TimeSteppingD1Minus::initOSNS, TimeSteppingD1Minus simulation must have two OneStepNonsmoothProblems.");

    //update all index sets
    updateIndexSets();

    // set evaluation levels (first is of velocity, second of acceleration type)
    (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->setLevels(1, 1);
    (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->initialize(shared_from_this());
    (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY + 1]->setLevels(2, 2);
    (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(shared_from_this());

    // update output
    for (unsigned int level = _levelMinForOutput; level < _levelMaxForOutput; level++)
      updateOutput(level);
  }
}

void TimeSteppingD1Minus::initializeInteraction(SP::Interaction inter)
{
  for (DSIterator it = inter->dynamicalSystemsBegin(); it != inter->dynamicalSystemsEnd(); ++it)
    inter->workZ()->insertPtr((*it)->z());

  RELATION::TYPES pbType = inter->relation()->getType();
  if (pbType == Lagrangian)
  {
    for (DSIterator it = inter->dynamicalSystemsBegin(); it != inter->dynamicalSystemsEnd(); ++it)
    {
      inter->workX()->insertPtr((boost::static_pointer_cast<LagrangianDS>(*it))->velocity());
      inter->workFree()->insertPtr((boost::static_pointer_cast<LagrangianDS>(*it))->workFree());
    }
  }
  else
    RuntimeException::selfThrow("TimeSteppingD1Minus::initializeInteractions - not implemented for Relation of type " + pbType);

}

TimeSteppingD1Minus::TimeSteppingD1Minus(SP::TimeDiscretisation td, int nb) : Simulation(td), impactOccuredLastTimeStep(false)
{
  (*_allNSProblems).resize(nb);
}

TimeSteppingD1Minus::~TimeSteppingD1Minus()
{
}

void TimeSteppingD1Minus::updateIndexSet(unsigned int i)
{
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeSteppingD1Minus::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0 &&
         "TimeSteppingD1Minus::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::InteractionsGraph indexSet0 = topo->indexSet(0); // ALL Interactions : formula (8.30) of Acary2008
  SP::InteractionsGraph indexSet1 = topo->indexSet(1); // ACTIVE Interactions for IMPACTS
  SP::InteractionsGraph indexSet2 = topo->indexSet(2); // ACTIVE Interactions for CONTACTS
  assert(indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  topo->setHasChanged(false); // only with changed topology, OSNS will be forced to update themselves

  DEBUG_PRINT("\nINDEXSETS BEFORE UPDATE\n");
  DEBUG_EXPR(indexSet0->display());
  DEBUG_EXPR(indexSet1->display());
  DEBUG_EXPR(indexSet2->display());

  InteractionsGraph::VIterator uipend, uip;

  for (tie(uip, uipend) = indexSet0->vertices(); uip != uipend; ++uip) // loop over ALL
  {
    SP::Interaction inter = indexSet0->bundle(*uip);

    if (i == 1) // ACTIVE FOR IMPACT CALCULATIONS? Contacts which have been closing in the last time step
    {
      impactOccuredLastTimeStep = false;
      DEBUG_PRINT("\nUPDATE INDEXSET 1\n");

      double y = (*(inter->y(0)))(0); // current position
      double yOld = (*(inter->yOld(0)))(0); // old position

      DEBUG_PRINTF("y= %f\n", y);
      DEBUG_PRINTF("yOld= %f\n", yOld);

      if (!indexSet1->is_vertex(inter))
      {
        if (y <= DEFAULT_TOL_D1MINUS && yOld > DEFAULT_TOL_D1MINUS)
        {
          // if Interaction has not been active in the previous calculation and now becomes active
          indexSet1->copy_vertex(inter, *indexSet0);
          topo->setHasChanged(true);
        }
      }
      else
      {
        indexSet1->remove_vertex(inter);
        topo->setHasChanged(true);
        impactOccuredLastTimeStep = true;
        inter->lambda(1)->zero(); // impuls is zero
      }
    }
    else if (i == 2) // ACTIVE FOR CONTACT CALCULATIONS? Contacts which are closed but have not been closing in the last time step
    {
      DEBUG_PRINT("\nUPDATE INDEXSET 2\n");

      double y = (*(inter->y(0)))(0); // current position

      DEBUG_PRINTF("y= %f\n", y);

      if (indexSet2->is_vertex(inter))
      {
        if (y > DEFAULT_TOL_D1MINUS)
        {
          // if Interaction has been active in the previous calculation and now becomes in-active
          indexSet2->remove_vertex(inter);
          topo->setHasChanged(true);
          inter->lambda(2)->zero(); // force is zero
        }
      }
      else
      {
        if (y <= DEFAULT_TOL_D1MINUS && !indexSet1->is_vertex(inter) && !impactOccuredLastTimeStep)
        {
          // if Interaction has is active but has not become active recently
          indexSet2->copy_vertex(inter, *indexSet0);
          topo->setHasChanged(true);
        }
      }
    }
    else
      RuntimeException::selfThrow("TimeSteppingD1Minus::updateIndexSet, IndexSet[i > 2] does not exist.");
  }

  DEBUG_PRINT("\nINDEXSETS AFTER UPDATE\n");
  DEBUG_EXPR(indexSet0->display());
  DEBUG_EXPR(indexSet1->display());
  DEBUG_EXPR(indexSet2->display());
}

void TimeSteppingD1Minus::update(unsigned int levelInput)
{
  // compute input (lambda -> r)
  if (!_allNSProblems->empty())
    updateInput(levelInput);

  // compute state for each dynamical system
  for (OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end(); ++itOSI)
    (*itOSI)->updateState(levelInput);

  // compute output (x -> y)
  if (!_allNSProblems->empty())
  {
    for (unsigned int level = _levelMinForOutput; level < _levelMaxForOutput; level++)
      updateOutput(level);
  }
}

void TimeSteppingD1Minus::run()
{
  unsigned int count = 0;
  cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" << endl;
  while (_eventsManager->hasNextEvent())
  {
    advanceToEvent();

    _eventsManager->processEvents();
    count++;
  }
  cout << "===== End of " << Type::name(*this) << "simulation. " << count << " events have been processed. ==== " << endl;
}

void TimeSteppingD1Minus::advanceToEvent()
{
  // we start after initialization (initOSNS) with
  // * initial state (q_0, v_0^+)
  // * updated indexset (I_0^+)
  // * updated  gaps and gap velocities (g_0^+)
  //
  // hence we end this procedure with
  // * state (q_{k+1}, v_{k+1}^+)
  // * updated gaps and gap velocities (g_{k+1}^+)
  // * indexset (I_{k+1}^+)

  // calculate residu without nonsmooth event with OSI
  // * calculate position q_{k+1} in ds->q()
  // * calculate velocity v_{k+1}^- and not free velocity in ds->velocity()
  // * calculate free residu in ds->freeResidu()
  computeResidu();

  // calculate state without nonsmooth event with OSI
  // * calculate free velocity and not v_{k+1}^- in ds->velocity
  computeFreeState();

  // event (impulse) calculation only when there has been a topology change (here: closing contact)
  // * calculate gap velocity using free velocity with OSI
  // * calculate local impulse (Lambda_{k+1}^+)
  //
  // Maurice Bremond: indices must be recomputed
  // as we deal with dynamic graphs, vertices and edges are stored
  // in lists for fast add/remove during updateIndexSet(i)
  // we need indices of list elements to build the OSNS Matrix so we
  // need an update if graph has changed
  // this should be done in updateIndexSet(i) for all integrators only
  // if a graph has changed
  //updateIndexSet(1);
  //model()->nonSmoothDynamicalSystem()->topology()->indexSet(1)->update_vertices_indices();
  //model()->nonSmoothDynamicalSystem()->topology()->indexSet(1)->update_edges_indices();

  //if(model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
  //{
  //  for(OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
  //  {
  //    (*itOsns)->setHasBeenUpdated(false);
  //  }
  //}

  if (!_allNSProblems->empty())
    computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);

  // update on impulse level
  // * calculate global impulse (p_{k+1}^+)
  // * update velocity (v_{k+1}^+) with OSI in ds->velocity
  // * calculate local gaps (g_{k+1}^+)
  update(1);

  // indexset (I_{k+1}^+) is calculated in Simulation::processEvent
}

void TimeSteppingD1Minus::updateInput(unsigned int level)
{
  //  assert(level>=0);

  double time = model()->currentTime();
  SP::Topology topology = model()->nonSmoothDynamicalSystem()->topology();
  InteractionsIterator it;

  // set dynamical systems non-smooth part to zero.
  for (OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end(); ++itOSI)
  {
    for (DSIterator itDS = (*itOSI)->dynamicalSystems()->begin(); itDS != (*itOSI)->dynamicalSystems()->end(); ++itDS)
    {
      Type::Siconos dsType = Type::value(**itDS);
      if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
        RuntimeException::selfThrow("TimeSteppingD1Minus::updateInput - not implemented for Dynamical system type: " + dsType);
      else
      {
        SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (*itDS);
        if (d->p(level)) d->p(level)->zero();
      }
    }
  }

  // we compute input using lambda(level).
  for (it = topology->interactions()->begin(); it != topology->interactions()->end(); it++)
  {
    assert((*it)->lowerLevelForInput() <= level);
    assert((*it)->upperLevelForInput() >= level);
    (*it)->computeInput(time, level);
  }
}

void TimeSteppingD1Minus::computeResidu()
{
  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->computeResidu();
}

void TimeSteppingD1Minus::computeFreeState()
{
  std::for_each(_allOSI->begin(), _allOSI->end(), boost::bind(&OneStepIntegrator::computeFreeState, _1));
}

TimeSteppingD1Minus* TimeSteppingD1Minus::convert(Simulation *str)
{
  return dynamic_cast<TimeSteppingD1Minus*>(str);
}
