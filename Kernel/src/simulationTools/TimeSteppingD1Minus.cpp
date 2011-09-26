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
#include "UnitaryRelation.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

#include <debug.h>

using namespace std;
using namespace RELATION;

void TimeSteppingD1Minus::initOSNS()
{
  // initialize OSNS for UnitaryRelationsGraph from Topology
  assert(model()->nonSmoothDynamicalSystem()->topology()->isUpToDate());
  SP::Topology topo =  model()->nonSmoothDynamicalSystem()->topology();
  SP::UnitaryRelationsGraph indexSet0 = topo->indexSet(0);

  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    initializeInteraction(indexSet0->bundle(*ui)->interaction());
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

TimeSteppingD1Minus::TimeSteppingD1Minus(SP::TimeDiscretisation td, int nb) : Simulation(td)
{
  (*_allNSProblems).resize(nb);
}

TimeSteppingD1Minus::~TimeSteppingD1Minus()
{
}

void TimeSteppingD1Minus::updateIndexSet(unsigned int i)
{
  // To update IndexSet i: add or remove UnitaryRelations from
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

  // For all Unitary Relations in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::UnitaryRelationsGraph indexSet0 = topo->indexSet(0); // ALL UnitaryRelations : formula (8.30) of Acary2008
  SP::UnitaryRelationsGraph indexSet1 = topo->indexSet(1); // ACTIVE UnitaryRelations : formula (8.31) of Acary2008
  SP::UnitaryRelationsGraph indexSet2 = topo->indexSet(2); // STAYING ACTIVE UnitaryRelations : formula (8.32) of Acary2008
  assert(indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  topo->setHasChanged(false); // only with changed topology, OSNS will be forced to update themselves

  DEBUG_PRINTF("update indexSets start : indexSet0 size : %d\n", indexSet0->size());
  DEBUG_PRINTF("update IndexSets start : indexSet1 size : %d\n", indexSet1->size());
  DEBUG_PRINTF("update IndexSets start : indexSet2 size : %d\n", indexSet2->size());

  UnitaryRelationsGraph::VIterator uipend, uip;

  for (tie(uip, uipend) = indexSet0->vertices(); uip != uipend; ++uip) // loop over ALL
  {
    SP::UnitaryRelation urp = indexSet0->bundle(*uip);
    if (i == 1) // ACTIVE?
    {
      double y = urp->getYRef(0); // position

      if (!indexSet1->is_vertex(urp))
      {
        if (y <= 0.)
        {
          // if UnitaryRelation has not been active in the previous calculation and now becomes active
          indexSet1->copy_vertex(urp, *indexSet0);
          topo->setHasChanged(true);
        }
      }
      else
      {
        if (y > 0.)
        {
          // if UnitaryRelation has been active in the previous calculation and now becomes in-active
          indexSet1->remove_vertex(urp);
          topo->setHasChanged(true);
          urp->lambda(1)->zero(); // impuls is zero
        }
      }
    }
    else if (i == 2) // STAYING ACTIVE
    {
      if (indexSet1->is_vertex(urp)) // if UnitaryRelation is active
      {
        double y = urp->getYRef(1); // velocity

        if (!indexSet2->is_vertex(urp))
        {
          if (y <= 0.)
          {
            // if UnitaryRelation has not been staying active in the previous calculation and now becomes staying active
            indexSet2->copy_vertex(urp, *indexSet0);
            topo->setHasChanged(true);
          }
        }
        else
        {
          if (y > 0.)
          {
            // if UnitaryRelation has been staying active in the previous calculation and now does not stay active
            indexSet2->remove_vertex(urp);
            topo->setHasChanged(true);
            urp->lambda(2)->zero(); // force is zero
          }
        }
      }
      else
      {
        if (indexSet2->is_vertex(urp))
        {
          // if UnitaryRelation is in-active and has been staying active in previous calculation
          indexSet2->remove_vertex(urp);
          topo->setHasChanged(true);
          urp->lambda(2)->zero(); // force is zero
        }
      }
    }
    else
      RuntimeException::selfThrow("TimeSteppingD1Minus::updateIndexSet, IndexSet[i > 2] does not exist.");
  }
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
  // * calculate position (q_{k+1})
  // * calculate velocity (v_{k+1}^-)
  computeResidu();

  // calculate state without nonsmooth event with OSI
  // * calculate free velocity (!=v_{k+1}^-)
  computeFreeState();

  // event (impulse) calculation
  // * calculate gap velocity (g_{k+1}^-) with OSI
  // * calculate local impulse (Lambda_{k+1}^+)
  if (!_allNSProblems->empty())
    computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);

  // update on impulse level
  // * calculate global impulse (p_{k+1}^+)
  // * update velocity (v_{k+1}^+) with OSI
  // * calculate local gaps (g_{k+1}^+)
  update(1);

  // indexset (I_{k+1}^+) is calculated in Simulation::processEvent
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
