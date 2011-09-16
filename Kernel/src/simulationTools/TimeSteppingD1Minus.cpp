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
#include "SimulationXML.hpp"
#include "Topology.hpp"
#include "LCP.hpp"
#include "Relay.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "UnitaryRelation.hpp"
#include "OneStepIntegrator.hpp"
#include "Interaction.hpp"
#include "EventsManager.hpp"
#include "FrictionContact.hpp"
#include "Moreau.hpp"
#include "D1MinusLinear.hpp"

#include <debug.h>

using namespace std;

TimeSteppingD1Minus::TimeSteppingD1Minus(SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb)
  : Simulation(td), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonOptions(SICONOS_TS_NONLINEAR)
{

  _computeResiduY = false;
  _computeResiduR = false;

  if (osi) insertIntegrator(osi);
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TS);
  if (osnspb) insertNonSmoothProblem(osnspb, SICONOS_OSNSP_TS_VELOCITY);
}

TimeSteppingD1Minus::TimeSteppingD1Minus(SP::TimeDiscretisation td, int nb)
  : Simulation(td), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonOptions(SICONOS_TS_NONLINEAR)
{
  _computeResiduY = false;
  _computeResiduR = false;

  (*_allNSProblems).resize(nb);
}

// --- Destructor ---
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

  topo->setHasChanged(false); // TODO NOETIG?

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
          urp->lambda(1)->zero(); // TODO WAS ZU NULL SETZEN?
        }
      }
    }
    //else if(i == 2) // STAYING ACTIVE
    //{
    //  if(indexSet1->is_vertex(urp)) // if UnitaryRelation is active
    //  {
    //    double y = urp->getYRef(1); // velocity

    //    if(!indexSet2->is_vertex(urp))
    //    {
    //      if(y <= 0.)
    //      { // if UnitaryRelation has not been staying active in the previous calculation and now becomes staying active
    //        indexSet2->copy_vertex(urp, *indexSet0);
    //        topo->setHasChanged(true);
    //      }
    //    }
    //    else
    //    {
    //      if(y > 0.)
    //      { // if UnitaryRelation has been staying active in the previous calculation and now does not stay active
    //        indexSet2->remove_vertex(urp);
    //        topo->setHasChanged(true);
    //        urp->lambda(2)->zero(); // TODO WAS ZU NULL SETZEN?
    //      }
    //    }
    //  }
    //  else
    //  {
    //    if(indexSet2->is_vertex(urp))
    //    { // if UnitaryRelation is in-active and has been staying active in previous calculation
    //      indexSet2->remove_vertex(urp);
    //      topo->setHasChanged(true);
    //      urp->lambda(2)->zero(); // TODO WAS ZU NULL SETZEN?
    //    }
    //  }
    //}
    //else
    //  RuntimeException::selfThrow("TimeSteppingD1Minus::updateIndexSet, IndexSet[i > 2] does not exist.");
  }
}

void TimeSteppingD1Minus::initOSNS()
{
  // === creates links between work vector in OSI and work vector in
  // Unitary Relations
  SP::OneStepIntegrator  osi;

  ConstDSIterator itDS;

  SP::Topology topo =  model()->nonSmoothDynamicalSystem()->topology();
  SP::UnitaryRelationsGraph indexSet0 = topo->indexSet(0);

  UnitaryRelationsGraph::VIterator ui, uiend;

  // For each Unitary relation in I0 ...
  for (boost::tie(ui, uiend) = indexSet0->vertices();
       ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet0->bundle(*ui);
    indexSet0->bundle(*ui)->initialize("TimeSteppingD1Minus");
    // creates a POINTER link between workX[ds] (xfree) and the
    // corresponding unitaryBlock in each UR for each ds of the
    // current UR.
    for (itDS = ur->interaction()->dynamicalSystemsBegin();
         itDS != ur->interaction()->dynamicalSystemsEnd(); ++itDS)
    {
      //osi = osiMap[*itDS];
      ur->insertInWorkFree((*itDS)->workFree()); // osi->getWorkX(*itDS));
    }
  }

  if (!_allNSProblems->empty()) // ie if some Interactions have been
    // declared and a Non smooth problem
    // built.
  {
    assert(model()->nonSmoothDynamicalSystem()->topology()->isUpToDate());

    initLevelMin();

    // === update all index sets ===
    updateIndexSets();

    // initialization of  OneStepNonSmoothProblem
    for (OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    {
      (*itOsns)->setLevels(_levelMin, _levelMax);
      (*itOsns)->initialize(shared_from_this());
    }
  }
}

void TimeSteppingD1Minus::initLevelMin()
{
  assert(model()->nonSmoothDynamicalSystem()->topology()->minRelativeDegree() >= 0);

  _levelMin = model()->nonSmoothDynamicalSystem()->topology()->minRelativeDegree();

  if (_levelMin != 0)
    _levelMin--;
}

void TimeSteppingD1Minus::initLevelMax()
{
  _levelMax = model()->nonSmoothDynamicalSystem()->topology()->maxRelativeDegree();
  // Interactions initialization (here, since level depends on the
  // type of simulation) level corresponds to the number of Y and
  // Lambda derivatives computed.

  if (_levelMax == 0)
    _levelMax++;
  // like event driven scheme
}

void TimeSteppingD1Minus::nextStep()
{
  processEvents();
}


void TimeSteppingD1Minus::update(unsigned int levelInput)
{
  // 1 - compute input (lambda -> r)
  if (!_allNSProblems->empty())
    updateInput(levelInput);

  // 2 - compute state for each dynamical system
  OSIIterator itOSI;
  for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    (*itOSI)->updateState(levelInput);
  /*Because the dof of DS have been updated,
    the world (CAO for example) must be updated.*/
  updateWorldFromDS();

  // 3 - compute output ( x ... -> y)
  if (!_allNSProblems->empty())
  {
    updateOutput(0, _levelMax);
  }
}

void TimeSteppingD1Minus::computeFreeState()
{
  std::for_each(_allOSI->begin(), _allOSI->end(), boost::bind(&OneStepIntegrator::computeFreeState, _1));
}

// compute simulation between current and next event.  Initial
// DS/interaction state is given by memory vectors and final state is
// the one saved in DS/Interaction at the end of this function
void TimeSteppingD1Minus::computeOneStep()
{
  advanceToEvent();
}

void TimeSteppingD1Minus::computeInitialResidu()
{
  //  cout<<"BEGIN computeInitialResidu"<<endl;
  double tkp1 = getTkp1();

  double time = model()->currentTime();
  assert(abs(time - tkp1) < 1e-14);

  updateOutput(0, _levelMax);
  updateInput(_levelMin);

  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    dsGraph->bundle(*vi)->updatePlugins(tkp1);
  }

  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();

  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->computeResidu();
  if (_computeResiduY)
    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->relation()->computeResiduY(tkp1);
    }
}

void TimeSteppingD1Minus::run()
{
  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of
  // events manager.
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
  newtonSolve(_newtonTolerance, _newtonMaxIteration);
}

/*update of the nabla */
/*discretisation of the Interactions */
void   TimeSteppingD1Minus::prepareNewtonIteration()
{
  DSOSIConstIterator it = _osiMap.begin();
  while (it != _osiMap.end())
  {
    if ((it->second)->getType() == OSI::MOREAU)
    {
      Moreau::convert(&(*(it->second)))->computeW(getTkp1(), it->first);
    }
    if ((it->second)->getType() == OSI::D1MINUSLINEAR)
    {
      D1MinusLinear::convert(&(*(it->second)))->computeW(getTkp1(), it->first);
    }
    ++it;
  }

  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
  for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  {
    (*it)->relation()->computeJach(getTkp1());
    (*it)->relation()->computeJacg(getTkp1());
  }

  /*reset to zero the ds buffers*/
  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();


  /* should be evaluated only if needed */
  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    dsGraph->bundle(*vi)->preparStep();
    //     (*itds)->xp()->zero();
    //     (*itds)->R()->zero();
  }
  /**/
  for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  {
    (*it)->relation()->preparNewtonIteration();
  }
  bool topoHasChanged = model()->nonSmoothDynamicalSystem()->topology()->hasChanged();
  if (topoHasChanged)
    for (OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    {
      (*itOsns)->setHasBeUpdated(false);
    }

}
void TimeSteppingD1Minus::saveYandLambdaInMemory()
{
  // Save OSNS state (Interactions) in Memory.
  OSNSIterator itOsns;
  for (itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    (*itOsns)->saveInMemory();

}
void TimeSteppingD1Minus::newtonSolve(double criterion, unsigned int maxStep)
{
  bool isNewtonConverge = false;
  _newtonNbSteps = 0; // number of Newton iterations
  int info = 0;
  bool isLinear  = (_model.lock())->nonSmoothDynamicalSystem()->isLinear();
  computeInitialResidu();

  if ((_newtonOptions == SICONOS_TS_LINEAR || _newtonOptions == SICONOS_TS_LINEAR_IMPLICIT)
      || isLinear)
  {
    _newtonNbSteps++;
    prepareNewtonIteration();
    computeFreeState();
    if (!_allNSProblems->empty())
      info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
    // Check output from solver (convergence or not ...)
    DefaultCheckSolverOutput(info);

    update(_levelMin);

    saveYandLambdaInMemory();
  }

  else if (_newtonOptions == SICONOS_TS_NONLINEAR)
  {
    while ((!isNewtonConverge) && (_newtonNbSteps < maxStep) && (!info))
    {
      _newtonNbSteps++;
      prepareNewtonIteration();
      computeFreeState();
      if (info)
        cout << "new loop because of info\n" << endl;
      if (!_allNSProblems->empty())
        info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
      if (info)
        cout << "info!" << endl;
      // Check output from solver (convergence or not ...)
      DefaultCheckSolverOutput(info);

      update(_levelMin);
      isNewtonConverge = newtonCheckConvergence(criterion);
      if (!isNewtonConverge && !info)
      {
        saveYandLambdaInMemory();
      }
    }
    if (!isNewtonConverge)
      cout << "TimeSteppingD1Minus::newtonSolve -- Newton process stopped: max. number of steps  reached." << endl ;
    else if (info)
      cout << "TimeSteppingD1Minus::newtonSolve -- Newton process stopped: solver failed." << endl ;
    //    else
    //      cout << "TimeSteppingD1Minus::newtonSolve succed nbit="<<_newtonNbSteps<<"maxStep="<<maxStep<<endl;
  }
  else
    RuntimeException::selfThrow("TimeSteppingD1Minus::NewtonSolve failed. Unknow newtonOptions: " + _newtonOptions);
}

bool TimeSteppingD1Minus::newtonCheckConvergence(double criterion)
{
  bool checkConvergence = true;
  //_relativeConvergenceCriterionHeld is true means that the RCC is
  //activated, and the relative criteron helds.  In this case the
  //newtonCheckConvergence has to return true. End of the Newton
  //iterations
  if (_relativeConvergenceCriterionHeld)
  {
    return true;
  }
  // get the nsds indicator of convergence
  // We compute cvg = abs(xi+1 - xi)/xi and if cvg < criterion
  //  if (nsdsConverge < criterion )

  double residu = 0.0;
  _newtonResiduDSMax = 0.0;
  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
  {
    residu = (*it)->computeResidu();

    if (residu > _newtonResiduDSMax) _newtonResiduDSMax = residu;
    if (residu > criterion)
    {
      checkConvergence = false;
      //break;
    }
  }

  if (_computeResiduY)
  {
    //check residuy.
    _newtonResiduYMax = 0.0;
    residu = 0.0;
    SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->relation()->computeResiduY(getTkp1());
      residu = (*it)->relation()->residuY()->norm2();
      if (residu > _newtonResiduYMax) _newtonResiduYMax = residu;
      if (residu > criterion)
      {
        //      cout<<"residuY > criteron"<<residu<<">"<<criterion<<endl;
        checkConvergence = false;
        //break;
      }
    }
  }
  if (_computeResiduR)
  {
    //check residur.
    _newtonResiduRMax = 0.0;
    residu = 0.0;
    SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->relation()->computeResiduR(getTkp1());
      residu = (*it)->relation()->residuR()->norm2();
      cout << "residuR =" << residu << ">" << criterion << endl;
      if (residu > _newtonResiduRMax) _newtonResiduRMax = residu;
      if (residu > criterion)
      {
        cout << "residuR > criteron" << residu << ">" << criterion << endl;
        checkConvergence = false;
        //break;
      }
    }
  }

  return(checkConvergence);
}

TimeSteppingD1Minus* TimeSteppingD1Minus::convert(Simulation *str)
{
  TimeSteppingD1Minus* ts = dynamic_cast<TimeSteppingD1Minus*>(str);
  return ts;
}

void TimeSteppingD1Minus::DefaultCheckSolverOutput(int info)
{
  // info = 0 => ok
  // else: depend on solver
  if (info != 0)
  {
    cout << "TimeSteppingD1Minus::check non smooth solver output warning: output message from solver is equal to " << info << endl;
  }
}
