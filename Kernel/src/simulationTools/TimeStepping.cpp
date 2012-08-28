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

#include "TimeStepping.hpp"
#include "SimulationXML.hpp"
#include "OneStepNSProblemXML.hpp"
#include "Topology.hpp"
#include "LCP.hpp"
#include "Relay.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
//#include "Interaction.hpp"
#include "OneStepIntegrator.hpp"
#include "Interaction.hpp"
#include "EventsManager.hpp"
#include "FrictionContact.hpp"
#include "FirstOrderNonLinearDS.hpp"

//#define DEBUG_MESSAGES 1

#include <debug.h>

using namespace std;
using namespace RELATION;

/** Pointer to function, used to set the behavior of simulation when
    ns solver failed.  If equal to null, use DefaultCheckSolverOutput
    else (set with setCheckSolverFunction) call the pointer below).
    Note FP: (temporary) bad method to set checkSolverOutput but it
    works ... It may be better to use plug-in?
*/
static CheckSolverFPtr checkSolverOutput = NULL;

TimeStepping::TimeStepping(SP::TimeDiscretisation td,
                           SP::OneStepIntegrator osi,
                           SP::OneStepNSProblem osnspb)
  : Simulation(td), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonOptions(SICONOS_TS_NONLINEAR), _newtonResiduDSMax(0.0), _newtonResiduYMax(0.0), _newtonResiduRMax(0.0)
{

  _computeResiduY = false;
  _computeResiduR = false;


  if (osi) insertIntegrator(osi);
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TS);
  if (osnspb) insertNonSmoothProblem(osnspb, SICONOS_OSNSP_TS_VELOCITY);

}

TimeStepping::TimeStepping(SP::TimeDiscretisation td, int nb)
  : Simulation(td), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonOptions(SICONOS_TS_NONLINEAR), _newtonResiduDSMax(0.0), _newtonResiduYMax(0.0), _newtonResiduRMax(0.0)
{
  _computeResiduY = false;
  _computeResiduR = false;

  (*_allNSProblems).resize(nb);
}

// --- XML constructor ---
TimeStepping::TimeStepping(SP::SimulationXML strxml, double t0,
                           double T, SP::DynamicalSystemsSet dsList,
                           SP::InteractionsSet interList):
  Simulation(strxml, t0, T, dsList, interList), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonOptions(SICONOS_TS_NONLINEAR), _newtonResiduDSMax(0.0), _newtonResiduYMax(0.0), _newtonResiduRMax(0.0)
{
  _computeResiduY = false;
  _computeResiduR = false;
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TS);
  // === One Step NS Problem === For time stepping, only one non
  // smooth problem is built.
  if (_simulationxml->hasOneStepNSProblemXML())  // ie if OSNSList is
    // not empty
  {
    SetOfOSNSPBXML OSNSList = _simulationxml->getOneStepNSProblemsXML();
    if (OSNSList.size() != 1)
      RuntimeException::selfThrow("TimeStepping::xml constructor - Two many inputs for OSNS problems (only one problem is required).");
    SP::OneStepNSProblemXML osnsXML = *(OSNSList.begin());
    // OneStepNSProblem - Memory allocation/construction
    string type = osnsXML->getNSProblemType();
    if (type == LCP_TAG)  // LCP
    {
      (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY].reset(new LCP(osnsXML));
    }
    else if (type == FRICTIONCONTACT_TAG)
    {
      (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY].reset(new FrictionContact(osnsXML));
    }
    else RuntimeException::selfThrow("TimeStepping::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");

    //      (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->setId("timeStepping");

    // Add QP and Relay cases when these classes will be fully
    // implemented.
  }
}

// --- Destructor ---
TimeStepping::~TimeStepping()
{
}

// bool TimeStepping::predictorDeactivate(SP::Interaction inter, unsigned int i)
// {
//   double h = timeStep();
//   double y = inter->getYRef(i-1); // for i=1 it is the position -> historic notation y
//   double yDot = inter->getYRef(i); // for i=1 it is the velocity -> historic notation yDot
//   DEBUG_PRINTF("TS::predictorDeactivate yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y+0.5*h*yDot);
//   y += 0.5*h*yDot;
//   assert(!isnan(y));
//   if(y>0)
//     DEBUG_PRINTF("TS::predictorDeactivate DEACTIVATE.\n");
//   return (y>0);
// }

// bool TimeStepping::predictorActivate(SP::Interaction inter, unsigned int i)
// {
//   double h = timeStep();
//   double y = inter->getYRef(i-1); // for i=1 it is the position -> historic notation y
//   double yDot = inter->getYRef(i); // for i=1 it is the velocity -> historic notation yDot
//   DEBUG_PRINTF("TS::predictorActivate yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y+0.5*h*yDot);
//   y += 0.5*h*yDot;
//   assert(!isnan(y));
//   if(y<=0)
//     DEBUG_PRINTF("TS::predictorActivate ACTIVATE.\n");
//   return (y<=0);
// }

void TimeStepping::updateIndexSet(unsigned int i)
{
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.
  // boost::default_color_type is used to organize update in InteractionsGraph:
  // - white_color : undiscovered vertex (Interaction)
  // - gray_color : discovered vertex (Interaction) but searching descendants
  // - black_color : discovered vertex (Interaction) together with the descendants

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0 &&
         "TimeStepping::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  assert(indexSet0);
  assert(indexSet1);

  topo->setHasChanged(false);

  DEBUG_PRINTF("update indexSets start : indexSet0 size : %ld\n", indexSet0->size());
  DEBUG_PRINTF("update IndexSets start : indexSet1 size : %ld\n", indexSet1->size());

  // Check indexSet1
  InteractionsGraph::VIterator ui1, ui1end, v1next;
  boost::tie(ui1, ui1end) = indexSet1->vertices();

  //Remove interactions from the indexSet1
  for (v1next = ui1; ui1 != ui1end; ui1 = v1next)
  {
    ++v1next;

    SP::Interaction inter1 = indexSet1->bundle(*ui1);
    if (indexSet0->is_vertex(inter1))
    {
      InteractionsGraph::VDescriptor inter1_descr0 = indexSet0->descriptor(inter1);

      assert((indexSet0->color(inter1_descr0) == boost::white_color));

      indexSet0->color(inter1_descr0) = boost::gray_color;
      if (Type::value(*(inter1->nonSmoothLaw())) != Type::EqualityConditionNSL)
      {
        SP::OneStepIntegrator Osi = integratorOfInteraction(inter1);
        //if(predictorDeactivate(inter1,i))
        if (Osi->removeInteractionInIndexSet(inter1, i))
        {
          // Interaction is not active
          // ui1 becomes invalid
          indexSet0->color(inter1_descr0) = boost::black_color;

          indexSet1->eraseProperties(*ui1);

          InteractionsGraph::OEIterator oei, oeiend;
          for (boost::tie(oei, oeiend) = indexSet1->out_edges(*ui1);
               oei != oeiend; ++oei)
          {
            InteractionsGraph::EDescriptor ed1, ed2;
            boost::tie(ed1, ed2) = indexSet1->edges(indexSet1->source(*oei), indexSet1->target(*oei));
            if (ed2 != ed1)
            {
              indexSet1->eraseProperties(ed1);
              indexSet1->eraseProperties(ed2);
            }
            else
            {
              indexSet1->eraseProperties(ed1);
            }
          }


          indexSet1->remove_vertex(inter1);
          /* \warning V.A. 25/05/2012 : Multiplier lambda are only set to zero if they are removed from the IndexSet*/
          inter1->lambda(1)->zero();
          topo->setHasChanged(true);
        }
      }
    }
    else
    {
      // Interaction is not in indexSet0 anymore.
      // ui1 becomes invalid
      indexSet1->eraseProperties(*ui1);
      InteractionsGraph::OEIterator oei, oeiend;
      for (boost::tie(oei, oeiend) = indexSet1->out_edges(*ui1);
           oei != oeiend; ++oei)
      {
        InteractionsGraph::EDescriptor ed1, ed2;
        boost::tie(ed1, ed2) = indexSet1->edges(indexSet1->source(*oei), indexSet1->target(*oei));
        if (ed2 != ed1)
        {
          indexSet1->eraseProperties(ed1);
          indexSet1->eraseProperties(ed2);
        }
        else
        {
          indexSet1->eraseProperties(ed1);
        }
      }

      indexSet1->remove_vertex(inter1);
      topo->setHasChanged(true);
    }
  }

  // indexSet0\indexSet1 scan
  InteractionsGraph::VIterator ui0, ui0end;
  //Add interaction in indexSet1
  for (boost::tie(ui0, ui0end) = indexSet0->vertices(); ui0 != ui0end; ++ui0)
  {
    if (indexSet0->color(*ui0) == boost::black_color)
    {
      // reset
      indexSet0->color(*ui0) = boost::white_color ;
    }
    else
    {
      if (indexSet0->color(*ui0) == boost::gray_color)
      {
        // reset
        indexSet0->color(*ui0) = boost::white_color;

        assert(indexSet1->is_vertex(indexSet0->bundle(*ui0)));
        /*assert( { !predictorDeactivate(indexSet0->bundle(*ui0),i) ||
          Type::value(*(indexSet0->bundle(*ui0)->nonSmoothLaw())) == Type::EqualityConditionNSL ;
          });*/
      }
      else
      {
        assert(indexSet0->color(*ui0) == boost::white_color);

        SP::Interaction inter0 = indexSet0->bundle(*ui0);
        assert(!indexSet1->is_vertex(inter0));
        bool activate = true;
        if (Type::value(*(inter0->nonSmoothLaw())) != Type::EqualityConditionNSL)
        {
          SP::OneStepIntegrator Osi = integratorOfInteraction(inter0);
          activate = Osi->addInteractionInIndexSet(inter0, i);
        }
        if (activate)
        {
          assert(!indexSet1->is_vertex(inter0));

          // vertex and edges insertion in indexSet1
          indexSet1->copy_vertex(inter0, *indexSet0);
          topo->setHasChanged(true);
          assert(indexSet1->is_vertex(inter0));
        }
      }
    }
  }

  assert(indexSet1->size() <= indexSet0->size());

  DEBUG_PRINTF("update indexSets end : indexSet0 size : %ld\n", indexSet0->size());
  DEBUG_PRINTF("update IndexSets end : indexSet1 size : %ld\n", indexSet1->size());
}

// void TimeStepping::insertNonSmoothProblem(SP::OneStepNSProblem osns)
// {
//   // A the time, a time stepping simulation can only have one non
//   // smooth problem.
//   if((*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY])
//      RuntimeException::selfThrow
//        ("TimeStepping,  insertNonSmoothProblem - A non smooth problem already exist. You can not have more than one.");

//   (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY] = osns;
// }

void TimeStepping::initializeInteraction(SP::Interaction inter)
{


  RELATION::TYPES pbType = inter->relation()->getType();
  if (pbType == FirstOrder)
  {
  }
}



void TimeStepping::initOSNS()
{
  // === creates links between work vector in OSI and work vector in
  // Interactions
  SP::OneStepIntegrator  osi;

  ConstDSIterator itDS;

  SP::Topology topo =  model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);

  InteractionsGraph::VIterator ui, uiend;

  // For each Interaction in I0 ...
  for (boost::tie(ui, uiend) = indexSet0->vertices();
       ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui);
    // indexSet0->bundle(*ui)->initialize("TimeStepping");
    initializeInteraction(indexSet0->bundle(*ui));



  }


  if (!_allNSProblems->empty()) // ie if some Interactions have been
    // declared and a Non smooth problem
    // built.
  {
    //if (_allNSProblems->size()>1)
    //  RuntimeException::selfThrow("TimeStepping::initialize, at the time, a time stepping simulation can not have more than one non smooth problem.");

    // At the time, we consider that for all systems, levelMin is
    // equal to the minimum value of the relative degree - 1 except
    // for degree 0 case where we keep 0.

    assert(model()->nonSmoothDynamicalSystem()->topology()->isUpToDate());


    // === update all index sets ===
    updateIndexSets();

    // initialization of  OneStepNonSmoothProblem
    for (OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    {
      (*itOsns)->setLevels(_levelMinForInput, _levelMaxForInput);
      (*itOsns)->initialize(shared_from_this());
    }
  }
}

// void TimeStepping::initLevelMin()
// {
//   assert(model()->nonSmoothDynamicalSystem()->topology()->minRelativeDegree()>=0);

//   _levelMin = model()->nonSmoothDynamicalSystem()->topology()->minRelativeDegree();

//   if(_levelMin!=0)
//     _levelMin--;
// }

// void TimeStepping::initLevelMax()
// {
//   _levelMax = model()->nonSmoothDynamicalSystem()->topology()->maxRelativeDegree();
//   // Interactions initialization (here, since level depends on the
//   // type of simulation) level corresponds to the number of Y and
//   // Lambda derivatives computed.

//   if(_levelMax!=0)
//     _levelMax--;
//   // level max is equal to relative degree-1. But for relative degree 0 case, we keep 0 value for _levelMax
// }

void TimeStepping::nextStep()
{
  processEvents();
}

void TimeStepping::update(unsigned int levelInput)
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
    for (unsigned int level = _levelMinForOutput;
         level < _levelMaxForOutput + 1;
         level++)
      updateOutput(level);
  }
}

void TimeStepping::computeFreeState()
{
  std::for_each(_allOSI->begin(), _allOSI->end(), boost::bind(&OneStepIntegrator::computeFreeState, _1));
}

// compute simulation between current and next event.  Initial
// DS/interaction state is given by memory vectors and final state is
// the one saved in DS/Interaction at the end of this function
void TimeStepping::computeOneStep()
{
  advanceToEvent();
}


void TimeStepping::computeInitialResidu()
{
  //  cout<<"BEGIN computeInitialResidu"<<endl;
  double tkp1 = getTkp1();

  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();

  if (!allInteractions->isEmpty())
  {
    //    assert(_levelMinForOutput >=0);
    assert(_levelMaxForOutput >= _levelMinForOutput);
    //    assert(_levelMinForInput >=0);
    assert(_levelMaxForInput >= _levelMinForInput);


    updateOutput(_levelMinForOutput);
    updateInput(_levelMaxForInput);
  }

  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    dsGraph->bundle(*vi)->updatePlugins(tkp1);
  }

  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->computeResidu();

  if (_computeResiduY)
    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->computeResiduY(tkp1);
    }

  //  cout<<"END computeInitialResidu"<<endl;
}

void TimeStepping::run()
{
  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of
  // events manager.
  cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" << endl;
  while (_eventsManager->hasNextEvent())
  {

    //  if (_newtonOptions=="linear" || )
    //       advanceToEvent();

    //     else if (_newtonOptions=="Newton")
    //       newtonSolve(criterion,maxIter);

    //     else
    //       RuntimeException::selfThrow("TimeStepping::run(opt) failed. Unknow simulation option: "+opt);

    advanceToEvent();

    _eventsManager->processEvents();
    count++;
  }
  cout << "===== End of " << Type::name(*this) << "simulation. " << count << " events have been processed. ==== " << endl;
}

void TimeStepping::advanceToEvent()
{
  //   computeInitialResidu();
  //   /*advance To Event consists of one Newton iteration, here the jacobians are updated.*/
  //   prepareNewtonIteration();
  //   // solve ...
  //   computeFreeState();
  //   int info = 0;
  //   if (!_allNSProblems->empty())
  //     info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
  //   // Check output from solver (convergence or not ...)
  //   if (!checkSolverOutput)
  //     DefaultCheckSolverOutput(info);
  //   else
  //     checkSolverOutput(info, this);
  //   // Update
  //   update(_levelMin);

  newtonSolve(_newtonTolerance, _newtonMaxIteration);



}

/*update of the nabla */
/*discretisation of the Interactions */
void   TimeStepping::prepareNewtonIteration()
{
  //  cout << "update the operators" <<endl ;


  for (OSIIterator itosi = _allOSI->begin();
       itosi != _allOSI->end(); ++itosi)
  {
    (*itosi)->prepareNewtonIteration(getTkp1());
  }


  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
  for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  {
    (*it)->relation()->computeJach(getTkp1(), **it);
    (*it)->relation()->computeJacg(getTkp1(), **it);
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
    (*it)->preparNewtonIteration();
  }
  bool topoHasChanged = model()->nonSmoothDynamicalSystem()->topology()->hasChanged();
  if (topoHasChanged)
  {
    for (OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    {
      (*itOsns)->setHasBeenUpdated(false);
    }
  }
}
void TimeStepping::saveYandLambdaInMemory()
{
  // Save OSNS state (Interactions) in Memory.
  OSNSIterator itOsns;
  for (itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    (*itOsns)->saveInMemory();

}
void TimeStepping::newtonSolve(double criterion, unsigned int maxStep)
{
  bool isNewtonConverge = false;
  _newtonNbSteps = 0; // number of Newton iterations
  int info = 0;
  //cout<<"||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||| BEGIN NEWTON IT "<<endl;
  bool isLinear  = (_model.lock())->nonSmoothDynamicalSystem()->isLinear();
  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();

  computeInitialResidu();

  if ((_newtonOptions == SICONOS_TS_LINEAR || _newtonOptions == SICONOS_TS_LINEAR_IMPLICIT)
      || isLinear)
  {
    _newtonNbSteps++;
    prepareNewtonIteration();
    computeFreeState();
    if (!_allNSProblems->empty() &&  !allInteractions->isEmpty())
      info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
    // Check output from solver (convergence or not ...)
    if (!checkSolverOutput)
      DefaultCheckSolverOutput(info);
    else
      checkSolverOutput(info, this);

    update(_levelMaxForInput);

    //isNewtonConverge = newtonCheckConvergence(criterion);
    if (!_allNSProblems->empty() &&  !allInteractions->isEmpty())
      saveYandLambdaInMemory();
  }

  else if (_newtonOptions == SICONOS_TS_NONLINEAR)
  {
    //  while((!isNewtonConverge)&&(_newtonNbSteps < maxStep)&&(!info))
    while ((!isNewtonConverge) && (_newtonNbSteps < maxStep))
    {
      _newtonNbSteps++;
      prepareNewtonIteration();
      computeFreeState();
      if (info)
        cout << "New Newton loop because of nonsmooth solver failed\n" << endl;

      // if there is not any Interaction at
      // the beginning of the simulation _allNSProblems may not be
      // empty here (check with SpaceFilter and one disk not on
      // the ground : MultiBodyTest::t2)

      // if((*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->simulation())
      // is also relevant here.
      if (!_allNSProblems->empty() && !allInteractions->isEmpty())
      {
        info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
      }
      //if(info)
      //  cout<<"info!"<<endl;
      // Check output from solver (convergence or not ...)
      if (!checkSolverOutput)
        DefaultCheckSolverOutput(info);
      else
        checkSolverOutput(info, this);

      update(_levelMaxForInput);
      isNewtonConverge = newtonCheckConvergence(criterion);

      if (!isNewtonConverge && !info)
      {
        if (!_allNSProblems->empty() &&  !allInteractions->isEmpty())
          saveYandLambdaInMemory();
      }
    }
    if (!isNewtonConverge)
    {
      cout << "TimeStepping::newtonSolve -- Newton process stopped: max. number of steps (" << maxStep << ") reached." << endl ;
      if (info)
        cout << "TimeStepping::newtonSolve -- nonsmooth solver failed." << endl ;
    }
    else
    {
      //      cout << "TimeStepping::newtonSolve succeed nbit="<<_newtonNbSteps<<"maxStep="<<maxStep<<endl;
    }
  }
  else
    RuntimeException::selfThrow("TimeStepping::NewtonSolve failed. Unknow newtonOptions: " + _newtonOptions);
}

bool TimeStepping::newtonCheckConvergence(double criterion)
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
      (*it)->computeResiduY(getTkp1());
      residu = (*it)->residuY()->norm2();
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
      (*it)->computeResiduR(getTkp1());
      residu = (*it)->residuR()->norm2();

      if (residu > _newtonResiduRMax) _newtonResiduRMax = residu;
      if (residu > criterion)
      {
        //cout<<"residuR > criteron"<<residu<<">"<<criterion<<endl;
        checkConvergence = false;
        //break;
      }
      //else
      //  cout<<"residuR ="<<residu<<"<"<<criterion<<endl;
    }

  }

  return(checkConvergence);
}

TimeStepping* TimeStepping::convert(Simulation *str)
{
  TimeStepping* ts = dynamic_cast<TimeStepping*>(str);
  return ts;
}

void TimeStepping::DefaultCheckSolverOutput(int info)
{
  // info = 0 => ok
  // else: depend on solver
  if (info != 0)
  {
    std::cout << "TimeStepping::DefaultCheckSolverOutput:" << std::endl;
    std::cout << "Non smooth solver warning/error: output message from solver is equal to " << info << std::endl;
    //       cout << "=> may have failed? (See Numerics solver documentation for details on the message meaning)." << endl;
    //      cout << "=> may have failed? (See Numerics solver documentation for details on the message meaning)." << endl;
    //     RuntimeException::selfThrow(" Non smooth problem, solver convergence failed ");
    /*      if(info == 1)
            cout <<" reach max iterations number with solver " << solverName << endl;
            else if (info == 2)
            {
            if (solverName == "LexicoLemke" || solverName == "CPG" || solverName == "NLGS")
            RuntimeException::selfThrow(" negative diagonal term with solver "+solverName);
            else if (solverName == "QP" || solverName == "NSQP" )
            RuntimeException::selfThrow(" can not satisfy convergence criteria for solver "+solverName);
            else if (solverName == "Latin")
            RuntimeException::selfThrow(" Choleski factorisation failed with solver Latin");
            }
            else if (info == 3 && solverName == "CPG")
            cout << "pWp null in solver CPG" << endl;
            else if (info == 3 && solverName == "Latin")
            RuntimeException::selfThrow("Null diagonal term with solver Latin");
            else if (info == 5 && (solverName == "QP" || solverName == "NSQP"))
            RuntimeException::selfThrow("Length of working array insufficient in solver "+solverName);
            else
            RuntimeException::selfThrow("Unknown error type in solver "+ solverName);
    */
  }
}

void TimeStepping::setCheckSolverFunction(CheckSolverFPtr newF)
{
  checkSolverOutput = newF;
}
