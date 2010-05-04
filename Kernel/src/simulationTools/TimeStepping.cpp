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

#include "TimeStepping.hpp"
#include "SimulationXML.hpp"
#include "OneStepNSProblemXML.hpp"
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

using namespace std;

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
  : Simulation(td)
{

  mComputeResiduY = false;

  if (osi) insertIntegrator(osi);
  (*_allNSProblems).resize(SICONOS_OSNSP_TS_NUMBER);
  if (osnspb) insertNonSmoothProblem(osnspb);

}

// --- XML constructor ---
TimeStepping::TimeStepping(SP::SimulationXML strxml, double t0,
                           double T, SP::DynamicalSystemsSet dsList,
                           SP::InteractionsSet interList):
  Simulation(strxml, t0, T, dsList, interList)
{
  mComputeResiduY = false;
  (*_allNSProblems).resize(SICONOS_OSNSP_TS_NUMBER);
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


// i=1
void TimeStepping::updateIndexSet(unsigned int i)
{
  // To update IndexSet number i: add or remove UnitaryRelations from
  // this set, depending on y values.

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");

  assert(i != 0 &&   // IndexSets[0] must not be updated in simulation,
         // since it belongs to the Topology.
         "TimeStepping::updateIndexSet(i=0), indexSets[0] can not be updated.");

  assert(i == 1);  // yes

  assert(topo->indexSet(0));
  assert(topo->indexSet(1));

  SP::UnitaryRelationsGraph indexSet0 = topo->indexSet(0);
  SP::UnitaryRelationsGraph indexSet1 = topo->indexSet(1);



  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]
  double y, yDot;
  bool inserted;

  // For all UR in Index[i-1] ...
  double h = timeStep();

  // indexSet1 scan
  UnitaryRelationsGraph::VIterator ui1, ui1end, v1next;
  boost::tie(ui1, ui1end) = indexSet1->vertices();
  //Remove interactions from the indexSet1
  for (v1next = ui1 ;
       ui1 != ui1end; ui1 = v1next)
  {
    ++v1next;

    SP::UnitaryRelation ur1 = indexSet1->bundle(*ui1);
    if (indexSet0->is_vertex(ur1))
    {
      UnitaryRelationsGraph::VDescriptor ur1_descr0 = indexSet0->descriptor(ur1);

      assert((indexSet0->color(ur1_descr0)
              == boost::white_color));

      indexSet0->color(ur1_descr0) = boost::gray_color;
      if (Type::value(*(ur1->interaction()->nonSmoothLaw())) != Type::EqualityConditionNSL)
      {
        y = ur1->getYRef(i - 1);
        yDot = ur1->getYRef(1);
        y += 0.5 * h * yDot;
        if (y > 0)
        {
          // Unitary relation is not active
          // ui1 becomes invalid
          indexSet0->color(ur1_descr0) = boost::black_color;
          indexSet1->remove_vertex(ur1);
          ur1->lambda(1)->zero();
        }
      }
    }
    else
    {
      // Unitary relation is not in indexSet0 anymore.
      // ui1 becomes invalid
      indexSet1->remove_vertex(ur1);
    }
  }

  // indexSet0\indexSet1 scan
  UnitaryRelationsGraph::VIterator ui0, ui0end;
  //Add interaction in indexSet1
  for (boost::tie(ui0, ui0end) = indexSet0->vertices();
       ui0 != ui0end; ++ui0)
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
        indexSet0->color(*ui0) = boost::white_color ;

        assert(indexSet1->is_vertex(indexSet0->bundle(*ui0)));
        assert( { y = indexSet0->bundle(*ui0)->getYRef(i - 1);
                  yDot = indexSet0->bundle(*ui0)->getYRef(1);
                  y += 0.5 * h*yDot;
                  y <= 0 || Type::value(*(indexSet0->bundle(*ui0)->interaction()->nonSmoothLaw())) == Type::EqualityConditionNSL ;
                });
      }

      else
      {
        assert(indexSet0->color(*ui0) == boost::white_color);

        SP::UnitaryRelation ur0 = indexSet0->bundle(*ui0);
        assert(!indexSet1->is_vertex(ur0));

        if (Type::value(*(ur0->interaction()->nonSmoothLaw())) != Type::EqualityConditionNSL)
        {
          y = ur0->getYRef(i - 1);
          yDot = ur0->getYRef(1);
          y += 0.5 * h * yDot;
        }
        else
          y = -1;
        if (y <= 0)
        {
          assert(!indexSet1->is_vertex(ur0));

          // vertex and edges insertion in indexSet1
          indexSet1->copy_vertex(ur0, *indexSet0);

          assert(indexSet1->is_vertex(ur0));
        }
      }
    }
  }

  indexSet1->update_vertices_indices();
  indexSet1->update_edges_indices();

  assert(indexSet1->size() <= indexSet0->size());

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

void TimeStepping::initOSNS()
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
    indexSet0->bundle(*ui)->initialize("TimeStepping");
    // creates a POINTER link between workX[ds] (xfree) and the
    // corresponding unitaryBlock in each UR for each ds of the
    // current UR.
    for (itDS = ur->dynamicalSystemsBegin();
         itDS != ur->dynamicalSystemsEnd(); ++itDS)
    {
      //osi = osiMap[*itDS];
      ur->insertInWorkFree((*itDS)->workFree()); // osi->getWorkX(*itDS));
    }
  }


  if (!_allNSProblems->empty()) // ie if some Interactions have been
    // declared and a Non smooth problem
    // built.
  {
    if (_allNSProblems->size() > 1)
      RuntimeException::selfThrow("TimeStepping::initialize, at the time, a time stepping simulation can not have more than one non smooth problem.");

    // At the time, we consider that for all systems, levelMin is
    // equal to the minimum value of the relative degree - 1 except
    // for degree 0 case where we keep 0.

    assert(model()->nonSmoothDynamicalSystem()->topology()->isUpToDate());
    assert(model()->nonSmoothDynamicalSystem()->topology()->minRelativeDegree() >= 0);

    _levelMin = model()->nonSmoothDynamicalSystem()->topology()->minRelativeDegree();

    if (_levelMin != 0)
      _levelMin--;

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

void TimeStepping::initLevelMax()
{
  _levelMax = model()->nonSmoothDynamicalSystem()->topology()->maxRelativeDegree();
  // Interactions initialization (here, since level depends on the
  // type of simulation) level corresponds to the number of Y and
  // Lambda derivatives computed.

  if (_levelMax != 0) // level max is equal to relative degree-1. But for
    // relative degree 0 case, we keep 0 value for
    // _levelMax
    _levelMax--;
}

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

  // 3 - compute output ( x ... -> y)
  if (!_allNSProblems->empty())
  {
    updateOutput(0, _levelMax);
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
  for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  {
    (*it)->relation()->computeg(tkp1);
    (*it)->relation()->computeh(tkp1);
  }

  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    dsGraph->bundle(*vi)->updatePlugins(tkp1);
  }


  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->computeResidu();
  if (mComputeResiduY)
    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->relation()->computeResiduY(tkp1);
    }

  //  cout<<"END computeInitialResidu"<<endl;
}

void TimeStepping::advanceToEvent()
{
  computeInitialResidu();

  // solve ...
  computeFreeState();
  int info = 0;
  if (!_allNSProblems->empty())
    info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
  // Check output from solver (convergence or not ...)
  if (!checkSolverOutput)
    DefaultCheckSolverOutput(info);
  else
    checkSolverOutput(info, this);
  // Update
  update(_levelMin);
}

/*update of the nabla */
/*discretisation of the Interactions */
void   TimeStepping::prepareNewtonIteration()
{
  //  cout << "update the operators" <<endl ;

  DSOSIConstIterator it = _osiMap.begin();
  while (it != _osiMap.end())
  {
    if ((it->second)->getType() == OSI::MOREAU)
    {
      Moreau::convert(&(*(it->second)))->computeW(getTkp1(), it->first);
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
  mNbNewtonSteps = 0; // number of Newton iterations
  //double residu = 0;
  int info = 0;
  //  cout<<"||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||| BEGIN NEWTON IT"<<endl;
  computeInitialResidu();

  while ((!isNewtonConverge || info) && (mNbNewtonSteps <= maxStep))
  {
    mNbNewtonSteps++;
    prepareNewtonIteration();
    computeFreeState();
    if (info)
      cout << "new loop because of info\n" << endl;
    if (!_allNSProblems->empty())
      info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
    if (info)
      cout << "info!" << endl;
    // Check output from solver (convergence or not ...)
    if (!checkSolverOutput)
      DefaultCheckSolverOutput(info);
    else
      checkSolverOutput(info, this);

    update(_levelMin);
    isNewtonConverge = newtonCheckConvergence(criterion);
    if (!isNewtonConverge && !info)
    {
      saveYandLambdaInMemory();
    }
  }
  if (!isNewtonConverge && !info)
    cout << "Newton process stopped: max. steps number reached." << endl ;
  //  cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  END NEWTON IT"<<endl;
}

bool TimeStepping::newtonCheckConvergence(double criterion)
{
  bool checkConvergence = true;
  //mRelativeConvergenceCriterionHeld is true means that the RCC is activated, and the relative criteron helds.
  //In this case the newtonCheckConvergence has to return true. End of the Newton iterations
  if (mRelativeConvergenceCriterionHeld)
  {
    return true;
  }
  // get the nsds indicator of convergence
  // We compute cvg = abs(xi+1 - xi)/xi and if cvg < criterion
  //  if (nsdsConverge < criterion )

  double residu = 0.0;
  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
  {
    residu = (*it)->computeResidu();
    if (residu > criterion)
    {
      checkConvergence = false;
      //break;
    }
  }
  if (!mComputeResiduY)
    return(checkConvergence);


  //check residuy.
  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
  for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  {
    (*it)->relation()->computeResiduY(getTkp1());
    residu = (*it)->relation()->getResiduY()->norm2();
    if (residu > criterion)
    {
      //      cout<<"residuY > criteron"<<residu<<">"<<criterion<<endl;
      checkConvergence = false;
      //break;
    }
  }
  return(checkConvergence);
}

void TimeStepping::run(const std::string& opt, double criterion, unsigned int maxIter)
{
  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of
  // events manager.
  cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" << endl;
  while (_eventsManager->hasNextEvent())
  {
    if (opt == "linear")
      advanceToEvent();

    else if (opt == "Newton")
      newtonSolve(criterion, maxIter);

    else
      RuntimeException::selfThrow("TimeStepping::run(opt) failed. Unknow simulation option: " + opt);

    _eventsManager->processEvents();
    count++;
  }
  cout << "===== End of " << Type::name(*this) << "simulation. " << count << " events have been processed. ==== " << endl;
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
    cout << "TimeStepping::check non smooth solver output warning: output message from solver is equal to " << info << endl;
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
