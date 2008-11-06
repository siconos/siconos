/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "Model.h"
#include "TimeDiscretisation.h"
#include "NonSmoothDynamicalSystem.h"
#include "UnitaryRelation.h"
#include "OneStepIntegrator.h"
#include "Interaction.h"
#include "EventsManager.h"
#include "FrictionContact.h"

using namespace std;

/** Pointer to function, used to set the behavior of simulation when
    ns solver failed.  If equal to null, use DefaultCheckSolverOutput
    else (set with setCheckSolverFunction) call the pointer below).
    Note FP: (temporary) bad method to set checkSolverOutput but it
    works ... It may be better to use plug-in?
 */
static CheckSolverFPtr checkSolverOutput = NULL;

TimeStepping::TimeStepping(SP::TimeDiscretisation td): Simulation(td, "TimeStepping")
{}

// --- XML constructor ---
TimeStepping::TimeStepping(SP::SimulationXML strxml, double t0, double T, SP::DynamicalSystemsSet dsList, SP::InteractionsSet interList): Simulation(strxml, t0, T, dsList, interList, "TimeStepping")
{
  // === One Step NS Problem === For time stepping, only one non
  // smooth problem is built.
  if (simulationxml->hasOneStepNSProblemXML())  // ie if OSNSList is
    // not empty
  {
    SetOfOSNSPBXML OSNSList = simulationxml->getOneStepNSProblemsXML();
    if (OSNSList.size() != 1)
      RuntimeException::selfThrow("TimeStepping::xml constructor - Two many inputs for OSNS problems (only one problem is required).");
    SP::OneStepNSProblemXML osnsXML = *(OSNSList.begin());
    // OneStepNSProblem - Memory allocation/construction
    string type = osnsXML->getNSProblemType();
    if (type == LCP_TAG)  // LCP
    {
      (*allNSProblems)["timeStepping"].reset(new LCP(osnsXML));
    }
    else if (type == FRICTIONCONTACT_TAG)
    {
      (*allNSProblems)["timeStepping"].reset(new FrictionContact(osnsXML));
    }
    else RuntimeException::selfThrow("TimeStepping::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");

    (*allNSProblems)["timeStepping"]->setId("timeStepping");

    // Add QP and Relay cases when these classes will be fully
    // implemented.
  }
}

// --- Destructor ---
TimeStepping::~TimeStepping()
{
}

void TimeStepping::updateIndexSet(unsigned int i)
{
  // To update IndexSet number i: add or remove UnitaryRelations from
  // this set, depending on y values.

  if (i > indexSets.size())
    RuntimeException::selfThrow("TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");

  if (i == 0) // IndexSets[0] must not be updated in simulation, since it
    // belongs to the Topology.
    RuntimeException::selfThrow("TimeStepping::updateIndexSet(i=0), indexSets[0] can not be updated.");

  // for all Unitary Relations in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]
  UnitaryRelationsIterator it, itForFind;

  double y;

  // For all UR in Index[i-1] ...
  for (it = (indexSets[i - 1])->begin(); it != (indexSets[i - 1])->end(); ++it)
  {
    // itForFind: indicator to check if current Unitary Relation (ie
    // *it) is already in indexSets[1] (if not itForFind will be
    // equal to indexSets.end())
    itForFind = (indexSets[i])->find(*it);

    // Get y values for the considered UnitaryRelation
    y = (*it)->getYRef(i - 1); // y[i-1](0)
    if (i == 1)
    {
      double h = getTimeStep(); // Current time step
      double yDot = (*it)->getYRef(1); // y[1](0)
      y += 0.5 * h * yDot; // y_"prediction" = y[0](0) + 0.5*h*y[1](0)
    }

    // if y <=0, then the unitary relation is added in indexSets[1]
    // (if it was not already there) else if y > 0 and if the
    // unitary relation was in the set, it is removed.
    if (y <= 0 && itForFind == (indexSets[i])->end())
      (indexSets[i])->insert(*it);

    else if (y > 0 && itForFind != (indexSets[i])->end())
    {
      (indexSets[i])->erase(*it);
      (*it)->getLambdaPtr(i)->zero();
    }
  }
}

void TimeStepping::recordNonSmoothProblem(SP::OneStepNSProblem osns)
{
  // A the time, a time stepping simulation can only have one non
  // smooth problem.
  if (!allNSProblems->empty())
    RuntimeException::selfThrow("TimeStepping,  recordNonSmoothProblem - A non smooth problem already exist. You can not have more than one.");
  string name = "timeStepping";
  osns->setId(name);
  (*allNSProblems)[name] = osns;
}

void TimeStepping::initOSNS()
{
  // === creates links between work vector in OSI and work vector in
  // Unitary Relations
  SP::OneStepIntegrator  osi;

  ConstDSIterator itDS;
  UnitaryRelationsIterator it;
  // For each Unitary relation in I0 ...
  for (it = indexSets[0]->begin(); it != indexSets[0]->end(); ++it)
  {
    (*it)->initialize("TimeStepping");
    // creates a POINTER link between workX[ds] (xfree) and the
    // corresponding unitaryBlock in each UR for each ds of the
    // current UR.
    for (itDS = (*it)->dynamicalSystemsBegin();
         itDS != (*it)->dynamicalSystemsEnd(); ++itDS)
    {
      osi = osiMap[*itDS];
      (*it)->insertInWorkX(osi->getWorkX(*itDS));
    }
  }


  if (!allNSProblems->empty()) // ie if some Interactions have been
    // declared and a Non smooth problem
    // built.
  {
    if (allNSProblems->size() > 1)
      RuntimeException::selfThrow("TimeStepping::initialize, at the time, a time stepping simulation can not have more than one non smooth problem.");

    // At the time, we consider that for all systems, levelMin is
    // equal to the minimum value of the relative degree - 1 except
    // for degree 0 case where we keep 0.

    levelMin = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getMinRelativeDegree();

    if (levelMin != 0)
      levelMin--;

    // === update all index sets ===
    updateIndexSets();

    // initialization of  OneStepNonSmoothProblem
    for (OSNSIterator itOsns = allNSProblems->begin(); itOsns != allNSProblems->end(); ++itOsns)
    {
      (itOsns->second)->setLevels(levelMin, levelMax);
      (itOsns->second)->initialize(shared_from_this());
    }
  }
}

void TimeStepping::initLevelMax()
{
  levelMax = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getMaxRelativeDegree();
  // Interactions initialization (here, since level depends on the
  // type of simulation) level corresponds to the number of Y and
  // Lambda derivatives computed.

  if (levelMax != 0) // level max is equal to relative degree-1. But for
    // relative degree 0 case, we keep 0 value for
    // levelMax
    levelMax--;
}

void TimeStepping::nextStep()
{
  eventsManager->processEvents();
}


void TimeStepping::update(unsigned int levelInput)
{
  // 1 - compute input (lambda -> r)
  updateInput(levelInput);

  // 2 - compute state for each dynamical system

  OSIIterator itOSI;
  for (itOSI = allOSI->begin(); itOSI != allOSI->end() ; ++itOSI)
    (*itOSI)->updateState(levelInput);

  // 3 - compute output ( x ... -> y)
  if (!allNSProblems->empty())
  {
    updateOutput(0, levelMax);
  }
}

void TimeStepping::computeFreeState()
{
  std::for_each(allOSI->begin(), allOSI->end(), boost::bind(&OneStepIntegrator::computeFreeState, _1));
}

// compute simulation between current and next event.  Initial
// DS/interaction state is given by memory vectors and final state is
// the one saved in DS/Interaction at the end of this function
void TimeStepping::computeOneStep()
{
  advanceToEvent();
}

void TimeStepping::advanceToEvent()
{
  // solve ...
  computeFreeState();
  int info = 0;
  if (!allNSProblems->empty())
    info = computeOneStepNSProblem("timeStepping");
  // Check output from solver (convergence or not ...)
  if (!checkSolverOutput)
    DefaultCheckSolverOutput(info);
  else
    checkSolverOutput(info, this);
  // Update
  update(levelMin);
}

void TimeStepping::newtonSolve(double criterion, unsigned int maxStep)
{
  bool isNewtonConverge = false;
  unsigned int nbNewtonStep = 0; // number of Newton iterations
  //double residu = 0;
  int info = 0;
  while ((!isNewtonConverge) && (nbNewtonStep <= maxStep))
  {
    nbNewtonStep++;
    computeFreeState();
    if (!allNSProblems->empty())
      info = computeOneStepNSProblem("timeStepping");
    // Check output from solver (convergence or not ...)
    if (!checkSolverOutput)
      DefaultCheckSolverOutput(info);
    else
      checkSolverOutput(info, this);

    update(levelMin);
    isNewtonConverge = newtonCheckConvergence(criterion);
  }

  if (!isNewtonConverge)
    cout << "Newton process stopped: max. steps number reached." << endl ;
}

bool TimeStepping::newtonCheckConvergence(double criterion)
{
  bool checkConvergence = false;
  // get the nsds indicator of convergence
  double nsdsConverge = model-> getNonSmoothDynamicalSystemPtr()->nsdsConvergenceIndicator();
  // We compute cvg = abs(xi+1 - xi)/xi and if cvg < criterion
  if (nsdsConverge < criterion)
  {
    checkConvergence = true ;
    double residu = 0.0;
    for (OSIIterator it = allOSI->begin(); it != allOSI->end() ; ++it)
    {
      residu = (*it)->computeResidu();
      //    cout << residu << endl;
      if (residu > criterion)
      {
        checkConvergence = false;
        break;
      }
    }
  }

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
      advanceToEvent();

    else if (opt == "Newton")
      newtonSolve(criterion, maxIter);

    else
      RuntimeException::selfThrow("TimeStepping::run(opt) failed. Unknow simulation option: " + opt);

    eventsManager->processEvents();
    count++;
  }
  cout << "===== End of " << simulationType << "simulation. " << count << " events have been processed. ==== " << endl;
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
