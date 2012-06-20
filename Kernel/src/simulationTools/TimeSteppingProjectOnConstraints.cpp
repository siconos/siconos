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

#include "TimeSteppingProjectOnConstraints.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerFrom1DLocalFrameR.hpp"
#include "OneStepIntegrator.hpp"
using namespace std;
static CheckSolverFPtr checkSolverOutputProjectOnConstraints = NULL;
//#define TSPROJ_DEBUG
//#define CORRECTIONSVELOCITIES
TimeSteppingProjectOnConstraints::TimeSteppingProjectOnConstraints(SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo,
    SP::OneStepNSProblem osnspb_pos,
    unsigned int level)
  : TimeStepping(td, osi, osnspb_velo)
{

  //if (Type::value(osi) != Type::MoreauProjectOnConstraintsOSI)
  OSI::TYPES typeOSI;
  typeOSI = (osi)->getType();
  if (typeOSI != OSI::MOREAUPROJECTONCONSTRAINTSOSI)
    RuntimeException::selfThrow("TimeSteppingProjectOnConstraints::TimeSteppingProjectOnConstraints.  wrong type of OneStepIntegrator");

  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, SICONOS_OSNSP_TS_POS);

  _indexSetLevelForProjection = level;
  _constraintTol = 1e-04;
  _constraintTolUnilateral = 1e-08;
  _projectionMaxIteration = 10;
  _doProj = 1;
  _doOnlyProj = 0;
  _maxViolationUnilateral = 0.0;
  _maxViolationEquality = 0.0;

}

// --- Destructor ---
TimeSteppingProjectOnConstraints::~TimeSteppingProjectOnConstraints()
{
}

void TimeSteppingProjectOnConstraints::initOSNS()
{
  TimeStepping::initOSNS();

  (*_allNSProblems)[SICONOS_OSNSP_TS_POS]->setLevelMin(_indexSetLevelForProjection);
  (*_allNSProblems)[SICONOS_OSNSP_TS_POS]->setLevelMax(_indexSetLevelForProjection);

  (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->setLevelMin(_levelMaxForInput);
  (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->setLevelMax(_levelMaxForInput);
}

void TimeSteppingProjectOnConstraints::nextStep()
{
  TimeStepping::nextStep();


  // Zeroing Lambda Muliplier of indexSet()

  SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    inter->lambda(0)->zero();
  }

}

void TimeSteppingProjectOnConstraints::advanceToEvent()
{
  /** First step, Solve the standard velocity formulation.*/
#ifdef TSPROJ_DEBUG
  cout << "TimeStepping::newtonSolve begin :\n";
#endif
  if (!_doOnlyProj)
    TimeStepping::newtonSolve(_newtonTolerance, _newtonMaxIteration);
#ifdef TSPROJ_DEBUG
  cout << "TimeStepping::newtonSolve end : Number of iterations=" << getNewtonNbSteps() << "\n";
  cout << "                              : newtonResiduDSMax=" << newtonResiduDSMax() << "\n";
  cout << "                              : newtonResiduYMax=" << newtonResiduYMax() << "\n";
  cout << "                              : newtonResiduRMax=" << newtonResiduRMax() << "\n";
#endif

  if (!_doProj)
    return;
  int info = 0;




  /** Second step, Perform the projection on constraints.*/
#ifdef TSPROJ_DEBUG
  cout << "TimeSteppingProjectOnConstraints::newtonSolve begin projection:\n";
#endif
  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();


#ifdef TSPROJ_CORRECTIONVELOCITIES
  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    SP::DynamicalSystem ds = dsGraph->bundle(*vi);
    Type::Siconos dsType = Type::value(*ds);
    if (dsType != Type::NewtonEulerDS)
      RuntimeException::selfThrow("TS:: - ds is not from NewtonEulerDS.");
    SP::NewtonEulerDS neds = boost::static_pointer_cast<NewtonEulerDS>(ds);
    *(neds->deltaq()) = *(neds->q());
  }
#endif

  bool runningProjection = false;
  _nbProjectionIteration = 0;
  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
  // for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++){
  //   double criteria = (*it)->relation()->y(0)->getValue(0);
  //   if (Type::value(*((*it)->nonSmoothLaw())) ==  Type::NewtonImpactFrictionNSL ||
  //  Type::value(*((*it)->nonSmoothLaw())) == Type::NewtonImpactNSL){
  //     SP::NewtonEulerFrom1DLocalFrameR ri = boost::static_pointer_cast<NewtonEulerFrom1DLocalFrameR> ((*it)->relation());
  //     if (criteria < -1e-7){
  //  ri->_isOnContact=true;
  //     }else{
  //  ri->_isOnContact=false;
  //     }
  //   }
  //   if (criteria < -_constraintTol)
  //     runningNewton=true;
  // }
  if (model()->nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() > _indexSetLevelForProjection)
    computeCriteria(&runningProjection);
  // Zeroing Lambda Muliplier of indexSet()

  SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    inter->lambda(0)->zero();
  }
  updateInput(0);

  //Store the q vector of each DS.

  for (DynamicalSystemsGraph::VIterator aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2)
  {
    SP::DynamicalSystem ds = dsGraph->bundle(*aVi2);
    Type::Siconos dsType = Type::value(*ds);
    if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS neds = boost::static_pointer_cast<NewtonEulerDS>(ds);
      neds->addWorkVector(neds->q(), DynamicalSystem::qtmp);
    }
    else if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
      d->addWorkVector(d->q(), DynamicalSystem::qtmp);
    }
    else
      RuntimeException::selfThrow("TimeSteppingProjectOnConstraints::advanceToEvent() :: - Ds is not from NewtonEulerDS neither from LagrangianDS.");
  }

  while (runningProjection && _nbProjectionIteration < _projectionMaxIteration)
  {
    _nbProjectionIteration++;
#ifdef TSPROJ_DEBUG
    printf("TimeSteppingProjectOnConstraints projection step = %d\n", _nbProjectionIteration);
#endif
    SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
    InteractionsGraph::VIterator ui, uiend;
    for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet->bundle(*ui);
      inter->lambda(0)->zero();
    }
    updateInput(0);
    info = 0;
#ifdef TSPROJ_DEBUG
    cout << "TimeSteppingProjectOnConstraint compute OSNSP." << endl ;
#endif
    info = computeOneStepNSProblem(SICONOS_OSNSP_TS_POS);

    if (info)
    {
      cout << " TimeSteppingProjectOnConstraints::advanceToEvent() project on constraints. solver failed." << endl ;
      return;
    }
    updateInput(0);
#ifdef TSPROJ_DEBUG

    std ::cout << "After update input" << std::endl;
    SP::InteractionsGraph indexSet1 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    std ::cout << "lamda(1) in IndexSet1" << std::endl;
    for (boost::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet1->bundle(*ui);
      inter->lambda(1)->display();
    }
    SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
    std ::cout << "lamda(0) in indexSet0" << std::endl;
    for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet0->bundle(*ui);
      inter->lambda(0)->display();
    }

#endif
    // This part should be in MoreauProjectOnConstraintsOS::updateState(level =0)
    for (DynamicalSystemsGraph::VIterator aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2)
    {
      SP::DynamicalSystem ds = dsGraph->bundle(*aVi2);
      Type::Siconos dsType = Type::value(*ds);
      if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS neds = boost::static_pointer_cast<NewtonEulerDS>(ds);
        SP::SiconosVector q = neds->q();
        SP::SiconosVector qtmp = neds->getWorkVector(DynamicalSystem::qtmp);
#ifdef TSPROJ_DEBUG
        std ::cout << "qtmp before  update " << std::endl;
        qtmp->display();
        std ::cout << "p(0) before  update " << std::endl;
        neds->p(0)->display();

#endif
        if (neds->p(0))
        {
          //*q = * qtmp +  *neds->p(0);
          *q += *neds->p(0); // Why it works like that and not with the previous line ?
        }
#ifdef TSPROJ_DEBUG
        std ::cout << "q after  update " << std::endl;
        q->display();
#endif



        neds->normalizeq();
        neds->updateT();
      }
      else if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
      {
        SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
        SP::SiconosVector q = d->q();
        SP::SiconosVector qtmp = d->getWorkVector(DynamicalSystem::qtmp);

        if (d->p(0))
        {
          //*q = * qtmp +  *d->p(0);
          *q += *d->p(0);
        }
      }
      else
        RuntimeException::selfThrow("TimeSteppingProjectOnConstraints::advanceToEvent() :: - Ds is not from NewtonEulerDS neither from LagrangianDS.");

    }

    updateWorldFromDS();

    computeCriteria(&runningProjection);

    //cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  Z:"<<endl;
    //(*_allNSProblems)[SICONOS_OSNSP_TS_POS]->display();
    //(boost::static_pointer_cast<LinearOSNS>((*_allNSProblems)[SICONOS_OSNSP_TS_POS]))->z()->display();

#ifdef TSPROJ_DEBUG
    cout << "TimeSteppingProjectOnConstraints::Projection end : Number of iterations=" << _nbProjectionIteration << "\n";
    std ::cout << "After update state in position" << std::endl;
    std ::cout << "lamda(1) in IndexSet1" << std::endl;
    for (boost::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet1->bundle(*ui);
      inter->lambda(1)->display();
    }
    std ::cout << "lamda(0) in indexSet0" << std::endl;
    for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet0->bundle(*ui);
      inter->lambda(0)->display();
    }
    std ::cout << "y(1) in IndexSet1" << std::endl;
    for (boost::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet1->bundle(*ui);
      inter->y(1)->display();
    }
    std ::cout << "y(0) in indexSet0" << std::endl;
    for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet0->bundle(*ui);
      inter->y(0)->display();
    }

#endif

    //cout<<"during projection before normalizing of q:\n";
    //for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    //{
    //  (*it)->relation()->computeh(getTkp1());
    //}
  }// end while(runningProjection && _nbProjectionIteration < _projectionMaxIteration)
  if (_nbProjectionIteration == _projectionMaxIteration)
  {
    cout << "TimeSteppingProjectOnConstraints::advanceToEvent() Max number of projection iterations reached (" << _nbProjectionIteration << ")"  << endl ;
    printf("              max criteria equality =  %e.\n", _maxViolationEquality);
    printf("              max criteria unilateral =  %e.\n", _maxViolationUnilateral);
    RuntimeException::selfThrow("youyou");
  }



#ifdef TSPROJ_DEBUG
  cout << "TimeSteppingProjectOnConstraints::newtonSolve end projection:\n";
#endif

  return;
  //#ifdef TSPROJ_CORRECTIONVELOCITIES
  //   /*The following reduces the velocity because the position step increase the energy of the system. This formulation works only with simple systems.To activate it, comment the next line.*/

  //   for(DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  //   {
  //     SP::DynamicalSystem ds = dsGraph->bundle(*vi);
  //     Type::Siconos dsType = Type::value(*ds);
  //     if(dsType !=Type::NewtonEulerDS)
  //       RuntimeException::selfThrow("TS:: - ds is not from NewtonEulerDS.");
  //     // SP::SiconosVector dotq = neds->dotq();
  //     // SP::SiconosVector q = neds->q();
  //     // SP::SiconosVector qold = neds->qMemory()->getSiconosVector(0);
  //     // double h = timeStep();
  //     // dotq->setValue(3,(q->getValue(3)-qold->getValue(3))/h);
  //     // dotq->setValue(4,(q->getValue(4)-qold->getValue(4))/h);
  //     // dotq->setValue(5,(q->getValue(5)-qold->getValue(5))/h);
  //     // dotq->setValue(6,(q->getValue(6)-qold->getValue(6))/h);

  //     /*compute the new velocity seeing the work of fext*/
  //     SP::NewtonEulerDS neds = boost::static_pointer_cast<NewtonEulerDS>(ds);
  //     *(neds->deltaq())-=*(neds->q());
  // #ifdef TSPROJ_DEBUG
  //     printf("TSProj NewtonSolve :deltaq:");
  //     (neds->deltaq())->display();
  // #endif
  //     //continue;
  //     double  n2q=neds->deltaq()->norm2();
  //     double n2=0.0;
  //     if(neds->fExt())
  //       n2=neds->fExt()->norm2();
  //     if(n2 > 1e-7 && n2q > 1e-14)
  //     {
  //       //if (n2q < 1e-14)
  //       // continue;

  //       SP::SiconosVector FextNorm(new SiconosVector(3));
  //       FextNorm->setValue(0,neds->fExt()->getValue(0));
  //       FextNorm->setValue(1,neds->fExt()->getValue(1));
  //       FextNorm->setValue(2,neds->fExt()->getValue(2));
  // #ifdef TSPROJ_DEBUG
  //       cout<<"TimeSteppingProjectOnConstraints::newtonSolve deltaQ :\n";
  //       neds->deltaq()->display();
  //       cout<<"TimeSteppingProjectOnConstraints::newtonSolve Fext :\n";
  //       FextNorm->display();
  // #endif

  //       (*FextNorm)*=(1./n2);
  //       /*work of external forces.*/
  //       double workFext= neds->fExt()->getValue(0) * neds->deltaq()->getValue(0)+
  //                        neds->fExt()->getValue(1) * neds->deltaq()->getValue(1)+
  //                        neds->fExt()->getValue(2) * neds->deltaq()->getValue(2);
  //       //workFext*=2.0;
  //       double VkFNorm=FextNorm->getValue(0)*neds->velocity()->getValue(0)+
  //                      FextNorm->getValue(1)*neds->velocity()->getValue(1)+
  //                      FextNorm->getValue(2)*neds->velocity()->getValue(2);
  //       double VkcFNorm =VkFNorm;
  //       VkcFNorm= VkFNorm*VkFNorm - 2*fabs(workFext)/(neds->massValue());
  //       if(VkcFNorm >0)
  //       {
  //         if(VkFNorm>0)
  //           VkcFNorm=sqrt(VkcFNorm);
  //         else
  //           VkcFNorm=-sqrt(VkcFNorm);
  //       }
  //       else
  //         VkcFNorm=0;
  //       // if (VkFNorm >= 0 && workFext >0){
  //       //   ;//VkcFNorm=sqrt (2*workFext/(neds->massValue())+VkFNorm*VkFNorm);
  //       // }else if (VkFNorm <= 0 && workFext < 0){
  //       //   ;//VkcFNorm=-sqrt (fabs(2*workFext/(neds->massValue())+VkFNorm*VkFNorm));
  //       // }else if (VkFNorm > 0 && workFext <0){
  //       //   VkcFNorm= VkFNorm*VkFNorm + 2*workFext/(neds->massValue());
  //       //   if (VkcFNorm >0)
  //       //     VkcFNorm=sqrt(VkcFNorm);
  //       //   else
  //       //     VkcFNorm=0;
  //       // }else if (VkFNorm < 0 && workFext > 0){
  //       //   VkcFNorm= VkFNorm*VkFNorm - 2*workFext/(neds->massValue());
  //       //   if (VkcFNorm >0)
  //       //     VkcFNorm=-sqrt(VkcFNorm);
  //       //   else
  //       //     VkcFNorm=0;
  //       // }
  // #ifdef TSPROJ_DEBUG
  //       printf("TimeSteppingProjectOnConstraints::newtonSolve velocity before update(prevComp=%e, newComp=%e)\n",VkFNorm,VkcFNorm);
  //       printf("VELOCITY1 ");
  //       neds->velocity()->display();
  // #endif
  //       neds->velocity()->setValue(0,neds->velocity()->getValue(0)+(VkcFNorm - VkFNorm)*FextNorm->getValue(0));
  //       neds->velocity()->setValue(1,neds->velocity()->getValue(1)+(VkcFNorm - VkFNorm)*FextNorm->getValue(1));
  //       neds->velocity()->setValue(2,neds->velocity()->getValue(2)+(VkcFNorm - VkFNorm)*FextNorm->getValue(2));
  // #ifdef TSPROJ_DEBUG
  //       cout<<"TimeSteppingProjectOnConstraints::newtonSolve velocity updated\n";
  //       printf("VELOCITY2 ");
  //       neds->velocity()->display();
  // #endif
  //     }
  //     SP::SiconosMatrix T = neds->T();
  //     SP::SiconosVector dotq = neds->dotq();
  //     prod(*T,*neds->velocity(),*dotq,true);
  //     if(!_allNSProblems->empty())
  //     {
  //       for(unsigned int level = _levelMinForOutput;
  //           level < _levelMaxForOutput;
  //           level++)
  //         updateOutput(level);
  //     }
  //   }
  //#endif
#ifdef TSPROJ_DEBUG
  cout << "TimeSteppingProjectOnConstraints::newtonSolve end projection:\n";
#endif

}

void TimeSteppingProjectOnConstraints::computeCriteria(bool * runningProjection)
{

  SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevelForProjection);
  InteractionsGraph::VIterator aVi, viend;

  double maxViolationEquality = -1e24;
  double minViolationEquality = +1e24;
  double maxViolationUnilateral = -1e24;
  //double minViolationUnilateral = +1e24;

  *runningProjection = false;

  for (boost::tie(aVi, viend) = indexSet->vertices();
       aVi != viend; ++aVi)
  {
    SP::Interaction inter = indexSet->bundle(*aVi);
    SP::Interaction interac = inter;
    interac->relation()->computeOutput(getTkp1(), 0);
    interac->relation()->computeJach(getTkp1());
    if (Type::value(*(interac->nonSmoothLaw())) ==  Type::NewtonImpactFrictionNSL ||
        Type::value(*(interac->nonSmoothLaw())) == Type::NewtonImpactNSL)
    {
      double criteria = std::max(0.0, - interac->y(0)->getValue(0));
#ifdef TSPROJ_DEBUG
      printf("Unilateral interac->y(0)->getValue(0) %e.\n", interac->y(0)->getValue(0));
#endif

      if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;
      //if (criteria < minViolationUnilateral) minViolationUnilateral=criteria;
      if (maxViolationUnilateral > _constraintTolUnilateral)
      {
        *runningProjection = true;
#ifdef TSPROJ_DEBUG
        printf("TSProj newton criteria unilateral true %e.\n", criteria);
#endif
      }
    }
    else
    {
#ifdef TSPROJ_DEBUG
      printf("Equality interac->y(0)->normInf() %e.\n", interac->y(0)->normInf());
#endif
      if (interac->y(0)->normInf()  > maxViolationEquality) maxViolationEquality = interac->y(0)->normInf() ;
      if (interac->y(0)->normInf()  < minViolationEquality) minViolationEquality = interac->y(0)->normInf() ;
      if (interac->y(0)->normInf() > _constraintTol)
      {
        *runningProjection = true;
#ifdef TSPROJ_DEBUG
        printf("TSProj  newton criteria equality true %e.\n", interac->y(0)->normInf());
#endif
      }
    }
    _maxViolationUnilateral = maxViolationUnilateral;
    _maxViolationEquality = maxViolationEquality;
  }
#ifdef TSPROJ_DEBUG
  printf("TSProj newton min/max criteria projection\n");
  std::cout << "             runningProjection "  << std::boolalpha << *runningProjection << std::endl;
  printf("              min criteria equality =  %e.\n", minViolationEquality);
  printf("              max criteria equality =  %e.\n", maxViolationEquality);
  printf("              max criteria unilateral =  %e.\n", maxViolationUnilateral);
  //printf("              min criteria unilateral =  %e.\n",minViolationUnilateral);
#endif
}



void TimeSteppingProjectOnConstraints::newtonSolve(double criterion, unsigned int maxStep)
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
    // updateOutput(0);
    // updateIndexSets();
    if (!_allNSProblems->empty() &&  !allInteractions->isEmpty())
      info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
    // Check output from solver (convergence or not ...)
    if (!checkSolverOutputProjectOnConstraints)
      DefaultCheckSolverOutput(info);
    else
      checkSolverOutputProjectOnConstraints(info, this);

    update(_levelMaxForInput);

    //isNewtonConverge = newtonCheckConvergence(criterion);
    if (!_allNSProblems->empty() &&  !allInteractions->isEmpty())
      saveYandLambdaInMemory();
  }

  else if (_newtonOptions == SICONOS_TS_NONLINEAR)
  {
    while ((!isNewtonConverge) && (_newtonNbSteps < maxStep) && (!info))
    {
      _newtonNbSteps++;
      prepareNewtonIteration();
      computeFreeState();
      // updateOutput(0);
      // updateIndexSets();
      if (info)
        cout << "new loop because of info\n" << endl;

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
      if (info)
        cout << "info!" << endl;
      // Check output from solver (convergence or not ...)
      if (!checkSolverOutputProjectOnConstraints)
        DefaultCheckSolverOutput(info);
      else
        checkSolverOutputProjectOnConstraints(info, this);

      update(_levelMaxForInput);
      isNewtonConverge = newtonCheckConvergence(criterion);
      if (!isNewtonConverge && !info)
      {
        if (!_allNSProblems->empty() &&  !allInteractions->isEmpty())
          saveYandLambdaInMemory();
      }
    }
    if (!isNewtonConverge)
      cout << "TimeStepping::newtonSolve -- Newton process stopped: max. number of steps (" << maxStep << ") reached." << endl ;
    else if (info)
      cout << "TimeStepping::newtonSolve -- Newton process stopped: solver failed." << endl ;
    //    else
    //      cout << "TimeStepping::newtonSolve succed nbit="<<_newtonNbSteps<<"maxStep="<<maxStep<<endl;
  }
  else
    RuntimeException::selfThrow("TimeStepping::NewtonSolve failed. Unknow newtonOptions: " + _newtonOptions);
}
