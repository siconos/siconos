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
#include "NewtonEulerDS.hpp"
#include "NewtonEulerFrom1DLocalFrameR.hpp"
using namespace std;
//#define TSPROJ_DEBUG
//#define CORRECTIONSVELOCITIES
TimeSteppingProjectOnConstraints::TimeSteppingProjectOnConstraints(SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo,
    SP::OneStepNSProblem osnspb_pos,
    unsigned int level)
  : TimeStepping(td, osi, osnspb_velo)
{
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, SICONOS_OSNSP_TS_POS);

  _indexSetLevelForProjection = level;
  _constraintTol = 1e-04;
  _constraintTolUnilateral = 1e-08;
  _projectionMaxIteration = 10;
  _doProj = 1;
  _doOnlyProj = 0;


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
  unsigned int cmp = 0;
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

  while (runningProjection && cmp < _projectionMaxIteration)
  {
    cmp++;
#ifdef TSPROJ_DEBUG
    printf("TimeSteppingProjectOnConstraints projection step = %d\n", cmp);
#endif
    info = 0;
#ifdef TSPROJ_DEBUG
    cout << "TimeSteppingProjectOnConstraint compute OSNSP." << endl ;
#endif
    info = computeOneStepNSProblem(SICONOS_OSNSP_TS_POS);

    if (info)
    {
      cout << "TimeSteppingProjectOnConstraints1 project on constraints failed." << endl ;
      return;
    }

    for (DynamicalSystemsGraph::VIterator aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2)
    {
      SP::DynamicalSystem ds = dsGraph->bundle(*aVi2);
      Type::Siconos dsType = Type::value(*ds);
      if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS neds = boost::static_pointer_cast<NewtonEulerDS>(ds);
        neds->normalizeq();
        neds->updateT();
      }
      else if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
      {

      }
      else
        RuntimeException::selfThrow("TimeSteppingProjectOnConstraints :: - Ds is not from NewtonEulerDS neither from LagrangianDS.");

    }

    updateWorldFromDS();

    computeCriteria(&runningProjection);

    //cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  Z:"<<endl;
    //(*_allNSProblems)[SICONOS_OSNSP_TS_POS]->display();
    //(boost::static_pointer_cast<LinearOSNS>((*_allNSProblems)[SICONOS_OSNSP_TS_POS]))->z()->display();
    if (info)
      cout << "TimeSteppingProjectOnConstraints2 project on constraints failed." << endl ;
#ifdef TSPROJ_DEBUG
    cout << "TimeSteppingProjectOnConstraints::Projection end : Number of iterations=" << cmp << "\n";
#endif

    //cout<<"during projection before normalizing of q:\n";
    //for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    //{
    //  (*it)->relation()->computeh(getTkp1());
    //}
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

  //       SP::SimpleVector FextNorm(new SimpleVector(3));
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
  double minViolationUnilateral = +1e24;

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
      double criteria = interac->y(0)->getValue(0);
#ifdef TSPROJ_DEBUG
      printf("unilatreal interac->y(0)->getValue(0) %e.\n", interac->y(0)->getValue(0));
#endif

      if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;
      if (criteria < minViolationUnilateral) minViolationUnilateral = criteria;
      if (criteria < - _constraintTolUnilateral)
      {
        *runningProjection = true;
#ifdef TSPROJ_DEBUG
        printf("TSProj newton criteria unilateral true %e.\n", criteria);
#endif
      }
    }
    else
    {
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

  }
#ifdef TSPROJ_DEBUG
  printf("TSProj newton min/max criteria projection\n");
  std::cout << "                 runningProjection"  << *runningProjection << std::endl;
  printf("              min criteria equality =  %e.\n", minViolationEquality);
  printf("              max criteria equality =  %e.\n", maxViolationEquality);
  printf("              max criteria unilateral =  %e.\n", maxViolationUnilateral);
  printf("              min criteria unilateral =  %e.\n", minViolationUnilateral);
#endif
}
