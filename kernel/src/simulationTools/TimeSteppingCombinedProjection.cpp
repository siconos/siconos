/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


#include "TimeStepping.hpp"

#include "TimeSteppingCombinedProjection.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerFrom1DLocalFrameR.hpp"
#include "OneStepIntegrator.hpp"
#include "MLCPProjectOnConstraints.hpp"
#include "NonSmoothLaw.hpp"
#include "NewtonEulerR.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"




//#define TSPROJ_DEBUG_LEVEL1
//#define TSPROJ_WITHOUT_PROJECTION
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>




TimeSteppingCombinedProjection::TimeSteppingCombinedProjection(SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo,
    SP::OneStepNSProblem osnspb_pos,
    unsigned int level)
  : TimeStepping(td, osi, osnspb_velo)
{
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, SICONOS_OSNSP_TS_POS);

  _indexSetLevelForProjection = level;
  if (_indexSetLevelForProjection != 2)
  {
    RuntimeException::selfThrow("TimeSteppingCombinedProjection::TimeSteppingCombinedProjection level not equal to 2 is not yet implemented.  ");

  }
  _constraintTol = 1e-08;
  _constraintTolUnilateral = 1e-08;
  _projectionMaxIteration = 50;
  _kIndexSetMax = 50;
  _doCombinedProj = true;
  _doCombinedProjOnEquality = true;
  _isIndexSetsStable = false;
  _maxViolationUnilateral = 0.0;
  _maxViolationEquality = 0.0;
}

// --- Destructor ---
TimeSteppingCombinedProjection::~TimeSteppingCombinedProjection()
{
}



struct TimeSteppingCombinedProjection::_SimulationEffectOnOSNSP : public SiconosVisitor
{

  using SiconosVisitor::visit;

  TimeSteppingCombinedProjection * _parent;

  _SimulationEffectOnOSNSP(TimeSteppingCombinedProjection * p) :
    _parent(p)
  {
    std::cout << "hello" << std::endl;
  };

  void visit(MLCPProjectOnConstraints &   onsnsp)
  {

    bool toto = (bool)_parent->doCombinedProjOnEquality();
    onsnsp.setDoProjOnEquality(toto);
  }
  void visit(MLCPProjectOnConstraints *   onsnsp)
  {
    std::cout << "hello" << std::endl;
  }
  void visit(MLCPProjectOnConstraints   onsnsp)
  {
    std::cout << "hello" << std::endl;
  }

};


void TimeSteppingCombinedProjection::initOSNS()
{
  TimeStepping::initOSNS();


  SP::OneStepNSProblem osnspb_pos = (*_allNSProblems)[SICONOS_OSNSP_TS_POS];

  osnspb_pos->setIndexSetLevel(_indexSetLevelForProjection);
  osnspb_pos->setInputOutputLevel(0);

  (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
  (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);


  // better with visitor but I am not able to fix it.
  //_SimulationEffectOnOSNSP simulationEffectOnOSNSP(this);
  //osnspb_pos->accept(simulationEffectOnOSNSP);
  if (Type::value(*osnspb_pos) ==     Type::MLCPProjectOnConstraints)
  {
    // std::cout << "Type::name(*osnspb_pos) "<< Type::name(*osnspb_pos) <<std::endl;
    MLCPProjectOnConstraints * toto = static_cast<MLCPProjectOnConstraints*>(osnspb_pos.get());
    bool tutu = (bool)_doCombinedProjOnEquality;
    toto -> setDoProjOnEquality(tutu);
  }

}



void TimeSteppingCombinedProjection::advanceToEventOLD()
{


  newtonSolve(_newtonTolerance, _newtonMaxIteration);

  assert(0);
}




void TimeSteppingCombinedProjection::advanceToEvent()
{
  DEBUG_PRINT("================================================");
  DEBUG_PRINT("TimeSteppingCombinedProjection::advanceToEvent()");
  DEBUG_PRINT("================================================\n");
  _isIndexSetsStable = false;
  _maxViolationUnilateral = 0.0;
  _maxViolationEquality = 0.0;

  // Update interactions if a manager was provided
  updateInteractions();

  SP::Topology topo = _nsds->topology();
  if (topo->numberOfIndexSet() > _indexSetLevelForProjection)
  {
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);
    SP::InteractionsGraph indexSet2 = topo->indexSet(2);
    assert(indexSet1);
    assert(indexSet2);

    InteractionsGraph::VIterator ui, uiend, vnext;

    // zeroing the lambda of indexSet1
    std11::tie(ui, uiend) = indexSet1->vertices();
    for (vnext = ui; ui != uiend; ui = vnext)
    {
      ++vnext;
      SP::Interaction inter1 = indexSet1->bundle(*ui);
      inter1->lambda(1)->zero();
      indexSet1->eraseProperties(*ui);
      InteractionsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = indexSet1->out_edges(*ui);
           oei != oeiend; ++oei)
      {
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet1->edges(indexSet1->source(*oei), indexSet1->target(*oei));
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
    }

    indexSet1->clear();

    std11::tie(ui, uiend) = indexSet2->vertices();
    for (vnext = ui; ui != uiend; ui = vnext)
    {
      ++vnext;
      indexSet2->eraseProperties(*ui);
      InteractionsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = indexSet2->out_edges(*ui);
           oei != oeiend; ++oei)
      {
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet2->edges(indexSet2->source(*oei), indexSet2->target(*oei));
        if (ed2 != ed1)
        {
          indexSet2->eraseProperties(ed1);
          indexSet2->eraseProperties(ed2);
        }
        else
        {
          indexSet2->eraseProperties(ed1);
        }
      }
    }



    indexSet2->clear();
  }



  _nbIndexSetsIteration = 0 ;
  _cumulatedNewtonNbIterations = 0 ;
  _nbCumulatedProjectionIteration = 0;

  while (!_isIndexSetsStable)
  {
    _nbIndexSetsIteration++ ;
    InteractionsGraph::VIterator ui, uiend;
    DEBUG_PRINTF( "====================== _nbIndexSetsIteration = %i \n ", _nbIndexSetsIteration );

#ifdef TSPROJ_DEBUG_LEVEL1
    SP::InteractionsGraph indexSet0 = topo->indexSet(0);
    std::cout << "indexSet0->size() " << indexSet0->size()   <<std::endl;
    unsigned int level;
    SP::Interaction inter;
#endif

#ifdef TSPROJ_DEBUG_LEVEL2

    if (topo->numberOfIndexSet() > _indexSetLevelForProjection)
    {
      SP::InteractionsGraph indexSet1 = topo->indexSet(1);
      SP::InteractionsGraph indexSet2 = topo->indexSet(2);
      std::cout << "indexSet1->size() " << indexSet1->size()   <<std::endl;
      std::cout << "indexSet2->size() " << indexSet2->size()   <<std::endl;
    }


    for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      inter = indexSet0->bundle(*ui);

      std::cout << "inter->number()" << inter->number()   << std::endl;

      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 0);
      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 1);

      //  inter->swapInMemory();

      level = 0;

      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);

      //std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     << std::endl;
      std::cout << "inter->y(" << level << ")"   << std::endl;
      inter->y(level)->display();

      std::cout << "inter->y_k(" << level << ")"   << std::endl;
      inter->y_k(level)->display();

      level = 1;
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      //std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     << std::endl;
      std::cout << "inter->y(" << level << ")"   << std::endl;
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")"   << std::endl;
      inter->y_k(level)->display();
    }

#endif



    if (_nbIndexSetsIteration > _kIndexSetMax)
    {
      RuntimeException::selfThrow("TimeSteppingCombinedProjection::TimeSteppingCombinedProjection _nbIndexSetsIteration >  _kIndexSetMax ");
    }


    /** First step, Solve the standard velocity formulation.*/
    TimeStepping::newtonSolve(_newtonTolerance, _newtonMaxIteration);
    _cumulatedNewtonNbIterations += getNewtonNbIterations();

#ifdef TSPROJ_DEBUG_LEVEL1

  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    std::cout << "inter->number()" << inter->number()   << std::endl;

    level = 0;

    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);

    //   std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     << std::endl;
    std::cout << "inter->y(" << level << ")"   << std::endl;
    inter->y(level)->display();
    std::cout << "inter->y_k(" << level << ")"   << std::endl;
    inter->y_k(level)->display();


    level = 1;
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    //std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     << std::endl;
    std::cout << "inter->y(" << level << ")"   << std::endl;
    inter->y(level)->display();
    std::cout << "inter->y_k(" << level << ")"   << std::endl;
    inter->y_k(level)->display();
  }
#endif

    int info = 0;

    // Zeroing Lambda Muliplier of indexSet()

    SP::InteractionsGraph indexSet = _nsds->topology()->indexSet(0);
    for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet->bundle(*ui);
      inter->lambda(0)->zero();
    }
    _nsds->updateInput(nextTime(),0);

#ifdef TSPROJ_WITHOUT_PROJECTION

#else
    /** Second step, Perform the projection on constraints.*/
    DEBUG_PRINT( "TimeSteppingCombinedProjection::newtonSolve begin projection:\n");

    SP::DynamicalSystemsGraph dsGraph = _nsds->dynamicalSystems();


    bool runningProjection = false;
    _nbProjectionIteration = 0;

    if (_nsds->topology()->numberOfIndexSet() > _indexSetLevelForProjection)
    {
      updateIndexSet(2);
      computeCriteria(&runningProjection);

#ifdef TSPROJ_DEBUG_LEVEL1

      SP::InteractionsGraph indexSet2 = topo->indexSet(2);
      SP::InteractionsGraph indexSet1 = topo->indexSet(1);
      if (indexSet2->size() > 1)
      {
        printf("indexSet2->size() = %i >1 \n", (int)indexSet2->size());
      }
      if (indexSet1->size() > 0)
      {
        printf("indexSet1->size() = %i >0 \n", (int)indexSet1->size());
      }

#endif


    }


    if (!runningProjection)
    {
      // Zeroing Lambda Muliplier of indexSet()

      SP::InteractionsGraph indexSet = _nsds->topology()->indexSet(0);
      InteractionsGraph::VIterator ui, uiend;
      for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter = indexSet->bundle(*ui);
        inter->lambda(0)->zero();
      }
      _nsds->updateInput(nextTime(),0);
    }

    //Store the q vector of each DS.

    for (DynamicalSystemsGraph::VIterator aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2)
    {
      SP::DynamicalSystem ds = dsGraph->bundle(*aVi2);
      Type::Siconos dsType = Type::value(*ds);
      VectorOfVectors& workVectors = *dsGraph->properties(*aVi2).workVectors;

      if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS neds = std11::static_pointer_cast<NewtonEulerDS>(ds);
        *workVectors[OneStepIntegrator::qtmp] = *neds->q();
      }
      else if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
        *workVectors[OneStepIntegrator::qtmp] = * d->q();
      }
      else
        RuntimeException::selfThrow("TimeSteppingCombinedProjection::advanceToEvent() :: - Ds is not from NewtonEulerDS neither from LagrangianDS.");
    }



    _nbProjectionIteration = 0;

    while ((runningProjection && _nbProjectionIteration < _projectionMaxIteration) && _doCombinedProj)
    {
      _nbProjectionIteration++;

      DEBUG_PRINTF("Projection iteration number   %d\t", _nbProjectionIteration);
      DEBUG_PRINT("================================================\n");


      // Zeroing Lambda Muliplier of indexSet()

      SP::InteractionsGraph indexSet = _nsds->topology()->indexSet(0);
      InteractionsGraph::VIterator ui, uiend;
      for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter = indexSet->bundle(*ui);
        inter->lambda(0)->zero();
      }




      info = 0;
#ifdef TSPROJ_DEBUG_LEVEL1
      std::cout << "TimeSteppingCombinedProjection compute OSNSP POS." <<std::endl ;
#endif
      info = computeOneStepNSProblem(SICONOS_OSNSP_TS_POS);


      if (info)
      {
        std::cout << " TimeSteppingCombinedProjection::advanceToEvent() project on constraints. solver failed." <<std::endl ;
        return;
      }



      _nsds->updateInput(nextTime(),0);



      for (DynamicalSystemsGraph::VIterator aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2)
      {
        SP::DynamicalSystem ds = dsGraph->bundle(*aVi2);
        Type::Siconos dsType = Type::value(*ds);
        VectorOfVectors& workVectors = *dsGraph->properties(*aVi2).workVectors;
           
        if (dsType == Type::NewtonEulerDS)
        {
          SP::NewtonEulerDS neds = std11::static_pointer_cast<NewtonEulerDS>(ds);
          SP::SiconosVector q = neds->q();
          
          
          SP::SiconosVector qtmp = workVectors[OneStepIntegrator::qtmp];
          if (neds->p(0))
          {
            //*q = * qtmp +  *neds->p(0);
            *q +=  *neds->p(0);
          }
          neds->normalizeq();
          //neds->computeT();

#ifdef TSPROJ_DEBUG_LEVEL1
          neds->display();
#endif


        }
        else if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
        {
          SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
          SP::SiconosVector q = d->q();
          SP::SiconosVector qtmp = workVectors[OneStepIntegrator::qtmp];
          if (d->p(0))
          {
            //*q = * qtmp +  *d->p(0);
            *q +=  *d->p(0);
          }
#ifdef TSPROJ_DEBUG_LEVEL1
          std::cout << " q=" << std::endl;
          q->display();
          std::cout << " p(0)=" << std::endl;
          d->p(0)->display();
          std::cout << " p(1)=" << std::endl;
          d->p(1)->display();
#endif
        }
        else
          RuntimeException::selfThrow("TimeSteppingCombinedProjection::advanceToEvent() - Ds is not from NewtonEulerDS neither from LagrangianDS.");

      }

      updateWorldFromDS();

      computeCriteria(&runningProjection);


      //cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  Z:"<<endl;
      //(*_allNSProblems)[SICONOS_OSNSP_TS_POS]->display();
      //(std11::static_pointer_cast<LinearOSNS>((*_allNSProblems)[SICONOS_OSNSP_TS_POS]))->z()->display();

#ifdef TSPROJ_DEBUG_LEVEL1

      // SP::InteractionsGraph indexSet1 = _nsds->topology()->indexSet(1);
      // std ::cout << "lambda(1) in IndexSet1" << std::endl;
      // for (std11::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
      // {
      //   SP::Interaction inter = indexSet1->bundle(*ui);
      //   inter->lambda(1)->display();
      // }
      SP::InteractionsGraph indexSet2 = _nsds->topology()->indexSet(2);
      std ::cout << "lambda(0) in indexSet2" << std::endl;
      for (std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter = indexSet2->bundle(*ui);
        inter->lambda(0)->display();
      }



#endif

      //cout<<"during projection before normalizing of q:\n";
      //for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
      //{
      //  (*it)->relation()->computeh(getTkp1());
      //}
    } // end while ((runningProjection && _nbProjectionIteration < _projectionMaxIteration) && _doCombinedProj)

    DEBUG_PRINTF( "TimeSteppingCombinedProjection::Projection end : Number of iterations= %i\n", _nbProjectionIteration);

    _nbCumulatedProjectionIteration += _nbProjectionIteration ;
    if (_nbProjectionIteration == _projectionMaxIteration)
    {
      std::cout << "TimeSteppingCombinedProjection::advanceToEvent() Max number of projection iterations reached (" << _nbProjectionIteration << ")"  <<std::endl ;
      printf("              max criteria equality =  %e.\n", _maxViolationEquality);
      printf("              max criteria unilateral =  %e.\n", _maxViolationUnilateral);
    }

#endif // TSPROJ_WITHOUT_PROJECTION


    DEBUG_PRINT( "TimeSteppingCombinedProjection::newtonSolve end projection:\n");

    // We update forces to start the Newton Loop the next tiem step with a correct value in swap
    for (DynamicalSystemsGraph::VIterator aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2)
      {
        SP::DynamicalSystem ds = dsGraph->bundle(*aVi2);
        Type::Siconos dsType = Type::value(*ds);
        if (dsType == Type::NewtonEulerDS)
        {
          SP::NewtonEulerDS neds = std11::static_pointer_cast<NewtonEulerDS>(ds);
          double time = nextTime();
          neds->computeForces(time, neds->q(), neds->twist());
        }
        else if (dsType == Type::LagrangianDS)
        {
          SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
          double time = nextTime();
          d->computeForces(time, d->q(),d->velocity());
        }
        else if (dsType == Type::LagrangianLinearTIDS)
        {
        }
        else
          RuntimeException::selfThrow("TimeSteppingCombinedProjection::advanceToEvent() - Ds is not from NewtonEulerDS neither from LagrangianDS.");
      }


    if (_nsds->topology()->numberOfIndexSet() > _indexSetLevelForProjection)
    {
      updateIndexSet(1);
    }
    else
    {
      _isIndexSetsStable = true;
    }
#ifdef TSPROJ_DEBUG_LEVEL1

    if (topo->numberOfIndexSet() > _indexSetLevelForProjection)
    {
      SP::InteractionsGraph indexSet1 = topo->indexSet(1);
      SP::InteractionsGraph indexSet2 = topo->indexSet(2);
      std::cout << "indexSet1->size() " << indexSet1->size()   <<std::endl;
      std::cout << "indexSet2->size() " << indexSet2->size()   <<std::endl;
    }

    level = 0;
    for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      inter =indexSet0->bundle(*ui);
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 0);
      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 1);

      std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     << std::endl;
      std::cout << "inter->y(" << level << ")"   << std::endl;
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")"   << std::endl;
      inter->y_k(level)->display();
    }
    level = 1;
    for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      inter =indexSet0->bundle(*ui);
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     << std::endl;
      std::cout << "inter->y(" << level << ")"   << std::endl;
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")"   << std::endl;
      inter->y_k(level)->display();
    }

#endif

  }// end  while (!_isIndexSetsStable)

  DEBUG_PRINTF( "TimeSteppingCombinedProjection::indexset stable end : Number of iterations= %i \n ",_nbIndexSetsIteration);


  return;
}

void TimeSteppingCombinedProjection::computeCriteria(bool * runningProjection)
{
  DEBUG_PRINT("TimeSteppingCombinedProjection::computeCriteria(bool * runningProjection)\n");
  // SP::InteractionsGraph indexSet = _nsds->topology()->indexSet(_indexSetLevelForProjection);
  SP::InteractionsGraph indexSet = _nsds->topology()->indexSet(_indexSetLevelForProjection);

  InteractionsGraph::VIterator aVi, viend;

  double maxViolationEquality = -1e24;
  double maxViolationUnilateral = -1e24;

  *runningProjection = false;

  for (std11::tie(aVi, viend) = indexSet->vertices();
       aVi != viend; ++aVi)
  {
    SP::Interaction interac = indexSet->bundle(*aVi);

    interac->computeOutput(getTkp1(), indexSet->properties(*aVi), 0);
    interac->relation()->computeJach(getTkp1(), *interac, indexSet->properties(*aVi));

    if (Type::value(*(interac->nonSmoothLaw())) ==  Type::NewtonImpactFrictionNSL ||
        Type::value(*(interac->nonSmoothLaw())) == Type::NewtonImpactNSL)
    {



#ifdef TSPROJ_DEBUG_LEVEL1
      printf("  TimeSteppingCombinedProjection::computeCriteria  Unilateral interac->y(0)->getValue(0) %e.\n", interac->y(0)->getValue(0));
#endif
      if (!_doCombinedProjOnEquality)
      {
        if (maxViolationUnilateral > _constraintTolUnilateral)
        {
          double criteria = std::max(0.0, - interac->y(0)->getValue(0));
          if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;

          *runningProjection = true;
#ifdef TSPROJ_DEBUG_LEVEL1
          printf("TSProj newton criteria unilateral true %e.\n", criteria);
#endif
        }
      }
      else
      {
        double criteria = interac->y(0)->getValue(0);
        if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;

        if (std::abs(criteria) >=  _constraintTolUnilateral)
        {
          *runningProjection = true;
#ifdef TSPROJ_DEBUG_LEVEL1
          printf("TSProj newton criteria unilateral true %e.\n", criteria);
#endif
        }
      }

    }
    else
    {
      DEBUG_EXPR(interac->y(0)->display(););
      if (interac->y(0)->normInf()  > maxViolationEquality)
      {
        DEBUG_EXPR(interac->y(0)->display(););
        maxViolationEquality = interac->y(0)->normInf() ;
      }
      if (interac->y(0)->normInf() > _constraintTol)
      {
        *runningProjection = true;
#ifdef TSPROJ_DEBUG_LEVEL1
        printf("TSProj  newton criteria equality true %e.\n", interac->y(0)->normInf());
#endif
      }
    }

  }

  _maxViolationUnilateral = maxViolationUnilateral;
  _maxViolationEquality = maxViolationEquality;



  DEBUG_PRINTF("              max criteria equality =  %e.\n", _maxViolationEquality);
  DEBUG_PRINTF("              max criteria unilateral =  %e.\n", _maxViolationUnilateral);



#ifdef TSPROJ_DEBUG_LEVEL1
  printf("TSProj newton min/max criteria projection\n");
  std::cout << "                 runningProjection  "  << *runningProjection << std::endl;
  printf("              max criteria equality =  %e.\n", maxViolationEquality);
  printf("              max criteria unilateral =  %e.\n", maxViolationUnilateral);
  //  printf("              min criteria unilateral =  %e.\n",minViolationUnilateral);
#endif
}

void TimeSteppingCombinedProjection::updateIndexSet(unsigned int i)
{
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.
  // boost::default_color_type is used to organize update in InteractionsGraph:
  // - white_color : undiscovered vertex (Interaction)
  // - gray_color : discovered vertex (Interactions) but searching descendants
  // - black_color : discovered vertex (Interaction) together with the descendants

  assert(_nsds);
  assert(_nsds->topology());

  SP::Topology topo = _nsds->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0 &&
         "TimeStepping::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].

  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);
  DynamicalSystemsGraph& DSG0= *nonSmoothDynamicalSystem()->dynamicalSystems();
  assert(indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  topo->setHasChanged(false);

  DEBUG_PRINTF("update indexSets start : indexSet0 size : %i\n", (int)(indexSet0->size()));

  // Check indexSet1

  if (i == 1)
  {
    InteractionsGraph::VIterator ui1, ui1end, v1next;

    std11::tie(ui1, ui1end) = indexSet1->vertices();
    _isIndexSetsStable = true ;

    DEBUG_PRINTF("update IndexSets start : indexSet1 size : %i\n", (int)(indexSet1->size()));
    // indexSet1->display();
    //Remove interactions from the indexSet1
    for (v1next = ui1; ui1 != ui1end; ui1 = v1next)
    {
      ++v1next;
      SP::Interaction inter1 = indexSet1->bundle(*ui1);
      if (indexSet0->is_vertex(inter1))
      {
        InteractionsGraph::VDescriptor ur1_descr0 = indexSet0->descriptor(inter1);
        assert((indexSet0->color(ur1_descr0) == boost::white_color));
        indexSet0->color(ur1_descr0) = boost::gray_color;
      }
      // else
      // {
      //   // Interactions is not in indexSet0 anymore.
      //   // ui1 becomes invalid
      //   indexSet1->remove_vertex(inter1);
      //   topo->setHasChanged(true);
      //   _isIndexSetsStable=false;
      // }
    }

    // indexSet0\indexSet1 scan
    InteractionsGraph::VIterator ui0, ui0end;
    //Add interaction in indexSet1
    for (std11::tie(ui0, ui0end) = indexSet0->vertices(); ui0 != ui0end; ++ui0)
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
            Type::value(*(indexSet0->bundle(*ui0)->interaction()->nonSmoothLaw())) == Type::EqualityConditionNSL ;
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
            //SP::OneStepIntegrator Osi = indexSet0->properties(*ui0).osi;
            // We assume that the integrator of the ds1 drive the update of the index set
            SP::DynamicalSystem ds1 = indexSet1->properties(*ui0).source;
            OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;
            
            activate = osi.addInteractionInIndexSet(inter0, i);
          }
          if (activate)
          {
            assert(!indexSet1->is_vertex(inter0));

            // vertex and edges insertion in indexSet1
            indexSet1->copy_vertex(inter0, *indexSet0);
            topo->setHasChanged(true);
            _isIndexSetsStable = false ;
            assert(indexSet1->is_vertex(inter0));
          }
        }
      }
    }
    indexSet1->update_vertices_indices();
    indexSet1->update_edges_indices();
    assert(indexSet1->size() <= indexSet0->size());
    DEBUG_PRINTF("update indexSets end : indexSet0 size : %i\n", (int)(indexSet0->size()));
    DEBUG_PRINTF("update IndexSets end : indexSet1 size : %i\n", (int)(indexSet1->size()));
  } // i==1

  if (i == 2)
  {
    InteractionsGraph::VIterator ui1, ui1end, v1next;
    std11::tie(ui1, ui1end) = indexSet2->vertices();

    for (v1next = ui1; ui1 != ui1end; ui1 = v1next)
    {
      ++v1next;
      indexSet2->eraseProperties(*ui1);
      InteractionsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = indexSet2->out_edges(*ui1);
           oei != oeiend; ++oei)
      {
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet2->edges(indexSet2->source(*oei), indexSet2->target(*oei));
        if (ed2 != ed1)
        {
          indexSet2->eraseProperties(ed1);
          indexSet2->eraseProperties(ed2);
        }
        else
        {
          indexSet2->eraseProperties(ed1);
        }
      }
    }




    indexSet2->clear();
    DEBUG_PRINTF("update IndexSets start : indexSet2 size : %i\n", (int)(indexSet2->size()));

    // Scan indexSet1
    std11::tie(ui1, ui1end) = indexSet1->vertices();
    for (v1next = ui1; ui1 != ui1end; ui1 = v1next)
    {
      ++v1next;
      SP::Interaction inter1 = indexSet1->bundle(*ui1);
      bool activate = true;
      if (Type::value(*(inter1->nonSmoothLaw())) != Type::EqualityConditionNSL)
      {
        //SP::OneStepIntegrator Osi = indexSet1->properties(*ui1).osi;
        // We assume that the integrator of the ds1 drive the update of the index set
        SP::DynamicalSystem ds1 = indexSet1->properties(*ui1).source;
        OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;

        activate = osi.addInteractionInIndexSet(inter1, i);
      }
      if (activate)
      {
        assert(!indexSet2->is_vertex(inter1));

        // vertex and edges insertion in indexSet2
        indexSet2->copy_vertex(inter1, *indexSet1);
        topo->setHasChanged(true);
        assert(indexSet2->is_vertex(inter1));
      }
    }
    DEBUG_PRINTF("update IndexSets end : indexSet0 size : %i\n", (int)(indexSet0->size()));
    DEBUG_PRINTF("update IndexSets end : indexSet2 size : %i\n", (int)(indexSet2->size()));
    indexSet2->update_vertices_indices();
    indexSet2->update_edges_indices();

  }



}
