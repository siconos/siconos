/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "ZeroOrderHoldOSI.hpp"
#include "SiconosAlgebraProd.hpp"
#include "EventDriven.hpp"
#include "EventsManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderLinearR.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "SubPluggedObject.hpp"
#include "FirstOrderType2R.hpp"
#include "MatrixIntegrator.hpp"
#include "ExtraAdditionalTerms.hpp"
#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"

//#define DEBUG_WHERE_MESSAGES

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"


using namespace RELATION;
// --- constructor from a set of data ---
//ZeroOrderHoldOSI::ZeroOrderHoldOSI():
//  OneStepIntegrator(OSI::ZOHOSI)
//{
//}

// --- constructor from a minimum set of data ---
ZeroOrderHoldOSI::ZeroOrderHoldOSI():
  OneStepIntegrator(OSI::ZOHOSI), _useGammaForRelation(false)
{
  _steps = 1;
  _levelMinForOutput= 0;
  _levelMaxForOutput =0;
  _levelMinForInput =0;
  _levelMaxForInput =0;
}

void ZeroOrderHoldOSI::initializeWorkVectorsForDS(double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  DEBUG_BEGIN("void ZeroOrderHoldOSI::initializeWorkVectorsForDS( double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(ZeroOrderHoldOSI::WORK_LENGTH);

  DynamicalSystemsGraph& DSG0 = *_dynamicalSystemsGraph;
  InteractionsGraph& IG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

  Type::Siconos dsType = Type::value(*ds);

  if((dsType != Type::FirstOrderLinearDS) && (dsType != Type::FirstOrderLinearTIDS))
    THROW_EXCEPTION("ZeroOrderHoldOSI::initialize - the DynamicalSystem does not have the right type");
  unsigned int indxIter = 0;
  siconos::graphs::DynamicalSystemsGraph::AVIterator avi, aviend;
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  if(!DSG0.Ad.hasKey(dsgVD))
  {
    DSG0.Ad[dsgVD] = std::make_shared<siconos::simulation::MatrixIntegrator>(*ds, *_simulation->nonSmoothDynamicalSystem(), _simulation->eventsManager()->timeDiscretisation()));
    if(DSG0.Ad.at(dsgVD)->isConst())
      DSG0.Ad.at(dsgVD)->integrate();
  }
  else
    THROW_EXCEPTION("ZeroOrderHoldOSI::initialize - Ad MatrixIntegrator is already initialized for ds the DS");

  if((static_cast<const FirstOrderLinearDS&>(*ds)).b())
  {
    auto E = std::make_shared<siconos::algebra::SimpleMatrix>(ds->n(), ds->n(), 0));
    E->eye();
    DSG0.AdInt.insert(dsgVD, SP::MatrixIntegrator(new MatrixIntegrator(*ds,* _simulation->nonSmoothDynamicalSystem(),_simulation->eventsManager()->timeDiscretisation(), E)));
    if(DSG0.AdInt.at(dsgVD)->isConst())
      DSG0.AdInt.at(dsgVD)->integrate();
  }

  // init extra term, usually to add control terms
  if(_extraAdditionalTerms)
    _extraAdditionalTerms->init(DSG0, *_simulation->nonSmoothDynamicalSystem(), _simulation->eventsManager()->timeDiscretisation());

  // Now we search for an Interaction dedicated to control
  for(std::tie(avi, aviend) = DSG0.adjacent_vertices(dsgVD);
      avi != aviend; ++avi)
  {
    siconos::graphs::DynamicalSystemsGraph::EDescriptor ed1, ed2;
    std::tie(ed1, ed2) = DSG0.edges(dsgVD, *avi);

    if(IG0.properties(IG0.descriptor(DSG0.bundle(ed1))).forControl)
    {
      auto& inter = *DSG0.bundle(ed1);
      Relation& rel = *inter.relation();
      if(rel.getType() != RELATION::FirstOrder)
        THROW_EXCEPTION("ZeroOrderHoldOSI::initialize - the Integrator can only deal with FirstOrder Relation");
      FirstOrderR& relR = static_cast<FirstOrderR&>(rel);

      if(indxIter == 0)
      {
        indxIter++;
        if(!relR.getPluginJacLg()->isPlugged())
        {
          DSG0.Bd[dsgVD] = std::make_shared<siconos::simulation::MatrixIntegrator>(*ds,*_simulation->nonSmoothDynamicalSystem(),_simulation->eventsManager()->timeDiscretisation(), relR.B()));
          if(DSG0.Bd.at(dsgVD)->isConst())
            DSG0.Bd.at(dsgVD)->integrate();
        }
        else
        {
          DSG0.Bd[dsgVD] = std::make_shared<siconos::simulation::MatrixIntegrator>(*ds, *_simulation->nonSmoothDynamicalSystem(),_simulation->eventsManager()->timeDiscretisation(), relR.getPluging(), inter.dimension()));
        }
      }
      else
      {
        //        THROW_EXCEPTION("ZeroOrderHoldOSI::initialize - DS linked with more that one iteraction");
        DEBUG_PRINTF("number of iteraction attached to the process : %d\n", indxIter);
      }
    }
  }

  // Get work buffers from the graph
  ds_work_vectors[ZeroOrderHoldOSI::FREE] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  ds_work_vectors[ZeroOrderHoldOSI::DELTA_X_FOR_RELATION] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  DEBUG_END("void ZeroOrderHoldOSI::initializeWorkVectorsForDS( double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
}

void ZeroOrderHoldOSI::initializeWorkVectorsForInteraction(siconos::modeling::Interaction&inter,
    siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph & DSG)
{
  auto ds1= interProp.source;
  auto ds2= interProp.target;
  assert(ds1);
  assert(ds2);

  if(!interProp.workVectors)
  {
    interProp.workVectors = std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>
    interProp.workVectors->resize(ZeroOrderHoldOSI::WORK_INTERACTION_LENGTH);
  }

  if(!interProp.workBlockVectors)
  {
    interProp.workBlockVectors = std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>();
    interProp.workBlockVectors->resize(ZeroOrderHoldOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  Relation &relation =  *inter.relation();
  auto relationType = relation.getType();
  auto relationSubType = inter.relation()->getSubType();

  if(!inter_work[ZeroOrderHoldOSI::OSNSP_RHS])
    inter_work[ZeroOrderHoldOSI::OSNSP_RHS] = std::make_shared<siconos::algebra::SiconosVector>(inter.dimension()));

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  if(!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  {
    THROW_EXCEPTION("ZeroOrderHoldOSI::initializeWorkVectorsForInteraction. The implementation is not correct for two different OSI for one interaction");
  }

  unsigned int xfree = ZeroOrderHoldOSI::xfree;

  auto &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if(relationType == FirstOrder)
  {

    if(relationSubType == NonLinearR || relationSubType == Type2R)
    {
      inter_work[ZeroOrderHoldOSI::H_ALPHA] = std::make_shared<siconos::algebra::SiconosVector>(inter.dimension()));
    }
    inter_work_block[xfree] = std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[xfree]->insertPtr(workVds1[ZeroOrderHoldOSI::FREE]);
  }


  if(ds1 != ds2)
  {
    auto &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if(relationType == siconos::modeling::RelationType::Lagrangian)
    {
      inter_work_block[ZeroOrderHoldOSI::xfree]->insertPtr(workVds2[ZeroOrderHoldOSI::FREE]);
    }
  }

  if(!inter_work_block[ZeroOrderHoldOSI::DELTA_X])
  {
    inter_work_block[ZeroOrderHoldOSI::DELTA_X] = std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[ZeroOrderHoldOSI::DELTA_X]->insertPtr(workVds1[ZeroOrderHoldOSI::DELTA_X_FOR_RELATION]);
  }
  else
    inter_work_block[ZeroOrderHoldOSI::DELTA_X]->setVectorPtr(0,workVds1[ZeroOrderHoldOSI::DELTA_X_FOR_RELATION]);
}

double ZeroOrderHoldOSI::computeResidu()
{
  DEBUG_BEGIN("double ZeroOrderHoldOSI::computeResidu()\n");
  // This function is used to compute the residu for each "MoreauJeanOSI-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    dsType = Type::value(*ds); // Its type
    // 1 - First Order Linear Systems
    if(dsType == Type::FirstOrderLinearDS)
    {
      // No residu with ZOH ...
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if(dsType == Type::FirstOrderLinearTIDS)
    {
      // No residu with ZOH ...
    }
    else
      THROW_EXCEPTION("ZeroOrderHoldOSI::computeResidu - not yet implemented for Dynamical system type: " + std::to_string(dsType));
  }
  DEBUG_END("double ZeroOrderHoldOSI::computeResidu()\n");
  return maxResidu;
}

void ZeroOrderHoldOSI::computeFreeState()
{
  DEBUG_BEGIN("void ZeroOrderHoldOSI::computeFreeState()\n");
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  // Operators computed at told have index i, and (i+1) at t.
  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  DynamicalSystemsGraph& DSG0 = *_dynamicalSystemsGraph;
  Type::Siconos dsType ; // Type of the current DS.

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  DEBUG_EXPR(display(););
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds); // Its type
    DEBUG_EXPR(ds->display(););
    DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
    auto& ds_work_vectors = *DSG0.properties(dsgVD).workVectors;
//    updateMatrices(dsDescr);
    if(dsType == Type::FirstOrderLinearTIDS || dsType == Type::FirstOrderLinearDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(*ds);
      // Check whether we have to recompute things
      if(!DSG0.Ad.at(dsgVD)->isConst())
        DSG0.Ad.at(dsgVD)->integrate();
      if(d.b() && !DSG0.AdInt.at(dsgVD)->isConst())
        DSG0.AdInt.at(dsgVD)->integrate();

      siconos::algebra::SiconosVector& xfree = *ds_work_vectors[ZeroOrderHoldOSI::FREE];
      DEBUG_EXPR(xfree.display(););

      siconos::algebra::prod(DSG0.Ad.at(dsgVD)->mat(), *d.x(), xfree); // xfree = Ad*xold
      DEBUG_EXPR(xfree.display(););
      if(d.b())
      {
        assert(DSG0.AdInt.hasKey(dsgVD));
        siconos::algebra::prod(DSG0.AdInt.at(dsgVD)->mat(), *d.b(), xfree, false); // xfree += AdInt*b
        DEBUG_EXPR(xfree.display(););
      }

      // add extra term, possible control terms
      if(_extraAdditionalTerms)
      {
        DEBUG_PRINT("add extra additional terms\n");
        _extraAdditionalTerms->addSmoothTerms(DSG0, dsgVD, h, xfree);
      }
      DEBUG_EXPR(xfree.display(););
    }
    else
      THROW_EXCEPTION("ZeroOrderHoldOSI::computeFreeState - not yet implemented for Dynamical system type: " + std::to_string(dsType));

  }
  DEBUG_END("void ZeroOrderHoldOSI::computeFreeState()\n");
}

void ZeroOrderHoldOSI::prepareNewtonIteration(double time)
{

  // siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  // for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  // {
  //   //if (!checkOSI(dsi)) continue;
  //   //auto ds = _dynamicalSystemsGraph->bundle(*dsi);
  //   //    computeMatrices(time, *ds);
  // }

  if(!_explicitJacobiansOfRelation)
  {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);
  }
}


struct ZeroOrderHoldOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  std::shared_ptr<siconos::modeling::Interaction> _inter;
  siconos::graphs::InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(OneStepNSProblem *p, std::shared_ptr<siconos::modeling::Interaction> inter, siconos::graphs::InteractionProperties& interProp) :
    _osnsp(p), _inter(inter), _interProp(interProp) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    siconos::algebra::SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[ZeroOrderHoldOSI::OSNSP_RHS];
    siconos::algebra::subscal(e, _inter->y_k(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord, false);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    siconos::algebra::SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[ZeroOrderHoldOSI::OSNSP_RHS];
    osnsp_rhs(0) +=  e * _inter->y_k(_osnsp->inputOutputLevel())(0);

  }
  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
};


void ZeroOrderHoldOSI::computeFreeOutput(siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem * osnsp)
{
  DEBUG_BEGIN("void ZeroOrderHoldOSI::computeFreeOutput(...)\n");
  /** \warning: ensures that it can also work with two different osi for two different ds ?
  */
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(vertex_inter);
  auto& DSlink = inter->linkToDSVariables();
  auto& inter_work = *indexSet->properties(vertex_inter).workVectors;
  auto& inter_work_block = *indexSet->properties(vertex_inter).workBlockVectors;



  // Get relation and non smooth law types
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();

  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;



  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;


  std::shared_ptr<siconos::algebra::BlockVector> deltax;
  deltax = inter_work_block[ZeroOrderHoldOSI::DELTA_X];

  siconos::algebra::SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[ZeroOrderHoldOSI::OSNSP_RHS];

  std::shared_ptr<siconos::algebra::BlockVector> Xfree;
  if(relationType == FirstOrder)
  {
    Xfree = inter_work_block[ZeroOrderHoldOSI::xfree];
  }
  assert(Xfree);

  SP::Relation rel = inter->relation();
  assert(inter);
  assert(rel);

  //  if (!IG0.properties(IG0.descriptor(inter)).forControl) // the integration is not for control
  {
    if(relationType == FirstOrder && relationSubType == Type2R)
    {
      std::shared_ptr<siconos::algebra::SiconosVector> lambda = inter->lambda(0);
      auto C = rel->C();
      auto D = std::static_pointer_cast<FirstOrderType2R>(rel)->D();
      assert(lambda);

      if(D)
      {
        coord[3] = D->size(1);
        coord[5] = D->size(1);
        siconos::algebra::subprod(*D, *lambda, osnsp_rhs, coord, true);

        osnsp_rhs *= -1.0;
      }
      if(C)
      {
        coord[3] = C->size(1);
        coord[5] = C->size(1);
        siconos::algebra::subprod(*C, *deltax, osnsp_rhs, coord, false);

      }

      if(_useGammaForRelation)
      {
        THROW_EXCEPTION("ZeroOrderHoldOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
      }

      siconos::algebra::SiconosVector& hAlpha= *inter_work[ZeroOrderHoldOSI::H_ALPHA];
      osnsp_rhs += hAlpha;
    }

    else
    {
      auto C = rel->C();

      if(C)
      {

        assert(Xfree);
        assert(deltax);

        coord[3] = C->size(1);
        coord[5] = C->size(1);
        // creates a POINTER link between workX[ds] (xfree) and the
        // corresponding interactionBlock in each Interactionfor each ds of the
        // current Interaction.

        if(_useGammaForRelation)
        {
          siconos::algebra::subprod(*C, *deltax, osnsp_rhs, coord, true);
        }
        else
        {
          siconos::algebra::subprod(*C, *Xfree, osnsp_rhs, coord, true);
        }
      }

      if(relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
      {
        // In the first order linear case it may be required to add e + FZ to q.
        // q = HXfree + e + FZ
        std::shared_ptr<siconos::algebra::SiconosVector> e;
        std::shared_ptr<siconos::algebra::SiconosMatrix> F;
        if(relationSubType == LinearTIR)
        {
          e = std::static_pointer_cast<FirstOrderLinearTIR>(rel)->e();
          F = std::static_pointer_cast<FirstOrderLinearTIR>(rel)->F();
        }
        else
        {
          e = std::static_pointer_cast<FirstOrderLinearR>(rel)->e();
          F = std::static_pointer_cast<FirstOrderLinearR>(rel)->F();
        }

        if(e)
          osnsp_rhs += *e;

        if(F)
        {
          coord[3] = F->size(1);
          coord[5] = F->size(1);
          siconos::algebra::subprod(*F, *DSlink[FirstOrderR::z], osnsp_rhs, coord, false);
        }
      }

    }
  }

  DEBUG_END("void ZeroOrderHoldOSI::computeFreeOutput(...)\n");
}
void ZeroOrderHoldOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  // This function should not be used
  THROW_EXCEPTION("ZeroOrderHoldOSI::integrate - should not be used");
}

void ZeroOrderHoldOSI::updateState(const unsigned int level)
{
  DEBUG_BEGIN("ZeroOrderHoldOSI::updateState(const unsigned int level)\n");
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph& DSG0 = *_dynamicalSystemsGraph;
  siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    DEBUG_EXPR(ds->display(););
    Type::Siconos dsType = Type::value(*ds);

    DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
    auto& ds_work_vectors = *DSG0.properties(dsgVD).workVectors;
    // 1 - First Order Linear Systems
    if(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(*ds);
      siconos::algebra::SiconosVector& x = *d.x();
      // 1 - First Order Linear Time Invariant Systems
      // \Phi is already computed
      x = *ds_work_vectors[ZeroOrderHoldOSI::FREE]; // x = xfree = Phi*xold (+ Bd*u ) (+  Ld*e)
      if(level != simulation::internal::LEVELMAX)
      {
        std::shared_ptr<siconos::modeling::Interaction> interC;
        // we have to find the control interaction
        for(std::tie(oei, oeiend) = DSG0.out_edges(dsgVD); oei != oeiend; ++oei)
        {
          if(DSG0.properties(*oei).forControl)
          {
            interC = DSG0.bundle(*oei);
            break;
          }
        }
        if(interC)
        {
          DEBUG_PRINT("A control interaction is found\n");
          MatrixIntegrator& Bd = *DSG0.Bd[dsgVD];
          if(!Bd.isConst())
          {
            Bd.integrate();
          }
          siconos::algebra::prod(Bd.mat(), *interC->lambda(0), x, false); // x += Bd*\lambda
        }
      }
      DEBUG_EXPR(ds->display(););
    }
    else
      THROW_EXCEPTION("ZeroOrderHoldOSI::updateState - not yet implemented for Dynamical system type: " + std::to_string(dsType));
  }
  DEBUG_END("ZeroOrderHoldOSI::updateState(const unsigned int level)\n");
}


bool ZeroOrderHoldOSI::addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = .5;
  DEBUG_PRINTF("ZeroOrderHoldOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!std::isnan(y));
  if(y <= 0)
  {
    DEBUG_PRINT("ZeroOrderHoldOSI::addInteractionInIndexSet ACTIVATE.\n");
  }
  return (y <= 0);
}


bool ZeroOrderHoldOSI::removeInteractionFromIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  return !(addInteractionInIndexSet(inter,i));
}

const SiconosMatrix& ZeroOrderHoldOSI::Ad(std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  DynamicalSystemsGraph& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  return DSG0.Ad.at(dsgVD)->mat();
}

const SiconosMatrix& ZeroOrderHoldOSI::Bd(std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  DynamicalSystemsGraph& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  return DSG0.Bd.at(dsgVD)->mat();
}


void ZeroOrderHoldOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== ZOH OSI display ======" <<std::endl;
  std::cout << "--------------------------------" << std::endl;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--> Phi of dynamical system number " <<  ": " <<    std::endl;
    Ad(ds).display();
    std::cout << "--> Psi of dynamical system number " <<  ": " <<    std::endl;
    Bd(ds).display();
  }
  std::cout << "================================" <<std::endl;
}

void ZeroOrderHoldOSI::updateMatrices(std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
//  DynamicalSystemsGraph& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  if (!DSG0.Ad[dsgVD]->isConst())
//    computeAd(dsgVD);
//  if (DSG0.Bd.hasKey(dsgVD) && !DSG->Bd[dsgVD]->isConst())
//    computeBd(ds);
}
