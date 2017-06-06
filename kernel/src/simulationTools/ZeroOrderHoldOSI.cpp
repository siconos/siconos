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

#include "ZeroOrderHoldOSI.hpp"
#include "EventDriven.hpp"
#include "Model.hpp"
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
#include "CxxStd.hpp"
#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"

//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>


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
  _steps = 0;
  _levelMinForOutput= 0;
  _levelMaxForOutput =0;
  _levelMinForInput =0;
  _levelMaxForInput =0;
}

void ZeroOrderHoldOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);

  DynamicalSystemsGraph& DSG0 = *_dynamicalSystemsGraph;
  InteractionsGraph& IG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

  Type::Siconos dsType = Type::value(*ds);

  if((dsType != Type::FirstOrderLinearDS) && (dsType != Type::FirstOrderLinearTIDS))
    RuntimeException::selfThrow("ZeroOrderHoldOSI::initialize - the DynamicalSystem does not have the right type");
  unsigned int indxIter = 0;
  DynamicalSystemsGraph::AVIterator avi, aviend;
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  if(!DSG0.Ad.hasKey(dsgVD))
  {
    DSG0.Ad[dsgVD].reset(new MatrixIntegrator(*ds, m));
    if(DSG0.Ad.at(dsgVD)->isConst())
      DSG0.Ad.at(dsgVD)->integrate();
  }
  else
    RuntimeException::selfThrow("ZeroOrderHoldOSI::initialize - Ad MatrixIntegrator is already initialized for ds the DS");

  if((static_cast<const FirstOrderLinearDS&>(*ds)).b())
  {
    SP::SiconosMatrix E(new SimpleMatrix(ds->n(), ds->n(), 0));
    E->eye();
    DSG0.AdInt.insert(dsgVD, SP::MatrixIntegrator(new MatrixIntegrator(*ds, m, E)));
    if(DSG0.AdInt.at(dsgVD)->isConst())
      DSG0.AdInt.at(dsgVD)->integrate();
  }

  // init extra term, usually to add control terms
  if(_extraAdditionalTerms)
    _extraAdditionalTerms->init(DSG0, m);

  // Now we search for an Interaction dedicated to control
  for(std11::tie(avi, aviend) = DSG0.adjacent_vertices(dsgVD);
      avi != aviend; ++avi)
  {
    DynamicalSystemsGraph::EDescriptor ed1, ed2;
    std11::tie(ed1, ed2) = DSG0.edges(dsgVD, *avi);

    if(IG0.properties(IG0.descriptor(DSG0.bundle(ed1))).forControl)
    {
      Interaction& inter = *DSG0.bundle(ed1);
      Relation& rel = *inter.relation();
      if(rel.getType() != RELATION::FirstOrder)
        RuntimeException::selfThrow("ZeroOrderHoldOSI::initialize - the Integrator can only deal with FirstOrder Relation");
      FirstOrderR& relR = static_cast<FirstOrderR&>(rel);

      if(indxIter == 0)
      {
        indxIter++;
        if(!relR.getPluginJacLg()->isPlugged())
        {
          DSG0.Bd[dsgVD].reset(new MatrixIntegrator(*ds, m, relR.B()));
          if(DSG0.Bd.at(dsgVD)->isConst())
            DSG0.Bd.at(dsgVD)->integrate();
        }
        else
        {
          DSG0.Bd[dsgVD].reset(new MatrixIntegrator(*ds, m, relR.getPluging(), inter.getSizeOfY()));
        }
      }
      else
      {
        //        RuntimeException::selfThrow("ZeroOrderHoldOSI::initialize - DS linked with more that one iteraction");
        DEBUG_PRINTF("number of iteraction attached to the process : %d\n", indxIter);
      }
    }
  }

  // Get work buffers from the graph
  workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
  workVectors[OneStepIntegrator::free].reset(new SiconosVector(ds->dimension()));
  workVectors[OneStepIntegrator::delta_x_for_relation].reset(new SiconosVector(ds->dimension()));
}

void ZeroOrderHoldOSI::fillDSLinks(Interaction &inter,
				     InteractionProperties& interProp,
				     DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);

  VectorOfVectors& workV = *interProp.workVectors;
  workV[FirstOrderR::osnsp_rhs].reset(new SiconosVector(inter.getSizeOfY()));

  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  Relation &relation =  *inter.relation();
  RELATION::TYPES relationType = relation.getType();

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  /* allocate and set work vectors for the osi */
  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  {
    RuntimeException::selfThrow("ZeroOrderHoldOSI::fillDSLinks. The implementation is not correct for two different OSI for one interaction");
  }

  VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == FirstOrder)
    {
      DSlink[FirstOrderR::xfree].reset(new BlockVector());
      DSlink[FirstOrderR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
    }


  if (ds1 != ds2)
    {
      VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      if (relationType == Lagrangian)
      {
        DSlink[FirstOrderR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
      }
    }

  if (!DSlink[FirstOrderR::deltax])
  {
    DSlink[FirstOrderR::deltax].reset(new BlockVector());
    DSlink[FirstOrderR::deltax]->insertPtr(workVds1[OneStepIntegrator::delta_x_for_relation]);
  }
  else
    DSlink[FirstOrderR::deltax]->setVectorPtr(0,workVds1[OneStepIntegrator::delta_x_for_relation]);
}

double ZeroOrderHoldOSI::computeResidu()
{

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
  DynamicalSystemsGraph::VIterator dsi, dsend;

  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

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
      RuntimeException::selfThrow("ZeroOrderHoldOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
  }

  return maxResidu;
}

void ZeroOrderHoldOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  // Operators computed at told have index i, and (i+1) at t.
  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  DynamicalSystemsGraph& DSG0 = *_dynamicalSystemsGraph;
  Type::Siconos dsType ; // Type of the current DS.

  DynamicalSystemsGraph::VIterator dsi, dsend;

  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds); // Its type

    DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
    VectorOfVectors& workVectors = *DSG0.properties(dsgVD).workVectors;
//    updateMatrices(dsDescr);
    if(dsType == Type::FirstOrderLinearTIDS || dsType == Type::FirstOrderLinearDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(*ds);
      // Check whether we have to recompute things
      if(!DSG0.Ad.at(dsgVD)->isConst())
        DSG0.Ad.at(dsgVD)->integrate();
      if(d.b() && !DSG0.AdInt.at(dsgVD)->isConst())
        DSG0.AdInt.at(dsgVD)->integrate();

      SiconosVector& xfree = *workVectors[OneStepIntegrator::free];
      prod(DSG0.Ad.at(dsgVD)->mat(), *d.x(), xfree); // xfree = Ad*xold
      if(d.b())
      {
        assert(DSG0.AdInt.hasKey(dsgVD));
        prod(DSG0.AdInt.at(dsgVD)->mat(), *d.b(), xfree, false); // xfree += AdInt*b
      }

      // add extra term, possible control terms
      if(_extraAdditionalTerms)
        _extraAdditionalTerms->addSmoothTerms(DSG0, dsgVD, h, xfree);
    }
    else
      RuntimeException::selfThrow("ZeroOrderHoldOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }

}

void ZeroOrderHoldOSI::prepareNewtonIteration(double time)
{

  // DynamicalSystemsGraph::VIterator dsi, dsend;

  // for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  // {
  //   //if (!checkOSI(dsi)) continue;
  //   //SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
  //   //    computeMatrices(time, *ds);
  // }
}


struct ZeroOrderHoldOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;
  InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter, InteractionProperties& interProp) :
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
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[FirstOrderR::osnsp_rhs];
    subscal(e, *_inter->y_k(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord, false);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[FirstOrderR::osnsp_rhs];
    osnsp_rhs(0) +=  e * (*_inter->y_k(_osnsp->inputOutputLevel()))(0);

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


void ZeroOrderHoldOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem * osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
  */
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
  VectorOfVectors& workV = *indexSet->properties(vertex_inter).workVectors;
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();

  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;



  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;


  SP::BlockVector deltax;
  deltax = DSlink[FirstOrderR::deltax];

  SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[FirstOrderR::osnsp_rhs];

  SP::BlockVector Xfree;
  if(relationType == FirstOrder)
  {
    Xfree = DSlink[FirstOrderR::xfree];
  }
  assert(Xfree);

  SP::Relation rel = inter->relation();
  assert(inter);
  assert(rel);

  //  if (!IG0.properties(IG0.descriptor(inter)).forControl) // the integration is not for control
  {
    if(relationType == FirstOrder && relationSubType == Type2R)
    {
      SP::SiconosVector lambda = inter->lambda(0);
      SP::SiconosMatrix C = rel->C();
      SP::SiconosMatrix D = std11::static_pointer_cast<FirstOrderType2R>(rel)->D();
      assert(lambda);

      if(D)
      {
        coord[3] = D->size(1);
        coord[5] = D->size(1);
        subprod(*D, *lambda, osnsp_rhs, coord, true);

        osnsp_rhs *= -1.0;
      }
      if(C)
      {
        coord[3] = C->size(1);
        coord[5] = C->size(1);
        subprod(*C, *deltax, osnsp_rhs, coord, false);

      }

      if(_useGammaForRelation)
      {
        RuntimeException::selfThrow("ZeroOrderHoldOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
      }

      SiconosVector& hAlpha= *workV[FirstOrderR::h_alpha];
      osnsp_rhs += hAlpha;
    }

    else
    {
      SP::SiconosMatrix C = rel->C();

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
          subprod(*C, *deltax, osnsp_rhs, coord, true);
        }
        else
        {
          subprod(*C, *Xfree, osnsp_rhs, coord, true);
        }
      }

      if(relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
      {
        // In the first order linear case it may be required to add e + FZ to q.
        // q = HXfree + e + FZ
        SP::SiconosVector e;
        SP::SiconosMatrix F;
        if(relationSubType == LinearTIR)
        {
          e = std11::static_pointer_cast<FirstOrderLinearTIR>(rel)->e();
          F = std11::static_pointer_cast<FirstOrderLinearTIR>(rel)->F();
        }
        else
        {
          e = std11::static_pointer_cast<FirstOrderLinearR>(rel)->e();
          F = std11::static_pointer_cast<FirstOrderLinearR>(rel)->F();
        }

        if(e)
          osnsp_rhs += *e;

        if(F)
        {
          coord[3] = F->size(1);
          coord[5] = F->size(1);
          subprod(*F, *DSlink[FirstOrderR::z], osnsp_rhs, coord, false);
        }
      }

    }
  }


}
void ZeroOrderHoldOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  // This function should not be used
  RuntimeException::selfThrow("ZeroOrderHoldOSI::integrate - should not be used");
}

void ZeroOrderHoldOSI::updateState(const unsigned int level)
{
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph& DSG0 = *_dynamicalSystemsGraph;
  DynamicalSystemsGraph::OEIterator oei, oeiend;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    Type::Siconos dsType = Type::value(*ds);

    DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
    VectorOfVectors& workVectors = *DSG0.properties(dsgVD).workVectors;
    // 1 - First Order Linear Systems
    if(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(*ds);
      SiconosVector& x = *d.x();
      // 1 - First Order Linear Time Invariant Systems
      // \Phi is already computed
      x = *workVectors[OneStepIntegrator::free]; // x = xfree = Phi*xold (+ Bd*u ) (+  Ld*e)
      if(level != LEVELMAX)
      {
        SP::Interaction interC;
        // we have to find the control interaction
        for(std11::tie(oei, oeiend) = DSG0.out_edges(dsgVD); oei != oeiend; ++oei)
        {
          if(DSG0.properties(*oei).forControl)
          {
            interC = DSG0.bundle(*oei);
            break;
          }
        }
        if(interC)
        {
          MatrixIntegrator& Bd = *DSG0.Bd[dsgVD];
          if(!Bd.isConst())
          {
            Bd.integrate();
          }
          prod(Bd.mat(), *interC->lambda(0), x, false); // x += Bd*\lambda
        }
      }
    }
    else
      RuntimeException::selfThrow("ZeroOrderHoldOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}


bool ZeroOrderHoldOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = .5;
  DEBUG_PRINTF("ZeroOrderHoldOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));
  if(y <= 0)
  {
    DEBUG_PRINT("ZeroOrderHoldOSI::addInteractionInIndexSet ACTIVATE.\n");
  }
  return (y <= 0);
}


bool ZeroOrderHoldOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = .5;
  DEBUG_PRINTF("ZeroOrderHoldOSI::removeInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));
  if(y > 0)
  {
    DEBUG_PRINT("ZeroOrderHoldOSI::removeInteractionInIndexSet DEACTIVATE.\n");
  }
  return (y > 0);
}

const SiconosMatrix& ZeroOrderHoldOSI::Ad(SP::DynamicalSystem ds)
{
  DynamicalSystemsGraph& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  return DSG0.Ad.at(dsgVD)->mat();
}

const SiconosMatrix& ZeroOrderHoldOSI::Bd(SP::DynamicalSystem ds)
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
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
//    cout << "--> Phi of dynamical system number " << itN << ": " << endl;
//    if (Ad(ds)) Ad(ds)->display();
//    else cout << "-> NULL" << endl;
//    cout << "--> Psi of dynamical system number " << itN << ": " << endl;
//    if (Bd(ds)) Bd(ds)->display();
//    else cout << "-> NULL" << endl;
  }
  std::cout << "================================" <<std::endl;
}

void ZeroOrderHoldOSI::updateMatrices(SP::DynamicalSystem ds)
{
//  DynamicalSystemsGraph& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  if (!DSG0.Ad[dsgVD]->isConst())
//    computeAd(dsgVD);
//  if (DSG0.Bd.hasKey(dsgVD) && !DSG->Bd[dsgVD]->isConst())
//    computeBd(ds);
}
