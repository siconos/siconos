/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "EulerMoreauOSI.hpp"
#include "Simulation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderNonLinearR.hpp"
#include "FirstOrderType2R.hpp"
#include "FirstOrderType1R.hpp"
#include "NonSmoothLaw.hpp"
#include "CxxStd.hpp"
#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include <debug.h>

using namespace RELATION;

// --- constructor with theta parameter value  ---
EulerMoreauOSI::EulerMoreauOSI(double theta):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =0;
  _levelMinForInput =0;
  _levelMaxForInput =0;
  _steps=1;
  _theta = theta;
}

// --- constructor from a set of data ---
EulerMoreauOSI::EulerMoreauOSI(double theta, double gamma):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _useGammaForRelation(false)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =0;
  _levelMinForInput =0;
  _levelMaxForInput =0;
  _steps=1;
  _theta = theta;
  _gamma = gamma;
  _useGamma = true;
}

const SimpleMatrix EulerMoreauOSI::getW(SP::DynamicalSystem ds)
{
  assert(ds &&
         "EulerMoreauOSI::getW(ds): ds == NULL.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "EulerMoreauOSI::getW(ds): W[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W); // Copy !!
}

SP::SimpleMatrix EulerMoreauOSI::W(SP::DynamicalSystem ds)
{
  assert(ds && "EulerMoreauOSI::W(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
}


const SimpleMatrix EulerMoreauOSI::getWBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds &&
         "EulerMoreauOSI::getWBoundaryConditions(ds): ds == NULL.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions
         && "EulerMoreauOSI::getWBoundaryConditions(ds): WBoundaryConditions[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions); // Copy !!
}

SP::SiconosMatrix EulerMoreauOSI::WBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds && "EulerMoreauOSI::WBoundaryConditions(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions;
}

void EulerMoreauOSI::initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds)
{
  VectorOfVectors& ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(EulerMoreauOSI::WORK_LENGTH);

  // Check dynamical system type
  SP::FirstOrderNonLinearDS fods = std11::static_pointer_cast<FirstOrderNonLinearDS> (ds);
  assert (Type::value(*ds) == Type::FirstOrderNonLinearDS ||
    Type::value(*ds) == Type::FirstOrderLinearDS ||
    Type::value(*ds) == Type::FirstOrderLinearTIDS);
  // Compute W (iteration matrix)
  initializeIterationMatrixW(t, ds);

  // buffers allocation (into the graph)
  ds_work_vectors[EulerMoreauOSI::RESIDU].reset(new SiconosVector(ds->dimension()));
  ds_work_vectors[EulerMoreauOSI::RESIDU_FREE].reset(new SiconosVector(ds->dimension()));
  ds_work_vectors[EulerMoreauOSI::FREE].reset(new SiconosVector(ds->dimension()));
  ds_work_vectors[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION].reset(new SiconosVector(ds->dimension()));
  ds_work_vectors[EulerMoreauOSI::DELTA_X_FOR_RELATION].reset(new SiconosVector(ds->dimension()));
  ds_work_vectors[EulerMoreauOSI::LOCAL_BUFFER].reset(new SiconosVector(ds->dimension()));

  // Update dynamical system components (for memory swap).
  fods->computef(t, fods->x()); // Only fold is concerned, for FirstOrderNonLinearDS.
  // Update memory buffers
  ds->swapInMemory();

}

void EulerMoreauOSI::initializeWorkVectorsForInteraction(Interaction &inter,
                                 InteractionProperties& interProp,
                                 DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);

  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();

  if (!interProp.workVectors)
  {
    interProp.workVectors.reset(new VectorOfVectors);
    interProp.workVectors->resize(EulerMoreauOSI::WORK_INTERACTION_LENGTH);
  }
  if (!interProp.workMatrices)
  {
    interProp.workMatrices.reset(new VectorOfSMatrices);
    interProp.workMatrices->resize(EulerMoreauOSI::MAT_WORK_LENGTH);
  }
  if (!interProp.workBlockVectors)
  {
    interProp.workBlockVectors.reset(new VectorOfBlockVectors);
    interProp.workBlockVectors->resize(EulerMoreauOSI::BLOCK_WORK_LENGTH);
  }

  VectorOfVectors& inter_work = *interProp.workVectors;
  VectorOfSMatrices& inter_work_mat = *interProp.workMatrices;
  VectorOfBlockVectors& inter_work_block = *interProp.workBlockVectors;

  Relation &relation =  *inter.relation();

  unsigned int sizeY = inter.dimension();
  unsigned int sizeOfDS = inter.getSizeOfDS();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();

  RELATION::TYPES relationType = relation.getType();
  RELATION::SUBTYPES relationSubType = relation.getSubType();
  inter_work[EulerMoreauOSI::OSNSP_RHS].reset(new SiconosVector(inter.dimension()));

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  if(checkOSI(DSG.descriptor(ds1)))
  {
    DEBUG_PRINTF("ds1->number() %i is taken in to account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;


    if (relationType == FirstOrder)
    {
      inter_work[EulerMoreauOSI::VEC_Z].reset(new SiconosVector(sizeZ));
      inter_work[EulerMoreauOSI::VEC_X].reset(new SiconosVector(sizeOfDS));

      if (relationSubType == NonLinearR || relationSubType == Type2R )
      {
        inter_work[EulerMoreauOSI::H_ALPHA].reset(new SiconosVector(sizeY));
        inter_work[EulerMoreauOSI::VEC_RESIDU_Y].reset(new SiconosVector(sizeY));
        inter_work[EulerMoreauOSI::G_ALPHA].reset(new SiconosVector(sizeOfDS));
        inter_work[EulerMoreauOSI::VEC_RESIDU_R].reset(new SiconosVector(sizeOfDS));
        inter_work_mat[EulerMoreauOSI::MAT_KHAT].reset(new SimpleMatrix(sizeOfDS, sizeY));
        inter_work_mat[EulerMoreauOSI::MAT_KTILDE].reset(new SimpleMatrix(sizeOfDS, sizeY));
      }


      if (!inter_work_block[EulerMoreauOSI::XFREE])
      {
        inter_work_block[EulerMoreauOSI::XFREE].reset(new BlockVector());
        inter_work_block[EulerMoreauOSI::XFREE]->insertPtr(workVds1[EulerMoreauOSI::FREE]);
      }
      else
        inter_work_block[EulerMoreauOSI::XFREE]->setVectorPtr(0,workVds1[EulerMoreauOSI::FREE]);

      if (!inter_work_block[EulerMoreauOSI::X_PARTIAL_NS])
      {
        inter_work_block[EulerMoreauOSI::X_PARTIAL_NS].reset(new BlockVector());
        inter_work_block[EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(workVds1[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
      }
      else
        inter_work_block[EulerMoreauOSI::X_PARTIAL_NS]->setVectorPtr(0,workVds1[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
      if (!inter_work_block[EulerMoreauOSI::DELTA_X])
      {
        inter_work_block[EulerMoreauOSI::DELTA_X].reset(new BlockVector());
        inter_work_block[EulerMoreauOSI::DELTA_X]->insertPtr(workVds1[EulerMoreauOSI::DELTA_X_FOR_RELATION]);
      }
      else
        inter_work_block[EulerMoreauOSI::DELTA_X]->setVectorPtr(0,workVds1[EulerMoreauOSI::DELTA_X_FOR_RELATION]);
    }
  }
  DEBUG_PRINTF("ds1->number() %i\n",ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n",ds2->number());

  if (ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");

    if(checkOSI(DSG.descriptor(ds2)))
    {
      DEBUG_PRINTF("ds2->number() %i is taken in to account\n",ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      if (relationType == FirstOrder)
      {
        if (!inter_work_block[EulerMoreauOSI::XFREE])
        {
          inter_work_block[EulerMoreauOSI::XFREE].reset(new BlockVector());
          //dummy insertion to reserve first vector for ds1
          inter_work_block[EulerMoreauOSI::XFREE]->insertPtr(workVds2[EulerMoreauOSI::FREE]);
          inter_work_block[EulerMoreauOSI::XFREE]->insertPtr(workVds2[EulerMoreauOSI::FREE]);
        }
        else
          inter_work_block[EulerMoreauOSI::XFREE]->insertPtr(workVds2[EulerMoreauOSI::FREE]);

        if (!inter_work_block[EulerMoreauOSI::X_PARTIAL_NS])
        {
          inter_work_block[EulerMoreauOSI::X_PARTIAL_NS].reset(new BlockVector());
          //dummy insertion to reserve first vector for ds1
          inter_work_block[EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(workVds2[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
          inter_work_block[EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(workVds2[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
        }
        else
          inter_work_block[EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(workVds2[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);


        if (!inter_work_block[EulerMoreauOSI::DELTA_X])
        {
          inter_work_block[EulerMoreauOSI::DELTA_X].reset(new BlockVector());
          //dummy insertion to reserve first vector for ds1
          inter_work_block[EulerMoreauOSI::DELTA_X]->insertPtr(workVds2[EulerMoreauOSI::DELTA_X_FOR_RELATION]);
          inter_work_block[EulerMoreauOSI::DELTA_X]->insertPtr(workVds2[EulerMoreauOSI::DELTA_X_FOR_RELATION]);
        }
        else
          inter_work_block[EulerMoreauOSI::DELTA_X]->insertPtr(workVds2[EulerMoreauOSI::DELTA_X_FOR_RELATION]);
      }
    }
  }
}

void EulerMoreauOSI::initializeIterationMatrixW(double time, SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending on its type.

  if(!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);

  if(_dynamicalSystemsGraph->properties(dsv).W)
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - W(ds) is already in the map and has been initialized.");

  unsigned int sizeW = ds->dimension(); // n for first order systems, ndof for lagrangian.
  // Memory allocation for W
  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  // 1 - All 'First order' systems
  if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = std11::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Memory allocation for W property of the graph
    if(d->M()) // W = M
    {
      d->computeM(time);
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->M()));
    }
    else // W = I
    {
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(sizeW, sizeW));
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }

    SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(dsv).W;
    // Add -h*_theta*jacobian_XF to W
    if(d->jacobianfx())
    {
      d->computeJacobianfx(time, ds->x());
      scal(-h * _theta, *d->jacobianfx(), *W, false);
    }
  }
  else RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
}


void EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if(!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  Type::Siconos dsType = Type::value(*ds);


  RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void EulerMoreauOSI::computeWBoundaryConditions(SP::DynamicalSystem ds)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  assert(ds &&
         "EulerMoreauOSI::computeWBoundaryConditions(t,ds) - ds == NULL");

  Type::Siconos dsType = Type::value(*ds);
  //unsigned int dsN = ds->number();
  RuntimeException::selfThrow("EulerMoreauOSI::computeWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void EulerMoreauOSI::computeW(double time, DynamicalSystem& ds,
                              DynamicalSystemsGraph::VDescriptor& dsv,
                              SiconosMatrix& W)
{
  DEBUG_BEGIN("EulerMoreauOSI::computeW(...)\n");
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, W is supposed to exist and not to be null
  // Memory allocation has been done during initializeIterationMatrixW.

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(ds);

  // 1 - First order non linear systems
  if(dsType == Type::FirstOrderNonLinearDS)
  {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    FirstOrderNonLinearDS& d = static_cast<FirstOrderNonLinearDS&>(ds);

    // Copy M or I if M is Null into W
    if(d.M())
    {
      d.computeM(time);
      W = *d.M();
    }
    else
      W.eye();

    if(d.jacobianfx())
    {
      d.computeJacobianfx(time, d.x());
      // Add -h*_theta*jacobianfx to W
      scal(-h * _theta, *d.jacobianfx(), W, false);
    }

    DEBUG_EXPR(W.display(););

  }
  // 2 - First order linear systems
  else if(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    FirstOrderLinearDS& fods = static_cast<FirstOrderLinearDS&>(ds);
    if(dsType == Type::FirstOrderLinearDS )
    {
      fods.computeA(time);
      fods.computeM(time);
    }

    if(fods.M())
      W = *fods.M();
    else
      W.eye();

    if(fods.A())
      scal(-h * _theta, *fods.A(), W, false);
  }
  else RuntimeException::selfThrow("EulerMoreauOSI::computeW - not yet implemented for Dynamical system type :" + dsType);

  //  if (_useGamma)
  {

//    InteractionsGraph& indexSet = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(0);

    DynamicalSystemsGraph::OEIterator oei, oeiend;
//    InteractionsGraph::VDescriptor ivd;
    SP::SiconosMatrix K;
    SP::Interaction inter;
    for(std11::tie(oei, oeiend) = _dynamicalSystemsGraph->out_edges(dsv); oei != oeiend; ++oei)
    {
      inter = _dynamicalSystemsGraph->bundle(*oei);
      VectorOfSMatrices& relationMat = inter->relationMatrices();
//      ivd = indexSet.descriptor(inter);
      FirstOrderR& rel = static_cast<FirstOrderR&>(*inter->relation());
      K = rel.K();
      if(!K) K = relationMat[FirstOrderR::mat_K];
      if(K)
      {
        scal(-h * _gamma, *K, W, false);
      }
    }
  }
  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
  DEBUG_EXPR(W.display());
  DEBUG_END("EulerMoreauOSI::computeW(...)\n");
}

void EulerMoreauOSI::computeKhat(Interaction& inter, SiconosMatrix& m,
                                 VectorOfSMatrices& workM, double h) const
{
  RELATION::TYPES relationType = inter.relation()->getType();

  if ((relationType == FirstOrder) && (workM[EulerMoreauOSI::MAT_KHAT]))
  {
    SP::SiconosMatrix K = std11::static_pointer_cast<FirstOrderR>(inter.relation())->K();
    if (!K) K = inter.relationMatrices()[FirstOrderR::mat_K];
    prod(*K, m, *workM[EulerMoreauOSI::MAT_KHAT], true);
    *workM[EulerMoreauOSI::MAT_KHAT] *= h;
  }
}



double EulerMoreauOSI::computeResidu()
{
  DEBUG_BEGIN("EulerMoreauOSI::computeResidu()\n");
  // This function is used to compute the residu for each "EulerMoreauOSI-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double time = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = time - told; // time step length

  DEBUG_PRINTF("nextTime %f\n", time);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);


  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds); // Its type

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
    // Maurice will do that with subgraph :)
    SiconosVector& residuFree = *ds_work_vectors[EulerMoreauOSI::RESIDU_FREE];
    SiconosVector& residu = *ds_work_vectors[EulerMoreauOSI::RESIDU];

    // 1 - First Order Non Linear Systems AND First Order Linear DS
    if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS)
    {
      // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
      //  $\mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
      //  $\mathcal R_{free}(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

      // Note: indices k/k+1 corresponds to value at the beginning/end of the time step.
      // Newton iterate are x and r

      FirstOrderNonLinearDS& fonlds = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      FirstOrderLinearDS& folds = *std11::static_pointer_cast<FirstOrderLinearDS>(ds);

      // 1 - Compute the free residu (purely on the "smooth" dynamics)

      residuFree = *(fonlds.x()); // last saved value for x: could be x_k or x_{k+1}^alpha
      const SiconosVector& xold = fonlds.xMemory().getSiconosVector(0);
      residuFree -= xold; // state x_k (at previous time step)

      SP::SiconosMatrix M = fonlds.M();
      if(M)
      {
        fonlds.computeM(time);
        prod(*M, residuFree, residuFree, true);
      }
      // at this step, we have residuFree = M(x - x_k)
      DEBUG_PRINT("EulerMoreauOSI::computeResidu residuFree = M(x - x_k)\n");
      DEBUG_EXPR(residuFree.display());
      double coef = -h * (1 - _theta);
      if(dsType == Type::FirstOrderLinearDS)
      {
         // computes f(t_k,x_k)
         // No fold in FirstOrderLinearDS.
         // residu is used as a tmp buffer to compute Ax + b
         residu.zero();
         if(folds.A())
         {
            folds.computeA(told);
            prod(*folds.A(), xold, residu);
         }

         if(folds.b())
         {
          folds.computeb(told);
          residu += *folds.b();
        }
        DEBUG_EXPR(residuFree.display());
        // residuFree += -h * (1 - _theta) * f(t_k,x_k)
        scal(coef, residu, residuFree, false);
        residu.zero();
        if(folds.A())
        {
          folds.computeA(time);
          prod(*folds.A(), *folds.x(), residu);
        }
        if(folds.b())
        {
          folds.computeb(time);
          residu += *folds.b();
        }
        // residuFree += -h * _theta * f(t_{x+1}, x_{k+1}^alpha)
        coef = -h * _theta;
        scal(coef, residu, residuFree, false);
        DEBUG_PRINT("- 3 -\n");
        DEBUG_EXPR(residuFree.display());
        DEBUG_EXPR(xold.display());
        DEBUG_EXPR(folds.x()->display());
      }
      else if(dsType == Type::FirstOrderNonLinearDS) // FirstOrderNonLinearDS
      {
        DEBUG_PRINT("dsType == Type::FirstOrderNonLinearDS\n");
        DEBUG_EXPR(fonlds.f()->display(););
        if(fonlds.f())
        {
          coef = -h * (1 - _theta);
          // for these systems, fold is available
          // residuFree += -h * (1 - _theta) * f(t_k,x_k)
          scal(coef, *fonlds.fold(), residuFree, false);

          // computes f(t_{x+1}, x_{k+1}^alpha)
          fonlds.computef(time, fonlds.x());
          coef = -h * _theta;
          // residuFree += -h * _theta * f(t_{x+1}, x_{k+1}^alpha)
          scal(coef, *(fonlds.f()), residuFree, false);
        }
      }

      // now we compute the residu = residuFree - h*gamma*r - h*(1-gamma)r_k
      residu = residuFree;

      if(!_useGamma)  // no gamma
      {
        DEBUG_EXPR(fonlds.r()->display(););
        DEBUG_EXPR(residu.display());
        scal(-h, *fonlds.r(), residu, false); // residu = residu - h*r
      } else {
        scal(-h*_gamma, *fonlds.r(), residu, false);
        scal(-h*(1-_gamma), fonlds.rMemory().getSiconosVector(0), residu, false);
      }

      normResidu = residu.norm2();
      DEBUG_EXPR(residu.display());

    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if(dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearTIDS& foltids = *std11::static_pointer_cast<FirstOrderLinearTIDS>(ds);
      //Don't use W because it is LU factorized
      //Residu : R_{free} = M(x^{\alpha}_{k+1} - x_{k}) -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b)

      // 1. R_{free} = -h * b
      if(foltids.b())
        scal(-h, *(foltids.b()), residuFree, true);
      else
        residuFree.zero();

      // 2. residuFree += -h * A (\theta x_{k+1}^{\alpha} + (1-\theta) x_k)
      // residu is used as a temp buffer
      if(foltids.A()) {
        SP::SiconosMatrix A = foltids.A();
        prod(*A, foltids.xMemory().getSiconosVector(0), residu, true);
        double coef = -h * (1 - _theta);
        scal(coef, residu, residuFree, false);

        prod(*A, *(foltids.x()), residu, true);
        coef = -h * _theta;
        scal(coef, residu, residuFree, false);
      }

      // 3. residuFree += M(x_{k+1}^{\alpha} - x_k)
      residu = *(foltids.x()) - foltids.xMemory().getSiconosVector(0);
      SP::SiconosMatrix M = foltids.M();
      if(M) {
         prod(*M, residu, residuFree, false);
      } else {
        residuFree += residu;
      }
    }
    else
      RuntimeException::selfThrow("EulerMoreauOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);


    DEBUG_PRINT("EulerMoreauOSI::computeResidu final residuFree\n");
    DEBUG_EXPR(residuFree.display());

    if(normResidu > maxResidu) maxResidu = normResidu;

  }
  DEBUG_END("EulerMoreauOSI::computeResidu()\n");
  return maxResidu;
}

void EulerMoreauOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.
  DEBUG_BEGIN("EulerMoreauOSI::computeFreeState()\n");

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d.rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.



  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
    // Maurice will do that with subgraph :)

    VectorOfVectors& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds); // Its type
    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W; // Its W EulerMoreauOSI matrix of iteration.

    // 1 - First Order Non Linear Systems
    if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      // xfree =  x - W^{-1} (ResiduFree - h(1-gamma)*rold)
      // with ResiduFree = = M(x - x_k) - h*theta*f(t_{k+1}, x) - h*(1-theta)*f(t_k, x_k)

      // to be updated at current time: W, f
      // fold is f at t_k
      // not time dependant: M
      FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //    SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

      SiconosVector& x = *d.x(); // x = x_k or x = x_{k+1}^{\alpha}
      // xfree gets ResiduFree at first
      SiconosVector& xfree = *ds_work_vectors[EulerMoreauOSI::FREE];
      xfree = *ds_work_vectors[EulerMoreauOSI::RESIDU_FREE];

      DEBUG_PRINT("EulerMoreauOSI::computeFreeState xfree <- residuFree\n");
      DEBUG_EXPR(xfree.display());

      if (_useGamma)
      {
        const SiconosVector& rold = d.rMemory().getSiconosVector(0);
        double coeff = -h * (1 - _gamma);
        scal(coeff, rold, xfree, false); //  xfree += -h(1-gamma)*rold
      }


      // At this point xfree = (ResiduFree - h(1-gamma)*rold)
      // -> Solve WX = xfree and set xfree = X
      W.PLUForwardBackwardInPlace(xfree);

      // at this point, xfree = W^{-1} (ResiduFree - h(1-gamma)*rold)
      // -> compute real xfree = x - W^{-1} (ResiduFree - h(1-gamma)*rold)
      xfree *= -1.0;
      xfree += x;

      DEBUG_EXPR(xfree.display());

      // now the crazy intermediate variables
      // xPartialNS was updated before this fonction call
      // It constains either 0 (first Newton iterate)
      // or g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
      SiconosVector& xPartialNS = *ds_work_vectors[EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION];
      DEBUG_PRINT("EulerMoreauOSI::computeFreeState xPartialNS from Interaction\n");
      DEBUG_EXPR(xPartialNS.display());

      // -> Solve WX = g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
      // and set xPartialNS = X
      W.PLUForwardBackwardInPlace(xPartialNS);
      scal(h, xPartialNS, xPartialNS);

      // compute real xPartialNS = xfree + ...
      xPartialNS += xfree;
      DEBUG_PRINT("EulerMoreauOSI::computeFreeState xPartialNS real value\n");
      DEBUG_EXPR(xPartialNS.display());

      // deltaxForRelation = (\widetilde{K}_{k+1}^{\alpha})^{-1} xPartialNS - x
      SiconosVector& deltaxForRelation = *ds_work_vectors[EulerMoreauOSI::DELTA_X_FOR_RELATION];
      deltaxForRelation = xPartialNS;

      deltaxForRelation -= x;

      DEBUG_EXPR(deltaxForRelation.display());

      // have a look at the end of the DevNotes for this part
      if(_useGammaForRelation)
      {
        if(!(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS))
          RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - _useGammaForRelation == true is only implemented for FirstOrderLinearDS or FirstOrderLinearTIDS");

        deltaxForRelation = xfree;

        scal(_gamma, deltaxForRelation, deltaxForRelation);
        const SiconosVector& xold = d.xMemory().getSiconosVector(0);

        scal(1.0 - _gamma, xold, deltaxForRelation, false);
      }

      // some output
      DEBUG_EXPR(xfree.display(););
      DEBUG_EXPR(xPartialNS.display(););
      DEBUG_EXPR(deltaxForRelation.display(););

    }
    else
      RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }
  DEBUG_END("EulerMoreauOSI::computeFreeState()\n");
}


void EulerMoreauOSI::prepareNewtonIteration(double time)
{
  // XXX TMP hack -- xhub
  // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
  // Maurice will do that with subgraph :)
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    DynamicalSystemsGraph::VDescriptor dsv = _dynamicalSystemsGraph->descriptor(ds);
    SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(*dsi).W;
    computeW(time, *ds, dsv, *W);
  }

  if(!_explicitJacobiansOfRelation)
  {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);

    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

    for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      Interaction & inter = *indexSet0->bundle(*ui);
      InteractionProperties& interProp = indexSet0->properties(*ui);

      VectorOfVectors& inter_work = *interProp.workVectors;
      VectorOfSMatrices& relationMat = inter.relationMatrices();
      VectorOfBlockVectors& inter_work_block = *interProp.workBlockVectors;


      RELATION::TYPES relationType = inter.relation()->getType();
      RELATION::SUBTYPES relationSubType = inter.relation()->getSubType();
      if(relationType == FirstOrder)
      {
        FirstOrderR& relation = static_cast<FirstOrderR&> (*inter.relation());
        BlockVector& xPartialNS = *inter_work_block[EulerMoreauOSI::X_PARTIAL_NS];

        if (relationSubType == NonLinearR || relationSubType == Type2R)
        {
          if (relation.B())
            prod(*relation.B(), *inter.lambda(0), *inter_work[EulerMoreauOSI::VEC_X], true);
          else
            prod(*relationMat[FirstOrderR::mat_B], *inter.lambda(0), *inter_work[EulerMoreauOSI::VEC_X], true);

          xPartialNS = *inter_work[EulerMoreauOSI::G_ALPHA];
          xPartialNS -= *inter_work[EulerMoreauOSI::VEC_X];
        }
      }
    }
  }





}

/// @cond

struct EulerMoreauOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter) :
    _osnsp(p), _inter(inter) {};

  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
};

/// @endcond

void EulerMoreauOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */
  DEBUG_BEGIN("EulerMoreauOSI::computeFreeOutput(...)\n");
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = inter->linkToDSVariables();
  VectorOfSMatrices& relationMat = inter->relationMatrices();
  VectorOfVectors & relationVec= inter->relationVectors();

  VectorOfVectors& inter_work = *indexSet->properties(vertex_inter).workVectors;
  VectorOfBlockVectors& inter_work_block = *(indexSet->properties(vertex_inter)).workBlockVectors;
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
  SP::SiconosMatrix  C;
  SP::SiconosMatrix  D;
  SP::SiconosMatrix  F;
  SP::BlockVector deltax;
  SP::BlockVector Xfree;

  SP::SiconosVector H_alpha;

  deltax = inter_work_block[EulerMoreauOSI::DELTA_X];
  DEBUG_EXPR(deltax->display(););
  SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[EulerMoreauOSI::OSNSP_RHS];

  Xfree = inter_work_block[EulerMoreauOSI::XFREE];
  DEBUG_EXPR(Xfree->display(););
  assert(Xfree);


  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if(relationType == FirstOrder && (relationSubType == Type2R || relationSubType == NonLinearR))
  {
    DEBUG_PRINT("relationType == FirstOrder && (relationSubType == Type2R || relationSubType == NonLinearR)\n")
    SiconosVector& lambda = *inter->lambda(0);
    FirstOrderR& rel = *std11::static_pointer_cast<FirstOrderR>(mainInteraction->relation());
    C = rel.C();
    if(!C) C = relationMat[FirstOrderR::mat_C];
    D = rel.D();
    if(!D) D = relationMat[FirstOrderR::mat_D];

    if(D)
    {
      coord[3] = D->size(1);
      coord[5] = D->size(1);
      subprod(*D, lambda, osnsp_rhs, coord, true);

      osnsp_rhs *= -1.0;
    }
    else
    {
      subscal(0, osnsp_rhs, osnsp_rhs, coord, true);
    }

    if(C)
    {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *deltax, osnsp_rhs, coord, false);

    }

    if(_useGammaForRelation)
    {
      RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Type2R and H_alpha->getValue() should return the mid-point value");
    }
    SiconosVector& hAlpha= *inter_work[EulerMoreauOSI::H_ALPHA];
    DEBUG_EXPR(hAlpha.display());
    osnsp_rhs += hAlpha;
    DEBUG_EXPR(osnsp_rhs.display(););
  }
  else if(relationType == FirstOrder && relationSubType == Type1R)
  {
    DEBUG_PRINT("relationType == FirstOrder && relationSubType == Type1R\n");
    FirstOrderType1R& rel = *std11::static_pointer_cast<FirstOrderType1R>(mainInteraction->relation());
    C = rel.C();
    if(!C) C = relationMat[FirstOrderR::mat_C];
    F = rel.F();
    if(!F) F = relationMat[FirstOrderR::mat_F];
    assert(Xfree);
    assert(deltax);

    if(F)
    {
      coord[3] = F->size(1);
      coord[5] = F->size(1);
      subprod(*F, *DSlink[FirstOrderR::z], osnsp_rhs, coord, true);

    }
    if(C)
    {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, osnsp_rhs, coord, false);

    }

    if(_useGammaForRelation)
    {
      RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
    }
    if(inter_work[EulerMoreauOSI::H_ALPHA])
    {
      osnsp_rhs += *inter_work[EulerMoreauOSI::H_ALPHA];
    }
  }
  else // First Order Linear Relation
  {
    DEBUG_PRINT("relationType == FirstOrder\n");
    C = mainInteraction->relation()->C();
    if(!C) C = relationMat[FirstOrderR::mat_C];

    if(C)
    {

      assert(Xfree);
      assert(deltax);

      coord[3] = C->size(1);
      coord[5] = C->size(1);

      if(_useGammaForRelation)
      {
        subprod(*C, *deltax, osnsp_rhs, coord, true);
      }
      else
      {
        subprod(*C, *Xfree, osnsp_rhs, coord, true);
      }
    }
    DEBUG_EXPR(osnsp_rhs.display(););
    if(relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
    {
      // In the first order linear case it may be required to add e + FZ to y.
      // y = CXfree + e + FZ
      SP::SiconosVector e;
      if(relationSubType == LinearTIR)
      {
        e = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->e();
        F = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->F();
      }
      else
      {
        e = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->e();
        if(!e) e = relationVec[FirstOrderR::e];
        F = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->F();
        if(!F) F = relationMat[FirstOrderR::mat_F];
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
    DEBUG_EXPR(osnsp_rhs.display(););

  }
  DEBUG_END("EulerMoreauOSI::computeFreeOutput(...)\n");
}

void EulerMoreauOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  // Last parameter is not used (required for LsodarOSI but not for EulerMoreauOSI).

  //double h = tend - tinit;
  tout = tend;


  SP::SiconosMatrix W;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    RuntimeException::selfThrow("EulerMoreauOSI::integrate - not yet implemented for Dynamical system type :" + dsType);
  }
}

void EulerMoreauOSI::updateState(const unsigned int)
{

  DEBUG_PRINT("EulerMoreauOSI::updateState\n");

  double h = _simulation->timeStep();

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph::VIterator dsi, dsend;

  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    // Get the DS type
    Type::Siconos dsType = Type::value(*ds);
    VectorOfVectors& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    SimpleMatrix&  W = *_dynamicalSystemsGraph->properties(*dsi).W;

    if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      SiconosVector& x = *ds->x();
      DEBUG_PRINT("EulerMoreauOSI::updateState Old value of x\n");
      DEBUG_EXPR(x.display());
      DEBUG_PRINT("EulerMoreauOSI::updateState residu value\n");
      DEBUG_EXPR(d.r()->display());

      // TODO ???
      bool baux = (useRCC && dsType == Type::FirstOrderNonLinearDS && _simulation->relativeConvergenceCriterionHeld());

      //    SP::SiconosVector xFree = d->xFree();

      // Save value of q in local_buffer for relative convergence computation
      if(baux)
        *ds_work_vectors[EulerMoreauOSI::LOCAL_BUFFER] = x;

      if(_useGamma)
      {
        // XXX UseGamma broken ? -- xhub
        scal(_gamma * h, *d.r(), x); // x = gamma*h*r
      }
      else
      {
        scal(h, *d.r(), x); // x = h*r
      }

      W.PLUForwardBackwardInPlace(x); // x = h* W^{-1} *r

      x += *ds_work_vectors[EulerMoreauOSI::FREE]; // x+=xfree

      if(baux)
      {
        double ds_norm_ref = 1. + ds->x0()->norm2(); // Should we save this in the graph?
        *ds_work_vectors[EulerMoreauOSI::LOCAL_BUFFER] -= x;
        double aux = (ds_work_vectors[EulerMoreauOSI::LOCAL_BUFFER]->norm2()) / (ds_norm_ref);
        if(aux > RelativeTol)
          _simulation->setRelativeConvergenceCriterionHeld(false);
      }
      DEBUG_PRINT("EulerMoreauOSI::updateState New value of x\n");
      DEBUG_EXPR(x.display());
    }
    else RuntimeException::selfThrow("EulerMoreauOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}

void EulerMoreauOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== EulerMoreauOSI OSI display ======" <<std::endl;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--------------------------------" <<std::endl;
    std::cout << "--> W of dynamical system number " << ds->number() << ": " <<std::endl;
    if(_dynamicalSystemsGraph->properties(*dsi).W) _dynamicalSystemsGraph->properties(*dsi).W->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "--> and corresponding theta is: " << _theta <<std::endl;
  }
  std::cout << "================================" <<std::endl;
}
void EulerMoreauOSI::updateOutput(double time)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  for (unsigned int level = _levelMinForOutput;
       level < _levelMaxForOutput + 1;
       level++)
    updateOutput(time,level);
}

void EulerMoreauOSI::updateInput(double time)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  for (unsigned int level = _levelMinForInput;
       level < _levelMaxForInput + 1;
       level++)
    updateInput(time,level);
}

void EulerMoreauOSI::updateOutput(double time, unsigned int level)
{
  DEBUG_BEGIN("EulerMoreauOSI::updateOutput(double time, unsigned int level)\n");
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  //_simulation->nonSmoothDynamicalSystem()->updateOutput(time,level);
  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet0->bundle(*ui);
    assert(inter.lowerLevelForOutput() <= level);
    assert(inter.upperLevelForOutput() >= level);


    VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
    VectorOfSMatrices& relationMat = inter.relationMatrices();

    InteractionProperties& interProp = indexSet0->properties(*ui);
    VectorOfVectors& inter_work = *interProp.workVectors;
    VectorOfBlockVectors& inter_work_block = *interProp.workBlockVectors;
    RELATION::SUBTYPES relationSubType = inter.relation()->getSubType();
    if (relationSubType == Type2R)
    {
      FirstOrderType2R & r = static_cast<FirstOrderType2R&>(*inter.relation());
            // compute the new y obtained by linearisation (see DevNotes)
      // y_{alpha+1}_{k+1} = h(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      // or equivalently
      // y_{alpha+1}_{k+1} = y_{alpha}_{k+1} - ResiduY_{k+1}^{alpha}
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      SiconosVector& y = *inter.y(level);
      DEBUG_EXPR(y.display());



      if (r.D())
        prod(*r.D(), *(inter.lambdaOld(level)), y, true);
      else
        prod(*relationMat[FirstOrderR::mat_D], *(inter.lambdaOld(level)), y, true);

      y *= -1.0;
      //SiconosVector yOld = *inter.yOld(0); // Retrieve  y_{alpha}_{k+1}
      DEBUG_PRINT("FirstOrderType2R::computeOutput : yOld(level) \n");
      DEBUG_EXPR(inter.yOld(level)->display());

      y += *inter.yOld(level);

      DEBUG_PRINT("EulerMoreauOSI::updateOutput : ResiduY() \n");
      SiconosVector& residuY = *inter_work[EulerMoreauOSI::VEC_RESIDU_Y];
      DEBUG_EXPR(residuY.display());

      y -= residuY;
      DEBUG_PRINT("EulerMoreauOSI::updateOutput : y(level) \n");
      DEBUG_EXPR(y.display());

      BlockVector& deltax = *inter_work_block[EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("EulerMoreauOSI::updateOutput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (r.C())
        prod(*r.C(), deltax, y, false);
      else
        prod(*relationMat[FirstOrderR::mat_C], deltax, y, false);


      DEBUG_PRINT("EulerMoreauOSI::updateOutput : y before osnsM\n");
      DEBUG_EXPR(y.display());
      if (interProp.block)
      {
        SiconosMatrix& osnsM = *interProp.block;
        prod(osnsM, *inter.lambda(level), y, false);
        DEBUG_EXPR(inter.lambda(level)->display());
        DEBUG_EXPR(osnsM.display());
        DEBUG_PRINT("EulerMoreauOSI::updateOutput : new linearized y \n");
        DEBUG_EXPR(y.display());
      }

      SiconosVector& x = *inter_work[EulerMoreauOSI::VEC_X];
      x = *DSlink[FirstOrderR::x];


      SiconosVector& hAlpha= *inter_work[EulerMoreauOSI::H_ALPHA];

      r.computeh(time, x, *inter.lambda(level), hAlpha);
      DEBUG_PRINT("EulerMoreauOSI::updateOutput : new Halpha \n");
      DEBUG_EXPR(hAlpha.display());
    }
    else if (relationSubType == NonLinearR )
    {
      FirstOrderNonLinearR & r = static_cast<FirstOrderNonLinearR&>(*inter.relation());
      // compute the new y  obtained by linearisation (see DevNotes)
      // y_{alpha+1}_{k+1} = h(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      // or equivalently
      // y_{alpha+1}_{k+1} = y_{alpha}_{k+1} - ResiduY_{k+1}^{alpha}
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      SiconosVector& y = *inter.y(level);
      DEBUG_EXPR(y.display());


      if (r.D())
        prod(*r.D(), *(inter.lambdaOld(level)), y, true);
      else
        prod(*relationMat[FirstOrderR::mat_D], *(inter.lambdaOld(level)), y, true);

      y *= -1.0;
      //SiconosVector yOld = *inter.yOld(0); // Retrieve  y_{alpha}_{k+1}
      DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : yOld(level) \n");
      DEBUG_EXPR(inter.yOld(level)->display());

      y += *inter.yOld(level);

      DEBUG_PRINT("EulerMoreauOSI::updateOutput : ResiduY() \n");
      SiconosVector& residuY = *inter_work[EulerMoreauOSI::VEC_RESIDU_Y];
      DEBUG_EXPR(residuY.display());

      y -= residuY;

      DEBUG_PRINT("EulerMoreauOSI::updateOutput : y(level) \n");
      DEBUG_EXPR(y.display());

      BlockVector& deltax = *inter_work_block[EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("EulerMoreauOSI::updateOutput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (r.C())
        prod(*r.C(), deltax, y, false);
      else
        prod(*relationMat[FirstOrderR::mat_C], deltax, y, false);

      if (interProp.block)
      {
        SiconosMatrix& osnsM = *interProp.block;
        // osnsM = h * C * W^-1 * B + D
        DEBUG_EXPR(osnsM.display(););
        prod(osnsM, *inter.lambda(level), y, false);
      }
      DEBUG_PRINT("EulerMoreauOSI::updateOutput : new linearized y \n");
      DEBUG_EXPR(y.display());

      SiconosVector& x = *inter_work[EulerMoreauOSI::VEC_X];
      x = *DSlink[FirstOrderR::x];
      SiconosVector& z = *inter_work[EulerMoreauOSI::VEC_Z];
      z = *DSlink[FirstOrderR::z];

      SiconosVector& hAlpha =  *inter_work[EulerMoreauOSI::H_ALPHA];
      r.computeh(time, x, *inter.lambda(level), z, hAlpha);
      DEBUG_EXPR(x.display(););
      DEBUG_PRINT("EulerMoreauOSI::updateOutput : new Halpha \n");
      DEBUG_EXPR(hAlpha.display());
      *DSlink[FirstOrderR::z] = z;
    }
    else
      inter.computeOutput(time, level);
  }
  DEBUG_END("EulerMoreauOSI::updateOutput(double time, unsigned int level)\n");
}

void EulerMoreauOSI::updateInput(double time, unsigned int level)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  //_simulation->nonSmoothDynamicalSystem()->updateInput(time,level);

  // Set dynamical systems non-smooth part to zero.
  _simulation->nonSmoothDynamicalSystem()->reset(level);


  InteractionsGraph::VIterator ui, uiend;

  SP::InteractionsGraph indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = * indexSet0->bundle(*ui);
    assert(inter.lowerLevelForInput() <= level);
    assert(inter.upperLevelForInput() >= level);

    VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
    VectorOfSMatrices& relationMat = inter.relationMatrices();

    InteractionProperties& interProp = indexSet0->properties(*ui);
    VectorOfVectors& inter_work = *interProp.workVectors;
    VectorOfSMatrices& inter_work_mat = *interProp.workMatrices;
    VectorOfBlockVectors& inter_work_block = *interProp.workBlockVectors;

    RELATION::SUBTYPES relationSubType = inter.relation()->getSubType();
    if (relationSubType == Type2R)
    {
      FirstOrderType2R & r = static_cast<FirstOrderType2R&>(*inter.relation());
      SiconosVector lambda = *inter.lambda(level);
      lambda -= *(inter.lambdaOld(level));

      if (r.B())
        prod(*r.B(), lambda, *inter_work[EulerMoreauOSI::G_ALPHA], false);
      else
        prod(*relationMat[FirstOrderR::mat_B], lambda, *inter_work[EulerMoreauOSI::G_ALPHA], false);


      *DSlink[FirstOrderR::r] += *inter_work[EulerMoreauOSI::G_ALPHA];
      DEBUG_EXPR(DSlink[FirstOrderR::r]->display(););
      //compute the new g_alpha
      r.computeg(time, *inter.lambda(level), *inter_work[EulerMoreauOSI::G_ALPHA]);
      DEBUG_EXPR(inter_work[EulerMoreauOSI::G_ALPHA]->display(););
    }
    else if (relationSubType == NonLinearR )
    {
      FirstOrderNonLinearR & r = static_cast<FirstOrderNonLinearR&>(*inter.relation());
      // compute the new r  obtained by linearisation
      // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
      //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )


      SiconosVector lambda = *inter.lambda(level);
      lambda -= *(inter.lambdaOld(level));

      SiconosVector& g_alpha = *inter_work[EulerMoreauOSI::G_ALPHA];

      if (r.B())
        prod(*r.B(), lambda, g_alpha, false);
      else
        prod(*relationMat[FirstOrderR::mat_B], lambda, g_alpha, false);

      BlockVector& deltax = *inter_work_block[EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("FirstOrderNonLinearR::computeInput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (r.K())
        prod(*r.K(), deltax, g_alpha, false);
      else
        prod(*relationMat[FirstOrderR::mat_K], deltax, g_alpha, false);


      // Khat = h * K * W^-1 * B
      prod(*inter_work_mat[EulerMoreauOSI::MAT_KHAT], *inter.lambda(level), g_alpha, false);

      *DSlink[FirstOrderR::r] += g_alpha;

      //compute the new g_alpha
      SiconosVector& x = *inter_work[EulerMoreauOSI::VEC_X];
      x = *DSlink[FirstOrderR::x];
      SiconosVector& z = *inter_work[EulerMoreauOSI::VEC_Z];
      z = *DSlink[FirstOrderR::z];
      r.computeg(time, x, *inter.lambda(level), z, g_alpha);
      *DSlink[FirstOrderR::z] = z;

    }
    else
    {
      inter.computeInput(time, level);
    }
  }

}


double EulerMoreauOSI::computeResiduOutput(double time, SP::InteractionsGraph indexSet)
{

  double residu =0.0;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    VectorOfVectors& inter_work = *indexSet->properties(*ui).workVectors;
    SiconosVector&  residuY = *inter_work[EulerMoreauOSI::VEC_RESIDU_Y];
    Interaction & inter = *indexSet->bundle(*ui);
    residuY = *inter_work[EulerMoreauOSI::H_ALPHA];
    scal(-1, residuY, residuY);
    residuY += *(inter.y(0));
    DEBUG_EXPR(residuY.display(););
    residu = std::max(residu,residuY.norm2());
  }
  return residu;
}
double EulerMoreauOSI::computeResiduInput(double time, SP::InteractionsGraph indexSet)
{
  double residu =0.0;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    InteractionProperties& interProp = indexSet->properties(*ui);
    VectorOfVectors& inter_work = *interProp.workVectors;
    SP::Interaction inter = indexSet->bundle(*ui);
    VectorOfBlockVectors& DSlink = inter->linkToDSVariables();
    SiconosVector&  residuR = *inter_work[EulerMoreauOSI::VEC_RESIDU_R];
    //Residu_r = r_alpha_k+1 - g_alpha;
    residuR = *DSlink[FirstOrderR::r];
    residuR -= *inter_work[EulerMoreauOSI::G_ALPHA];
    DEBUG_EXPR(residuR.display(););
    residu = std::max(residu,residuR.norm2());
  }
  return residu;
}
