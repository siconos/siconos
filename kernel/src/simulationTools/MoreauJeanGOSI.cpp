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
#include "MoreauJeanGOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "CxxStd.hpp"

#include "TypeName.hpp"

#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"


// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>


using namespace RELATION;

// --- constructor from a set of data ---
MoreauJeanGOSI::MoreauJeanGOSI(double theta, double gamma):
  OneStepIntegrator(OSI::MOREAUJEANGOSI), _useGammaForRelation(false), _explicitNewtonEulerDSOperators(false)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =1;
  _levelMinForInput =0;
  _levelMaxForInput =1;
  _steps=1;
  _theta = theta;
  if(!isnan(gamma))
  {
    _gamma = gamma;
    _useGamma = true;
  }
  else
  {
    _gamma = 1.0;
    _useGamma = false;
  }
}
void MoreauJeanGOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type
  Type::Siconos dsType = Type::value(*ds);

  // Compute W (iteration matrix)
  initializeIterationMatrixW(t, ds);
  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
    {
      SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS> (ds);

      workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
      workVectors[OneStepIntegrator::residu_free].reset(new SiconosVector(lds->dimension()));
      workVectors[OneStepIntegrator::free].reset(new SiconosVector(lds->dimension()));
      workVectors[OneStepIntegrator::local_buffer].reset(new SiconosVector(lds->dimension()));

      lds->computeForces(m.t0(), lds->q(), lds->velocity());
      lds->swapInMemory();
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS neds = std11::static_pointer_cast<NewtonEulerDS> (ds);
      workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
      workVectors[OneStepIntegrator::residu_free].reset(new SiconosVector(neds->dimension()));
      workVectors[OneStepIntegrator::free].reset(new SiconosVector(neds->dimension()));
      //Compute a first value of the dotq  to store it in  _dotqMemory
      SP::SiconosMatrix T = neds->T();
      SP::SiconosVector dotq = neds->dotq();
      SP::SiconosVector v = neds->twist();
      prod(*T, *v, *dotq, true);

      //Compute a first value of the forces to store it in _forcesMemory
      neds->computeForces(m.t0(), neds->q(), v);
      neds->swapInMemory();
    }
}

void MoreauJeanGOSI::fillDSLinks(Interaction &inter,
				   InteractionProperties& interProp,
				   DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);
  VectorOfVectors& workV = *interProp.workVectors;
  workV.resize(MoreauJeanGOSI::WORK_INTERACTION_LENGTH);
  workV[MoreauJeanGOSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));
  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  Relation &relation =  *inter.relation();  
  RELATION::TYPES relationType = relation.getType();

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  // if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  // {
  //   RuntimeException::selfThrow("MoreauJeanGOSI::fillDSLinks. The implementation is not correct for two different OSI for one interaction");
  // }



  
  /* allocate and set work vectors for the osi */
  unsigned int xfree =0;
  if (relationType == Lagrangian)
  {
    xfree = LagrangianR::xfree;
  }
  else if (relationType == NewtonEuler)
  {
    xfree = NewtonEulerR::xfree;
  }
  if (ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!DSlink[xfree]) || (DSlink[xfree]->numberOfBlocks() !=2 ))
      DSlink[xfree].reset(new BlockVector(2));
  }
  else
  {
    if ((!DSlink[xfree]) || (DSlink[xfree]->numberOfBlocks() !=1 ))
      DSlink[xfree].reset(new BlockVector(1));
  }

  if(checkOSI(DSG.descriptor(ds1)))
  {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    DSlink[xfree]->setVectorPtr(0,workVds1[OneStepIntegrator::free]);
  }
  DEBUG_PRINTF("ds1->number() %i\n",ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n",ds2->number());


  if (ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");
    if(checkOSI(DSG.descriptor(ds2)))
    {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n",ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      DSlink[xfree]->setVectorPtr(1,workVds2[OneStepIntegrator::free]);
    }
  }
  
}

void MoreauJeanGOSI::initialize_nonsmooth_problems()
{
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->setIndexSetLevel(1);
  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->setInputOutputLevel(1);
  //  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->initialize(_simulation);
}

void MoreauJeanGOSI::initializeIterationMatrixW(double t, SP::DynamicalSystem ds)
{
  DEBUG_PRINT("MoreauJeanGOSI::initializeIterationMatrixW starts\n");
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending on its type.

  if(!ds)
    RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixW(t, ds, dsv) - ds == NULL or dsv == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  
  if(_dynamicalSystemsGraph->properties(dsv).W)
    RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixW(t,ds) - W(ds) is already in the map and has been initialized.");

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);
  unsigned int sizeW = ds->dimension();
  if(dsType == Type::LagrangianDS)
  {
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
    if(d->mass())
    {
      d->computeMass(d->q());
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->mass())); //*W = *d->mass();
    }
    else
    {
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(sizeW, sizeW));
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }
    // Compute the W matrix
    computeW(t, ds, *_dynamicalSystemsGraph->properties(dsv).W);
    // WBoundaryConditions initialization
    if(d->boundaryConditions())
      initializeIterationMatrixWBoundaryConditions(d);
  }
  // 2 - Lagrangian linear systems
  else if(dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);
    _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->mass()));
    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(dsv).W;

    SP::SiconosMatrix K = d->K();
    SP::SiconosMatrix C = d->C();

    if(C)
      scal(h * _theta, *C, W, false); // W += h*_theta *C
    if(K)
      scal(h * h * _theta * _theta, *K, W, false); // W = h*h*_theta*_theta*K

    // WBoundaryConditions initialization
    if(d->boundaryConditions())
      initializeIterationMatrixWBoundaryConditions(d);
  }

  // === ===
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
    _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->mass()));

    computeW(t,ds, *_dynamicalSystemsGraph->properties(dsv).W);
    // WBoundaryConditions initialization
    if(d->boundaryConditions())
      initializeIterationMatrixWBoundaryConditions(d);

  }
  else RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixW - not yet implemented for Dynamical system of type : " + Type::name(*ds));

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
  DEBUG_PRINT("MoreauJeanGOSI::initializeIterationMatrixW ends\n");


}


void MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  DEBUG_PRINT("MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds) starts\n");
  if(!ds)
    RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - ds does not belong to the OSI.");

  Type::Siconos dsType = Type::value(*ds);
  unsigned int dsN = ds->number();

  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS || dsType == Type::NewtonEulerDS)
  {
    if(_WBoundaryConditionsMap.find(dsN) != _WBoundaryConditionsMap.end())
      RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - WBoundaryConditions(ds) is already in the map and has been initialized.");
    
    // Memory allocation for WBoundaryConditions
    unsigned int sizeWBoundaryConditions = ds->dimension(); // n for first order systems, ndof for lagrangian.

    SP::BoundaryCondition bc;
    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      bc = d->boundaryConditions();
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      bc = d->boundaryConditions();
    }
    unsigned int numberBoundaryConditions = bc->velocityIndices()->size();
    _WBoundaryConditionsMap[dsN].reset(new SimpleMatrix(sizeWBoundaryConditions, numberBoundaryConditions));
    computeWBoundaryConditions(ds);
  }
  else
    RuntimeException::selfThrow("MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions - not yet implemented for Dynamical system of type :" +  Type::name(*ds));
  DEBUG_PRINT("MoreauJeanGOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds) ends \n");
}


void MoreauJeanGOSI::computeWBoundaryConditions(SP::DynamicalSystem ds)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  assert(ds &&
         "MoreauJeanGOSI::computeWBoundaryConditions(t,ds) - ds == NULL");

  Type::Siconos dsType = Type::value(*ds);
  unsigned int dsN = ds->number();
  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS ||  dsType == Type::NewtonEulerDS)
  {
    assert((_WBoundaryConditionsMap.find(dsN) != _WBoundaryConditionsMap.end()) &&
           "MoreauJeanGOSI::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

    SP::SimpleMatrix WBoundaryConditions = _WBoundaryConditionsMap[dsN];

    SP::SiconosVector columntmp(new SiconosVector(ds->dimension()));

    int columnindex = 0;

    std::vector<unsigned int>::iterator itindex;

    SP::BoundaryCondition bc;
    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&>(*ds);
      bc = d.boundaryConditions();
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      NewtonEulerDS& d = static_cast<NewtonEulerDS&>(*ds);
      bc = d.boundaryConditions();
    }
    SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;

    for(itindex = bc->velocityIndices()->begin() ;
        itindex != bc->velocityIndices()->end();
        ++itindex)
    {

      W->getCol(*itindex, *columntmp);
      /*\warning we assume that W is symmetric in the Lagrangian case
        we store only the column and not the row */

      WBoundaryConditions->setCol(columnindex, *columntmp);
      double diag = (*columntmp)(*itindex);
      columntmp->zero();
      (*columntmp)(*itindex) = diag;

      W->setCol(*itindex, *columntmp);
      W->setRow(*itindex, *columntmp);


      columnindex ++;
    }
    DEBUG_EXPR(W->display());
  }
  else
    RuntimeException::selfThrow("MoreauJeanGOSI::computeWBoundaryConditions - not yet implemented for Dynamical system type : " +  Type::name(*ds));
}


void MoreauJeanGOSI::computeW(double t, SP::DynamicalSystem ds, SiconosMatrix& W)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.
  DEBUG_PRINT("MoreauJeanGOSI::computeW starts\n");

  assert(ds &&
         "MoreauJeanGOSI::computeW(t,ds) - ds == NULL");

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  if(dsType == Type::LagrangianLinearTIDS)
  {
    // Nothing: W does not depend on time.
  }
  else if(dsType == Type::LagrangianDS)
  {

    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
    SP::SiconosMatrix K = d->jacobianqForces(); // jacobian according to q
    SP::SiconosMatrix C = d->jacobianqDotForces(); // jacobian according to velocity

    if(d->mass())
      {
	d->computeMass(d->q());
	W = *d->mass();
      }
    else
      W.zero();

    if(C)
    {
      d->computeJacobianqDotForces(t);
      scal(-h * _theta, *C, W, false); // W -= h*_theta*C
    }

    if(K)
    {
      d->computeJacobianqForces(t);
      scal(-h * h * _theta * _theta, *K, W, false); //*W -= h*h*_theta*_theta**K;
    }
  }
  // === ===
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
    W = *(d->mass());

    SP::SiconosMatrix K = d->jacobianqForces(); // jacobian according to q
    SP::SiconosMatrix C = d->jacobianvForces(); // jacobian according to velocity

    if(C)
    {
      d->computeJacobianvForces(t);
      scal(-h * _theta, *C, W, false); // W -= h*_theta*C
    }
    if(K)
    {
      d->computeJacobianqForces(t);
      SP::SiconosMatrix T = d->T();
      DEBUG_EXPR(T->display(););
      DEBUG_EXPR(K->display(););
      SP::SimpleMatrix  buffer(new SimpleMatrix(*(d->mass())));
      prod(*K, *T, *buffer, true);
      scal(-h * h * _theta * _theta, *buffer, W, false);
      //*W -= h*h*_theta*_theta**K;
    }
  }
  else RuntimeException::selfThrow("MoreauJeanGOSI::computeW - not yet implemented for Dynamical system of type : " +Type::name(*ds));
  DEBUG_PRINT("MoreauJeanGOSI::computeW ends\n");
  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}

void MoreauJeanGOSI::computeInitialNewtonState()
{
  DEBUG_PRINT("MoreauJeanGOSI::computeInitialNewtonState() starts\n");
  // Compute the position value giving the initial velocity.
  // The goal of to save one newton iteration for nearly linear system
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem &ds = *_dynamicalSystemsGraph->bundle(*dsi);

    if(_explicitNewtonEulerDSOperators)
    {
      if(Type::value(ds) == Type::NewtonEulerDS)
      {
        // The goal is to update T() one time at the beginning of the Newton Loop
        // We want to be explicit on this function since we do not compute their Jacobians.
        NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
        SP::SiconosVector qold = d.qMemory()->getSiconosVector(0);
        //SP::SiconosVector q = d.q();
        computeT(qold,d.T());
      }
    }
    // The goal is to converge in one iteration of the system is almost linear
    // we start the Newton loop q = q0+hv0
    updatePosition(ds);


  }
  DEBUG_PRINT("MoreauJeanGOSI::computeInitialNewtonState() ends\n");
}



double MoreauJeanGOSI::computeResidu()
{
  DEBUG_PRINT("\nMoreauJeanGOSI::computeResidu(), start\n");
  // This function is used to compute the residu for each "MoreauJeanGOSI-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
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
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds); // Its type

    // 3 - Lagrangian Non Linear Systems
    if(dsType == Type::LagrangianDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::computeResidu(), dsType == Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1) - h*(1-theta)*forces(ti,vi,qi) - p_i+1
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      SiconosVector& free = *workVectors[OneStepIntegrator::free];

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);
      SP::SiconosVector q = d->q();

      SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      //residuFree.zero();
      DEBUG_EXPR(residuFree.display());

      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());
      DEBUG_EXPR(q->display());
      DEBUG_EXPR(v->display());

      residuFree = *v;
      sub(residuFree, *vold, residuFree);
      if(d->mass())
      {
        d->computeMass(d->q());
        prod(*(d->mass()), residuFree, residuFree); // residuFree = M(v - vold)
      }

      if(d->forces())
      {
        // Cheaper version: get forces(ti,vi,qi) from memory
        SP::SiconosVector fold = d->forcesMemory()->getSiconosVector(0);
        double coef = -h * (1 - _theta);
        scal(coef, *fold, residuFree, false);

        // Expensive computes forces(ti,vi,qi)
        // d->computeForces(told, qold, vold);
        // double coef = -h * (1 - _theta);
        // // residuFree += coef * fL_i
        // scal(coef, *d->forces(), *residuFree, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        d->computeForces(t,q,v);
        coef = -h * _theta;
        scal(coef, *d->forces(), residuFree, false);

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        //SP::SiconosVector qbasedonv(new SiconosVector(*qold));
        //*qbasedonv +=  h * ((1 - _theta)* *vold + _theta * *v);
        //d->computeForces(t, qbasedonv, v);
        //coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        //scal(coef, *d->forces(), residuFree, false);


      }

      if(d->boundaryConditions())
      {
        d->boundaryConditions()->computePrescribedVelocity(t);

        unsigned int columnindex = 0;
        SP::SimpleMatrix WBoundaryConditions = _WBoundaryConditionsMap[ds->number()];
        SP::SiconosVector columntmp(new SiconosVector(ds->dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d->boundaryConditions()->velocityIndices()->begin() ;
            itindex != d->boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          double DeltaPrescribedVelocity =
            d->boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - v->getValue(*itindex);

          WBoundaryConditions->getCol(columnindex, *columntmp);
          residuFree -= *columntmp * (DeltaPrescribedVelocity);

          residuFree.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;
        }
      }

      free = residuFree; // copy residuFree in Workfree

      //       std::cout << "MoreauJeanGOSI::ComputeResidu LagrangianDS residufree :"  << std::endl;
      DEBUG_EXPR(residuFree.display());

      if(d->p(1))
        free -= *d->p(1); // Compute Residu in Workfree Notation !!
      // We use DynamicalSystem::free as tmp buffer

      if(d->boundaryConditions())
      {
        unsigned int columnindex = 0;
        SP::SimpleMatrix WBoundaryConditions = _WBoundaryConditionsMap[ds->number()];
        SP::SiconosVector columntmp(new SiconosVector(ds->dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d->boundaryConditions()->velocityIndices()->begin() ;
            itindex != d->boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          double DeltaPrescribedVelocity =
            d->boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - v->getValue(*itindex);

          WBoundaryConditions->getCol(columnindex, *columntmp);

          free.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;
        }
      }


      DEBUG_EXPR(free.display());
      normResidu = free.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }
    // 4 - Lagrangian Linear Systems
    else if(dsType == Type::LagrangianLinearTIDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::computeResidu(), dsType == Type::LagrangianLinearTIDS\n");
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual for v = v_i
      // otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); //vi

      DEBUG_EXPR(qold->display(););
      DEBUG_EXPR(vold->display(););
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););

      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      SiconosVector& free_rhs = *workVectors[OneStepIntegrator::free];
      // --- ResiduFree computation Equation (1) ---
      residuFree.zero();
      SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      prod(W, *vold, residuFree);

      double coeff;
      // -- No need to update W --

      SP::SiconosVector v = d->velocity(); // v = v_k,i+1

      SP::SiconosMatrix C = d->C();
      if(C)
        prod(h, *C, *vold, residuFree, false); // vfree += -h*C*vi

      SP::SiconosMatrix K = d->K();
      if(K)
      {
        coeff = -h * h * _theta;
        prod(coeff, *K, *vold, residuFree, false); // vfree += -h^2*_theta*K*vi
        prod(-h, *K, *qold, residuFree, false); // vfree += -h*K*qi
      }

      SP::SiconosVector Fext = d->fExt();
      if(Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - _theta);
        scal(coeff, *(d->fExt()), residuFree, false); // vfree += h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * _theta;
        scal(coeff, *(d->fExt()), residuFree, false); // vfree += h*_theta * fext(ti+1)
      }


      // Computation of the complete residual Equation (2)
      //   ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      //       SP::SiconosMatrix M = d->mass();
      //       SP::SiconosVector realresiduFree (new SiconosVector(*residuFree));
      //       realresiduFree->zero();
      //       prod(*M, (*v-*vold), *realresiduFree); // residuFree = M(v - vold)
      //       SP::SiconosVector qkplustheta (new SiconosVector(*qold));
      //       qkplustheta->zero();
      //       *qkplustheta = *qold + h *((1-_theta)* *vold + _theta* *v);
      //       if (C){
      //         double coef = h*(1-_theta);
      //         prod(coef, *C, *vold , *realresiduFree, false);
      //         coef = h*(_theta);
      //         prod(coef,*C, *v , *realresiduFree, false);
      //       }
      //       if (K){
      //         double coef = h*(1-_theta);
      //         prod(coef,*K , *qold , *realresiduFree, false);
      //         coef = h*(_theta);
      //         prod(coef,*K , *qkplustheta , *realresiduFree, false);
      //       }

      //       if (Fext)
      //       {
      //         // computes Fext(ti)
      //         d->computeFExt(told);
      //         coeff = -h*(1-_theta);
      //         scal(coeff, *Fext, *realresiduFree, false); // vfree -= h*(1-_theta) * fext(ti)
      //         // computes Fext(ti+1)
      //         d->computeFExt(t);
      //         coeff = -h*_theta;
      //         scal(coeff, *Fext, *realresiduFree, false); // vfree -= h*_theta * fext(ti+1)
      //       }




      if(d->boundaryConditions())
      {
         RuntimeException::selfThrow("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }

      free_rhs = residuFree; // copy residuFree in Workfree
      DEBUG_EXPR(free_rhs.display());
      // if(d->p(1))
      //   free -= *d->p(1); // Compute Residu in Workfree Notation !!
      // We use DynamicalSystem::free as tmp buffer

      normResidu = 0.0; // we assume that v = vfree + W^(-1) p
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::computeResidu(), dsType == Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) - h*(1-_theta)*forces(ti,vi,qi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      SiconosVector& free = *workVectors[OneStepIntegrator::free];
      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0);


      // Get the current state vector
      SP::SiconosVector q = d->q();
      SP::SiconosVector v = d->twist(); // v = v_k,i+1

      // Get the (constant mass matrix)
      SP::SiconosMatrix massMatrix = d->mass();
      prod(*massMatrix, (*v - *vold), residuFree, true); // residuFree = M(v - vold)
      DEBUG_EXPR(residuFree.display(););

      if(d->forces())   // if fL exists
      {
        DEBUG_PRINTF("MoreauJeanGOSI:: _theta = %e\n",_theta);
        DEBUG_PRINTF("MoreauJeanGOSI:: h = %e\n",h);

        // Cheaper version: get forces(ti,vi,qi) from memory
        SP::SiconosVector fold = d->forcesMemory()->getSiconosVector(0);
        double coef = -h * (1 - _theta);
        scal(coef, *fold, residuFree, false);

        // Expensive version to check ...
        //d->computeForces(told,qold,vold);
        //double coef = -h * (1.0 - _theta);
        //scal(coef, *d->forces(), residuFree, false);

        DEBUG_PRINT("MoreauJeanGOSI:: old forces :\n");
        DEBUG_EXPR(d->forces()->display(););
        DEBUG_EXPR(residuFree.display(););

        // computes forces(ti,v,q)
        d->computeForces(t,q,v);
        coef = -h * _theta;
        scal(coef, *d->forces(), residuFree, false);
        DEBUG_PRINT("MoreauJeanGOSI:: new forces :\n");
        DEBUG_EXPR(d->forces()->display(););
        DEBUG_EXPR(residuFree.display(););

      }


      if(d->boundaryConditions())
      {
         RuntimeException::selfThrow("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }

      free = residuFree;
      if(d->p(1))
        free -= *d->p(1);// We use DynamicalSystem::free as tmp buffer


      if(d->boundaryConditions())
      {
        RuntimeException::selfThrow("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }

      DEBUG_PRINT("MoreauJeanGOSI::computeResidu :\n");
      DEBUG_EXPR(residuFree.display(););
      DEBUG_EXPR(if(d->p(1)) d->p(1)->display(););
      DEBUG_EXPR(free.display(););

      normResidu = free.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }
    else
      RuntimeException::selfThrow("MoreauJeanGOSI::computeResidu - not yet implemented for Dynamical system of type: " + Type::name(*ds));

    if(normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void MoreauJeanGOSI::computeFreeState()
{
  DEBUG_PRINT("\nMoreauJeanGOSI::computeFreeState() starts\n");

  DEBUG_PRINT("MoreauJeanGOSI::computeFreeState() ends\n");
}

void MoreauJeanGOSI::prepareNewtonIteration(double time)
{
  DEBUG_BEGIN(" MoreauJeanOSI::prepareNewtonIteration(double time)\n");
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
    computeW(time, ds, W);
  }

  if(!_explicitNewtonEulerDSOperators)
  {
    DynamicalSystemsGraph::VIterator dsi, dsend;

    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;

      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

      //  VA <2016-04-19 Tue> We compute T to be consistent with the Jacobian at the beginning of the Newton iteration and not at the end
      Type::Siconos dsType = Type::value(*ds);
      if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        computeT(d->q(),d->T());
      }
    }

  }


  DEBUG_END(" MoreauJeanOSI::prepareNewtonIteration(double time)\n");

}


struct MoreauJeanGOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem& _osnsp;
  Interaction& _inter;
  InteractionProperties& _interProp;

  _NSLEffectOnFreeOutput(OneStepNSProblem& p, Interaction& inter, InteractionProperties& interProp) :
    _osnsp(p), _inter(inter), _interProp(interProp) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter.nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[MoreauJeanGOSI::OSNSP_RHS];
    subscal(e, *_inter.y_k(_osnsp.inputOutputLevel()), osnsp_rhs, subCoord, true);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    DEBUG_PRINTF("e= %e\n", e)
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[MoreauJeanGOSI::OSNSP_RHS];
    DEBUG_PRINTF("y_k = %e\n", (*_inter.y_k(_osnsp.inputOutputLevel()))(0));
    DEBUG_PRINTF("level = %i\n", _osnsp.inputOutputLevel() );

    osnsp_rhs(0) =  e * (*_inter.y_k(_osnsp.inputOutputLevel()))(0);

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

void MoreauJeanGOSI::NSLcontrib(SP::Interaction inter, OneStepNSProblem& osnsp)
{
  if(inter->relation()->getType() == Lagrangian || inter->relation()->getType() == NewtonEuler)
  {
    InteractionsGraph& indexSet = *osnsp.simulation()->indexSet(osnsp.indexSetLevel());
    InteractionsGraph::VDescriptor ivd = indexSet.descriptor(inter);
    _NSLEffectOnFreeOutput nslEffectOnFreeOutput = _NSLEffectOnFreeOutput(osnsp, *inter, indexSet.properties(ivd) );
    inter->nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
}

void MoreauJeanGOSI::integrate(double& tinit, double& tend, double& tout, int& notUsed)
{
}

void MoreauJeanGOSI::updatePosition(DynamicalSystem& ds)
{
  DEBUG_END("MoreauJeanGOSI::updatePosition(const unsigned int )\n");
  double h = _simulation->timeStep();

  Type::Siconos dsType = Type::value(ds);

  // 1 - Lagrangian Systems
  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    // get dynamical system

    LagrangianDS& d = static_cast<LagrangianDS&> (ds);

    // Compute q
    SP::SiconosVector v = d.velocity();
    SP::SiconosVector q = d.q();
    DEBUG_EXPR(v->display());
    DEBUG_EXPR(q->display());

    //  -> get previous time step state
    SP::SiconosVector vold = d.velocityMemory()->getSiconosVector(0);
    SP::SiconosVector qold = d.qMemory()->getSiconosVector(0);
    // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    double coeff = h * _theta;
    scal(coeff, *v, *q) ; // q = h*theta*v
    coeff = h * (1 - _theta);
    scal(coeff, *vold, *q, false); // q += h(1-theta)*vold
    *q += *qold;
    DEBUG_EXPR(v->display());
    DEBUG_EXPR(q->display());

  }
  else if(dsType == Type::NewtonEulerDS)
  {
    // get dynamical system
    NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
    SP::SiconosVector v = d.twist();
    DEBUG_PRINT("MoreauJeanGOSI::updateState()\n ")
    DEBUG_EXPR(d.display());
    DEBUG_PRINT("MoreauJeanGOSI::updateState() prev v\n")
    DEBUG_EXPR(v->display());

    //compute q
    //first step consists in computing  \dot q.
    //second step consists in updating q.
    //
    SP::SiconosMatrix T = d.T();
    SP::SiconosVector dotq = d.dotq();
    prod(*T, *v, *dotq, true);

    DEBUG_PRINT("MoreauJeanGOSI::updateState v\n");
    DEBUG_EXPR(v->display());
    DEBUG_EXPR(dotq->display());

    SP::SiconosVector q = d.q();

    //  -> get previous time step state
    SP::SiconosVector dotqold = d.dotqMemory()->getSiconosVector(0);
    SP::SiconosVector qold = d.qMemory()->getSiconosVector(0);
    // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    double coeff = h * _theta;
    scal(coeff, *dotq, *q) ; // q = h*theta*v
    coeff = h * (1 - _theta);
    scal(coeff, *dotqold, *q, false); // q += h(1-theta)*vold
    *q += *qold;
    DEBUG_PRINT("new q before normalizing\n");
    DEBUG_EXPR(q->display());

    //q[3:6] must be normalized
    d.normalizeq();

    /* \warning VA 02/06/2013.
     * What is the reason of doing the following computation ?
     */
    // dotq->setValue(3, (q->getValue(3) - qold->getValue(3)) / h);
    // dotq->setValue(4, (q->getValue(4) - qold->getValue(4)) / h);
    // dotq->setValue(5, (q->getValue(5) - qold->getValue(5)) / h);
    // dotq->setValue(6, (q->getValue(6) - qold->getValue(6)) / h);

    // d->computeT(); //  VA 09/06/2015. We prefer only compute T() every time--step for Newton convergence reasons.

  }
  DEBUG_END("MoreauJeanGOSI::updatePosition(const unsigned int )\n");

}

void MoreauJeanGOSI::updateState(const unsigned int )
{

  DEBUG_PRINT("MoreauJeanGOSI::updateState(const unsigned int )\n");

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem& ds = *_dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // Get the DS type

    Type::Siconos dsType = Type::value(ds);

    // 3 - Lagrangian Systems
    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);
      bool baux = dsType == Type::LagrangianDS && useRCC && _simulation->relativeConvergenceCriterionHeld();

      SiconosVector &q = *d.q();
      SiconosVector& local_buffer = *workVectors[OneStepIntegrator::local_buffer];

      // Save value of q in stateTmp for future convergence computation
      if(baux)
        local_buffer = q;

      updatePosition(ds);

      if(baux)
      {
        double ds_norm_ref = 1. + ds.x0()->norm2(); // Should we save this in the graph?
        local_buffer -= q;
        double aux = (local_buffer.norm2()) / ds_norm_ref;
        if(aux > RelativeTol)
          _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::updateState(const unsigned int ), dsType == Type::NewtonEulerDS \n");
      updatePosition(ds);
    }
    else RuntimeException::selfThrow("MoreauJeanGOSI::updateState - not yet implemented for Dynamical system of type: " +  Type::name(ds));

  }
}


bool MoreauJeanGOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("addInteractionInIndexSet(SP::Interaction inter, unsigned int i)\n");

  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity

  double gamma = 1.0 / 2.0;
  if(_useGamma)
  {
    gamma = _gamma;
  }
  DEBUG_PRINTF("MoreauJeanGOSI::addInteractionInIndexSet of level = %i y=%e, yDot=%e, y_estimated=%e.\n", i,  y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;

  DEBUG_PRINTF("y = %e\n", y);
  assert(!isnan(y));
  DEBUG_EXPR_WE(
    if(y <= 0)
    {
      DEBUG_PRINT("MoreauJeanGOSI::addInteractionInIndexSet ACTIVATED.\n");
    }
    else
    {
      DEBUG_PRINT("MoreauJeanGOSI::addInteractionInIndexSet NOT ACTIVATED.\n");
    }
    );
  return (y <= 0.0);
}


bool MoreauJeanGOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  if(_useGamma)
  {
    gamma = _gamma;
  }
  DEBUG_PRINTF("MoreauJeanGOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));

  DEBUG_EXPR(
    if(y > 0)
    DEBUG_PRINT("MoreauJeanGOSI::removeInteractionInIndexSet DEACTIVATE.\n");
  );
  return (y > 0.0);
}



void MoreauJeanGOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanOSI OSI display ======" <<std::endl;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  if(_dynamicalSystemsGraph)
  {
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
  }
  std::cout << "================================" <<std::endl;
}
