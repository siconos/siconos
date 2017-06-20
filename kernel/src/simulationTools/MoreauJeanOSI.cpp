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
#include "MoreauJeanOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianLinearDiagonalDS.hpp"

#include "FirstOrderR.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianCompliantLinearTIR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "CxxStd.hpp"

#include "TypeName.hpp"

#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include <debug.h>


using namespace RELATION;

// --- constructor from a set of data ---
MoreauJeanOSI::MoreauJeanOSI(double theta, double gamma):
  OneStepIntegrator(OSI::MOREAUJEANOSI), _useGammaForRelation(false),_explicitNewtonEulerDSOperators(false)
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

const SimpleMatrix MoreauJeanOSI::getW(SP::DynamicalSystem ds)
{
  assert(ds &&
         "MoreauJeanOSI::getW(ds): ds == NULL.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "MoreauJeanOSI::getW(ds): W[ds] == NULL.");
  return *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W; // Copy !!
}

SP::SimpleMatrix MoreauJeanOSI::W(SP::DynamicalSystem ds)
{
  assert(ds && "MoreauJeanOSI::W(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
}


const SimpleMatrix MoreauJeanOSI::getWBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds &&
         "MoreauJeanOSI::getWBoundaryConditions(ds): ds == NULL.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions &&
         "MoreauJeanOSI::getWBoundaryConditions(ds): WBoundaryConditions[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions); // Copy !!
}

SP::SiconosMatrix MoreauJeanOSI::WBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds && "MoreauJeanOSI::WBoundaryConditions(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W ;
}


void MoreauJeanOSI::initializeDynamicalSystem(Model&, double t, SP::DynamicalSystem ds)
{
  DEBUG_BEGIN("MoreauJeanOSI::initializeDynamicalSystem(Model&, double t, SP::DynamicalSystem ds)\n");
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);
  workVectors.resize(MoreauJeanOSI::WORK_LENGTH);

  // Check dynamical system type
  Type::Siconos dsType = Type::value(*ds);
  assert(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS || Type::NewtonEulerDS
     || dsType == Type::LagrangianLinearDiagonalDS);

  // Compute W (iteration matrix)
  initializeIterationMatrixW(t, ds);

  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearDiagonalDS)
  {
    // buffers allocation (inside the graph)
    SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS> (ds);
    workVectors[MoreauJeanOSI::RESIDU_FREE].reset(new SiconosVector(lds->dimension()));
    workVectors[MoreauJeanOSI::VFREE].reset(new SiconosVector(lds->dimension()));
    workVectors[MoreauJeanOSI::BUFFER].reset(new SiconosVector(lds->dimension()));

    // Update dynamical system components (for memory swap).
    lds->computeForces(t, lds->q(), lds->velocity());
    lds->swapInMemory();
  }
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS neds = std11::static_pointer_cast<NewtonEulerDS> (ds);
    DEBUG_PRINTF("neds->number() %i \n",neds->number());
    // buffers allocation (into the graph)
    workVectors[MoreauJeanOSI::RESIDU_FREE].reset(new SiconosVector(neds->dimension()));
    workVectors[MoreauJeanOSI::VFREE].reset(new SiconosVector(neds->dimension()));

    //Compute a first value of the dotq  to store it in  _dotqMemory
    SP::SiconosMatrix T = neds->T();
    SP::SiconosVector dotq = neds->dotq();
    SP::SiconosVector v = neds->twist();
    prod(*T, *v, *dotq, true);

    //Compute a first value of the forces to store it in _forcesMemory
    neds->computeForces(t, neds->q(), v);
    neds->swapInMemory();
  }
  DEBUG_END("MoreauJeanOSI::initializeDynamicalSystem(Model&, double t, SP::DynamicalSystem ds)\n");

}
void MoreauJeanOSI::fillDSLinks(Interaction &inter, InteractionProperties& interProp,
                                  DynamicalSystemsGraph & DSG)
{
  DEBUG_BEGIN("MoreauJeanOSI::fillDSLinks(Interaction &inter, InteractionProperties& interProp, DynamicalSystemsGraph & DSG)\n");
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);
  assert(interProp.workVectors);
  
  VectorOfVectors& workV = *interProp.workVectors;
  workV.resize(MoreauJeanOSI::WORK_INTERACTION_LENGTH);
  workV[MoreauJeanOSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));

  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  Relation &relation =  *inter.relation();

  RELATION::TYPES relationType = relation.getType();

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

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
    DSlink[xfree]->setVectorPtr(0,workVds1[MoreauJeanOSI::VFREE]);
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
      DSlink[xfree]->setVectorPtr(1,workVds2[MoreauJeanOSI::VFREE]);
    }
  }
  DEBUG_END("MoreauJeanOSI::fillDSLinks(Interaction &inter, InteractionProperties& interProp, DynamicalSystemsGraph & DSG)\n");

}

void MoreauJeanOSI::initialize_nonsmooth_problems()
{
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->setIndexSetLevel(1);
  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->setInputOutputLevel(1);
  //  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->initialize(_simulation);
}

void MoreauJeanOSI::initializeIterationMatrixW(double time, SP::DynamicalSystem ds)
{
  DEBUG_BEGIN("MoreauJeanOSI::initializeIterationMatrixW\n");
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending on its type.

  if(!ds)
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixW(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);

  if(_dynamicalSystemsGraph->properties(dsv).W)
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixW(t,ds) - W(ds) is already in the map and has been initialized.");

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);
  unsigned int sizeW = ds->dimension();
  if(dsType == Type::LagrangianDS)
  {
    LagrangianDS& d = static_cast<LagrangianDS&> (*ds);
    // Memory allocation for W property of the grap
    if(d.mass())
    {
      d.computeMass(d.q());
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d.mass())); //*W = *d->mass();
    }
    else
    {
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(sizeW, sizeW));
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }
    // Compute the W matrix
    computeW(time,d, *_dynamicalSystemsGraph->properties(dsv).W);
    // WBoundaryConditions initialization
    if(d.boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(*ds, dsv);
  }
  // 2 - Lagrangian linear systems
  else if(dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);
    if(d->mass())
    {
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->mass())); //*W = *d->mass();
    }
    else
    {
      _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(sizeW, sizeW));
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }

    SP::SiconosMatrix K = d->K();
    SP::SiconosMatrix C = d->C();
    SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(dsv).W;
    if(C)
      scal(h * _theta, *C, *W, false); // W += h*_theta *C
    if(K)
      scal(h * h * _theta * _theta, *K, *W, false); // W = h*h*_theta*_theta*K

    // WBoundaryConditions initialization
    if(d->boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(*d,dsv);
  }
  else if(dsType == Type::LagrangianLinearDiagonalDS)
  {
    LagrangianLinearDiagonalDS& lldds = static_cast<LagrangianLinearDiagonalDS&> (*ds);
    unsigned int ndof = lldds.dimension();
    _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(ndof, ndof, Siconos::BANDED, 0, 0));
    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(dsv).W;

    if(lldds.mass())
      W = *lldds.mass();
    else
      W.eye();

    double htheta = h * _theta;
    double h2theta2 = h * h * _theta * _theta;
    if(lldds.damping())
    {
      SiconosVector& C = *lldds.damping();
      for(unsigned int i=0;i<ndof;++i)
      {
        W(i, i) += htheta * C(i);
      }
    }

    if(lldds.stiffness())
    {
      SiconosVector& K = *lldds.stiffness();
      for(unsigned int i=0;i<ndof;++i)
      {
        W(i, i) += h2theta2 * K(i);
      }
    }
    // WBoundaryConditions initialization
    if(lldds.boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(lldds, dsv);

    for(unsigned int i=0;i<ndof;++i)
    {
      W(i, i) = 1. / W(i, i);
    }
  }

  // === ===
  else if(dsType == Type::NewtonEulerDS)
  {
    NewtonEulerDS& d = static_cast<NewtonEulerDS&> (*ds);
    _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d.mass()));

    computeW(time, d, *_dynamicalSystemsGraph->properties(dsv).W);

    // WBoundaryConditions initialization
    if(d.boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(*ds,dsv);

  }
  else RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixW - not yet implemented for Dynamical system of type : " + Type::name(*ds));

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
  DEBUG_END("MoreauJeanOSI::initializeIterationMatrixW\n");
}


void MoreauJeanOSI::_initializeIterationMatrixWBoundaryConditions(DynamicalSystem& ds, const DynamicalSystemsGraph::VDescriptor& dsv)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  DEBUG_BEGIN("MoreauJeanOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds)\n");

  if(!(checkOSI(dsv)))
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - ds does not belong to the OSI.");

  if(_dynamicalSystemsGraph->properties(dsv).WBoundaryConditions)
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - WBoundaryConditions(ds) is already in the map and has been initialized.");


  Type::Siconos dsType = Type::value(ds);
  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS || dsType == Type::NewtonEulerDS
     || dsType == Type::LagrangianLinearDiagonalDS)
  {
    // Memory allocation for WBoundaryConditions
    unsigned int sizeWBoundaryConditions = ds.dimension(); // n for first order systems, ndof for lagrangian.

    SP::BoundaryCondition bc;
    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);
      bc = d.boundaryConditions();
    }
    else if(dsType == Type::NewtonEulerDS)
    {
     NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
      bc = d.boundaryConditions();
    }
    unsigned int numberBoundaryConditions = bc->velocityIndices()->size();
    _dynamicalSystemsGraph->properties(dsv).WBoundaryConditions.reset(new SimpleMatrix(sizeWBoundaryConditions, numberBoundaryConditions));
    _computeWBoundaryConditions(ds,*_dynamicalSystemsGraph->properties(dsv).WBoundaryConditions,
                                *_dynamicalSystemsGraph->properties(dsv).W);
  }
  else
    RuntimeException::selfThrow("MoreauJeanOSI::initializeIterationMatrixWBoundaryConditions - not yet implemented for Dynamical system of type :" +  Type::name(ds));
  DEBUG_END("MoreauJeanOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds) \n");
}


void MoreauJeanOSI::_computeWBoundaryConditions(DynamicalSystem& ds, SiconosMatrix& WBoundaryConditions, SiconosMatrix& iteration_matrix)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  Type::Siconos dsType = Type::value(ds);
  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS ||  dsType == Type::NewtonEulerDS || dsType == Type::LagrangianLinearDiagonalDS)
  {

    SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

    int columnindex = 0;

    std::vector<unsigned int>::iterator itindex;

    SP::BoundaryCondition bc;

    //SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;

    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);
      bc = d.boundaryConditions();
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
      bc = d.boundaryConditions();
    }

    for(itindex = bc->velocityIndices()->begin() ;
        itindex != bc->velocityIndices()->end();
        ++itindex)
    {

      iteration_matrix.getCol(*itindex, *columntmp);
      /*\warning we assume that W is symmetric in the Lagrangian case
        we store only the column and not the row */

      WBoundaryConditions.setCol(columnindex, *columntmp);
      double diag = (*columntmp)(*itindex);
      columntmp->zero();
      (*columntmp)(*itindex) = diag;

      iteration_matrix.setCol(*itindex, *columntmp);
      iteration_matrix.setRow(*itindex, *columntmp);


      columnindex ++;
    }
  }
  else
    RuntimeException::selfThrow("MoreauJeanOSI::computeWBoundaryConditions - not yet implemented for Dynamical system type : " +  Type::name(ds));
}


void MoreauJeanOSI::computeW(double t, DynamicalSystem& ds, SiconosMatrix& W)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.
  DEBUG_PRINT("MoreauJeanOSI::computeW starts\n");

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(ds);

  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
  {
    // Nothing: W does not depend on time.
  }
  else if(dsType == Type::LagrangianDS)
  {


    LagrangianDS& d = static_cast<LagrangianDS&> (ds);
    if(d.mass())
    {
      d.computeMass();
      W = *d.mass();
    }
    else
      W.eye();

    if(d.jacobianqDotForces())
    {
      SiconosMatrix& C = *d.jacobianqDotForces(); // jacobian according to velocity
      d.computeJacobianqDotForces(t);
      scal(-h * _theta, C, W, false); // W -= h*_theta*C
    }

    if(d.jacobianqForces())
    {
      SiconosMatrix& K = *d.jacobianqForces(); // jacobian according to q
      d.computeJacobianqForces(t);
      scal(-h * h * _theta * _theta, K, W, false); //*W -= h*h*_theta*_theta**K;
    }
  }
  // === ===
  else if(dsType == Type::NewtonEulerDS)
  {
    NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
    W = *(d.mass());

    if(d.jacobianvForces())
    {
      SiconosMatrix& C = *d.jacobianvForces(); // jacobian according to velocity

      d.computeJacobianvForces(t);
      scal(-h * _theta, C, W, false); // W -= h*_theta*C
    }
    if(d.jacobianqForces())
    {
      SiconosMatrix& K = *d.jacobianqForces(); // jacobian according to q
      d.computeJacobianqForces(t);
      SiconosMatrix& T = *d.T();
      DEBUG_EXPR(T.display(););
      DEBUG_EXPR(K.display(););
      SP::SimpleMatrix  buffer(new SimpleMatrix(*(d.mass())));
      prod(K, T, *buffer, true);
      scal(-h * h * _theta * _theta, *buffer, W, false);
      //*W -= h*h*_theta*_theta**K;
    }
    DEBUG_EXPR(W.display(););
    DEBUG_EXPR_WE(std::cout <<  std::boolalpha << "W.isPLUFactorized() = "<< W.isPLUFactorized() << std::endl;);

  }
  else RuntimeException::selfThrow("MoreauJeanOSI::computeW - not yet implemented for Dynamical system of type : " +Type::name(ds));
  DEBUG_PRINT("MoreauJeanOSI::computeW ends\n");
  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}

void MoreauJeanOSI::computeInitialNewtonState()
{
  DEBUG_BEGIN("MoreauJeanOSI::computeInitialNewtonState()\n");
  // Compute the position value giving the initial velocity.
  // The goal is to save one newton iteration for nearly linear system
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem&  ds = *_dynamicalSystemsGraph->bundle(*dsi);

    if(_explicitNewtonEulerDSOperators)
    {
      if(Type::value(ds) == Type::NewtonEulerDS)
      {
        // The goal is to update T() one time at the beginning of the Newton Loop
        // We want to be explicit on this function since we do not compute their Jacobians.
        NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
        SP::SiconosVector qold = d.qMemory()->getSiconosVector(0);
        //SP::SiconosVector q = d->q();
        computeT(qold,d.T());
      }
    }
    // The goal is to converge in one iteration if the system is almost linear
    // we start the Newton loop q = q0+hv0
    updatePosition(ds);
  }
  DEBUG_END("MoreauJeanOSI::computeInitialNewtonState()\n");
}



double MoreauJeanOSI::computeResidu()
{
  DEBUG_BEGIN("MoreauJeanOSI::computeResidu()\n");
  // This function is used to compute the residu for each "MoreauJeanOSI-discretized" dynamical system.
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
  //SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem& ds = *_dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(ds); // Its type

    // 3 - Lagrangian Non Linear Systems
    if(dsType == Type::LagrangianDS)
    {
      DEBUG_PRINT("MoreauJeanOSI::computeResidu(), dsType == Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1) - h*(1-theta)*forces(ti,vi,qi) - p_i+1
      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& free = *workVectors[MoreauJeanOSI::VFREE];

      // -- Convert the DS into a Lagrangian one.
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SiconosVector &vold = *d.velocityMemory()->getSiconosVector(0);

      SiconosVector &v = *d.velocity(); // v = v_k,i+1
      //residuFree.zero();
      DEBUG_EXPR(residuFree.display());
      DEBUG_EXPR(vold.display());
      DEBUG_EXPR(v.display());

      residuFree = v;
      sub(residuFree, vold, residuFree);
      if(d.mass())
      {
        d.computeMass(d.q());
        prod(*(d.mass()), residuFree, residuFree); // residuFree = M(v - vold)
      }

      if(d.forces())
      {
        // Cheaper version: get forces(ti,vi,qi) from memory
        SiconosVector& fold = *d.forcesMemory()->getSiconosVector(0);
        double coef = -h * (1 - _theta);
        scal(coef, fold, residuFree, false);

        // Expensive computes forces(ti,vi,qi)
        // d.computeForces(told, qold, vold);
        // double coef = -h * (1 - _theta);
        // // residuFree += coef * fL_i
        // scal(coef, *d.forces(), residuFree, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        d.computeForces(t,d.q(),d.velocity());
        coef = -h * _theta;
        scal(coef, *d.forces(), residuFree, false);

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        //SP::SiconosVector qbasedonv(new SiconosVector(*qold));
        //*qbasedonv +=  h * ((1 - _theta)* *vold + _theta * *v);
        //d.computeForces(t, qbasedonv, v);
        //coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        //scal(coef, *d.forces(), *residuFree, false);


      }

      if(d.boundaryConditions())
      {
        d.boundaryConditions()->computePrescribedVelocity(t);

        unsigned int columnindex = 0;
        SimpleMatrix & WBoundaryConditions  = *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions ;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          double DeltaPrescribedVelocity =
            d.boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - v.getValue(*itindex);

          WBoundaryConditions.getCol(columnindex, *columntmp);
          residuFree -= *columntmp * (DeltaPrescribedVelocity);

          residuFree.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;
        }
      }

      free = residuFree; // copy residuFree into Workfree
      DEBUG_EXPR(residuFree.display());

      if(d.p(1))
        free -= *d.p(1); // Compute Residu in Workfree Notation !!
      // We use DynamicalSystem::free as tmp buffer

      if(d.boundaryConditions())
      {
        unsigned int columnindex = 0;
        SimpleMatrix& WBoundaryConditions = *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions ;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          double DeltaPrescribedVelocity =
            d.boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - v.getValue(*itindex);

          WBoundaryConditions.getCol(columnindex, *columntmp);

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
      DEBUG_PRINT("MoreauJeanOSI::computeResidu(), dsType == Type::LagrangianLinearTIDS\n");
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual for v = v_i
      // otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      LagrangianLinearTIDS& d = static_cast<LagrangianLinearTIDS&> (ds);

      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& free = *workVectors[MoreauJeanOSI::VFREE];


      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SiconosVector& qold = *d.qMemory()->getSiconosVector(0); // qi
      SiconosVector& vold = *d.velocityMemory()->getSiconosVector(0); //vi

      DEBUG_EXPR(qold.display(););
      DEBUG_EXPR(vold.display(););
      DEBUG_EXPR(d.q()->display(););
      DEBUG_EXPR(d.velocity()->display(););

      // --- ResiduFree computation Equation (1) ---
      residuFree.zero();
      double coeff;
      // -- No need to update W --

      if(d.C())
      {
        prod(h, *d.C() , vold, residuFree, false); // vfree += h*C*vi
      }
      if(d.K())
      {
        coeff = h * h * _theta;
        prod(coeff, *d.K(), vold, residuFree, false); // vfree += h^2*_theta*K*vi
        prod(h, *d.K(), qold, residuFree, false); // vfree += h*K*qi
      }

      if(d.fExt())
      {
        // computes Fext(ti)
        d.computeFExt(told);
        coeff = -h * (1 - _theta);
        scal(coeff, *(d.fExt()), residuFree, false); // vfree -= h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d.computeFExt(t);
        coeff = -h * _theta;
        scal(coeff, *(d.fExt()), residuFree, false); // vfree -= h*_theta * fext(ti+1)
      }


      // Computation of the complete residual Equation (2)
      //   ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      //       SP::SiconosMatrix M = d.mass();
      //       SP::SiconosVector realresiduFree (new SiconosVector(residuFree));
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
      //         d.computeFExt(told);
      //         coeff = -h*(1-_theta);
      //         scal(coeff, *Fext, *realresiduFree, false); // vfree -= h*(1-_theta) * fext(ti)
      //         // computes Fext(ti+1)
      //         d.computeFExt(t);
      //         coeff = -h*_theta;
      //         scal(coeff, *Fext, *realresiduFree, false); // vfree -= h*_theta * fext(ti+1)
      //       }


      if(d.boundaryConditions())
      {
        d.boundaryConditions()->computePrescribedVelocity(t);

        unsigned int columnindex = 0;
        SimpleMatrix& WBoundaryConditions = *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {

          double DeltaPrescribedVelocity =
            d.boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - vold.getValue(*itindex);

          WBoundaryConditions.getCol(columnindex, *columntmp);
          residuFree += *columntmp * (DeltaPrescribedVelocity);

          residuFree.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;

        }
      }

      free = residuFree; // copy residuFree into free
      if(d.p(1))
        free-= *d.p(1); // Compute Residu in Workfree Notation !!
      // We use free as tmp buffer
      DEBUG_EXPR(free.display());
      DEBUG_EXPR(residuFree.display());

      normResidu = 0.0; // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();

    }

    else if(dsType == Type::LagrangianLinearDiagonalDS)
    {
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual for v = v_i
      // otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      LagrangianLinearDiagonalDS& d = static_cast<LagrangianLinearDiagonalDS&> (ds);

      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& free = *workVectors[MoreauJeanOSI::VFREE];


      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SiconosVector& qold = *d.qMemory()->getSiconosVector(0); // qi
      SiconosVector& vold = *d.velocityMemory()->getSiconosVector(0); //vi
      // --- ResiduFree computation Equation (1) ---
      residuFree.zero();
      double coeff;
      // -- No need to update W --
      if(d.damping())
      {
        SiconosVector & sigma = *d.damping();
        for(unsigned int i=0;i<d.dimension();++i)
          residuFree(i) += h * sigma(i) * vold(i);
      }
      if(d.stiffness())
      {
        coeff = h * h * _theta;
        SiconosVector & omega = *d.stiffness();
        for(unsigned int i=0;i<d.dimension();++i)
          residuFree(i) += coeff * omega(i) * vold(i) + h * omega(i) * qold(i);
      }

      if(d.fExt())
      {
        // computes Fext(ti)
        d.computeFExt(told);
        coeff = -h * (1 - _theta);
        scal(coeff, *(d.fExt()), residuFree, false); // vfree -= h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d.computeFExt(t);
        coeff = -h * _theta;
        scal(coeff, *(d.fExt()), residuFree, false); // vfree -= h*_theta * fext(ti+1)
      }

      if(d.boundaryConditions())
      {
        d.boundaryConditions()->computePrescribedVelocity(t);

        unsigned int columnindex = 0;
        SimpleMatrix& WBoundaryConditions = *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {

          double DeltaPrescribedVelocity =
            d.boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - vold.getValue(*itindex);

          WBoundaryConditions.getCol(columnindex, *columntmp);
          residuFree += *columntmp * (DeltaPrescribedVelocity);

          residuFree.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;

        }
      }

      free = residuFree; // copy residuFree into free
      if(d.p(1))
        free-= *d.p(1); // Compute Residu in Workfree Notation !!
      // We use DynamicalSystem::free as tmp buffer

      normResidu = 0.0; // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();

    }


    else if(dsType == Type::NewtonEulerDS)
    {
      DEBUG_PRINT("MoreauJeanOSI::computeResidu(), dsType == Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) - h*(1-_theta)*forces(ti,vi,qi) - pi+1

      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& free = *workVectors[MoreauJeanOSI::VFREE];


      // -- Convert the DS into a Lagrangian one.
      NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);

      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"
      SiconosVector& vold = *d.twistMemory()->getSiconosVector(0);

      // Get the current state vector
      //SiconosVector& q = *d.q();
      SiconosVector& v = *d.twist(); // v = v_k,i+1

      // Get the (constant mass matrix)
      SiconosMatrix &massMatrix = *d.mass();
      prod(massMatrix, (v - vold), residuFree, true); // residuFree = M(v - vold)
      DEBUG_EXPR(residuFree.display(););

      if(d.forces())   // if fL exists
      {
        DEBUG_PRINTF("MoreauJeanOSI:: _theta = %e\n",_theta);
        DEBUG_PRINTF("MoreauJeanOSI:: h = %e\n",h);

        // Cheaper version: get forces(ti,vi,qi) from memory
        SiconosVector& fold = *d.forcesMemory()->getSiconosVector(0);
        DEBUG_PRINT("MoreauJeanOSI:: old forces :\n");
        DEBUG_EXPR(fold.display(););

        double coef = -h * (1 - _theta);
        scal(coef, fold, residuFree, false);

        //Expensive version to check ...
        //SP::SiconosVector qold = d.qMemory()->getSiconosVector(0);
        //SP::SiconosVector vold = d.twistMemory()->getSiconosVector(0);
        // d.computeForces(told,qold,vold);
        // DEBUG_EXPR(d.forces()->display(););
        //double coef = -h * (1.0 - _theta);
        //scal(coef, *d.forces(), *residuFree, false);

        DEBUG_EXPR(residuFree.display(););

        // computes forces(ti,v,q)
        d.computeForces(t,d.q(),d.twist());
        coef = -h * _theta;
        scal(coef, *d.forces(), residuFree, false);
        DEBUG_PRINT("MoreauJeanOSI:: new forces :\n");
        DEBUG_EXPR(d.forces()->display(););
        DEBUG_EXPR(residuFree.display(););

      }


      if(d.boundaryConditions())
      {
        d.boundaryConditions()->computePrescribedVelocity(t);

        unsigned int columnindex = 0;
        SimpleMatrix& WBoundaryConditions = *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {

          DEBUG_PRINTF("columnindex = %i\n",columnindex);
          DEBUG_PRINTF("*itindex = %i\n",*itindex);
          double DeltaPrescribedVelocity =
            d.boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - v.getValue(*itindex);

          DEBUG_EXPR(d.boundaryConditions()->prescribedVelocity()->display());

          WBoundaryConditions.getCol(columnindex, *columntmp);
          residuFree -= *columntmp * (DeltaPrescribedVelocity);


          residuFree.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;
        }
      }

      free = residuFree;

      if(d.p(1))
        free -= *d.p(1);// We use DynamicalSystem::free as tmp buffer


      if(d.boundaryConditions())
      {
        unsigned int columnindex = 0;
        SimpleMatrix &  WBoundaryConditions = *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          double DeltaPrescribedVelocity =
            d.boundaryConditions()->prescribedVelocity()->getValue(columnindex)
            - v.getValue(*itindex);

          WBoundaryConditions.getCol(columnindex, *columntmp);

          free.setValue(*itindex, - columntmp->getValue(*itindex)   * (DeltaPrescribedVelocity));

          columnindex ++;
        }
      }

      DEBUG_PRINT("MoreauJeanOSI::computeResidu :\n");
      DEBUG_EXPR(residuFree.display(););
      DEBUG_EXPR(if(d.p(1)) d.p(1)->display(););
      DEBUG_EXPR(free.display(););

      normResidu =free.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }
    else
      RuntimeException::selfThrow("MoreauJeanOSI::computeResidu - not yet implemented for Dynamical system of type: " + Type::name(ds));

    if(normResidu > maxResidu) maxResidu = normResidu;

  }
  DEBUG_END("MoreauJeanOSI::computeResidu()\n");
  return maxResidu;


}

void MoreauJeanOSI::computeFreeState()
{
  DEBUG_BEGIN("MoreauJeanOSI::computeFreeState()\n");
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  double t = _simulation->nextTime(); // End of the time step

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //


  //SP::DynamicalSystem ds; // Current Dynamical System.
  //SP::SiconosMatrix W; // W MoreauJeanOSI matrix of the current DS.
  Type::Siconos dsType ; // Type of the current DS.

  DynamicalSystemsGraph::VIterator dsi, dsend;

  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem & ds = *_dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(ds); // Its type
    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W; // Its W MoreauJeanOSI matrix of iteration.
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    // 3 - Lagrangian Non Linear Systems
    if(dsType == Type::LagrangianDS)
    {
      DEBUG_PRINT("MoreauJeanOSI::computeFreeState(), dsType == Type::LagrangianDS\n");
      // IN to be updated at current time: W, M, q, v, fL
      // IN at told: qi,vi, fLi

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      // vFree = v_k,i+1 - W^{-1} ResiduFree
      // with
      // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*forces(t,v_k,i+1, q_k,i+1) - h*(1-theta)*forces(ti,vi,qi)

      // -- Convert the DS into a Lagrangian one.
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SiconosVector &v = *d.velocity(); // v = v_k,i+1
      DEBUG_EXPR(v.display());

      // --- ResiduFree computation ---
      // ResFree = M(v-vold) - h*[theta*forces(t) + (1-theta)*forces(told)]
      //
      // vFree pointer is used to compute and save ResiduFree in this first step.
      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& vfree = *workVectors[MoreauJeanOSI::VFREE];

      vfree = residuFree;

      // -- Update W --
      // Note: during computeW, mass and jacobians of forces will be computed/
      computeW(t, d, W);

      // -- vfree =  v - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      W.PLUForwardBackwardInPlace(vfree);
      // -> compute real vfree
      vfree *= -1.0;
      vfree += v;
      DEBUG_EXPR(vfree.display());

    }
    // 4 - Lagrangian Linear Systems
    else if(dsType == Type::LagrangianLinearTIDS)
    {
      DEBUG_PRINT("MoreauJeanOSI::computeFreeState(), dsType == Type::LagrangianLinearTIDS\n");
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // vFree = v_i + W^{-1} ResiduFree    // with
      // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 + h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      LagrangianLinearTIDS& d = static_cast<LagrangianLinearTIDS&> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SiconosVector& vold = *d.velocityMemory()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---
      // vFree pointer is used to compute and save ResiduFree in this first step.

      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& vfree = *workVectors[MoreauJeanOSI::VFREE];

      vfree = residuFree;
      DEBUG_EXPR(vfree.display());
      W.PLUForwardBackwardInPlace(vfree);
      vfree *= -1.0;
      vfree += vold;

      DEBUG_EXPR(vfree.display());


    }
    // 4 - Lagrangian Linear Diagonal Systems
    else if(dsType == Type::LagrangianLinearDiagonalDS)
    {
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // vFree = v_i + W^{-1} ResiduFree    // with
      // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 + h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      LagrangianLinearDiagonalDS& d = static_cast<LagrangianLinearDiagonalDS&> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SiconosVector& vold = *d.velocityMemory()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---
      // vFree pointer is used to compute and save ResiduFree in this first step.

      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SiconosVector& vfree = *workVectors[MoreauJeanOSI::VFREE];
      // W is diagonal and contains the inverse of the iteration matrix!
      for(unsigned int i=0;i<d.dimension();++i)
        vfree(i) = -W(i, i) * vfree(i) + vold(i);

    }
    else if(dsType == Type::NewtonEulerDS)
    {
      // IN to be updated at current time: W, M, q, v, fL
      // IN at told: qi,vi,

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      // vFree = v_k,i+1 - W^{-1} ResiduFree
      // with
      // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*forces(t,v_k,i+1, q_k,i+1)
      //                                        - h*(1-theta)*forces(ti,vi,qi)

      // -- Convert the DS into a NewtonEuler one.
      NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);

      // --- ResiduFree computation ---
      // ResFree = M(v-vold) - h*[theta*forces(t) + (1-theta)*forces(told)]
      //
      // vFree pointer is used to compute and save ResiduFree in this first step.

      SiconosVector& residuFree = *workVectors[MoreauJeanOSI::RESIDU_FREE];
      SiconosVector& vfree = *workVectors[MoreauJeanOSI::VFREE];


      vfree = residuFree;

      // -- Update W --
      // Note: during computeW, mass and jacobians of forces will be computed/
      SimpleMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      computeW(t, d, W);
      SiconosVector& v = *d.twist(); // v = v_k,i+1

      // -- vfree =  v - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      //    std::cout<<"MoreauJeanOSI::computeFreeState residu free"<<endl;
      //    vfree->display();
      DEBUG_EXPR(residuFree.display(););

      W.PLUForwardBackwardInPlace(vfree);
      //    std::cout<<"MoreauJeanOSI::computeFreeState -WRfree"<<endl;
      //    vfree->display();
      //    scal(h,*vfree,*vfree);
      // -> compute real vfree
      vfree *= -1.0;
      DEBUG_EXPR(vfree.display(););
      vfree += v;
      DEBUG_EXPR(vfree.display(););
    }
    else
      RuntimeException::selfThrow("MoreauJeanOSI::computeFreeState - not yet implemented for Dynamical system of type: " +  Type::name(ds));

  }
  DEBUG_END("MoreauJeanOSI::computeFreeState()\n");
}

void MoreauJeanOSI::prepareNewtonIteration(double time)
{
  DEBUG_BEGIN(" MoreauJeanOSI::prepareNewtonIteration(double time)\n");
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem& ds = *_dynamicalSystemsGraph->bundle(*dsi);
    computeW(time, ds, *_dynamicalSystemsGraph->properties(*dsi).W);
  }

  if(!_explicitNewtonEulerDSOperators)
  {
    DynamicalSystemsGraph::VIterator dsi, dsend;

    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;

      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

      //  VA <2016-04-19 Tue> We compute T to be consitent with the Jacobian at the beginning of the Newton iteration and not at the end
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


struct MoreauJeanOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  Interaction& _inter;
  InteractionProperties& _interProp;
  
  _NSLEffectOnFreeOutput(OneStepNSProblem *p, Interaction& inter, InteractionProperties& interProp) :
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
    //VectorOfVectors &workV = *_interProp.workVectors;
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[MoreauJeanOSI::OSNSP_RHS];
    subscal(e, *_inter.y_k(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord, false);
  }

  void visit(const RelayNSL& nslaw)
  {
    // since velocity lower-/upper-bounds are fully specified in NSL,
    // nothing to do here
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[MoreauJeanOSI::OSNSP_RHS];

    // The normal part is multiplied depends on en
    if (nslaw.en() > 0.0)
    {
      osnsp_rhs (0) +=  nslaw.en() * (*_inter.y_k(_osnsp->inputOutputLevel()))(0);
    }
    // The tangential part is multiplied depends on et
    if (nslaw.et() > 0.0)
    {
      osnsp_rhs (1) +=  nslaw.et()  * (*_inter.y_k(_osnsp->inputOutputLevel()))(1);
      if (_inter.nonSmoothLaw()->size()>=2)
      {
        osnsp_rhs (2) +=  nslaw.et()  * (*_inter.y_k(_osnsp->inputOutputLevel()))(2);
      }
    }
  }
  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const ComplementarityConditionNSL& nslaw)
  {
    ;
  }
};


void MoreauJeanOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */
  DEBUG_BEGIN("MoreauJeanOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)\n");
  SP::OneStepNSProblems allOSNS  = _simulation->oneStepNSProblems();
  InteractionsGraph& indexSet = *osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  assert(indexSet.bundle(vertex_inter));
  Interaction& inter = *indexSet.bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet.properties(vertex_inter).DSlink;
  // Get relation and non smooth law types
  assert(inter.relation());
  RELATION::TYPES relationType = inter.relation()->getType();
  RELATION::SUBTYPES relationSubType = inter.relation()->getSubType();

  
  unsigned int sizeY = inter.nonSmoothLaw()->size();

  unsigned int relativePosition = 0;

  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  SP::SiconosMatrix  F;
  //  SP::BlockVector deltax;
  
  //SiconosVector& yForNSsolver = *inter.yForNSsolver()

  SiconosVector& osnsp_rhs = *(*indexSet.properties(vertex_inter).workVectors)[MoreauJeanOSI::OSNSP_RHS];
  
  SP::BlockVector Xfree;

  /** \todo VA. All of these values should be stored in a node in the interactionGraph
   * corresponding to the Interaction
   * when a MoreauJeanOSI scheme is used.
   */

//  deltax = DSlink[FirstOrderR::deltax];;

  if(relationType == NewtonEuler)
  {
    Xfree = DSlink[NewtonEulerR::xfree];
  }
  else if(relationType == Lagrangian)
  {
    Xfree = DSlink[LagrangianR::xfree];
  }

  assert(Xfree);
  DEBUG_EXPR(Xfree->display(););

  Interaction& mainInteraction = inter;
  assert(mainInteraction.relation());

  if(relationType == NewtonEuler)
  {
    if(std11::static_pointer_cast<NewtonEulerR>(mainInteraction.relation())->jachqT())
    {
      SiconosMatrix& CT =  *static_cast<NewtonEulerR&>(*mainInteraction.relation()).jachqT();
      coord[3] = CT.size(1);
      coord[5] = CT.size(1);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interaction for each ds of the
      // current Interaction.
      // XXX Big quirks !!! -- xhub
      subprod(CT, *Xfree, osnsp_rhs, coord, true);
    }

  }
  else
  {

    if(mainInteraction.relation()->C())
    {
      SiconosMatrix&  C = *mainInteraction.relation()->C() ;

      coord[3] = C.size(1);
      coord[5] = C.size(1);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interactionfor each ds of the
      // current Interaction.
      if(_useGammaForRelation)
      {
        RuntimeException::selfThrow("MoreauJeanOSI::computeFreeOutput Configuration not possible");
//        subprod(C, *deltax, osnsp_rhs, coord, true);
      }
      else
      {
        subprod(C, *Xfree, osnsp_rhs, coord, true);
      }
    }

    if(relationType == Lagrangian)
    {
      SP::SiconosMatrix ID(new SimpleMatrix(sizeY, sizeY));
      ID->eye();

      Index xcoord(8);
      xcoord[0] = 0;
      xcoord[1] = sizeY;
      xcoord[2] = 0;
      xcoord[3] = sizeY;
      xcoord[4] = 0;
      xcoord[5] = sizeY;
      xcoord[6] = 0;
      xcoord[7] = sizeY;

      // For the relation of type LagrangianRheonomousR
      if(relationSubType == RheonomousR)
      {
        if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
        {
          SiconosVector q = *DSlink[LagrangianR::q0];
          SiconosVector z = *DSlink[LagrangianR::z];

          std11::static_pointer_cast<LagrangianRheonomousR>(inter.relation())->computehDot(simulation()->getTkp1(), q, z);
          *DSlink[LagrangianR::z] = z;
          subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter.relation())->hDot()), osnsp_rhs, xcoord, false); // y += hDot
        }
        else
          RuntimeException::selfThrow("MoreauJeanOSI::computeFreeOutput not yet implemented for SICONOS_OSNSP ");
      }
      if(relationSubType == CompliantLinearTIR)
      {
        if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
        {

          SiconosMatrix&  C = *mainInteraction.relation()->C() ;
          double h = _simulation->timeStep();
          osnsp_rhs *= h * _theta ;

          /* we have to check that the value are at the beginnning of the time step */
          SiconosVector q = *DSlink[LagrangianR::q0];
          SiconosVector z = *DSlink[LagrangianR::z];
          SiconosVector v = *DSlink[LagrangianR::q1];


          // + C q_k
          subprod(C, q, osnsp_rhs, coord, false);
          // + h(1-_theta)v_k

          v *= (1-_theta)* h ;
          subprod(C, v, osnsp_rhs, coord, false);


          if (std11::static_pointer_cast<LagrangianCompliantLinearTIR>(inter.relation())->e())
          {
            SiconosVector& e = *std11::static_pointer_cast<LagrangianCompliantLinearTIR>(inter.relation())->e();
            osnsp_rhs += e;
          }
        }
        else
          RuntimeException::selfThrow("MoreauJeanOSI::computeFreeOutput not yet implemented for SICONOS_OSNSP ");
      }
      DEBUG_EXPR(osnsp_rhs.display(););


      // For the relation of type LagrangianScleronomousR
      if(relationSubType == ScleronomousR)
      {

      }


    }

  }

  if(inter.relation()->getType() == Lagrangian || inter.relation()->getType() == NewtonEuler)
  {
    _NSLEffectOnFreeOutput nslEffectOnFreeOutput = _NSLEffectOnFreeOutput(osnsp, inter,
                                                                          indexSet.properties(vertex_inter));
    inter.nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
  DEBUG_EXPR(osnsp_rhs.display(););


  DEBUG_END("MoreauJeanOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)\n");
}

void MoreauJeanOSI::integrate(double& tinit, double& tend, double& tout, int& notUsed)
{
  // Last parameter is not used (required for LsodarOSI but not for MoreauJeanOSI).

  double h = tend - tinit;
  tout = tend;


  SP::SiconosMatrix W;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    W =  _dynamicalSystemsGraph->properties(*dsi).W;
    Type::Siconos dsType = Type::value(*ds);

    if(dsType == Type::LagrangianLinearTIDS)
    {
      // get the ds
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);
      // get velocity pointers for current time step
      SP::SiconosVector v = d->velocity();
      // get q and velocity pointers for previous time step
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      // get p pointer

      SP::SiconosVector p = d->p(1);

      // velocity computation :
      //
      // v = vi + W^{-1}[ -h*C*vi - h*h*theta*K*vi - h*K*qi + h*theta*Fext(t) + h*(1-theta) * Fext(ti) ] + W^{-1}*pi+1
      //

      *v = *p; // v = p

      double coeff;
      // -- No need to update W --
      SP::SiconosMatrix C = d->C();
      if(C)
        prod(-h, *C, *vold, *v, false); // v += -h*C*vi

      SP::SiconosMatrix K = d->K();
      if(K)
      {
        coeff = -h * h * _theta;
        prod(coeff, *K, *vold, *v, false); // v += -h^2*theta*K*vi
        prod(-h, *K, *qold, *v, false); // v += -h*K*qi
      }

      SP::SiconosVector Fext = d->fExt();
      if(Fext)
      {
        // computes Fext(ti)
        d->computeFExt(tinit);
        coeff = h * (1 - _theta);
        scal(coeff, *Fext, *v, false); // v += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(tout);
        coeff = h * _theta;
        scal(coeff, *Fext, *v, false); // v += h*theta * fext(ti+1)
      }
      // -> Solve WX = v and set v = X
      W->PLUForwardBackwardInPlace(*v);
      *v += *vold;
    }
    else RuntimeException::selfThrow("MoreauJeanOSI::integrate - not yet implemented for Dynamical system of type :" +  Type::name(*ds));
  }
}

void MoreauJeanOSI::updatePosition(DynamicalSystem& ds)
{
  DEBUG_BEGIN("MoreauJeanOSI::updatePosition(SP::DynamicalSystem ds)\n");

  double h = _simulation->timeStep();

  Type::Siconos dsType = Type::value(ds);

  // 1 - Lagrangian Systems
  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
  {
    // get dynamical system
    LagrangianDS& d = static_cast<LagrangianDS&> (ds);

    // Compute q
    SiconosVector& v = *d.velocity();
    SiconosVector& q = *d.q();
    //  -> get previous time step state
    SiconosVector& vold = *d.velocityMemory()->getSiconosVector(0);
    SiconosVector& qold = *d.qMemory()->getSiconosVector(0);
    // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    double coeff = h * _theta;
    scal(coeff, v, q) ; // q = h*theta*v
    coeff = h * (1 - _theta);
    scal(coeff, vold, q, false); // q += h(1-theta)*vold
    q += qold;
  }
  else if(dsType == Type::NewtonEulerDS)
  {
    // get dynamical system
    NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
    SiconosVector &v = *d.twist();
    //DEBUG_EXPR(d.display());

    //compute q
    //first step consists in computing  \dot q.
    //second step consists in updating q.
    //
    SiconosMatrix& T = *d.T();
    SiconosVector& dotq = *d.dotq();
    DEBUG_EXPR(v.display());
    prod(T, v, dotq, true);
    DEBUG_EXPR(dotq.display());

    SiconosVector& q = *d.q();
    //  -> get previous time step state
    SiconosVector& dotqold = *d.dotqMemory()->getSiconosVector(0);
    DEBUG_EXPR(dotqold.display());


    // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    double coeff = h * _theta;
    scal(coeff, dotq, q) ; // q = h*theta*v

    coeff = h * (1 - _theta);
    scal(coeff, dotqold, q, false); // q += h(1-theta)*vold

    SiconosVector& qold = *d.qMemory()->getSiconosVector(0);
    DEBUG_EXPR(qold.display());
    q += qold;   // q += qold

    DEBUG_PRINT("new q before normalizing\n");
    DEBUG_EXPR(q.display());

    //q[3:6] must be normalized
    d.normalizeq();
    DEBUG_PRINT("new q after normalizing\n");
    DEBUG_EXPR(q.display());

  }
  DEBUG_END("MoreauJeanOSI::updatePosition(SP::DynamicalSystem ds)\n");

}

void MoreauJeanOSI::updateState(const unsigned int )
{

  DEBUG_BEGIN("MoreauJeanOSI::updateState(const unsigned int )\n");

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

    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
    // Get the DS type

    Type::Siconos dsType = Type::value(ds);

    // 3 - Lagrangian Systems
    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
    {
      DEBUG_PRINT("MoreauJeanOSI::updateState(const unsigned int ), dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS \n");
      // get dynamical system
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);
      SiconosVector& vfree = *workVectors[MoreauJeanOSI::VFREE];

      //    SiconosVector *vfree = d.velocityFree();
      SiconosVector &v = *d.velocity();
      bool baux = dsType == Type::LagrangianDS && useRCC && _simulation->relativeConvergenceCriterionHeld();

      if(d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0)
      {

        assert(((d.p(_levelMaxForInput)).get()) &&
               " MoreauJeanOSI::updateState() *d.p(_levelMaxForInput) == NULL.");
        v = *d.p(_levelMaxForInput); // v = p
        if(d.boundaryConditions())
          for(std::vector<unsigned int>::iterator
                itindex = d.boundaryConditions()->velocityIndices()->begin() ;
              itindex != d.boundaryConditions()->velocityIndices()->end();
              ++itindex)
            v.setValue(*itindex, 0.0);
        if(dsType == Type::LagrangianLinearDiagonalDS)
        {
          for(unsigned int i=0;i<d.dimension();++i)
            v(i) = vfree(i) + W(i, i) * v(i);
        }
        else
        {
          W.PLUForwardBackwardInPlace(v);
          v +=  vfree;
        }
      }
      else
      {
        v =  vfree;
      }
      DEBUG_EXPR(v.display());



      if(d.boundaryConditions())
      {
        int bc = 0;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          _dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions->getCol(bc, *columntmp);
          /*\warning we assume that W is symmetric in the Lagrangian case*/
          if (!_dynamicalSystemsGraph->properties(*dsi).W->isSymmetric(1e-10))
            std::cout <<"Warning, we apply boundary conditions assuming W symmetric" << std::endl;
          double value = - inner_prod(*columntmp, v);
          if( d.p(_levelMaxForInput)&& d.p(_levelMaxForInput)->size() > 0)
          {
            value += (d.p(_levelMaxForInput))->getValue(*itindex);
          }
          /* \warning the computation of reactionToBoundaryConditions take into
             account the contact impulse but not the external and internal forces.
             A complete computation of the residu should be better */
          d.reactionToBoundaryConditions()->setValue(bc, value) ;
          bc++;
        }
      }

      SiconosVector& q = *d.q();
      SiconosVector& local_buffer = *workVectors[MoreauJeanOSI::BUFFER];
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
      DEBUG_PRINT("MoreauJeanOSI::updateState(const unsigned int), dsType == Type::NewtonEulerDS \n");

      // get dynamical system
      NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
      SiconosVector& v = *d.twist();
      // DEBUG_PRINT("MoreauJeanOSI::updateState()\n ")
      // DEBUG_EXPR(d.display());
      DEBUG_PRINT("MoreauJeanOSI::updateState() prev v\n")
        DEBUG_EXPR(v.display());

      // failure on bullet sims
      // d.p(_levelMaxForInput) is checked in next condition
      // assert(((d.p(_levelMaxForInput)).get()) &&
      //       " MoreauJeanOSI::updateState() *d.p(_levelMaxForInput) == NULL.");

      SiconosVector& vfree = *workVectors[MoreauJeanOSI::VFREE];


      if( d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0)
      {
        /*d.p has been fill by the Relation->computeInput, it contains
          B \lambda _{k+1}*/
        v = *d.p(_levelMaxForInput); // v = p
        if(d.boundaryConditions())
          for(std::vector<unsigned int>::iterator
                itindex = d.boundaryConditions()->velocityIndices()->begin() ;
              itindex != d.boundaryConditions()->velocityIndices()->end();
              ++itindex)
            v.setValue(*itindex, 0.0);

        _dynamicalSystemsGraph->properties(*dsi).W->PLUForwardBackwardInPlace(v);

        DEBUG_EXPR(d.p(_levelMaxForInput)->display());
        DEBUG_PRINT("MoreauJeanOSI::updatestate W CT lambda\n");
        DEBUG_EXPR(v.display());
        v +=  vfree;
      }
      else
        v =  vfree;

      DEBUG_PRINT("MoreauJeanOSI::updatestate work free\n");
      DEBUG_EXPR(vfree.display());
      DEBUG_PRINT("MoreauJeanOSI::updatestate new v\n");
      DEBUG_EXPR(v.display());

      if(d.boundaryConditions())
      {
        int bc = 0;
        SP::SiconosVector columntmp(new SiconosVector(ds.dimension()));

        for(std::vector<unsigned int>::iterator  itindex = d.boundaryConditions()->velocityIndices()->begin() ;
            itindex != d.boundaryConditions()->velocityIndices()->end();
            ++itindex)
        {
          _dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions->getCol(bc, *columntmp);
          /*\warning we assume that W is symmetric in the Lagrangian case*/
          double value = - inner_prod(*columntmp, v);
          if( d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0)
          {
            value += (d.p(_levelMaxForInput))->getValue(*itindex);
          }
          /* \warning the computation of reactionToBoundaryConditions take into
             account the contact impulse but not the external and internal forces.
             A complete computation of the residu should be better */
          d.reactionToBoundaryConditions()->setValue(bc, value) ;
          bc++;
        }
      }

      updatePosition(ds);

    }
    else RuntimeException::selfThrow("MoreauJeanOSI::updateState - not yet implemented for Dynamical system of type: " +  Type::name(ds));

  }
  DEBUG_END("MoreauJeanOSI::updateState(const unsigned int)\n");
}


bool MoreauJeanOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
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
  DEBUG_PRINTF("MoreauJeanOSI::addInteractionInIndexSet of level = %i yref=%e, yDot=%e, y_estimated=%e.\n", i,  y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));
  DEBUG_EXPR(
    if(y <= 0)
      DEBUG_PRINT("MoreauJeanOSI::addInteractionInIndexSet ACTIVATE.\n");
    );
  return (y <= 0.0);
}


bool MoreauJeanOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
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
  DEBUG_PRINTF("MoreauJeanOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));

  DEBUG_EXPR(
    if(y > 0)
      DEBUG_PRINT("MoreauJeanOSI::removeInteractionInIndexSet DEACTIVATE.\n");
    );
  return (y > 0.0);
}



void MoreauJeanOSI::display()
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
