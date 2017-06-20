
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
#include "NewMarkAlphaOSI.hpp"
#include "Simulation.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianScleronomousR.hpp"
#include "LagrangianR.hpp"
#include "NonSmoothLaw.hpp"
#include "NewtonEulerR.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"
using namespace RELATION;
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
NewMarkAlphaOSI::NewMarkAlphaOSI(double new_beta, double new_gamma, double new_alpha_m, double new_alpha_f, bool flag = false):
  OneStepIntegrator(OSI::NEWMARKALPHAOSI)
{
  _steps=1;
  _beta = new_beta;
  _gamma = new_gamma;
  _alpha_m = new_alpha_m;
  _alpha_f = new_alpha_f;
  _orderDenseOutput = 5.0;
  _IsVelocityLevel = flag;
  // Set levels. This may depend on the nonsmooth law and will be updated during fillDSLinks(...) call.
  _levelMinForOutput=0;
  _levelMaxForOutput=2;
  _levelMinForInput=1;
  _levelMaxForInput=2;
}

NewMarkAlphaOSI::NewMarkAlphaOSI(double _rho_infty, bool flag = false):
  OneStepIntegrator(OSI::NEWMARKALPHAOSI)
{
  _steps=1;
  _alpha_m = (2 * _rho_infty - 1) / (_rho_infty + 1);
  _alpha_f = _rho_infty / (_rho_infty + 1);
  _gamma = 0.5 + _alpha_f - _alpha_m;
  _beta = 0.25 * std::pow((_gamma + 0.5), 2);
  _orderDenseOutput = 5.0;
  _IsVelocityLevel = flag;
  // Set levels. This may depend on the nonsmooth law and will be updated during fillDSLinks(...) call.
  _levelMinForOutput=0;
  _levelMaxForOutput=2;
  _levelMinForInput=1;
  _levelMaxForInput=2;
}

const SimpleMatrix NewMarkAlphaOSI::getW(SP::DynamicalSystem ds)
{
  assert(ds && "NewMarkAlphaOSI::getW(ds): ds == NULL.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W && "NewMarkAlphaOSI::getW(ds): W[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W); // Copy !!
}

SP::SimpleMatrix NewMarkAlphaOSI::W(SP::DynamicalSystem ds)
{
  assert(ds && "NewMarkAlphaOSI::W(ds): ds == NULL.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W && "NewMarkAlphaOSI::W(ds): W[ds] == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
}

void NewMarkAlphaOSI::initializeIterationMatrixW(SP::DynamicalSystem ds)
{

  if(!ds)
    RuntimeException::selfThrow("NewMarkAlphaOSI::initializeIterationMatrixW(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("NewMarkAlphaOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  if(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W)
    RuntimeException::selfThrow("NewMarkAlphaOSI::initializeIterationMatrixW(t,ds) - W(ds) is already in the map and has been initialized.");

  _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W.reset(
    new SimpleMatrix(ds->dimension(), ds->dimension())); // allocate memory
  computeW(ds,*_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W);
}


void NewMarkAlphaOSI::computeW(SP::DynamicalSystem ds, SiconosMatrix& W)
{
  double beta_prime = (1 - _alpha_m) / ((1 - _alpha_f) * _beta);
  double gamma_prime = _gamma / _beta;
  double h = _simulation->nextTime() - _simulation->startingTime(); // step size
  if(h < 100 * MACHINE_PREC)
    RuntimeException::selfThrow("In NewMarkAlphaOSI::initializeIterationMatrixW(t,ds), time integration is too small");
  // make sure that W is initialized before computing
  Type::Siconos dsType = Type::value(*ds);
  SP::SiconosMatrix M;
  SP::SiconosMatrix K;
  SP::SiconosMatrix C;
  if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
  {
    SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS>(ds);
    if(dsType == Type::LagrangianDS)
    {
      K = lds->jacobianqForces(); // jacobian according to q
      C = lds->jacobianqDotForces(); // jacobian according to velocity
      lds->computeMass(lds->q());
    }
    else // LagrangianLinearTIDS
    {
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS>(ds);
      K = d->K();       // matrix K
      if(K)
        *K *= -1.0;     // K = -K
      C = d->C();       // matrix C
      if(C)
        *C *= -1.0;     // C = -C
    }
    // Compute W = (beta_prime/h^2)*M - (gamma_prime/h)*C - K
    if(lds->mass())
      scal(beta_prime / (h * h), *(lds->mass()), W, true);
    else
    {
      W.eye();
      W *= beta_prime/( h * h);
    }
    if(C)
      scal(-gamma_prime / h, *C, W, false);
    if(K)
      scal(-1.0, *K, W, false);
    //
#ifdef DEBUG_NEWMARK
    std::cout.precision(15);
    std::cout << "Iteration matrix W: ";
    W->display();
#endif
  }
  else
  {
    RuntimeException::selfThrow("In NewMarkAlphaOSI::initializeIterationMatrixW(t,ds), this type of Dynamical System is not yet implemented");
  }
}


double NewMarkAlphaOSI::computeResidu()
{
  DEBUG_BEGIN("NewMarkAlphaOSI::computeResidu() \n");
  // Compute the residual for each Dynamical system at step n and at iteration k
  // R_{n,k} = M_{n,k} ddotq_{n,k} - F_{n,k} - p_{n,k}
  // R_free = M_{n,k} ddotq_{n,k} - F_{n,k};
  // Compute norm of R_{n,k} for each DS
  // Take maximum norm of R_{n,k} over all DS
  double t = _simulation->nextTime(); // End of the time step
  // Iteration through the set of Dynamical Systems.
  //
  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.
  double maxResidu = 0.0;
  double normResidu = 0.0;
  SP::SiconosVector _residu;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds); // Its type
    //    SP::SiconosVector freeR = ds->workspace(DynamicalSystem::freeresidu);
    //    freeR->zero();
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {

      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
      // get position, velocity and acceleration
      SP::SiconosVector q = d->q();
      SP::SiconosVector v = d->velocity();
      SP::SiconosVector a = d->acceleration();
      SP::SiconosVector F;
      // get the reaction force p
      SP::SiconosVector p = d->p(2);
      // Compute free residual
      residuFree.zero();
      if(d->mass())
        prod(*d->mass(), *a, residuFree, true); // freeR = M*a
      else
        residuFree = *a;
      // For LagrangianDS (non linear Lagrangian DS)
      if(dsType == Type::LagrangianDS)
	    {
	      // Update mass matrix
	      d->computeMass();
	      F = d->forces();
	      if(F)
          // Compute F = F_ext - F_int - F_Gyr
          d->computeForces(t, d->q(), d->velocity());
	    }
      // For LagrangianLinearTIDS
      if(dsType == Type::LagrangianLinearTIDS)
	    {
	      // We need to add F_int = Cv + Kq to freeR
	      SP::LagrangianLinearTIDS dtids = std11::static_pointer_cast<LagrangianLinearTIDS>(ds);
	      SP::SiconosMatrix K = dtids->K();
	      SP::SiconosMatrix C = dtids->C();
	      F = dtids->fExt(); // Note that for LagrangianLinearTIDS, F = F_ext
	      if(F)
        {
          dtids->computeFExt(t);
        }
	      if(K)
          prod(*K, *q, residuFree, false); // f = M*a + C*v + K*q
	      if(C)
          prod(*C, *v, residuFree, false); // freeR = M*a + C*v
	    }

      residuFree -= *F;            // freeR = _workspace[freeresidu] - F
      // Compute residual
      _residu.reset(new SiconosVector(residuFree)); // _residu = freeR
      *_residu -= *p;                                 // _residu = _workspace[freeresidu] - p
      // Compute Euclidean norm of the residual
      normResidu = _residu->norm2();
      // Take maximum value of norm over all DS
      if(normResidu > maxResidu)
	    {
	      maxResidu = normResidu;
	    }
      //

      DEBUG_EXPR(residuFree.display(););
      DEBUG_EXPR(_residu->display(););

    }
    else
    {
      RuntimeException::selfThrow("In NewMarkAlphaOSI::computeResidu(t,ds), this type of Dynamical System is not yet implemented");
    }
  }
  DEBUG_END("NewMarkAlphaOSI::computeResidu() \n");
  return maxResidu;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::computeFreeState()
{
  DEBUG_BEGIN("NewMarkAlphaOSI::computeFreeState()\n");
  // Compute delta q_free = - ((W_{n,k})^-1)*R_free
  //Loop through the set of DS
  SP::DynamicalSystem ds;   // Current Dynamical System.
  SP::SiconosMatrix W;      // W matrix of the current DS.
  Type::Siconos dsType ;    // Type of the current DS.

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds); // Its type
    // Get iteration matrix W, make sure that W was updated before
    W = _dynamicalSystemsGraph->properties(*dsi).W; // Its W matrix of iteration.


    SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
    SiconosVector& qfree = *workVectors[OneStepIntegrator::free];

    // -- Convert the DS into a Lagrangian one.
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
      qfree = residuFree;
      W->PLUForwardBackwardInPlace(qfree); //_qfree = (W^-1)*R_free
      qfree *= -1.0; //_qfree = -(W^-1)*R_free
      //
      DEBUG_EXPR(qfree.display(););
    }
    else
    {
      RuntimeException::selfThrow("In NewMarkAlphaOSI::computeResidu(t,ds), this type of Dynamical System is not yet implemented");
    }
  }
  DEBUG_END("NewMarkAlphaOSI::computeFreeState()\n");
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  DEBUG_BEGIN("NewMarkAlphaOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)\n");
  double t = _simulation->nextTime();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;

  // Get the type of relation
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();
  // Get the set of OSNSPs
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  // get the size of the interaction
  unsigned int sizeY = inter->nonSmoothLaw()->size();
  // get pointer to delta q_free of Dynamical Systems concerned with the interaction

  SP::BlockVector q_free;
  if(relationType == Lagrangian)
  {
    q_free = DSlink[LagrangianR::xfree];
  }
  assert(q_free);
  DEBUG_EXPR(q_free->display(););

  // get pointer to yForNSsolver vector
  

  SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[NewMarkAlphaOSI::OSNSP_RHS];

  assert(q_free && "In NewMarkAlphaOSI::computeFreeOutput: pointer q_free has not initialized yet");
  assert(inter->relation() && "In NewMarkAlphaOSI::computeFreeOutput, relation associated with the interaction does not exist.");
  SP::SiconosMatrix C = inter->relation()->C();
  assert(C && "In NewMarkAlphaOSI::computeFreeOutput: Jacobian matrix does not exist");
  if(relationType == Lagrangian)
  {
    if(relationSubType == RheonomousR)
    {
      RuntimeException::selfThrow("NewMarkAlphaOSI::computeFreeOutput  not yet implemented with LagrangianRheonomousR");
    }

    if(relationSubType == ScleronomousR)
    {
      Index coord(8);
      coord[0] = 0;
      coord[1] = sizeY;
      coord[2] = 0;
      coord[3] = C->size(1);
      coord[4] = 0;
      coord[5] = C->size(1);
      coord[6] = 0;
      coord[7] = sizeY;
      if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)  // LCP at acceleration level
	    {
	      subprod(*C, *q_free, osnsp_rhs, coord, true);
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
	      SP::LagrangianScleronomousR _SclerR = std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation());
	      _SclerR->computedotjacqhXqdot(t, *inter, DSlink);
	      subprod(*ID, *(_SclerR->dotjacqhXqdot()), osnsp_rhs, xcoord, false); // y += NonLinearPart
	    }
      else if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_POS]).get() == osnsp)  // LCP at position level
	    {
	      // Update Jacobian matrix
	      inter->relation()->computeJach(t, *inter, indexSet->properties(vertex_inter));
	      // compute osnsp_rhs = y_{n,k} + G*q_free
	      if(!_IsVelocityLevel)  // output at the position level y_{n,k} = g_{n,k}
        {
          inter->computeOutput(t, indexSet->properties(vertex_inter), 0); // Update output of level 0
          osnsp_rhs = *(inter->y(0)); //g_{n,k}
        }
	      else                  // output at the velocity level y_{n,k} = (h/gamma_prime)*dotg_{n,k}
        {
          double h = _simulation->nextTime() - _simulation->startingTime();
          double gamma_prime = _gamma / _beta;
          inter->computeOutput(t, indexSet->properties(vertex_inter), 1); // Update output of level 1
          osnsp_rhs = (h / gamma_prime) * (*(inter->y(1))); //(h/gamma_prime)*dotg_{n,k}
        }
	      subprod(*C, *q_free, osnsp_rhs, coord, false);
	    }
      else
	    {
	      RuntimeException::selfThrow("NewMarkAlphaOSI::computeFreeOutput, this OSNSP does not exist");
	    }
    }

    DEBUG_EXPR(osnsp_rhs.display(););
  }
  else
  {
    RuntimeException::selfThrow("In NewMarkAlphaOSI::computeFreeOutput, this type of relation is not yet implemented");
  }
  DEBUG_END("NewMarkAlphaOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)\n");
}

void NewMarkAlphaOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  DEBUG_BEGIN("NewMarkAlphaOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)\n")

  // Get work buffers from the graph
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);

  // Get work matrices from the graph
  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  VectorOfMatrices& workMatrices = *_dynamicalSystemsGraph->properties(dsv).workMatrices;
  // Initialize memory buffers
  // _dynamicalSystemsGraph->bundle(dsv)->initMemory(getSizeMem());
  // // Force dynamical system to its initial state
  // _dynamicalSystemsGraph->bundle(dsv)->resetToInitialState();
  // Check dynamical system type
  Type::Siconos dsType = Type::value(*ds);
  assert(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS);

  // Compute W (iteration matrix)
  initializeIterationMatrixW(ds);
  // allocate memory for work space for Newton iteration procedure
  assert(_dynamicalSystemsGraph->properties(dsv).W && "W is NULL");
  // ds->allocateWorkVector(DynamicalSystem::local_buffer,   _dynamicalSystemsGraph->properties(dsv).W->size(0));
  //Allocate the memory to stock the acceleration-like variable
  if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
  {

    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
    d->initRhs(t);

    workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
    workVectors[OneStepIntegrator::residu_free].reset(new SiconosVector(d->dimension()));
    workVectors[OneStepIntegrator::free].reset(new SiconosVector(d->dimension()));
    workVectors[OneStepIntegrator::acce_like].reset(new SiconosVector(*(d->acceleration()))); // set a0 = ddotq0
    workVectors[OneStepIntegrator::acce_memory].reset(new SiconosVector(*(d->acceleration()))); // set a0 = ddotq0

    // Allocate the memory to stock coefficients of the polynomial for the dense output
    workMatrices.resize(OneStepIntegrator::work_vector_of_matrix_size);
    workMatrices[OneStepIntegrator::dense_output_coefficients].reset(new SimpleMatrix(ds->dimension(), (getOrderDenseOutput() + 1)));

    //*(d->workspace(DynamicalSystem::acce_like)) = *(d->acceleration());
    // ds->allocateWorkVector(DynamicalSystem::acce_like, ds->dimension()); // allocate memory for the acceleration-like of DS
    // ds->allocateWorkVector(DynamicalSystem::acce_memory, ds->dimension()); // allocate memory to stock acceleration

    //          d->allocateWorkMatrix(OneStepIntegrator::dense_output_coefficients, ds->dimension(), (osi_NewMark->getOrderDenseOutput() + 1));
    DEBUG_EXPR(d->display());

    d->swapInMemory();
  }
  else
  {
    RuntimeException::selfThrow("In NewMarkAlphaOSI::initialize: this type of DS is not yet implemented");
  }

  DEBUG_END("NewMarkAlphaOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)\n")


    }
void NewMarkAlphaOSI::fillDSLinks(Interaction &inter,
                                    InteractionProperties& interProp,
                                    DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);
  
  VectorOfVectors& workV = *interProp.workVectors;
  workV.resize(NewMarkAlphaOSI::WORK_INTERACTION_LENGTH);
  workV[NewMarkAlphaOSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));

  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  Relation &relation =  *inter.relation();  
  RELATION::TYPES relationType = relation.getType();

  NonSmoothLaw & nslaw = *inter.nonSmoothLaw();
  Type::Siconos nslType = Type::value(nslaw);
 
  if (nslType == Type::NewtonImpactNSL || nslType == Type::MultipleImpactNSL)
  {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 2 ;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
  }
  else if (nslType ==  Type::NewtonImpactFrictionNSL)
  {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 4;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
    RuntimeException::selfThrow("NewMarkAlphaOSI::fillDSLinks  not yet implemented for nonsmooth law of type NewtonImpactFrictionNSL");
  }
  else
    RuntimeException::selfThrow("NewMarkAlphaOSI::fillDSLinks not yet implemented  for nonsmooth of type");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  /* allocate and set work vectors for the osi */
  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  {
    RuntimeException::selfThrow("NewMarkAlphaOSI::fillDSLinks. The implementation is not correct for two different OSI for one interaction");
  }
  VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == Lagrangian)
  {
    DSlink[LagrangianR::xfree].reset(new BlockVector());
    DSlink[LagrangianR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  }
  // else if (relationType == NewtonEuler)
  // {
  //   DSlink[NewtonEulerR::xfree].reset(new BlockVector());
  //   DSlink[NewtonEulerR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  // }

  if (ds1 != ds2)
  {
    VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == Lagrangian)
    {
      DSlink[LagrangianR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
    }
    // else if (relationType == NewtonEuler)
    // {
    //   DSlink[NewtonEulerR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
    // }
  }
}

void NewMarkAlphaOSI::prepareNewtonIteration(double time)
{
  // Compute matrix W for all Dynamical Systems
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
    computeW(ds, W);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::prediction()
{
  DEBUG_BEGIN("NewMarkAlphaOSI::prediction()\n");
  // Step size
  double h = _simulation->nextTime() - _simulation->startingTime();
  if(h < 100 * MACHINE_PREC)
    RuntimeException::selfThrow("In NewMarkAlphaOSI::prediction, time integration is too small");
  // Loop over all DS
  Type::Siconos dsType ;    // Type of the current DS.
  SP::SiconosVector _q, _dotq, _ddotq;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;


    dsType = Type::value(*ds); // Its type
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
      DEBUG_EXPR(d->display(););
      _q = d->q();                // generalized coordinate
      _dotq = d->velocity();      // generalized velocity
      _ddotq = d->acceleration(); // generalized acceleration
      SiconosVector& acce_like = *workVectors[OneStepIntegrator::acce_like];
      // Save the acceleration before the prediction
      acce_like = *(_ddotq);

      DEBUG_PRINT("Before prediction :\n")
        DEBUG_EXPR(_q->display(););
      DEBUG_EXPR(_dotq->display(););
      DEBUG_EXPR(_ddotq->display(););
      DEBUG_EXPR(acce_like.display(););

      *_q = *_q + (*_dotq) * h + acce_like * (h*h * (0.5 - _beta)); // q_{n+1} = q_n + (h^2)*(0.5 - beta)*a_n
      *_dotq = *_dotq + acce_like * (h * (1 - _gamma));              // dotq_{n+1} = dotq_n + h*(1 - gamma)*a_n
      acce_like = (_alpha_f / (1 - _alpha_m)) * (*_ddotq) - (_alpha_m / (1 - _alpha_m)) * acce_like; // a_{n+1} = (alpha_f*ddotq_n - alpha_m*a_n)/(1 - alpha_m)
      *_q = *_q + (h*h * _beta) * acce_like;
      *_dotq = *_dotq + (h * _gamma) * acce_like;
      _ddotq->zero();

      DEBUG_PRINT("After prediction :\n")
        DEBUG_EXPR(_q->display(););
      DEBUG_EXPR(_dotq->display(););
      DEBUG_EXPR(_ddotq->display(););
      DEBUG_EXPR(acce_like.display(););
    }
    else
    {
      RuntimeException::selfThrow("In NewMarkAlphaOSI::prediction: this type of DS is not yet implemented");
    }
  }
  DEBUG_END("NewMarkAlphaOSI::prediction()\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::correction()
{
  double h = _simulation->nextTime() - _simulation->startingTime();
  double beta_prime = (1 - _alpha_m) / ((1 - _alpha_f) * _beta);
  double gamma_prime = _gamma / _beta;
  //Make sure that the input of the concerned Dynamical Systems is updated after solving LCP
  Type::Siconos dsType ;    // Type of the current DS
  SP::SiconosVector delta_q;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    SP::SimpleMatrix W = _dynamicalSystemsGraph->properties(*dsi).W; // Its W matrix of iteration.
    ; // Iteration matrix W_{n+1,k} computed at kth iteration
    SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
    SiconosVector& acce_like = *workVectors[OneStepIntegrator::acce_like];

    
    dsType = Type::value(*ds); // Its type
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
      SP::SiconosVector _p = d->p(2); // resultant force p_{n+1,k+1} of DS at (k+1)th iteration
      // Compute delta_q = W_{n+1,k}^{-1}(p_{n+1,k+1} - r_{n+1,k})
      delta_q.reset(new SiconosVector(*_p - residuFree)); // copy (p_{n+1,k+1} - r_{n+1,k}) to delta_q
      W->PLUForwardBackwardInPlace(*delta_q);
      // Correction q_{n+1,k+1}, dotq_{n+1,k+1}, ddotq_{n+1,k+1}
      *(d->q()) += *delta_q; // q_{n+1,k+1} = q_{n+1,k} + delta_q
      *(d->velocity()) += (gamma_prime / h) * (*delta_q); // dotq_{n+1,k+1} = dotq_{n+1,k} + (gamma_prime/h)*delta_q
      *(d->acceleration()) += beta_prime / (h*h) * (*delta_q); // ddotq_{n+1,k+1} = ddotq_{n+1,k} + (beta_prime/h^2)*delta_q
      //a_{n+1,k+1} = a_{n+1,k} + ((1-alpha_f)/(1-alpha_m))*(beta_prime/h^2)*delta_q
      
      acce_like += ((1 - _alpha_f) / (1 - _alpha_m)) * ((beta_prime / (h*h)) * (*delta_q));

      DEBUG_PRINT("After correction : \n");
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););
      DEBUG_EXPR(d->acceleration()->display(););
      DEBUG_EXPR(acce_like.display(););

    }
    else
    {
      RuntimeException::selfThrow("In NewMarkAlphaOSI::updateState: this type of DS is not yet implemented");
    }
  }

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::integrate(double& t_ini, double& t_end, double& t_out, int& flag)
{
  RuntimeException::selfThrow("In NewMarkAlphaOSI::integrate, this method does nothing in the NewMarkAlpha OneStepIntegrator!!!");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::updateState(const unsigned int level)
{
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  if(level == 1)  // ie impact case: compute velocity
  {
    DynamicalSystemsGraph::VIterator dsi, dsend;
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS>(ds);
      lds->computePostImpactVelocity();
    }
  }
  else if(level == 2)
  {
    double time = _simulation->nextTime();
    DynamicalSystemsGraph::VIterator dsi, dsend;
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      ds->update(time);
    }
  }
  else RuntimeException::selfThrow("In NewMarkAlphaOSI::updateState, index is out of range. Index = " + level);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::computeCoefsDenseOutput(SP::DynamicalSystem ds)
{
  double h = _simulation->nextTime() - _simulation->startingTime();
  Type::Siconos dsType = Type::value(*ds);    // Type of the current DS
  SP::SiconosVector q_n, dotq_n, ddotq_n, q_np1, dotq_np1, ddotq_np1;
  SP::SiconosVector _vec(new SiconosVector(ds->dimension()));
  VectorOfMatrices& workMatrices = *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).workMatrices;
  VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).workVectors;
  if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
  {
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
    q_n = d->qMemory()->getSiconosVector(0); // q_n
    dotq_n = d->velocityMemory()->getSiconosVector(0); // dotq_n
    //ddotq_n = d->workspace(DynamicalSystem::acce_memory); // ddotq_n
    SiconosVector& ddotq_n = *workVectors[OneStepIntegrator::acce_memory];

    q_np1 = d->q(); // q_{n+1}
    dotq_np1 = d->velocity(); // dotq_{n+1}
    ddotq_np1 = d->acceleration(); // ddotq_{n+1}
    SP::SiconosMatrix _CoeffsDense = workMatrices[OneStepIntegrator::dense_output_coefficients];
    //d->workMatrix(OneStepIntegrator::dense_output_coefficients); // matrix of coefficients [a0 a1 a2 a3 a4 a5]
    if(_CoeffsDense->size(1) != 6)
    {
      RuntimeException::selfThrow("In NewMarkAlphaOSI::computeCoefsDenseOutput: the number of polynomial coeffcients considered here must equal to 6 (dense output polynomial of order 5)");
    }
    //a0 = q_n
    (*_vec) = (*q_n);
    _CoeffsDense->setCol(0, (*_vec));
    std::cout << "a0: ";
    _vec->display();
    //a1 = h*dotq_n
    (*_vec) = h * (*dotq_n);
    _CoeffsDense->setCol(1, (*_vec));
    std::cout << "a1: ";
    _vec->display();
    //a2 = 0.5*h^2*ddotq_n
    (*_vec) = (0.5 * h * h) * (ddotq_n);
    _CoeffsDense->setCol(2, (*_vec));
    std::cout << "a2: ";
    _vec->display();
    //a3 = -10*q_n - 6*h*dotq_n - 1.5*h^2*ddotq_n + 10*q_{n+1} - 4*h*dotq_{n+1} + 0.5*h^2*ddotq_{n+1}
    (*_vec) = (-10.0) * (*q_n) - (6.0 * h) * (*dotq_n) - (1.5 * h * h) * (ddotq_n) + 10.0 * (*q_np1) - (4.0 * h) * (*dotq_np1) + (0.5 * h *h) * (*ddotq_np1);
    _CoeffsDense->setCol(3, (*_vec));
    std::cout << "a3: ";
    _vec->display();
    //a4 = 15*q_n + 8*h*dotq_n + 1.5*h^2*ddotq_n - 15*q_{n+1} + 7*h*dotq_{n+1} - h^2*ddotq_{n+1}
    (*_vec) = 15.0 * (*q_n) + (8.0 * h) * (*dotq_n) + (1.5 * h *h) * (ddotq_n) - 15.0 * (*q_np1) + (7.0 * h) * (*dotq_np1) - h*h * (*ddotq_np1);
    _CoeffsDense->setCol(4, (*_vec));
    std::cout << "a4: ";
    _vec->display();
    //a5 = -6*q_n - 3*h*dotq_n - 0.5*h^2*ddotq_n + 6*q_{n+1} - 3*h*dotq_{n+1} + 0.5*h^2*ddotq_{n+1}
    (*_vec) = (-6.0) * (*q_n) - (3.0 * h) * (*dotq_n) - (0.5 * h*h) * (ddotq_n) + 6.0 * (*q_np1) - (3.0 * h) * (*dotq_np1) + (0.5 * h*h) * (*ddotq_np1);
    _CoeffsDense->setCol(5, (*_vec));
    std::cout << "a5: ";
    _vec->display();
    //
#ifdef DEBUG_NEWMARK
    std::cout << "==================== In NewMarkAlphaOSI::computeCoefsDenseOutput ================" <<std::endl;
    std::cout << "DS number: " << ds->number() <<std::endl;
    std::cout << "q_n: ";
    q_n->display();
    std::cout << "dotq_n: ";
    dotq_n->display();
    std::cout << "ddotq_n: ";
    ddotq_n->display();
    std::cout << "q_n+1: ";
    q_np1->display();
    std::cout << "dotq_n+1: ";
    dotq_np1->display();
    std::cout << "ddotq_n+1: ";
    ddotq_np1->display();
    std::cout << "Dense output coefficient matrix: " <<std::endl;
    _CoeffsDense->display();
#endif
  }
  else
  {
    RuntimeException::selfThrow("In NewMarkAlphaOSI::computeCoefsDenseOutput: this type of DS has not been implemented yet");
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::prepareEventLocalization()
{

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    // Compute coefficients of the dense output polynomial for all Dynamical Systems
    computeCoefsDenseOutput(ds);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::DenseOutputallDSs(double t)
{
  // Make sure that all coefficients of the dense output polynomial for all DSs has been computed before
  double t_n = _simulation->startingTime();
  double t_np1 = _simulation->nextTime();
  double h = t_np1 - t_n;
  double theta = (t - t_n) / h;
  SP::SiconosVector _vec1(new SiconosVector(_orderDenseOutput + 1));
  assert((_vec1->size() == 6) && "There are six coefficients of the dense output polynomial");
  (*_vec1)(0) = 1.0;
  (*_vec1)(1) = theta;
  (*_vec1)(2) = std::pow(theta, 2);
  (*_vec1)(3) = std::pow(theta, 3);
  (*_vec1)(4) = std::pow(theta, 4);
  (*_vec1)(5) = std::pow(theta, 5);
  //
  SP::SiconosVector _vec2(new SiconosVector(_orderDenseOutput + 1));
  assert((_vec2->size() == 6) && "There are six coefficients of the dense output polynomial");
  (*_vec2)(0) = 0.0;
  (*_vec2)(1) = 1.0 / h;
  (*_vec2)(2) = (2.0 * theta) / h;
  (*_vec2)(3) = (3.0 * std::pow(theta, 2)) / h;
  (*_vec2)(4) = (4.0 * std::pow(theta, 3)) / h;
  (*_vec2)(5) = (5.0 * std::pow(theta, 4)) / h;
  //
  SP::SiconosVector _vec3(new SiconosVector(_orderDenseOutput + 1));
  assert((_vec3->size() == 6) && "There are six coefficients of the dense output polynomial");
  (*_vec3)(0) = 0.0;
  (*_vec3)(1) = 0.0;
  (*_vec3)(2) = 2.0 / std::pow(h, 2);
  (*_vec3)(3) = (6.0 * theta) / std::pow(h, 2);
  (*_vec3)(4) = (12.0 * std::pow(theta, 2)) / std::pow(h, 2);
  (*_vec3)(5) = (20.0 * std::pow(theta, 3)) / std::pow(h, 2);
  //
  Type::Siconos dsType;    // Type of the current DS

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds);
    VectorOfMatrices& workMatrices = *_dynamicalSystemsGraph->properties(*dsi).workMatrices;

    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(ds);
      SiconosMatrix &coeffsDense = *workMatrices[OneStepIntegrator::dense_output_coefficients];
      prod(coeffsDense, *_vec1, *(d->q()), true); // q = Matrix_coeffs*_vec1
      prod(coeffsDense, *_vec2, *(d->velocity()), true); // dotq = Matrix_coeffs*_vec2
      prod(coeffsDense, *_vec3, *(d->acceleration()), true); // ddotq = Matrix_coeffs*_vec3
    }
    else
    {
      RuntimeException::selfThrow("In NewMarkAlphaOSI::DenseOutputallDSs: this type of DS has not been implemented yet");
    }

  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NewMarkAlphaOSI::display()
{

}
