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
#include "MoreauJeanGOSI.hpp"
#include "SiconosAlgebraProd.hpp"
#include "SiconosAlgebraScal.hpp"
//#include "SiconosVectorFriends.hpp"
#include "Simulation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "TypeName.hpp"

#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"


// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"


using namespace RELATION;

/// for non-owned shared pointers (passing const SiconosVector into
/// functions that take SP::SiconosVector without copy -- warning
/// const abuse!)
//static void null_deleter(const SiconosVector *) {}
// template <typename T> static std::shared_ptr<T> ptr(const T& a)
// {
//   return std::shared_ptr<SiconosVector>(&*(T*)&a, null_deleter);
// }


void MoreauJeanGOSI::initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  VectorOfVectors& ds_work_vectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type
  Type::Siconos dsType = Type::value(*ds);
  SP::SecondOrderDS sods = std::static_pointer_cast<SecondOrderDS> (ds);
  // Compute W (iteration matrix)
  initializeIterationMatrixW(t, sods);
  if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
  {
    SP::LagrangianDS lds = std::static_pointer_cast<LagrangianDS> (ds);

    ds_work_vectors.resize(MoreauJeanGOSI::WORK_LENGTH);
    ds_work_vectors[MoreauJeanGOSI::RESIDU_FREE].reset(new SiconosVector(lds->dimension()));
    ds_work_vectors[MoreauJeanGOSI::FREE].reset(new SiconosVector(lds->dimension()));
    ds_work_vectors[MoreauJeanGOSI::LOCAL_BUFFER].reset(new SiconosVector(lds->dimension()));

    lds->computeForces(t, lds->q(), lds->velocity());
    lds->swapInMemory();
  }
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS neds = std::static_pointer_cast<NewtonEulerDS> (ds);
    ds_work_vectors.resize(MoreauJeanGOSI::WORK_LENGTH);
    ds_work_vectors[MoreauJeanGOSI::RESIDU_FREE].reset(new SiconosVector(neds->dimension()));
    ds_work_vectors[MoreauJeanGOSI::FREE].reset(new SiconosVector(neds->dimension()));
    //Compute a first value of the dotq  to store it in  _dotqMemory
    SP::SiconosMatrix T = neds->T();
    SP::SiconosVector dotq = neds->dotq();
    SP::SiconosVector v = neds->twist();
    prod(*T, *v, *dotq, true);

    //Compute a first value of the forces to store it in _forcesMemory
    neds->computeForces(t, neds->q(), v);
    neds->swapInMemory();
  }
}

void MoreauJeanGOSI::initializeWorkVectorsForInteraction(Interaction &inter,
    InteractionProperties& interProp,
    DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);

  if(!interProp.workVectors)
  {
    interProp.workVectors.reset(new VectorOfVectors);
    interProp.workVectors->resize(MoreauJeanGOSI::WORK_INTERACTION_LENGTH);
  }

  if(!interProp.workBlockVectors)
  {
    interProp.workBlockVectors.reset(new VectorOfBlockVectors);
    interProp.workBlockVectors->resize(MoreauJeanGOSI::BLOCK_WORK_LENGTH);
  }

  VectorOfVectors& inter_work = *interProp.workVectors;
  VectorOfBlockVectors& inter_block_work = *interProp.workBlockVectors;


  if(!inter_work[MoreauJeanGOSI::OSNSP_RHS])
    inter_work[MoreauJeanGOSI::OSNSP_RHS].reset(new SiconosVector(inter.dimension()));

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);


  /* allocate and set work vectors for the osi */
  unsigned int xfree = MoreauJeanGOSI::xfree;

  if(ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");
    if((!inter_block_work[xfree]) || (inter_block_work[xfree]->numberOfBlocks() !=2))
      inter_block_work[xfree].reset(new BlockVector(2));
  }
  else
  {
    if((!inter_block_work[xfree]) || (inter_block_work[xfree]->numberOfBlocks() !=1))
      inter_block_work[xfree].reset(new BlockVector(1));
  }

  if(checkOSI(DSG.descriptor(ds1)))
  {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    inter_block_work[xfree]->setVectorPtr(0,workVds1[MoreauJeanGOSI::FREE]);
  }
  DEBUG_PRINTF("ds1->number() %i\n",ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n",ds2->number());


  if(ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");
    if(checkOSI(DSG.descriptor(ds2)))
    {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n",ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      inter_block_work[xfree]->setVectorPtr(1,workVds2[MoreauJeanGOSI::FREE]);
    }
  }
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
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds); // Its type

    // 3 - Lagrangian Non Linear Systems
    if(dsType == Type::LagrangianDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::computeResidu(), dsType == Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1) - h*(1-theta)*forces(ti,vi,qi) - p_i+1
      SiconosVector& residu = *ds_work_vectors[MoreauJeanGOSI::RESIDU_FREE];
      SiconosVector& free_rhs = *ds_work_vectors[MoreauJeanGOSI::FREE];

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = std::static_pointer_cast<LagrangianDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      // const SiconosVector& qold = d->qMemory().getSiconosVector(0);
      const SiconosVector& vold = d->velocityMemory().getSiconosVector(0);
      SP::SiconosVector q = d->q();

      SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      //residu.zero();
      DEBUG_EXPR(residu.display());

      // DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());
      DEBUG_EXPR(q->display());
      DEBUG_EXPR(v->display());
      free_rhs.zero();
      SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      prod(W, vold, free_rhs);



      if(d->forces())
      {
        // Cheaper version: get forces(ti,vi,qi) from memory
        const SiconosVector& fold = d->forcesMemory().getSiconosVector(0);
        double coef = h * (1 - _theta);
        scal(coef, fold, free_rhs, false);

        // Expensive computes forces(ti,vi,qi)
        // d->computeForces(told, qold, vold);
        // double coef = -h * (1 - _theta);
        // // residu += coef * fL_i
        // scal(coef, *d->forces(), *residu, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        d->computeForces(t,q,v);
        coef = h * _theta;
        scal(coef, *d->forces(), free_rhs, false);

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        //SP::SiconosVector qbasedonv(new SiconosVector(*qold));
        //*qbasedonv +=  h * ((1 - _theta)* *vold + _theta * *v);
        //d->computeForces(t, qbasedonv, v);
        //coef = -h * _theta;
        // residu += coef * fL_k,i+1
        //scal(coef, *d->forces(), residu, false);


      }


      residu =  -1.0* free_rhs;
      prod(1.0, W, *v, residu, false);

      DEBUG_EXPR(residu.display());

      if(d->p(1))
        residu -= *d->p(1); // Compute Residu in Workfree Notation !!

      if(d->boundaryConditions())
      {
        THROW_EXCEPTION("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }


      DEBUG_EXPR(residu.display());
      normResidu = residu.norm2();
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
      SP::LagrangianLinearTIDS d = std::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      const SiconosVector& qold = d->qMemory().getSiconosVector(0); // qi
      const SiconosVector& vold = d->velocityMemory().getSiconosVector(0); //vi

      DEBUG_EXPR(qold.display(););
      DEBUG_EXPR(vold.display(););
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););

      SiconosVector& residu = *ds_work_vectors[MoreauJeanGOSI::RESIDU_FREE];
      SiconosVector& free_rhs = *ds_work_vectors[MoreauJeanGOSI::FREE];
      // --- ResiduFree computation Equation (1) ---
      residu.zero();
      SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      prod(W, vold, free_rhs);

      double coeff;
      // -- No need to update W --

      SP::SiconosVector v = d->velocity(); // v = v_k,i+1

      SP::SiconosMatrix C = d->C();
      if(C)
        prod(h, *C, vold, free_rhs, false); // free_rhs += h*C*vi

      SP::SiconosMatrix K = d->K();
      if(K)
      {
        coeff = -h * h * _theta;
        prod(coeff, *K, vold, residu, false); // free_rhs += -h^2*_theta*K*vi
        prod(-h, *K, qold, free_rhs, false); // free_rhs += -h*K*qi
      }

      SP::SiconosVector Fext = d->fExt();
      if(Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - _theta);
        scal(coeff, *(d->fExt()), free_rhs, false); // free_rhs += h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * _theta;
        scal(coeff, *(d->fExt()), free_rhs, false); // free_rhs += h*_theta * fext(ti+1)
      }
      DEBUG_EXPR(free_rhs.display());

      if(d->boundaryConditions())
      {
        THROW_EXCEPTION("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }

      // residu = -1.0*free_rhs;
      // prod(1.0, W, *v, residu, false);
      // DEBUG_EXPR(free_rhs.display());
      // if(d->p(1))
      //   residu -= *d->p(1); // Compute Residu in Workfree Notation !!

      normResidu = 0.0; // we assume that v = vfree + W^(-1) p
    }



    else if(dsType == Type::NewtonEulerDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::computeResidu(), dsType == Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) - h*(1-_theta)*forces(ti,vi,qi) - pi+1


      // -- Convert the DS into a Lagrangian one.
      SP::NewtonEulerDS d = std::static_pointer_cast<NewtonEulerDS> (ds);
      DEBUG_EXPR(d->display());
      SiconosVector& residu = *ds_work_vectors[MoreauJeanGOSI::RESIDU_FREE];
      SiconosVector& free_rhs = *ds_work_vectors[MoreauJeanGOSI::FREE];
      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"
      const SiconosVector& vold = d->twistMemory().getSiconosVector(0);


      // Get the current state vector
      SP::SiconosVector q = d->q();
      SP::SiconosVector v = d->twist(); // v = v_k,i+1

      DEBUG_EXPR(vold.display());
      DEBUG_EXPR(q->display());
      DEBUG_EXPR(v->display());


      residu.zero();
      // Get the (constant mass matrix)
      // SP::SiconosMatrix massMatrix = d->mass();
      // prod(*massMatrix, (*v - vold), residu, true); // residu = M(v - vold)
      // DEBUG_EXPR(residu.display(););
      free_rhs.zero();
      SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      prod(W, vold, free_rhs);



      if(d->forces())   // if fL exists
      {
        DEBUG_PRINTF("MoreauJeanGOSI:: _theta = %e\n",_theta);
        DEBUG_PRINTF("MoreauJeanGOSI:: h = %e\n",h);

        // Cheaper version: get forces(ti,vi,qi) from memory
        const SiconosVector& fold = d->forcesMemory().getSiconosVector(0);
        double coef = h * (1 - _theta);
        scal(coef, fold, free_rhs, false);

        // Expensive version to check ...
        //d->computeForces(told,qold,vold);
        //double coef = -h * (1.0 - _theta);
        //scal(coef, *d->forces(), residu, false);

        DEBUG_PRINT("MoreauJeanGOSI:: old forces :\n");
        DEBUG_EXPR(d->forces()->display(););
        DEBUG_EXPR(residu.display(););

        // computes forces(ti,v,q)
        d->computeForces(t,q,v);
        coef = h * _theta;
        scal(coef, *d->forces(), free_rhs, false);
        DEBUG_PRINT("MoreauJeanGOSI:: new forces :\n");
        DEBUG_EXPR(d->forces()->display(););
        DEBUG_EXPR(residu.display(););

      }


      if(d->boundaryConditions())
      {
        THROW_EXCEPTION("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }

      residu =  -1.0* free_rhs;

      prod(1.0, W, *v, residu, false);


      if(d->p(1))
        residu -= *d->p(1);


      if(d->boundaryConditions())
      {
        THROW_EXCEPTION("MoreauJeanGOSI::computeResidu - boundary conditions not yet implemented for Dynamical system of type: " + Type::name(*ds));
      }

      DEBUG_PRINT("MoreauJeanGOSI::computeResidu :\n");
      DEBUG_EXPR(residu.display(););
      DEBUG_EXPR(if(d->p(1)) d->p(1)->display(););
      DEBUG_EXPR(free_rhs.display(););

      normResidu = residu.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }
    else
      THROW_EXCEPTION("MoreauJeanGOSI::computeResidu - not yet implemented for Dynamical system of type: " + Type::name(*ds));

    if(normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void MoreauJeanGOSI::computeFreeState()
{
  DEBUG_BEGIN("MoreauJeanGOSI::computeFreeState()\n");
  DEBUG_END("MoreauJeanGOSI::computeFreeState()\n");
}


void MoreauJeanGOSI::NonSmoothLawContributionToOutput(SP::Interaction inter, OneStepNSProblem& osnsp)
{
  if(inter->relation()->getType() == Lagrangian || inter->relation()->getType() == NewtonEuler)
  {
    InteractionsGraph& indexSet = *osnsp.simulation()->indexSet(osnsp.indexSetLevel());
    InteractionsGraph::VDescriptor ivd = indexSet.descriptor(inter);
    struct _NSLEffectOnFreeOutput nslEffectOnFreeOutput =  _NSLEffectOnFreeOutput(osnsp, *inter, indexSet.properties(ivd));
    SiconosVector & osnsp_rhs = *(*indexSet.properties(ivd).workVectors)[MoreauJeanOSI::OSNSP_RHS];
    osnsp_rhs.zero();
    inter->nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
}

void MoreauJeanGOSI::integrate(double& tinit, double& tend, double& tout, int& notUsed)
{
}

void MoreauJeanGOSI::updateState(const unsigned int)
{

  DEBUG_BEGIN("MoreauJeanGOSI::updateState(const unsigned int )\n");

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    DynamicalSystem& ds = *_dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // Get the DS type
    Type::Siconos dsType = Type::value(ds);

    // 3 - Lagrangian Systems
    if(dsType == Type::LagrangianDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&>(ds);
      bool baux = dsType == Type::LagrangianDS && useRCC && _simulation->relativeConvergenceCriterionHeld();

      SiconosVector &q = *d.q();
      SiconosVector& local_buffer = *ds_work_vectors[MoreauJeanGOSI::LOCAL_BUFFER];

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
    else if (dsType == Type::LagrangianLinearTIDS)
    {
      //LagrangianDS& d = static_cast<LagrangianDS&>(ds);
      // bool baux = dsType == Type::LagrangianDS && useRCC && _simulation->relativeConvergenceCriterionHeld();

      // SiconosVector &q = *d.q();
      // SiconosVector& local_buffer = *ds_work_vectors[MoreauJeanGOSI::LOCAL_BUFFER];

      // // Save value of q in stateTmp for future convergence computation
      // if(baux)
      //   local_buffer = q;

      updatePosition(ds);

      // if(baux)
      // {
      //   double ds_norm_ref = 1. + ds.x0()->norm2(); // Should we save this in the graph?
      //   local_buffer -= q;
      //   double aux = (local_buffer.norm2()) / ds_norm_ref;
      //   if(aux > RelativeTol)
      //     _simulation->setRelativeConvergenceCriterionHeld(false);
      // }
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      DEBUG_PRINT("MoreauJeanGOSI::updateState(const unsigned int ), dsType == Type::NewtonEulerDS \n");
      updatePosition(ds);
    }
    else THROW_EXCEPTION("MoreauJeanGOSI::updateState - not yet implemented for Dynamical system of type: " +  Type::name(ds));

  }
  DEBUG_END("MoreauJeanGOSI::updateState(const unsigned int )\n");
}


void MoreauJeanGOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanOSI OSI display ======" <<std::endl;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  if(_dynamicalSystemsGraph)
  {
    for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------" <<std::endl;
      std::cout << "--> W of dynamical system number " << ds->number() << ": " <<std::endl;
      if(_dynamicalSystemsGraph->properties(*dsi).W) _dynamicalSystemsGraph->properties(*dsi).W->display();
      else std::cout << "-> nullptr" <<std::endl;
      std::cout << "--> and corresponding theta is: " << _theta <<std::endl;
    }
  }
  std::cout << "================================" <<std::endl;
}
