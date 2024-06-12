/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "BlockVector.hpp"
#include "BoundaryCondition.hpp"
#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for prod, subprod ...
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "../../mechanics/src/fem/native/SolidLinearTIDS.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

void siconos::integrators::MoreauJeanGOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type
  auto sods = std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds);
  // Compute W (iteration matrix)
  initializeIterationMatrixW(t, sods);
  //  )dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
  if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    ds_work_vectors.resize(siconos::integrators::MoreauJeanGOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    ds_work_vectors[siconos::integrators::MoreauJeanGOSI::LOCAL_BUFFER] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    auto dsType = siconos::types::type_value(*ds);
    if (dsType == siconos::modeling::Type::SolidLinearTIDS) {
        auto solidds = std::static_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
//        ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE]->resize(solidds->velocityDimension()+solidds->stressDimension());
        ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_SIGMAFREE] =
                std::make_shared<siconos::algebra::SiconosVector>(solidds->stressDimension());
        ds_work_vectors[siconos::integrators::MoreauJeanOSI::SIGMAFREE] =
                std::make_shared<siconos::algebra::SiconosVector>(solidds->stressDimension());
        ds_work_vectors[siconos::integrators::MoreauJeanOSI::Q_SIGMAFREE] =
                std::make_shared<siconos::algebra::SiconosVector>(solidds->dimension()+solidds->stressDimension());
    }
    lds->computeForces(t, lds->q(), lds->velocity());
    lds->swapInMemory();
  } else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
    ds_work_vectors.resize(siconos::integrators::MoreauJeanGOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(neds->dimension());
    ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(neds->dimension());
    // Compute a first value of the dotq  to store it in  _dotqMemory
    auto T = neds->T();
    auto dotq = neds->dotq();
    auto v = neds->twist();
    siconos::algebra::prod(*T, *v, *dotq, true);

    // Compute a first value of the forces to store it in _forcesMemory
    neds->computeForces(t, neds->q(), v);
    neds->swapInMemory();
  }
}

void siconos::integrators::MoreauJeanGOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            MoreauJeanGOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            MoreauJeanGOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_block_work = *interProp.workBlockVectors;

  if (!inter_work[siconos::integrators::MoreauJeanGOSI::OSNSP_RHS])
    inter_work[siconos::integrators::MoreauJeanGOSI::OSNSP_RHS] =
        std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current
  // osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  auto xfree = siconos::integrators::MoreauJeanGOSI::xfree;

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!inter_block_work[xfree]) || (inter_block_work[xfree]->numberOfBlocks() != 2))
      inter_block_work[xfree] = std::make_shared<siconos::algebra::BlockVector>(2);
  } else {
    if ((!inter_block_work[xfree]) || (inter_block_work[xfree]->numberOfBlocks() != 1))
      inter_block_work[xfree] = std::make_shared<siconos::algebra::BlockVector>(1);
  }

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    inter_block_work[xfree]->setVectorPtr(
        0, workVds1[siconos::integrators::MoreauJeanGOSI::FREE]);
  }
  DEBUG_PRINTF("ds1->number() %i\n", ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n", ds2->number());

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      inter_block_work[xfree]->setVectorPtr(
          1, workVds2[siconos::integrators::MoreauJeanGOSI::FREE]);
    }
  }
}

double siconos::integrators::MoreauJeanGOSI::computeResidu() {
  DEBUG_PRINT("\nsiconos::integrators::MoreauJeanGOSI::computeResidu(), start\n");
  // This function is used to compute the residu for each "MoreauJeanGOSI-discretized"
  // dynamical system. It then computes the norm of each of them and finally return the
  // maximum value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  auto t = _simulation->nextTime();         // End of the time step
  auto told = _simulation->startingTime();  // Beginning of the time step
  auto h = t - told;                        // time step length
  siconos::modeling::Type dsType;  // Type of the current DS.

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds;  // Current Dynamical System.

  double maxResidu = 0;
  double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    dsType = siconos::types::type_value(*ds);  // Its type
    if (dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
        auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
          "Type::LagrangianLinearTIDS\n");
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual for v = v_i
      // otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v +
      // K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      const auto& qold = d->qMemory().getSiconosVector(0);         // qi
      const auto& vold = d->velocityMemory().getSiconosVector(0);  // vi

      DEBUG_EXPR(qold.display(););
      DEBUG_EXPR(vold.display(););
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););

      auto& residu = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_FREE];
      auto& free_rhs = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE];
      // --- ResiduFree computation Equation (1) ---
      residu.zero();
      auto& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      siconos::algebra::prod(W, vold, free_rhs);

      double coeff;
      // -- No need to update W --

      auto v = d->velocity();  // v = v_k,i+1

      auto C = d->C();
      if (C) siconos::algebra::prod(h, *C, vold, free_rhs, false);  // free_rhs += h*C*vi

      auto K = d->K();
      if (K) {
        coeff = -h * h * _theta;
        siconos::algebra::prod(coeff, *K, vold, free_rhs,
                               false);                          // free_rhs += -h^2*_theta*K*vi
        siconos::algebra::prod(-h, *K, qold, free_rhs, false);  // free_rhs += -h*K*qi
      }

      auto Fext = d->fExt();
      if (Fext) {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - _theta);
        scal(coeff, *(d->fExt()), free_rhs, false);  // free_rhs += h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * _theta;
        scal(coeff, *(d->fExt()), free_rhs, false);  // free_rhs += h*_theta * fext(ti+1)
      }
      DEBUG_EXPR(free_rhs.display());
      applyBoundaryConditions(*d, free_rhs, dsi, t);

//      if (d->boundaryConditions()) {
//        THROW_EXCEPTION(
//            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
//            "yet implemented for this type of Dynamical system\n");
//      }

      // residu = -1.0*free_rhs;
      // siconos::algebra::prod(1.0, W, *v, residu, false);
      // DEBUG_EXPR(free_rhs.display());
      // if(d->p(1))
      //   residu -= *d->p(1); // Compute Residu in Workfree Notation !!

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
    }
    else     if (dsType == siconos::modeling::Type::SolidLinearTIDS) {
        DEBUG_PRINT(
            "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
            "Type::SolidLinearTIDS\n");
        // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
        // This formulae is only valid for the first computation of the residual for v = v_i
        // otherwise the complete formulae must be applied, that is
        // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v +
        // K(q_i+h(1-theta)v_i+h theta v)))
        //                     +hFext_theta     (2)
        // for v != vi, the formulae (1) is wrong.
        // in the sequel, only the equation (1) is implemented

        // Get state i (previous time step) from Memories -> var. indexed with "Old"
        auto &solid = static_cast<siconos::mechanics::fem::SolidLinearTIDS &>(*ds);

        const auto& qold = solid.qMemory().getSiconosVector(0);         // qi
        const auto& vold = solid.velocityMemory().getSiconosVector(0);  // vi
        const auto &sigmaold = solid.stressMemory().getSiconosVector(0);  // sigmai
        const auto &sigma = solid.stress(); //sigmaf


        DEBUG_EXPR(qold.display(););
        DEBUG_EXPR(vold.display(););
        DEBUG_EXPR(d->q()->display(););
        DEBUG_EXPR(d->velocity()->display(););

        auto& residu = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_FREE];
        auto& free_rhs = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE];
        auto &residusigmafree = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_SIGMAFREE];
        auto &sigmafree = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::SIGMAFREE];

        // --- ResiduFree computation Equation (1) ---
        residu.zero();
        auto& W = *_dynamicalSystemsGraph->properties(*dsi).W;
        auto qSigmaold = siconos::algebra::SiconosVector(solid.velocityDimension()+solid.stressDimension());
        auto full_free_rhs = siconos::algebra::SiconosVector(solid.velocityDimension()+solid.stressDimension());


        vold.toBlock(qSigmaold, solid.velocityDimension(), 0, 0);
        sigmaold.toBlock(qSigmaold, solid.stressDimension(), 0, solid.velocityDimension());  // q_sigma = [v; sigma]

        siconos::algebra::prod(W, qSigmaold, full_free_rhs);
        std::vector<std::size_t> subCoord(4);
        subCoord[0] = 0;
        subCoord[1] = solid.velocityDimension();
        subCoord[2] = 0;
        subCoord[3] = solid.velocityDimension();
        siconos::algebra::subscal(1.0, full_free_rhs, free_rhs,
                                  subCoord, true);
        subCoord[0] = solid.velocityDimension();
        subCoord[1] = qSigmaold.size();
        subCoord[2] = 0;
        subCoord[3] = solid.stressDimension();
        siconos::algebra::subscal(1.0, full_free_rhs, residusigmafree,
                                  subCoord, true);


        if (solid.B()) {
            siconos::algebra::prod( sigmaold, *(solid.B()), free_rhs, true); // residufree = B^T sigma_i
            siconos::algebra::scal( -h, free_rhs, free_rhs, true); // residufree = -h*B^T sigma_i
        }

        if (solid.fExt()) {
            // computes Fext(ti)
            solid.computeFExt(told);
            auto fextTheta = std::make_shared<siconos::algebra::SiconosVector>(solid.fExt()->size());
            double coeff = (1 - _theta);
            siconos::algebra::scal(coeff, *(solid.fExt()), *fextTheta,
                                   true);  // fext_k+theta = (1-_theta) * fext(ti)
            // computes Fext(ti+1)
            solid.computeFExt(t);
            coeff = _theta;
            siconos::algebra::scal(coeff, *(solid.fExt()), *fextTheta,
                                   false);  // fext_k+theta += _theta * fext(ti+1)
            double conditionningMagicCoeff = 1/solid.S()->normInf();

            siconos::algebra::scal(h/conditionningMagicCoeff, *fextTheta, free_rhs ,
                                   false);  // residufree += h*fext_k+theta
        }

        siconos::algebra::prod(-h, *solid.B(), vold, residusigmafree, true); // residuSigmaFree = -h*B*v_i

        applyBoundaryConditions(solid, free_rhs, dsi, t);
//        free = residuFree;            // copy residuFree into free
  //      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!
//        if (solid.p(1)) free += *solid.p(1);  // Compute Residu in Workfree Notation !!

//        sigfreed = residuSigfreed;
//        if (d.epsilonp(1)) {
//            sigfreed -= *d.epsilonp(1)*h;
//        }



        DEBUG_EXPR(free_rhs.display());
//        applyBoundaryConditions(*d, free_rhs, dsi, t);

  //      if (d->boundaryConditions()) {
  //        THROW_EXCEPTION(
  //            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
  //            "yet implemented for this type of Dynamical system\n");
  //      }

        // residu = -1.0*free_rhs;
        // siconos::algebra::prod(1.0, W, *v, residu, false);
        // DEBUG_EXPR(free_rhs.display());
        // if(d->p(1))
        //   residu -= *d->p(1); // Compute Residu in Workfree Notation !!

        normResidu = 0.0;  // we assume that v_sigma = vfree_sigma + W^(-1)  [ p ; z]
      }


    else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
          "Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1) -
      // h*(1-theta)*forces(ti,vi,qi) - p_i+1
      auto& residu = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_FREE];
      auto& free_rhs = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE];

      // -- Convert the DS into a Lagrangian one.

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      // const auto& qold =
      // d->qMemory().getsiconos::algebra::SiconosVector(0);
      const auto& vold = lds->velocityMemory().getSiconosVector(0);
      auto q = lds->q();

      auto v = lds->velocity();  // v = v_k,i+1
      // residu.zero();
      DEBUG_EXPR(residu.display());

      // DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());
      DEBUG_EXPR(q->display());
      DEBUG_EXPR(v->display());
      free_rhs.zero();
      auto& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      siconos::algebra::prod(W, vold, free_rhs);

      if (lds->forces()) {
        // Cheaper version: get forces(ti,vi,qi) from memory
        const auto& fold = lds->forcesMemory().getSiconosVector(0);
        double coef = h * (1 - _theta);
        scal(coef, fold, free_rhs, false);

        // Expensive computes forces(ti,vi,qi)
        // d->computeForces(told, qold, vold);
        // double coef = -h * (1 - _theta);
        // // residu += coef * fL_i
        // scal(coef, *d->forces(), *residu, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        lds->computeForces(t, q, v);
        coef = h * _theta;
        scal(coef, *lds->forces(), free_rhs, false);

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        // auto qbasedonv =
        // std::make_shared<siconos::algebra::SiconosVector>(*qold)); *qbasedonv +=  h * ((1
        // -
        //_theta)* *vold + _theta * *v); d->computeForces(t, qbasedonv, v); coef = -h *
        //_theta;
        // residu += coef * fL_k,i+1
        // scal(coef, *d->forces(), residu, false);
      }

      residu = -1.0 * free_rhs;
      siconos::algebra::prod(1.0, W, *v, residu, false);

      DEBUG_EXPR(residu.display());

      if (lds->p(1)) residu -= *lds->p(1);  // Compute Residu in Workfree Notation !!

      if (lds->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      DEBUG_EXPR(residu.display());
      normResidu = residu.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
          "Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) -
      // h*(1-_theta)*forces(ti,vi,qi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      DEBUG_EXPR(d->display());
      auto& residu = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::RESIDU_FREE];
      auto& free_rhs = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE];
      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"
      const auto& vold = d->twistMemory().getSiconosVector(0);

      // Get the current state vector
      auto q = d->q();
      auto v = d->twist();  // v = v_k,i+1

      DEBUG_EXPR(vold.display());
      DEBUG_EXPR(q->display());
      DEBUG_EXPR(v->display());

      residu.zero();
      // Get the (constant mass matrix)
      // auto massMatrix = d->mass();
      // siconos::algebra::prod(*massMatrix, (*v - vold), residu, true); // residu = M(v -
      // vold) DEBUG_EXPR(residu.display(););
      free_rhs.zero();
      auto& W = *_dynamicalSystemsGraph->properties(*dsi).W;
      siconos::algebra::prod(W, vold, free_rhs);

      if (d->forces())  // if fL exists
      {
        DEBUG_PRINTF("siconos::integrators::MoreauJeanGOSI:: _theta = %e\n", _theta);
        DEBUG_PRINTF("siconos::integrators::MoreauJeanGOSI:: h = %e\n", h);

        // Cheaper version: get forces(ti,vi,qi) from memory
        const auto& fold = d->forcesMemory().getSiconosVector(0);
        auto coef = h * (1 - _theta);
        scal(coef, fold, free_rhs, false);

        // Expensive version to check ...
        // d->computeForces(told,qold,vold);
        // double coef = -h * (1.0 - _theta);
        // scal(coef, *d->forces(), residu, false);

        DEBUG_PRINT("siconos::integrators::MoreauJeanGOSI:: old forces :\n");
        DEBUG_EXPR(d->forces()->display(););
        DEBUG_EXPR(residu.display(););

        // computes forces(ti,v,q)
        d->computeForces(t, q, v);
        coef = h * _theta;
        scal(coef, *d->forces(), free_rhs, false);
        DEBUG_PRINT("siconos::integrators::MoreauJeanGOSI:: new forces :\n");
        DEBUG_EXPR(d->forces()->display(););
        DEBUG_EXPR(residu.display(););
      }

      if (d->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      residu = -1.0 * free_rhs;

      siconos::algebra::prod(1.0, W, *v, residu, false);

      if (d->p(1)) residu -= *d->p(1);

      if (d->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      DEBUG_PRINT("siconos::integrators::MoreauJeanGOSI::computeResidu :\n");
      DEBUG_EXPR(residu.display(););
      DEBUG_EXPR(if (d->p(1)) d->p(1)->display(););
      DEBUG_EXPR(free_rhs.display(););

      normResidu = residu.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanGOSI::computeResidu - not yet implemented for this "
          "type of DS\n");

    if (normResidu > maxResidu) maxResidu = normResidu;
  }
  return maxResidu;
}

void siconos::integrators::MoreauJeanGOSI::applyBoundaryConditions(
    siconos::modeling::SecondOrderDS &d,
    siconos::algebra::SiconosVector &residu,
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, double t) {

  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(...)\n");
  if (d.boundaryConditions()) {
    d.boundaryConditions()->computePrescribedVelocity(t);

    unsigned int columnindex = 0;
    auto &WBoundaryConditions =
        *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions;
    auto columntmp =
        std::make_shared<siconos::algebra::SiconosVector>(d.dimension());
    auto dsType = siconos::types::type_value(d);
    if (dsType == siconos::modeling::Type::SolidLinearTIDS){
        auto &solid = static_cast<siconos::mechanics::fem::SolidLinearTIDS &>(d);
        columntmp->resize(solid.dimension() + solid.stressDimension());
    }
    for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
        std::cout << "Iterating BCs with GOSI: " << itindex << " val: " << d.boundaryConditions()->prescribedVelocity()->getValue(columnindex) <<std::endl;
      double DeltaPrescribedVelocity =
          d.boundaryConditions()->prescribedVelocity()->getValue(columnindex);
      DEBUG_PRINTF(
          "index  = %i, value = %e\n", *itindex,
          d.boundaryConditions()->prescribedVelocity()->getValue(columnindex));
      DEBUG_PRINTF("DeltaPrescribedVelocity = %e\n", DeltaPrescribedVelocity);
      std::cout << "columnindex: " << columnindex << " columntmp size: " << columntmp->size() << " WBoundaryConditions size: " << WBoundaryConditions.size(0) << " " << WBoundaryConditions.size(1) << std::endl;
      WBoundaryConditions.getCol(columnindex, *columntmp);
      residu -= *columntmp * (DeltaPrescribedVelocity);

      residu.setValue(
          itindex, -columntmp->getValue(itindex) * (DeltaPrescribedVelocity));

      columnindex++;
    }
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(...)\n");
}

void siconos::integrators::MoreauJeanGOSI::computeFreeState() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanGOSI::computeFreeState()\n");
  DEBUG_END("siconos::integrators::MoreauJeanGOSI::computeFreeState()\n");
}

void siconos::integrators::MoreauJeanGOSI::NonSmoothLawContributionToOutput(
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::nonsmooth_formulations::OneStepNSProblem& osnsp) {
  if (inter->relation()->getType() == siconos::modeling::RelationType::Lagrangian ||
      inter->relation()->getType() == siconos::modeling::RelationType::NewtonEuler) {
    auto& indexSet = *osnsp.simulation()->indexSet(osnsp.indexSetLevel());
    auto ivd = indexSet.descriptor(inter);
    struct MoreauJeanOSI::_NSLEffectOnFreeOutput nslEffectOnFreeOutput =
        _NSLEffectOnFreeOutput(osnsp, *inter, indexSet.properties(ivd));
    auto& osnsp_rhs = *(*indexSet.properties(ivd).workVectors)[MoreauJeanOSI::OSNSP_RHS];
    osnsp_rhs.zero();
    inter->nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
}

void siconos::integrators::MoreauJeanGOSI::integrate(double& tinit, double& tend, double& tout,
                                                     int& notUsed) {}

void siconos::integrators::MoreauJeanGOSI::updateState(const unsigned int) {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanGOSI::updateState(const unsigned int )\n");

  auto RelativeTol = _simulation->relativeConvergenceTol();
  auto useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    if (std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
      // LagrangianDS& d = static_cast<LagrangianDS&>(ds);
      //  bool baux = dsType == Type::LagrangianDS && useRCC &&
      //  _simulation->relativeConvergenceCriterionHeld();

      // auto &q = *d.q();
      // auto& local_buffer =
      // *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::LOCAL_BUFFER];

      // // Save value of q in stateTmp for future convergence computation
      // if(baux)
      //   local_buffer = q;

      updatePosition(*ds);

      // if(baux)
      // {
      //   double ds_norm_ref = 1. + ds.x0()->norm2(); // Should we save this in the graph?
      //   local_buffer -= q;
      //   double aux = (local_buffer.norm2()) / ds_norm_ref;
      //   if(aux > RelativeTol)
      //     _simulation->setRelativeConvergenceCriterionHeld(false);
      // }
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      bool baux = useRCC && _simulation->relativeConvergenceCriterionHeld();

      auto& q = *d->q();
      auto& local_buffer =
          *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::LOCAL_BUFFER];

      // Save value of q in stateTmp for future convergence computation
      if (baux) local_buffer = q;

      updatePosition(*ds);

      if (baux) {
        double ds_norm_ref = 1. + ds->x0()->norm2();  // Should we save this in the graph?
        local_buffer -= q;
        double aux = (local_buffer.norm2()) / ds_norm_ref;
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    } else if (std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::updateState(const unsigned int ), dsType "
          "== "
          "Type::NewtonEulerDS \n");
      updatePosition(*ds);
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanGOSI::updateState - not yet implemented for this "
          "kind of ds.")
  }
  DEBUG_END("siconos::integrators::MoreauJeanGOSI::updateState(const unsigned int )\n");
}

void siconos::integrators::MoreauJeanGOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanOSI OSI display ======\n";
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  if (_dynamicalSystemsGraph) {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------\n";
      std::cout << "--> W of dynamical system number " << ds->number() << ": "
                << "\n";
      if (_dynamicalSystemsGraph->properties(*dsi).W)
        _dynamicalSystemsGraph->properties(*dsi).W->display();
      else
        std::cout << "-> nullptr"
                  << "\n";
      std::cout << "--> and corresponding theta is: " << _theta << "\n";
    }
  }
  std::cout << "================================\n";
}
