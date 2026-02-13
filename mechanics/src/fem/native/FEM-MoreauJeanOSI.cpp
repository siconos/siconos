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
#include "FEM-MoreauJeanOSI.hpp"

#include <memory>

#include "BoundaryCondition.hpp"
#include "MoreauJeanOSI.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "SolidLinearTIDS.hpp"
#include "StressLinearTIR.hpp"
// // #define DEBUG_NOCOLOR
// // #define DEBUG_STDOUT
// // #define DEBUG_MESSAGES
// // #define DEBUG_BEGIN_END_ONLY
// // #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::integrators::MoreauJeanOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // Check dynamical system type
  auto solid = std::dynamic_pointer_cast<SolidLinearTIDS>(ds);
  assert(solid);
  initializeIterationMatrix(t, solid);

  // Initialize work vectors
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(tools::enum_to_index(wk_ds::size));

  ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());
  ds_work_vectors[tools::enum_to_index(wk_ds::vfree)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());
  ds_work_vectors[tools::enum_to_index(wk_ds::buffer)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());
  ds_work_vectors[tools::enum_to_index(wk_ds::v_iter)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());
  ds_work_vectors[tools::enum_to_index(wk_ds::sigma_iter)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->stressDimension());

  ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->stressDimension());
  ds_work_vectors[tools::enum_to_index(wk_ds::sigma_free)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->stressDimension());
  ds_work_vectors[tools::enum_to_index(wk_ds::q_sigma_free)] =
      std::make_shared<siconos::algebra::SiconosVector>(solid->dimension() +
                                                        solid->stressDimension());
  // Update dynamical system components (for memory swap).

  solid->computeTotalForces(solid->velocity_read(), solid->q_read(), t);
  solid->swapInMemory();
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  assert(ds1);
  auto is_ds1_integrated_by_this_osi = checkOSI(DSG.descriptor(ds1));

  auto ds2 = interProp.target;
  assert(ds2);
  auto use_two_ds = ds1 != ds2;
  bool is_ds2_integrated_by_this_osi = use_two_ds && checkOSI(DSG.descriptor(ds2));

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            tools::enum_to_index(siconos::integrators::MoreauJeanOSI::wk_inter::size));
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            tools::enum_to_index(siconos::integrators::MoreauJeanOSI::wkb_inter::size));
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_block_work = *interProp.workBlockVectors;

  inter.relation()->checkSize(inter);

  if (!inter_work[tools::enum_to_index(
          siconos::integrators::MoreauJeanOSI::wk_inter::osnsp_rhs)])
    inter_work[tools::enum_to_index(
        siconos::integrators::MoreauJeanOSI::wk_inter::osnsp_rhs)] =
        std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current
  // osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  auto label_xfree =
      tools::enum_to_index(siconos::integrators::MoreauJeanOSI::wkb_inter::xfree);

  if (use_two_ds) {
    if ((!inter_block_work[label_xfree]) ||
        (inter_block_work[label_xfree]->numberOfBlocks() != 2))
      inter_block_work[label_xfree] = std::make_shared<siconos::algebra::BlockVector>(2);
  } else {
    if ((!inter_block_work[label_xfree]) ||
        (inter_block_work[label_xfree]->numberOfBlocks() != 1))
      inter_block_work[label_xfree] = std::make_shared<siconos::algebra::BlockVector>(1);
  }
  auto relation =
      std::dynamic_pointer_cast<siconos::mechanics::fem::StressLinearTIR>(inter.relation());

  if (is_ds1_integrated_by_this_osi) {
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    if (relation) {
      // Notice: this is the only place in this method where there is a difference with kernel
      // MJGOSI
      std::cout << "StressLinearTIR !!! " << std::endl;
      inter_block_work[label_xfree]->setVectorPtr(
          0, workVds1[tools::enum_to_index(wk_ds::q_sigma_free)]);

    } else {
      inter_block_work[label_xfree]->setVectorPtr(
          0, workVds1[tools::enum_to_index(wk_ds::vfree)]);
    }
  }
  if (is_ds2_integrated_by_this_osi) {
    DEBUG_PRINTF("ds2->number() %i is taken into account\n", ds2->number());
    assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
    auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relation) {
      // Notice: this is the only place in this method where there is a difference with kernel
      // MJGOSI
      std::cout << "StressLinearTIR !!! " << std::endl;
      inter_block_work[label_xfree]->setVectorPtr(
          1, workVds2[tools::enum_to_index(wk_ds::q_sigma_free)]);

    } else {
      inter_block_work[label_xfree]->setVectorPtr(
          1, workVds2[tools::enum_to_index(wk_ds::vfree)]);
    }
  }
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::initializeIterationMatrix(
    double time, std::shared_ptr<siconos::modeling::SecondOrderDS> ds) const {
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical
  // system, depending on its type.
  assert(ds);
  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "MoreauJeanOSI::initializeIterationMatrix(t,ds) "
        "- ds does not "
        "belong to the OSI.");

  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);

  /* For SolidLinearTIDS systems, we define W by blocks: W = [M h*theta*B^T; h*theta*B -S]
   so that the discretised system is W v_{\sigma} = RHS:
   $ \begin{bmatrix}
        \boldsymbol M & h\theta {\boldsymbol B}^T \\
        h\theta \boldsymbol B & - \boldsymbol S
        \end{bmatrix}  \begin{bmatrix} v_{k+\theta} \\
        \sigma_{k+\theta}
        \end{bmatrix}}
        = \begin{bmatrix}
        \boldsymbol M v_{k} + h\theta F_{ext,k+\theta} + \theta \boldsymbol H p_{N,k+1}\\
        -\boldsymbol S \sigma_{k} -h\theta z_{k+\theta}
        \end{bmatrix}$
  */

  auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
  assert(solid);

  // Allocate iteration matrix (in the graph)
  assert(!_dynamicalSystemsGraph->properties(dsv).iterationMatrix);
  auto size_iteration_mat = solid->real_size();
  _dynamicalSystemsGraph->properties(dsv).iterationMatrix =
      std::make_shared<siconos::algebra::SiconosMatrix>(size_iteration_mat,
                                                        size_iteration_mat);
  // Get a link to the iteration matrix
  auto iterationMat = _dynamicalSystemsGraph->properties(dsv).iterationMatrix;

  // Start computation
  double timeStep = _simulation->timeStep();
  double htheta = timeStep * _theta;
  double conditionningMagicCoeff = 1.0;

  // W = [M 0; 0 0]
  iterationMat->topLeftCorner(solid->dimension(), solid->dimension()) =
      conditionningMagicCoeff * solid->mass();

  // W = [M 0; 0 -S]
  iterationMat->bottomRightCorner(solid->stressDimension(), solid->stressDimension()) =
      -conditionningMagicCoeff * solid->elasticityMatrix();

  iterationMat->bottomLeftCorner(solid->stressDimension(), solid->dimension()) =
      htheta * solid->BMatrix();
  // W = [M h*theta*B^T; h*theta*B -S]
  iterationMat->topRightCorner(solid->dimension(), solid->stressDimension()) =
      htheta * solid->BMatrix().transpose();
  // LU factorization
  _dynamicalSystemsGraph->properties(dsv).LUW =
      std::make_shared<siconos::algebra::SiconosDenseLUMatrix>(*iterationMat);

  // Apply BC
  if (solid->boundaryConditions()) initializeIterationMatrixBoundaryConditions(*solid, dsv);
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::computeInitialNewtonState() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeInitialNewtonState()\n");
  // Compute the position value giving the initial velocity.
  // The goal is to save one newton iteration for nearly linear system
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    // Copy current velocity in V_ITER
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto& v_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::v_iter)];
    auto& sigma_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::sigma_iter)];

    auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
    assert(solid);
    v_iter = solid->velocity_read();
    sigma_iter = solid->stress_read();
  }
}

double siconos::mechanics::fem::integrators::MoreauJeanOSI::computeResidu() {
  DEBUG_BEGIN("MoreauJeanOSI::computeResidu()\n");
  // This function is used to compute the residu for each
  // "MoreauJeanOSI-discretized" dynamical system. It then computes the norm of
  // each of them and finally return the maximum value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) -
  //  h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) -
  //  h(1-\theta)f(x_k,t_k) $

  double t = _simulation->nextTime();         // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double time_step = t - told;                // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", time_step);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.

  // double maxResidu = 0;
  // double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    DEBUG_PRINT(
        "MoreauJeanOSI::computeResidu(), dsType == "
        "siconos::modeling::Type::LagrangianLinearTIDS\n");
    // Here, for TIDS, computeResidu computes the Right hand side.
    // We compute the top block of the RHS: ResiduFree = -h B^T sigma_i + h Fext_theta
    // and residuSigmaFree the bottom block: ResiduSigmaFree = -h B v_i
    // This way we have W [v_{i+1} - v_i; sigma_{i+1} - sigma_i] = [ResiduFree;
    // ResiduSigmaFree]

    // -- Convert the DS into a SolidLinearTIDS one.
    auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
    assert(solid);

    auto& residuFree = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
    residuFree.setZero();
    if (solid->hasBMatrix()) {
      const auto& sigmaold = solid->stressMemory().getSiconosVector(0);  // sigmai
      residuFree.noalias() = -time_step * solid->BMatrix().transpose() * sigmaold;
    }

    if (solid->hasExternalForces()) {
      double conditionningMagicCoeff =
          time_step / siconos::algebra::normInf(solid->elasticityMatrix());

      // fext_k+theta += (1-theta) * fext(ti) + theta * fext(ti+1)

      // computes Fext(ti)
      solid->computeFext(told);
      auto coeff = conditionningMagicCoeff * (1 - _theta);
      residuFree += coeff * solid->fext();
      // computes Fext(ti+1)
      solid->computeFext(t);
      coeff = conditionningMagicCoeff * _theta;
      residuFree += coeff * solid->fext();
    }

    auto& residuSigfreed = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)];
    const auto& vold = solid->velocityMemory().getSiconosVector(0);  // vi
    residuSigfreed = -time_step * solid->BMatrix() * vold;

    applyBoundaryConditions(*solid, residuFree, dsi, t, vold);

    auto& free = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
    free = residuFree;  // copy residuFree into free

    //      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!
    if (solid->p(1)) free += solid->p_read(1);  // Compute Residu in Workfree Notation !!

    auto& sigfreed = *ds_work_vectors[tools::enum_to_index(wk_ds::sigma_free)];
    sigfreed = residuSigfreed;

    sigfreed -= solid->epsilon_p(1) * time_step;

    // We use free as tmp buffer
    DEBUG_EXPR(free.display());
    DEBUG_EXPR(residuFree.display());

    // normResidu = 0.0;  // we assume that v_sigma = vfree_sigma + W^(-1)  [ p ; z]
    //      normResidu = realresiduFree->norm2();
  }
  DEBUG_END("MoreauJeanOSI::computeResidu()\n");
  return 0.;  // maxResidu;
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::computeFreeState() {
  DEBUG_BEGIN("MoreauJeanOSI::computeFreeState()\n");
  // This function computes "free" states of the DS belonging to this
  // Integrator. "Free" means without taking non-smooth effects into account.

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  auto *rold =
  //  static_cast<SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  // std::shared_ptr<siconos::modeling::DynamicalSystem> ds; // Current
  // Dynamical System. std::shared_ptr<siconos::algebra::SiconosMatrix> W; // W
  // MoreauJeanOSI matrix of the current DS.

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    // IN to be updated at current time: W, M, q, v, fL
    // IN at told: qi,vi, fLi

    // Note: indices i/i+1 corresponds to value at the beginning/end of the time
    // step. Index k stands for Newton iteration and thus corresponds to the
    // last computed value, ie the one saved in the DynamicalSystem. "i" values
    // are saved in memory vectors.

    // vFree = v_k,i+1 - W^{-1} ResiduFree
    // with
    // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - time_step*theta*forces(t,v_k,i+1,
    // q_k,i+1) - time_step*(1-theta)*forces(ti,vi,qi)

    // Get storage for W from the graph
    // Get workVectors (rfree, vfree ...) from the graph
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // Current dynamical system
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    auto solid = std::dynamic_pointer_cast<SolidLinearTIDS>(ds);
    assert(solid);

    //        auto &residuFree =
    //            *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
    //        auto &vfree = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];

    auto& qSigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::q_sigma_free)];
    auto& residuFree = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
    qSigmafree.head(solid->dimension()) = residuFree;
    auto& residusigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)];
    qSigmafree.tail(solid->stressDimension()) = residusigmafree;
    // q_sigma_free = [residuFree; residuSigmaFree]

    // -- sigmafree =   W^{-1} q_sigma --

    // -> Solve WX = sigmafree and set sigmafree = X
    qSigmafree = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(qSigmafree);

    // vfree += vold + [W^{-1} * qsigmafree](0:size V)
    auto& vfree = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
    vfree = solid->velocityMemory().getSiconosVector(0) + qSigmafree.head(solid->dimension());

    // sigmafree += sigmaold + [W^{-1} * qsigmafree](sizeV:end)
    auto& sigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::sigma_free)];
    sigmafree =
        solid->stressMemory().getSiconosVector(0) + qSigmafree.tail(solid->stressDimension());

    // We don't need to do vfree = -vfree: we computed the
    // real right hand side in computeResidu(), and not its opposite like it is done for
    // generanl lagrangianDS case
  }
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::integrate(double& tinit,
                                                                    double& tend, double& tout,
                                                                    int& notUsed) {
  // Last parameter is not used (required for LsodarOSI but not for
  // MoreauJeanOSI).

  double time_step = tend - tinit;
  tout = tend;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    // get the ds
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto solid = std::dynamic_pointer_cast<SolidLinearTIDS>(ds);

    // [velocity ; stress] computation :
    //
    // [v; sigma] = W^{-1} [-h B^T \sigma_i + h F_{ext,i+\theta} ; -h B v_i] +
    // W^{-1}*[p_{i+1}; h \dot \epsilon_{i+\theta}]
    //
    // get velocity and stress pointers for current time step
    auto v = solid->velocity();
    auto sigma = solid->stress();

    *v = solid->p_read(1);
    *sigma = -time_step * solid->plasticRate_read();
    const auto& vold = solid->velocityMemory().getSiconosVector(0);
    const auto& sigmaOld = solid->stressMemory().getSiconosVector(0);

    if (solid->hasBMatrix()) {
      // v += -h*B^T*sigma_i
      *v -= time_step * solid->BMatrix().transpose() * sigmaOld;

      // sigma += -h*B*vi
      *sigma -= time_step * solid->BMatrix() * vold;
    }

    // -- No need to update W --
    if (solid->hasExternalForces()) {
      // computes Fext(ti)
      solid->computeFext(tinit);
      auto coeff = time_step * (1 - _theta);
      *v += coeff * solid->fext();
      // computes Fext(ti+1)
      solid->computeFext(tout);
      coeff = time_step * _theta;
      *v += coeff * solid->fext();
    }

    siconos::algebra::SiconosVector qSigma{v->size() + sigma->size()};
    qSigma.head(solid->dimension()) = *v;
    qSigma.tail(solid->stressDimension()) = *sigma;
    // q_sigma = [v; sigma]

    // -> Solve WX = qSigma and set qSigma = X
    qSigma = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(qSigma);
    *v += vold + qSigma.head(solid->dimension());
    *sigma += sigmaOld + qSigma.tail(solid->stressDimension());
  }
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::updateState(const unsigned int) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto& v_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::v_iter)];
    auto& sigma_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::sigma_iter)];

    auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(
        _dynamicalSystemsGraph->bundle(*dsi));
    assert(solid);
    // Update ds variables with the current iteration state
    *solid->velocity() = v_iter;
    *solid->stress() = sigma_iter;
    auto time_step = _simulation->timeStep();
    *solid->q() = time_step * _theta * solid->velocity_read() +
                  time_step * (1. - _theta) * solid->velocityMemory().getSiconosVector(0) +
                  solid->qMemory().getSiconosVector(0);
  }
}

void siconos::mechanics::fem::integrators::MoreauJeanOSI::computeIteration() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeIteration()\n");

  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);
  double time_step = _simulation->timeStep();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
    assert(solid);

    auto& v_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::v_iter)];
    auto& sigma_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::sigma_iter)];

    bool hasContact = (solid->p(_levelMaxForInput) && solid->p(_levelMaxForInput)->size() > 0);
    bool hasPlasticity = (solid->has_epsilon_p(_levelMaxForInput) &&
                          solid->epsilon_p(_levelMaxForInput).size() > 0);

    if (hasContact) {
      assert(((solid->p(_levelMaxForInput))) &&
             " siconos::integrators::MoreauJeanOSI::ComputeIteration() "
             "*solid->p(_levelMaxForInput) == nullptr.");
      v_iter = solid->p_read(_levelMaxForInput);  // v = p
      if (solid->boundaryConditions()) {
        for (auto itindex : solid->boundaryConditions()->velocityIndices()) {
          (v_iter)(itindex) = 0.0;
        }
      }

    } else {
      v_iter.setZero();
    }

    if (hasPlasticity) {
      sigma_iter = -time_step * solid->epsilon_p(_levelMaxForInput);
    } else {
      sigma_iter.setZero();
    }

    siconos::algebra::SiconosVector qSigma;
    if (hasContact || hasPlasticity) {
      qSigma.resize(solid->real_size());
      qSigma.head(solid->dimension()) = v_iter;
      qSigma.tail(solid->stressDimension()) = sigma_iter;
      // q_sigma = [v; sigma]
      // -> Solve WX = qSigma and set qSigma = X
      qSigma = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(qSigma);
      v_iter = qSigma.head(solid->dimension());
      sigma_iter = qSigma.tail(solid->stressDimension());
    }

    auto& vfree = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
    v_iter += vfree;

    auto& sigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::sigma_free)];
    sigma_iter += sigmafree;

    if (solid->boundaryConditions()) {
      int bc = 0;
      siconos::algebra::SiconosVector columntmp22{22};

      for (auto itindex : solid->boundaryConditions()->velocityIndices()) {
        auto columntmp =
            _dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions->col(
                bc);
        /*\warning we assume that W is symmetric in the Lagrangian case*/

        double value = 0.;
        if (hasContact || hasPlasticity)  // qSigma has been computed
          value = -1. * columntmp.dot(qSigma);
        if (hasContact) {
          value += solid->p_read(_levelMaxForInput)(itindex);
        }
        /* \warning the computation of reactionToBoundaryConditions take into
           account the contact impulse but not the external and internal
           forces. A complete computation of the residu should be better */
        (*solid->reactionToBoundaryConditions())(bc) = value;
        bc++;
      }
    }
  }
}

bool siconos::mechanics::fem::integrators::MoreauJeanOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction>"
      " inter, "
      "unsigned int i)\n");

  assert(i == 1);
  auto relation =
      std::dynamic_pointer_cast<siconos::mechanics::fem::StressLinearTIR>(inter->relation());
  if (relation) {
    std::cout << "Stress interaction !!!  " << std::endl;
    return true;
  } else {
    return siconos::integrators::MoreauJeanOSI::addInteractionInIndexSet(inter, i);
  }
}
