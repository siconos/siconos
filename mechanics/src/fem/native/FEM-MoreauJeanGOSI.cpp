/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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
#include "FEM-MoreauJeanGOSI.hpp"

#include <memory>

#include "BlockVector.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "MoreauJeanGOSI.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "SolidLinearTIDS.hpp"
#include "StressLinearTIR.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "Tools.hpp"
#include "siconos_debug.h"

void siconos::mechanics::fem::integrators::MoreauJeanGOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type
  if (auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds)) {
    // Compute W (iteration matrix)
    initializeIterationMatrix(t, solid);
    ds_work_vectors.resize(tools::enum_to_index(wk_ds::size));
    ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)] =
        std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());
    ds_work_vectors[tools::enum_to_index(wk_ds::free)] =
        std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());
    ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)] =
        std::make_shared<siconos::algebra::SiconosVector>(solid->stressDimension());
    // ds_work_vectors[tools::enum_to_index(wk_ds::sigma_free)] =
    //     std::make_shared<siconos::algebra::SiconosVector>(solid->stressDimension());
    //   // Not   used
    ds_work_vectors[tools::enum_to_index(wk_ds::buffer)] =
        std::make_shared<siconos::algebra::SiconosVector>(solid->dimension());

    ds_work_vectors[tools::enum_to_index(wk_ds::q_sigma_free)] =
        std::make_shared<siconos::algebra::SiconosVector>(solid->dimension() +
                                                          solid->stressDimension());
    solid->computeTotalForces(solid->velocity_read(), solid->q_read(), t);

    solid->swapInMemory();
  } else {
    siconos::integrators::MoreauJeanGOSI::initializeWorkVectorsForDS(t, ds);
  }
}

void siconos::mechanics::fem::integrators::MoreauJeanGOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

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

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!inter_block_work[label_xfree]) ||
        (inter_block_work[label_xfree]->numberOfBlocks() != 2))
      inter_block_work[label_xfree] = std::make_shared<siconos::algebra::BlockVector>(2);
  } else {
    if ((!inter_block_work[label_xfree]) ||
        (inter_block_work[label_xfree]->numberOfBlocks() != 1))
      inter_block_work[label_xfree] = std::make_shared<siconos::algebra::BlockVector>(1);
  }

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    auto relation =
        std::dynamic_pointer_cast<siconos::mechanics::fem::StressLinearTIR>(inter.relation());
    if (relation) {
      // Notice: this is the only place in this method where there is a difference with kernel
      // MJGOSI
      std::cout << "StressLinearTIR !!! " << std::endl;
      inter_block_work[label_xfree]->setVectorPtr(
          0, workVds1[tools::enum_to_index(wk_ds::q_sigma_free)]);

    } else {
      inter_block_work[label_xfree]->setVectorPtr(0,
                                                  workVds1[tools::enum_to_index(wk_ds::free)]);
    }
  }
  DEBUG_PRINTF("ds1->number() %i\n", ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n", ds2->number());

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      inter_block_work[label_xfree]->setVectorPtr(1,
                                                  workVds2[tools::enum_to_index(wk_ds::free)]);
    }
  }
}

double siconos::mechanics::fem::integrators::MoreauJeanGOSI::computeResidu() {
  DEBUG_PRINT("\nMoreauJeanGOSI::computeResidu(), start\n");
  // This function is used to compute the residu for each "MoreauJeanGOSI-discretized"
  // dynamical system. It then computes the norm of each of them and finally return the
  // maximum value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  auto tend = _simulation->nextTime();      // End of the time step
  auto told = _simulation->startingTime();  // Beginning of the time step
  auto time_step = tend - told;             // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  // double maxResidu = 0;
  // double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
    assert(solid);
    DEBUG_PRINT(
        "siconos::mechanics::fem::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
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

    const auto& sigma = solid->stress();  // sigmaf

    DEBUG_EXPR(qold.display(););
    DEBUG_EXPR(vold.display(););
    DEBUG_EXPR(d->q()->display(););
    DEBUG_EXPR(d->velocity()->display(););

    // --- ResiduFree computation Equation (1) ---
    auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
    residu.setZero();

    auto& residusigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)];

    siconos::algebra::SiconosVector qSigmaold{solid->dimension() + solid->stressDimension()};
    auto start = std::chrono::high_resolution_clock::now();
    const auto& vold = solid->velocityMemory().getSiconosVector(0);    // vi
    const auto& sigmaold = solid->stressMemory().getSiconosVector(0);  // sigmai
    // q_sigma = [v; sigma] at the beginning of the time step
    qSigmaold.head(solid->dimension()) = vold;
    qSigmaold.tail(solid->stressDimension()) = sigmaold;

    // Get buffer and compute full-free rhs
    auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::free)];
    auto& full_free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::q_sigma_free)];
    auto& iteration_matrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    full_free_rhs.noalias() = iteration_matrix * qSigmaold;
    free_rhs.noalias() = full_free_rhs.head(solid->dimension());
    residusigmafree.noalias() = full_free_rhs.tail(solid->stressDimension());
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout << "W*[v;sigma] time:" << float_ms.count() << " ms" << std::endl;

    if (solid->hasBMatrix()) {
      // residufree += -h*B^T sigma_i
      free_rhs.noalias() -= time_step * solid->BMatrix().transpose() * sigmaold;
      // residuSigmaFree = -h*B*v_i
      residusigmafree.noalias() -= time_step * solid->BMatrix() * vold;
    }
    if (solid->hasExternalForces()) {
      // double conditionningMagicCoeff = (time_step/solid->S()->normInf())/100000;
      double conditionningMagicCoeff = time_step;
      // computes Fext(ti)
      solid->computeFext(told);
      auto coeff = conditionningMagicCoeff * (1 - _theta);
      // free_rhs += time_step*(1-_theta) * fext(ti)
      free_rhs.noalias() += coeff * solid->fext();
      // computes Fext(ti+1)
      solid->computeFext(tend);
      coeff = conditionningMagicCoeff * _theta;
      // free_rhs += time_step*_theta * fext(ti+1)
      free_rhs.noalias() += coeff * solid->fext();
    }

    applyBoundaryConditions(*solid, free_rhs, dsi, tend, vold);
    qSigmaold.head(solid->dimension()) = free_rhs;
    qSigmaold.tail(solid->stressDimension()) = qSigmaold;
    // q_sigma = [v; sigma]

    DEBUG_EXPR(free_rhs.display());
    //        applyBoundaryConditions(*d, free_rhs, dsi, t);

    //      if (d->boundaryConditions()) {
    //        THROW_EXCEPTION(
    //            "siconos::mechanics::fem::integrators::MoreauJeanGOSI::computeResidu -
    //            boundary conditions not " "yet implemented for this type of Dynamical
    //            system\n");
    //      }

    // normResidu = 0.0;  // we assume that v_sigma = vfree_sigma + W^(-1)  [ p ; z]

    // if (normResidu > maxResidu) maxResidu = normResidu;
  }
  return 0.;  // maxResidu;
}

void siconos::mechanics::fem::integrators::MoreauJeanGOSI::updateState(const unsigned int) {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanGOSI::updateState(const unsigned int )\n");

  auto RelativeTol = _simulation->relativeConvergenceTol();
  auto useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;

    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    if (auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(
            _dynamicalSystemsGraph->bundle(*dsi))) {
      siconos::integrators::moreau_jean::updatePosition(_simulation->timeStep(), _theta,
                                                        *lltids);
    } else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(
                   _dynamicalSystemsGraph->bundle(*dsi))) {
      bool baux = useRCC && _simulation->relativeConvergenceCriterionHeld();

      auto& local_buffer = *ds_work_vectors[tools::enum_to_index(wk_ds::buffer)];

      // Save value of q in stateTmp for future convergence computation
      if (baux) local_buffer = lds->q_read();

      siconos::integrators::moreau_jean::updatePosition(_simulation->timeStep(), _theta, *lds);
      if (baux) {
        double ds_norm_ref = 1. + lds->x0().norm();  // Should we save this in the graph?
        local_buffer -= lds->q_read();
        double aux = (local_buffer.norm()) / ds_norm_ref;
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanGOSI::updateState - not yet implemented for this "
          "kind of ds.")
  }
  DEBUG_END("siconos::integrators::MoreauJeanGOSI::updateState(const unsigned int )\n");
}
