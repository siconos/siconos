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

#include "OneStepNSProblem.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "SolidLinearTIDS.hpp"
#include "StressLinearTIR.hpp"  // IWYU pragma: keep
#include "Tools.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

double siconos::mechanics::fem::integrators::MoreauJeanGOSI::computeResidu() {
  DEBUG_PRINT("\nMoreauJeanGOSI::computeResidu(), start\n");
  // This function is used to compute the residu for each "MoreauJeanGOSI-discretized"
  // dynamical system. It then computes the norm of each of them and finally return the
  // maximum value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $
  auto sim = simulation();
  auto tend = sim->nextTime();      // End of the time step
  auto told = sim->startingTime();  // Beginning of the time step
  auto time_step = tend - told;     // time step length

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
    // auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
    // residu.setZero();

    auto& residusigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)];

    siconos::algebra::SiconosVector qSigmaold{solid->dimension() + solid->stressDimension()};
    const auto& vold = solid->velocityMemory().getSiconosVector(0);    // vi
    const auto& sigmaold = solid->stressMemory().getSiconosVector(0);  // sigmai
    // q_sigma = [v; sigma] at the beginning of the time step
    qSigmaold.head(solid->dimension()) = vold;
    qSigmaold.tail(solid->stressDimension()) = sigmaold;

            // Get buffer and compute full-free rhs
    auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
    auto& full_free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::q_sigma_free)];
    auto& iteration_matrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    full_free_rhs.noalias() = iteration_matrix * qSigmaold;
    free_rhs.noalias() = full_free_rhs.head(solid->dimension());
    residusigmafree.noalias() = full_free_rhs.tail(solid->stressDimension());

    if (solid->hasBMatrix()) {
      // residufree += -h*B^T sigma_i
      free_rhs.noalias() -= time_step * solid->BMatrix().transpose() * sigmaold;
      // residuSigmaFree = -h*B*v_i
      residusigmafree.noalias() -= time_step * solid->BMatrix() * vold;
    }
    if (solid->hasExternalForces()) {
      // double conditionningMagicCoeff = (time_step/solid->S()->normInf())/100000;
      // computes Fext(ti)
      solid->computeFext(told);
      auto coeff =  time_step * (1 - _theta);
      // free_rhs += time_step*(1-_theta) * fext(ti)
      free_rhs.noalias() += coeff * solid->fext();
      // computes Fext(ti+1)
      solid->computeFext(tend);
      coeff =  time_step * _theta;
      // free_rhs += time_step*_theta * fext(ti+1)
      free_rhs.noalias() += coeff * solid->fext();
    }

    applyBoundaryConditions(*solid, free_rhs, dsi, tend, vold);
    qSigmaold.head(solid->dimension()) = free_rhs;
    qSigmaold.tail(solid->stressDimension()) = residusigmafree;
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

void siconos::mechanics::fem::integrators::MoreauJeanGOSI::NonSmoothLawContributionToOutput(
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::nonsmooth_formulations::OneStepNSProblem& osnsp) {
  if (inter->relation()->getType() == siconos::modeling::RelationType::Lagrangian ||
      inter->relation()->getType() == siconos::modeling::RelationType::NewtonEuler) {
    auto& indexSet = *osnsp.simulation()->indexSet(osnsp.indexSetLevel());
    auto ivd = indexSet.descriptor(inter);
    MoreauJeanOSI::_NSLEffectOnFreeOutput nslEffectOnFreeOutput{
        _NSLEffectOnFreeOutput(osnsp, *inter, indexSet.properties(ivd), _theta)};
    auto& osnsp_rhs =
        *(*indexSet.properties(ivd).workVectors)[tools::enum_to_index(wk_inter::osnsp_rhs)];
    osnsp_rhs.setZero();
    inter->nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
}

void siconos::mechanics::fem::integrators::MoreauJeanGOSI::integrate(double& tinit,
                                                                     double& tend,
                                                                     double& tout,
                                                                     int& notUsed) {}

void siconos::mechanics::fem::integrators::MoreauJeanGOSI::computeIteration() {}
