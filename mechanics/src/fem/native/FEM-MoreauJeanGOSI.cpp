// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "FEM-MoreauJeanGOSI.hpp"

#include <memory>

#include "BoundaryCondition.hpp"
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
    auto& residusigmafree = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_sigma_free)];
    const auto& vold = solid->velocityMemory().getSiconosVector(0);    // vi
    const auto& sigmaold = solid->stressMemory().getSiconosVector(0);  // sigmai

    // Right-hand side q of W x = q + H r, x = [v_{k+1}; sigma_{k+1}],
    // W = [M h theta B^T; -h theta B S] (unconstrained blocks):
    //   q_v     = M v_i - h (1-theta) B^T sigma_i + h F_ext,theta
    //   q_sigma = S sigma_i + h (1-theta) B v_i
    // (the plastic term -h epsdot_p is the reaction H r of the stress interactions).
    auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
    solid->use_mass([&](const auto& M) { free_rhs.noalias() = M * vold; });
    residusigmafree.noalias() = solid->elasticityMatrix() * sigmaold;
    if (solid->hasBMatrix()) {
      free_rhs.noalias() -= time_step * (1 - _theta) * solid->BMatrix().transpose() * sigmaold;
      residusigmafree.noalias() += time_step * (1 - _theta) * solid->BMatrix() * vold;
    }
    if (solid->hasExternalForces()) {
      // computes Fext(ti)
      solid->computeFext(told);
      auto coeff = time_step * (1 - _theta);
      // free_rhs += time_step*(1-_theta) * fext(ti)
      free_rhs.noalias() += coeff * solid->fext();
      // computes Fext(ti+1)
      solid->computeFext(tend);
      coeff = time_step * _theta;
      // free_rhs += time_step*_theta * fext(ti+1)
      free_rhs.noalias() += coeff * solid->fext();
    }

    // Prescribed velocities v_j = g_j (t_{k+1}): the rows and columns j of W have been
    // zeroed (diagonal kept) and the original columns stored in
    // iterationMatrixBoundaryConditions. The known unknowns are moved to the right-hand
    // side, in the velocity AND stress rows (-h theta B(:, j) g_j), and row j gives
    // W_jj v_j = W_jj g_j. Unlike MoreauJeanOSI::applyBoundaryConditions (increment
    // formulation, Delta v = g - v_i), the global formulation solves for v_{k+1} itself.
    if (solid->boundaryConditions()) {
      auto& bc = *solid->boundaryConditions();
      bc.computePrescribedVelocity(tend);
      auto g = bc.prescribedVelocity();
      const auto& Wbc =
          *_dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions;
      auto ndof = solid->dimension();
      auto nsigma = solid->stressDimension();
      siconos::algebra::Index j = 0;
      for ([[maybe_unused]] auto idx : bc.velocityIndices()) {
        free_rhs.noalias() -= g(j) * Wbc.col(j).head(ndof);
        residusigmafree.noalias() -= g(j) * Wbc.col(j).segment(ndof, nsigma);
        j++;
      }
      j = 0;
      for (auto idx : bc.velocityIndices()) {
        free_rhs(idx) = Wbc(idx, j) * g(j);
        j++;
      }
    }

    auto& full_free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::q_sigma_free)];
    full_free_rhs.head(solid->dimension()) = free_rhs;
    full_free_rhs.tail(solid->stressDimension()) = residusigmafree;

    DEBUG_EXPR(free_rhs.display());
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
    double h = simulation()->timeStep();
    MoreauJeanOSI::_NSLEffectOnFreeOutput nslEffectOnFreeOutput{
        _NSLEffectOnFreeOutput(osnsp, *inter, indexSet.properties(ivd), _theta, h)};
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
