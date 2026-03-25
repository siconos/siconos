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

#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianSparseLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "RotationQuaternion.hpp"  // for quaternionFromTwistVector and compositionLawLieGroup
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "TypeName.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

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
  auto time_step = t - told;                // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  double maxResidu = 0;
  double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto& v_iter = *ds_work_vectors[tools::enum_to_index(wk_ds::v_iter)];

    if (auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
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

      auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
      auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
      // --- ResiduFree computation Equation (1) ---
      residu.setZero();
      auto& iterationMatrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;

      const auto& vold = lltids->velocityMemory().getSiconosVector(0);  // vi
      free_rhs = iterationMatrix * vold;

      // -- No need to update W --
      //    auto v = d->velocity();  // v = v_k,i+1
      // free_rhs += h*C*vi
      if (lltids->hasDampingMatrix()) free_rhs -= time_step * lltids->dampingMatrix() * vold;
      if (lltids->hasStiffnessMatrix()) {
        auto coeff = time_step * time_step * _theta;
        free_rhs -=
            coeff * lltids->stiffnessMatrix() * vold;  // vfree += time_step^2*_theta*K*vi
        free_rhs -= time_step * lltids->stiffnessMatrix() *
                    lltids->qMemory().getSiconosVector(0);  // vfree += time_step*K*qi
      }

      if (lltids->hasExternalForces()) {
        // computes Fext(ti)
        lltids->computeFext(told);
        auto coeff = time_step * (1 - _theta);
        // vfree -= time_step*(1-_theta) * fext(ti)
        free_rhs += coeff * lltids->fext();
        // computes Fext(ti+1)
        lltids->computeFext(t);
        coeff = time_step * _theta;
        // vfree -= time_step*_theta * fext(ti+1)
        free_rhs += coeff * lltids->fext();
      }

      DEBUG_EXPR(siconos::algebra::print(free_rhs));

      if (lltids->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      // residu = -1.0*free_rhs;
      // residu.noalias() += W **v;
      // DEBUG_EXPR(siconos::algebra::print(free_rhs));
      // if(d->p(1))
      //   residu -= *d->p(1); // Compute Residu in Workfree Notation !!

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
    }

    else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
          "Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1) -
      // h*(1-theta)*forces(ti,vi,qi) - p_i+1
      auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
      auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];

      // -- Convert the DS into a Lagrangian one.

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      // residu.setZero();

      auto& iterationMatrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
      const auto& vold = lltids->velocityMemory().getSiconosVector(0);  // vi
      free_rhs = iterationMatrix * vold;

      if (lds->hasTotalForces()) {
        // Cheaper version: get forces(ti,vi,qi) from memory
        free_rhs += time_step * (1. - _theta) * lds->forcesMemory().getSiconosVector(0);

        // Expensive computes forces(ti,vi,qi)
        // d->computeTotalForces(vold, qold, told);
        // double coef = -h * (1 - _theta);
        // // residu += coef * fL_i
        // scal(coef, *d->totalForces(), *residu, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)

        // we use residu as a local buffer to compute the current iterate in position
        auto& qold = lds->qMemory().getSiconosVector(0);
        residu = qold + time_step * (1. - _theta) * vold + time_step * _theta * v_iter;

        lds->computeTotalForces(v_iter, residu, t);
        free_rhs += time_step * _theta * lds->totalForces();

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        // auto qbasedonv =
        // std::make_shared<siconos::algebra::SiconosVector>(*qold)); *qbasedonv +=  h * ((1
        // -
        //_theta)* *vold + _theta * *v);
        // d->computeTotalForces(v, qbasedonv, t); coef = -h *
        //_theta;
        // residu += coef * fL_k,i+1
        // scal(coef, *d->totalForces(), residu, false);
      }

      residu = iterationMatrix * v_iter - free_rhs;

      DEBUG_EXPR(siconos::algebra::print(residu));

      if (lds->p(1)) residu -= lds->p_read(1);  // Compute Residu in Workfree Notation !!

      if (lds->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      DEBUG_EXPR(siconos::algebra::print(residu));
      normResidu = residu.norm();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    } else if (auto lltids =
                   std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseLinearTIDS>(
                       ds)) {
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual for v = v_i
      // otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v +
      // K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
      auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
      // --- ResiduFree computation Equation (1) ---
      residu.setZero();
      auto& iterationMatrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;

      const auto& vold = lltids->velocityMemory().getSiconosVector(0);  // vi
      free_rhs = iterationMatrix * vold;

      // -- No need to update W --
      //    auto v = d->velocity();  // v = v_k,i+1
      // free_rhs += h*C*vi
      if (lltids->hasDampingMatrix()) free_rhs -= time_step * lltids->dampingMatrix() * vold;
      if (lltids->hasStiffnessMatrix()) {
        auto coeff = time_step * time_step * _theta;
        free_rhs -=
            coeff * lltids->stiffnessMatrix() * vold;  // vfree += time_step^2*_theta*K*vi
        free_rhs -= time_step * lltids->stiffnessMatrix() *
                    lltids->qMemory().getSiconosVector(0);  // vfree += time_step*K*qi
      }

      if (lltids->hasExternalForces()) {
        // computes Fext(ti)
        lltids->computeFext(told);
        auto coeff = time_step * (1 - _theta);
        // vfree -= time_step*(1-_theta) * fext(ti)
        free_rhs += coeff * lltids->fext();
        // computes Fext(ti+1)
        lltids->computeFext(t);
        coeff = time_step * _theta;
        // vfree -= time_step*_theta * fext(ti+1)
        free_rhs += coeff * lltids->fext();
      }

      DEBUG_EXPR(siconos::algebra::print(free_rhs));

      if (lltids->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      // residu = -1.0*free_rhs;
      // residu.noalias() += W **v;
      // DEBUG_EXPR(siconos::algebra::print(free_rhs));
      // if(d->p(1))
      //   residu -= *d->p(1); // Compute Residu in Workfree Notation !!

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
    }

    else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
          "Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1) -
      // h*(1-theta)*forces(ti,vi,qi) - p_i+1
      auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
      auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];

      // -- Convert the DS into a Lagrangian one.

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      // residu.setZero();

      auto& iterationMatrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
      const auto& vold = lltids->velocityMemory().getSiconosVector(0);  // vi
      free_rhs = iterationMatrix * vold;

      if (lds->hasTotalForces()) {
        // Cheaper version: get forces(ti,vi,qi) from memory
        free_rhs += time_step * (1. - _theta) * lds->forcesMemory().getSiconosVector(0);

        // Expensive computes forces(ti,vi,qi)
        // d->computeTotalForces(vold, qold, told);
        // double coef = -h * (1 - _theta);
        // // residu += coef * fL_i
        // scal(coef, *d->totalForces(), *residu, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)

        // we use residu as a local buffer to compute the current iterate in position
        auto& qold = lds->qMemory().getSiconosVector(0);
        residu = qold + time_step * (1. - _theta) * vold + time_step * _theta * v_iter;

        lds->computeTotalForces(v_iter, residu, t);
        free_rhs += time_step * _theta * lds->totalForces();

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        // auto qbasedonv =
        // std::make_shared<siconos::algebra::SiconosVector>(*qold)); *qbasedonv +=  h * ((1
        // -
        //_theta)* *vold + _theta * *v);
        // d->computeTotalForces(v, qbasedonv, t); coef = -h *
        //_theta;
        // residu += coef * fL_k,i+1
        // scal(coef, *d->totalForces(), residu, false);
      }

      residu = iterationMatrix * v_iter - free_rhs;

      DEBUG_EXPR(siconos::algebra::print(residu));

      if (lds->p(1)) residu -= lds->p_read(1);  // Compute Residu in Workfree Notation !!

      if (lds->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      DEBUG_EXPR(siconos::algebra::print(residu));
      normResidu = residu.norm();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }

    else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeResidu(), dsType == "
          "Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) -
      // h*(1-_theta)*forces(ti,vi,qi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      DEBUG_EXPR(d->display());
      auto& residu = *ds_work_vectors[tools::enum_to_index(wk_ds::residu_free)];
      auto& free_rhs = *ds_work_vectors[tools::enum_to_index(wk_ds::vfree)];
      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"

      // Get the (constant mass matrix)
      // auto massMatrix = d->mass();
      // residu = *massMatrix * (*v - vold);
      //  // residu = M(v -
      // vold) DEBUG_EXPR(siconos::algebra::print(residu););

      auto& iterationMatrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
      const auto& vold = neds->twistMemory().getSiconosVector(0);
      // Cheaper version: get forces(ti,vi,qi) from memory
      free_rhs = iterationMatrix * vold +
                 time_step * (1 - _theta) * neds->forcesMemory().getSiconosVector(0);

      // Expensive version to check ...
      // d->computeTotalForces(vold,qold,told);
      // double coef = -h * (1.0 - _theta);
      // scal(coef, *d->totalForces(), residu, false);

      // computes forces(ti,v,q)
      // we use residu as a local buffer to compute the current iterate in position
      siconos::algebra::SiconosVector6 velocityIncrement;

      velocityIncrement = time_step * _theta * v_iter +
                          time_step * (1. - _theta) * neds->twistMemory().getSiconosVector(0);

      siconos::algebra::SiconosVector7 qtmp{7};
      qtmp.head(3) = velocityIncrement.head(3);
      siconos::geometry::quaternionFromTwistVector(velocityIncrement, qtmp);

      siconos::geometry::compositionLawLieGroup(neds->qMemory().getSiconosVector(0), qtmp);

      neds->computeWrench(v_iter, qtmp, t);
      free_rhs += time_step * _theta * neds->wrench();
      DEBUG_PRINT("siconos::integrators::MoreauJeanGOSI:: new forces :\n");
      DEBUG_EXPR(siconos::algebra::print(*d->totalForces()););
      DEBUG_EXPR(siconos::algebra::print(residu););

      if (neds->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      residu = iterationMatrix * v_iter - free_rhs;
      if (neds->p(1)) residu -= neds->p_read(1);

      if (neds->boundaryConditions()) {
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanGOSI::computeResidu - boundary conditions not "
            "yet implemented for this type of Dynamical system\n");
      }

      DEBUG_PRINT("siconos::integrators::MoreauJeanGOSI::computeResidu :\n");
      DEBUG_EXPR(siconos::algebra::print(residu););
      DEBUG_EXPR(if (neds->p(1)) siconos::algebra::print(*neds->p(1)););
      DEBUG_EXPR(siconos::algebra::print(free_rhs););

      normResidu = residu.norm();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanGOSI::computeResidu - not yet implemented for this "
          "type of DS\n");

    if (normResidu > maxResidu) maxResidu = normResidu;
  }
  return maxResidu;
}

void siconos::integrators::MoreauJeanGOSI::computeFreeState() {}

void siconos::integrators::MoreauJeanGOSI::NonSmoothLawContributionToOutput(
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::nonsmooth_formulations::OneStepNSProblem& osnsp) {
  if (inter->relation()->getType() == siconos::modeling::RelationType::Lagrangian ||
      inter->relation()->getType() == siconos::modeling::RelationType::NewtonEuler) {
    auto& indexSet = *osnsp.simulation()->indexSet(osnsp.indexSetLevel());
    auto ivd = indexSet.descriptor(inter);
    struct MoreauJeanOSI::_NSLEffectOnFreeOutput nslEffectOnFreeOutput =
        _NSLEffectOnFreeOutput(osnsp, *inter, indexSet.properties(ivd), _theta);
    auto& osnsp_rhs =
        *(*indexSet.properties(ivd).workVectors)[tools::enum_to_index(wk_inter::osnsp_rhs)];
    osnsp_rhs.setZero();
    inter->nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
}

void siconos::integrators::MoreauJeanGOSI::integrate(double& tinit, double& tend, double& tout,
                                                     int& notUsed) {}

void siconos::integrators::MoreauJeanGOSI::computeIteration() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanGOSI::computeIteration(const unsigned int )\n");

  auto RelativeTol = _simulation->relativeConvergenceTol();
  auto useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;

    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    if (auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(
            _dynamicalSystemsGraph->bundle(*dsi))) {
    } else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(
                   _dynamicalSystemsGraph->bundle(*dsi))) {
      bool baux = useRCC && _simulation->relativeConvergenceCriterionHeld();

      auto& local_buffer = *ds_work_vectors[tools::enum_to_index(wk_ds::buffer)];

      // Save value of q in stateTmp for future convergence computation
      if (baux) local_buffer = lds->q_read();

      if (baux) {
        double ds_norm_ref = 1. + lds->x0().norm();  // Should we save this in the graph?
        local_buffer -= lds->q_read();
        double aux = (local_buffer.norm()) / ds_norm_ref;
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    } else if (auto lltids =
                   std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseLinearTIDS>(
                       _dynamicalSystemsGraph->bundle(*dsi))) {
    } else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(
                   _dynamicalSystemsGraph->bundle(*dsi))) {
      bool baux = useRCC && _simulation->relativeConvergenceCriterionHeld();

      auto& local_buffer = *ds_work_vectors[tools::enum_to_index(wk_ds::buffer)];

      // Save value of q in stateTmp for future convergence computation
      if (baux) local_buffer = lds->q_read();

      if (baux) {
        double ds_norm_ref = 1. + lds->x0().norm();  // Should we save this in the graph?
        local_buffer -= lds->q_read();
        double aux = (local_buffer.norm()) / ds_norm_ref;
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    } else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(
                   _dynamicalSystemsGraph->bundle(*dsi))) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanGOSI::computeIteration(const unsigned int ), "
          "dsType "
          "== "
          "Type::NewtonEulerDS \n");
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanGOSI::computeIteration - not yet implemented for "
          "this "
          "kind of ds.")
  }
  DEBUG_END("siconos::integrators::MoreauJeanGOSI::computeIteration(const unsigned int )\n");
}

void siconos::integrators::MoreauJeanGOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanGOSI OSI display ======\n";
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  if (_dynamicalSystemsGraph) {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------\n";
      std::cout << "--> W of dynamical system number " << ds->number() << ": " << "\n";
      if (_dynamicalSystemsGraph->properties(*dsi).iterationMatrix)
        siconos::algebra::print(*_dynamicalSystemsGraph->properties(*dsi).iterationMatrix);
      else
        std::cout << "-> nullptr" << "\n";
      std::cout << "--> and corresponding theta is: " << _theta << "\n";
    }
  }
  std::cout << "================================\n";
}
