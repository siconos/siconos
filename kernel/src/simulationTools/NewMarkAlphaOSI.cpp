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
#include "NewMarkAlphaOSI.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosConst.hpp"  // For MACHINE_PREC
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
siconos::integrators::NewMarkAlphaOSI::NewMarkAlphaOSI(double new_beta, double new_gamma,
                                                       double new_alpha_m, double new_alpha_f,
                                                       bool flag = false)
    : OneStepIntegrator(siconos::integrators::IntegratorType::NEWMARKALPHAOSI, 1, 0, 2, 1, 2),
      _alpha_m{new_alpha_m},
      _alpha_f{new_alpha_f},
      _gamma{new_gamma},
      _beta{new_beta},
      _IsVelocityLevel{flag} {}

siconos::integrators::NewMarkAlphaOSI::NewMarkAlphaOSI(double rho_infty, bool flag = false)
    : OneStepIntegrator(siconos::integrators::IntegratorType::NEWMARKALPHAOSI, 1, 0, 2, 1, 2) {
  _alpha_m = (2 * rho_infty - 1) / (rho_infty + 1);
  _alpha_f = rho_infty / (rho_infty + 1);
  _gamma = 0.5 + _alpha_f - _alpha_m;
  _beta = 0.25 * std::pow((_gamma + 0.5), 2);
  _IsVelocityLevel = flag;
}

void siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrix(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds);

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrix(t,ds) - ds does "
        "not belong to the OSI.");

  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
  assert(!_dynamicalSystemsGraph->properties(dsv).iterationMatrix);

  // Allocate storage for W in the graph
  _dynamicalSystemsGraph->properties(dsv).iterationMatrix =
      std::make_shared<siconos::algebra::SiconosMatrix>(ds->dimension(), ds->dimension());

  if (auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    auto gamma_prime = _gamma / _beta;
    auto inv_time_step =
        1. / (_simulation->nextTime() - _simulation->startingTime());  // 1. / step size
    auto beta_prime =
        (1. - _alpha_m) / ((1. - _alpha_f) * _beta) * inv_time_step * inv_time_step;

    if (lltids->hasMass())
      *_dynamicalSystemsGraph->properties(dsv).iterationMatrix = lltids->mass();
    else
      _dynamicalSystemsGraph->properties(dsv).iterationMatrix->setIdentity();

    *_dynamicalSystemsGraph->properties(dsv).iterationMatrix *= beta_prime;

    if (lltids->hasDampingMatrix())
      *_dynamicalSystemsGraph->properties(dsv).iterationMatrix +=
          gamma_prime * inv_time_step * lltids->dampingMatrix();
    if (lltids->hasStiffnessMatrix())
      *_dynamicalSystemsGraph->properties(dsv).iterationMatrix += lltids->stiffnessMatrix();

    // LU Factorisation
    _dynamicalSystemsGraph->properties(dsv).LUW =
        std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(
            *_dynamicalSystemsGraph->properties(dsv).iterationMatrix);

  } else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    computeIterationMatrix(
        ds,
        *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
             .iterationMatrix,
        _dynamicalSystemsGraph->properties(dsv).LUW);
  } else {
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrix(t,ds), this "
        "type of Dynamical System is not yet implemented");
  }
}

void siconos::integrators::NewMarkAlphaOSI::computeIterationMatrix(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    siconos::algebra::SiconosMatrix& iterationMatrix,
    std::shared_ptr<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>> luw) {
  // We assume that the iteration matrix is properly allocated (call to
  // initializeIterationMatrix)

  if (auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    return;  // Time invariant system, W constant already up to date.
  }

  auto gamma_prime = _gamma / _beta;
  auto time_step = _simulation->nextTime() - _simulation->startingTime();  // step size
  auto inv_time_step = 1. / time_step;
  auto beta_prime =
      (1. - _alpha_m) / ((1. - _alpha_f) * _beta) * inv_time_step * inv_time_step;

  if (time_step < 100 * siconos::internal::MACHINE_PREC)
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::computeIterationMatrix(t,ds), "
        "integration time step is too small");

  if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    if (lds->hasMass()) {
      lds->computeMass(lds->q_read());
      iterationMatrix = lds->mass();
    } else
      iterationMatrix.setIdentity();

    iterationMatrix *= beta_prime;

    if (lds->hasJacobianTotalForces()) {
      iterationMatrix -= gamma_prime * inv_time_step * lds->jacobianTotalForcesOver_velocity();
      iterationMatrix -= lds->jacobianTotalForcesOver_q();
    }

    // LU Factorisation
    luw = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(iterationMatrix);

    //
#ifdef DEBUG_NEWMARK
    std::cout.precision(15);
    std::cout << "Iteration matrix W: ";
    W->display();
#endif
  } else {
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrix(t,ds), this "
        "type of Dynamical System is not yet implemented");
  }
}

double siconos::integrators::NewMarkAlphaOSI::computeResidu() {
  DEBUG_BEGIN("siconos::integrators::NewMarkAlphaOSI::computeResidu() \n");
  // Compute the residual for each Dynamical system at step n and at iteration k
  // R_{n,k} = M_{n,k} ddotq_{n,k} - F_{n,k} - p_{n,k}
  // R_free = M_{n,k} ddotq_{n,k} - F_{n,k};
  // Compute norm of R_{n,k} for each DS
  // Take maximum norm of R_{n,k} over all DS
  double t = _simulation->nextTime();  // End of the time step
  // Iteration through the set of Dynamical Systems.
  //
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds;  // Current Dynamical System.
  double maxResidu = 0.0;
  double normResidu = 0.0;
  std::shared_ptr<siconos::algebra::SiconosVector> _residu;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto& residuFree = *workVectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE];
      // -- Convert the DS into a Lagrangian one.
      // Compute free residual
      residuFree.setZero();
      if (d->hasMass())
        residuFree = d->mass() * d->acceleration_read();
      else
        residuFree = d->acceleration_read();
      // LinearTIDS case
      if (auto lltids =
              std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
        // f = M*a + C*v + K*q
        // freeR = M*a + C*v

        if (lltids->hasExternalForces()) {
          lltids->computeFext(t);
          residuFree -= lltids->fext();  // freeR = _workspace[freeresidu] - F
        }
        if (lltids->hasDampingMatrix())
          residuFree += lltids->dampingMatrix() * lltids->velocity_read();
        if (lltids->hasStiffnessMatrix())
          residuFree += lltids->hasStiffnessMatrix() * lltids->q_read();

      } else  // LagrangianDS
      {
        // Update mass matrix
        d->computeMass(d->q_read());
        if (d->hasTotalForces()) {
          // Compute F = F_ext - F_int - F_Gyr
          d->computeTotalForces(d->velocity_read(), d->q_read(), t);
          residuFree -= d->totalForces();  // freeR = _workspace[freeresidu] - F
        }
      }
      // Compute residual
      _residu =
          std::make_shared<siconos::algebra::SiconosVector>(residuFree);  // _residu = freeR
      *_residu -= d->p_read(2);  // _residu = _workspace[freeresidu] - p
      // Compute Euclidean norm of the residual
      normResidu = _residu->norm2();
      // Take maximum value of norm over all DS
      if (normResidu > maxResidu) {
        maxResidu = normResidu;
      }
      //

      DEBUG_EXPR(residuFree.display(););
      DEBUG_EXPR(_residu->display(););
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::computeResidu(t,ds), this type of "
          "Dynamical System is not yet implemented");
    }
  }
  DEBUG_END("siconos::integrators::NewMarkAlphaOSI::computeResidu() \n");
  return maxResidu;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::computeFreeState() {
  DEBUG_BEGIN("siconos::integrators::NewMarkAlphaOSI::computeFreeState()\n");
  // Compute delta q_free = - ((W_{n,k})^-1)*R_free
  // Loop through the set of DS
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    auto& residuFree = *workVectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE];
    auto& qfree = *workVectors[siconos::integrators::NewMarkAlphaOSI::FREE];

    // -- Convert the DS into a Lagrangian one.
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      qfree = residuFree;
      //  assume that iterationMatrix and LU are uptodate
      qfree =
          _dynamicalSystemsGraph->properties(*dsi).LUW->solve(qfree);  //_qfree = (W^-1)*R_free
      qfree *= -1.0;  //_qfree = -(W^-1)*R_free
      //
      DEBUG_EXPR(qfree.display(););
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::computeResidu(t,ds), this type of "
          "Dynamical System is not yet implemented");
    }
  }
  DEBUG_END("siconos::integrators::NewMarkAlphaOSI::computeFreeState()\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  DEBUG_BEGIN(
      "siconos::integrators::NewMarkAlphaOSI::computeFreeOutput(siconos::graphs::"
      "InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)\n");
  auto t = _simulation->nextTime();
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);

  auto& workBlockV = *indexSet->properties(vertex_inter).workBlockVectors;
  // Get the type of relation
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();
  // Get the set of OSNSPs
  auto allOSNS = _simulation->oneStepNSProblems();
  // get the size of the interaction
  auto sizeY = inter->nonSmoothLaw()->size();
  // get pointer to delta q_free of Dynamical Systems concerned with the interaction

  std::shared_ptr<siconos::algebra::BlockVector> q_free;
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    q_free = workBlockV[siconos::integrators::NewMarkAlphaOSI::xfree];
  }
  assert(q_free);
  DEBUG_EXPR(q_free->display(););

  // get pointer to yForNSsolver vector

  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::NewMarkAlphaOSI::OSNSP_RHS];

  assert(q_free &&
         "In siconos::integrators::NewMarkAlphaOSI::computeFreeOutput: pointer q_free has not "
         "initialized yet");
  assert(inter->relation() &&
         "In siconos::integrators::NewMarkAlphaOSI::computeFreeOutput, relation associated "
         "with the interaction does not exist.");
  auto C = inter->relation()->C();
  assert(C &&
         "In siconos::integrators::NewMarkAlphaOSI::computeFreeOutput: Jacobian matrix does "
         "not exist");
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    if (relationSubType == siconos::modeling::RelationSubType::RheonomousR) {
      THROW_EXCEPTION(
          "siconos::integrators::NewMarkAlphaOSI::computeFreeOutput  not yet implemented "
          "with "
          "LagrangianRheonomousR");
    }
    DEBUG_EXPR(
        std::cout << "((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS]).get()"
                  << ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS]).get()
                  << std::endl;
        std::cout << "((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get()"
                  << ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get()
                  << std::endl;
        std::cout << "osnsp" << osnsp << std::endl;);
    if (relationSubType == siconos::modeling::RelationSubType::ScleronomousR) {
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get() ==
          osnsp)  // LCP at acceleration level
      {
        siconos::algebra::matrixBlockVector_prod(*C, *q_free, osnsp_rhs, true);
        auto ID = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
        ID->setIdentity();
        auto _SclerR = std::static_pointer_cast<siconos::modeling::LagrangianScleronomousR>(
            inter->relation());
        _SclerR->computedotjacqhXqdot(t, *inter);
        osnsp_rhs += *ID * *(_SclerR->dotjacqhXqdot());
        // y += NonLinearPart
      } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS]).get() ==
                 osnsp)  // LCP at position level
      {
        // Update Jacobian matrix
        inter->relation()->computeJach(t, *inter);
        // compute osnsp_rhs = y_{n,k} + G*q_free
        if (!_IsVelocityLevel)  // output at the position level y_{n,k} = g_{n,k}
        {
          inter->computeOutput(t, 0);  // Update output of level 0
          osnsp_rhs = *(inter->y(0));  // g_{n,k}
        }
        siconos::algebra::matrixBlockVector_prod(*C, *q_free, osnsp_rhs, false);
      } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]).get() ==
                 osnsp)  // output at the velocity level y_{n,k} = (h/gamma_prime)*dotg_{n,k}
      {
        auto h = _simulation->nextTime() - _simulation->startingTime();
        auto gamma_prime = _gamma / _beta;
        inter->computeOutput(t, 1);                        // Update output of level 1
        osnsp_rhs = (h / gamma_prime) * (*(inter->y(1)));  //(h/gamma_prime)*dotg_{n,k}
        siconos::algebra::matrixBlockVector_prod(*C, *q_free, osnsp_rhs, false);
      } else {
        osnsp->display();
        THROW_EXCEPTION(
            "siconos::integrators::NewMarkAlphaOSI::computeFreeOutput, this OSNSP does not "
            "exist");
      }
    }  // endif(relationSubType == ScleronomousR)

    DEBUG_EXPR(osnsp_rhs.display(););
  } else {
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::computeFreeOutput, this type of relation "
        "is not yet implemented");
  }
  DEBUG_END(
      "siconos::integrators::NewMarkAlphaOSI::computeFreeOutput(siconos::graphs::"
      "InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)\n");
}

void siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  DEBUG_BEGIN(
      "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n")

  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);

  // Get work matrices from the graph
  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
  auto& workMatrices = *_dynamicalSystemsGraph->properties(dsv).workMatrices;
  // Initialize memory buffers
  // _dynamicalSystemsGraph->bundle(dsv)->initMemory(getSizeMem());
  // // Force dynamical system to its initial state
  // _dynamicalSystemsGraph->bundle(dsv)->resetToInitialState();
  // Check dynamical system type
  auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
  assert(lds);
  // Compute W (iteration matrix)
  initializeIterationMatrix(ds);
  // allocate memory for work space for Newton iteration procedure
  assert(_dynamicalSystemsGraph->properties(dsv).iterationMatrix && "W is nullptr");
  // Allocate the memory to stock the acceleration-like variable
  if (lds) {
    lds->initRhs(t);

    ds_work_vectors.resize(siconos::integrators::NewMarkAlphaOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    // ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::FREE] =
    // std::make_shared<siconos::algebra::SiconosVector>(d->dimension()));
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->acceleration_read());
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_LIKE] =
        std::make_shared<siconos::algebra::SiconosVector>(
            lds->acceleration_read());  // set a0 = ddotq0
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_MEMORY] =
        std::make_shared<siconos::algebra::SiconosVector>(
            lds->acceleration_read());  // set a0 = ddotq0

    // Allocate the memory to stock coefficients of the polynomial for the dense output
    workMatrices.resize(siconos::integrators::NewMarkAlphaOSI::MAT_WORK_LENGTH);
    workMatrices[siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS] =
        std::make_shared<siconos::algebra::SiconosMatrix>(ds->dimension(),
                                                          (getOrderDenseOutput() + 1));

    //*(lds->workspace(DynamicalSystem::acce_like)) = *(lds->acceleration());
    // ds->allocateWorkVector(DynamicalSystem::acce_like, ds->dimension()); // allocate
    // memory for the acceleration-like of DS
    // ds->allocateWorkVector(DynamicalSystem::acce_memory, ds->dimension()); // allocate
    // memory to stock acceleration

    //          lds->allocateWorkMatrix(siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS,
    //          ds->dimension(), (osi_NewMark->getOrderDenseOutput() + 1));
    DEBUG_EXPR(lds->display());

    lds->swapInMemory();
  } else {
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::initialize: this type of DS is not yet "
        "implemented");
  }

  DEBUG_END(
      "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n")
}
void siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  DEBUG_BEGIN(
      "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction(...)\n")
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  auto& DSlink = inter.linkToDSVariables();

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::NewMarkAlphaOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::NewMarkAlphaOSI::BLOCK_WORK_LENGTH);
  }

  if (!interProp.workMatrices) {
    interProp.workMatrices =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>>>(
            siconos::integrators::NewMarkAlphaOSI::MAT_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();
  auto relationType = relation.getType();

  inter_work[siconos::integrators::NewMarkAlphaOSI::OSNSP_RHS] =
      std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  auto& nslaw = *inter.nonSmoothLaw();
  auto nslType = siconos::types::type_value(nslaw);

  if (nslType == siconos::modeling::Type::NewtonImpactNSL ||
      nslType == siconos::modeling::Type::MultipleImpactNSL) {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 2;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
  } else if (nslType == siconos::modeling::Type::NewtonImpactFrictionNSL) {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 4;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction  not "
        "yet "
        "implemented for nonsmooth law of type NewtonImpactFrictionNSL");
  } else
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction not yet "
        "implemented  for nonsmooth of type");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current
  // osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2)))) {
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction. The "
        "implementation is not correct for two different OSI for one interaction");
  }
  auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1);
    inter_work_block[siconos::integrators::NewMarkAlphaOSI::xfree] =
        std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[siconos::integrators::NewMarkAlphaOSI::xfree]->insertPtr(
        workVds1[siconos::integrators::NewMarkAlphaOSI::FREE]);
    DSlink[siconos::modeling::LagrangianR::p2] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::p2]->insertPtr(lds.p(2));
    DSlink[siconos::modeling::LagrangianR::q2] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::q2]->insertPtr(lds.acceleration());
  }
  // else if (relationType == NewtonEuler)
  // {
  //   inter_work_block[siconos::integrators::NewMarkAlphaOSI::xfree] =
  //   std::make_shared<siconos::algebra::BlockVector>();
  //   inter_work_block[siconos::integrators::NewMarkAlphaOSI::xfree]->insertPtr(workVds1[siconos::integrators::NewMarkAlphaOSI::FREE]);
  // }

  if (ds1 != ds2) {
    auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2);
      inter_work_block[siconos::integrators::NewMarkAlphaOSI::xfree]->insertPtr(
          workVds2[siconos::integrators::NewMarkAlphaOSI::FREE]);
      DSlink[siconos::modeling::LagrangianR::p2]->insertPtr(lds.p(2));
      DSlink[siconos::modeling::LagrangianR::q2]->insertPtr(lds.acceleration());
    }
    // else if (relationType == NewtonEuler)
    // {
    //   inter_work_block[siconos::integrators::NewMarkAlphaOSI::xfree]->insertPtr(workVds2[siconos::integrators::NewMarkAlphaOSI::FREE]);
    // }
  }

  DEBUG_END(
      "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction(...)\n")
}

void siconos::integrators::NewMarkAlphaOSI::prepareNewtonIteration(double time) {
  DEBUG_BEGIN("siconos::integrators::NewMarkAlphaOSI::prepareNewtonIteration(double time)\n");
  // Compute matrix W for all Dynamical Systems
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    computeIterationMatrix(_dynamicalSystemsGraph->bundle(*dsi),
                           *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix,
                           _dynamicalSystemsGraph->properties(*dsi).LUW);
  }
  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);
  }
  DEBUG_END("siconos::integrators::NewMarkAlphaOSI::prepareNewtonIteration(double time)\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::prediction() {
  DEBUG_BEGIN("siconos::integrators::NewMarkAlphaOSI::prediction()\n");
  // Step size
  double h = _simulation->nextTime() - _simulation->startingTime();
  if (h < 100 * siconos::internal::MACHINE_PREC)
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::prediction, time integration is too "
        "small");
  // Loop over all DS
  std::shared_ptr<siconos::algebra::SiconosVector> _q, _dotq, _ddotq;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      DEBUG_EXPR(d->display(););
      siconos::algebra::SiconosVector& acce_like =
          *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_LIKE];
      // Save the acceleration before the prediction
      acce_like = d->acceleration_read();

      DEBUG_PRINT("Before prediction :\n")
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););
      DEBUG_EXPR(d->acceleration()->display(););
      DEBUG_EXPR(acce_like.display(););

      *d->q() +=
          d->velocity_read() * h +
          acce_like * (h * h * (0.5 - _beta));  // q_{n+1} = q_n + (h^2)*(0.5 - beta)*a_n
      *d->velocity() +=
          acce_like * (h * (1 - _gamma));  // dotq_{n+1} = dotq_n + h*(1 - gamma)*a_n
      acce_like = (_alpha_f / (1 - _alpha_m)) * d->acceleration_read() -
                  (_alpha_m / (1 - _alpha_m)) *
                      acce_like;  // a_{n+1} = (alpha_f*ddotq_n - alpha_m*a_n)/(1 - alpha_m)
      *d->q() += (h * h * _beta) * acce_like;
      *d->velocity() += (h * _gamma) * acce_like;
      d->acceleration()->setZero();

      DEBUG_PRINT("After prediction :\n")
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););
      DEBUG_EXPR(d->acceleration()->display(););
      DEBUG_EXPR(acce_like.display(););
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::prediction: this type of DS is not yet "
          "implemented");
    }
  }
  DEBUG_END("siconos::integrators::NewMarkAlphaOSI::prediction()\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::correction() {
  auto h = _simulation->nextTime() - _simulation->startingTime();
  auto beta_prime = (1. - _alpha_m) / ((1. - _alpha_f) * _beta);
  auto gamma_prime = _gamma / _beta;
  // Make sure that the input of the concerned Dynamical Systems is updated after solving LCP
  std::shared_ptr<siconos::algebra::SiconosVector> delta_q;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    auto W = _dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    // Its W matrix of iteration.
    // Iteration matrix W_{n+1,k} computed at kth iteration
    auto& residuFree = *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE];
    auto& acce_like = *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_LIKE];

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // resultant force p_{n+1,k+1} of DS at (k+1)th iteration
      std::shared_ptr<siconos::algebra::SiconosVector> _p = d->p(2);
      // Compute delta_q = W_{n+1,k}^{-1}(p_{n+1,k+1} - r_{n+1,k})
      siconos::algebra::SiconosVector delta_q{*_p - residuFree};
      // copy (p_{n+1,k+1} - r_{n+1,k}) to delta_q

      delta_q = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(delta_q);

      // Correction q_{n+1,k+1}, dotq_{n+1,k+1}, ddotq_{n+1,k+1}
      *(d->q()) += delta_q;  // q_{n+1,k+1} = q_{n+1,k} + delta_q
      *(d->velocity()) += (gamma_prime / h) *
                          delta_q;  // dotq_{n+1,k+1} = dotq_{n+1,k} + (gamma_prime/h)*delta_q
      *(d->acceleration()) +=
          beta_prime / (h * h) *
          delta_q;  // ddotq_{n+1,k+1} = ddotq_{n+1,k} + (beta_prime/h^2)*delta_q
      // a_{n+1,k+1} = a_{n+1,k} + ((1-alpha_f)/(1-alpha_m))*(beta_prime/h^2)*delta_q

      acce_like += ((1 - _alpha_f) / (1 - _alpha_m)) * ((beta_prime / (h * h)) * delta_q);

      DEBUG_PRINT("After correction : \n");
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););
      DEBUG_EXPR(d->acceleration()->display(););
      DEBUG_EXPR(acce_like.display(););
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::updateState: this type of DS is not "
          "yet "
          "implemented");
    }
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::integrate(double& t_ini, double& t_end,
                                                      double& t_out, int& flag) {
  THROW_EXCEPTION(
      "In siconos::integrators::NewMarkAlphaOSI::integrate, this method does nothing in the "
      "NewMarkAlpha OneStepIntegrator!!!");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::updateState(const unsigned int level) {
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  if (level == 1)  // ie impact case: compute velocity
  {
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto lds = std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds);
      lds->computePostImpactVelocity();
    }
  } else if (level == 2) {
    auto time = _simulation->nextTime();
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      ds->update(time);
    }
  } else
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::updateState, index is out of range. "
        "Index "
        "= " +
        std::to_string(level));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::computeCoefsDenseOutput(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  auto h = _simulation->nextTime() - _simulation->startingTime();

  auto _vec = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  auto& workMatrices =
      *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).workMatrices;
  auto& ds_work_vectors =
      *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).workVectors;
  if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    const siconos::algebra::SiconosVector &q_n(d->qMemory().getSiconosVector(0)),  // q_n
        dotq_n(d->velocityMemory().getSiconosVector(0)),                           // dotq_n
        // ddotq_n = d->workspace(DynamicalSystem::acce_memory); // ddotq_n
        q_np1(d->q_read()),                 // q_{n+1}
        dotq_np1(d->velocity_read()),       // dotq_{n+1}
        ddotq_np1(d->acceleration_read());  // ddotq_{n+1}

    auto& ddotq_n = *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_MEMORY];

    auto _CoeffsDense =
        workMatrices[siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS];
    // d->workMatrix(siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS); //
    // matrix of coefficients [a0 a1 a2 a3 a4 a5]
    if (_CoeffsDense->size(1) != 6) {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::computeCoefsDenseOutput: the number of "
          "polynomial coeffcients considered here must equal to 6 (dense output polynomial "
          "of "
          "order 5)");
    }
    // a0 = q_n
    (*_vec) = q_n;
    _CoeffsDense->col(0) = (*_vec);
    DEBUG_EXPR(std::cout << "a0: "; _vec->display(););
    // a1 = h*dotq_n
    (*_vec) = h * dotq_n;
    _CoeffsDense->col(1) = (*_vec);
    DEBUG_EXPR(std::cout << "a1: "; _vec->display(););
    // a2 = 0.5*h^2*ddotq_n
    (*_vec) = (0.5 * h * h) * ddotq_n;
    _CoeffsDense->col(2) = (*_vec);
    DEBUG_EXPR(std::cout << "a2: "; _vec->display(););
    // a3 = -10*q_n - 6*h*dotq_n - 1.5*h^2*ddotq_n + 10*q_{n+1} - 4*h*dotq_{n+1} +
    // 0.5*h^2*ddotq_{n+1}
    (*_vec) = (-10.0) * q_n - (6.0 * h) * dotq_n - (1.5 * h * h) * ddotq_n + 10.0 * q_np1 -
              (4.0 * h) * dotq_np1 + (0.5 * h * h) * ddotq_np1;
    _CoeffsDense->col(3) = (*_vec);
    DEBUG_EXPR(std::cout << "a3: "; _vec->display(););
    // a4 = 15*q_n + 8*h*dotq_n + 1.5*h^2*ddotq_n - 15*q_{n+1} + 7*h*dotq_{n+1} -
    // h^2*ddotq_{n+1}
    (*_vec) = 15.0 * q_n + (8.0 * h) * dotq_n + (1.5 * h * h) * ddotq_n - 15.0 * q_np1 +
              (7.0 * h) * dotq_np1 - h * h * ddotq_np1;
    _CoeffsDense->col(4) = (*_vec);
    DEBUG_EXPR(std::cout << "a4: "; _vec->display(););
    // a5 = -6*q_n - 3*h*dotq_n - 0.5*h^2*ddotq_n + 6*q_{n+1} - 3*h*dotq_{n+1} +
    // 0.5*h^2*ddotq_{n+1}
    (*_vec) = (-6.0) * q_n - (3.0 * h) * dotq_n - (0.5 * h * h) * ddotq_n + 6.0 * q_np1 -
              (3.0 * h) * dotq_np1 + (0.5 * h * h) * ddotq_np1;
    _CoeffsDense->col(5) = (*_vec);
    DEBUG_EXPR(std::cout << "a5: "; _vec->display(););
    //
#ifdef DEBUG_NEWMARK
    std::cout
        << "==================== In "
           "siconos::integrators::NewMarkAlphaOSI::computeCoefsDenseOutput ================"
        << std::endl;
    std::cout << "DS number: " << ds->number() << std::endl;
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
    std::cout << "Dense output coefficient matrix: " << std::endl;
    _CoeffsDense->display();
#endif
  } else {
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::computeCoefsDenseOutput: this type of DS "
        "has not been implemented yet");
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::prepareEventLocalization() {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    // Compute coefficients of the dense output polynomial for all Dynamical Systems
    computeCoefsDenseOutput(ds);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::DenseOutputallDSs(double t) {
  // Make sure that all coefficients of the dense output polynomial for all DSs has been
  // computed before
  auto t_n = _simulation->startingTime();
  auto t_np1 = _simulation->nextTime();
  auto h = t_np1 - t_n;
  auto theta = (t - t_n) / h;
  siconos::algebra::SiconosVector vec1{6};
  vec1 << 1., theta, std::pow(theta, 2), std::pow(theta, 3), std::pow(theta, 4),
      std::pow(theta, 5);

  siconos::algebra::SiconosVector vec2{6};
  vec2 << 0., 1. / h, (2.0 * theta) / h, (3.0 * std::pow(theta, 2)) / h,
      (4.0 * std::pow(theta, 3)) / h, (5.0 * std::pow(theta, 4)) / h;

  siconos::algebra::SiconosVector vec3{6};
  vec3 << 0., 0., 2.0 / std::pow(h, 2), (6.0 * theta) / std::pow(h, 2),
      (12.0 * std::pow(theta, 2)) / std::pow(h, 2),
      (20.0 * std::pow(theta, 3)) / std::pow(h, 2);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workMatrices = *_dynamicalSystemsGraph->properties(*dsi).workMatrices;

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto& coeffsDense =
          *workMatrices[siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS];
      *d->q() = coeffsDense * vec1;
      *d->velocity() = coeffsDense * vec2;      // vel = Matrix_coeffs*_vec2
      *d->acceleration() = coeffsDense * vec3;  // ddotq = Matrix_coeffs*_vec3
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::DenseOutputallDSs: this type of DS has "
          "not been implemented yet");
    }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::display() const {}
