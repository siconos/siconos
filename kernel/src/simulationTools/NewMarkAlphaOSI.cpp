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
#include "SiconosMatrixOp.hpp"        // scal
#include "SiconosMatrixVectorOp.hpp"  // mat-vec prod
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
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

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::integrators::NewMarkAlphaOSI::W(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds && "siconos::integrators::NewMarkAlphaOSI::W(ds): ds == nullptr.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "siconos::integrators::NewMarkAlphaOSI::W(ds): W[ds] == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
}

void siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrixW(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrixW(t,ds) - ds == "
        "nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrixW(t,ds) - ds does "
        "not belong to the OSI.");

  if (_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W)
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrixW(t,ds) - W(ds) is "
        "already in the map and has been initialized.");

  _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W =
      std::make_shared<siconos::algebra::SimpleMatrix>(ds->dimension(), ds->dimension());

  computeW(ds, *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W);
}

void siconos::integrators::NewMarkAlphaOSI::computeW(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    siconos::algebra::SiconosMatrix& W) {
  auto beta_prime = (1. - _alpha_m) / ((1. - _alpha_f) * _beta);
  auto gamma_prime = _gamma / _beta;
  auto h = _simulation->nextTime() - _simulation->startingTime();  // step size
  if (h < 100 * siconos::internal::MACHINE_PREC)
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrixW(t,ds), time "
        "integration is too small");
  // make sure that W is initialized before computing
  std::shared_ptr<siconos::algebra::SiconosMatrix> M;
  std::shared_ptr<siconos::algebra::SiconosMatrix> K;
  std::shared_ptr<siconos::algebra::SiconosMatrix> C;
  if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    if (auto ltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
      K = ltids->K();     // matrix K
      if (K) *K *= -1.0;  // K = -K
      C = ltids->C();     // matrix C
      if (C) *C *= -1.0;  // C = -C
    } else                // LagrangianLinearDS
    {
      K = lds->jacobianqForces();  // jacobian according to q
      C = lds->jacobianvForces();  // jacobian according to velocity
      lds->computeMass(lds->q());
    }

    // Compute W = (beta_prime/h^2)*M - (gamma_prime/h)*C - K
    if (lds->mass())
      siconos::algebra::scal(beta_prime / (h * h), *(lds->mass()), W, true);
    else {
      W.eye();
      W *= beta_prime / (h * h);
    }
    if (C) siconos::algebra::scal(-gamma_prime / h, *C, W, false);
    if (K) siconos::algebra::scal(-1.0, *K, W, false);
      //
#ifdef DEBUG_NEWMARK
    std::cout.precision(15);
    std::cout << "Iteration matrix W: ";
    W->display();
#endif
  } else {
    THROW_EXCEPTION(
        "In siconos::integrators::NewMarkAlphaOSI::initializeIterationMatrixW(t,ds), this "
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
      // get position, velocity and acceleration
      auto q = d->q();
      auto v = d->velocity();
      auto a = d->acceleration();
      std::shared_ptr<siconos::algebra::SiconosVector> F;
      // get the reaction force p
      auto p = d->p(2);
      // Compute free residual
      residuFree.zero();
      if (d->mass())
        siconos::algebra::prod(*d->mass(), *a, residuFree, true);  // freeR = M*a
      else
        residuFree = *a;
      // LinearTIDS case
      if (auto ltids =
              std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
        // We need to add F_int = Cv + Kq to freeR
        auto K = ltids->K();
        auto C = ltids->C();
        F = ltids->fExt();  // Note that for LagrangianLinearTIDS, F = F_ext
        if (F) {
          ltids->computeFExt(t);
        }
        if (K) siconos::algebra::prod(*K, *q, residuFree, false);  // f = M*a + C*v + K*q
        if (C) siconos::algebra::prod(*C, *v, residuFree, false);  // freeR = M*a + C*v
      } else                                                       // LagrangianDS
      {
        // Update mass matrix
        d->computeMass();
        F = d->forces();
        if (F)
          // Compute F = F_ext - F_int - F_Gyr
          d->computeForces(t, d->q(), d->velocity());
      }

      residuFree -= *F;  // freeR = _workspace[freeresidu] - F
      // Compute residual
      _residu =
          std::make_shared<siconos::algebra::SiconosVector>(residuFree);  // _residu = freeR
      *_residu -= *p;  // _residu = _workspace[freeresidu] - p
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

    // Get iteration matrix W, make sure that W was updated before
    auto W = _dynamicalSystemsGraph->properties(*dsi).W;  // Its W matrix of iteration.

    auto& residuFree = *workVectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE];
    auto& qfree = *workVectors[siconos::integrators::NewMarkAlphaOSI::FREE];

    // -- Convert the DS into a Lagrangian one.
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      qfree = residuFree;
      W->Solve(qfree);  //_qfree = (W^-1)*R_free
      qfree *= -1.0;    //_qfree = -(W^-1)*R_free
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
          "siconos::integrators::NewMarkAlphaOSI::computeFreeOutput  not yet implemented with "
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
      std::vector<std::size_t> coord(8);
      coord[0] = 0;
      coord[1] = sizeY;
      coord[2] = 0;
      coord[3] = C->size(1);
      coord[4] = 0;
      coord[5] = C->size(1);
      coord[6] = 0;
      coord[7] = sizeY;
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get() ==
          osnsp)  // LCP at acceleration level
      {
        siconos::algebra::subprod(*C, *q_free, osnsp_rhs, coord, true);
        auto ID = std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeY);
        ID->eye();
        std::vector<std::size_t> xcoord(8);
        xcoord[0] = 0;
        xcoord[1] = sizeY;
        xcoord[2] = 0;
        xcoord[3] = sizeY;
        xcoord[4] = 0;
        xcoord[5] = sizeY;
        xcoord[6] = 0;
        xcoord[7] = sizeY;
        auto _SclerR = std::static_pointer_cast<siconos::modeling::LagrangianScleronomousR>(
            inter->relation());
        _SclerR->computedotjacqhXqdot(t, *inter);
        siconos::algebra::subprod(*ID, *(_SclerR->dotjacqhXqdot()), osnsp_rhs, xcoord,
                                  false);  // y += NonLinearPart
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
        siconos::algebra::subprod(*C, *q_free, osnsp_rhs, coord, false);
      } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]).get() ==
                 osnsp)  // output at the velocity level y_{n,k} = (h/gamma_prime)*dotg_{n,k}
      {
        auto h = _simulation->nextTime() - _simulation->startingTime();
        auto gamma_prime = _gamma / _beta;
        inter->computeOutput(t, 1);                        // Update output of level 1
        osnsp_rhs = (h / gamma_prime) * (*(inter->y(1)));  //(h/gamma_prime)*dotg_{n,k}
        siconos::algebra::subprod(*C, *q_free, osnsp_rhs, coord, false);
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
  initializeIterationMatrixW(ds);
  // allocate memory for work space for Newton iteration procedure
  assert(_dynamicalSystemsGraph->properties(dsv).W && "W is nullptr");
  // Allocate the memory to stock the acceleration-like variable
  if (lds) {
    lds->initRhs(t);

    ds_work_vectors.resize(siconos::integrators::NewMarkAlphaOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    // ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::FREE] =
    // std::make_shared<siconos::algebra::SiconosVector>(d->dimension()));
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(*(lds->acceleration()));
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_LIKE] =
        std::make_shared<siconos::algebra::SiconosVector>(
            *(lds->acceleration()));  // set a0 = ddotq0
    ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_MEMORY] =
        std::make_shared<siconos::algebra::SiconosVector>(
            *(lds->acceleration()));  // set a0 = ddotq0

    // Allocate the memory to stock coefficients of the polynomial for the dense output
    workMatrices.resize(siconos::integrators::NewMarkAlphaOSI::MAT_WORK_LENGTH);
    workMatrices[siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS] =
        std::make_shared<siconos::algebra::SimpleMatrix>(ds->dimension(),
                                                         (getOrderDenseOutput() + 1));

    //*(lds->workspace(DynamicalSystem::acce_like)) = *(lds->acceleration());
    // ds->allocateWorkVector(DynamicalSystem::acce_like, ds->dimension()); // allocate memory
    // for the acceleration-like of DS ds->allocateWorkVector(DynamicalSystem::acce_memory,
    // ds->dimension()); // allocate memory to stock acceleration

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
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SimpleMatrix>>>(
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
        "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction  not yet "
        "implemented for nonsmooth law of type NewtonImpactFrictionNSL");
  } else
    THROW_EXCEPTION(
        "siconos::integrators::NewMarkAlphaOSI::initializeWorkVectorsForInteraction not yet "
        "implemented  for nonsmooth of type");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
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
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& W = *_dynamicalSystemsGraph->properties(*dsi).W;
    computeW(ds, W);
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
        "In siconos::integrators::NewMarkAlphaOSI::prediction, time integration is too small");
  // Loop over all DS
  std::shared_ptr<siconos::algebra::SiconosVector> _q, _dotq, _ddotq;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      DEBUG_EXPR(d->display(););
      _q = d->q();                 // generalized coordinate
      _dotq = d->velocity();       // generalized velocity
      _ddotq = d->acceleration();  // generalized acceleration
      siconos::algebra::SiconosVector& acce_like =
          *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_LIKE];
      // Save the acceleration before the prediction
      acce_like = *(_ddotq);

      DEBUG_PRINT("Before prediction :\n")
      DEBUG_EXPR(_q->display(););
      DEBUG_EXPR(_dotq->display(););
      DEBUG_EXPR(_ddotq->display(););
      DEBUG_EXPR(acce_like.display(););

      *_q = *_q + (*_dotq) * h +
            acce_like * (h * h * (0.5 - _beta));  // q_{n+1} = q_n + (h^2)*(0.5 - beta)*a_n
      *_dotq =
          *_dotq + acce_like * (h * (1 - _gamma));  // dotq_{n+1} = dotq_n + h*(1 - gamma)*a_n
      acce_like = (_alpha_f / (1 - _alpha_m)) * (*_ddotq) -
                  (_alpha_m / (1 - _alpha_m)) *
                      acce_like;  // a_{n+1} = (alpha_f*ddotq_n - alpha_m*a_n)/(1 - alpha_m)
      *_q = *_q + (h * h * _beta) * acce_like;
      *_dotq = *_dotq + (h * _gamma) * acce_like;
      _ddotq->zero();

      DEBUG_PRINT("After prediction :\n")
      DEBUG_EXPR(_q->display(););
      DEBUG_EXPR(_dotq->display(););
      DEBUG_EXPR(_ddotq->display(););
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

    auto W = _dynamicalSystemsGraph->properties(*dsi).W;
    // Its W matrix of iteration.
    // Iteration matrix W_{n+1,k} computed at kth iteration
    auto& residuFree = *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::RESIDU_FREE];
    auto& acce_like = *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_LIKE];

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // resultant force p_{n+1,k+1} of DS at (k+1)th iteration
      std::shared_ptr<siconos::algebra::SiconosVector> _p = d->p(2);
      // Compute delta_q = W_{n+1,k}^{-1}(p_{n+1,k+1} - r_{n+1,k})
      delta_q = std::make_shared<siconos::algebra::SiconosVector>(
          *_p - residuFree);  // copy (p_{n+1,k+1} - r_{n+1,k}) to delta_q
      W->Solve(*delta_q);
      // Correction q_{n+1,k+1}, dotq_{n+1,k+1}, ddotq_{n+1,k+1}
      *(d->q()) += *delta_q;  // q_{n+1,k+1} = q_{n+1,k} + delta_q
      *(d->velocity()) +=
          (gamma_prime / h) *
          (*delta_q);  // dotq_{n+1,k+1} = dotq_{n+1,k} + (gamma_prime/h)*delta_q
      *(d->acceleration()) +=
          beta_prime / (h * h) *
          (*delta_q);  // ddotq_{n+1,k+1} = ddotq_{n+1,k} + (beta_prime/h^2)*delta_q
      // a_{n+1,k+1} = a_{n+1,k} + ((1-alpha_f)/(1-alpha_m))*(beta_prime/h^2)*delta_q

      acce_like += ((1 - _alpha_f) / (1 - _alpha_m)) * ((beta_prime / (h * h)) * (*delta_q));

      DEBUG_PRINT("After correction : \n");
      DEBUG_EXPR(d->q()->display(););
      DEBUG_EXPR(d->velocity()->display(););
      DEBUG_EXPR(d->acceleration()->display(););
      DEBUG_EXPR(acce_like.display(););
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::updateState: this type of DS is not yet "
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
        "In siconos::integrators::NewMarkAlphaOSI::updateState, index is out of range. Index "
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
        q_np1(*d->q()),                 // q_{n+1}
        dotq_np1(*d->velocity()),       // dotq_{n+1}
        ddotq_np1(*d->acceleration());  // ddotq_{n+1}

    auto& ddotq_n = *ds_work_vectors[siconos::integrators::NewMarkAlphaOSI::ACCE_MEMORY];

    auto _CoeffsDense =
        workMatrices[siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS];
    // d->workMatrix(siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS); //
    // matrix of coefficients [a0 a1 a2 a3 a4 a5]
    if (_CoeffsDense->size(1) != 6) {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::computeCoefsDenseOutput: the number of "
          "polynomial coeffcients considered here must equal to 6 (dense output polynomial of "
          "order 5)");
    }
    // a0 = q_n
    (*_vec) = q_n;
    _CoeffsDense->setCol(0, (*_vec));
    DEBUG_EXPR(std::cout << "a0: "; _vec->display(););
    // a1 = h*dotq_n
    (*_vec) = h * dotq_n;
    _CoeffsDense->setCol(1, (*_vec));
    DEBUG_EXPR(std::cout << "a1: "; _vec->display(););
    // a2 = 0.5*h^2*ddotq_n
    (*_vec) = (0.5 * h * h) * ddotq_n;
    _CoeffsDense->setCol(2, (*_vec));
    DEBUG_EXPR(std::cout << "a2: "; _vec->display(););
    // a3 = -10*q_n - 6*h*dotq_n - 1.5*h^2*ddotq_n + 10*q_{n+1} - 4*h*dotq_{n+1} +
    // 0.5*h^2*ddotq_{n+1}
    (*_vec) = (-10.0) * q_n - (6.0 * h) * dotq_n - (1.5 * h * h) * ddotq_n + 10.0 * q_np1 -
              (4.0 * h) * dotq_np1 + (0.5 * h * h) * ddotq_np1;
    _CoeffsDense->setCol(3, (*_vec));
    DEBUG_EXPR(std::cout << "a3: "; _vec->display(););
    // a4 = 15*q_n + 8*h*dotq_n + 1.5*h^2*ddotq_n - 15*q_{n+1} + 7*h*dotq_{n+1} -
    // h^2*ddotq_{n+1}
    (*_vec) = 15.0 * q_n + (8.0 * h) * dotq_n + (1.5 * h * h) * ddotq_n - 15.0 * q_np1 +
              (7.0 * h) * dotq_np1 - h * h * ddotq_np1;
    _CoeffsDense->setCol(4, (*_vec));
    DEBUG_EXPR(std::cout << "a4: "; _vec->display(););
    // a5 = -6*q_n - 3*h*dotq_n - 0.5*h^2*ddotq_n + 6*q_{n+1} - 3*h*dotq_{n+1} +
    // 0.5*h^2*ddotq_{n+1}
    (*_vec) = (-6.0) * q_n - (3.0 * h) * dotq_n - (0.5 * h * h) * ddotq_n + 6.0 * q_np1 -
              (3.0 * h) * dotq_np1 + (0.5 * h * h) * ddotq_np1;
    _CoeffsDense->setCol(5, (*_vec));
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
  auto _vec1 = std::make_shared<siconos::algebra::SiconosVector>(_orderDenseOutput + 1);
  assert((_vec1->size() == 6) && "There are six coefficients of the dense output polynomial");
  (*_vec1)(0) = 1.0;
  (*_vec1)(1) = theta;
  (*_vec1)(2) = std::pow(theta, 2);
  (*_vec1)(3) = std::pow(theta, 3);
  (*_vec1)(4) = std::pow(theta, 4);
  (*_vec1)(5) = std::pow(theta, 5);
  //
  auto _vec2 = std::make_shared<siconos::algebra::SiconosVector>(_orderDenseOutput + 1);
  assert((_vec2->size() == 6) && "There are six coefficients of the dense output polynomial");
  (*_vec2)(0) = 0.0;
  (*_vec2)(1) = 1.0 / h;
  (*_vec2)(2) = (2.0 * theta) / h;
  (*_vec2)(3) = (3.0 * std::pow(theta, 2)) / h;
  (*_vec2)(4) = (4.0 * std::pow(theta, 3)) / h;
  (*_vec2)(5) = (5.0 * std::pow(theta, 4)) / h;
  //
  auto _vec3 = std::make_shared<siconos::algebra::SiconosVector>(_orderDenseOutput + 1);
  assert((_vec3->size() == 6) && "There are six coefficients of the dense output polynomial");
  (*_vec3)(0) = 0.0;
  (*_vec3)(1) = 0.0;
  (*_vec3)(2) = 2.0 / std::pow(h, 2);
  (*_vec3)(3) = (6.0 * theta) / std::pow(h, 2);
  (*_vec3)(4) = (12.0 * std::pow(theta, 2)) / std::pow(h, 2);
  (*_vec3)(5) = (20.0 * std::pow(theta, 3)) / std::pow(h, 2);
  //

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workMatrices = *_dynamicalSystemsGraph->properties(*dsi).workMatrices;

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto& coeffsDense =
          *workMatrices[siconos::integrators::NewMarkAlphaOSI::DENSE_OUTPUT_COEFFICIENTS];
      siconos::algebra::prod(coeffsDense, *_vec1, *(d->q()), true);  // q = Matrix_coeffs*_vec1
      siconos::algebra::prod(coeffsDense, *_vec2, *(d->velocity()),
                             true);  // dotq = Matrix_coeffs*_vec2
      siconos::algebra::prod(coeffsDense, *_vec3, *(d->acceleration()),
                             true);  // ddotq = Matrix_coeffs*_vec3
    } else {
      THROW_EXCEPTION(
          "In siconos::integrators::NewMarkAlphaOSI::DenseOutputallDSs: this type of DS has "
          "not been implemented yet");
    }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::integrators::NewMarkAlphaOSI::display() const {}
