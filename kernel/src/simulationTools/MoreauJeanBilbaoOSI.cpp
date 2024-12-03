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
#include "MoreauJeanBilbaoOSI.hpp"

#include "BlockVector.hpp"
#include "BoundaryCondition.hpp"
#include "Interaction.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonImpactNSL.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for subprod
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
#include "Simulation.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
#include "siconos_debug.h"

// --- constructor from a set of data ---
siconos::integrators::MoreauJeanBilbaoOSI::MoreauJeanBilbaoOSI()
    : OneStepIntegrator(IntegratorType::MOREAUJEANBILBAOOSI, 1, 0, 1, 1, 1) {}

void siconos::integrators::MoreauJeanBilbaoOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // Get work buffers from the graph and initialize ds state and memory
  auto& work_ds = *_initializeDSWorkVectors(ds);

  DEBUG_PRINT("initializeWorkVectorsForDS() \n");
  // Check consistency between OSI type and DS type
  if (auto lldds =
          std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearDiagonalDS>(ds)) {
    // Prepare work buffers
    work_ds.resize(siconos::integrators::MoreauJeanBilbaoOSI::WORK_LENGTH);
    // - work buffer, used to save vfree.
    work_ds[siconos::integrators::MoreauJeanBilbaoOSI::VFREE] =
        std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
    // Other buffers required to compute iteration matrix
    work_ds[siconos::integrators::MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR] =
        std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
    work_ds[siconos::integrators::MoreauJeanBilbaoOSI::ONE_MINUS_THETA] =
        std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());

    // Initialize the iteration matrix and other parameters of the osi (all
    // ds-dependent)
    _initialize_iteration_matrix(ds);
    // Update dynamical system components (for memory swap).
    lldds->computeTotalForces(lldds->velocity_read(), lldds->q_read(), t);
    lldds->swapInMemory();
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanBilbaoOSI::initializeWorkVectorsForDS "
        "- only "
        "implemented for LagrangianLinearDiagonalDS");
}

void siconos::integrators::MoreauJeanBilbaoOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  // Get the dynamical systems linked by inter
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::MoreauJeanBilbaoOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::MoreauJeanBilbaoOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();
  auto relationType = relation.getType();
  auto relationSubType = inter.relation()->getSubType();

  inter_work[siconos::integrators::MoreauJeanBilbaoOSI::OSNSP_RHS] =
      std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  if (relationType != siconos::modeling::RelationType::Lagrangian ||
      relationSubType != siconos::modeling::RelationSubType::LinearTIR)
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanBilbaoOSI::computeFreeOutput only "
        "Lagrangian Linear "
        "Relations are allowed.");

  // // Check if interations levels (i.e. y and lambda sizes) are compliant with
  // the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  // allocate and set work vectors for the osi

  if (checkOSI(DSG.descriptor(ds1))) {
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    if (!inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]) {
      inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree] =
          std::make_shared<siconos::algebra::BlockVector>();
      inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]->insertPtr(
          workVds1[siconos::integrators::MoreauJeanBilbaoOSI::VFREE]);
    } else {
      inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]->setVectorPtr(
          0, workVds1[siconos::integrators::MoreauJeanBilbaoOSI::VFREE]);
    }
  }
  if (ds1 != ds2) {
    if (checkOSI(DSG.descriptor(ds2))) {
      auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      if (!inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]) {
        inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree] =
            std::make_shared<siconos::algebra::BlockVector>();
        // dummy insertion to reserve first vector for ds1
        inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]->insertPtr(
            workVds2[siconos::integrators::MoreauJeanBilbaoOSI::VFREE]);
        inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]->insertPtr(
            workVds2[siconos::integrators::MoreauJeanBilbaoOSI::VFREE]);
      } else {
        inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree]->insertPtr(
            workVds2[siconos::integrators::MoreauJeanBilbaoOSI::VFREE]);
      }
    }
  }
}

void siconos::integrators::MoreauJeanBilbaoOSI::computeInitialNewtonState() {
  DEBUG_BEGIN("MoreauJeanOSI::computeInitialNewtonState()\n");
  // Compute the position value giving the initial velocity.
  // The goal is to save one newton iteration for nearly linear system
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto& ds = *_dynamicalSystemsGraph->bundle(*dsi);

    // The goal is to converge in one iteration if the system is almost linear
    // we start the Newton loop q = q0+hv0
    updatePosition(ds);
  }
  DEBUG_END("MoreauJeanBilbaoOSI::computeInitialNewtonState()\n");
}

void siconos::integrators::MoreauJeanBilbaoOSI::initialize_nonsmooth_problems() {
  auto allOSNS = _simulation->oneStepNSProblems();
  ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY])->setIndexSetLevel(1);
  ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY])->setInputOutputLevel(1);
}

void siconos::integrators::MoreauJeanBilbaoOSI::_initialize_iteration_matrix(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // This function:
  // - allocates memory for the matrix W
  // - updates its content for the current (initial) state of the dynamical
  // system, depending on its type.
  //
  // W = mass + time_step**2/2(I - Theta)Stiffness + time_step * Sigma*
  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
  if (_dynamicalSystemsGraph->properties(dsv).iterationMatrix)
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanBilbaoOSI::_initialize_iteration_"
        "matrix(ds) - W has "
        "already been initialized by another osi");

  auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  if (auto lldds =
          std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearDiagonalDS>(ds)) {
    auto ndof = lldds->dimension();
    // Allocate work buffers for:
    // - Iteration matrix
    _dynamicalSystemsGraph->properties(dsv).iterationMatrix =
        std::make_shared<siconos::algebra::SiconosMatrix>(
            ndof, ndof);  // WARNING : Use bandmatrix instead ?

    // - I - theta
    // ds_work_vectors[siconos::integrators::MoreauJeanBilbaoOSI::ONE_MINUS_THETA]
    // = std::make_shared<siconos::algebra::SiconosVector>(ndof));
    // - dt * sigma*
    // ds_work_vectors[siconos::integrators::MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR]
    // = std::make_shared<siconos::algebra::SiconosVector>(ndof));
    auto& iteration_matrix = *_dynamicalSystemsGraph->properties(dsv).iterationMatrix;
    auto omega2 = lldds->stiffnessMatrix();
    auto damp = lldds->dampingMatrix();
    auto mass = lldds->massMatrix();

    double one_minus_theta, dt_sigma_star;
    auto time_step = _simulation->timeStep();
    double coeff = 0.5 * time_step * time_step;
    double one = 1.;
    double two = 2.;
    double omega2k, sigmak, massk;
    for (unsigned int k = 0; k < ndof; ++k) {
      massk = 1.;
      sigmak = 0.;
      omega2k = 0.;
      if (lldds->hasMassMatrix()) massk = mass(k, k);
      if (lldds->hasDampingMatrix()) sigmak = 0.5 * damp(k);
      if (lldds->hasStiffnessMatrix()) omega2k = omega2(k);
      compute_parameters(time_step, omega2k, sigmak, one_minus_theta, dt_sigma_star);
      iteration_matrix(k, k) =
          one / (massk + coeff * one_minus_theta * omega2k + dt_sigma_star);
      (*ds_work_vectors[siconos::integrators::MoreauJeanBilbaoOSI::ONE_MINUS_THETA])(k) =
          one_minus_theta;
      (*ds_work_vectors[siconos::integrators::MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR])(k) =
          two * dt_sigma_star;
    }
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanBilbaoOSI::initialize_iteration_"
        "matrix - Only "
        "implemented for LagrangianLinearDiagonalDS");
}

void siconos::integrators::MoreauJeanBilbaoOSI::compute_parameters(double time_step,
                                                                   double omega2, double sigma,
                                                                   double& one_minus_theta,
                                                                   double& dt_sigma_star) {
  // Computes:
  // 1 - theta
  // sigma_star = dt * sigma^*
  double ek = std::exp(-sigma * time_step);
  std::complex<double> buff =
      std::sqrt(std::complex<double>(sigma * sigma - omega2)) * time_step;
  std::complex<double> cAk = ek * (std::exp(buff) + std::exp(-buff));
  assert(std::imag(cAk) < 1e-12);
  double Ak = std::real(cAk);
  double one = 1.;
  double two = 2.;
  double one_over_two = 0.5;
  ek = ek * ek;
  double res = omega2 * time_step * time_step;
  one_minus_theta = one - two / res + Ak / (one + ek - Ak);
  dt_sigma_star = (one - ek) / (one + ek) * (one + one_over_two * res * one_minus_theta);
}

double siconos::integrators::MoreauJeanBilbaoOSI::computeResidu() { return 0.; }

void siconos::integrators::MoreauJeanBilbaoOSI::computeFreeState() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanBilbaoOSI::computeFreeState()")
  // This function computes "free" states of the DS belonging to this
  // Integrator. "Free" means without taking non-smooth effects into account.

  double time_step = _simulation->timeStep();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    // Nothing to be done if the osi is not linked to the ds
    if (!checkOSI(dsi)) continue;
    //
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& lldds = static_cast<siconos::modeling::LagrangianLinearDiagonalDS&>(*ds);
    auto& work_ds = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // Get velocity computed at the beginning of the time step.
    const auto& v_i = lldds.velocityMemory().getSiconosVector(0);
    const auto& q_i = lldds.qMemory().getSiconosVector(0);
    DEBUG_EXPR(q_i.display(););
    DEBUG_EXPR(v_i.display(););
    auto stiffness = lldds.stiffnessMatrix();
    // Get iteration matrix
    const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
    auto& inv_iteration_matrix = *_dynamicalSystemsGraph->properties(dsv).iterationMatrix;
    DEBUG_EXPR(inv_iteration_matrix.display(););
    // Get 2.*dt*sigma^*
    auto& two_dt_sigma_star =
        *work_ds[siconos::integrators::MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR];
    // buffer for vfree
    auto& vfree = *work_ds[siconos::integrators::MoreauJeanBilbaoOSI::VFREE];
    // Compute vfree
    auto dimension = lldds.dimension();
    for (decltype(dimension) k = 0; k < dimension; ++k)
      vfree(k) = v_i(k) - inv_iteration_matrix(k, k) * (time_step * stiffness(k) * q_i(k) +
                                                        two_dt_sigma_star(k) * v_i(k));
    DEBUG_EXPR(vfree.display(););
  }
  DEBUG_END("siconos::integrators::MoreauJeanBilbaoOSI::computeFreeState()")
}

struct siconos::integrators::MoreauJeanBilbaoOSI::_NSLEffectOnFreeOutput
    : public siconos::internal::SiconosVisitor {
  using SiconosVisitor::visit;

  siconos::nonsmooth_formulations::OneStepNSProblem* _osnsp{nullptr};
  siconos::modeling::Interaction& _inter;
  siconos::graphs::InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(const _NSLEffectOnFreeOutput&) = delete;
  _NSLEffectOnFreeOutput(siconos::nonsmooth_formulations::OneStepNSProblem* p,
                         siconos::modeling::Interaction& inter,
                         siconos::graphs::InteractionProperties& interProp)
      : _osnsp(p), _inter(inter), _interProp(interProp) {};

  void visit(const siconos::modeling::NewtonImpactNSL& nslaw) const override {
    auto e = nslaw.e();
    auto& osnsp_rhs =
        *(*_interProp.workVectors)[siconos::integrators::MoreauJeanBilbaoOSI::OSNSP_RHS];
    osnsp_rhs += e * _inter.y_k(_osnsp->inputOutputLevel());
  }
};

void siconos::integrators::MoreauJeanBilbaoOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  auto& indexSet = *osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  assert(indexSet.bundle(vertex_inter));
  // current interaction
  auto& inter = *indexSet.bundle(vertex_inter);
  // auto& DSlink = inter.linkToDSVariables();
  assert(inter.relation());
  auto& inter_work_block = *indexSet.properties(vertex_inter).workBlockVectors;
  // Get relation and non smooth law types
  // auto relationType = inter.relation()->getType();
  // auto relationSubType = inter.relation()->getSubType();
  // check relation type: done in initializeWorkVectorsForInteraction.
  // if(relationType != siconos::modeling::RelationType::Lagrangian ||
  // relationSubType != siconos::modeling::RelationSubType::LinearR)
  //   THROW_EXCEPTION("siconos::integrators::MoreauJeanBilbaoOSI::computeFreeOutput
  //   only Lagrangian Linear Relations are allowed.");

  auto sizeY = inter.nonSmoothLaw()->size();
  // buffer used to save output
  auto& x_free = *inter_work_block[siconos::integrators::MoreauJeanBilbaoOSI::xfree];
  auto& osnsp_rhs = *(*indexSet.properties(vertex_inter)
                           .workVectors)[siconos::integrators::MoreauJeanBilbaoOSI::OSNSP_RHS];

  auto lagr = std::dynamic_pointer_cast<siconos::modeling::LagrangianR>(inter.relation());
  assert(lagr && "MoreauJeanBilbaoOSI only implemented for Lagrangian systems.");

  auto H = lagr->jacobianhOver_q();
  siconos::algebra::matrixBlockVector_prod(H, x_free, osnsp_rhs, true);
  _NSLEffectOnFreeOutput nslEffectOnFreeOutput{osnsp, inter,
                                               indexSet.properties(vertex_inter)};
  inter.nonSmoothLaw()->accept(nslEffectOnFreeOutput);
}

void siconos::integrators::MoreauJeanBilbaoOSI::integrate(double& tinit, double& tend,
                                                          double& tout, int& notUsed) {
  THROW_EXCEPTION(
      "siconos::integrators::MoreauJeanBilbaoOSI::integrate - Not yet "
      "implemented for "
      "MoreauJeanBilbaoOSI.");
}

void siconos::integrators::MoreauJeanBilbaoOSI::updatePosition(
    siconos::modeling::DynamicalSystem& ds) {
  // --  Update current position for the last computed velocities --
  // q(i+1) = q(i) + dt v(i+1)
  double time_step = _simulation->timeStep();
  // get dynamical system
  auto& d = static_cast<siconos::modeling::LagrangianLinearDiagonalDS&>(ds);
  // Compute q
  //  -> get positions at the beginning of the time step
  const auto& qold = d.qMemory().getSiconosVector(0);
  // update positions
  *d.q() = time_step * d.velocity_read();
  *d.q() += qold;
  DEBUG_EXPR(q.display(););
}

void siconos::integrators::MoreauJeanBilbaoOSI::updateState(const unsigned int) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanBilbaoOSI::updateState(const unsigned "
      "int)");
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto& inv_iteration_matrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    // get dynamical system and work vector
    auto& lldds = static_cast<siconos::modeling::LagrangianLinearDiagonalDS&>(
        *_dynamicalSystemsGraph->bundle(*dsi));
    auto& work_ds = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto& vfree = *work_ds[siconos::integrators::MoreauJeanBilbaoOSI::VFREE];

    auto& v = *lldds.velocity();
    if (lldds.p(_levelMaxForInput) && lldds.p(_levelMaxForInput)->size() > 0) {
      v = *lldds.p(_levelMaxForInput);  // v = p
      if (lldds.boundaryConditions())
        for (const auto itindex : lldds.boundaryConditions()->velocityIndices()) {
          v.setValue(itindex, 0.0);
        }
      auto ndof = lldds.dimension();
      for (unsigned int k = 0; k < ndof; ++k)
        v(k) = vfree(k) + v(k) * inv_iteration_matrix(k, k);
    } else
      v = vfree;
    DEBUG_EXPR(v.display(););
    // Update positions with the last computed velocities.
    updatePosition(lldds);
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanBilbaoOSI::updateState(const unsigned "
      "int)");
}

void siconos::integrators::MoreauJeanBilbaoOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanBilbaoOSI OSI display ======\n";
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  if (_dynamicalSystemsGraph) {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------\n";
      std::cout << "--> W of dynamical system number " << ds->number() << ": \n";
      if (_dynamicalSystemsGraph->properties(*dsi).iterationMatrix)
        _dynamicalSystemsGraph->properties(*dsi).iterationMatrix->display();
      else
        std::cout << "-> nullptr\n";
    }
  }
  std::cout << "================================\n ";
}

void siconos::integrators::MoreauJeanBilbaoOSI::prepareNewtonIteration(double time) {
  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);
  }
}

bool siconos::integrators::MoreauJeanBilbaoOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction>"
      " inter, "
      "unsigned int i)\n");

  assert(i == 1);
  double h = _simulation->timeStep();
  auto y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  auto yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity
  bool _useGamma = false;
  double _gamma = 1.0 / 2.0;
  double _constraintActivationThreshold = 0.0;

  double gamma = 1.0 / 2.0;
  if (_useGamma) {
    gamma = _gamma;
  }
  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanBilbaoOSI::addInteractionInIndexSet of "
      "level = %i "
      "yref=%e, yDot=%e, y_estimated=%e.,  _constraintActivationThreshold "
      "=%e\n",
      i, y, yDot, y + gamma * h * yDot, _constraintActivationThreshold);
  y += gamma * h * yDot;
  assert(!std::isnan(y));
  DEBUG_EXPR_WE(if (y <= _constraintActivationThreshold) std::cout
                    << "siconos::integrators::MoreauJeanBilbaoOSI::"
                       "addInteractionInIndexSet ACTIVATE."
                    << y << "<= " << _constraintActivationThreshold << "\n";
                ;);
  return (y <= _constraintActivationThreshold);
}

bool siconos::integrators::MoreauJeanBilbaoOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  return !(addInteractionInIndexSet(inter, i));
}
