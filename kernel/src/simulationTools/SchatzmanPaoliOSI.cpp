
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
#include "SchatzmanPaoliOSI.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianLinearTIR.hpp"
#include "LagrangianR.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosMatrixOp.hpp"        // for scal
#include "SiconosMatrixVectorOp.hpp"  // for mat-vecprod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for subscal
#include "SiconosVisitor.hpp"
#include "SiconosMatrix.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// --- constructor from a set of data ---
siconos::integrators::SchatzmanPaoliOSI::SchatzmanPaoliOSI(double theta)
    : OneStepIntegrator{IntegratorType::SCHATZMANPAOLIOSI, 2, 0, 0, 0, 0}, _theta{theta} {
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY;
}

// --- constructor from a set of data ---
siconos::integrators::SchatzmanPaoliOSI::SchatzmanPaoliOSI(double theta, double gamma)
    : OneStepIntegrator{IntegratorType::SCHATZMANPAOLIOSI, 2, 0, 0, 0, 0},
      _theta{theta},
      _gamma{gamma},
      _useGamma{true} {
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY;
}

std::shared_ptr<siconos::algebra::SiconosMatrix> siconos::integrators::SchatzmanPaoliOSI::W(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds && "siconos::integrators::SchatzmanPaoliOSI::W(ds): ds == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
  ;
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::SchatzmanPaoliOSI::WBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds &&
         "siconos::integrators::SchatzmanPaoliOSI::WBoundaryConditions(ds): ds == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
      .WBoundaryConditions;
}

void siconos::integrators::SchatzmanPaoliOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  DEBUG_BEGIN(
      "siconos::integrators::SchatzmanPaoliOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");

  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type
  assert(std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds));
  if (auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    // buffers allocation (inside the graph)

    ds_work_vectors.resize(siconos::integrators::SchatzmanPaoliOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lltids->dimension());
    ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lltids->dimension());
    ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::LOCAL_BUFFER] =
        std::make_shared<siconos::algebra::SiconosVector>(lltids->dimension());
    auto q0 = lltids->q0();
    auto q = lltids->q();
    auto v0 = lltids->velocity0();
    auto velocity = lltids->velocity();

    // We first swap the initial value contained in q and v after initialization.
    lltids->swapInMemory();
    // we compute the new state values
    auto h = _simulation->timeStep();
    *q = *q0 + h * *v0;

    //*velocity=*velocity; we do nothing for the velocity
    lltids->swapInMemory();
  }

  // W initialization
  initializeIterationMatrixW(t, ds);

  for (auto k = _levelMinForInput; k < _levelMaxForInput + 1; k++) {
    ds->initializeNonSmoothInput(k);
  }

  DEBUG_EXPR(ds->display());
  DEBUG_END(
      "siconos::integrators::SchatzmanPaoliOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
}

void siconos::integrators::SchatzmanPaoliOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  auto& DSlink = inter.linkToDSVariables();

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::SchatzmanPaoliOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::SchatzmanPaoliOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();

  inter_work[siconos::integrators::SchatzmanPaoliOSI::OSNSP_RHS] =
      std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  auto relationType = relation.getType();

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2)))) {
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::initializeWorkVectorsForInteraction. The "
        "implementation is not correct for two different OSI for one interaction");
  }

  /* allocate and set work vectors for the osi */
  auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1);
    DSlink[siconos::modeling::LagrangianR::p0] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::p0]->insertPtr(lds.p(0));

    inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree] =
        std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree]->insertPtr(
        workVds1[siconos::integrators::SchatzmanPaoliOSI::FREE]);
  } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree] =
        std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree]->insertPtr(
        workVds1[siconos::integrators::SchatzmanPaoliOSI::FREE]);
  }

  if (ds1 != ds2) {
    auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree]->insertPtr(
          workVds2[siconos::integrators::SchatzmanPaoliOSI::FREE]);
      auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2);
      DSlink[siconos::modeling::LagrangianR::p0]->insertPtr(lds.p(0));
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree]->insertPtr(
          workVds2[siconos::integrators::SchatzmanPaoliOSI::FREE]);
    }
  }
}

void siconos::integrators::SchatzmanPaoliOSI::initializeIterationMatrixW(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending
  // on its type.

  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::initializeIterationMatrixW(t,ds) - ds == "
        "nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::initializeIterationMatrixW(t,ds) - ds "
        "does "
        "not belong to the OSI.");

  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);

  if (_dynamicalSystemsGraph->properties(dsv).W)
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::initializeIterationMatrixW(t,ds) - W(ds) "
        "is "
        "already in the map and has been initialized.");

  // Memory allocation for W
  auto h = _simulation->timeStep();
  auto sizeW = ds->dimension();
  if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    if (d->mass()) {
      _dynamicalSystemsGraph->properties(dsv).W =
          std::make_shared<siconos::algebra::SiconosMatrix>(*d->mass());  //*W = *d->mass();
    } else {
      _dynamicalSystemsGraph->properties(dsv).W =
          std::make_shared<siconos::algebra::SiconosMatrix>(sizeW, sizeW);
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }

    auto K = d->K();
    auto C = d->C();
    auto W = _dynamicalSystemsGraph->properties(dsv).W;
    if (C)  // W += 1/2.0*h*_theta *C
      siconos::algebra::scal(1 / 2.0 * h * _theta, *C, *W, false);

    if (K)  // W = h*h*_theta*_theta*K
      siconos::algebra::scal(h * h * _theta * _theta, *K, *W, false);

    // WBoundaryConditions initialization
    if (d->boundaryConditions()) initializeIterationMatrixWBoundaryConditions(d, dsv);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::initializeIterationMatrixW - only "
        "implemented for LagrangianLinearTIDS");

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
}

void siconos::integrators::SchatzmanPaoliOSI::initializeIterationMatrixWBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsv) {
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::"
        "initializeIterationMatrixWBoundaryConditions(t,ds) - ds == nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::"
        "initializeIterationMatrixWBoundaryConditions(t,ds) - ds does not belong to the "
        "OSI.");

  THROW_EXCEPTION(
      "siconos::integrators::SchatzmanPaoliOSI::"
      "initializeIterationMatrixWBoundaryConditions "
      "- not yet implemented .")
}

void siconos::integrators::SchatzmanPaoliOSI::computeWBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    siconos::algebra::SiconosMatrix& WBoundaryConditions) {
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  assert(ds &&
         "siconos::integrators::SchatzmanPaoliOSI::computeWBoundaryConditions(t,ds) - ds == "
         "nullptr");

  THROW_EXCEPTION(
      "siconos::integrators::SchatzmanPaoliOSI::computeWBoundaryConditions - not yet "
      "implemented ");
}

void siconos::integrators::SchatzmanPaoliOSI::computeW(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    siconos::algebra::SiconosMatrix& W) {
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  assert(ds && "siconos::integrators::SchatzmanPaoliOSI::computeW(t,ds) - ds == nullptr");

  // double h = _simulation->timeStep();

  if (std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    // Nothing: W does not depend on time.
  } else
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::computeW - only implemented for "
        "LagrangianLinearTIDS.");

  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}

double siconos::integrators::SchatzmanPaoliOSI::computeResidu() {
  DEBUG_BEGIN("siconos::integrators::SchatzmanPaoliOSI::computeResidu()\n");
  // This function is used to compute the residu for each "SchatzmanPaoliOSI-discretized"
  // dynamical system. It then computes the norm of each of them and finally return the
  // maximum value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  auto t = _simulation->nextTime();         // End of the time step
  auto told = _simulation->startingTime();  // Beginning of the time step
  auto h = t - told;                        // time step length

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

    // 1 - Lagrangian Non Linear Systems
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
      // ResiduFree =  M(-q_{k}+q_{k-1})  + h^2 (K q_k)+  h^2 C (\theta
      // \Frac{q_k-q_{k-1}}{2h}+ (1-\theta) v_k))  (1) This formulae is only valid for the
      // first computation of the residual for q = q_k otherwise the complete formulae must
      // be applied, that is ResiduFree   M(q-2q_{k}+q_{k-1})  + h^2 (K(\theta q+ (1-\theta)
      // q_k)))+  h^2 C (\theta \Frac{q-q_{k-1}}{2h}+ (1-\theta) v_k))  (2) for q != q_k, the
      // formulae (1) is wrong. in the sequel, only the equation (1) is implemented

      DEBUG_EXPR(d->display());
      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      const auto& q_k = d->qMemory().getSiconosVector(0);         // q_k
      const auto& q_k_1 = d->qMemory().getSiconosVector(1);       // q_{k-1}
      const auto& v_k = d->velocityMemory().getSiconosVector(0);  // v_k

      // --- ResiduFree computation Equation (1) ---
      auto& residuFree =
          *ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::RESIDU_FREE];
      auto& free = *ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::FREE];

      DEBUG_EXPR(free.display());
      DEBUG_EXPR(residuFree.display());
      residuFree.zero();
      double coeff;
      // -- No need to update W --

      residuFree = q_k_1;
      siconos::algebra::sub(residuFree, q_k, residuFree);
      if (d->mass())
        siconos::algebra::prod(*(d->mass()), residuFree,
                               residuFree);  // residuFree = M(-q_{k}+q_{k-1})

      auto K = d->K();
      if (K) {
        siconos::algebra::prod(h * h, *K, q_k, residuFree, false);  // residuFree += h^2*K*qi
      }

      auto C = d->C();
      if (C)
        siconos::algebra::prod(
            h * h, *C, (1.0 / (2.0 * h) * _theta * (q_k - q_k_1) + (1.0 - _theta) * v_k),
            residuFree, false);
      // residufree += h^2 C (\theta \Frac{q-q_{k-1}}{2h}+ (1-\theta) v_k))

      auto Fext = d->fExt();
      if (Fext) {
        // computes Fext(ti)
        d->computeFext(told);
        coeff = -h * h * (1 - _theta);
        siconos::algebra::scal(coeff, *Fext, residuFree,
                               false);  // residufree -= h^2*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFext(t);
        coeff = -h * h * _theta;
        siconos::algebra::scal(coeff, *Fext, residuFree,
                               false);  // residufree -= h^2*_theta * fext(ti+1)
      }

      DEBUG_EXPR(free.display());
      DEBUG_EXPR(residuFree.display());

      //  std::cout << "siconos::integrators::SchatzmanPaoliOSI::ComputeResidu
      //  LagrangianLinearTIDS residufree :"  << std::endl;
      // residuFree->display();

      free = residuFree;              // copy residuFree in Workfree
      if (d->p(0)) free -= *d->p(0);  // Compute Residu in Workfree Notation !!
      DEBUG_EXPR(free.display());
      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree.norm2();
    } else
      THROW_EXCEPTION(
          "siconos::integrators::SchatzmanPaoliOSI::computeResidu - only implemented for "
          "LagrangianLinearTIDS systems.");

    if (normResidu > maxResidu) maxResidu = normResidu;
  }
  DEBUG_END("siconos::integrators::SchatzmanPaoliOSI::computeResidu()\n");
  return maxResidu;
}

void siconos::integrators::SchatzmanPaoliOSI::computeFreeState() {
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  // double t = _simulation->nextTime(); // End of the time step
  // double told = _simulation->startingTime(); // Beginning of the time step
  // double h = t-told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  siconos::algebra::SiconosVector *rold =
  //  static_cast<siconos::algebra::SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  // W SchatzmanPaoliOSI matrix of the current DS.
  std::shared_ptr<siconos::algebra::SiconosMatrix> W;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;

    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    W = _dynamicalSystemsGraph->properties(*dsi)
            .W;  // Its W SchatzmanPaoliOSI matrix of iteration.

    // 1 - Lagrangian Non Linear Systemsv
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // vFree = v_i + W^{-1} ResiduFree    // with
      // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 +
      // h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      const auto& qold = d->qMemory().getSiconosVector(0);  // q_k
      //   auto vold =
      //   d->velocityMemory()->getSiconosVector(0); //v_k

      // --- ResiduFree computation ---

      // vFree pointer is used to compute and save ResiduFree in this first step.
      auto& residuFree =
          *ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::RESIDU_FREE];
      auto& qfree = *ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::FREE];

      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      qfree = residuFree;

      siconos::algebra::solveInPlace(*W, qfree);
      qfree *= -1.0;
      qfree += qold;
    } else
      THROW_EXCEPTION(
          "siconos::integrators::SchatzmanPaoliOSI::computeFreeState - Only implemented for "
          "LagrangianLinearTIDS systems.");
  }
}

void siconos::integrators::SchatzmanPaoliOSI::prepareNewtonIteration(double time) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    computeW(time, ds, *_dynamicalSystemsGraph->properties(*dsi).W);
  }
  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);
  }
}

struct siconos::integrators::SchatzmanPaoliOSI::_NSLEffectOnFreeOutput
    : public siconos::internal::SiconosVisitor {
  using siconos::internal::SiconosVisitor::visit;

  siconos::nonsmooth_formulations::OneStepNSProblem* _osnsp{nullptr};
  std::shared_ptr<siconos::modeling::Interaction> _inter;
  siconos::graphs::InteractionProperties& _interProp;

  _NSLEffectOnFreeOutput(siconos::nonsmooth_formulations::OneStepNSProblem* p,
                         std::shared_ptr<siconos::modeling::Interaction> inter,
                         siconos::graphs::InteractionProperties& interProp)
      : _osnsp(p), _inter(inter), _interProp(interProp){};

  void visit(const siconos::modeling::NewtonImpactNSL& nslaw) const override {
    double e;
    e = nslaw.e();
    auto sizeY = _inter->nonSmoothLaw()->size();
    std::vector<std::size_t> subCoord{0, sizeY, 0, sizeY};
    // Only the normal part is multiplied by e
    const auto& y_k_1(_inter->yMemory(_osnsp->inputOutputLevel()).getSiconosVector(1));

    DEBUG_PRINTF("_osnsp->inputOutputLevel() = %i \n ", _osnsp->inputOutputLevel());
    DEBUG_EXPR(y_k_1.display());
    ;
    auto& osnsp_rhs =
        *(*_interProp.workVectors)[siconos::integrators::SchatzmanPaoliOSI::OSNSP_RHS];
    siconos::algebra::subscal(e, y_k_1, osnsp_rhs, subCoord, false);
  }

  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nslaw) const override {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    const auto& y_k_1(_inter->yMemory(_osnsp->inputOutputLevel()).getSiconosVector(1));
    auto& osnsp_rhs =
        *(*_interProp.workVectors)[siconos::integrators::SchatzmanPaoliOSI::OSNSP_RHS];
    osnsp_rhs(0) += e * (y_k_1)(0);
  }
  void visit(const siconos::modeling::EqualityConditionNSL& nslaw) const override { ; }
  void visit(const siconos::modeling::MixedComplementarityConditionNSL& nslaw) const override {
    ;
  }
};

void siconos::integrators::SchatzmanPaoliOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  DEBUG_BEGIN(
      "siconos::integrators::SchatzmanPaoliOSI::computeFreeOutput(siconos::graphs::"
      "InteractionsGraph::VDescriptor& vertex_inter, "
      "siconos::nonsmooth_formulations::OneStepNSProblem* "
      "osnsp)\n");
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */

  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);
  auto allOSNS = _simulation->oneStepNSProblems();
  auto& inter_work_block = *indexSet->properties(vertex_inter).workBlockVectors;

  // Get relation and non smooth law types
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();
  auto sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;

  std::vector<std::size_t> coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  std::shared_ptr<siconos::algebra::SiconosMatrix> C;
  std::shared_ptr<siconos::algebra::SiconosMatrix> D;
  std::shared_ptr<siconos::algebra::SiconosMatrix> F;
  std::shared_ptr<siconos::algebra::BlockVector> deltax;
  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::SchatzmanPaoliOSI::OSNSP_RHS];

  std::shared_ptr<siconos::algebra::SiconosVector> e;
  auto Xfree = inter_work_block[siconos::integrators::SchatzmanPaoliOSI::xfree];
  ;
  assert(Xfree);

  auto mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());
  DEBUG_EXPR(inter->display(););
  if (relationSubType == siconos::modeling::RelationSubType::LinearTIR) {
    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() != osnsp)
      THROW_EXCEPTION(
          "siconos::integrators::SchatzmanPaoliOSI::computeFreeOutput not yet implemented "
          "for "
          "siconos::simulation::SICONOS_OSNSP ");

    C = mainInteraction->relation()->C();

    if (C) {
      assert(Xfree);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interactionfor each ds of the
      // current Interaction.

      if (_useGammaForRelation) {
        assert(deltax);
        siconos::algebra::subprod(*C, *deltax, osnsp_rhs, coord, true);
      } else {
        siconos::algebra::subprod(*C, *Xfree, osnsp_rhs, coord, true);
      }
    }
    auto ltir = std::static_pointer_cast<siconos::modeling::LagrangianLinearTIR>(
        mainInteraction->relation());
    e = ltir->e();
    if (e) {
      osnsp_rhs += *e;
    }
  } else
    THROW_EXCEPTION(
        "siconos::integrators::SchatzmanPaoliOSI::ComputeFreeOutput not yet implemented  "
        "for "
        "relation of Type : " +
        siconos::tools::enum_to_string(relationType));

  if (inter->relation()->getSubType() == siconos::modeling::RelationSubType::LinearTIR) {
    auto nslEffectOnFreeOutput = std::make_shared<_NSLEffectOnFreeOutput>(
        osnsp, inter, indexSet->properties(vertex_inter));
    inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
  }

  DEBUG_END(
      "siconos::integrators::SchatzmanPaoliOSI::computeFreeOutput(siconos::graphs::"
      "InteractionsGraph::VDescriptor& vertex_inter, "
      "siconos::nonsmooth_formulations::OneStepNSProblem* "
      "osnsp)\n");
}
void siconos::integrators::SchatzmanPaoliOSI::integrate(double& tinit, double& tend,
                                                        double& tout, int&) {
  THROW_EXCEPTION(
      "siconos::integrators::SchatzmanPaoliOSI::integrate - not yet implemented :");
}

void siconos::integrators::SchatzmanPaoliOSI::updateState(const unsigned int) {
  DEBUG_BEGIN("siconos::integrators::SchatzmanPaoliOSI::updateState(const unsigned int )\n");

  double h = _simulation->timeStep();

  auto RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  std::shared_ptr<siconos::algebra::SiconosMatrix> W;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    W = _dynamicalSystemsGraph->properties(*dsi).W;

    // 1 - Lagrangian Systems
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // get dynamical system
      auto& qfree = *ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::FREE];

      //    siconos::algebra::SiconosVector *vfree = d->velocityFree();
      auto& q = *d->q();
      bool baux =
          ((not(std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds))) &&
           useRCC && _simulation->relativeConvergenceCriterionHeld());

      // To compute q, we solve W(q - qfree) = p
      if (d->p(_levelMaxForInput)) {
        q = *d->p(_levelMaxForInput);  // q = p
        siconos::algebra::solveInPlace(*W, q);
      } else
        q.zero();

      q += qfree;

      // Computation of the velocity

      auto& v = *d->velocity();
      const auto& q_k_1 = d->qMemory().getSiconosVector(1);  // q_{k-1}

      //  std::cout << "siconos::integrators::SchatzmanPaoliOSI::updateState - q_k_1 ="
      //  <<std::endl;
      // q_k_1->display();
      //  std::cout << "siconos::integrators::SchatzmanPaoliOSI::updateState - q ="
      //  <<std::endl;
      // q->display();

      v = 1.0 / (2.0 * h) * (q - q_k_1);
      //  std::cout << "siconos::integrators::SchatzmanPaoliOSI::updateState - v ="
      //  <<std::endl;
      // v->display();

      // int bc=0;
      // auto columntmp =
      // std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));

      // if (d->boundaryConditions())
      // {
      //   for (const auto itindex : d->boundaryConditions()->velocityIndices()){
      //   {
      //     _WBoundaryConditionsMap[ds]->getCol(bc,*columntmp);
      //     /*\warning we assume that W is symmetric in the Lagrangian case*/
      //     double value = - siconos::algebra::inner_prod(*columntmp, *v);
      //     value += (d->p(level))->getValue(itindex);
      //     /* \warning the computation of reactionToBoundaryConditions take into
      //        account the contact impulse but not the external and internal forces.
      //        A complete computation of the residue should be better */
      //     d->reactionToBoundaryConditions()->setValue(bc,value) ;
      //     bc++;
      //   }

      if (baux) {
        double ds_norm_ref = 1. + ds->x0()->norm2();  // Should we save this in the graph?
        *ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::LOCAL_BUFFER] -= q;
        auto aux =
            (ds_work_vectors[siconos::integrators::SchatzmanPaoliOSI::LOCAL_BUFFER]->norm2()) /
            ds_norm_ref;
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    }
    // 2 - Newton Euler Systems
    else {  // if (dsType == Type::NewtonEulerDS) {
      //  // get dynamical system
      //       auto d = std::static_pointer_cast<NewtonEulerDS> (ds);
      //       auto v = d->velocity();
      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"siconos::integrators::SchatzmanPaoliOSI::updatestate prev
      //       v"<<endl; v->display();
      // #endif

      //       /*d->p has been fill by the Relation->computeInput, it contains
      //            B \lambda _{k+1}*/
      //       *v = *d->p(level); // v = p
      //       d->luW()->Solve(*v);

      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"siconos::integrators::SchatzmanPaoliOSI::updatestate hWB
      //       lambda"<<endl; v->display();
      // #endif

      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"siconos::integrators::SchatzmanPaoliOSI::updatestate work
      //       free"<<endl; std::cout<<"siconos::integrators::SchatzmanPaoliOSI::updatestate
      //       new v"<<endl; v->display();
      // #endif
      //       //compute q
      //       //first step consists in computing  \dot q.
      //       //second step consists in updating q.
      //       //
      //       auto T = d->T();
      //       auto dotq = d->dotq();
      //       siconos::algebra::prod(*T,*v,*dotq,true);
      //       // std::cout<<"siconos::integrators::SchatzmanPaoliOSI::updateState v"<<endl;
      //       // v->display();
      //       // std::cout<<"siconos::integrators::SchatzmanPaoliOSI::updateState
      //       dotq"<<endl;
      //       // dotq->display();

      //       auto q = d->q();

      //       //  -> get previous time step state
      //       auto dotqold =
      //       d->dotqMemory()->getSiconosVector(0);
      //       auto qold =
      //       d->qMemory()->getSiconosVector(0);
      //       // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
      //       double coeff = h*_theta;
      //       scal(coeff, *dotq, *q) ; // q = h*theta*v
      //       coeff = h*(1-_theta);
      //       scal(coeff,*dotqold,*q,false); // q += h(1-theta)*vold
      //       *q += *qold;
      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"new q before normalizing"<<endl;
      //       q->display();
      // #endif

      //       //q[3:6] must be normalized
      //       d->normalizeq();
      //       dotq->setValue(3,(q->getValue(3)-qold->getValue(3))/h);
      //       dotq->setValue(4,(q->getValue(4)-qold->getValue(4))/h);
      //       dotq->setValue(5,(q->getValue(5)-qold->getValue(5))/h);
      //       dotq->setValue(6,(q->getValue(6)-qold->getValue(6))/h);
      //       d->updateT();
      THROW_EXCEPTION(
          "siconos::integrators::SchatzmanPaoliOSI::updateState - only implemented for "
          "LagrangianLinearTIDS");
    }
  }
  DEBUG_END("siconos::integrators::SchatzmanPaoliOSI::updateState(const unsigned int)\n");
}

void siconos::integrators::SchatzmanPaoliOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== SchatzmanPaoliOSI OSI display ======\n";

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--------------------------------\n";
    std::cout << "--> W of dynamical system number " << ds->number() << ":\n";
    if (_dynamicalSystemsGraph->properties(*dsi).W)
      _dynamicalSystemsGraph->properties(*dsi).W->display();
    else
      std::cout << "-> nullptr"
                << "\n";
    std::cout << "--> and corresponding theta is: " << _theta << "\n";
  }
  std::cout << "================================\n";
}
