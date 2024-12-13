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
#include "EulerMoreauOSI.hpp"

#include "BlockVector.hpp"
#include "FirstOrderLinearDS.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderNonLinearR.hpp"
#include "FirstOrderType1R.hpp"
#include "FirstOrderType2R.hpp"
#include "Interaction.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for prod and subprod
#include "Simulation.hpp"
#include "Topology.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

// --- constructor with theta parameter value  ---
siconos::integrators::EulerMoreauOSI::EulerMoreauOSI(double theta)
    : OneStepIntegrator(IntegratorType::EULERMOREAUOSI, 1, 0, 0, 0, 0), _theta(theta) {}

// --- constructor from a set of data ---
siconos::integrators::EulerMoreauOSI::EulerMoreauOSI(double theta, double gamma)
    : OneStepIntegrator(IntegratorType::EULERMOREAUOSI, 1, 0, 0, 0, 0),
      _theta{theta},
      _gamma{gamma},
      _useGamma{true} {}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::EulerMoreauOSI::IterationMatrixBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds &&
         "siconos::integrators::EulerMoreauOSI::IterationMatrixBoundaryConditions(ds): ds == "
         "nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
      .iterationMatrixBoundaryConditions;
}

void siconos::integrators::EulerMoreauOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(siconos::integrators::EulerMoreauOSI::WORK_LENGTH);

  // Check dynamical system type
  auto fods = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds);
  assert(fods);

  // Compute W (iteration matrix)
  initializeIterationMatrix(t, ds);

  // buffers allocation (into the graph)
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU_FREE] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::FREE] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());

  // Update dynamical system components (for memory swap).
  fods->computefVector(*fods->x(), t);  // Only fold is concerned, for FirstOrderNonLinearDS.
  // Update memory buffers
  ds->swapInMemory();
}

void siconos::integrators::EulerMoreauOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::EulerMoreauOSI::WORK_INTERACTION_LENGTH);
  }
  if (!interProp.workMatrices) {
    interProp.workMatrices =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>>>(
            siconos::integrators::EulerMoreauOSI::MAT_WORK_LENGTH);
  }
  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::EulerMoreauOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_mat = *interProp.workMatrices;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();

  auto sizeY = inter.dimension();
  auto sizeOfDS = inter.getSizeOfDS();
  auto relationType = relation.getType();
  auto relationSubType = relation.getSubType();
  inter_work[siconos::integrators::EulerMoreauOSI::OSNSP_RHS] =
      std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  assert((relationType == siconos::modeling::RelationType::FirstOrder) &&
         "EulerMoreauOSI only implemented for first-order type relations");

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken in to account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;

    inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y] =
        std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA] =
        std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    inter_work[siconos::integrators::EulerMoreauOSI::YOLD] =
        std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD] =
        std::make_shared<siconos::algebra::SiconosVector>(sizeY);

    if (relationSubType == siconos::modeling::RelationSubType::NonLinearR ||
        relationSubType == siconos::modeling::RelationSubType::Type2R) {
      inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA] =
          std::make_shared<siconos::algebra::BlockVector>(1, sizeOfDS);

      // Work vector used during computeResiduInput and updateInput
      inter_work[siconos::integrators::EulerMoreauOSI::WORK_DS] =
          std::make_shared<siconos::algebra::SiconosVector>(sizeOfDS);

      inter_work_mat[siconos::integrators::EulerMoreauOSI::MAT_KHAT] =
          std::make_shared<siconos::algebra::SiconosMatrix>(sizeOfDS, sizeY);
    }

    if (!inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]) {
      inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE] =
          std::make_shared<siconos::algebra::BlockVector>();
      inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
          workVds1[siconos::integrators::EulerMoreauOSI::FREE]);
    } else
      inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->setVectorPtr(
          0, workVds1[siconos::integrators::EulerMoreauOSI::FREE]);

    if (!inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]) {
      inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS] =
          std::make_shared<siconos::algebra::BlockVector>();
      inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
          workVds1[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
    } else
      inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->setVectorPtr(
          0, workVds1[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
    if (!inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]) {
      inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X] =
          std::make_shared<siconos::algebra::BlockVector>();
      inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
          workVds1[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
    } else
      inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->setVectorPtr(
          0, workVds1[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
  }
  DEBUG_PRINTF("ds1->number() %i\n", ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n", ds2->number());

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");

    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken in to account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      if (relationType == siconos::modeling::RelationType::FirstOrder) {
        if (!inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]) {
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE] =
              std::make_shared<siconos::algebra::BlockVector>();
          // dummy insertion to reserve first vector for ds1
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::FREE]);
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::FREE]);
        } else
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::FREE]);

        if (!inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]) {
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS] =
              std::make_shared<siconos::algebra::BlockVector>();
          // dummy insertion to reserve first vector for ds1
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
        } else
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);

        if (!inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]) {
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X] =
              std::make_shared<siconos::algebra::BlockVector>();
          // dummy insertion to reserve first vector for ds1
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
        } else
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
      }
    }
  }
}

void siconos::integrators::EulerMoreauOSI::initializeIterationMatrix(
    double time, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending on
  // its type.
  assert(ds);
  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrix(t,ds) - ds does not "
        "belong to the OSI.");

  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);

  assert(!_dynamicalSystemsGraph->properties(dsv).iterationMatrix);

  double timeStep = _simulation->timeStep();

  // Memory allocation for W
  _dynamicalSystemsGraph->properties(dsv).iterationMatrix =
      std::make_shared<siconos::algebra::SiconosMatrix>(ds->dimension(), ds->dimension());

  auto iterationMat = _dynamicalSystemsGraph->properties(dsv).iterationMatrix;
  auto fods = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds);
  assert(fods &&
         "EulerMoreauOSI::initializeIterationMatrix implemented only for first order ds");

  // Linear time-invariant system: W constant
  // if (auto fods = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds)) {
  if (fods->isTimeInvariant()) {
    auto foltids = std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds);
    if (fods->hasMMatrix())
      *iterationMat = fods->MMatrix();
    else
      iterationMat->setIdentity();

    if (foltids->hasA()) *iterationMat -= _simulation->timeStep() * _theta * foltids->A();

    //  if (_useGamma)
    {
      siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std::tie(oei, oeiend) = _dynamicalSystemsGraph->out_edges(dsv); oei != oeiend;
           ++oei) {
        auto inter = _dynamicalSystemsGraph->bundle(*oei);
        //      ivd = indexSet.descriptor(inter);
        auto& rel = static_cast<siconos::modeling::FirstOrderR&>(*inter->relation());
        if (rel.hasJacobiangOver_state()) {
          auto K = rel.jacobiangOver_state();
          *iterationMat -= _simulation->timeStep() * _gamma * K;
        }
      }
    }
    // LU Factorisation
    _dynamicalSystemsGraph->properties(dsv).LUW =
        std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*iterationMat);
  }
  //}
  // In all other cases, we use computeIterationMatrix
  else
    // if  (auto d = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds)) {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    computeIterationMatrix(time, *ds, dsv,
                           *_dynamicalSystemsGraph->properties(dsv).iterationMatrix);
}

void siconos::integrators::EulerMoreauOSI::initializeIterationMatrixBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // This function:
  // - allocate memory for a matrix IterationMatrixBoundaryConditions
  // - insert this matrix into IterationMatrixBoundaryConditionsMap with ds as a key

  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixBoundaryConditions("
        "t,"
        "ds) - ds == nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrix(t,ds) - ds does "
        "not "
        "belong to the OSI.");

  THROW_EXCEPTION(
      "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixBoundaryConditions - "
      "not yet implemented.");
}

void siconos::integrators::EulerMoreauOSI::computeIterationMatrixBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // Compute IterationMatrixBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, IterationMatrixBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixBoundaryConditions.

  assert(ds &&
         "siconos::integrators::EulerMoreauOSI::computeIterationMatrixBoundaryConditions(t,ds)"
         " - ds == "
         "nullptr");

  // unsigned int dsN = ds->number();
  THROW_EXCEPTION(
      "siconos::integrators::EulerMoreauOSI::computeIterationMatrixBoundaryConditions - not "
      "yet "
      "implemented.");
}

void siconos::integrators::EulerMoreauOSI::computeIterationMatrix(
    double time, siconos::modeling::DynamicalSystem& ds,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsv,
    siconos::algebra::SiconosMatrix& iterationMatrix) {
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeIterationMatrix(...)\n");
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, we assume that memory has been allocated for iterationMatrix
  // (call to initializeIterationMatrix)

  auto h = _simulation->timeStep();
  auto fods = dynamic_cast<siconos::modeling::FirstOrderNonLinearDS*>(&ds);
  // No need to check if the system is first-order. This must have been done during init.

  // 1 - First order linear and time-invariant coeff systems, W constant, nothing to be done
  if (fods->isTimeInvariant()) return;

  // W =  M - h*_theta* [jacobian_x f(t,x,z)]
  // Copy M or I if M is Null into W

  if (fods->hasMMatrix()) {
    fods->computeMMatrix(time);
    iterationMatrix = fods->MMatrix();
  } else
    iterationMatrix.setIdentity();

  // 2 - First order linear systems with coeff depending on time, update A matrix
  if (auto folds = dynamic_cast<siconos::modeling::FirstOrderLinearDS*>(&ds)) {
    if (folds->hasA()) {
      folds->computeA(time);
      iterationMatrix -= h * _theta * folds->A();
    }
  }
  // 3 - // nonlinear general case
  else {
    if (fods->hasJacobianfOver_x()) {
      fods->computeJacobianfOver_x(*fods->x(), time);
      // Add -h*_theta*jacobianfx to W
      iterationMatrix -= h * _theta * fods->jacobianfOver_x();
    }
  }

  DEBUG_EXPR(W.display(););

  //  if (_useGamma)
  {
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = _dynamicalSystemsGraph->out_edges(dsv); oei != oeiend;
         ++oei) {
      auto inter = _dynamicalSystemsGraph->bundle(*oei);
      auto& rel = static_cast<siconos::modeling::FirstOrderR&>(*inter->relation());
      if (rel.hasJacobiangOver_state()) {
        iterationMatrix -= h * _gamma * rel.jacobiangOver_state();
      }
    }
  }
  // LU Factorisation
  _dynamicalSystemsGraph->properties(dsv).LUW =
      std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(iterationMatrix);

  DEBUG_EXPR(W.display());
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeIterationMatrix(...)\n");
}

void siconos::integrators::EulerMoreauOSI::computeKhat(
    siconos::modeling::Interaction& inter, siconos::algebra::SiconosMatrix& m,
    std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>>& workM, double h) const {
  auto relationType = inter.relation()->getType();
  if ((relationType == siconos::modeling::RelationType::FirstOrder) &&
      (workM[siconos::integrators::EulerMoreauOSI::MAT_KHAT])) {
    auto& rel = static_cast<siconos::modeling::FirstOrderR&>(*inter.relation());
    if (rel.hasJacobiangOver_state()) {
      *workM[siconos::integrators::EulerMoreauOSI::MAT_KHAT] = rel.jacobiangOver_state() * m;
      *workM[siconos::integrators::EulerMoreauOSI::MAT_KHAT] *= h;
    }
  }
}

double siconos::integrators::EulerMoreauOSI::computeResidu() {
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeResidu()\n");
  // This function is used to compute the residu for each "EulerMoreauOSI-discretized"
  // dynamical system. It then computes the norm of each of them and finally return the
  // maximum value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double time = _simulation->nextTime();      // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double h = time - told;                     // time step length

  DEBUG_PRINTF("nextTime %f\n", time);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  double maxResidu = 0;
  double normResidu = maxResidu;

  auto htheta = h * _theta;
  auto h_one_minus_theta = h * (1. - _theta);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary
    // anymore Maurice will do that with subgraph :)
    auto& residuFree = *ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU_FREE];
    auto& residu = *ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU];

    // Check if ds are first-order
    auto fods = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds);
    assert(fods && "EulerMoreau OSI only implemented for 1st order ds.");

    // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
    // Residu = ResiduFree

    // Initialize rFree with M(x_{k+1}^{\alpha} - x_k)
    if (fods->hasMMatrix()) {
      fods->computeMMatrix(time);  // does nothing if time-invariant
      residuFree = fods->MMatrix() * (*(fods->x()) - fods->xMemory().getSiconosVector(0));
    } else
      residuFree = *(fods->x()) - fods->xMemory().getSiconosVector(0);

    if (auto folds = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds)) {
      if (folds->isTimeInvariant()) {  // linear Systems with time invariant coefficients
        // 1. R_{free} -= h * b
        if (folds->hasbVector()) residuFree -= h * folds->bVector();

        // 2. residuFree -= h * A (\theta x_{k+1}^{\alpha} + (1-\theta) x_k)
        if (folds->hasA()) {
          residuFree -= h_one_minus_theta * folds->A() * folds->xMemory().getSiconosVector(0);
          residuFree -= htheta * folds->A() * *(folds->x());
        }
      } else {  // First Order Linear ds with time-dependant coeff
        // Note: indices k/k+1 corresponds to value at the beginning/end of the time step.
        // Newton iterate are x and r

        // computes f(t_k,x_k)
        // No fold in FirstOrderLinearDS.

        if (folds->hasA()) {
          folds->computeA(told);
          residuFree -= h_one_minus_theta * folds->A() * folds->xMemory().getSiconosVector(0);
          folds->computeA(time);
          residuFree -= htheta * folds->A() * *(folds->x());
        }
        if (folds->hasbVector()) {
          folds->computeb(told);
          residuFree -= h_one_minus_theta * folds->bVector();
          folds->computeb(time);
          residuFree -= htheta * folds->bVector();
        }
      }
    } else {  // First order non Linear Systems, general case
      if (fods->hasfVector()) {
        // fods->computefVector(*fods->xMemory().getSiconosVector(0), time);
        //  Not required, we have fold
        residuFree -= h_one_minus_theta * fods->fold();
        fods->computefVector(*fods->x(), time);
        residuFree -= htheta * fods->fVector();
      }
    }

    // now we compute the residu = residuFree - h*gamma*r - h*(1-gamma)r_k
    residu = residuFree;

    if (!fods->isTimeInvariant()) {
      if (!_useGamma)  // no gamma
      {
        DEBUG_EXPR(fonlds->r()->display(););
        DEBUG_EXPR(residu.display());
        residu -= h * *fods->r();
      } else {
        residu -= h * _gamma * *fods->r();
        residu -= h * (1. - _gamma) * fods->rMemory().getSiconosVector(0);
      }
    }

    normResidu = residu.norm2();
    DEBUG_EXPR(residu.display());
    DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::computeResidu final residuFree\n");
    DEBUG_EXPR(residuFree.display());

    if (normResidu > maxResidu) maxResidu = normResidu;
  }

  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeResidu()\n");
  return maxResidu;
}

void siconos::integrators::EulerMoreauOSI::computeFreeState() {
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeFreeState()\n");

  double t = _simulation->nextTime();         // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double h = t - told;                        // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  siconos::algebra::SiconosVector *rold =
  //  static_cast<siconos::algebra::SiconosVector*>(d.rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary
    // anymore Maurice will do that with subgraph :)

    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    auto& W = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;

    // No need to check if the DS is a first order: it must have been done during init steps.
    auto d = std::static_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds);

    // xfree =  x - W^{-1} (ResiduFree - h(1-gamma)*rold)
    // with ResiduFree = = M(x - x_k) - h*theta*f(t_{k+1}, x) - h*(1-theta)*f(t_k, x_k)

    // to be updated at current time: W, f
    // fold is f at t_k
    // not time dependant: M
    // Get state i (previous time step) from Memories -> var. indexed with "Old"
    //    std::shared_ptr<siconos::algebra::SiconosVector> xold =
    //    d->xMemory()->getSiconosVector(0); // xi

    auto& x = *d->x();  // x = x_k or x = x_{k+1}^{\alpha}
    // xfree gets ResiduFree at first
    auto& xfree = *ds_work_vectors[siconos::integrators::EulerMoreauOSI::FREE];
    xfree = *ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU_FREE];

    DEBUG_PRINT(
        "siconos::integrators::EulerMoreauOSI::computeFreeState xfree <- residuFree\n");
    DEBUG_EXPR(xfree.display());

    if (_useGamma) {
      const auto& rold = d->rMemory().getSiconosVector(0);
      xfree -= h * (1 - _gamma) * rold;  //  xfree += -h(1-gamma)*rold
    }

    // At this point xfree = (ResiduFree - h(1-gamma)*rold)
    // -> Solve WX = xfree and set xfree = X
    // LUW must be uptodate in the graph (update during computeIterationMatrix call)
    xfree = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(xfree);

    // at this point, xfree = W^{-1} (ResiduFree - h(1-gamma)*rold)
    // -> compute real xfree = x - W^{-1} (ResiduFree - h(1-gamma)*rold)
    xfree *= -1.0;
    xfree += x;

    DEBUG_EXPR(xfree.display());

    // now the crazy intermediate variables
    // xPartialNS was updated before this fonction call
    // It constains either 0 (first Newton iterate)
    // or g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
    auto& xPartialNS =
        *ds_work_vectors[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION];
    DEBUG_PRINT(
        "siconos::integrators::EulerMoreauOSI::computeFreeState xPartialNS from "
        "Interaction\n");
    DEBUG_EXPR(xPartialNS.display());

    // -> Solve WX = g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
    // and set xPartialNS = X
    xPartialNS = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(xPartialNS);
    xPartialNS *= h;

    // compute real xPartialNS = xfree + ...
    xPartialNS += xfree;
    DEBUG_PRINT(
        "siconos::integrators::EulerMoreauOSI::computeFreeState xPartialNS real value\n");
    DEBUG_EXPR(xPartialNS.display());

    // deltaxForRelation = (\widetilde{K}_{k+1}^{\alpha})^{-1} xPartialNS - x
    auto& deltaxForRelation =
        *ds_work_vectors[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION];
    deltaxForRelation = xPartialNS;

    deltaxForRelation -= x;

    DEBUG_EXPR(deltaxForRelation.display());

    // have a look at the end of the DevNotes for this part
    if (_useGammaForRelation) {
      if (not(std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds))) {
        THROW_EXCEPTION(
            "siconos::integrators::EulerMoreauOSI::computeFreeState - "
            "_useGammaForRelation "
            "== true is only implemented for FirstOrderLinearDS");
      }
      deltaxForRelation = _gamma * xfree;
      const auto& xold = d->xMemory().getSiconosVector(0);
      deltaxForRelation += (1. - _gamma) * xold;
    }

    // some output
    DEBUG_EXPR(xfree.display(););
    DEBUG_EXPR(xPartialNS.display(););
    DEBUG_EXPR(deltaxForRelation.display(););
  }
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeFreeState()\n");
}

void siconos::integrators::EulerMoreauOSI::prepareNewtonIteration(double time) {
  // XXX TMP hack -- xhub
  // we have to iterate over the edges of the DSG0 -> the following won't be necessary
  // anymore Maurice will do that with subgraph :)
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto dsv = _dynamicalSystemsGraph->descriptor(ds);
    auto W = _dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    computeIterationMatrix(time, *ds, dsv, *W);
  }

  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      auto& inter = *indexSet0->bundle(*ui);
      auto& interProp = indexSet0->properties(*ui);

      auto& inter_work = *interProp.workVectors;
      auto& inter_work_block = *interProp.workBlockVectors;

      auto relationType = inter.relation()->getType();
      auto relationSubType = inter.relation()->getSubType();
      if (relationSubType == siconos::modeling::RelationSubType::NonLinearR ||
          relationSubType == siconos::modeling::RelationSubType::Type2R) {
        auto& relation = static_cast<siconos::modeling::FirstOrderR&>(*inter.relation());
        auto& xPartialNS =
            *inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS];

        if (relation.hasJacobiangOver_lambda())
          *inter_work[siconos::integrators::EulerMoreauOSI::WORK_DS] =
              relation.jacobiangOver_lambda() * *inter.lambda(0);
        xPartialNS = *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA];
        xPartialNS -= *inter_work[siconos::integrators::EulerMoreauOSI::WORK_DS];
      }
    }
  }
  // update alpha variables
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet0->bundle(*ui);
    auto& interProp = indexSet0->properties(*ui);

    auto& inter_work = *interProp.workVectors;
    *inter_work[siconos::integrators::EulerMoreauOSI::YOLD] = *(inter.y(0));
    *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD] = *(inter.lambda(0));
  }
}

/// @cond

// Visitor is not required for EulerMoreauOSI ?
// struct siconos::integrators::EulerMoreauOSI::_NSLEffectOnFreeOutput : public
// SiconosVisitor
// {
//   using SiconosVisitor::visit;

//   siconos::nonsmooth_formulations::OneStepNSProblem* _osnsp{nullptr};
//   std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};

//   _NSLEffectOnFreeOutput(siconos::nonsmooth_formulations::OneStepNSProblem* p,
//                          std::shared_ptr<siconos::modeling::Interaction> inter)
//       : _osnsp(p), _inter(inter){};

//   void visit(const EqualityConditionNSL& nslaw) const override { ; }
//   void visit(const MixedComplementarityConditionNSL& nslaw) const override { ; }
// };

/// @endcond

void siconos::integrators::EulerMoreauOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeFreeOutput(...)\n");
  auto allOSNS = _simulation->oneStepNSProblems();
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);

  auto& DSlink = inter->linkToDSVariables();
  auto& relationVec = inter->relationVectors();

  auto& inter_work = *indexSet->properties(vertex_inter).workVectors;
  auto& inter_work_block = *(indexSet->properties(vertex_inter)).workBlockVectors;
  // Get relation and non smooth law types
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();

  auto sizeY = static_cast<std::size_t>(inter->nonSmoothLaw()->size());

  std::shared_ptr<siconos::algebra::SiconosVector> H_alpha{nullptr};

  auto deltax = inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
  DEBUG_EXPR(deltax->display(););
  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::EulerMoreauOSI::OSNSP_RHS];

  auto Xfree = inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE];
  DEBUG_EXPR(Xfree->display(););
  assert(Xfree);

  auto mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());
  assert(relationType == siconos::modeling::RelationType::FirstOrder);
  auto forel =
      std::static_pointer_cast<siconos::modeling::FirstOrderR>(mainInteraction->relation());
  if (relationSubType == siconos::modeling::RelationSubType::Type2R ||
      relationSubType == siconos::modeling::RelationSubType::NonLinearR) {
    DEBUG_PRINT(
        "relationType == siconos::modeling::RelationType::FirstOrder && (relationSubType == "
        "Type2R || relationSubType == "
        "NonLinearR)\n")
    auto& lambda = *inter->lambda(0);
    if (forel->hasJacobianhOver_lambda()) {
      auto D = forel->jacobianhOver_lambda();  // read only view
      osnsp_rhs = D * lambda;
      osnsp_rhs *= -1.0;
    } else {
      osnsp_rhs.setZero();
    }

    if (forel->hasJacobianhOver_state()) {
      auto C = forel->jacobianhOver_state();  // read-only view
      siconos::algebra::matrixBlockVector_prod(C, *deltax, osnsp_rhs, false);
    }

    if (_useGammaForRelation) {
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::ComputeFreeOutput not yet implemented with "
          "useGammaForRelation() for FirstOrderR and Type2R and H_alpha->getValue() should "
          "return the mid-point value");
    }
    auto& hAlpha = *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
    DEBUG_EXPR(hAlpha.display());
    osnsp_rhs += hAlpha;
    DEBUG_EXPR(osnsp_rhs.display(););
  } else if (relationSubType == siconos::modeling::RelationSubType::Type1R) {
    DEBUG_PRINT(
        "relationType == siconos::modeling::RelationType::FirstOrder && relationSubType == "
        "Type1R\n");
    assert(Xfree);
    assert(deltax);

    if (forel->hasJacobianhOver_state()) {
      auto C = forel->jacobianhOver_state();  // read-only view
      siconos::algebra::matrixBlockVector_prod(C, *Xfree, osnsp_rhs, true);
    } else
      osnsp_rhs.setZero();

    if (_useGammaForRelation) {
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::ComputeFreeOutput not yet implemented with "
          "useGammaForRelation() for FirstOrderR and Typ2R and H_alpha->getValue() should "
          "return the mid-point value");
    }
    if (inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA]) {
      osnsp_rhs += *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
    }
  } else  // First Order Linear Relation
  {
    DEBUG_PRINT("relationType == siconos::modeling::RelationType::FirstOrder\n");
    if (forel->hasJacobianhOver_state()) {
      assert(Xfree);
      assert(deltax);
      auto C = forel->jacobianhOver_state();
      if (_useGammaForRelation) {
        siconos::algebra::matrixBlockVector_prod(C, *deltax, osnsp_rhs, true);
      } else {
        siconos::algebra::matrixBlockVector_prod(C, *Xfree, osnsp_rhs, true);
      }
    } else
      osnsp_rhs.setZero();

    DEBUG_EXPR(osnsp_rhs.display(););
    //    if (relationSubType == siconos::modeling::RelationSubType::LinearTIR ||
    //  relationSubType == siconos::modeling::RelationSubType::LinearR) {
    // In the first order linear case it may be required to add e to y.
    // y = CXfree + e
    if (relationSubType == siconos::modeling::RelationSubType::LinearTIR) {
      auto linrel = std::static_pointer_cast<siconos::modeling::FirstOrderLinearTIR>(forel);
      if (linrel->haseVector()) {
        auto e = linrel->eVector();
        osnsp_rhs += e;
      }
    } else {
      auto linrel = std::static_pointer_cast<siconos::modeling::FirstOrderLinearR>(forel);
      if (linrel->haseVector()) {
        auto e = linrel->eVector();
        osnsp_rhs += e;
      }
    }
    DEBUG_EXPR(osnsp_rhs.display(););
    //}
  }
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeFreeOutput(...)\n");
}

void siconos::integrators::EulerMoreauOSI::integrate(double& tinit, double& tend, double& tout,
                                                     int&) {
  // Last parameter is not used (required for LsodarOSI but not for EulerMoreauOSI).

  // double h = tend - tinit;
  tout = tend;

  std::shared_ptr<siconos::algebra::SiconosMatrix> W{nullptr};
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    // auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::integrate - not yet implemented for any "
        "kind "
        "of dynamical system.");
  }
}

void siconos::integrators::EulerMoreauOSI::updateState(const unsigned int) {
  DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState\n");

  double h = _simulation->timeStep();

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    // Get the DS type
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    auto fods = std::static_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds);
    auto& x = *ds->x();
    DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState Old value of x\n");
    DEBUG_EXPR(x.display());
    DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState residu value\n");
    DEBUG_EXPR(d->r()->display());

    // TODO ???
    bool baux = (useRCC &&
                 (not std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds)) &&
                 _simulation->relativeConvergenceCriterionHeld());

    //    std::shared_ptr<siconos::algebra::SiconosVector> xFree = d->xFree();

    // Save value of q in local_buffer for relative convergence computation
    if (baux) *ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER] = x;

    if (_useGamma) {
      // XXX UseGamma broken ? -- xhub
      x = h * _gamma * *fods->r();  // x = gamma*h*r
    } else {
      x = h * *fods->r();
    }

    x = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(x);  // x = h* W^{-1} *r

    x += *ds_work_vectors[siconos::integrators::EulerMoreauOSI::FREE];  // x+=xfree

    if (baux) {
      auto ds_norm_ref = 1. + fods->x0().norm();  // Should we save this in the graph?
      *ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER] -= x;
      auto aux =
          (ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER]->norm2()) /
          (ds_norm_ref);
      if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
    }
    DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState New value of x\n");
    DEBUG_EXPR(x.display());
  }
}

void siconos::integrators::EulerMoreauOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== EulerMoreauOSI OSI display ======"
            << "\n";

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--------------------------------\n";
    std::cout << "--> W of dynamical system number " << ds->number() << ": "
              << "\n";
    if (_dynamicalSystemsGraph->properties(*dsi).iterationMatrix)
      _dynamicalSystemsGraph->properties(*dsi).iterationMatrix->display();
    else
      std::cout << "-> nullptr\n";
    std::cout << "--> and corresponding theta is: " << _theta << "\n";
  }
  std::cout << "================================\n";
}
void siconos::integrators::EulerMoreauOSI::updateOutput(double time) {
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi
   */
  for (auto level = _levelMinForOutput; level < _levelMaxForOutput + 1; level++)
    updateOutput(time, level);
}

void siconos::integrators::EulerMoreauOSI::updateInput(double time) {
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi
   */
  for (auto level = _levelMinForInput; level < _levelMaxForInput + 1; level++)
    updateInput(time, level);
}

void siconos::integrators::EulerMoreauOSI::updateOutput(double time, unsigned int level) {
  DEBUG_BEGIN(
      "siconos::integrators::EulerMoreauOSI::updateOutput(double time, unsigned int "
      "level)\n");
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  //_simulation->nonSmoothDynamicalSystem()->updateOutput(time,level);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet0->bundle(*ui);
    assert(inter.lowerLevelForOutput() <= level);
    assert(inter.upperLevelForOutput() >= level);

    auto& DSlink = inter.linkToDSVariables();

    auto& interProp = indexSet0->properties(*ui);
    auto& inter_work = *interProp.workVectors;
    auto& inter_work_block = *interProp.workBlockVectors;
    auto relationSubType = inter.relation()->getSubType();
    if (relationSubType == siconos::modeling::RelationSubType::Type2R ||
        relationSubType == siconos::modeling::RelationSubType::NonLinearR) {
      auto& rel = static_cast<siconos::modeling::FirstOrderNonLinearR&>(*inter.relation());
      // compute the new y obtained by linearisation (see DevNotes)
      // y_{alpha+1}_{k+1} = h(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha}
      //                     )
      // or equivalently
      // y_{alpha+1}_{k+1} = y_{alpha}_{k+1} - ResiduY_{k+1}^{alpha}
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha}
      //                     )
      auto& y = *inter.y(level);
      DEBUG_EXPR(y.display());

      auto& residuY = *inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y];
      DEBUG_EXPR(residuY.display());
      y = *inter_work[siconos::integrators::EulerMoreauOSI::YOLD] - residuY;

      if (rel.hasJacobianhOver_lambda())
        y -= rel.jacobianhOver_lambda() *
             *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD];

      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : y(level) \n");
      DEBUG_EXPR(y.display());

      auto& deltax = *inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (rel.hasJacobianhOver_state()) {
        auto C = rel.jacobianhOver_state();
        siconos::algebra::matrixBlockVector_prod(C, deltax, y, false);
      }

      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : y before osnsM\n");
      DEBUG_EXPR(y.display());
      if (interProp.block) {
        y += *interProp.block * *inter.lambda(level);
        // block = h* C* W ^ -1 * B + D
        // y += osnsM * *inter.lambda(level);
        DEBUG_EXPR(inter.lambda(level)->display());
        DEBUG_EXPR(osnsM.display());
        DEBUG_PRINT(
            "siconos::integrators::EulerMoreauOSI::updateOutput : new linearized y \n");
        DEBUG_EXPR(y.display());
      }

      auto& hAlpha = *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];

      rel.computeh(*DSlink[siconos::modeling::FirstOrderR::Xxx], time, *inter.lambda(level),
                   hAlpha);
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : new Halpha \n");
      DEBUG_EXPR(hAlpha.display());
    } else
      inter.computeOutput(time, level);
  }
  DEBUG_END(
      "siconos::integrators::EulerMoreauOSI::updateOutput(double time, unsigned int "
      "level)\n");
}

void siconos::integrators::EulerMoreauOSI::updateInput(double time, unsigned int level) {
  // Set dynamical systems non-smooth part to zero.
  // Warning: This reset may be prone to issue with multiple osis.
  // _simulation->nonSmoothDynamicalSystem()->resetNonSmoothPart(level);

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;

  auto& indexSet0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui) {
    if (!checkInteractionOSI(indexSet0, ui)) continue;
    auto& inter = *indexSet0.bundle(*ui);

    assert(inter.lowerLevelForInput() <= level);
    assert(inter.upperLevelForInput() >= level);

    auto& DSlink = inter.linkToDSVariables();

    auto& interProp = indexSet0.properties(*ui);
    auto& inter_work = *interProp.workVectors;
    auto& inter_work_mat = *interProp.workMatrices;
    auto& inter_work_block = *interProp.workBlockVectors;

    auto relationSubType = inter.relation()->getSubType();
    if (relationSubType == siconos::modeling::RelationSubType::Type2R) {
      auto& r = static_cast<siconos::modeling::FirstOrderType2R&>(*inter.relation());
      auto lambda = *inter.lambda(level);
      lambda -= *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD];

      if (r.hasJacobianhOver_lambda())
        *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA] +=
            r.jacobianhOver_lambda() * lambda;
      *DSlink[siconos::modeling::FirstOrderR::Rrr] +=
          *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA];
      DEBUG_EXPR(DSlink[siconos::modeling::FirstOrderR::Rrr]->display(););
      // compute the new g_alpha

      r.computeg(*inter.lambda(level),
                 *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]);
      DEBUG_EXPR(inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]->display(););
    } else if (relationSubType == siconos::modeling::RelationSubType::NonLinearR) {
      auto forel =
          std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearR>(inter.relation());
      // compute the new r  obtained by linearisation
      // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
      //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}-
      //                     lambda_{k+1}^{alpha}
      //                     )

      auto lambda = *inter.lambda(level);
      lambda -= *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD];

      // Remind that g_alpha has only one block
      auto& g_alpha =
          *(*inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]).vector(0);

      if (forel->hasJacobiangOver_lambda()) g_alpha += forel->jacobiangOver_lambda() * lambda;

      auto& deltax = *inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("siconos::modeling::FirstOrderNonLinearR::computeInput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (forel->hasJacobiangOver_state())
        siconos::algebra::matrixBlockVector_prod(forel->jacobiangOver_lambda(), deltax,
                                                 g_alpha, false);
      // Khat = h * K * W^-1 * B
      g_alpha += *inter_work_mat[siconos::integrators::EulerMoreauOSI::MAT_KHAT] *
                 *inter.lambda(level);

      *DSlink[siconos::modeling::FirstOrderR::Rrr] += g_alpha;

      // compute the new g_alpha
      forel->computeg(*DSlink[siconos::modeling::FirstOrderR::Xxx], time, *inter.lambda(level),
                      *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]);
    } else {
      inter.computeInput(time, level);
    }
  }
}

double siconos::integrators::EulerMoreauOSI::computeResiduOutput(
    double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet) {
  double residu = 0.0;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto& inter_work = *indexSet->properties(*ui).workVectors;
    auto& residuY = *inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y];
    auto& inter = *indexSet->bundle(*ui);
    residuY = *(inter.y(0)) - *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
    residu = std::max(residu, residuY.norm2());
  }
  return residu;
}
double siconos::integrators::EulerMoreauOSI::computeResiduInput(
    double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet) {
  double residu = 0.0;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto& interProp = indexSet->properties(*ui);
    auto& inter_work = *interProp.workVectors;
    auto& inter_work_block = *interProp.workBlockVectors;
    auto inter = indexSet->bundle(*ui);
    auto& DSlink = inter->linkToDSVariables();
    auto& residuR = *inter_work[siconos::integrators::EulerMoreauOSI::WORK_DS];
    // Residu_r = r_alpha_k+1 - g_alpha;

    auto r = DSlink[siconos::modeling::FirstOrderR::Rrr];
    auto galpha = inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA];
    double norm2 = 0;
    assert(r->size() == galpha->size());

    // Compute euclidian norm of the difference between r and galpha
    for (size_t i = 0; i < r->size(); ++i) {
      const siconos::algebra::SiconosVector& v = *r->vector(i);
      const siconos::algebra::SiconosVector& w = *galpha->vector(i);
      norm2 += (v - w).squaredNorm();
    }
    norm2 = std::sqrt(norm2);
    DEBUG_EXPR(norm2.display(););
    residu = std::max(residu, norm2);
  }
  return residu;
}
