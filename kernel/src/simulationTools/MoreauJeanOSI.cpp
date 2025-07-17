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
#include "MoreauJeanOSI.hpp"

#include "BlockVector.hpp"
#include "BoundaryCondition.hpp"
#include "DynamicalSystem.hpp"
#include "FremondImpactFrictionNSL.hpp"
#include "Interaction.hpp"
#include "LagrangianCompliantLinearTIR.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianRheonomousR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactFrictionNSL.hpp"         // for nslaw visitor
#include "NewtonImpactNSL.hpp"                 // for nslaw visitor
#include "NewtonImpactRollingFrictionNSL.hpp"  // for nslaw visitor
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "RotationQuaternion.hpp"  // for quaternionFromTwistVector and compositionLawLieGroup
#include "SiconosAlgebraAddons.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosPointers.hpp"  // For createSPtr
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::_NSLEffectOnFreeOutput(
    siconos::nonsmooth_formulations::OneStepNSProblem &p,
    siconos::modeling::Interaction &inter, siconos::graphs::InteractionProperties &interProp,
    double theta)
    : _osnsp{p}, _inter{inter}, _interProp{interProp}, _theta{theta} {};

void siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::NewtonImpactNSL &nslaw) {
  auto e = nslaw.e();
  auto &osnsp_rhs = *(*_interProp.workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];
  osnsp_rhs += e * _inter.y_k(_osnsp.inputOutputLevel());
}

void siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::NewtonImpactFrictionNSL &nslaw) {
  auto &osnsp_rhs = *(*_interProp.workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];

  // The normal part is multiplied depends on en
  if (nslaw.en() > 0.0) {
    osnsp_rhs(0) += nslaw.en() * _inter.y_k(_osnsp.inputOutputLevel())(0);
  }
  // The tangential part is multiplied depends on et
  if (nslaw.et() > 0.0) {
    osnsp_rhs(1) += nslaw.et() * _inter.y_k(_osnsp.inputOutputLevel())(1);
    if (_inter.nonSmoothLaw()->size() > 2) {
      osnsp_rhs(2) += nslaw.et() * _inter.y_k(_osnsp.inputOutputLevel())(2);
    }
  }
}

void siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::FremondImpactFrictionNSL &nslaw) {
  auto &osnsp_rhs = *(*_interProp.workVectors)[MoreauJeanOSI::OSNSP_RHS];

  // compute the local tangential velocity at t_{k+theta}
  osnsp_rhs = _theta * osnsp_rhs + (1. - _theta) * _inter.y_k(_osnsp.inputOutputLevel());

  // The normal part is multiplied depends on en
  osnsp_rhs(0) += (_theta * (1. + nslaw.en()) - 1.) * _inter.y_k(_osnsp.inputOutputLevel())(0);

  // The tangential part is multiplied depends on et

  if (nslaw.et() > 0.0) {
    THROW_EXCEPTION(
        "MoreauJeanOSI::computeFreeOutput : do  not know how to take into account a "
        "tangential coefficient of restitution with Fremon NSL");
    // osnsp_rhs(1) +=  nslaw.et()  * _inter.y_k(_osnsp.inputOutputLevel())(1);
    // if(_inter.nonSmoothLaw()->size()>2)
    // 	{
    // 	  osnsp_rhs(2) +=  nslaw.et()  * _inter.y_k(_osnsp.inputOutputLevel())(2);    }
  }
}

void siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::NewtonImpactRollingFrictionNSL &nslaw) {
  auto &osnsp_rhs = *(*_interProp.workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];

  // The normal part is multiplied depends on en
  if (nslaw.en() > 0.0) {
    osnsp_rhs(0) += nslaw.en() * _inter.y_k(_osnsp.inputOutputLevel())(0);
  }
  // The tangential part is multiplied depends on et
  if (nslaw.et() > 0.0) {
    osnsp_rhs(1) += nslaw.et() * _inter.y_k(_osnsp.inputOutputLevel())(1);
    if (_inter.nonSmoothLaw()->size() > 2) {
      osnsp_rhs(2) += nslaw.et() * _inter.y_k(_osnsp.inputOutputLevel())(2);
    }
  }
}

namespace siconos::integrators::moreau_jean {
template <typename MatrixType, typename LUMatrixType, typename LagrangianType>
void computeIterationMatrix_Lagrangian_T(double time, double time_step, double theta,
                                         LagrangianType &lds, MatrixType &iterationMatrix,
                                         std::shared_ptr<LUMatrixType> &luw) {
  if (lds.isLinear()) return;

  auto ndof = lds.dimension();
  if (lds.hasMass()) {
    lds.computeMass(lds.q_read());
    iterationMatrix = lds.mass();  // copy
  } else {
    if constexpr (std::is_base_of_v<Eigen::SparseMatrixBase<MatrixType>, MatrixType>) {
      iterationMatrix.resize(ndof, ndof);
      for (auto i = 0; i < ndof; ++i) iterationMatrix.insert(i, i) = 1.0;
      iterationMatrix.makeCompressed();
    } else {
      iterationMatrix.setIdentity();
    }
  }

  lds.computeJacobianTotalForcesOver_velocity(lds.velocity_read(), lds.q_read(), time);
  if (lds.hasJacobianTotalForcesOver_velocity()) {
    iterationMatrix -= time_step * theta * lds.jacobianTotalForcesOver_velocity();
  }

  lds.computeJacobianTotalForcesOver_q(lds.velocity_read(), lds.q_read(), time);
  if (lds.hasJacobianTotalForcesOver_q()) {
    iterationMatrix -= time_step * time_step * theta * theta * lds.jacobianTotalForcesOver_q();
  }

  // LU factorization
  luw = std::make_shared<LUMatrixType>();
  if constexpr (std::is_same_v<LUMatrixType, Eigen::SparseLU<MatrixType>>) {
    luw->analyzePattern(iterationMatrix);
    luw->factorize(iterationMatrix);
  } else {
    *luw = LUMatrixType(iterationMatrix);
  }
}
}  // namespace siconos::integrators::moreau_jean

// --- constructor from a set of data ---
siconos::integrators::MoreauJeanOSI::MoreauJeanOSI(double theta, double gamma)
    : OneStepIntegrator(IntegratorType::MOREAUJEANOSI, 1, 0, 1, 1, 1) {
  _theta = theta;
  if (!std::isnan(gamma)) {
    _gamma = gamma;
    _useGamma = true;
  } else {
    _gamma = 1.0 / 2.0;
    _useGamma = false;
  }
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::MoreauJeanOSI::IterationMatrixBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds &&
         "siconos::integrators::MoreauJeanOSI::IterationMatrixBoundaryConditions(ds): ds == "
         "nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
      .iterationMatrix;
}

void siconos::integrators::MoreauJeanOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::initializeWorkVectorsForDS(Model&, "
      "double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");

  // Check dynamical system type
  auto dsType = siconos::types::type_value(*ds);

  assert(dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
         dsType == siconos::modeling::Type::LagrangianDS ||
         dsType == siconos::modeling::Type::LagrangianSparseDS ||
         dsType == siconos::modeling::Type::NewtonEulerDS ||
         dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS);

  // Compute W (iteration matrix)
  auto sods = std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds);
  initializeIterationMatrix(t, sods);

  // Initialize work vectors
  auto &ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(siconos::integrators::MoreauJeanOSI::WORK_LENGTH);

  ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE] =
      std::make_shared<siconos::algebra::SiconosVector>(sods->dimension());
  ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE] =
      std::make_shared<siconos::algebra::SiconosVector>(sods->dimension());

  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    // buffers allocation (inside the graph)
    auto lds = std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds);
    ds_work_vectors[siconos::integrators::MoreauJeanOSI::BUFFER] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    // Update dynamical system components (for memory swap).
    lds->computeTotalForces(lds->velocity_read(), lds->q_read(), t);
    lds->swapInMemory();

  } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
    auto neds = std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds);
    // Compute a first value of the dotq  to store it in  _dotqMemory
    *neds->dotq() = neds->T() * neds->twist_read();
    neds->computeWrench(neds->twist_read(), neds->q_read(), t);
    neds->swapInMemory();
  }
  DEBUG_END("siconos::integrators::MoreauJeanOSI::initializeWorkVectorsForDS()\n");
}

void siconos::integrators::MoreauJeanOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction &inter, siconos::graphs::InteractionProperties &interProp,
    siconos::graphs::DynamicalSystemsGraph &DSG) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::"
      "initializeWorkVectorsForInteraction(Interaction "
      "&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds1 = interProp.source;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  DEBUG_PRINTF("interaction number %i\n", inter.number());

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::MoreauJeanOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::MoreauJeanOSI::BLOCK_WORK_LENGTH);
  }

  auto &inter_work = *interProp.workVectors;
  auto &inter_work_block = *interProp.workBlockVectors;

  auto &relation = *inter.relation();
  relation.checkSize(inter);

  if (!inter_work[siconos::integrators::MoreauJeanOSI::OSNSP_RHS])
    inter_work[siconos::integrators::MoreauJeanOSI::OSNSP_RHS] =
        std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  // Check if interations levels (i.e. y and lambda sizes) are compliant with
  // the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  auto xfree = siconos::integrators::MoreauJeanOSI::xfree;
  DEBUG_PRINTF("ds1->number() %i\n", ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n", ds2->number());

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!inter_work_block[xfree]) || (inter_work_block[xfree]->numberOfBlocks() != 2))
      inter_work_block[xfree] = std::make_shared<siconos::algebra::BlockVector>(2);
  } else {
    if ((!inter_work_block[xfree]) || (inter_work_block[xfree]->numberOfBlocks() != 1))
      inter_work_block[xfree] = std::make_shared<siconos::algebra::BlockVector>(1);
  }

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    auto &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    inter_work_block[xfree]->setVectorPtr(
        0, workVds1[siconos::integrators::MoreauJeanOSI::VFREE]);
  }

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      auto &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      inter_work_block[xfree]->setVectorPtr(
          1, workVds2[siconos::integrators::MoreauJeanOSI::VFREE]);
    }
  }

  DEBUG_EXPR(siconos::algebra::print(*inter_work_block[xfree]););
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::"
      "initializeWorkVectorsForInteraction(Interaction "
      "&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");
}

void siconos::integrators::MoreauJeanOSI::initialize_nonsmooth_problems() {
  auto allOSNS = _simulation->oneStepNSProblems();
  ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY])->setIndexSetLevel(1);
  ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY])->setInputOutputLevel(1);
  //  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->initialize(_simulation);
}

void siconos::integrators::MoreauJeanOSI::initializeIterationMatrix(
    double time, std::shared_ptr<siconos::modeling::SecondOrderDS> ds) {
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical
  // system, depending on its type.
  assert(ds);
  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrix(t,ds) "
        "- ds does not "
        "belong to the OSI.");

  const auto &dsv = _dynamicalSystemsGraph->descriptor(ds);

  assert(!_dynamicalSystemsGraph->properties(dsv).iterationMatrix);
  double timeStep = _simulation->timeStep();
  auto ndof = ds->dimension();

  // Allocate storage for W in the graph
  _dynamicalSystemsGraph->properties(dsv).iterationMatrix =
      std::make_shared<siconos::algebra::SiconosMatrix>(ndof, ndof);

  auto iterationMat = _dynamicalSystemsGraph->properties(dsv).iterationMatrix;
  double htheta = timeStep * _theta;
  double h2theta2 = timeStep * timeStep * _theta * _theta;

  // Linear Lagrangian sytems, where W is constant must be computed only once
  if (auto lldds =
          std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearDiagonalDS>(ds)) {
    // For these systems, W is a diagonal matrix. W^(-1) is straightforward and saved in W.
    if (lldds->hasMassMatrix())
      iterationMat->diagonal() = lldds->massMatrix();
    else
      iterationMat->setIdentity();

    if (lldds->hasDampingMatrix()) {
      iterationMat->diagonal() += htheta * lldds->dampingMatrix();
    }

    if (lldds->hasStiffnessMatrix()) {
      iterationMat->diagonal() += h2theta2 * lldds->stiffnessMatrix();
    }
    if (lldds->boundaryConditions()) _initializeIterationMatrixBoundaryConditions(*lldds, dsv);

    // W inverse in place
    iterationMat->diagonal() = iterationMat->diagonal().array().inverse();
  }

  // Linear Lagrangian time-invariant sytems, where W is constant must be computed only once
  else if (auto ltids =
               std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    if (ltids->hasMass())
      *iterationMat = ltids->mass();
    else
      iterationMat->setIdentity();

    if (ltids->hasDampingMatrix()) *iterationMat += htheta * ltids->dampingMatrix();
    if (ltids->hasStiffnessMatrix()) *iterationMat += h2theta2 * ltids->stiffnessMatrix();
    if (ltids->boundaryConditions()) _initializeIterationMatrixBoundaryConditions(*ltids, dsv);
    // LU factorization
    _dynamicalSystemsGraph->properties(dsv).LUW =
        std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*iterationMat);
  }
  // General Lagrangian DS
  else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    siconos::integrators::moreau_jean::computeIterationMatrix_Lagrangian_T(
        time, _simulation->timeStep(), _theta, *lds, *iterationMat,
        _dynamicalSystemsGraph->properties(dsv).LUW);
    if (lds->boundaryConditions()) _initializeIterationMatrixBoundaryConditions(*lds, dsv);
  }

  else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
    siconos::integrators::moreau_jean::computeIterationMatrix_NewtonEuler(
        time, _simulation->timeStep(), _theta, *neds, *iterationMat,
        _dynamicalSystemsGraph->properties(dsv).LUW);
    if (neds->boundaryConditions()) _initializeIterationMatrixBoundaryConditions(*neds, dsv);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrix - not "
        "yet implemented for Dynamical system of type : " +
        siconos::types::type_name(*ds));
  // Remark: except for LagrangianLinearDiagonalDS and LagrangianLinearTIDS, W is not
  // LU-factorized nor inversed here.

  // if (isWSymmetricDefinitePositive()) {
  //   auto &W = *_dynamicalSystemsGraph->properties(dsv).iterationMatrix;
  //   // W.setIsSymmetric(true);
  //   // W.setIsPositiveDefinite(true); // TODO : make sure it is not needed!
}

void siconos::integrators::MoreauJeanOSI::_initializeIterationMatrixBoundaryConditions(
    siconos::modeling::SecondOrderDS &ds,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor &dsv) {
  // This function:
  // - allocate memory for a matrix IterationMatrixBoundaryConditions
  // - insert this matrix into IterationMatrixBoundaryConditionsMap with ds as a key

  assert(!_dynamicalSystemsGraph->properties(dsv).iterationMatrixBoundaryConditions);

  auto dsType = siconos::types::type_value(ds);
  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::NewtonEulerDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    // Memory allocation for IterationMatrixBoundaryConditions
    auto sizeIterationMatrixBoundaryConditions =
        ds.dimension();  // n for first order systems, ndof for lagrangian.

    auto &d = static_cast<siconos::modeling::SecondOrderDS &>(ds);
    _dynamicalSystemsGraph->properties(dsv).iterationMatrixBoundaryConditions =
        std::make_shared<siconos::algebra::SiconosMatrix>(
            sizeIterationMatrixBoundaryConditions, d.boundaryConditions()->size());
    _computeIterationMatrixBoundaryConditions(
        ds, *_dynamicalSystemsGraph->properties(dsv).iterationMatrixBoundaryConditions,
        *_dynamicalSystemsGraph->properties(dsv).iterationMatrix);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixBoundaryConditions - "
        "not yet implemented for Dynamical system of type :" +
        siconos::types::type_name(ds));
}

void siconos::integrators::MoreauJeanOSI::_computeIterationMatrixBoundaryConditions(
    siconos::modeling::SecondOrderDS &ds,
    siconos::algebra::SiconosMatrix &IterationMatrixBoundaryConditions,
    siconos::algebra::SiconosMatrix &iteration_matrix) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::_computeIterationMatrixBoundaryConditions\n");
  // Compute IterationMatrixBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, IterationMatrixBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixBoundaryConditions.
  auto dsType = siconos::types::type_value(ds);

  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::NewtonEulerDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    auto columntmp = std::make_shared<siconos::algebra::SiconosVector>(ds.dimension());

    int columnindex = 0;

    auto &d = static_cast<siconos::modeling::SecondOrderDS &>(ds);

    if (!siconos::algebra::isSymmetric(
            iteration_matrix,
            1e-10))  // Warning this operation could be quite expensive
    {
      // siconos::algebra::print(iteration_matrix);
      std::cout << "Warning, we apply boundary conditions assuming W symmetric" << std::endl;
    }

    for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
      *columntmp = iteration_matrix.col(itindex);
      /*\warning we assume that W is symmetric
        we store only the column and not the row */
      IterationMatrixBoundaryConditions.col(columnindex) = *columntmp;
      columnindex++;
    }

    columnindex = 0;
    for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
      double diag = iteration_matrix(itindex, itindex);
      columntmp->setZero();
      (*columntmp)(itindex) = diag;
      iteration_matrix.col(itindex) = *columntmp;
      iteration_matrix.row(itindex) = *columntmp;
      columnindex++;
    }
    DEBUG_EXPR(siconos::algebra::print(iteration_matrix););
    DEBUG_EXPR(siconos::algebra::print(IterationMatrixBoundaryConditions););
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::computeIterationMatrixBoundaryConditions - not "
        "yet "
        "implemented for Dynamical system type : " +
        siconos::types::type_name(ds));
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::_computeIterationMatrixBoundaryConditions\n");
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::MoreauJeanOSI::iterationMatrixInverse(
    std::shared_ptr<siconos::modeling::SecondOrderDS> ds,
    const Eigen::FullPivLU<siconos::algebra::SiconosMatrix> &LUW) {
  /* We compute and return the current inverse of W matrix */

  // WARNING FP: we assume that LUW is up to date and contains the LU factorization of W

  const auto &dsv = _dynamicalSystemsGraph->descriptor(ds);
  if (!_dynamicalSystemsGraph->properties(dsv).iterationMatrixInverse) {
    // Allocate iterationMatrixInverse at first call
    auto sizeW = ds->dimension();
    _dynamicalSystemsGraph->properties(dsv).iterationMatrixInverse =
        std::make_shared<siconos::algebra::SiconosMatrix>(sizeW, sizeW);
    *_dynamicalSystemsGraph->properties(dsv).iterationMatrixInverse = LUW.inverse();
  } else {
    auto dsType = siconos::types::type_value(*ds);
    if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
        dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      // nothing, Winverse is constant and has been created at first call
    } else {
      *_dynamicalSystemsGraph->properties(dsv).iterationMatrixInverse = LUW.inverse();
    }
  }
  return _dynamicalSystemsGraph->properties(dsv).iterationMatrixInverse;
}

void siconos::integrators::MoreauJeanOSI::computeInitialNewtonState() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeInitialNewtonState()\n");
  // Compute the position value giving the initial velocity.
  // The goal is to save one newton iteration for nearly linear system
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);

    if (_explicitNewtonEulerDSOperators) {
      if (siconos::types::type_value(ds) == siconos::modeling::Type::NewtonEulerDS) {
        // The goal is to update T() one time at the beginning of the Newton
        // Loop We want to be explicit on this function since we do not compute
        // their Jacobians.
        auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);
        auto qold = d.qMemory().getSiconosVector(0);
        d.computeT(qold);
      }
    }
    // The goal is to converge in one iteration if the system is almost linear
    // we start the Newton loop q = q0+hv0
    moreau_jean::updatePosition(_simulation->timeStep(), _theta, ds);
  }
  DEBUG_END("siconos::integrators::MoreauJeanOSI::computeInitialNewtonState()\n");
}

void siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(
    siconos::modeling::SecondOrderDS &d, siconos::algebra::SiconosVector &residu,
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, double t,
    const siconos::algebra::SiconosVector &v) {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(...)\n");
  if (d.boundaryConditions()) {
    d.boundaryConditions()->computePrescribedVelocity(t);

    unsigned int columnindex = 0;
    auto &IterationMatrixBoundaryConditions =
        *_dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions;
    auto columntmp = std::make_shared<siconos::algebra::SiconosVector>(d.dimension());

    auto prescribedVelocity = d.boundaryConditions()->prescribedVelocity();
    for (const auto &itindex : d.boundaryConditions()->velocityIndices()) {
      double DeltaPrescribedVelocity = prescribedVelocity(columnindex) - v(itindex);
      DEBUG_PRINTF("index  = %i, value = %e\n", *itindex, prescribedVelocity(columnindex));
      DEBUG_PRINTF("DeltaPrescribedVelocity = %e\n", DeltaPrescribedVelocity);
      *columntmp = IterationMatrixBoundaryConditions.col(columnindex);
      residu -= *columntmp * (DeltaPrescribedVelocity);

      residu(itindex) = -(*columntmp)(itindex) * (DeltaPrescribedVelocity);

      columnindex++;
    }
  }
  DEBUG_END("siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(...)\n");
}

double siconos::integrators::MoreauJeanOSI::computeResidu() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeResidu()\n");
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
  //
  // std::shared_ptr<siconos::modeling::DynamicalSystem> ds; // Current
  // Dynamical System.
  siconos::modeling::Type dsType;  // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);
    auto &ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = siconos::types::type_value(ds);  // Its type

    // 3 - Lagrangian Non Linear Systems
    if (dsType == siconos::modeling::Type::LagrangianDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1)
      // - h*(1-theta)*forces(ti,vi,qi) - p_i+1
      auto &residuFree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::modeling::LagrangianDS &>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      const auto &vold = d.velocityMemory().getSiconosVector(0);

      const auto &v = *d.velocity();  // v = v_k,i+1
      // residuFree.setZero();
      DEBUG_EXPR(siconos::algebra::print(residuFree));

      DEBUG_EXPR(siconos::algebra::print(vold));
      DEBUG_EXPR(siconos::algebra::print(v));

      residuFree = v - vold;
      if (d.hasMass()) {
        d.computeMass(*d.q());
        residuFree = d.mass() * residuFree;  // residuFree = M(v - vold)
      }

      if (d.hasTotalForces()) {
        // Cheaper version: get forces(ti,vi,qi) from memory
        residuFree -= time_step * (1. - _theta) * d.forcesMemory().getSiconosVector(0);

        // Expensive computes forces(ti,vi,qi)
        // auto &qold = *d.qMemory()->getSiconosVector(0);
        // auto &vold = *d.velocityMemory()->getSiconosVector(0);

        // d.computeTotalForces(vold, qold, told);
        // double coef = -h * (1 - _theta);
        // // residuFree += coef * fL_i
        // siconos::algebra::scal(coef, *d.totalForces(), residuFree, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        d.computeTotalForces(d.velocity_read(), d.q_read(), t);
        residuFree -= time_step * _theta * d.totalForces();
        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        // std::shared_ptr<siconos::algebra::SiconosVector> qbasedonv(new
        // SiconosVector(*qold)); *qbasedonv +=  time_step * ((1 - _theta)* *vold +
        // _theta * *v); d.computeTotalForces(v, qbasedonv, t); coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        // siconos::algebra::scal(coef, *d.totalForces(), *residuFree, false);
      }

      applyBoundaryConditions(d, residuFree, dsi, t, v);

      free = residuFree;  // copy residuFree into Workfree
      DEBUG_EXPR(siconos::algebra::print(residuFree));

      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!

      applyBoundaryConditions(d, free, dsi, t, v);

      DEBUG_EXPR(siconos::algebra::print(free));
      normResidu = free.norm();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::LagrangianLinearTIDS\n");
      // ResiduFree = time_step*C*v_i + time_step*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual
      // for v = v_i otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + time_step*((1-theta)*(C v_i + K q_i) +theta * ( C*v
      // + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      auto &lltids = static_cast<siconos::modeling::LagrangianLinearTIDS &>(ds);

      auto &residuFree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // --- ResiduFree computation Equation (1) ---
      residuFree.setZero();
      // -- No need to update W --
      const auto &vold = lltids.velocityMemory().getSiconosVector(0);  // vi
      if (lltids.hasDampingMatrix()) residuFree += time_step * lltids.dampingMatrix() * vold;
      if (lltids.hasStiffnessMatrix()) {
        auto coeff = time_step * time_step * _theta;
        residuFree +=
            coeff * lltids.stiffnessMatrix() * vold;  // vfree += time_step^2*_theta*K*vi
        residuFree += time_step * lltids.stiffnessMatrix() *
                      lltids.qMemory().getSiconosVector(0);  // vfree += time_step*K*qi
      }

      if (lltids.hasExternalForces()) {
        // computes Fext(ti)
        lltids.computeFext(told);
        auto coeff = time_step * (1 - _theta);
        // vfree -= time_step*(1-_theta) * fext(ti)
        residuFree -= coeff * lltids.fext();
        // computes Fext(ti+1)
        lltids.computeFext(t);
        coeff = time_step * _theta;
        // vfree -= time_step*_theta * fext(ti+1)
        residuFree -= coeff * lltids.fext();
      }

      applyBoundaryConditions(lltids, residuFree, dsi, t, vold);

      free = residuFree;                          // copy residuFree into free
      if (lltids.p(1)) free -= lltids.p_read(1);  // Compute Residu in Workfree Notation !!
      // We use free as tmp buffer
      DEBUG_EXPR(siconos::algebra::print(free));
      DEBUG_EXPR(siconos::algebra::print(residuFree));

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm();
    }

    else if (dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      // ResiduFree = time_step*C*v_i + time_step*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual
      // for v = v_i otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + time_step*((1-theta)*(C v_i + K q_i) +theta * ( C*v
      // + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::modeling::LagrangianLinearDiagonalDS &>(ds);

      auto &residuFree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      const auto &qold = d.qMemory().getSiconosVector(0);         // qi
      const auto &vold = d.velocityMemory().getSiconosVector(0);  // vi
      // --- ResiduFree computation Equation (1) ---
      residuFree.setZero();
      // -- No need to update W --
      if (d.hasDampingMatrix()) {
        residuFree += time_step * d.dampingMatrix().cwiseProduct(vold);
      }
      if (d.hasStiffnessMatrix()) {
        auto coeff = time_step * time_step * _theta;
        residuFree += coeff * d.stiffnessMatrix().cwiseProduct(vold);
        residuFree += time_step * d.stiffnessMatrix().cwiseProduct(qold);
      }

      if (d.hasExternalForces()) {
        // computes Fext(ti)
        d.computeFext(told);
        auto coeff = time_step * (1 - _theta);
        residuFree -= coeff * d.fext();  // vfree -= time_step*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d.computeFext(t);
        residuFree -= time_step * _theta * d.fext();  // vfree -= time_step*_theta * fext(ti+1)
      }

      applyBoundaryConditions(d, residuFree, dsi, t, vold);

      free = residuFree;                // copy residuFree into free
      if (d.p(1)) free -= d.p_read(1);  // Compute Residu in Workfree Notation !!

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm();
    }

    else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - time_step*_theta*forces(t,v_k,i+1, q_k,i+1) -
      // time_step*(1-_theta)*forces(ti,vi,qi) - pi+1

      auto &residuFree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);

      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"
      const auto &vold = d.twistMemory().getSiconosVector(0);

      residuFree = d.totalInertiaMatrix() * (d.twist_read() - vold);

      // Cheaper version: get forces(ti,vi,qi) from memory
      residuFree -= time_step * (1. - _theta) * d.forcesMemory().getSiconosVector(0);

      // Expensive version to check ...
      // std::shared_ptr<siconos::algebra::SiconosVector> qold =
      // d.qMemory()->getSiconosVector(0);
      // std::shared_ptr<siconos::algebra::SiconosVector> vold =
      // d.twistMemory()->getSiconosVector(0);
      //  d.computeTotalForces(vold,qold,told);
      //  DEBUG_EXPR(siconos::algebra::print(*d.totalForces()););
      // double coef = -h * (1.0 - _theta);
      // siconos::algebra::scal(coef, *d.totalForces(), *residuFree, false);

      DEBUG_EXPR(siconos::algebra::print(residuFree););

      // computes forces(ti,v,q)
      d.computeWrench(d.twist_read(), d.q_read(), t);
      residuFree -= time_step * _theta * d.wrench();
      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI:: new forces :\n");
      DEBUG_EXPR(siconos::algebra::print(*d.totalForces()););
      DEBUG_EXPR(siconos::algebra::print(residuFree););

      applyBoundaryConditions(d, residuFree, dsi, t, d.twist_read());

      free = residuFree;

      if (d.p(1)) free -= d.p_read(1);

      applyBoundaryConditions(d, free, dsi, t, d.twist_read());

      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::computeResidu :\n");
      DEBUG_EXPR(siconos::algebra::print(residuFree););
      DEBUG_EXPR(if (d.p(1)) siconos::algebra::print(*d.p(1)););
      DEBUG_EXPR(siconos::algebra::print(free););

      normResidu = free.norm();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanOSI::computeResidu - not yet "
          "implemented for "
          "Dynamical system of type: " +
          siconos::types::type_name(ds));

    if (normResidu > maxResidu) maxResidu = normResidu;
  }
  DEBUG_END("siconos::integrators::MoreauJeanOSI::computeResidu()\n");
  return maxResidu;
}

void siconos::integrators::MoreauJeanOSI::computeFreeState() {
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeFreeState()\n");
  // This function computes "free" states of the DS belonging to this
  // Integrator. "Free" means without taking non-smooth effects into account.

  double nextTime = _simulation->nextTime();  // End of the time step
  double timeStep = _simulation->timeStep();

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
    auto iterationMatrix = _dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    // Get workVectors (rfree, vfree ...) from the graph
    auto &ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    auto &residuFree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
    auto &vfree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

    // Current dynamical system
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    // ========= Lagrangian Systems ===========
    if (auto lldds =
            std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearDiagonalDS>(ds)) {
      vfree = residuFree;  // initialize vfree with the content of residuFree
                           // W is diagonal and contains the inverse of the iteration matrix!
      vfree -= iterationMatrix->asDiagonal() * vfree;
      vfree += lldds->velocityMemory().getSiconosVector(0);
    } else if (auto lltids =
                   std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
      vfree = residuFree;  // initialize vfree with the content of residuFree
      // -- vfree =  vold - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve in place WX = vfree
      vfree = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(
          vfree);  // LUW is already up to date
      vfree *= -1.0;
      vfree += lltids->velocityMemory().getSiconosVector(0);
    } else if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // --- ResiduFree computation ---
      // ResFree = M(v-vold) - time_step*[theta*forces(t) + (1-theta)*forces(told)]
      //
      vfree = residuFree;  // initialize vfree with the content of residuFree
      // -- Update W --
      // Note: during computeIterationMatrix, mass and jacobians of forces will be computed/
      moreau_jean::computeIterationMatrix_Lagrangian_T(
          nextTime, timeStep, _theta, *lds, *iterationMatrix,
          _dynamicalSystemsGraph->properties(*dsi).LUW);
      if (lds->boundaryConditions()) {
        _computeIterationMatrixBoundaryConditions(
            *lds, *_dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions,
            *iterationMatrix);
      }
      vfree = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(vfree);
      vfree *= -1.0;
      vfree += lds->velocity_read();  // v = v_k,i+1;
    }

    // ========= NewtonEuler Systems ===========
    else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      vfree = residuFree;
      // -- Update W --
      // Note: during computeIterationMatrix, mass and jacobians of forces will be computed/
      moreau_jean::computeIterationMatrix_NewtonEuler(
          nextTime, timeStep, _theta, *neds, *iterationMatrix,
          _dynamicalSystemsGraph->properties(*dsi).LUW);
      if (neds->boundaryConditions()) {
        _computeIterationMatrixBoundaryConditions(
            *neds, *_dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions,
            *iterationMatrix);
      }
      // -- vfree =  twist - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      vfree = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(vfree);
      vfree *= -1.0;
      vfree += neds->twist_read();
    } else {
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanOSI::initialize - Only implemented for "
          "Lagrangian or NewtonEuler dynamical systems.");
    }
  }
}

void siconos::integrators::MoreauJeanOSI::prepareNewtonIteration(double time) {
  DEBUG_BEGIN(
      " siconos::integrators::MoreauJeanOSI::prepareNewtonIteration(double "
      "time)\n");
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      moreau_jean::computeIterationMatrix_Lagrangian_T(
          time, _simulation->timeStep(), _theta, *lds,
          *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix,
          _dynamicalSystemsGraph->properties(*dsi).LUW);

    } else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      moreau_jean::computeIterationMatrix_NewtonEuler(
          time, _simulation->timeStep(), _theta, *neds,
          *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix,
          _dynamicalSystemsGraph->properties(*dsi).LUW);
    }
  }

  if (!_explicitNewtonEulerDSOperators) {
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;

      std::shared_ptr<siconos::modeling::DynamicalSystem> ds =
          _dynamicalSystemsGraph->bundle(*dsi);

      //  VA <2016-04-19 Tue> We compute T to be consistent with the Jacobian
      //   at the beginning of the Newton iteration and not at the end
      auto dsType = siconos::types::type_value(*ds);
      if (dsType == siconos::modeling::Type::NewtonEulerDS) {
        auto d = std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds);
        d->computeT(d->q_read());
      }
    }
  }
  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);
  }

  DEBUG_END(
      " siconos::integrators::MoreauJeanOSI::prepareNewtonIteration(double "
      "time)\n");
}

void siconos::integrators::MoreauJeanOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor &vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem *osnsp) {
  /** \warning: ensures that it can also work with two different osi for two
   * different ds ?
   */
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::computeFreeOutput("
      "InteractionsGraph::VDescriptor& "
      "vertex_inter, siconos::nonsmooth_formulations::OneStepNSProblem* "
      "osnsp)\n");
  auto allOSNS = _simulation->oneStepNSProblems();
  auto &indexSet = *osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  assert(indexSet.bundle(vertex_inter));

  auto &inter = *indexSet.bundle(vertex_inter);

  assert(inter.relation());
  auto relationType = inter.relation()->getType();
  assert(!(inter.relation()->getType() == siconos::modeling::RelationType::FirstOrder) &&
         "You can not use MoreauJeanOSI with first-order type relations");

  auto &inter_work_block = *indexSet.properties(vertex_inter).workBlockVectors;

  auto &osnsp_rhs = *(*indexSet.properties(vertex_inter)
                           .workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];

  auto xfree = inter_work_block[siconos::integrators::MoreauJeanOSI::xfree];
  assert(xfree);

  // 1 - product H Xfree{}
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    auto H = std::dynamic_pointer_cast<siconos::modeling::LagrangianR>(inter.relation())
                 ->jacobianhOver_q();
    siconos::algebra::matrixBlockVector_prod(H, *xfree, osnsp_rhs, true);
  } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    auto H = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerR>(inter.relation())
                 ->H_NE_prod_T();
    siconos::algebra::matrixBlockVector_prod(H, *xfree, osnsp_rhs, true);
  }

  auto relationSubType = inter.relation()->getSubType();

  // 2 -  compute additional terms for ScleronomousR and CompliantLinearTIR

  // Get relation and non smooth law types

  if ((relationType == siconos::modeling::RelationType::Lagrangian) &&
      (relationSubType != siconos::modeling::RelationSubType::ScleronomousR)) {
    auto &DSlink = inter.linkToDSVariables();

    // For the relation of type LagrangianRheonomousR
    if (relationSubType == siconos::modeling::RelationSubType::RheonomousR) {
      assert(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp);

      auto rheoR =
          std::static_pointer_cast<siconos::modeling::LagrangianRheonomousR>(inter.relation());
      rheoR->computehdot(*DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q0)],
                         simulation()->getTkp1());

      // siconos::algebra::SiconosMatrix ID =
      //     siconos::algebra::SiconosMatrix::Identity(sizeY, sizeY);
      // This should be optimized -- vacary
      // y += hDot
      // auto hDot = rheoR->hdot();
      // osnsp_rhs += ID * hDot;
      osnsp_rhs += rheoR->hdot();
    }
    if (relationSubType == siconos::modeling::RelationSubType::CompliantLinearTIR) {
      assert(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp);
      auto compR = std::static_pointer_cast<siconos::modeling::LagrangianCompliantLinearTIR>(
          inter.relation());

      auto C = compR->CMatrix();
      double h = _simulation->timeStep();
      osnsp_rhs *= h * _theta;

      /* we have to check that the value are at the beginnning of the time
       * step */
      // + C q_k
      siconos::algebra::matrixBlockVector_prod(
          C, *DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q0)], osnsp_rhs,
          false);

      // + h(1-_theta)v_k

      *DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q1)] *= (1 - _theta) * h;
      siconos::algebra::matrixBlockVector_prod(
          C, *DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q1)], osnsp_rhs,
          false);

      if (compR->haseVector()) {
        osnsp_rhs += compR->eVector();
      }
    }
    DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););
  }

  // 3 - add part due to NonSmoothLaw
  _NSLEffectOnFreeOutput nslEffectOnFreeOutput(*osnsp, inter,
                                               indexSet.properties(vertex_inter), _theta);
  inter.nonSmoothLaw()->accept(nslEffectOnFreeOutput);

  DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););

  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::computeFreeOutput("
      "InteractionsGraph::VDescriptor& "
      "vertex_inter, siconos::nonsmooth_formulations::OneStepNSProblem* "
      "osnsp)\n");
}

void siconos::integrators::MoreauJeanOSI::integrate(double &tinit, double &tend, double &tout,
                                                    int &notUsed) {
  // Last parameter is not used (required for LsodarOSI but not for
  // MoreauJeanOSI).

  double h = tend - tinit;
  tout = tend;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    // get the ds
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto lltids = std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds);

    if (!lltids)
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanOSI::integrate - not yet implemented ds of type:" +
          siconos::types::type_name(*ds));

    // get q and velocity pointers for previous time step
    const auto &vold = lltids->velocityMemory().getSiconosVector(0);
    const auto &qold = lltids->qMemory().getSiconosVector(0);

    // velocity computation :
    //
    // v = vi + W^{-1}[ -h*C*vi - h*h*theta*K*vi - h*K*qi + h*theta*Fext(t) +
    // h*(1-theta)
    // * Fext(ti) ] + W^{-1}*pi+1
    //
    // get velocity pointers for current time step
    auto &v = *lltids->velocity();
    v = lltids->p_read(1);

    // -- No need to update W --
    if (lltids->hasDampingMatrix()) v -= h * lltids->dampingMatrix() * vold;
    if (lltids->hasStiffnessMatrix()) {
      auto coeff = h * h * _theta;
      v -= coeff * lltids->stiffnessMatrix() * vold;  // v += -h^2*theta*K*vi
      v -= h * lltids->stiffnessMatrix() * qold;
    }
    if (lltids->hasExternalForces()) {
      // computes Fext(ti)
      lltids->computeFext(tinit);
      auto coeff = h * (1 - _theta);
      v += coeff * lltids->fext();
      // computes Fext(ti+1)
      lltids->computeFext(tout);
      coeff = h * _theta;
      v += coeff * lltids->fext();
    }
    // -> Solve in place WX = v
    v = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(v);
    v += vold;
  }
}

void siconos::integrators::MoreauJeanOSI::updateState(const unsigned int) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::updateState(const unsigned int "
      ")\n");

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);

    auto &ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    auto &iterationMatrix = *_dynamicalSystemsGraph->properties(*dsi).iterationMatrix;
    // Get the DS type

    auto dsType = siconos::types::type_value(ds);

    // 3 - Lagrangian Systems
    if (dsType == siconos::modeling::Type::LagrangianDS ||
        dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
        dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::updateState(const unsigned int "
          "), dsType == "
          "siconos::modeling::Type::LagrangianDS || dsType == "
          "siconos::modeling::Type::LagrangianLinearTIDS \n");
      // get dynamical system
      auto &d = static_cast<siconos::modeling::LagrangianDS &>(ds);
      auto &vfree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      //    auto *vfree = d.velocityFree();
      // auto &v = *d.velocity();

      auto baux = dsType == siconos::modeling::Type::LagrangianDS && useRCC &&
                  _simulation->relativeConvergenceCriterionHeld();

      if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
        assert(((d.p(_levelMaxForInput)).get()) &&
               " siconos::integrators::MoreauJeanOSI::updateState() "
               "*d.p(_levelMaxForInput) "
               "== nullptr.");
        *d.velocity() = d.p_read(_levelMaxForInput);  // v = p
        if (d.boundaryConditions()) {
          for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
            (*d.velocity())(itindex) = 0.0;
          }
        }

        if (dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
          *d.velocity() = vfree;
          *d.velocity() += iterationMatrix.diagonal().cwiseProduct(d.velocity_read());
        } else {
          *d.velocity() =
              _dynamicalSystemsGraph->properties(*dsi).LUW->solve(d.velocity_read());
          *d.velocity() += vfree;
        }
      } else {
        *d.velocity() = vfree;
      }
      DEBUG_EXPR(siconos::algebra::print(v));

      if (d.boundaryConditions()) {
        int bc = 0;
        auto columntmp = std::make_shared<siconos::algebra::SiconosVector>(ds.dimension());

        for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
          *columntmp =
              _dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions->col(
                  bc);
          /*\warning we assume that W is symmetric in the Lagrangian case*/

          double value = -1. * columntmp->dot(d.velocity_read());
          if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
            value += (d.p_read(_levelMaxForInput))(itindex);
          }
          /* \warning the computation of reactionToBoundaryConditions take into
             account the contact impulse but not the external and internal
             forces. A complete computation of the residu should be better */
          (*d.reactionToBoundaryConditions())(bc) = value;
          bc++;
        }
      }

      auto &local_buffer = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::BUFFER];
      // Save value of q in stateTmp for future convergence computation
      if (baux) local_buffer = d.q_read();

      moreau_jean::updatePosition(_simulation->timeStep(), _theta, ds);

      if (baux) {
        double ds_norm_ref = 1. + ds.x0().norm();  // Should we save this in the graph?
        local_buffer -= d.q_read();
        double aux = (local_buffer.norm()) / ds_norm_ref;
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::updateState(), "
          "dsType == "
          "siconos::modeling::Type::NewtonEulerDS \n");

      // get dynamical system
      auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);
      // auto &v = *d.twist();
      // DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updateState()\n ")
      // DEBUG_EXPR(siconos::algebra::print(d));
      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updateState() prev v\n")

      // failure on bullet sims
      // d.p(_levelMaxForInput) is checked in next condition
      // assert(((d.p(_levelMaxForInput)).get()) &&
      //       " siconos::integrators::MoreauJeanOSI::updateState()
      //       *d.p(_levelMaxForInput)
      //       == nullptr.");

      auto &vfree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
        /*d.p has been fill by the Relation->computeInput, it contains
          B \lambda _{k+1}*/
        *d.twist() = d.p_read(_levelMaxForInput);  // v = p
        if (d.boundaryConditions())
          for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
            (*d.twist())(itindex) = 0.;
          }
        *d.twist() = _dynamicalSystemsGraph->properties(*dsi).LUW->solve(*d.twist());

        DEBUG_EXPR(std::cout << d.p_read(_levelMaxForInput));
        DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updatestate W CT lambda\n");
        DEBUG_EXPR(siconos::algebra::print(v));
        *d.twist() += vfree;
      } else
        *d.twist() = vfree;

      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updatestate work free\n");
      DEBUG_EXPR(siconos::algebra::print(vfree));
      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updatestate new v\n");
      DEBUG_EXPR(siconos::algebra::print(v));

      if (d.boundaryConditions()) {
        int bc = 0;
        auto columntmp = std::make_shared<siconos::algebra::SiconosVector>(ds.dimension());

        for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
          *columntmp =
              _dynamicalSystemsGraph->properties(*dsi).iterationMatrixBoundaryConditions->col(
                  bc);
          /*\warning we assume that W is symmetric in the Lagrangian case*/
          double value = -1. * columntmp->dot(d.twist_read());
          if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
            value += (d.p_read(_levelMaxForInput))(itindex);
          }
          /* \warning the computation of reactionToBoundaryConditions take into
             account the contact impulse but not the external and internal
             forces. A complete computation of the residu should be better */
          (*d.reactionToBoundaryConditions())(bc) = value;
          bc++;
        }
      }
      moreau_jean::updatePosition(_simulation->timeStep(), _theta, ds);
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanOSI::updateState - not yet "
          "implemented for "
          "Dynamical system of type: " +
          siconos::types::type_name(ds));
  }
  DEBUG_END("siconos::integrators::MoreauJeanOSI::updateState(const unsigned int)\n");
}

bool siconos::integrators::MoreauJeanOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction>"
      " inter, "
      "unsigned int i)\n");

  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (*inter->y(i - 1))(0);   // for i=1 y(i-1) is the position
  double yDot = (*(inter->y(i)))(0);  // for i=1 y(i) is the velocity
#ifdef DEBUG_ACTIVATION
  double yDot_k = (inter->y_k(i))(0);
  // for i=1 y(i) is the velocity et the beginning of the time step
#endif
  double gamma = 1.0 / 2.0;
  if (_useGamma) {
    gamma = _gamma;
  }
  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanOSI::addInteractionInIndexSet of level "
      "= %i yref=%e, "
      "yDot=%e, y_estimated=%e.,  _constraintActivationThreshold=%e\n",
      i, y, yDot, y + gamma * h * yDot, _constraintActivationThreshold);
  y += gamma * h * yDot;
  assert(!std::isnan(y));

  if (_activateWithNegativeRelativeVelocity) {
#ifdef DEBUG_ACTIVATION
    if (fabs(yDot - yDot_k) > 1e-10) {
      std::cout << "MoreauJeanOSI::addInteractionInIndexSet ACTIVATE." << y
                << "<= " << _constraintActivationThreshold << std::endl;

      getchar();
    }

    if (y <= _constraintActivationThreshold) {
      if (not(yDot <= _constraintActivationThresholdVelocity)) {
        std::cout << "\n MoreauJeanOSI::addInteractionInIndexSet activation at the position "
                     "level but not at the velocity level."
                  << "\n number :" << inter->number() << " y=" << y
                  << "<= " << _constraintActivationThreshold << " yDot_k =" << yDot_k << " > "
                  << _constraintActivationThresholdVelocity << " yDot =" << yDot << " > "
                  << _constraintActivationThresholdVelocity << std::endl;
        // getchar();
      } else {
        std::cout << "\n MoreauJeanOSI::addInteractionInIndexSet activation at the position "
                     "level but not at the velocity level."
                  << "\n number :" << inter->number() << " y=" << y
                  << "<= " << _constraintActivationThreshold << " yDot_k =" << yDot_k
                  << "<= " << _constraintActivationThresholdVelocity << " yDot =" << yDot
                  << " <= " << _constraintActivationThresholdVelocity << std::endl;
      }
    }
    DEBUG_EXPR_WE(
        if ((y <= _constraintActivationThreshold) and
            (yDot_k <= _constraintActivationThreshold)) std::cout
            << "MoreauJeanOSI::addInteractionInIndexSet ACTIVATE with velocity threshold."
            << " gamma " << gamma << " y=" << y << "<= " << _constraintActivationThreshold
            << " yDot_k =" << yDot_k << "<= " << _constraintActivationThresholdVelocity
            << std::endl;);
#endif
    return ((y <= _constraintActivationThreshold) and
            (yDot <= _constraintActivationThresholdVelocity));
  } else {
#ifdef DEBUG_ACTIVATION
    if (fabs(yDot - yDot_k) > 1e-10) {
      if (y <= _constraintActivationThreshold)
        std::cout << "MoreauJeanOSI::addInteractionInIndexSet ACTIVATE." << y
                  << "<= " << _constraintActivationThreshold << std::endl;

      getchar();
    }
    if (y <= _constraintActivationThreshold) {
      std::cout << "ACTIVATED "
                << "number :" << inter->number() << " y=" << y
                << "<= " << _constraintActivationThreshold << " yDot_k =" << yDot_k
                << " yDot =" << yDot << std::endl;
    } else {
      std::cout << "NOT ACTIVATED "
                << " number :" << inter->number() << " y=" << y
                << "<= " << _constraintActivationThreshold << " yDot_k =" << yDot_k
                << " yDot =" << yDot << std::endl;
      // getchar();
    }

#endif
    DEBUG_EXPR_WE(if (y <= _constraintActivationThreshold) std::cout
                      << "MoreauJeanOSI::addInteractionInIndexSet ACTIVATE." << y
                      << "<= " << _constraintActivationThreshold << std::endl;);
    return (y <= _constraintActivationThreshold);
  }
}

bool siconos::integrators::MoreauJeanOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  return !(addInteractionInIndexSet(inter, i));
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::MoreauJeanOSI::computeWorkForces() {
  DEBUG_BEGIN("MoreauJeanOSI::computeWorkForces()\n");

  double t = _simulation->nextTime();         // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double h = t - told;                        // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //

  auto number_of_ds = _simulation->nonSmoothDynamicalSystem()->getNumberOfDS();
  auto workForces = std::make_shared<siconos::algebra::SiconosMatrix>(number_of_ds, 2);

  size_t cnt_ds = 0;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    auto dsType = siconos::types::type_value(*ds);
    // 3 - Lagrangian Non Linear Systems
    if (dsType == siconos::modeling::Type::LagrangianDS ||
        dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
      // -- Convert the DS into a Lagrangian one.
      auto lds = static_pointer_cast<siconos::modeling::LagrangianDS>(ds);
      siconos::algebra::SiconosVector f_k_theta{lds->dimension()};
      siconos::algebra::SiconosVector v_k_theta{lds->dimension()};
      f_k_theta.setZero();
      v_k_theta.setZero();

      if (lds->hasTotalForces()) {
        // Cheaper version: get forces(ti,vi,qi) from memory
        f_k_theta += (1. - _theta) * lds->forcesMemory().getSiconosVector(0);
        v_k_theta += (1. - _theta) * lds->velocityMemory().getSiconosVector(0);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        lds->computeTotalForces(lds->velocity_read(), lds->q_read(), t);
        f_k_theta += _theta * lds->totalForces();
        v_k_theta += _theta * lds->velocity_read();

        DEBUG_PRINT("MoreauJeanOSI:: new forces :\n");
        DEBUG_EXPR(siconos::algebra::print(*lds->totalForces()););
        DEBUG_EXPR(siconos::algebra::print(*f_k_theta););

        // scalar product
        (*workForces)(ds->number(), 0) = ds->number();
        (*workForces)(ds->number(), 1) = h * f_k_theta.dot(v_k_theta);

        DEBUG_EXPR(siconos::algebra::print(*workForces););
      }
    } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
      DEBUG_PRINT("MoreauJeanOSI::computeWorkForces(), dsType == Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) -
      // h*(1-_theta)*forces(ti,vi,qi) - pi+1

      // -- Convert the DS into a NewtonEulerDS one.
      auto neds = static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds);

      siconos::algebra::SiconosVector f_k_theta{neds->dimension()};
      siconos::algebra::SiconosVector v_k_theta{neds->dimension()};

      DEBUG_PRINTF("MoreauJeanOSI:: _theta = %e\n", _theta);
      DEBUG_PRINTF("MoreauJeanOSI:: h = %e\n", h);

      // Cheaper version: get forces(ti,vi,qi) from memory
      f_k_theta = (1. - _theta) * neds->forcesMemory().getSiconosVector(0);
      v_k_theta = (1. - _theta) * neds->twistMemory().getSiconosVector(0);

      // computes forces(ti,v,q)
      neds->computeWrench(neds->twist_read(), neds->q_read(), t);
      f_k_theta += _theta * neds->wrench();
      v_k_theta = _theta * neds->twist_read();

      DEBUG_PRINT("MoreauJeanOSI:: new forces :\n");
      DEBUG_EXPR(siconos::algebra::print(*neds->totalForces()););
      DEBUG_EXPR(siconos::algebra::print(*f_k_theta););

      // scalar product
      (*workForces)(cnt_ds, 0) = ds->number();
      (*workForces)(cnt_ds, 1) = h * f_k_theta.dot(v_k_theta);

      cnt_ds++;

      DEBUG_PRINT("MoreauJeanOSI::computeWorkForces :\n");
    } else
      THROW_EXCEPTION(
          "MoreauJeanOSI::computeWorkForces - not yet implemented this kind of Dynamical "
          "system.");
  }
  DEBUG_END("MoreauJeanOSI::computeWorkForces()\n");
  return workForces;
}

void siconos::integrators::MoreauJeanOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanOSI OSI display ======" << std::endl;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  if (_dynamicalSystemsGraph) {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds =
          _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------" << std::endl;
      std::cout << "--> W of dynamical system number " << ds->number() << ": " << std::endl;
      if (_dynamicalSystemsGraph->properties(*dsi).iterationMatrix)
        siconos::algebra::print(*_dynamicalSystemsGraph->properties(*dsi).iterationMatrix);
      else
        std::cout << "-> nullptr" << std::endl;
      std::cout << "--> and corresponding theta is: " << _theta << std::endl;
    }
  }
  std::cout << "================================" << std::endl;
}

// --------------- Free functions  ---------------
// void siconos::integrators::moreau_jean::computeIterationMatrix_Lagrangian(
//     double time, double time_step, double theta, siconos::modeling::LagrangianDS &lds,
//     Eigen::Ref<siconos::algebra::SiconosMatrix> iterationMatrix,
//     std::shared_ptr<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>> &luw) {
//   // Compute W matrix of the Dynamical System ds, at time t and for the current
//   // ds state.
//   DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeIterationMatrix\n");

//   if (lds.isLinear()) return;  // Nothing: W does not depend on time.

//   if (lds.hasMass()) {
//     lds.computeMass(lds.q_read());
//     iterationMatrix = lds.mass();  // copy
//   } else
//     iterationMatrix.setIdentity();

//   lds.computeJacobianTotalForcesOver_velocity(lds.velocity_read(), lds.q_read(), time);
//   if (lds.hasJacobianTotalForcesOver_velocity()) {
//     iterationMatrix -= time_step * theta * lds.jacobianTotalForcesOver_velocity();
//   }

//   lds.computeJacobianTotalForcesOver_q(lds.velocity_read(), lds.q_read(), time);
//   if (lds.hasJacobianTotalForcesOver_q()) {
//     iterationMatrix -= time_step * time_step * theta * theta *
//     lds.jacobianTotalForcesOver_q();
//   }

//   // LU Factorize the iteration matrix
//   luw =
//   std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(iterationMatrix);
// }

void siconos::integrators::moreau_jean::computeIterationMatrix_NewtonEuler(
    double time, double time_step, double theta, siconos::modeling::NewtonEulerDS &neds,
    Eigen::Ref<siconos::algebra::SiconosMatrix> iterationMatrix,
    std::shared_ptr<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>> &luw) {
  // Compute W matrix of the Dynamical System ds, at time t and for the current
  // ds state.
  iterationMatrix = neds.totalInertiaMatrix();
  auto coeff = time_step * theta;
  if (neds.hasJacobianWrenchOver_twist()) {
    neds.computeJacobianWrenchOver_twist(neds.twist_read(), neds.q_read(), time);
    iterationMatrix -= coeff * neds.jacobianWrenchOver_twist();
  }
  if (neds.hasJacobianWrenchOver_q()) {
    neds.computeJacobianWrenchOver_q(neds.twist_read(), neds.q_read(), time);

    iterationMatrix -= coeff * coeff * (neds.jacobianWrenchOver_q() * neds.T());
    //*W -= h*h*_theta*_theta**K;
  }

  // LU Factorize the iteration matrix
  luw = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(iterationMatrix);
}

void siconos::integrators::moreau_jean::updatePosition(
    double time_step, double theta, siconos::modeling::DynamicalSystem &ds) {
  DEBUG_BEGIN("siconos::integrators::moreau_jean::updatePosition(...)\n");

  auto dsType = siconos::types::type_value(ds);
  // 1 - Lagrangian Systems
  if (dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    siconos::modeling::LagrangianDS &lds = static_cast<siconos::modeling::LagrangianDS &>(ds);
    // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    *(lds.q()) = time_step * theta * lds.velocity_read() +
                 time_step * (1. - theta) * lds.velocityMemory().getSiconosVector(0) +
                 lds.qMemory().getSiconosVector(0);

  } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
    siconos::modeling::NewtonEulerDS &neds =
        static_cast<siconos::modeling::NewtonEulerDS &>(ds);

    siconos::algebra::SiconosVector6 velocityIncrement{neds.dimension()};
    velocityIncrement = time_step * theta * neds.twist_read() +
                        time_step * (1. - theta) * neds.twistMemory().getSiconosVector(0);
    DEBUG_EXPR(siconos::algebra::print(*velocityIncrement));

    neds.q()->head(3) = velocityIncrement.head(3);
    siconos::geometry::quaternionFromTwistVector(velocityIncrement, *neds.q());

    DEBUG_EXPR(siconos::algebra::print(q));
    siconos::geometry::compositionLawLieGroup(neds.qMemory().getSiconosVector(0), *neds.q());
    DEBUG_EXPR(siconos::algebra::print(q));
  }
  DEBUG_END("siconos::integrators::moreau_jean::updatePosition(...)");
}
