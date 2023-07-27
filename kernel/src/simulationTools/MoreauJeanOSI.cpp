/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "Interaction.hpp"
#include "LagrangianCompliantLinearTIR.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianRheonomousR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonImpactFrictionNSL.hpp"         // for nslaw visitor
#include "NewtonImpactNSL.hpp"                 // for nslaw visitor
#include "NewtonImpactRollingFrictionNSL.hpp"  // for nslaw visitor
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "RotationQuaternion.hpp"  // for quaternionFromTwistVector and compositionLawLieGroup
#include "SiconosException.hpp"
#include "SiconosMatrixOp.hpp"        // for prod, scal, ...
#include "SiconosMatrixVectorOp.hpp"  // for prod, subprod ...
#include "SiconosPointers.hpp"        // For createSPtr
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for prod, subprod ...
#include "SiconosVisitor.hpp"
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "../../mechanics/src/fem/native/SolidLinearTIDS.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::
    _NSLEffectOnFreeOutput(siconos::nonsmooth_formulations::OneStepNSProblem &p,
                           siconos::modeling::Interaction &inter,
                           siconos::graphs::InteractionProperties &interProp)
    : _osnsp(p), _inter(inter), _interProp(interProp){};

void siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::NewtonImpactNSL &nslaw) const {
  double e;
  e = nslaw.e();
  std::vector<std::size_t> subCoord(4);
  subCoord[0] = 0;
  subCoord[1] = _inter.nonSmoothLaw()->size();
  subCoord[2] = 0;
  subCoord[3] = subCoord[1];
  auto &osnsp_rhs = *(
      *_interProp.workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];
  siconos::algebra::subscal(e, _inter.y_k(_osnsp.inputOutputLevel()), osnsp_rhs,
                            subCoord, false);
}

void siconos::integrators::MoreauJeanOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::NewtonImpactFrictionNSL &nslaw) const {
  auto &osnsp_rhs = *(
      *_interProp.workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];

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
    const siconos::modeling::NewtonImpactRollingFrictionNSL &nslaw) const {
  auto &osnsp_rhs = *(
      *_interProp.workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];

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

const siconos::algebra::SimpleMatrix siconos::integrators::MoreauJeanOSI::getW(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds && "siconos::integrators::MoreauJeanOSI::getW(ds): ds == nullptr.");
  assert(
      _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
          .W &&
      "siconos::integrators::MoreauJeanOSI::getW(ds): W[ds] == nullptr.");
  // Copy !!
  return *_dynamicalSystemsGraph
              ->properties(_dynamicalSystemsGraph->descriptor(ds))
              .W;
}

std::shared_ptr<siconos::algebra::SimpleMatrix>
siconos::integrators::MoreauJeanOSI::W(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds && "siconos::integrators::MoreauJeanOSI::W(ds): ds == nullptr.");
  return _dynamicalSystemsGraph
      ->properties(_dynamicalSystemsGraph->descriptor(ds))
      .W;
}

const siconos::algebra::SimpleMatrix
siconos::integrators::MoreauJeanOSI::getWBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds &&
         "siconos::integrators::MoreauJeanOSI::getWBoundaryConditions(ds): ds "
         "== nullptr.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(
      _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
          .WBoundaryConditions &&
      "siconos::integrators::MoreauJeanOSI::getWBoundaryConditions(ds): "
      "WBoundaryConditions[ds] == nullptr.");
  return *(
      _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
          .WBoundaryConditions);  // Copy !!
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::MoreauJeanOSI::WBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds &&
         "siconos::integrators::MoreauJeanOSI::WBoundaryConditions(ds): ds == "
         "nullptr.");
  return _dynamicalSystemsGraph
      ->properties(_dynamicalSystemsGraph->descriptor(ds))
      .W;
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
         dsType == siconos::modeling::Type::NewtonEulerDS ||
         dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS);

  // Compute W (iteration matrix)
  auto sods = std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds);
  initializeIterationMatrixW(t, sods);

  // Initialize work vectors
  auto &ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(siconos::integrators::MoreauJeanOSI::WORK_LENGTH);

  ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE] =
      std::make_shared<siconos::algebra::SiconosVector>(sods->dimension());
  ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE] =
      std::make_shared<siconos::algebra::SiconosVector>(sods->dimension());

  if (dsType == siconos::modeling::Type::SolidLinearTIDS){
      auto solidds = std::static_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
      ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_SIGMAFREE] =
          std::make_shared<siconos::algebra::SiconosVector>(solidds->stressDimension());
      ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_SIGMAFREE] =
          std::make_shared<siconos::algebra::SiconosVector>(solidds->stressDimension());
  }

  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    // buffers allocation (inside the graph)
    auto lds = std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds);
    ds_work_vectors[siconos::integrators::MoreauJeanOSI::BUFFER] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
  } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
    auto neds = std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds);
    DEBUG_PRINTF("neds->number() %i \n", neds->number());
    // Compute a first value of the dotq  to store it in  _dotqMemory
    std::shared_ptr<siconos::algebra::SiconosMatrix> T = neds->T();
    std::shared_ptr<siconos::algebra::SiconosVector> dotq = neds->dotq();
    std::shared_ptr<siconos::algebra::SiconosVector> v = neds->twist();
    siconos::algebra::prod(*T, *v, *dotq, true);
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::initializeWorkVectorsForDS(Model&, "
      "double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
  // Update dynamical system components (for memory swap).
  sods->computeForces(t, sods->q(), sods->velocity());
  sods->swapInMemory();
}
void siconos::integrators::MoreauJeanOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction &inter,
    siconos::graphs::InteractionProperties &interProp,
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
    interProp.workVectors = std::make_shared<
        std::vector<std::shared_ptr<siconos::algebra::SiconosVector> > >(
        siconos::integrators::MoreauJeanOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors = std::make_shared<
        std::vector<std::shared_ptr<siconos::algebra::BlockVector> > >(
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
    if ((!inter_work_block[xfree]) ||
        (inter_work_block[xfree]->numberOfBlocks() != 2))
      inter_work_block[xfree] =
          std::make_shared<siconos::algebra::BlockVector>(2);
  } else {
    if ((!inter_work_block[xfree]) ||
        (inter_work_block[xfree]->numberOfBlocks() != 1))
      inter_work_block[xfree] =
          std::make_shared<siconos::algebra::BlockVector>(1);
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

  DEBUG_EXPR(inter_work_block[xfree]->display(););
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::"
      "initializeWorkVectorsForInteraction(Interaction "
      "&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");
}

void siconos::integrators::MoreauJeanOSI::initialize_nonsmooth_problems() {
  auto allOSNS = _simulation->oneStepNSProblems();
  ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY])
      ->setIndexSetLevel(1);
  ((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY])
      ->setInputOutputLevel(1);
  //  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->initialize(_simulation);
}

void siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW(
    double time, std::shared_ptr<siconos::modeling::SecondOrderDS> ds) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW\n");
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical
  // system, depending on its type.
  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW(t,ds) "
        "- ds == "
        "nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW(t,ds) "
        "- ds does not "
        "belong to the OSI.");

  const siconos::graphs::DynamicalSystemsGraph::VDescriptor &dsv =
      _dynamicalSystemsGraph->descriptor(ds);

  if (_dynamicalSystemsGraph->properties(dsv).W)
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW(t,ds) "
        "- W(ds) is "
        "already in the map and has been initialized.");

  double h = _simulation->timeStep();
  auto dsType = siconos::types::type_value(*ds);
  auto sizeW = ds->dimension();
  if (dsType == siconos::modeling::Type::LagrangianDS) {
    auto &d = static_cast<siconos::modeling::LagrangianDS &>(*ds);
    // Memory allocation for W property of the grap
    if (d.mass()) {
      d.computeMass(d.q());
      _dynamicalSystemsGraph->properties(dsv).W =
          std::make_shared<siconos::algebra::SimpleMatrix>(
              *d.mass());  //*W = *d->mass();
    } else {
      _dynamicalSystemsGraph->properties(dsv).W =
          std::make_shared<siconos::algebra::SimpleMatrix>(sizeW, sizeW);
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }
    // Compute the W matrix
    computeW(time, d, *_dynamicalSystemsGraph->properties(dsv).W);
    // WBoundaryConditions initialization
    if (d.boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(*ds, dsv);
  }
  // 2 - Lagrangian linear systems
  else if (dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
    auto d =
        std::static_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds);
    if (d->mass()) {
      _dynamicalSystemsGraph->properties(dsv).W =
          std::make_shared<siconos::algebra::SimpleMatrix>(
              *d->mass());  //*W = *d->mass();
    } else {
      _dynamicalSystemsGraph->properties(dsv).W =
          std::make_shared<siconos::algebra::SimpleMatrix>(sizeW, sizeW);
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }

    auto K = d->K();
    auto C = d->C();
    auto W = _dynamicalSystemsGraph->properties(dsv).W;
    if (C)
      siconos::algebra::scal(h * _theta, *C, *W, false);  // W += h*_theta *C
    if (K)
      siconos::algebra::scal(h * h * _theta * _theta, *K, *W,
                             false);  // W = h*h*_theta*_theta*K

    // WBoundaryConditions initialization
    if (d->boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(*d, dsv);
  }   else if (dsType == siconos::modeling::Type::SolidLinearTIDS) {
      auto d =
          std::static_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
      if (d->mass()) {
        _dynamicalSystemsGraph->properties(dsv).W =
            std::make_shared<siconos::algebra::SimpleMatrix>(
                *d->mass());  //*W = *d->mass();
      } else {
        _dynamicalSystemsGraph->properties(dsv).W =
            std::make_shared<siconos::algebra::SimpleMatrix>(sizeW, sizeW);
        _dynamicalSystemsGraph->properties(dsv).W->eye();
      }

      auto K = d->K();
      auto C = d->C();
      auto S = d->S();
      auto B = d->B();
      auto W = _dynamicalSystemsGraph->properties(dsv).W;
      if (d->B()) {
          auto Btrans = std::make_shared<siconos::algebra::SimpleMatrix>(B->size(1),B->size(0));
          double coeff = h * h * _theta;
          d->mass()->Solve(*Btrans); // Btrans = M^-1 B^T
          siconos::algebra::prod(*B, *Btrans, *W,
                                 false);  // W = B M^-1 B^T
          if (S)
              siconos::algebra::scal(1.0, *S, *W,
                                     false);  // W += S

      // WBoundaryConditions initialization
      if (d->boundaryConditions())
        _initializeIterationMatrixWBoundaryConditions(*d, dsv);
    }

  else if (dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    auto &lldds =
        static_cast<siconos::modeling::LagrangianLinearDiagonalDS &>(*ds);
    auto ndof = lldds.dimension();
    _dynamicalSystemsGraph->properties(dsv).W =
        std::make_shared<siconos::algebra::SimpleMatrix>(
            ndof, ndof, siconos::algebra::UblasType::BANDED, 0, 0);
    auto &W = *_dynamicalSystemsGraph->properties(dsv).W;

    if (lldds.mass())
      W = *lldds.mass();
    else
      W.eye();

    double htheta = h * _theta;
    double h2theta2 = h * h * _theta * _theta;
    if (lldds.damping()) {
      auto &C = *lldds.damping();
      for (unsigned int i = 0; i < ndof; ++i) {
        W(i, i) += htheta * C(i);
      }
    }

    if (lldds.stiffness()) {
      auto &K = *lldds.stiffness();
      for (unsigned int i = 0; i < ndof; ++i) {
        W(i, i) += h2theta2 * K(i);
      }
    }
    // WBoundaryConditions initialization
    if (lldds.boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(lldds, dsv);

    for (unsigned int i = 0; i < ndof; ++i) {
      W(i, i) = 1. / W(i, i);
    }
  }

  // === ===

  else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
    auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(*ds);
    _dynamicalSystemsGraph->properties(dsv).W =
        std::make_shared<siconos::algebra::SimpleMatrix>(*d.mass());

    computeW(time, d, *_dynamicalSystemsGraph->properties(dsv).W);

    // WBoundaryConditions initialization
    if (d.boundaryConditions())
      _initializeIterationMatrixWBoundaryConditions(*ds, dsv);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW - not "
        "yet "
        "implemented for Dynamical system of type : " +
        siconos::types::type_name(*ds));

  if (isWSymmetricDefinitePositive()) {
    auto &W = *_dynamicalSystemsGraph->properties(dsv).W;
    W.setIsSymmetric(true);
    W.setIsPositiveDefinite(true);
  }

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::initializeIterationMatrixW\n");
}

void siconos::integrators::MoreauJeanOSI::
    _initializeIterationMatrixWBoundaryConditions(
        siconos::modeling::SecondOrderDS &ds,
        const siconos::graphs::DynamicalSystemsGraph::VDescriptor &dsv) {
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::"
      "initializeIterationMatrixWBoundaryConditions(std::"
      "shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");

  if (!(checkOSI(dsv)))
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::"
        "initializeIterationMatrixWBoundaryConditions(t,"
        "ds) - ds does not belong to the OSI.");

  if (_dynamicalSystemsGraph->properties(dsv).WBoundaryConditions)
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::"
        "initializeIterationMatrixWBoundaryConditions(t,"
        "ds) - WBoundaryConditions(ds) is already in the map and has been "
        "initialized.");

  auto dsType = siconos::types::type_value(ds);
  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::NewtonEulerDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    // Memory allocation for WBoundaryConditions
    auto sizeWBoundaryConditions =
        ds.dimension();  // n for first order systems, ndof for lagrangian.

    auto &d = static_cast<siconos::modeling::SecondOrderDS &>(ds);
    _dynamicalSystemsGraph->properties(dsv).WBoundaryConditions =
        std::make_shared<siconos::algebra::SimpleMatrix>(
            sizeWBoundaryConditions, d.boundaryConditions()->size(),
            d.mass()->num());
    _computeWBoundaryConditions(
        ds, *_dynamicalSystemsGraph->properties(dsv).WBoundaryConditions,
        *_dynamicalSystemsGraph->properties(dsv).W);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::"
        "initializeIterationMatrixWBoundaryConditions - "
        "not yet implemented for Dynamical system of type :" +
        siconos::types::type_name(ds));
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::"
      "initializeIterationMatrixWBoundaryConditions(std::"
      "shared_ptr<siconos::modeling::DynamicalSystem> ds) \n");
}

void siconos::integrators::MoreauJeanOSI::_computeWBoundaryConditions(
    siconos::modeling::SecondOrderDS &ds,
    siconos::algebra::SiconosMatrix &WBoundaryConditions,
    siconos::algebra::SiconosMatrix &iteration_matrix) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::_computeWBoundaryConditions\n");
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  auto dsType = siconos::types::type_value(ds);
  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::NewtonEulerDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    auto columntmp = std::make_shared<siconos::algebra::SiconosVector>(
        ds.dimension(), WBoundaryConditions.num());

    int columnindex = 0;

    auto &d = static_cast<siconos::modeling::SecondOrderDS &>(ds);

    if (!iteration_matrix.checkSymmetry(
            1e-10))  // Warning this operation could be quite expensive
    {
      // iteration_matrix.display();
      std::cout << "Warning, we apply boundary conditions assuming W symmetric"
                << std::endl;
    }

    for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
      iteration_matrix.getCol(itindex, *columntmp);
      /*\warning we assume that W is symmetric
        we store only the column and not the row */
      WBoundaryConditions.setCol(columnindex, *columntmp);
      columnindex++;
    }

    columnindex = 0;
    for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
      double diag = iteration_matrix.getValue(itindex, itindex);
      columntmp->zero();
      (*columntmp)(itindex) = diag;
      iteration_matrix.setCol(itindex, *columntmp);
      iteration_matrix.setRow(itindex, *columntmp);
      columnindex++;
    }
    DEBUG_EXPR(iteration_matrix.display(););
    DEBUG_EXPR(WBoundaryConditions.display(););
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::computeWBoundaryConditions - not "
        "yet "
        "implemented for Dynamical system type : " +
        siconos::types::type_name(ds));
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::_computeWBoundaryConditions\n");
}

void siconos::integrators::MoreauJeanOSI::computeW(
    double t, siconos::modeling::SecondOrderDS &ds,
    siconos::algebra::SiconosMatrix &W) {
  // Compute W matrix of the Dynamical System ds, at time t and for the current
  // ds state.
  DEBUG_BEGIN("siconos::integrators::MoreauJeanOSI::computeW\n");

  double h = _simulation->timeStep();
  auto dsType = siconos::types::type_value(ds);

  if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    // Nothing: W does not depend on time.
  } else if (dsType == siconos::modeling::Type::LagrangianDS) {
    auto &d = static_cast<siconos::modeling::LagrangianDS &>(ds);
    if (d.mass()) {
      d.computeMass();
      W = *d.mass();
    } else
      W.eye();

    if (d.jacobianvForces()) {
      auto &C = *d.jacobianvForces();  // jacobian according to velocity
      d.computeJacobianqDotForces(t);
      siconos::algebra::scal(-h * _theta, C, W, false);  // W -= h*_theta*C
    }

    if (d.jacobianqForces()) {
      auto &K = *d.jacobianqForces();  // jacobian according to q
      d.computeJacobianqForces(t);
      siconos::algebra::scal(-h * h * _theta * _theta, K, W,
                             false);  //*W -= h*h*_theta*_theta**K;
    }
  }
  // === ===
  else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
    auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);
    W = *(d.mass());

    if (d.jacobianvForces()) {
      auto &C = *d.jacobianvForces();  // jacobian according to velocity

      d.computeJacobianvForces(t);
      siconos::algebra::scal(-h * _theta, C, W, false);  // W -= h*_theta*C
    }
    if (d.jacobianqForces()) {
      auto &K = *d.jacobianqForces();  // jacobian according to q
      d.computeJacobianqForces(t);
      auto &T = *d.T();
      DEBUG_EXPR(T.display(););
      DEBUG_EXPR(K.display(););
      auto buffer = std::make_shared<siconos::algebra::SimpleMatrix>(*d.mass());
      siconos::algebra::prod(K, T, *buffer, true);
      siconos::algebra::scal(-h * h * _theta * _theta, *buffer, W, false);
      //*W -= h*h*_theta*_theta**K;
    }
    DEBUG_EXPR(W.display(););
    DEBUG_EXPR_WE(std::cout << std::boolalpha << "W.isFactorized() = "
                            << W.isFactorized() << std::endl;);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanOSI::computeW - not yet implemented "
        "for Dynamical "
        "system of type : " +
        siconos::types::type_name(ds));

  DEBUG_END("siconos::integrators::MoreauJeanOSI::computeW\n");
  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}

std::shared_ptr<siconos::algebra::SimpleMatrix>
siconos::integrators::MoreauJeanOSI::Winverse(
    std::shared_ptr<siconos::modeling::SecondOrderDS> ds, bool keepW) {
  /* We compute and return the current inverse the W matrix */

  const auto &dsv = _dynamicalSystemsGraph->descriptor(ds);

  auto W = _dynamicalSystemsGraph->properties(dsv).W;
  auto Winverse = _dynamicalSystemsGraph->properties(dsv).Winverse;
  if (!Winverse) {
    auto sizeW = ds->dimension();
    _dynamicalSystemsGraph->properties(dsv).Winverse =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeW, sizeW);
    Winverse = _dynamicalSystemsGraph->properties(dsv).Winverse;
    Winverse->eye();
    if (keepW) {
      // std::cout << "MoreauJeanOSI keepW" << std::endl;
      auto Wtmp = std::make_shared<siconos::algebra::SimpleMatrix>(*W);
      Wtmp->Solve(*Winverse);
    } else {
      W->Solve(*Winverse);
    }
  } else {
    auto dsType = siconos::types::type_value(*ds);
    if (dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
        dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      // Nothing: W does not depend on time.
    } else {
      Winverse->eye();
      if (keepW) {
        auto Wtmp = std::make_shared<siconos::algebra::SimpleMatrix>(*W);
        Wtmp->Solve(*Winverse);
      } else {
        W->Solve(*Winverse);
      }
    }
  }

  return Winverse;
}

void siconos::integrators::MoreauJeanOSI::computeInitialNewtonState() {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::computeInitialNewtonState()\n");
  // Compute the position value giving the initial velocity.
  // The goal is to save one newton iteration for nearly linear system
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend;
       ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);

    if (_explicitNewtonEulerDSOperators) {
      if (siconos::types::type_value(ds) ==
          siconos::modeling::Type::NewtonEulerDS) {
        // The goal is to update T() one time at the beginning of the Newton
        // Loop We want to be explicit on this function since we do not compute
        // their Jacobians.
        auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);
        auto qold = d.qMemory().getSiconosVector(0);
        auto pqold =
            siconos::pointers::createSPtr<siconos::algebra::SiconosVector>(
                qold);

        // std::shared_ptr<siconos::algebra::SiconosVector> q = d->q();
        siconos::modeling::computeT(pqold, d.T());
      }
    }
    // The goal is to converge in one iteration if the system is almost linear
    // we start the Newton loop q = q0+hv0
    updatePosition(ds);
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::computeInitialNewtonState()\n");
}

void siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(
    siconos::modeling::SecondOrderDS &d,
    siconos::algebra::SiconosVector &residu,
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, double t,
    const siconos::algebra::SiconosVector &v) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(...)\n");
  if (d.boundaryConditions()) {
    d.boundaryConditions()->computePrescribedVelocity(t);

    unsigned int columnindex = 0;
    auto &WBoundaryConditions =
        *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions;
    auto columntmp =
        std::make_shared<siconos::algebra::SiconosVector>(d.dimension());

    for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
      double DeltaPrescribedVelocity =
          d.boundaryConditions()->prescribedVelocity()->getValue(columnindex) -
          v.getValue(itindex);
      DEBUG_PRINTF(
          "index  = %i, value = %e\n", *itindex,
          d.boundaryConditions()->prescribedVelocity()->getValue(columnindex));
      DEBUG_PRINTF("DeltaPrescribedVelocity = %e\n", DeltaPrescribedVelocity);
      WBoundaryConditions.getCol(columnindex, *columntmp);
      residu -= *columntmp * (DeltaPrescribedVelocity);

      residu.setValue(
          itindex, -columntmp->getValue(itindex) * (DeltaPrescribedVelocity));

      columnindex++;
    }
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::applyBoundaryConditions(...)\n");
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
  double h = t - told;                        // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  // std::shared_ptr<siconos::modeling::DynamicalSystem> ds; // Current
  // Dynamical System.
  siconos::modeling::Type dsType;  // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend;
       ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);
    auto &ds_work_vectors =
        *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = siconos::types::type_value(ds);  // Its type

    // 3 - Lagrangian Non Linear Systems
    if (dsType == siconos::modeling::Type::LagrangianDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::LagrangianDS\n");
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*forces(t_i+1,v_k,i+1, q_k,i+1)
      // - h*(1-theta)*forces(ti,vi,qi) - p_i+1
      auto &residuFree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::modeling::LagrangianDS &>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      const auto &vold = d.velocityMemory().getSiconosVector(0);

      const auto &v = *d.velocity();  // v = v_k,i+1
      // residuFree.zero();
      DEBUG_EXPR(residuFree.display());

      DEBUG_EXPR(vold.display());
      DEBUG_EXPR(v.display());

      residuFree = v;
      sub(residuFree, vold, residuFree);
      if (d.mass()) {
        d.computeMass(d.q());
        siconos::algebra::prod(*(d.mass()), residuFree,
                               residuFree);  // residuFree = M(v - vold)
      }

      if (d.forces()) {
        // Cheaper version: get forces(ti,vi,qi) from memory
        const auto &fold = d.forcesMemory().getSiconosVector(0);
        double coef = -h * (1 - _theta);
        siconos::algebra::scal(coef, fold, residuFree, false);

        // Expensive computes forces(ti,vi,qi)
        // auto &qold = *d.qMemory()->getSiconosVector(0);
        // auto &vold = *d.velocityMemory()->getSiconosVector(0);

        // d.computeForces(told, qold, vold);
        // double coef = -h * (1 - _theta);
        // // residuFree += coef * fL_i
        // siconos::algebra::scal(coef, *d.forces(), residuFree, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        d.computeForces(t, d.q(), d.velocity());
        coef = -h * _theta;
        siconos::algebra::scal(coef, *d.forces(), residuFree, false);

        // or  forces(ti+1, v_k,i+\theta, q(v_k,i+\theta))
        // std::shared_ptr<siconos::algebra::SiconosVector> qbasedonv(new
        // SiconosVector(*qold)); *qbasedonv +=  h * ((1 - _theta)* *vold +
        // _theta * *v); d.computeForces(t, qbasedonv, v); coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        // siconos::algebra::scal(coef, *d.forces(), *residuFree, false);
      }

      applyBoundaryConditions(d, residuFree, dsi, t, v);

      free = residuFree;  // copy residuFree into Workfree
      DEBUG_EXPR(residuFree.display());

      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!

      applyBoundaryConditions(d, free, dsi, t, v);

      DEBUG_EXPR(free.display());
      normResidu = free.norm2();
      DEBUG_PRINTF("normResidu= %e\n", normResidu);
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::LagrangianLinearTIDS\n");
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual
      // for v = v_i otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v
      // + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::modeling::LagrangianLinearTIDS &>(ds);

      auto &residuFree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      const auto &qold = d.qMemory().getSiconosVector(0);         // qi
      const auto &vold = d.velocityMemory().getSiconosVector(0);  // vi

      DEBUG_EXPR(qold.display(););
      DEBUG_EXPR(vold.display(););
      DEBUG_EXPR(d.q()->display(););
      DEBUG_EXPR(d.velocity()->display(););

      // --- ResiduFree computation Equation (1) ---
      residuFree.zero();
      double coeff;
      // -- No need to update W --

      if (d.C()) {
        siconos::algebra::prod(h, *d.C(), vold, residuFree,
                               false);  // vfree += h*C*vi
      }
      if (d.K()) {
        coeff = h * h * _theta;
        siconos::algebra::prod(coeff, *d.K(), vold, residuFree,
                               false);  // vfree += h^2*_theta*K*vi
        siconos::algebra::prod(h, *d.K(), qold, residuFree,
                               false);  // vfree += h*K*qi
      }

      if (d.fExt()) {
        // computes Fext(ti)
        d.computeFExt(told);
        coeff = -h * (1 - _theta);
        siconos::algebra::scal(coeff, *(d.fExt()), residuFree,
                               false);  // vfree -= h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d.computeFExt(t);
        coeff = -h * _theta;
        siconos::algebra::scal(coeff, *(d.fExt()), residuFree,
                               false);  // vfree -= h*_theta * fext(ti+1)
      }

      // Computation of the complete residual Equation (2)
      //   ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * (
      //   C*v + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      //       auto M = d.mass();
      //       std::shared_ptr<siconos::algebra::SiconosVector> realresiduFree
      //       (new SiconosVector(residuFree)); realresiduFree->zero();
      //       siconos::algebra::prod(*M,
      //       (*v-*vold), *realresiduFree); // residuFree = M(v - vold)
      //       std::shared_ptr<siconos::algebra::SiconosVector> qkplustheta (new
      //       SiconosVector(*qold)); qkplustheta->zero(); *qkplustheta = *qold
      //       + h
      //       *((1-_theta)* *vold + _theta* *v); if (C){
      //         double coef = h*(1-_theta);
      //         siconos::algebra::prod(coef, *C, *vold , *realresiduFree,
      //         false); coef = h*(_theta); siconos::algebra::prod(coef,*C, *v ,
      //         *realresiduFree, false);
      //       }
      //       if (K){
      //         double coef = h*(1-_theta);
      //         siconos::algebra::prod(coef,*K , *qold , *realresiduFree,
      //         false); coef = h*(_theta); siconos::algebra::prod(coef,*K ,
      //         *qkplustheta , *realresiduFree, false);
      //       }

      //       if (Fext)
      //       {
      //         // computes Fext(ti)
      //         d.computeFExt(told);
      //         coeff = -h*(1-_theta);
      //         siconos::algebra::scal(coeff, *Fext, *realresiduFree, false);
      //         // vfree -= h*(1-_theta) * fext(ti)
      //         // computes Fext(ti+1)
      //         d.computeFExt(t);
      //         coeff = -h*_theta;
      //         siconos::algebra::scal(coeff, *Fext, *realresiduFree, false);
      //         // vfree -= h*_theta * fext(ti+1)
      //       }

      applyBoundaryConditions(d, residuFree, dsi, t, vold);

      free = residuFree;            // copy residuFree into free
      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!
      // We use free as tmp buffer
      DEBUG_EXPR(free.display());
      DEBUG_EXPR(residuFree.display());

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();
    }
    else if (dsType == siconos::modeling::Type::SolidLinearTIDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::LagrangianLinearTIDS\n");
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i + h B sigma_{k+theta} + hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual
      // for v = v_i and sigma = sigma_i otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v
      // + K(q_i+h(1-theta)v_i+h theta v))) + h B sigma_{k+theta}
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::mechanics::fem::SolidLinearTIDS &>(ds);

      auto &residuFree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];
      auto &residuSigfreed = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_SIGMAFREE];
      auto &sigfreed = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::SIGMAFREE];

      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      const auto &qold = d.qMemory().getSiconosVector(0);         // qi
      const auto &vold = d.velocityMemory().getSiconosVector(0);  // vi
      const auto &sigmaold = d.stressMemory().getSiconosVector(0);  // vi

      DEBUG_EXPR(qold.display(););
      DEBUG_EXPR(vold.display(););
      DEBUG_EXPR(d.q()->display(););
      DEBUG_EXPR(d.velocity()->display(););

      // --- ResiduFree computation Equation (1) ---
      residuFree.zero();
      double coeff;
      // -- No need to update W --

      if (d.C()) {
        siconos::algebra::prod(h, *d.C(), vold, residuFree,
                               false);  // vfree += h*C*vi
      }
      if (d.K()) {
        coeff = h * h * _theta;
        siconos::algebra::prod(coeff, *d.K(), vold, residuFree,
                               false);  // vfree += h^2*_theta*K*vi
        siconos::algebra::prod(h, *d.K(), qold, residuFree,
                               false);  // vfree += h*K*qi
      }
      if (d.B()) {
          if (d.S()) {
              auto Btrans = std::make_shared<siconos::algebra::SimpleMatrix>(d.B()->size(1),d.B()->size(0));
              auto Bvold = std::make_shared<siconos::algebra::SiconosVector>(vold.size());
              coeff = h * h * _theta;
              siconos::algebra::prod(coeff, *d.B(), vold, *Bvold,
                                     true);  // Bvold = h*h*theta*B*v_{i}
              d.S()->Solve(*Bvold); // Bvold = S^-1 Bvold
              siconos::algebra::prod(h, *Btrans, *Bvold, residuFree,
                                     false);  // vfree += h*h*theta*B^T*S^-1*B*vold

              siconos::algebra::prod(h, *Btrans, sigmaold, residuFree,
                                     false);  // vfree += h*B*sigma_{i}
          }
      }

      if (d.fExt()) {
        // computes Fext(ti)
        d.computeFExt(told);
        coeff = -h * (1 - _theta);
        siconos::algebra::scal(coeff, *(d.fExt()), residuFree,
                               false);  // vfree -= h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d.computeFExt(t);
        coeff = -h * _theta;
        siconos::algebra::scal(coeff, *(d.fExt()), residuFree,
                               false);  // vfree -= h*_theta * fext(ti+1)
      }

      applyBoundaryConditions(d, residuFree, dsi, t, vold);

      free = residuFree;            // copy residuFree into free
      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!
      // We use free as tmp buffer
      DEBUG_EXPR(free.display());
      DEBUG_EXPR(residuFree.display());

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();
    }
    else if (dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      // ResiduFree = h*C*v_i + h*Kq_i +h*h*theta*Kv_i+hFext_theta     (1)
      // This formulae is only valid for the first computation of the residual
      // for v = v_i otherwise the complete formulae must be applied, that is
      // ResiduFree = M(v - vold) + h*((1-theta)*(C v_i + K q_i) +theta * ( C*v
      // + K(q_i+h(1-theta)v_i+h theta v)))
      //                     +hFext_theta     (2)
      // for v != vi, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      auto &d =
          static_cast<siconos::modeling::LagrangianLinearDiagonalDS &>(ds);

      auto &residuFree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      const auto &qold = d.qMemory().getSiconosVector(0);         // qi
      const auto &vold = d.velocityMemory().getSiconosVector(0);  // vi
      // --- ResiduFree computation Equation (1) ---
      residuFree.zero();
      double coeff;
      // -- No need to update W --
      if (d.damping()) {
        auto &sigma = *d.damping();
        for (unsigned int i = 0; i < d.dimension(); ++i)
          residuFree(i) += h * sigma(i) * vold(i);
      }
      if (d.stiffness()) {
        coeff = h * h * _theta;
        auto &omega = *d.stiffness();
        for (unsigned int i = 0; i < d.dimension(); ++i)
          residuFree(i) += coeff * omega(i) * vold(i) + h * omega(i) * qold(i);
      }

      if (d.fExt()) {
        // computes Fext(ti)
        d.computeFExt(told);
        coeff = -h * (1 - _theta);
        siconos::algebra::scal(coeff, *(d.fExt()), residuFree,
                               false);  // vfree -= h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d.computeFExt(t);
        coeff = -h * _theta;
        siconos::algebra::scal(coeff, *(d.fExt()), residuFree,
                               false);  // vfree -= h*_theta * fext(ti+1)
      }

      applyBoundaryConditions(d, residuFree, dsi, t, vold);

      free = residuFree;            // copy residuFree into free
      if (d.p(1)) free -= *d.p(1);  // Compute Residu in Workfree Notation !!

      normResidu = 0.0;  // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();
    }

    else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::computeResidu(), dsType == "
          "siconos::modeling::Type::NewtonEulerDS\n");
      // residu = M (v_k,i+1 - v_i) - h*_theta*forces(t,v_k,i+1, q_k,i+1) -
      // h*(1-_theta)*forces(ti,vi,qi) - pi+1

      auto &residuFree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
      auto &free = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      // -- Convert the DS into a Lagrangian one.
      auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);

      // Get the state  (previous time step) from memory vector
      // -> var. indexed with "Old"
      const auto &vold = d.twistMemory().getSiconosVector(0);

      // Get the current state vector
      // auto& q = *d.q();
      const auto &v = *d.twist();  // v = v_k,i+1

      // Get the (constant mass matrix)
      const auto &massMatrix = *d.mass();
      siconos::algebra::prod(massMatrix, (v - vold), residuFree,
                             true);  // residuFree = M(v - vold)
      DEBUG_EXPR(residuFree.display(););

      if (d.forces())  // if fL exists
      {
        DEBUG_PRINTF("siconos::integrators::MoreauJeanOSI:: _theta = %e\n",
                     _theta);
        DEBUG_PRINTF("siconos::integrators::MoreauJeanOSI:: h = %e\n", h);

        // Cheaper version: get forces(ti,vi,qi) from memory
        const auto &fold = d.forcesMemory().getSiconosVector(0);
        DEBUG_PRINT("siconos::integrators::MoreauJeanOSI:: old forces :\n");
        DEBUG_EXPR(fold.display(););

        double coef = -h * (1 - _theta);
        siconos::algebra::scal(coef, fold, residuFree, false);

        // Expensive version to check ...
        // std::shared_ptr<siconos::algebra::SiconosVector> qold =
        // d.qMemory()->getSiconosVector(0);
        // std::shared_ptr<siconos::algebra::SiconosVector> vold =
        // d.twistMemory()->getSiconosVector(0);
        //  d.computeForces(told,qold,vold);
        //  DEBUG_EXPR(d.forces()->display(););
        // double coef = -h * (1.0 - _theta);
        // siconos::algebra::scal(coef, *d.forces(), *residuFree, false);

        DEBUG_EXPR(residuFree.display(););

        // computes forces(ti,v,q)
        d.computeForces(t, d.q(), d.twist());
        coef = -h * _theta;
        siconos::algebra::scal(coef, *d.forces(), residuFree, false);
        DEBUG_PRINT("siconos::integrators::MoreauJeanOSI:: new forces :\n");
        DEBUG_EXPR(d.forces()->display(););
        DEBUG_EXPR(residuFree.display(););
      }

      applyBoundaryConditions(d, residuFree, dsi, t, v);

      free = residuFree;

      if (d.p(1)) free -= *d.p(1);

      applyBoundaryConditions(d, free, dsi, t, v);

      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::computeResidu :\n");
      DEBUG_EXPR(residuFree.display(););
      DEBUG_EXPR(if (d.p(1)) d.p(1)->display(););
      DEBUG_EXPR(free.display(););

      normResidu = free.norm2();
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

  double t = _simulation->nextTime();  // End of the time step

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  auto *rold =
  //  static_cast<SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  // std::shared_ptr<siconos::modeling::DynamicalSystem> ds; // Current
  // Dynamical System. std::shared_ptr<siconos::algebra::SiconosMatrix> W; // W
  // MoreauJeanOSI matrix of the current DS.
  siconos::modeling::Type dsType;  // Type of the current DS.

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend;
       ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);
    dsType = siconos::types::type_value(ds);  // Its type
    auto &W = *_dynamicalSystemsGraph->properties(*dsi)
                   .W;  // Its W MoreauJeanOSI matrix of iteration.
    auto &ds_work_vectors =
        *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    // // 3 - Lagrangian Non Linear Systems
    // if(dsType == siconos::modeling::Type::LagrangianDS ||
    //    dsType == siconos::modeling::Type::NewtonEulerDS)
    // {
    DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::computeFreeState()\n");
    // IN to be updated at current time: W, M, q, v, fL
    // IN at told: qi,vi, fLi

    // Note: indices i/i+1 corresponds to value at the beginning/end of the time
    // step. Index k stands for Newton iteration and thus corresponds to the
    // last computed value, ie the one saved in the DynamicalSystem. "i" values
    // are saved in memory vectors.

    // vFree = v_k,i+1 - W^{-1} ResiduFree
    // with
    // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*forces(t,v_k,i+1,
    // q_k,i+1) - h*(1-theta)*forces(ti,vi,qi)

    // -- Convert the DS into a Lagrangian one.
    auto &d = static_cast<siconos::modeling::SecondOrderDS &>(ds);
    const auto &vold = d.velocityMemory().getSiconosVector(0);  // vi (vold)
    const auto &v = *d.velocity();                              // v = v_k,i+1

    DEBUG_EXPR(v.display());
    DEBUG_EXPR(vold.display());

    // --- ResiduFree computation ---
    // ResFree = M(v-vold) - h*[theta*forces(t) + (1-theta)*forces(told)]
    //
    // vFree pointer is used to compute and save ResiduFree in this first step.
    auto &residuFree =
        *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
    auto &vfree = *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

    vfree = residuFree;
    DEBUG_EXPR(vfree.display());
    // -- Update W --
    // Note: during computeW, mass and jacobians of forces will be computed/
    if (dsType == siconos::modeling::Type::LagrangianDS ||
        dsType == siconos::modeling::Type::NewtonEulerDS) {
      computeW(t, d, W);
      if (d.boundaryConditions()) {
        _computeWBoundaryConditions(
            d, *_dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions,
            W);
      }
    }

    DEBUG_EXPR(W.display(););
    if (dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      // W is diagonal and contains the inverse of the iteration matrix!
      for (unsigned int i = 0; i < d.dimension(); ++i)
        vfree(i) = -W(i, i) * vfree(i) + vold(i);
    } else {
      // -- vfree =  v - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      W.Solve(vfree);
      // -> compute real vfree
      vfree *= -1.0;
      // Get state i (previous time step) from Memories -> var. indexed with
      // "Old"
      if (dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
        vfree += vold;
      } else {
        vfree += v;
      }
      DEBUG_EXPR(vfree.display());
    }
    // }
    // // 4 - Lagrangian Linear Systems
    // else if(dsType == siconos::modeling::Type::LagrangianLinearTIDS)
    // {
    //   DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::computeFreeState(),
    //   dsType == siconos::modeling::Type::LagrangianLinearTIDS\n");
    //   // IN to be updated at current time: Fext
    //   // IN at told: qi,vi, fext
    //   // IN constants: K,C

    //   // Note: indices i/i+1 corresponds to value at the beginning/end of the
    //   time step.
    //   // "i" values are saved in memory vectors.

    //   // vFree = v_i + W^{-1} ResiduFree    // with
    //   // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 +
    //   h*(1-theta)*Fext_i

    //   // -- Convert the DS into a Lagrangian one.
    //   LagrangianLinearTIDS& d = static_cast<LagrangianLinearTIDS&> (ds);

    //   // Get state i (previous time step) from Memories -> var. indexed with
    //   "Old" const auto& vold = d.velocityMemory().getSiconosVector(0); //vi

    //   // --- ResiduFree computation ---
    //   // vFree pointer is used to compute and save ResiduFree in this first
    //   step.

    //   // Velocity free and residu. vFree = RESfree (pointer equality !!).
    //   auto& residuFree =
    //   *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
    //   auto& vfree =
    //   *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

    //   vfree = residuFree;
    //   DEBUG_EXPR(vfree.display());
    //   W.Solve(vfree);
    //   vfree *= -1.0;
    //   vfree += vold;

    //   DEBUG_EXPR(vfree.display());

    // }
    // // 4 - Lagrangian Linear Diagonal Systems
    // else if(dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS)
    // {
    //   // IN to be updated at current time: Fext
    //   // IN at told: qi,vi, fext
    //   // IN constants: K,C

    //   // Note: indices i/i+1 corresponds to value at the beginning/end of the
    //   time step.
    //   // "i" values are saved in memory vectors.

    //   // vFree = v_i + W^{-1} ResiduFree    // with
    //   // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 +
    //   h*(1-theta)*Fext_i

    //   // -- Convert the DS into a Lagrangian one.
    //   LagrangianLinearDiagonalDS& d =
    //   static_cast<LagrangianLinearDiagonalDS&> (ds);

    //   // Get state i (previous time step) from Memories -> var. indexed with
    //   "Old" const auto& vold = d.velocityMemory().getSiconosVector(0); //vi

    //   // --- ResiduFree computation ---
    //   // vFree pointer is used to compute and save ResiduFree in this first
    //   step.

    //   // Velocity free and residu. vFree = RESfree (pointer equality !!).
    //   auto& vfree =
    //   *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];
    //   // W is diagonal and contains the inverse of the iteration matrix!
    //   for(unsigned int i=0;i<d.dimension();++i)
    //     vfree(i) = -W(i, i) * vfree(i) + vold(i);

    // }
    // // else if  (dsType == siconos::modeling::Type::NewtonEulerDS)
    // {
    //   // IN to be updated at current time: W, M, q, v, fL
    //   // IN at told: qi,vi,

    //   // Note: indices i/i+1 corresponds to value at the beginning/end of the
    //   time step.
    //   // Index k stands for Newton iteration and thus corresponds to the last
    //   computed
    //   // value, ie the one saved in the DynamicalSystem.
    //   // "i" values are saved in memory vectors.

    //   // vFree = v_k,i+1 - W^{-1} ResiduFree
    //   // with
    //   // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*forces(t,v_k,i+1,
    //   q_k,i+1)
    //   //                                        -
    //   h*(1-theta)*forces(ti,vi,qi)

    //   // -- Convert the DS into a NewtonEuler one.
    //   NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
    //   // --- ResiduFree computation ---
    //   // ResFree = M(v-vold) - h*[theta*forces(t) + (1-theta)*forces(told)]
    //   //
    //   // vFree pointer is used to compute and save ResiduFree in this first
    //   step. auto& residuFree =
    //   *ds_work_vectors[siconos::integrators::MoreauJeanOSI::RESIDU_FREE];
    //   auto& vfree =
    //   *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

    //   vfree = residuFree;

    //   // -- Update W --
    //   // Note: during computeW, mass and jacobians of forces will be
    //   computed/
    //   //SimpleMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;
    //   computeW(t, d, W);
    //   const auto& v = *d.twist(); // v = v_k,i+1

    //   // -- vfree =  v - W^{-1} ResiduFree --
    //   // At this point vfree = residuFree
    //   // -> Solve WX = vfree and set vfree = X
    //   //    std::cout<<"siconos::integrators::MoreauJeanOSI::computeFreeState
    //   residu free"<<endl;
    //   //    vfree->display();
    //   DEBUG_EXPR(residuFree.display(););

    //   W.Solve(vfree);
    //   //    std::cout<<"siconos::integrators::MoreauJeanOSI::computeFreeState
    //   -WRfree"<<endl;
    //   //    vfree->display();
    //   //    siconos::algebra::scal(h,*vfree,*vfree);
    //   // -> compute real vfree
    //   vfree *= -1.0;
    //   DEBUG_EXPR(vfree.display(););
    //   vfree += v;
    //   DEBUG_EXPR(vfree.display(););
    // }
    // else
    //   THROW_EXCEPTION("siconos::integrators::MoreauJeanOSI::computeFreeState
    //   - not yet implemented for Dynamical system of type: " +
    //   siconos::types::type_name(ds));
  }
  DEBUG_END("siconos::integrators::MoreauJeanOSI::computeFreeState()\n");
}

void siconos::integrators::MoreauJeanOSI::prepareNewtonIteration(double time) {
  DEBUG_BEGIN(
      " siconos::integrators::MoreauJeanOSI::prepareNewtonIteration(double "
      "time)\n");
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend;
       ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &sods = *(std::static_pointer_cast<siconos::modeling::SecondOrderDS>(
        _dynamicalSystemsGraph->bundle(*dsi)));
    computeW(time, sods, *_dynamicalSystemsGraph->properties(*dsi).W);
  }

  if (!_explicitNewtonEulerDSOperators) {
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices();
         dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;

      std::shared_ptr<siconos::modeling::DynamicalSystem> ds =
          _dynamicalSystemsGraph->bundle(*dsi);

      //  VA <2016-04-19 Tue> We compute T to be consistent with the Jacobian
      //   at the beginning of the Newton iteration and not at the end
      auto dsType = siconos::types::type_value(*ds);
      if (dsType == siconos::modeling::Type::NewtonEulerDS) {
        auto d = std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds);
        siconos::modeling::computeT(d->q(), d->T());
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
  auto &inter_work_block = *indexSet.properties(vertex_inter).workBlockVectors;

  auto &osnsp_rhs =
      *(*indexSet.properties(vertex_inter)
             .workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS];

  std::shared_ptr<siconos::algebra::BlockVector> Xfree =
      inter_work_block[siconos::integrators::MoreauJeanOSI::xfree];
  assert(Xfree);
  DEBUG_EXPR(Xfree->display(););

  // 1 - product H Xfree
  auto &H = *inter.relation()->H();
  siconos::algebra::prod(H, *Xfree, osnsp_rhs, true);

  // 2 -  compute additional terms for ScleronomousR and CompliantLinearTIR

  // Get relation and non smooth law types
  assert(inter.relation());
  auto relationType = inter.relation()->getType();
  auto relationSubType = inter.relation()->getSubType();

  if ((relationType == siconos::modeling::RelationType::Lagrangian) &&
      (relationSubType != siconos::modeling::RelationSubType::ScleronomousR)) {
    auto sizeY = inter.nonSmoothLaw()->size();

    // Index _selected_coordinates(8);
    _selected_coordinates[0] = 0;
    _selected_coordinates[1] = sizeY;
    _selected_coordinates[2] = 0;
    _selected_coordinates[3] = sizeY;
    _selected_coordinates[4] = 0;
    _selected_coordinates[5] = sizeY;
    _selected_coordinates[6] = 0;
    _selected_coordinates[7] = sizeY;

    auto &DSlink = inter.linkToDSVariables();
    // For the relation of type LagrangianRheonomousR
    if (relationSubType == siconos::modeling::RelationSubType::RheonomousR) {
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() ==
          osnsp) {
        std::static_pointer_cast<siconos::modeling::LagrangianRheonomousR>(
            inter.relation())
            ->computehDot(simulation()->getTkp1(),
                          *DSlink[siconos::modeling::LagrangianR::q0],
                          *DSlink[siconos::modeling::LagrangianR::z]);
        auto ID =
            std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeY);
        ID->eye();
        // This should be optimized -- vacary
        siconos::algebra::subprod(
            *ID,
            *(std::static_pointer_cast<
                  siconos::modeling::LagrangianRheonomousR>(inter.relation())
                  ->hDot()),
            osnsp_rhs, _selected_coordinates, false);  // y += hDot
      } else
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanOSI::computeFreeOutput not yet "
            "implemented for "
            "SICONOS_OSNSP ");
    }
    if (relationSubType ==
        siconos::modeling::RelationSubType::CompliantLinearTIR) {
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() ==
          osnsp) {
        auto &C = *inter.relation()->C();
        double h = _simulation->timeStep();
        osnsp_rhs *= h * _theta;

        /* we have to check that the value are at the beginnning of the time
         * step */
        // + C q_k
        siconos::algebra::subprod(C,
                                  *DSlink[siconos::modeling::LagrangianR::q0],
                                  osnsp_rhs, _selected_coordinates, false);
        // + h(1-_theta)v_k

        *DSlink[siconos::modeling::LagrangianR::q1] *= (1 - _theta) * h;
        siconos::algebra::subprod(C,
                                  *DSlink[siconos::modeling::LagrangianR::q1],
                                  osnsp_rhs, _selected_coordinates, false);

        if (std::static_pointer_cast<
                siconos::modeling::LagrangianCompliantLinearTIR>(
                inter.relation())
                ->e()) {
          auto &e = *std::static_pointer_cast<
                         siconos::modeling::LagrangianCompliantLinearTIR>(
                         inter.relation())
                         ->e();
          osnsp_rhs += e;
        }
      } else
        THROW_EXCEPTION(
            "siconos::integrators::MoreauJeanOSI::computeFreeOutput not yet "
            "implemented for "
            "SICONOS_OSNSP ");
    }
    DEBUG_EXPR(osnsp_rhs.display(););
  }

  // 3 - add part due to NonSmoothLaw
  if (inter.relation()->getType() ==
          siconos::modeling::RelationType::Lagrangian ||
      inter.relation()->getType() ==
          siconos::modeling::RelationType::NewtonEuler) {
    _NSLEffectOnFreeOutput nslEffectOnFreeOutput(
        *osnsp, inter, indexSet.properties(vertex_inter));
    inter.nonSmoothLaw()->accept(nslEffectOnFreeOutput);
  }
  DEBUG_EXPR(osnsp_rhs.display(););

  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::computeFreeOutput("
      "InteractionsGraph::VDescriptor& "
      "vertex_inter, siconos::nonsmooth_formulations::OneStepNSProblem* "
      "osnsp)\n");
}

void siconos::integrators::MoreauJeanOSI::integrate(double &tinit, double &tend,
                                                    double &tout,
                                                    int &notUsed) {
  // Last parameter is not used (required for LsodarOSI but not for
  // MoreauJeanOSI).

  double h = tend - tinit;
  tout = tend;

  std::shared_ptr<siconos::algebra::SiconosMatrix> W;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend;
       ++dsi) {
    if (!checkOSI(dsi)) continue;
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds =
        _dynamicalSystemsGraph->bundle(*dsi);

    W = _dynamicalSystemsGraph->properties(*dsi).W;
    auto dsType = siconos::types::type_value(*ds);

    if (dsType == siconos::modeling::Type::LagrangianLinearTIDS) {
      // get the ds
      auto d =
          std::static_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds);
      // get velocity pointers for current time step
      auto &v = *d->velocity();
      // get q and velocity pointers for previous time step
      const auto &vold = d->velocityMemory().getSiconosVector(0);
      const auto &qold = d->qMemory().getSiconosVector(0);
      // get p pointer

      auto &p = *d->p(1);

      // velocity computation :
      //
      // v = vi + W^{-1}[ -h*C*vi - h*h*theta*K*vi - h*K*qi + h*theta*Fext(t) +
      // h*(1-theta)
      // * Fext(ti) ] + W^{-1}*pi+1
      //

      v = p;

      double coeff;
      // -- No need to update W --
      auto C = d->C();
      if (C) siconos::algebra::prod(-h, *C, vold, v, false);  // v += -h*C*vi

      auto K = d->K();
      if (K) {
        coeff = -h * h * _theta;
        siconos::algebra::prod(coeff, *K, vold, v,
                               false);                   // v += -h^2*theta*K*vi
        siconos::algebra::prod(-h, *K, qold, v, false);  // v += -h*K*qi
      }

      std::shared_ptr<siconos::algebra::SiconosVector> Fext = d->fExt();
      if (Fext) {
        // computes Fext(ti)
        d->computeFExt(tinit);
        coeff = h * (1 - _theta);
        siconos::algebra::scal(coeff, *Fext, v,
                               false);  // v += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(tout);
        coeff = h * _theta;
        siconos::algebra::scal(coeff, *Fext, v,
                               false);  // v += h*theta * fext(ti+1)
      }
      // -> Solve WX = v and set v = X
      W->Solve(v);
      v += vold;
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanOSI::integrate - not yet "
          "implemented for Dynamical "
          "system of type :" +
          siconos::types::type_name(*ds));
  }
}

void siconos::integrators::MoreauJeanOSI::updatePosition(
    siconos::modeling::DynamicalSystem &ds) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::updatePosition(std::shared_ptr<"
      "siconos::modeling::"
      "DynamicalSystem> ds)\n");

  double h = _simulation->timeStep();

  auto dsType = siconos::types::type_value(ds);

  // 1 - Lagrangian Systems
  if (dsType == siconos::modeling::Type::LagrangianDS ||
      dsType == siconos::modeling::Type::LagrangianLinearTIDS ||
      dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
    // get dynamical system
    auto &d = static_cast<siconos::modeling::LagrangianDS &>(ds);

    // Compute q
    auto &v = *d.velocity();
    auto &q = *d.q();
    //  -> get previous time step state
    const auto &vold = d.velocityMemory().getSiconosVector(0);
    const auto &qold = d.qMemory().getSiconosVector(0);
    // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    double coeff = h * _theta;
    siconos::algebra::scal(coeff, v, q);  // q = h*theta*v
    coeff = h * (1 - _theta);
    siconos::algebra::scal(coeff, vold, q, false);  // q += h(1-theta)*vold
    q += qold;
  } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
    // Old Version with projection
    //  NewtonEulerDS& d = static_cast<NewtonEulerDS&> (ds);
    // auto &v = *d.twist();
    // DEBUG_EXPR(d.display());
    // //compute q
    // //first step consists in computing  \dot q.
    // //second step consists in updating q.
    // //
    // SiconosMatrix& T = *d.T();
    // auto& dotq = *d.dotq();
    // DEBUG_EXPR(v.display());
    // siconos::algebra::prod(T, v, dotq, true);
    // DEBUG_EXPR(dotq.display());
    // auto& q = *d.q();
    // //  -> get previous time step state
    // auto& dotqold = *d.dotqMemory()->getSiconosVector(0);
    // DEBUG_EXPR(dotqold.display());
    // // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
    // double coeff = h * _theta;
    // siconos::algebra::scal(coeff, dotq, q) ; // q = h*theta*v
    // coeff = h * (1 - _theta);
    // siconos::algebra::scal(coeff, dotqold, q, false); // q += h(1-theta)*vold
    // auto& qold = *d.qMemory()->getSiconosVector(0);
    // DEBUG_EXPR(qold.display());
    // q += qold;   // q += qold
    // DEBUG_PRINT("new q before normalizing\n");
    // DEBUG_EXPR(q.display());
    // //q[3:6] must be normalized
    // d.normalizeq();
    // DEBUG_PRINT("new q after normalizing\n");
    // DEBUG_EXPR(q.display());

    // get dynamical system
    auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);

    const auto &qold = d.qMemory().getSiconosVector(0);
    const auto &vold = d.twistMemory().getSiconosVector(0);
    auto &q = *d.q();
    auto &v = *d.twist();

    auto velocityIncrement =
        std::make_shared<siconos::algebra::SiconosVector>(v.size());
    double coeff = h * _theta;
    siconos::algebra::scal(
        coeff, v, *velocityIncrement);  //  velocityIncrement= h*theta*v
    coeff = h * (1 - _theta);
    siconos::algebra::scal(coeff, vold, *velocityIncrement,
                           false);  // velocityIncrement += h(1-theta)*vold
    DEBUG_EXPR(velocityIncrement->display());

    q.setValue(0, velocityIncrement->getValue(0));
    q.setValue(1, velocityIncrement->getValue(1));
    q.setValue(2, velocityIncrement->getValue(2));
    siconos::geometry::quaternionFromTwistVector(*velocityIncrement, q);

    DEBUG_EXPR(q.display());
    siconos::geometry::compositionLawLieGroup(qold, q);
    DEBUG_EXPR(q.display());
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::updatePosition(std::shared_ptr<"
      "siconos::modeling::"
      "DynamicalSystem> ds)\n");
}

void siconos::integrators::MoreauJeanOSI::updateState(const unsigned int) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanOSI::updateState(const unsigned int "
      ")\n");

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend;
       ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto &ds = *_dynamicalSystemsGraph->bundle(*dsi);

    auto &ds_work_vectors =
        *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    auto &W = *_dynamicalSystemsGraph->properties(*dsi).W;
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
      auto &vfree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      //    auto *vfree = d.velocityFree();
      auto &v = *d.velocity();
      auto baux = dsType == siconos::modeling::Type::LagrangianDS && useRCC &&
                  _simulation->relativeConvergenceCriterionHeld();

      if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
        assert(((d.p(_levelMaxForInput)).get()) &&
               " siconos::integrators::MoreauJeanOSI::updateState() "
               "*d.p(_levelMaxForInput) "
               "== nullptr.");
        v = *d.p(_levelMaxForInput);  // v = p
        if (d.boundaryConditions()) {
          for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
            v.setValue(itindex, 0.0);
          }
        }

        if (dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
          for (unsigned int i = 0; i < d.dimension(); ++i)
            v(i) = vfree(i) + W(i, i) * v(i);
        } else {
          W.Solve(v);
          v += vfree;
        }
      } else {
        v = vfree;
      }
      DEBUG_EXPR(v.display());

      if (d.boundaryConditions()) {
        int bc = 0;
        auto columntmp =
            std::make_shared<siconos::algebra::SiconosVector>(ds.dimension());

        for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
          _dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions->getCol(
              bc, *columntmp);
          /*\warning we assume that W is symmetric in the Lagrangian case*/

          double value = -siconos::algebra::inner_prod(*columntmp, v);
          if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
            value += (d.p(_levelMaxForInput))->getValue(itindex);
          }
          /* \warning the computation of reactionToBoundaryConditions take into
             account the contact impulse but not the external and internal
             forces. A complete computation of the residu should be better */
          d.reactionToBoundaryConditions()->setValue(bc, value);
          bc++;
        }
      }

      auto &q = *d.q();
      auto &local_buffer =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::BUFFER];
      // Save value of q in stateTmp for future convergence computation
      if (baux) local_buffer = q;

      updatePosition(ds);

      if (baux) {
        double ds_norm_ref =
            1. + ds.x0()->norm2();  // Should we save this in the graph?
        local_buffer -= q;
        double aux = (local_buffer.norm2()) / ds_norm_ref;
        if (aux > RelativeTol)
          _simulation->setRelativeConvergenceCriterionHeld(false);
      }
    } else if (dsType == siconos::modeling::Type::NewtonEulerDS) {
      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::updateState(const unsigned "
          "int), "
          "dsType == "
          "siconos::modeling::Type::NewtonEulerDS \n");

      // get dynamical system
      auto &d = static_cast<siconos::modeling::NewtonEulerDS &>(ds);
      auto &v = *d.twist();
      // DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updateState()\n ")
      // DEBUG_EXPR(d.display());
      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updateState() prev v\n")
      DEBUG_EXPR(v.display());

      // failure on bullet sims
      // d.p(_levelMaxForInput) is checked in next condition
      // assert(((d.p(_levelMaxForInput)).get()) &&
      //       " siconos::integrators::MoreauJeanOSI::updateState()
      //       *d.p(_levelMaxForInput)
      //       == nullptr.");

      auto &vfree =
          *ds_work_vectors[siconos::integrators::MoreauJeanOSI::VFREE];

      if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
        /*d.p has been fill by the Relation->computeInput, it contains
          B \lambda _{k+1}*/
        v = *d.p(_levelMaxForInput);  // v = p
        if (d.boundaryConditions())
          for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
            v.setValue(itindex, 0.0);
          }

        _dynamicalSystemsGraph->properties(*dsi).W->Solve(v);

        DEBUG_EXPR(d.p(_levelMaxForInput)->display());
        DEBUG_PRINT(
            "siconos::integrators::MoreauJeanOSI::updatestate W CT lambda\n");
        DEBUG_EXPR(v.display());
        v += vfree;
      } else
        v = vfree;

      DEBUG_PRINT(
          "siconos::integrators::MoreauJeanOSI::updatestate work free\n");
      DEBUG_EXPR(vfree.display());
      DEBUG_PRINT("siconos::integrators::MoreauJeanOSI::updatestate new v\n");
      DEBUG_EXPR(v.display());

      if (d.boundaryConditions()) {
        int bc = 0;
        auto columntmp =
            std::make_shared<siconos::algebra::SiconosVector>(ds.dimension());

        for (const auto itindex : d.boundaryConditions()->velocityIndices()) {
          _dynamicalSystemsGraph->properties(*dsi).WBoundaryConditions->getCol(
              bc, *columntmp);
          /*\warning we assume that W is symmetric in the Lagrangian case*/
          double value = -siconos::algebra::inner_prod(*columntmp, v);
          if (d.p(_levelMaxForInput) && d.p(_levelMaxForInput)->size() > 0) {
            value += (d.p(_levelMaxForInput))->getValue(itindex);
          }
          /* \warning the computation of reactionToBoundaryConditions take into
             account the contact impulse but not the external and internal
             forces. A complete computation of the residu should be better */
          d.reactionToBoundaryConditions()->setValue(bc, value);
          bc++;
        }
      }

      updatePosition(ds);
    } else
      THROW_EXCEPTION(
          "siconos::integrators::MoreauJeanOSI::updateState - not yet "
          "implemented for "
          "Dynamical system of type: " +
          siconos::types::type_name(ds));
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanOSI::updateState(const unsigned int)\n");
}

bool siconos::integrators::MoreauJeanOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction>"
      " inter, "
      "unsigned int i)\n");

  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity

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
  DEBUG_EXPR_WE(if (y <= _constraintActivationThreshold) std::cout
                    << "siconos::integrators::MoreauJeanOSI::"
                       "addInteractionInIndexSet ACTIVATE."
                    << y << "<= " << _constraintActivationThreshold
                    << std::endl;
                ;);
  return (y <= _constraintActivationThreshold);
}

bool siconos::integrators::MoreauJeanOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  return !(addInteractionInIndexSet(inter, i));
}

void siconos::integrators::MoreauJeanOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanOSI OSI display ======" << std::endl;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  if (_dynamicalSystemsGraph) {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices();
         dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds =
          _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------" << std::endl;
      std::cout << "--> W of dynamical system number " << ds->number() << ": "
                << std::endl;
      if (_dynamicalSystemsGraph->properties(*dsi).W)
        _dynamicalSystemsGraph->properties(*dsi).W->display();
      else
        std::cout << "-> nullptr" << std::endl;
      std::cout << "--> and corresponding theta is: " << _theta << std::endl;
    }
  }
  std::cout << "================================" << std::endl;
}
