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

#include "D1MinusLinearOSI.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactNSL.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosMatrixVectorOp.hpp"  // mat-vec prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for subscal
#include "SiconosVisitor.hpp"
#include "SiconosMatrix.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"  // For enum_to_string
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::integrators::D1MinusLinearOSI::_NSLEffectOnFreeOutput::_NSLEffectOnFreeOutput(
    siconos::nonsmooth_formulations::OneStepNSProblem* p,
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionProperties& interProp)
    : _osnsp(p), _inter(inter), _interProp(interProp){};

void siconos::integrators::D1MinusLinearOSI::_NSLEffectOnFreeOutput::visit(
    const siconos::modeling::NewtonImpactNSL& nslaw) const {
  double e = nslaw.e();
  std::size_t nsl_size = _inter->nonSmoothLaw()->size();
  std::vector<std::size_t> subCoord = {0, nsl_size, 0, nsl_size};
  siconos::algebra::SiconosVector& osnsp_rhs =
      *(*_interProp.workVectors)[siconos::integrators::D1MinusLinearOSI::OSNSP_RHS];
  siconos::algebra::subscal(e, osnsp_rhs, osnsp_rhs, subCoord, false);
}

siconos::integrators::D1MinusLinearOSI::D1MinusLinearOSI(Type type)
    : OneStepIntegrator{IntegratorType::D1MINUSLINEAROSI, 2, 0, 2, 1, 2},
      _typeOfD1MinusLinearOSI{type} {}

unsigned int siconos::integrators::D1MinusLinearOSI::numberOfIndexSets() const {
  switch (_typeOfD1MinusLinearOSI) {
    case Type::halfexplicit_acceleration_level:
      return 4;
    case Type::halfexplicit_acceleration_level_full:
      return 4;
    case Type::halfexplicit_velocity_level:
      return 3;
  }
  THROW_EXCEPTION(
      "siconos::integrators::D1MinusLinearOSI::numberOfIndexSet - not implemented for "
      "D1minusLinear of type: " +
      siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
  return 0;
}
void siconos::integrators::D1MinusLinearOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type

  if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    lds->init_generalized_coordinates(2);  // acceleration is required for the ds
    lds->init_inverse_mass();              // invMass required to update post-impact velocity

    ds_work_vectors.resize(siconos::integrators::D1MinusLinearOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    ds_work_vectors[siconos::integrators::D1MinusLinearOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    ds_work_vectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
    // Update dynamical system components (for memory swap).
    lds->computeForces(t, lds->q(), lds->velocity());
    lds->swapInMemory();
  }

  else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
    neds->init_inverse_mass();  // invMass required to update post-impact velocity
    ds_work_vectors.resize(siconos::integrators::D1MinusLinearOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(neds->dimension());
    ds_work_vectors[siconos::integrators::D1MinusLinearOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(neds->dimension());
    ds_work_vectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG] =
        std::make_shared<siconos::algebra::SiconosVector>(neds->dimension());
    // Compute a first value of the forces to store it in _forcesMemory
    neds->computeForces(t, neds->q(), neds->twist());
    neds->swapInMemory();
  } else {
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::initialize - Only implemented for Lagrangian "
        "or NewtonEuler dynamical systems.");
  }

  for (unsigned int k = _levelMinForInput; k < _levelMaxForInput + 1; k++) {
    ds->initializeNonSmoothInput(k);
  }
}

void siconos::integrators::D1MinusLinearOSI::initialize_nonsmooth_problems() {
  auto allOSNSP = _simulation->oneStepNSProblems();  // all OSNSP

  bool isOSNSPinitialized = false;
  switch (_typeOfD1MinusLinearOSI) {
    case Type::halfexplicit_acceleration_level:
      // set evaluation levels (first is of velocity, second of acceleration type)
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->initialize(_simulation);

      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->setIndexSetLevel(2);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->setInputOutputLevel(2);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(_simulation);
      isOSNSPinitialized = true;
      DEBUG_EXPR((*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->display());
      break;
    case Type::halfexplicit_acceleration_level_full:
      // set evaluation levels (first is of velocity, second of acceleration type)
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->initialize(_simulation);

      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->setIndexSetLevel(2);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->setInputOutputLevel(2);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(_simulation);
      isOSNSPinitialized = true;
      break;
    case Type::halfexplicit_velocity_level:
      // set evaluation levels (first is of velocity, second of acceleration type)
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->initialize(_simulation);

      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->setIndexSetLevel(
          1); /** !!! */
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->setInputOutputLevel(2);
      (*allOSNSP)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(_simulation);
      isOSNSPinitialized = true;
      break;
  }

  if (!isOSNSPinitialized) {
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::initialize() - not implemented for type of "
        "D1MinusLinearOSI: " +
        siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
  }
}

void siconos::integrators::D1MinusLinearOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  DEBUG_BEGIN(
      "siconos::integrators::D1MinusLinearOSI::initializeWorkVectorsForInteraction(siconos::"
      "modeling::Interaction&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);
  DEBUG_PRINTF("interaction number %i\n", inter.number());
  auto& DSlink = inter.linkToDSVariables();

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::D1MinusLinearOSI::WORK_INTERACTION_LENGTH);
  }
  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::D1MinusLinearOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();
  auto relationType = relation.getType();

  inter_work[siconos::integrators::D1MinusLinearOSI::OSNSP_RHS] =
      std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2)))) {
    std::cout << "checkOSI(DSG.descriptor(ds1)): " << std::boolalpha
              << checkOSI(DSG.descriptor(ds1)) << std::endl;
    std::cout << "checkOSI(DSG.descriptor(ds2)): " << std::boolalpha
              << checkOSI(DSG.descriptor(ds2)) << std::endl;

    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::initializeWorkVectorsForInteraction. The "
        "implementation is not correct for two different OSI for one interaction");
  }

  /* allocate and set work vectors for the osi */
  auto xfree = siconos::integrators::D1MinusLinearOSI::xfree;
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
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    inter_work_block[xfree]->setVectorPtr(
        0, workVds1[siconos::integrators::D1MinusLinearOSI::FREE]);
  }
  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      inter_work_block[xfree]->setVectorPtr(
          1, workVds2[siconos::integrators::D1MinusLinearOSI::FREE]);
    }
  }

  DEBUG_EXPR(inter_work_block[xfree]->display(););

  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1);
    DSlink[siconos::modeling::LagrangianR::p2] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::p2]->insertPtr(lds.p(2));
    DSlink[siconos::modeling::LagrangianR::q2] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::q2]->insertPtr(lds.acceleration());
  } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
  }

  if (ds1 != ds2) {
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2);
      DSlink[siconos::modeling::LagrangianR::p2]->insertPtr(lds.p(2));
      DSlink[siconos::modeling::LagrangianR::q2]->insertPtr(lds.acceleration());
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    }
  }

  DEBUG_END(
      "siconos::integrators::D1MinusLinearOSI::initializeWorkVectorsForInteraction(siconos::"
      "modeling::Interaction&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");
}

double siconos::integrators::D1MinusLinearOSI::computeResidu() {
  DEBUG_PRINT("\n ******************************************************************\n");
  DEBUG_PRINT(" ******************************************************************\n");
  DEBUG_PRINT(" ******************************************************************\n");

  DEBUG_BEGIN("siconos::integrators::D1MinusLinearOSI::computeResidu()\n");

  DEBUG_PRINTF("nextTime %f\n", _simulation->nextTime());
  DEBUG_PRINTF("startingTime %f\n", _simulation->startingTime());
  DEBUG_PRINTF("time step size %f\n", _simulation->timeStep());

  switch (_typeOfD1MinusLinearOSI) {
    case Type::halfexplicit_acceleration_level:
      DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeResidu()\n");
      return computeResiduHalfExplicitAccelerationLevel();
    case Type::halfexplicit_velocity_level:
      DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeResidu()\n");
      return computeResiduHalfExplicitVelocityLevel();
    case Type::halfexplicit_acceleration_level_full:
      DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeResidu()\n");
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeResidu() - not implemented for type "
          "of "
          "D1MinusLinearOSI: " +
          siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
      // return computeResiduHalfExplicitAccelerationLevelFull(); // obsolete
  }
  THROW_EXCEPTION(
      "siconos::integrators::D1MinusLinearOSI::computeResidu() - not implemented for type of "
      "D1MinusLinearOSI: " +
      siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
  DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeResidu()\n");
  return 1;
}

void siconos::integrators::D1MinusLinearOSI::computeFreeState() {
  DEBUG_BEGIN("siconos::integrators::D1MinusLinearOSI::computeFreeState()\n");
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    /* \warning the following conditional statement should be removed with a MechanicalDS class
     */
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // Lagrangian Systems

      // get left state from memory
      const auto& vold = d->velocityMemory().getSiconosVector(0);  // right limit
      DEBUG_EXPR(vold.display());
      auto& residuFree = *ds_work_vectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
      siconos::algebra::SiconosVector& vfree =
          *d->velocity();  // POINTER CONSTRUCTOR : contains free velocity

      // get right information
      vfree = residuFree;
      DEBUG_EXPR(residuFree.display());
      // d->computeMass();
      // M->resetFactorizationFlags();
      // M->Solve(vfree);
      // DEBUG_EXPR(M->display());

      vfree *= -1.;
      vfree += vold;
      DEBUG_EXPR(vfree.display());
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      // NewtonEuler Systems

      // get left state from memory
      const auto& vold = d->twistMemory().getSiconosVector(0);  // right limit
      DEBUG_EXPR(vold.display());

      // get right information
      // auto M = std::make_shared<siconos::algebra::SiconosMatrix>(*(d->mass()))); // we copy
      // the mass matrix to avoid its factorization;
      auto& vfree = *d->twist();  // POINTER CONSTRUCTOR : contains free velocity
      auto& residuFree = *ds_work_vectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];

      vfree = residuFree;
      DEBUG_EXPR(residuFree.display());

      vfree *= -1.;
      vfree += vold;
      DEBUG_EXPR(vfree.display());
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeResidu - only implemented for "
          "Lagrangian or Newton Euler dynamical systems.");
  }
  DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeFreeState()\n");
}

void siconos::integrators::D1MinusLinearOSI::updateState(const unsigned int) {
  DEBUG_BEGIN("siconos::integrators::D1MinusLinearOSI::updateState(const unsigned int)\n");

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;

    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    /* \warning the following conditional statement should be removed with a MechanicalDS
     * class
     */
    /* Lagrangian DS*/
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto v = d->velocity();

      DEBUG_PRINT("Position and velocity before update\n");
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());

      /* Add the contribution of the impulse if any */
      if (d->p(1)) {
        DEBUG_EXPR(d->p(1)->display());
        /* copy the value of the impulse */
        auto dummy = std::make_shared<siconos::algebra::SiconosVector>(*(d->p(1)));
        /* Compute the velocity jump due to the impulse */
        if (d->inverseMass()) {
          d->update_inverse_mass();
          siconos::algebra::solveInPlace(*(d->inverseMass()), *dummy);
        }
        /* Add the velocity jump to the free velocity */
        *v += *dummy;
      }

      DEBUG_PRINT("Position and velocity after update\n");
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    /*  NewtonEuler Systems */
    else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto v = d->twist();  // POINTER CONSTRUCTOR : contains new velocity
      if (d->p(1)) {
        // Update the velocity
        auto dummy = std::make_shared<siconos::algebra::SiconosVector>(
            *(d->p(1)));  // value = nonsmooth impulse
        if (d->inverseMass()) {
          d->update_inverse_mass();
          siconos::algebra::solveInPlace(*(d->inverseMass()), *dummy);
        }
        *v += *dummy;  // add free velocity

        // update \f$ \dot q \f$
        auto T = d->T();
        auto dotq = d->dotq();
        siconos::algebra::prod(*T, *v, *dotq, true);

        DEBUG_PRINT("\nRIGHT IMPULSE\n");
        DEBUG_EXPR(d->p(1)->display());
      }
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeResidu - only implemented for "
          "Lagrangian or Newton Euler dynamical systems.");
  }

  DEBUG_END("\n siconos::integrators::D1MinusLinearOSI::updateState(const unsigned int )\n");
}

void siconos::integrators::D1MinusLinearOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  DEBUG_PRINT("siconos::integrators::D1MinusLinearOSI::computeFreeOutput(), start\n");
  switch (_typeOfD1MinusLinearOSI) {
    case Type::halfexplicit_acceleration_level:
      computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter, osnsp);
      DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeFreeOutput()\n");
      return;
    case Type::halfexplicit_acceleration_level_full:
      computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter, osnsp);
      DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeFreeOutput()\n");
      return;
    case Type::halfexplicit_velocity_level:
      computeFreeOutputHalfExplicitVelocityLevel(vertex_inter, osnsp);
      DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeFreeOutput()\n");
      return;
  }
  THROW_EXCEPTION(
      "siconos::integrators::D1MinusLinearOSI::computeResidu() - not implemented for type "
      "of "
      "D1MinusLinearOSI: " +
      siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
  DEBUG_END("siconos::integrators::D1MinusLinearOSI::computeFreeOutput()\n");
}

bool siconos::integrators::D1MinusLinearOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_BEGIN("siconos::integrators::D1MinusLinearOSI::addInteractionInIndexSet.\n");
  DEBUG_END("siconos::integrators::D1MinusLinearOSI::addInteractionInIndexSet.\n");
  switch (_typeOfD1MinusLinearOSI) {
    case Type::halfexplicit_acceleration_level:
      return addInteractionInIndexSetHalfExplicitAccelerationLevel(inter, i);
    case Type::halfexplicit_velocity_level:
      return addInteractionInIndexSetHalfExplicitVelocityLevel(inter, i);
    default:
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::addInteractionInIndexSet() - not "
          "implemented "
          "for type of D1MinusLinearOSI: " +
          siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
  }
  return 0;
}

bool siconos::integrators::D1MinusLinearOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_BEGIN("siconos::integrators::D1MinusLinearOSI::removeInteractionFromIndexSet.\n");
  DEBUG_END("siconos::integrators::D1MinusLinearOSI::removeInteractionFromIndexSet.\n");
  switch (_typeOfD1MinusLinearOSI) {
    case Type::halfexplicit_acceleration_level:
      return removeInteractionFromIndexSetHalfExplicitAccelerationLevel(inter, i);
    case Type::halfexplicit_velocity_level:
      return removeInteractionFromIndexSetHalfExplicitVelocityLevel(inter, i);
    default:
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::removeInteractionFromIndexSet() - not "
          "implemented for type of D1MinusLinearOSI: " +
          siconos::tools::enum_to_string(_typeOfD1MinusLinearOSI));
  }
  return 0;
}
