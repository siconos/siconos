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
#include "MoreauJeanCombinedProjectionOSI.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "Relation.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"
// #include "Simulation.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::integrators::MoreauJeanCombinedProjectionOSI::MoreauJeanCombinedProjectionOSI(
    double theta)
    : MoreauJeanOSI(theta) {
  _levelMinForInput = 0;
  //_integratorType = IntegratorType::MOREAUDIRECTPROJECTIONOSI;
}

void siconos::integrators::MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForDS( "
      "double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds) \n");

  MoreauJeanOSI::initializeWorkVectorsForDS(t, ds);

  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
  auto& workVectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    workVectors[MoreauJeanOSI::QTMP] =
        std::make_shared<siconos::algebra::SiconosVector>(d->dimension());
  } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
    workVectors[MoreauJeanOSI::QTMP] =
        std::make_shared<siconos::algebra::SiconosVector>(d->getqDim());
  } else {
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanCombinedProjectionOSI::initialize() - DS not of the "
        "right type");
  }
  for (unsigned int k = _levelMinForInput; k < _levelMaxForInput + 1; k++) {
    ds->initializeNonSmoothInput(k);
  }

  DEBUG_END(
      "siconos::integrators::MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForDS( "
      "double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds) \n");
}

void siconos::integrators::MoreauJeanCombinedProjectionOSI::
    initializeWorkVectorsForInteraction(siconos::modeling::Interaction& inter,
                                        siconos::graphs::InteractionProperties& interProp,
                                        siconos::graphs::DynamicalSystemsGraph& DSG) {
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanCombinedProjectionOSI::"
      "initializeWorkVectorsForInteraction(siconos::modeling:"
      ":Interaction&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");

  MoreauJeanOSI::initializeWorkVectorsForInteraction(inter, interProp, DSG);

  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);
  auto& DSlink = inter.linkToDSVariables();

  auto& relation = *inter.relation();
  auto relationType = relation.getType();

  unsigned int p0 = 0;
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    p0 = tools::enum_to_index(modeling::LagrangianR::WorkDS::p0);
  } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    p0 = siconos::tools::enum_to_index(siconos::modeling::NewtonEulerR::WorkDS::p0);
  }
  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!DSlink[p0]) || (DSlink[p0]->numberOfBlocks() != 2))
      DSlink[p0] = std::make_shared<siconos::algebra::BlockVector>(2);
  } else {
    if ((!DSlink[p0]) || (DSlink[p0]->numberOfBlocks() != 1))
      DSlink[p0] = std::make_shared<siconos::algebra::BlockVector>(1);
  }

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1);
      DSlink[p0]->setVectorPtr(0, lds.p(0));
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      auto& neds = *std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds1);
      DSlink[p0]->setVectorPtr(0, neds.p(0));
    }
  }
  DEBUG_PRINTF("ds1->number() %i\n", ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n", ds2->number());

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      if (relationType == siconos::modeling::RelationType::Lagrangian) {
        auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2);
        DSlink[p0]->setVectorPtr(1, lds.p(0));
      } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
        auto& neds = *std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds2);
        DSlink[p0]->setVectorPtr(1, neds.p(0));
      }
    }
  }

  DEBUG_END(
      "siconos::integrators::MoreauJeanCombinedProjectionOSI::"
      "initializeWorkVectorsForInteraction(siconos::modeling:"
      ":Interaction&inter, siconos::graphs::InteractionProperties& interProp, "
      "siconos::graphs::DynamicalSystemsGraph & DSG)\n");
}

bool siconos::integrators::MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(i == 1 || i == 2);
  // double h = _simulation->timeStep();
  if (i == 1)  // index set for resolution at the velocity
  {
    auto y = (inter->y(0))->getValue(0);  // y(0) is the position
    DEBUG_PRINTF(
        "siconos::integrators::MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet "
        "yref=%e \n",
        y);
    DEBUG_EXPR(
        if (y <= 0) printf(
            "siconos::integrators::MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet "
            "ACTIVATE in indexSet level = %i.\n",
            i);)
    return (y <= 0);
  } else if (i == 2)  //  special index for the projection
  {
    DEBUG_EXPR(
        auto lambda = (inter->lambda(1))->getValue(0);
        // lambda(1) is the contact impulse for MoreauJeanOSI scheme
        printf("siconos::integrators::MoreauJeanCombinedProjectionOSI::"
               "addInteractionInIndexSet lambdaref=%e \n",
               lambda);
        if (lambda > 0) printf(
            "siconos::integrators::MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet "
            "ACTIVATE in indexSet level = %i.\n",
            i);)
    //    return (lambda > 0);
    return true;
  } else {
    return false;
  }
}

bool siconos::integrators::MoreauJeanCombinedProjectionOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(0);
  return (0);
}
