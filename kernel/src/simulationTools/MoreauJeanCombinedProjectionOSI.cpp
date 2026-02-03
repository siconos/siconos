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

#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "Relation.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosVector.hpp"
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

  // const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
  // auto& workVectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  // if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
  //   workVectors[MoreauJeanOSI::QTMP] =
  //       std::make_shared<siconos::algebra::SiconosVector>(d->dimension());
  // } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
  //   workVectors[MoreauJeanOSI::QTMP] =
  //       std::make_shared<siconos::algebra::SiconosVector>(d->getqDim());
  // } else {
  //   THROW_EXCEPTION(
  //       "siconos::integrators::MoreauJeanCombinedProjectionOSI::initialize() - DS not of the
  //       " "right type");
  // }
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
      "initializeWorkVectorsForInteraction(siconos::modeling::Interaction&inter, "
      "InteractionProperties& interProp, siconos::graphs::DynamicalSystemsGraph & DSG)\n");

  MoreauJeanOSI::initializeWorkVectorsForInteraction(inter, interProp, DSG);

  auto ds1 = interProp.source;
  assert(ds1);
  auto is_ds1_integrated_by_this_osi = checkOSI(DSG.descriptor(ds1));

  auto ds2 = interProp.target;
  assert(ds2);
  auto use_two_ds = ds1 != ds2;
  bool is_ds2_integrated_by_this_osi = use_two_ds && checkOSI(DSG.descriptor(ds2));

  auto relationType = inter.relation()->getType();

  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    if (is_ds1_integrated_by_this_osi) {
      // SecondOrderDS to include Lagrangian and LagrangianSparse
      auto& lds1 = *std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds1);
      // inter.link[p0][0] = ds1->p[0]
      inter.append_to_dynamical_systems_variables(siconos::modeling::LagrangianR::ds_var::p0,
                                                  0, lds1.p(0));
    }
    if (is_ds2_integrated_by_this_osi) {
      auto& lds2 = *std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds2);
      // inter.link[p0][1] = ds2->p[0]
      inter.append_to_dynamical_systems_variables(siconos::modeling::LagrangianR::ds_var::p0,
                                                  1, lds2.p(0));
    }
  }

  else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    if (is_ds1_integrated_by_this_osi) {
      auto& lds1 = *std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds1);
      // inter.link[p0][0] = ds1->p[0]
      inter.append_to_dynamical_systems_variables(siconos::modeling::NewtonEulerR::ds_var::p0,
                                                  0, lds1.p(0));
    }
    if (is_ds2_integrated_by_this_osi) {
      auto& lds2 = *std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds2);
      // inter.link[p0][1] = ds2->p[0]
      inter.append_to_dynamical_systems_variables(siconos::modeling::NewtonEulerR::ds_var::p0,
                                                  1, lds2.p(0));
    }
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::"
      "initializeWorkVectorsForInteraction(siconos::modeling::Interaction&inter, "
      "InteractionProperties& interProp, siconos::graphs::DynamicalSystemsGraph & DSG)\n");
}

bool siconos::integrators::MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(i == 1 || i == 2);
  // double h = _simulation->timeStep();
  if (i == 1)  // index set for resolution at the velocity
  {
    auto y = (*(inter->y(0)))(0);  // y(0) is the position
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
        auto lambda = (*(inter->lambda(1)))(0);
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
