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
#include "MoreauJeanDirectProjectionOSI.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
//#define STANDARD_ACTIVATION
#define FIRSTWAY_ACTIVATION
//#define SECONDWAY_ACTIVATION
//#define QFREE_ACTIVATION

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
//#define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::integrators::MoreauJeanDirectProjectionOSI::MoreauJeanDirectProjectionOSI(
    double theta)
    : MoreauJeanOSI(theta)
{
  _levelMinForInput = 0;
  _integratorType = IntegratorType::MOREAUDIRECTPROJECTIONOSI;
}

siconos::integrators::MoreauJeanDirectProjectionOSI::MoreauJeanDirectProjectionOSI(
    double theta, double gamma)
    : MoreauJeanOSI(theta, gamma)
{
  _levelMinForInput = 0;
  _integratorType = IntegratorType::MOREAUDIRECTPROJECTIONOSI;
}

void siconos::integrators::MoreauJeanDirectProjectionOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::initializeWorkVectorsForDS( "
      "double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) \n");
  MoreauJeanOSI::initializeWorkVectorsForDS(t, ds);
  const auto& dsv = _dynamicalSystemsGraph->descriptor(ds);
  auto& workVectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    workVectors[MoreauJeanOSI::QTMP] =
        std::make_shared<siconos::algebra::SiconosVector>(d->dimension());
  }
  else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
    workVectors[MoreauJeanOSI::QTMP] =
        std::make_shared<siconos::algebra::SiconosVector>(d->getqDim());
  }
  else {
    THROW_EXCEPTION(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::initialize() - DS not of the "
        "right type");
  }
  for (unsigned int k = _levelMinForInput; k < _levelMaxForInput + 1; k++) {
    DEBUG_PRINTF("ds->initializeNonSmoothInput(%i)\n", k);
    ds->initializeNonSmoothInput(k);
    DEBUG_EXPR_WE(auto d = std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds);
                  if (d->p(k)) std::cout << "d->p(" << k << " ) exists\n";);
  }
  DEBUG_END(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::initializeWorkVectorsForDS( "
      "double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) \n");
}

void siconos::integrators::MoreauJeanDirectProjectionOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG)
{
  DEBUG_BEGIN(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::"
      "initializeWorkVectorsForInteraction(siconos::modeling::Interaction&inter, "
      "InteractionProperties& interProp, siconos::graphs::DynamicalSystemsGraph & DSG)\n");

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
    p0 = siconos::modeling::LagrangianR::p0;
  }
  else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    p0 = siconos::modeling::NewtonEulerR::p0;
  }

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!DSlink[p0]) || (DSlink[p0]->numberOfBlocks() != 2))
      DSlink[p0] = std::make_shared<siconos::algebra::BlockVector>(2);
  }
  else {
    if ((!DSlink[p0]) || (DSlink[p0]->numberOfBlocks() != 1))
      DSlink[p0] = std::make_shared<siconos::algebra::BlockVector>(1);
  }

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1);
      DSlink[p0]->setVectorPtr(0, lds.p(0));
    }
    else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
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
      }
      else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
        auto& neds = *std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds2);
        DSlink[p0]->setVectorPtr(1, neds.p(0));
      }
    }
  }

  DEBUG_END(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::"
      "initializeWorkVectorsForInteraction(siconos::modeling::Interaction&inter, "
      "InteractionProperties& interProp, siconos::graphs::DynamicalSystemsGraph & DSG)\n");
}

void siconos::integrators::MoreauJeanDirectProjectionOSI::computeFreeState()
{
  MoreauJeanOSI::computeFreeState();
}

#ifdef STANDARD_ACTIVATION
bool siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  return MoreauJeanOSI::addInteractionInIndexSet(inter, i);
}

bool siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  return MoreauJeanOSI::removeInteractionFromIndexSet(inter, i);
}
#endif

#ifdef FIRSTWAY_ACTIVATION
bool siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  assert(i == 1);
  auto h = _simulation->timeStep();
  auto y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  auto yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  if (_useGamma) {
    gamma = _gamma;
  }
  DEBUG_PRINTF(
      "\nsiconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
      "inter->number() = %i\n",
      inter->number());
  DEBUG_EXPR(inter->display(););
  DEBUG_PRINTF("MoreauJeanOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n",
               y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
      "yref=%e, yDot=%e.\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet  "
      "_activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
      _activateYPosThreshold, _activateYVelThreshold);

  assert(!std::isnan(y));
#ifdef DEBUG_MESSAGES
  if (y <= _activateYPosThreshold)
    DEBUG_PRINT(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
        "ACTIVATE.\n");
#endif
  return (y <= _activateYPosThreshold);
}

bool siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)

{
  assert(i == 1);
  auto h = _simulation->timeStep();
  auto y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  auto yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  if (_useGamma) {
    gamma = _gamma;
  }
  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet "
      "yref=%e, yDot=%e .\n",
      y, yDot);
  y += gamma * h * yDot;

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet  "
      "_deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
      _deactivateYPosThreshold, _deactivateYVelThreshold);

  assert(!std::isnan(y));
#ifdef DEBUG_MESSAGES
  if (y > _deactivateYPosThreshold && yDot >= _deactivateYVelThreshold)
    DEBUG_PRINT(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet "
        "DEACTIVATE.\n");
#endif
  return (y > _deactivateYPosThreshold && yDot >= _deactivateYVelThreshold);
}

#endif

#ifdef SECONDWAY_ACTIVATION
bool siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  assert(i == 1);
  auto y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
#ifdef DEBUG_MESSAGES
  auto yDot = (inter->y(i))->getValue(0);  // for i=1 y(i) is the velocity
#endif
  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
      "inter->number() = %i\n",
      inter->number());
  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
      "yref=%e, yDot=%e.\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet  "
      "_activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
      _activateYPosThreshold, _activateYVelThreshold);

  assert(!std::isnan(y));

  if (y <= _activateYPosThreshold)
    DEBUG_PRINT(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
        "ACTIVATE.\n");
  return (y <= _activateYPosThreshold);
}

bool siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)

{
  assert(i == 1);
  auto y = (inter->y(i - 1))->getValue(0);        // for i=1 y(i-1) is the position
  auto yDot = (inter->y(i))->getValue(0);         // for i=1 y(i) is the velocity
  auto lambda = (inter->lambda(i))->getValue(0);  // for i=1 y(i) is the velocity

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet "
      "yref=%e, yDot=%e .\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet  "
      "_deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
      _deactivateYPosThreshold, _deactivateYVelThreshold);

  assert(!std::isnan(y));
  if (y > _deactivateYPosThreshold && lambda <= _deactivateYVelThreshold)
    DEBUG_PRINT(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet "
        "DEACTIVATE.\n");
  return (y > _deactivateYPosThreshold && lambda <= _deactivateYVelThreshold);
}

#endif

#ifdef QFREE_ACTIVATION
bool siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  assert(i == 1);
  auto y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
#ifdef DEBUG_MESSAGES
  auto yDot = (inter->y(i))->getValue(0);  // for i=1 y(i) is the velocity
#endif
  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
      "yref=%e, yDot=%e.\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet  "
      "_activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
      _activateYPosThreshold, _activateYVelThreshold);

  assert(!std::isnan(y));

  if (y <= _activateYPosThreshold)
    DEBUG_PRINT(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::addInteractionInIndexSet "
        "ACTIVATE.\n");
  return (y <= _activateYPosThreshold);
}

bool siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)

{
  assert(i == 1);
  auto y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  auto yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet "
      "yref=%e, yDot=%e .\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet  "
      "_deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
      _deactivateYPosThreshold, _deactivateYVelThreshold);

  assert(!std::isnan(y));
  if (y > _deactivateYPosThreshold)
    DEBUG_PRINT(
        "siconos::integrators::MoreauJeanDirectProjectionOSI::removeInteractionFromIndexSet "
        "DEACTIVATE.\n");
  return (y > _deactivateYPosThreshold);
}

#endif
