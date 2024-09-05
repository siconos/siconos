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

#include "MBTB_MoreauJeanOSI.hpp"

#include "Interaction.hpp"
#include "SiconosVector.hpp"
#define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

// #define STANDARD_ACTIVATION
#define FIRSTWAY_ACTIVATION

#ifdef STANDARD_ACTIVATION
bool siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(i == 1);
  double h = _simulation->timeStep();
  double y =
      (inter->y(i - 1))->getValxxue(0);      // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0);  // for i=1 y(i) is the velocity

  double gamma = 1.0 / 2.0;
  if (_useGamma) {
    gamma = _gamma;
  }
  DEBUG_PRINTF(
      "siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet "
      "yref=%e, yDot=%e, "
      "y_estimated=%e.\n",
      y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!std::isnan(y));
  if (y <= 0)
    DEBUG_PRINT(
        "siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet "
        "ACTIVATE.\n");
  return (y <= 0);
}
bool siconos::mechanisms::MBTB_MoreauJeanOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  if (_useGamma) {
    gamma = _gamma;
  }
  DEBUG_PRINTF(
      "siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet "
      "yref=%e, yDot=%e .\n",
      y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!std::isnan(y));
  if (y > 0)
    DEBUG_PRINT(
        "siconos::mechanisms::MBTB_MoreauJeanOSI::"
        "removeInteractionFromIndexSet DEACTIVATE.\n");
  return (y > 0);
}
#endif

#ifdef FIRSTWAY_ACTIVATION
bool siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(i == 1);
  double y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  // double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity

  DEBUG_PRINTF(
      "siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet "
      "yref=%e, yDot=%e.\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet  "
      "_activateYPosThreshold "
      "=%e, _activateYVelThreshold=%e\n",
      _activateYPosThreshold, _activateYVelThreshold);

  assert(!std::isnan(y));

  if (y <= _activateYPosThreshold)
    DEBUG_PRINT(
        "siconos::mechanisms::MBTB_MoreauJeanOSI::addInteractionInIndexSet "
        "ACTIVATE.\n");
  return (y <= _activateYPosThreshold);
}

bool siconos::mechanisms::MBTB_MoreauJeanOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)

{
  assert(i == 1);
  //  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0);  // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0);   // for i=1 y(i) is the velocity

  DEBUG_PRINTF(
      "siconos::mechanisms::MBTB_MoreauJeanOSI::removeInteractionFromIndexSet "
      "yref=%e, yDot=%e .\n",
      y, yDot);

  DEBUG_PRINTF(
      "siconos::mechanisms::MBTB_MoreauJeanOSI::removeInteractionFromIndexSet  "
      "_deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
      _deactivateYPosThreshold, _deactivateYVelThreshold);

  assert(!std::isnan(y));
  if (y > _deactivateYPosThreshold && yDot >= _deactivateYVelThreshold)
    DEBUG_PRINT(
        "siconos::mechanisms::MBTB_MoreauJeanOSI::"
        "removeInteractionFromIndexSet DEACTIVATE.\n");
  return (y > _deactivateYPosThreshold && yDot >= _deactivateYVelThreshold);
}

#endif
