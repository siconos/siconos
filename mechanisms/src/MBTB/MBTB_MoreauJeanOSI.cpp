// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionsGraph::size_type i) {
  assert(i == 1);
  double h = simulation()->timeStep();
  double y = (inter->y(i - 1))->getValxxue(0);  // for i=1 y(i-1) is the position
  double yDot = (*(inter->y(i)))(0);            // for i=1 y(i) is the velocity

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
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionsGraph::size_type i) {
  assert(i == 1);
  double h = simulation()->timeStep();
  double y = (*inter->y(i - 1))(0);   // for i=1 y(i-1) is the position
  double yDot = (*(inter->y(i)))(0);  // for i=1 y(i) is the velocity
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
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionsGraph::size_type i) {
  assert(i == 1);
  double y = (*inter->y(i - 1))(0);  // for i=1 y(i-1) is the position
  // double yDot = (*(inter->y(i)))(0); // for i=1 y(i) is the velocity

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
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionsGraph::size_type i)

{
  assert(i == 1);
  //  double h = simulation()->timeStep();
  double y = (*inter->y(i - 1))(0);   // for i=1 y(i-1) is the position
  double yDot = (*(inter->y(i)))(0);  // for i=1 y(i) is the velocity

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
