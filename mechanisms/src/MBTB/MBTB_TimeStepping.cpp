/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "MBTB_TimeStepping.hpp"

#include "MBTB_PYTHON_API.hpp"    // For MBTB_updateDSFromSiconos
#include "MBTB_internalTool.hpp"  // For  _MBTB_updateContactFromDS();
// #define TS_DEBUG

void siconos::mechanisms::MBTB_TimeStepping::updateWorldFromDS() {
#ifdef TS_DEBUG
  printf("MBTB_TimeStepping::updateWordFromDS \n");
#endif
  MBTB_updateDSFromSiconos();
  mbtb::internal::MBTB_updateContactFromDS();
}
