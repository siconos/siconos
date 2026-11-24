// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MBTB_TimeSteppingCombinedProj.hpp"

#include "MBTB_PYTHON_API.hpp"
#include "MBTB_internalTool.hpp"
// #define TSPROJ_DEBUG

void siconos::mechanisms::MBTB_TimeSteppingCombinedProj::updateWorldFromDS() {
#ifdef TSPROJ_DEBUG
  printf("MBTB_TimeSteppingCombinedProj::updateWordFromDS \n");
#endif
  MBTB_updateDSFromSiconos();
  mbtb::internal::MBTB_updateContactFromDS();
}
