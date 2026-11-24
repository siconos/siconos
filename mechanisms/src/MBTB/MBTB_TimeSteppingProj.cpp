// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MBTB_TimeSteppingProj.hpp"

#include "MBTB_PYTHON_API.hpp"
#include "MBTB_internalTool.hpp"
// #define TSPROJ_DEBUG

void siconos::mechanisms::MBTB_TimeSteppingProj::updateWorldFromDS() {
#ifdef TSPROJ_DEBUG
  printf("MBTB_TimeSteppingProj::updateWordFromDS \n");
#endif
  MBTB_updateDSFromSiconos();
  mbtb::internal::MBTB_updateContactFromDS();
}
