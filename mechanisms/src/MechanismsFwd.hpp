// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MechanismsFwd_hpp
#define MechanismsFwd_hpp
#include <MechanicsFwd.hpp>
#include <SiconosPointers.hpp>
#define MECHANISMS_CLASSES()         \
  REGISTER(MBTB_FC3DContactRelation) \
  REGISTER(MBTB_ContactRelation)

#include <SiconosVisitables.hpp>

#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES() \
  KERNEL_CLASSES()           \
  MECHANICS_CLASSES()

#undef REGISTER
#undef REGISTER_BASE
#define REGISTER(X) DEFINE_SPTR(X);
#define REGISTER_BASE(X, Y) DEFINE_SPTR(X);
MECHANISMS_CLASSES();
#undef REGISTER

#endif
