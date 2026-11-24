// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NumericsDataVersion.h"

#include "assert.h"
#include "numerics_errors.h"

void NDV_inc(NumericsDataVersion* v) {
  assert(v->number < UINT64_MAX);
  v->number += 1;
}
