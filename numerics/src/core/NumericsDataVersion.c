/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "NumericsDataVersion.h"

#include "assert.h"
#include "numerics_errors.h"

void NDV_inc(NumericsDataVersion* v) {
  assert(v->number < UINT64_MAX);
  v->number += 1;
}
