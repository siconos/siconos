/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <cmath>

// The scalar functions of the scene, independent of the algebra backend: a
// collision relation needs hypot and nothing else, and must not pull the
// assembled matrices of numerics_backend.hpp.
namespace siconos::vkernel::algebra {
using std::hypot;
}
