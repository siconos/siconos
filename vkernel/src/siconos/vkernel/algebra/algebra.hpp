// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

// The whole algebra layer: the contract (linear_algebra.hpp), the scalar
// functions (math.hpp), the per-system backend (eigen_backend.hpp) and the
// assembled backend (numerics_backend.hpp, with sparse_factor.hpp).  A file
// that needs only one of them includes that one: a collision relation takes
// math.hpp, the model takes eigen_backend.hpp, the simulation takes
// numerics_backend.hpp.

#include "siconos/vkernel/algebra/eigen_backend.hpp"
#include "siconos/vkernel/algebra/linear_algebra.hpp"
#include "siconos/vkernel/algebra/math.hpp"
#include "siconos/vkernel/algebra/numerics_backend.hpp"
#include "siconos/vkernel/algebra/sparse_factor.hpp"
