// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/storage/pattern/pattern.hpp"
#include "siconos/vkernel/storage/some/some.hpp"
#include "siconos/vkernel/storage/storage.hpp"

namespace siconos::vkernel::collision {
using namespace storage;
using namespace storage::pattern;

template <typename...>
constexpr bool always_false = false;

}  // namespace siconos::vkernel::collision
