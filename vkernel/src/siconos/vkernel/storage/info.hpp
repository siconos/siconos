// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/pattern/pattern.hpp"

namespace siconos::vkernel::storage {

using namespace pattern;
// the storage info key
struct info {};

template <typename D>
auto get_info(D&& data) {
  return mp::get<info>(data.store());
}

template <typename D>
using get_info_t = decltype(get_info(std::decay_t<D>{}));

}  // namespace siconos::vkernel::storage
