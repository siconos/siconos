/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "siconos/vkernel/storage/default_interface.hpp"
#include "siconos/vkernel/storage/some/some.hpp"

namespace siconos::vkernel::storage {

template <typename Struct>
struct data_holder : item {
  using attributes = gather<attribute<"instance", some::specific<pointer<Struct>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) instance() { return attr<"instance">(*self()); };
  };
};

}  // namespace siconos::vkernel::storage
