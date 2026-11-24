// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/storage/get.hpp"
#include "siconos/vkernel/storage/handle.hpp"
#include "siconos/vkernel/storage/info.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"

namespace siconos::vkernel::storage {

template <typename Handle>
struct default_interface {
  decltype(auto) self() {
    return static_cast<Handle*>(this);  // handle inherits from default_interface
  };

  auto env() {
    auto& data = self()->data();
    using info_t = get_info_t<decltype(data)>;
    return typename info_t::template env<typename Handle::type>{};
  }

  auto params() { return typename decltype(env())::params{}; }

  template <string_literal S>
  constexpr auto env_param() {
    return mp::get<pattern::param<S>>(self()->params()).value;
  }

  auto make_handle(auto& index) { return storage::make_handle(self()->data(), index); }
};

}  // namespace siconos::vkernel::storage
