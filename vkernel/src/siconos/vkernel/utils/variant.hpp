/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <siconos/vkernel/storage/storage.hpp>
#include <tuple>
#include <variant>

namespace siconos::vkernel::variant {

template <typename D, typename V, typename F>
decltype(auto) visit(D& data, V&& var, F&& fun) {
  return std::visit(
      [&data, &fun](auto& rvar) {
        auto h = storage::make_handle(data, rvar);
        return static_cast<F&&>(fun)(h);
      },
      static_cast<V&&>(var));
}
}  // namespace siconos::vkernel::variant
