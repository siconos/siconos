/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <siconos/storage/storage.hpp>
#include <tuple>
#include <variant>

namespace siconos::variant {

template <typename D, typename V, typename F>
decltype(auto) visit(D& data, V&& var, F&& fun) {
  return std::visit(
      [&data, &fun](auto& rvar) {
        auto h = storage::make_handle(data, rvar);
        return static_cast<F&&>(fun)(h);
      },
      static_cast<V&&>(var));
}
}  // namespace siconos::variant
