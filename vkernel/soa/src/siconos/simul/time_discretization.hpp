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

#include "siconos/simul/simul_head.hpp"

namespace siconos::vkernel::simul {

template <typename... Params>
struct time_discretization : item {
  struct attributes {
    some::scalar h;
    some::scalar t0;
    some::scalar tmax;
    some::indice step;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) t0() { return attr<"t0">(*self()); };
    decltype(auto) tmax() { return attr<"tmax">(*self()); };
    decltype(auto) h() { return attr<"h">(*self()); };
    decltype(auto) step() { return attr<"step">(*self()); };
  };
};
}  // namespace siconos::vkernel::simul
