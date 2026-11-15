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

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"

namespace siconos::storage {

using namespace pattern;
// the storage info key
struct info {};

template <typename D>
auto get_info(D&& data) {
  return mp::get<info>(data.store());
}

template <typename D>
using get_info_t = decltype(get_info(std::decay_t<D>{}));

}  // namespace siconos::storage
