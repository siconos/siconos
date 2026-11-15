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

#include "siconos/storage/storage.hpp"

namespace siconos::vkernel::config {
template <typename... Cfs>
using map = storage::mp::map<Cfs...>;

template <typename K, typename V>
using assoc = storage::mp::key_value<K, V>;

using storage::pattern::param;
using storage::pattern::param_type;
using storage::pattern::param_val;

template <storage::pattern::string_literal S, auto N>
using iparam = assoc<param<S>, param_val<N>>;

template <storage::pattern::string_literal S, typename T>
using tparam = assoc<param<S>, param_type<T>>;
}  // namespace siconos::vkernel::config
