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

#include "siconos/storage/add.hpp"
#include "siconos/storage/data_holder.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/permute.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/registry.hpp"
#include "siconos/storage/remove.hpp"
#include "siconos/storage/sparse_set.hpp"

namespace siconos::vkernel::storage {

// template <match::item T>
// static constexpr void for_each_attribute(T)
// {
//   return mp::compose(mp::for_each, attributes)(T{});
// };

using pattern::attr_t;
using pattern::wrap;

}  // namespace siconos::vkernel::storage
