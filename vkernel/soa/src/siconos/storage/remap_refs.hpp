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

#include <cstddef>

#include "siconos/storage/refs.hpp"

namespace siconos::vkernel::storage {

/**
 * The rows of the items of Item moved as r says (see ref_remap): rewrite the
 * refs to them held in the store, in every column and every kept step.
 * storage::remove, storage::remove_many and storage::permute call it.
 *
 * The columns registered as referring to Item, dynamic properties included.
 */
template <typename Item>
void remap_refs(auto& data, const ref_remap& r) {
  data.store()._registry.template remap_refs<Item>(r);
}

}  // namespace siconos::vkernel::storage
