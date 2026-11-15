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

#include <stdexcept>
#include <vector>

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/remap_refs.hpp"

namespace siconos::storage {

/** Reorder the items of type Item: the item of index k becomes the item of
    index perm[k] before the call, in every attribute, every attached storage
    and every dynamic property (re-keyed), in every kept step.

    perm is a permutation of [0, number of items). The refs to the items of
    Item held in the store follow them (see remap_refs); the handles and the
    refs held outside the store do not. */
template <typename Item>
void permute(auto& data, const std::vector<std::size_t>& perm) {
  std::vector<std::size_t> inverse(std::size(perm));
  for (std::size_t k = 0; k < std::size(perm); ++k) {
    if (perm[k] >= std::size(perm)) [[unlikely]] {
      throw std::invalid_argument("storage::permute: index out of the permutation");
    }
    inverse[perm[k]] = k;
  }

  // every column of Item at every kept step, the refs to the moved items
  // with them (remap_refs)
  data.store()._registry.template permute_item<Item>(perm, inverse);
}

}  // namespace siconos::storage
