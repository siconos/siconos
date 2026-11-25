/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <cstddef>

#include "siconos/vkernel/storage/refs.hpp"

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
