/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <siconos/vkernel/storage/pattern/pattern.hpp>

using namespace siconos::vkernel::storage::pattern;
using namespace siconos::vkernel::storage::some;

static_assert(
    std::is_same_v<decltype(cons(int{}, gather<char, float>{})), gather<int, char, float>>);

static_assert(std::is_same_v<decltype(append(gather<int>{}, gather<float>{})),
                             gather<int, gather<float>>>);

static_assert(
    std::is_same_v<decltype(concat(gather<int>{}, gather<float>{})), gather<int, float>>);

// static_assert(
//     attribute_name(siconos::vkernel::storage::pattern::attribute<
//                    "hello", siconos::vkernel::storage::some::scalar>{}) ==
//     "hello");

int main() {};
