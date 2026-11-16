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
