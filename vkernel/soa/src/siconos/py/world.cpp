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

// The one place where the disks world is constructed, hence where its columns
// are registered: kept in its own translation unit, compiled in parallel with
// the binding parts (see data_t in py_head.hpp).
#include "py_head.hpp"

namespace siconos::vkernel::python::disks {
data_t::data_t() : _data(idata_t{}) {}
}  // namespace siconos::vkernel::python::disks
