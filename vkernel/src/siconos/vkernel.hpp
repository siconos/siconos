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

#include <charconv>
#include <functional>
#include <tuple>
#include <typeinfo>

#include "siconos/vkernel/algebra/numerics.hpp"
#include "siconos/vkernel/collision/diskdisk_r.hpp"
#include "siconos/vkernel/collision/diskline_r.hpp"
#include "siconos/vkernel/collision/neighborhood.hpp"
#include "siconos/vkernel/collision/point.hpp"
#include "siconos/vkernel/collision/shape/disk.hpp"
#include "siconos/vkernel/collision/shape/line.hpp"
#include "siconos/vkernel/collision/shape/polytope.hpp"
#include "siconos/vkernel/config/config.hpp"
#include "siconos/vkernel/config/environment.hpp"
#include "siconos/vkernel/io/io.hpp"
#include "siconos/vkernel/model/lagrangian_ds.hpp"
#include "siconos/vkernel/model/lagrangian_r.hpp"
#include "siconos/vkernel/model/nslaws.hpp"
#include "siconos/vkernel/simul/interaction.hpp"
#include "siconos/vkernel/simul/interaction_manager.hpp"
#include "siconos/vkernel/simul/one_step_integrator.hpp"
#include "siconos/vkernel/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/vkernel/simul/time_discretization.hpp"
#include "siconos/vkernel/simul/time_stepping.hpp"
#include "siconos/vkernel/simul/topology.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/pattern/pattern.hpp"
#include "siconos/vkernel/storage/storage.hpp"
#if defined(__clang__)
#endif

namespace siconos::vkernel {
using siconos::vkernel::storage::access;
using siconos::vkernel::storage::default_interface;
using siconos::vkernel::storage::prop;

using siconos::vkernel::storage::pattern::collect;
using siconos::vkernel::storage::pattern::method;

using namespace boost::hana::literals;

}  // namespace siconos::vkernel
