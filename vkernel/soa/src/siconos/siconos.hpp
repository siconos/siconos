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

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskline_r.hpp"
#include "siconos/collision/neighborhood.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/line.hpp"
#include "siconos/collision/shape/polytope.hpp"
#include "siconos/config/config.hpp"
#include "siconos/config/environment.hpp"
#include "siconos/io/io.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/interaction.hpp"
#include "siconos/simul/interaction_manager.hpp"
#include "siconos/simul/one_step_integrator.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/time_discretization.hpp"
#include "siconos/simul/time_stepping.hpp"
#include "siconos/simul/topology.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/storage.hpp"
#if defined(__clang__)
#endif

namespace siconos {
using siconos::storage::access;
using siconos::storage::default_interface;
using siconos::storage::prop;

using siconos::storage::pattern::collect;
using siconos::storage::pattern::method;

using namespace boost::hana::literals;

}  // namespace siconos
