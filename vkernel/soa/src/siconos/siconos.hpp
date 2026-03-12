#pragma once

#include <charconv>
#include <functional>
#include <tuple>
#include <typeinfo>

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskline_r.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/line.hpp"
#include "siconos/collision/neighborhood.hpp"
#include "siconos/collision/space_filter.hpp"
#include "siconos/config/config.hpp"
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
#include "siconos/utils/environment.hpp"
#if defined(__clang__)
#include <boost/hana/experimental/type_name.hpp>
#endif

namespace siconos {
using siconos::storage::access;
using siconos::storage::default_interface;
using siconos::storage::prop;

using siconos::storage::pattern::collect;
using siconos::storage::pattern::method;

using namespace boost::hana::literals;

}  // namespace siconos
