#pragma once

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <typeinfo>

#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/diskfsegment_r.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/collision/neighborhood.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/chained_segment.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/mesh.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/collision/space_filter.hpp"
#include "siconos/collision/translated.hpp"
#include "siconos/io/io.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/range.hpp"

namespace py = pybind11;

namespace siconos::config::disks {

struct disk : model::lagrangian_ds {};
struct fem : model::elastic_lagrangian_ds {};

struct nslaw : model::newton_impact_friction {};
using diskdisk_r = collision::diskdisk_r;
using diskfdisk_r = collision::diskfdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
using diskmesh_r = collision::diskmesh_r;
using segment_shape = collision::shape::segment;
using disk_shape = collision::shape::disk;
using mesh_shape = collision::shape::mesh;
using chained_segment_shape = collision::shape::chained_segment;
using translated_disk_shape = collision::translated<disk_shape>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
using solver_options = simul::solver_options;
using trace_params = simul::trace_params;
struct ct_interaction
    : simul::interaction<nslaw, diskdisk_r, diskfdisk_r, diskfsegment_r> {};
struct rt_interaction : simul::rt_ct_interaction<nslaw, diskmesh_r> {};
struct topo : simul::topology<disk, ct_interaction, fem, rt_interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<disk, collision::empty_shape>;
using pointf = collision::point<fem, mesh_shape>;
using pointl = collision::point<storage::pattern::empty_item, segment_shape>;
using pointtds =
    collision::point<storage::pattern::empty_item, translated_disk_shape>;
struct neighborhood
    : collision::neighborhood<pointd, pointf, pointl, pointtds> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

struct io : siconos::io::io<osi, segment_shape, translated_disk_shape> {};
template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};
}  // namespace siconos::config::disks

using namespace siconos;
namespace some = siconos::storage::some;

namespace pattern = siconos::storage::pattern;

namespace siconos::python::disks {

namespace config = siconos::config::disks;

struct maker
    : storage::make<
          config::env, config::simulation, config::interaction_manager,
          config::io, config::segment_shape, config::disk_shape,
          config::mesh_shape,
          storage::with_properties<
              storage::wrapped<config::disk, some::unbounded_collection>,
              storage::wrapped<config::diskdisk_r,
                               some::unbounded_collection>,
              storage::wrapped<config::diskfsegment_r,
                               some::unbounded_collection>,
              storage::wrapped<config::diskmesh_r,
                               some::unbounded_collection>,
              storage::wrapped<config::diskfdisk_r,
                               some::unbounded_collection>,
              storage::wrapped<config::pointl, some::unbounded_collection>,
              storage::wrapped<config::pointd, some::unbounded_collection>,
              storage::wrapped<config::pointf, some::unbounded_collection>,
              storage::wrapped<config::pointtds, some::unbounded_collection>,
              storage::wrapped<config::ct_interaction,
                               some::unbounded_collection>,
              storage::wrapped<config::rt_interaction,
                               some::unbounded_collection>,
              storage::wrapped<config::segment_shape,
                               some::unbounded_collection>,
              storage::wrapped<config::mesh_shape,
                               some::unbounded_collection>,
              storage::wrapped<config::disk_shape,
                               some::unbounded_collection>,
              storage::wrapped<config::translated_disk_shape,
                               some::unbounded_collection>,
              storage::attached<config::disk,
                                storage::pattern::symbol<"shape">,
                                storage::some::item_ref<config::disk_shape>>,
              storage::attached<config::fem,
                                storage::pattern::symbol<"shape">,
                                storage::some::item_ref<config::mesh_shape>>,
              storage::sparse<
                  storage::pattern::attr_t<config::fem, "mass_matrix">>,
              storage::sparse<
                  storage::pattern::attr_t<config::fem, "k_matrix">>,
              storage::unbounded<storage::pattern::attr_t<config::fem, "q">>,
              storage::unbounded<
                  storage::pattern::attr_t<config::fem, "velocity">>,
              storage::unbounded<
                  storage::pattern::attr_t<config::fem, "fext">>,
              storage::time_invariant<
                  storage::pattern::attr_t<config::disk, "fext">>,
              storage::diagonal<
                  storage::pattern::attr_t<config::disk, "mass_matrix">>,
              storage::assembled_diagonal<storage::pattern::attr_t<
                  typename config::osi::assembled_osi_t,
                  "mass_matrix_assembled">>,
              storage::bind<config::disk, "disk">,
              storage::bind<config::fem, "fem">,
              storage::bind<config::nslaw, "nslaw">,
              storage::bind<config::diskdisk_r, "diskdisk_r">,
              storage::bind<config::diskfdisk_r, "diskfdisk_r">,
              storage::bind<config::diskfsegment_r, "diskfsegment_r">,
              storage::bind<config::diskmesh_r, "diskmesh_r">,
              storage::bind<config::neighborhood, "neighborhood">,
              storage::bind<config::space_filter, "space_filter">,
              storage::bind<config::segment_shape, "segment_shape">,
              storage::bind<config::chained_segment_shape, "chained_segment">,
              storage::bind<config::disk_shape, "disk_shape">,
              storage::bind<config::mesh_shape, "mesh_shape">,
              storage::bind<config::translated_disk_shape,
                            "translated_disk_shape">,
              storage::bind<config::interaction_manager,
                            "interaction_manager">,
              storage::bind<config::ct_interaction, "ct_interaction">,
              storage::bind<config::rt_interaction, "rt_interaction">,
              storage::bind<config::solver_options, "solver_options">,
              storage::bind<config::trace_params, "trace_params">,
              storage::bind<config::osi, "osi">,
              storage::bind<config::td, "time_discretization">,
              storage::bind<config::topo, "topology">,
              storage::bind<config::simulation, "simulation">,
              storage::bind<config::osnspb, "osnspb">,
              storage::bind<typename config::osi::assembled_osi_t,
                            "assembled_osi">,
              storage::bind<config::fc2d, "fc2d">,
              storage::bind<config::io, "io">>> {};

static auto imake_storage() { return maker(); }

static auto idata = imake_storage();

struct idata_t : std::decay_t<decltype(idata)> {};

// just hide idata_t to pybind11
struct data_t {
  data_t() : _data(idata_t{}) {};

  idata_t& operator()() { return _data; };

  idata_t _data;
};

}  // namespace siconos::python::disks

namespace ground = siconos::storage::mp;
namespace match = siconos::storage::pattern::match;
template <typename H, typename T>
static decltype(auto) out_formatter(H h, T&& out_value)
{
  using out_t = std::decay_t<T>;

  if constexpr (!match::diagonal_matrix<out_t>) {
    if constexpr (match::matrix<out_t>) {
      return algebra::matrix_ref<out_t>(static_cast<T&&>(out_value));
    }
    else if constexpr (match::index<out_t>) {
      auto ret_val = siconos::storage::handle<
          siconos::storage::handle_base, typename out_t::type,
          typename out_t::value_t, siconos::python::disks::idata_t>(
          h.data(), out_value);
      return ret_val;
    }
    else {
      return static_cast<T&&>(out_value);
    }
  }
  else {
    return out_value.diagonal();
  }
}

template <typename H, typename T>
static decltype(auto) in_formatter(H&& h, T&& in_value)
{
  using in_t = std::decay_t<T>;

  if constexpr (!match::diagonal_matrix<in_t>) {
    return static_cast<T&&>(in_value);
  }
  else {
    return in_value.diagonal();
  }
}

// Hide std::vector from pybind11 to avoid a list binding
template <typename T>
struct handles_wrap {
  T handles;

  handles_wrap(T&& hs) : handles(std::move(hs)) {}
};

using namespace boost::hana::literals;
