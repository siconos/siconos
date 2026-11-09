#pragma once

#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {

template <typename SpaceFilter>
struct interaction_manager : item {
  using space_filter = SpaceFilter;
  using nslaw = typename space_filter::nslaw;

  using ncgroups = some::indice_parameter<"ncgroups">;

  using attributes =
      gather<attribute<"space_filter", some::item_ref<space_filter>>,
             attribute<"nslaws", some::matrix<some::item_ref<nslaw>, ncgroups, ncgroups>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) insert_nonsmooth_law(auto nslaw, auto gid1, auto gid2) {
      attr<"nslaws"> (*self())(gid1, gid2) = nslaw;
      return *self();
    }

    decltype(auto) get_nonsmooth_law(auto gid1, auto gid2) {
      auto& data = self()->data();
      return storage::make_handle(data, attr<"nslaws">(*self())(gid1, gid2));
    }

    void update_interactions(auto step) {
      auto& data = self()->data();
      storage::make_handle(data, attr<"space_filter">(*self())).update_index_set0(step);
    }

    auto methods() {
      using data_t = std::decay_t<decltype(self()->data())>;
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      // a handle, not an index: python passes handles, and pybind would
      // reinterpret a handle (data reference first) as an index
      using ref_nslaw_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(), storage::index<nslaw, indice>{}))>;

      return collect(
          method("insert_nonsmooth_law",
                 &interface<Handle>::insert_nonsmooth_law<ref_nslaw_t, indice, indice>),
          method("update_interactions", &interface<Handle>::update_interactions<indice>));
    }
  };
};
}  // namespace siconos::simul
