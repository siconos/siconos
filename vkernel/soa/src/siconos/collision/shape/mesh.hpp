#pragma once

#include "siconos/collision/collision.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/chained_segment.hpp"

namespace siconos::collision::shape {

struct mesh : item {
  using items = gather<chained_segment>;

  struct attributes {
    some::item_ref<chained_segment> segments;
    some::unbounded_vector<some::indice> global_indices;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    void __init__()
    {
      // Create the associated chained_segment
      auto& data = self()->data();
      auto seg = storage::add<chained_segment>(data);
      attr<"segments">(*self()) = seg.index();
    }

    void __del__()
    {
      throw std::runtime_error("Stable pointers need to be implemented");
      // Create the associated chained_segment
      storage::remove(segments());
    }

    decltype(auto) segments()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"segments">(*self()));
    }

    decltype(auto) global_indices()
    {
      return storage::attr<"global_indices">(*self());
    }

    decltype(auto) set_nodes(auto& nodes)
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto segments = self()->segments();

      // segments.nodes().clear();

      for (indice i = 0; i < algebra::rows(nodes); ++i) {
        segments.insert(nodes(i, 0), nodes(i, 1), nodes(i, 2));
      }
    }

    decltype(auto) indices() { return storage::attr<"indices">(*self()); }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using unbounded_matrix =
          typename env_t::template unbounded_matrix<scalar>;

      return collect(method("set_nodes",
                            &interface<Handle>::set_nodes<unbounded_matrix>));
    }
  };
};

}  // namespace siconos::collision::shape
