#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/storage.hpp"

#define USE_DOUBLE
#ifdef WITH_GPU
#include "cuNSearch.h"
namespace CompactNSearch = cuNSearch;
#else
#include "CompactNSearch/CompactNSearch.h"
#endif

namespace siconos::collision {

/* the neighborhood engine */
template <typename... Points>
struct neighborhood
    : storage::data_holder<CompactNSearch::NeighborhoodSearch> {
  using items = gather<Points...>;
  using points_t = gather<Points...>;

  using attributes = storage::pattern::cons_x<
      attribute<"point_set_id",
                some::array<some::indice,
                            some::indice_value<ground::size(points_t{})>>>,
      typename storage::data_holder<
          CompactNSearch::NeighborhoodSearch>::attributes>;

  template <typename Handle>
  struct interface : storage::data_holder<
                         CompactNSearch::NeighborhoodSearch>::
                         template interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) point_set_id()
    {
      return storage::attr<"point_set_id">(*self());
    };

    void create(auto radius)
    {
      self()->instance().reset(
          new CompactNSearch::NeighborhoodSearch(radius));
    }

    void add_point_sets(auto step)
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto& data = self()->data();

      indice i = 0;

      auto& psid = storage::attr<"point_set_id">(*self());
      auto& instance = self()->instance();
      ground::for_each(points_t{}, [&data, &step, &i, &psid,
                                    &instance]<typename Point>(Point) {
        auto& coords = storage::attr_values<Point, "coord">(data, step);

        // std::vector assumed for coords
        psid[i++] =
            instance->add_point_set(coords.front().data(), coords.size());
      });
    }

    void set_active(auto ps1_id, auto ps2_id, auto value)
    {
      self()->instance()->set_active((unsigned int)ps1_id,
                                     (unsigned int)ps2_id, value);
    }

    bool is_active(auto i, auto j)
    {
      return self()->instance()->is_active(i, j);
    }

    void update(auto step)
    {
      auto& data = self()->data();
      ground::for_each(points_t{}, [&data, &step]<typename Point>(Point) {
        for (auto point : storage::handles<Point>(data, step)) {
          point.update();
        }
      });
    }

    void search() { self()->instance()->find_neighbors(); };

    void sort()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto& data = self()->data();
      auto& instance = self()->instance();

      instance->z_sort();

      indice i = 0;
      ground::for_each(
          points_t{}, [&instance, &data, &i]<typename Point>(Point p) {
            auto ps = instance->point_set(i++);
            // apply function only if some points exist
            if (ps.n_points() > 0) {
              storage::apply_fun(data, p, [&ps]<typename Array>(Array& a) {
                ps.sort_field(a.data());
              });
            }
          });
    }
    auto methods()
    {
      using env_t = decltype(self()->env());

      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(
          method("point_set_id", &interface<Handle>::point_set_id),
          method("create", &interface<Handle>::create<scalar>),
          method("add_point_sets",
                 &interface<Handle>::add_point_sets<indice>),
          method("update", &interface<Handle>::update<indice>),
          method("set_active",
                 &interface<Handle>::set_active<indice, indice, bool>),
          method("is_active", &interface<Handle>::is_active<indice, indice>),
          method("search", &interface<Handle>::search),
          method("sort", &interface<Handle>::sort));
    }
  };
};

}  // namespace siconos::collision
