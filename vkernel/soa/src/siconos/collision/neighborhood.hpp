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

#include <numeric>
#include <vector>

#include "siconos/collision/collision_head.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/permute.hpp"
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
struct neighborhood : storage::data_holder<CompactNSearch::NeighborhoodSearch> {
  using items = gather<Points...>;
  using points_t = gather<Points...>;

  using attributes = storage::pattern::cons_x<
      attribute<"point_set_id",
                some::array<some::indice, some::indice_value<mp::size(points_t{})>>>,
      typename storage::data_holder<CompactNSearch::NeighborhoodSearch>::attributes>;

  template <typename Handle>
  struct interface : storage::data_holder<
                         CompactNSearch::NeighborhoodSearch>::template interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) point_set_id() { return storage::attr<"point_set_id">(*self()); };

    void create(auto radius) {
      self()->instance().reset(new CompactNSearch::NeighborhoodSearch(radius));
    }

    void add_point_sets(auto step) {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto& data = self()->data();

      indice i = 0;

      auto& psid = storage::attr<"point_set_id">(*self());
      auto& instance = self()->instance();
      mp::for_each(points_t{}, [&data, &step, &i, &psid, &instance]<typename Point>(Point) {
        auto& coords = storage::attr_values<Point, "coord">(data, step);

        // std::vector assumed for coords
        if (coords.empty()) {
          // dummy
          float dummy[3] = {0, 0, 0};
          psid[i++] = instance->add_point_set(dummy, 0);
        } else {
          // fixed-shape point kinds (empty_item specializations: segments,
          // fixed disks) are kinematic: CompactNSearch must not recompute
          // their positions every search; only moving bodies (pointd) and,
          // when present, the FEM mesh probes (pointf) are dynamic.
          psid[i++] = instance->add_point_set(coords.front().data(), coords.size(),
                                              Point::is_dynamic_v);
        }
      });
    }

    void set_active(auto ps1_id, auto ps2_id, auto value) {
      self()->instance()->set_active((unsigned int)ps1_id, (unsigned int)ps2_id, value);
    }

    bool is_active(auto i, auto j) { return self()->instance()->is_active(i, j); }

    void update(auto step) {
      auto& data = self()->data();
      mp::for_each(points_t{}, [&data, &step]<typename Point>(Point) {
        for (auto point : storage::handles<Point>(data, step)) {
          point.update(step);
        }
      });
    }

    void search() { self()->instance()->find_neighbors(); };

    void sort() {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto& data = self()->data();
      auto& instance = self()->instance();

      instance->z_sort();

      indice i = 0;
      mp::for_each(points_t{}, [&instance, &data, &i]<typename Point>(Point) {
        auto&& ps = instance->point_set(i++);
        if (ps.n_points() > 0) {
          // the z-order of the points as a permutation (sort_field puts the
          // item perm[k] at k), applied by storage::permute to every storage
          // of the points, dynamic properties included. The columns keep
          // their buffers: the point set still refers to the coordinates
          std::vector<std::size_t> perm(ps.n_points());
          std::iota(perm.begin(), perm.end(), std::size_t{0});
          ps.sort_field(perm.data());
          storage::permute<Point>(data, perm);
        }
      });
    }
    auto methods() {
      using env_t = decltype(self()->env());

      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(
          method("create", &interface<Handle>::create<scalar>),
          method("add_point_sets", &interface<Handle>::add_point_sets<indice>),
          method("update", &interface<Handle>::update<indice>),
          method("set_active", &interface<Handle>::set_active<indice, indice, bool>),
          method("search", &interface<Handle>::search),
          method("sort", &interface<Handle>::sort));
    }
  };
};

}  // namespace siconos::collision
