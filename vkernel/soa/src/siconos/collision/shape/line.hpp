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

#include <math.h>

#include <cassert>

#include "siconos/collision/collision_head.hpp"

namespace siconos::collision::shape {

struct line : item {
  // a*x + b*y + c = 0
  using attributes =
      gather<attribute<"a", some::scalar>, attribute<"b", some::scalar>,
             attribute<"c", some::scalar>,
             attribute<"p0", some::vector<some::scalar, some::indice_value<3>>>,
             attribute<"direction", some::vector<some::scalar, some::indice_value<3>>>,
             attribute<"maxpoints", some::scalar>, attribute<"invsqrta2pb2", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;
    decltype(auto) a() { return storage::attr<"a">(*self()); }
    decltype(auto) b() { return storage::attr<"b">(*self()); }
    decltype(auto) c() { return storage::attr<"c">(*self()); }

    decltype(auto) p0() { return storage::attr<"p0">(*self()); }
    decltype(auto) direction() { return storage::attr<"direction">(*self()); }
    decltype(auto) maxpoints() { return storage::attr<"maxpoints">(*self()); };
    decltype(auto) invsqrta2pb2() { return storage::attr<"invsqrta2pb2">(*self()); }

    void initialize() {
      using vect_t = std::decay_t<decltype(p0())>;

      if (b()) {
        p0() = vect_t{{0., -c() / b(), 0.}};
      } else {
        assert(a());
        p0() = vect_t{{-c() / a(), 0., 0.}};
      }

      invsqrta2pb2() = 1. / (sqrt(a() * a() + b() * b()));

      direction() = vect_t{{-b(), a(), 0.}} * invsqrta2pb2();
    }
    decltype(auto) points_coords() {
      auto& p0 = self()->p0();
      auto& maxpoints = self()->maxpoints();
      auto& dir = self()->direction();
      return view::iota(0, self()->maxpoints()) |
             view::transform([&](auto i) { return p0 + ((i - maxpoints / 2) / 10.) * dir; });
    }
    decltype(auto) distance(auto& q) {
      auto& x = q(0);
      auto& y = q(1);

      return fabs(attr<"a">(*self()) * x + attr<"b">(*self()) * y + attr<"c">(*self())) *
             attr<"invsqrta2pb2">(*self());
    }

    auto methods() { return collect(method("initialize", &interface<Handle>::initialize)); }
  };
};
}  // namespace siconos::collision::shape
