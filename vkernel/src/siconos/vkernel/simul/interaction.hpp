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

#include "siconos/vkernel/simul/simul_head.hpp"

namespace siconos::vkernel::simul {

/* compile time interaction */
template <match::item Nslaw, match::item... Relations>
struct interaction : item {
  using nslaw = Nslaw;
  using relations = gather<Relations...>;

  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<nslaw::size>;

  struct attributes {
    some::polymorphic_attribute<some::item_ref<Relations>...> relation;
    some::item_ref<Nslaw> nslaw;
    some::matrix<some::scalar, nslaw_size, dof> h_matrix1;
    some::matrix<some::scalar, nslaw_size, dof> h_matrix2;
    some::vector<some::scalar, nslaw_size> lambda;
    some::vector<some::scalar, nslaw_size> y;
    some::vector<some::scalar, nslaw_size> ydot;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return attr<"nslaw">(*self()); }

    // read only: a relation is set by topology::bind, which keeps the
    // relation sets
    const auto& relation() { return attr<"relation">(*self()); }
    decltype(auto) h_matrix1() { return attr<"h_matrix1">(*self()); }
    decltype(auto) h_matrix2() { return attr<"h_matrix2">(*self()); }

    decltype(auto) lambda() { return attr<"lambda">(*self()); }
    decltype(auto) y() { return attr<"y">(*self()); }
    decltype(auto) ydot() { return attr<"ydot">(*self()); }
  };
};

/* runtime / compile_time interaction */
template <match::item Nslaw, match::item... Relations>
struct rt_ct_interaction : item {
  using nslaw = Nslaw;
  using relations = gather<Relations...>;

  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<nslaw::size>;

  using attributes = gather<
      attribute<"relation", some::polymorphic_attribute<some::item_ref<Relations>...>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"h_matrix1", some::matrix<some::scalar, nslaw_size, dof>>,
      attribute<"h_matrix2", some::matrix<some::scalar, nslaw_size, some::indice_value<4>>>,
      attribute<"lambda", some::vector<some::scalar, nslaw_size>>,
      attribute<"y", some::vector<some::scalar, nslaw_size>>,
      attribute<"ydot", some::vector<some::scalar, nslaw_size>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return attr<"nslaw">(*self()); }

    // read only: a relation is set by topology::bind, which keeps the
    // relation sets
    const auto& relation() { return attr<"relation">(*self()); }
    decltype(auto) h_matrix1() { return attr<"h_matrix1">(*self()); }
    decltype(auto) h_matrix2() { return attr<"h_matrix2">(*self()); }

    decltype(auto) lambda() { return attr<"lambda">(*self()); }
    decltype(auto) y() { return attr<"y">(*self()); }
    decltype(auto) ydot() { return attr<"ydot">(*self()); }
  };
};

/* runtime / runtime interaction: draft */
template <match::item Nslaw, match::item... Relations>
struct rt_rt_interaction : item {
  using nslaw = Nslaw;
  using relations = gather<Relations...>;

  using nslaw_size = some::indice_value<nslaw::size>;

  using attributes =
      gather<attribute<"relation", some::polymorphic_attribute<some::item_ref<Relations>...>>,
             attribute<"nslaw", some::item_ref<nslaw>>,
             attribute<"h_matrix1", some::unbounded_col_matrix<some::scalar, nslaw_size>>,
             attribute<"h_matrix2", some::unbounded_col_matrix<some::scalar, nslaw_size>>,
             attribute<"lambda", some::vector<some::scalar, nslaw_size>>,
             attribute<"y", some::vector<some::scalar, nslaw_size>>,
             attribute<"ydot", some::vector<some::scalar, nslaw_size>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return attr<"nslaw">(*self()); }

    // read only: a relation is set by topology::bind, which keeps the
    // relation sets
    const auto& relation() { return attr<"relation">(*self()); }
    decltype(auto) h_matrix1() { return attr<"h_matrix1">(*self()); }
    decltype(auto) h_matrix2() { return attr<"h_matrix2">(*self()); }

    decltype(auto) lambda() { return attr<"lambda">(*self()); }
    decltype(auto) y() { return attr<"y">(*self()); }
    decltype(auto) ydot() { return attr<"ydot">(*self()); }
  };
};

}  // namespace siconos::vkernel::simul
