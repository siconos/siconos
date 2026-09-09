#pragma once

#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {

/* compile time interaction */
template <match::item Nslaw, match::item... Relations>
struct interaction : item {
  using nslaw = Nslaw;
  using relations = gather<Relations...>;

  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<nslaw::size>;

  struct attributes {
    some::polymorphic_attribute<some::item_ref<Relations>...> relation;
    some::item_ref<nslaw> nslaw;
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

    decltype(auto) relation() { return attr<"relation">(*self()); };
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
      attribute<"relation",
                some::polymorphic_attribute<some::item_ref<Relations>...>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"h_matrix1", some::matrix<some::scalar, nslaw_size, dof>>,
      attribute<"h_matrix2", some::matrix<some::scalar, nslaw_size,
                                          some::indice_value<4>>>,
      attribute<"lambda", some::vector<some::scalar, nslaw_size>>,
      attribute<"y", some::vector<some::scalar, nslaw_size>>,
      attribute<"ydot", some::vector<some::scalar, nslaw_size>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return attr<"nslaw">(*self()); }

    decltype(auto) relation() { return attr<"relation">(*self()); };
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

  using attributes = gather<
      attribute<"relation",
                some::polymorphic_attribute<some::item_ref<Relations>...>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"h_matrix1",
                some::unbounded_col_matrix<some::scalar, nslaw_size>>,
      attribute<"h_matrix2",
                some::unbounded_col_matrix<some::scalar, nslaw_size>>,
      attribute<"lambda", some::vector<some::scalar, nslaw_size>>,
      attribute<"y", some::vector<some::scalar, nslaw_size>>,
      attribute<"ydot", some::vector<some::scalar, nslaw_size>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return attr<"nslaw">(*self()); }

    decltype(auto) relation() { return attr<"relation">(*self()); };
    decltype(auto) h_matrix1() { return attr<"h_matrix1">(*self()); }
    decltype(auto) h_matrix2() { return attr<"h_matrix2">(*self()); }

    decltype(auto) lambda() { return attr<"lambda">(*self()); }
    decltype(auto) y() { return attr<"y">(*self()); }
    decltype(auto) ydot() { return attr<"ydot">(*self()); }
  };
};

}  // namespace siconos::simul
