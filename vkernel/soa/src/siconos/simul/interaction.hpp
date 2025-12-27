#pragma once

#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {

template <match::item Nslaw, match::item... Relations>
struct interaction : item {
  using nslaw = Nslaw;
  using relations = gather<Relations...>;

  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<nslaw::size>;

  using attributes = gather<
      attribute<"relation",
                some::polymorphic_attribute<some::item_ref<Relations>...>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"h_matrix1", some::matrix<some::scalar, nslaw_size, dof>>,
      attribute<"h_matrix2", some::matrix<some::scalar, nslaw_size, dof>>,
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

/* runtime interaction */
template <match::item Nslaw, match::item... Relations>
struct rt_interaction : item {
  using nslaw = Nslaw;
  using relations = gather<Relations...>;

  using dof = some::indice_parameter<"dof">;
  using nslaw_size = some::indice_value<nslaw::size>;

  using attributes = gather<
      attribute<"relation",
                some::polymorphic_attribute<some::item_ref<Relations>...>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"h_matrix1", some::unbounded_col_matrix<some::scalar, nslaw_size>>,
      attribute<"h_matrix2", some::matrix<some::scalar, nslaw_size, dof>>,
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
