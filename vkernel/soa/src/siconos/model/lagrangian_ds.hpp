#pragma once

#include "siconos/model/model_head.hpp"

namespace siconos::model {

struct lagrangian_ds
    : item<description<"A lagrangian dynamical system [...]">> {
  using dof = some::indice_parameter<"dof">;

  using attributes =
      gather<attribute<"q", some::vector<some::scalar, dof>>,
             attribute<"velocity", some::vector<some::scalar, dof>>,
             attribute<"mass_matrix", some::matrix<some::scalar, dof, dof>>,
             attribute<"fext", some::vector<some::scalar, dof>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mass_matrix() { return attr<"mass_matrix">(*self()); }

    decltype(auto) velocity() { return attr<"velocity">(*self()); }

    decltype(auto) q() { return attr<"q">(*self()); }

    decltype(auto) fext() { return attr<"fext">(*self()); }
  };
};

  struct rt_lagrangian_ds : item<> {
  using without_attached_storages_bindings = void;

  using attributes =
      gather<attribute<"q", some::unbounded_vector<some::scalar>>,
             attribute<"velocity", some::unbounded_vector<some::scalar>>,
             attribute<"mass_matrix", some::unbounded_matrix<some::scalar>>,
             attribute<"fext", some::unbounded_vector<some::scalar>>>;
};
}  // namespace siconos::model
