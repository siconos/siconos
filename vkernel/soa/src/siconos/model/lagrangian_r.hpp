#pragma once

#include <functional>  // mem_fn

#include "siconos/algebra/algebra.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/model_head.hpp"
namespace siconos::storage::pattern::match {
template <typename T>
concept linear_relation = match::handle<T, model::linear>;

template <typename T>
concept relation1 = match::handle<T, model::relation1>;

template <typename T>
concept relation2 = match::handle<T, model::relation2>;

}  // namespace siconos::storage::pattern::match

namespace siconos::model {

template <auto NSLSize>
struct lagrangian_r : item,
                      linear,
                      relation1,
                      relation2,
                      any_lagrangian_relation {
  using nslaw_size = some::param_val<NSLSize>;
  using dof = some::indice_parameter<"dof">;

  using attributes = gather<
      attribute<"b", some::scalar>,
      attribute<"h_matrix", some::matrix<some::scalar, nslaw_size, dof>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) h_matrix() { return attr<"h_matrix">(*self()); }

    decltype(auto) b() { return attr<"b">(*self()); }

    decltype(auto) compute_jachq(auto step, auto& ds1, auto& ds2,
                                 auto& h_matrix1, auto& h_matrix2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }
    decltype(auto) compute_jachq(auto step, auto& ds, auto& h_matrix1)
    {
      h_matrix1 << -h_matrix();
    }
  };
};

}  // namespace siconos::model
