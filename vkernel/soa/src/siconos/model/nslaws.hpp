#pragma once

#include "siconos/model/model_head.hpp"

namespace siconos::model {
struct equality_condition {};
struct relay {};

struct newton_impact_friction : item<> {

  static constexpr auto size = 2;

  using attributes =
      gather<attribute<"e", some::scalar>, attribute<"mu", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;
    decltype(auto) e() { return attr<"e">(*self()); }
    decltype(auto) mu() { return attr<"mu">(*self()); }
  };
};

struct newton_impact : item<> {
  static constexpr auto size = 1;
  using attributes = gather<attribute<"e", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;
    decltype(auto) e() { return attr<"e">(*self()); };
  };

};
}  // namespace siconos::model
