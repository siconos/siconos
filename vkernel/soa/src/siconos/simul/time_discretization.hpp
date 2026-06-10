#pragma once

#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {

template <typename... Params>
struct time_discretization : item {

  struct attributes {
    some::scalar h;
    some::scalar t0;
    some::scalar tmax;
    some::indice step;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) t0() { return attr<"t0">(*self()); };
    decltype(auto) tmax() { return attr<"tmax">(*self()); };
    decltype(auto) h() { return attr<"h">(*self()); };
    decltype(auto) step() { return attr<"step">(*self()); };
  };
};
}  // namespace siconos::simul
