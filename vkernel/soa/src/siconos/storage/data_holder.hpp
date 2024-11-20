#pragma once

#include "siconos/storage/info.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/get.hpp"

namespace siconos::storage {

template <typename Handle>
struct default_interface {
  decltype(auto) self()
  {
    return static_cast<Handle*>(
        this);  // handle inherits from default_interface
  };

  auto env()
  {
    auto& data = self()->data();
    using info_t = get_info_t<decltype(data)>;
    return typename info_t::env{};
  }

  auto params() { return typename decltype(env())::params{}; }

  template <string_literal S>
  constexpr auto env_param()
  {
    return ground::get<pattern::param<S>>(self()->params()).value;
  }
};

template <typename Struct>
struct data_holder : item<> {
  using attributes =
      gather<attribute<"instance", some::specific<pointer<Struct>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) instance() { return attr<"instance">(*self()); };
  };
};

}  // namespace siconos::storage
