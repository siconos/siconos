#pragma once

#include "siconos/storage/some/some.hpp"
#include "siconos/storage/default_interface.hpp"

namespace siconos::storage {


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
