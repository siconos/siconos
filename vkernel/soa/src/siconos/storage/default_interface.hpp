#pragma once

#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/mp/mp.hpp"

namespace siconos::storage {

template <typename Handle>
struct default_interface {
  decltype(auto) self()
  {
    return static_cast<Handle*>(
        this);  // handle inherits from default_interface
  };

  // Helper trait to check if item has runtime_properties attachment
  template <typename Item, typename Data>
  static constexpr bool has_runtime_properties_v =
      mp::any_of(typename get_info_t<Data>::all_properties_t{},
                 []<typename X>(X) -> bool {
                   return match::attached_storage<X, Item> &&
                          match::tag<X, symbol<"runtime_properties">>;
                 });

  // Conditionally enable operator[] only for items with runtime_properties
  template <typename H = Handle, typename Data = typename H::data_t,
            typename Item = typename H::type,
            bool HasProps = has_runtime_properties_v<Item, Data>,
            typename = std::enable_if_t<HasProps>>
  decltype(auto) operator[](const char* key)
  {
    return prop<"runtime_properties">(*self())[key];
  }

  template <typename Value, typename H = Handle,
            typename Data = typename H::data_t,
            typename Item = typename H::type,
            bool HasProps = has_runtime_properties_v<Item, Data>,
            typename = std::enable_if_t<HasProps>>
  Value& get(const char* key)
  {
    return prop<"runtime_properties">(*self()).template get<Value>(key);
  }

  auto env()
  {
    auto& data = self()->data();
    using info_t = get_info_t<decltype(data)>;
    return typename info_t::template env<typename Handle::type>{};
  }

  auto params() { return typename decltype(env())::params{}; }

  template <string_literal S>
  constexpr auto env_param()
  {
    return mp::get<pattern::param<S>>(self()->params()).value;
  }

  auto make_handle(auto& index)
  {
    return storage::make_handle(self()->data(), index);
  }
};

}  // namespace siconos::storage
