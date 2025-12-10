#pragma once

#include <siconos/storage/storage.hpp>
#include <tuple>
#include <variant>

namespace siconos::variant {

template <typename Variant, typename V, typename F>
  requires requires {
             {
               std::variant_size_v<Variant> + 0
               } -> std::same_as<std::size_t>;
           }
V apply_if_valid(Variant& var, V&& default_value, F&& fun)
{
  return siconos::storage::ground::call_with_index<
      std::variant_size_v<Variant>, V>(
      var.index(), static_cast<V&&>(default_value), [&var, &fun](auto i) {
        return static_cast<F&&>(fun)(std::get<i>(var));
      });
}

  template <typename D, typename V, typename F>
  decltype(auto) visit(D& data, V&& var, F&& fun)
{
  return std::visit(
    [&data, &fun](auto& rvar) {
      auto h = storage::make_handle(data, rvar);
        return static_cast<F&&>(fun)(h);
      },
      static_cast<V&&>(var));
}
}  // namespace siconos::variant
