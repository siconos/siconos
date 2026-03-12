#pragma once

#include <siconos/storage/storage.hpp>
#include <tuple>
#include <variant>

namespace siconos::variant {

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
