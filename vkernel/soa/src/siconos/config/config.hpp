#pragma once

#include "siconos/storage/storage.hpp"

namespace siconos::config
{
  template<typename ...Cfs>
  using map = storage::mp::map<Cfs...>;

  template<typename K, typename V>
  using assoc = storage::mp::key_value<K, V>;

  using storage::pattern::param;
  using storage::pattern::param_val;
  using storage::pattern::param_type;

  template<storage::pattern::string_literal S, auto N>
  using iparam = assoc<param<S>, param_val<N>>;

  template<storage::pattern::string_literal S, typename T>
  using tparam = assoc<param<S>, param_type<T>>;
}
