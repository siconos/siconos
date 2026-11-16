/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "siconos/vkernel/config/config.hpp"
#include "siconos/vkernel/config/environment.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/pattern/base.hpp"
#include "siconos/vkernel/storage/properties.hpp"
#include "siconos/vkernel/storage/registry.hpp"

namespace siconos::vkernel::config {

using storage::pattern::string_literal;
using storage::pattern::string_view_to_fixed_string;
namespace mp = storage::mp;

// Converts a hana::string to string_literal
template <typename S>
constexpr auto to_string_literal(S s) {
  constexpr auto str = mp::to<char const*>(s);
  constexpr std::size_t len = mp::size(s);
  return string_view_to_fixed_string<len + 1>(std::string_view(str, len));
};

// core configuration builder
template <typename ParamsMap = boost::hana::map<>, typename ItemsList = mp::tuple<>,
          typename PropertiesList = mp::tuple<>>
class config_builder {
 public:
  // parameters
  template <typename Value>
  auto with_param(auto str, Value) const {
    constexpr auto str_li = to_string_literal(str);
    using key_value = assoc<param<str_li>, Value>;
    auto new_params = boost::hana::insert(ParamsMap{}, key_value{});

    return config_builder<decltype(new_params), ItemsList, PropertiesList>{};
  }

  template <string_literal str, auto Value>
  auto with_param() const {
    using key_value = assoc<param<str>, param_val<Value>>;
    auto new_params = boost::hana::insert(ParamsMap{}, key_value{});

    return config_builder<decltype(new_params), ItemsList, PropertiesList>{};
  }

  // dof
  template <auto Value>
  auto with_dof() const {
    return with_param<"dof", param_val<Value>>();
  }

  // items
  template <typename... ItemTypes>
  auto with_items() const {
    using NewItems = decltype(mp::concat(ItemsList{}, mp::tuple<ItemTypes...>{}));
    return config_builder<ParamsMap, NewItems, PropertiesList>{};
  }

  // properties
  template <typename Item, string_literal str>
  auto with_time_invariant() const {
    using Prop = storage::time_invariant<storage::attr_t<Item, str>>;
    using NewProps = decltype(mp::prepend(PropertiesList{}, Prop{}));
    return config_builder<ParamsMap, ItemsList, NewProps>{};
  }

  template <typename Item, string_literal str>
  auto with_diagonal() const {
    using Prop = storage::diagonal<storage::attr_t<Item, str>>;
    using NewProps = decltype(mp::prepend(PropertiesList{}, Prop{}));
    return config_builder<ParamsMap, ItemsList, NewProps>{};
  }

  template <typename OsiType, string_literal str>
  auto with_assembled_diagonal() const {
    using Prop = storage::assembled_diagonal<storage::attr_t<OsiType, str>>;
    using NewProps = decltype(mp::prepend(PropertiesList{}, Prop{}));
    return config_builder<ParamsMap, ItemsList, NewProps>{};
  }

  // general chaining
  template <typename OtherBuilder>
  auto then(OtherBuilder) const {
    // Merge two builders
    using MergedParams =
        decltype(boost::hana::union_(ParamsMap{}, typename OtherBuilder::params_map{}));
    using MergedItems = decltype(mp::concat(ItemsList{}, typename OtherBuilder::items_list{}));
    using MergedProps =
        decltype(mp::concat(PropertiesList{}, typename OtherBuilder::properties_list{}));
    return config_builder<MergedParams, MergedItems, MergedProps>{};
  }

  template <typename T>
  struct env : standard_environment<T> {
    using params = decltype(mp::unpack(
        ParamsMap{}, []<typename... Pairs>(Pairs... pairs) { return mp::map<Pairs...>{}; }));
  };

  // the runtime (registry-backed) database type (storage/registry.hpp),
  // computed without constructing one
  using world_type = typename decltype(mp::unpack(
      mp::concat(mp::concat(mp::tuple<>{}, ItemsList{}),
                 mp::unpack(PropertiesList{},
                            []<typename... Props>(Props...) {
                              return mp::tuple<storage::with_properties<Props...>>{};
                            })),
      []<typename... Ts>(Ts...) {
        return std::type_identity<storage::make_runtime_database<env, Ts...>>{};
      }))::type;

  // build the runtime database (compile-time typed per column, registered at
  // run time).  The return type is explicit: struct data_t :
  // decltype(storage()...build_runtime_database()) {} must not instantiate
  // this body, whose construction registers the database's columns and so
  // needs the whole item closure; only the TU that constructs data_t does
  // that.
  world_type build_runtime_database() const { return world_type{}; }

  // previous name of build_runtime_database, kept so existing call sites
  // compile
  world_type build_world() const { return build_runtime_database(); }

  // Type aliases for merging
  using params_map = ParamsMap;
  using items_list = ItemsList;
  using properties_list = PropertiesList;
};

// Factory function
inline auto storage() { return config_builder<>{}; }

}  // namespace siconos::vkernel::config
