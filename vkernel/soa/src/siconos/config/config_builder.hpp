#pragma once

#include "siconos/config/config.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/config/environment.hpp"

namespace siconos::config {

using storage::pattern::string_literal;
using storage::pattern::string_view_to_fixed_string;
namespace mp = storage::mp;

// Converts a hana::string to string_literal
template <typename S>
constexpr auto to_string_literal(S s)
{
  constexpr auto str = mp::to<char const*>(s);
  constexpr std::size_t len = mp::size(s);
  return string_view_to_fixed_string<len + 1>(std::string_view(str, len));
};

// core configuration builder
template <typename ParamsMap = boost::hana::map<>,
          typename ItemsList = mp::tuple<>,
          typename PropertiesList = mp::tuple<>>
class ConfigBuilder {
 public:
  // parameters
  template <typename Value>
  auto with_param(auto str, Value) const
  {
    constexpr auto str_li = to_string_literal(str);
    using KeyValue = assoc<param<str_li>, Value>;
    auto new_params = boost::hana::insert(ParamsMap{}, KeyValue{});

    return ConfigBuilder<decltype(new_params), ItemsList, PropertiesList>{};
  }

  template <string_literal str, auto Value>
  auto with_param() const
  {
    using KeyValue = assoc<param<str>, param_val<Value>>;
    auto new_params = boost::hana::insert(ParamsMap{}, KeyValue{});

    return ConfigBuilder<decltype(new_params), ItemsList, PropertiesList>{};
  }

  
  // dof
  template <auto Value>
  auto with_dof() const
  {
    return with_param<"dof", param_val<Value>>();
  }

  // items
  template <typename... ItemTypes>
  auto with_items() const
  {
    using NewItems =
        decltype(mp::concat(ItemsList{}, mp::tuple<ItemTypes...>{}));
    return ConfigBuilder<ParamsMap, NewItems, PropertiesList>{};
  }

  // properties
  template <typename Item, string_literal str>
  auto with_time_invariant() const
  {
    using Prop = storage::time_invariant<storage::attr_t<Item, str>>;
    using NewProps = decltype(mp::prepend(PropertiesList{}, Prop{}));
    return ConfigBuilder<ParamsMap, ItemsList, NewProps>{};
  }

  template <typename Item, string_literal str>
  auto with_diagonal() const
  {
    using Prop = storage::diagonal<storage::attr_t<Item, str>>;
    using NewProps = decltype(mp::prepend(PropertiesList{}, Prop{}));
    return ConfigBuilder<ParamsMap, ItemsList, NewProps>{};
  }

  template <typename OsiType, string_literal str>
  auto with_assembled_diagonal() const
  {
    using Prop = storage::assembled_diagonal<
        storage::attr_t<typename OsiType::assembled_osi_t, str>>;
    using NewProps = decltype(mp::prepend(PropertiesList{}, Prop{}));
    return ConfigBuilder<ParamsMap, ItemsList, NewProps>{};
  }

  // general chaining
  template <typename OtherBuilder>
  auto then(OtherBuilder) const
  {
    // Merge two builders
    using MergedParams = decltype(boost::hana::union_(
        ParamsMap{}, typename OtherBuilder::params_map{}));
    using MergedItems = decltype(mp::concat(
        ItemsList{}, typename OtherBuilder::items_list{}));
    using MergedProps = decltype(mp::concat(
        PropertiesList{}, typename OtherBuilder::properties_list{}));
    return ConfigBuilder<MergedParams, MergedItems, MergedProps>{};
  }

  // build
  auto build() const
  {
    // mp::map
    using ParamMap = decltype(mp::unpack(
        ParamsMap{}, []<typename... Pairs>(Pairs... pairs) {
          return mp::map<Pairs...>{};
        }));

    return mp::unpack(
        mp::concat(
            mp::concat(mp::tuple<standard_environment<ParamMap>>{},
                       ItemsList{}),
            mp::unpack(
                PropertiesList{},
                []<typename... Props>(Props...) {
                  return mp::tuple<storage::with_properties<Props...>>{};
                })),
        []<typename... Ts>(Ts... ts) { return storage::make<Ts...>{}; });
  }

  // Type aliases for merging
  using params_map = ParamsMap;
  using items_list = ItemsList;
  using properties_list = PropertiesList;
};

// Factory function
inline auto storage() { return ConfigBuilder<>{}; }

}  // namespace siconos::config
