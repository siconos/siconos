#pragma once

#include <any>

#include "siconos/storage/info.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/sparse_set.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::storage {

using namespace pattern;

template <typename A>
static auto get = mp::overload(
    // get<Attr>(data, step, handle)
    []<match::handle_attribute<A> Handle, match::store Data>(
        Data&& data, auto step, Handle&& handle) constexpr -> decltype(auto) {
      return memory(
          step,
          mp::get<A>(
              static_cast<Data&&>(data).store()))[handle.index().value()];
    },
    // get<Attr>(data, step, handle)
    []<match::handle_attribute<A> Handle, match::store Data>(
        Data& data, auto step, Handle& handle) constexpr -> decltype(auto) {
      return memory(step, mp::get<A>(data.store()))[handle.index().value()];
    },
    // get<Attr>(data, handle)
    []<match::handle_attribute<A> Handle, match::store Data>(
        Data& data, Handle& handle) constexpr -> decltype(auto) {
      return memory(0, mp::get<A>(data.store()))[handle.index().value()];
    });

template <typename T>
struct access {
  static constexpr auto at = mp::overload(
      []<typename Data, typename U = T,
         match::handle<decltype(item_attribute<U>(
             typename get_info_t<Data>::all_items_t{}))>
             Handle>(Handle h, Data& data) -> decltype(auto) {
        return siconos::storage::get<U>(h, data);
      },
      []<typename U = T, typename FullHandle>(FullHandle h)
          -> decltype(auto) { return siconos::storage::get<U>(h.data(), h); },
      []<typename U = T, typename FullHandle>(
          FullHandle h, typename FullHandle::indice step) -> decltype(auto) {
        return siconos::storage::get<U>(h.data(), step, h);
      });
};

template <string_literal S>
static auto param = [](auto h) constexpr -> decltype(auto) {
  return h.template env_param<S>();
};

template <string_literal S>
static auto prop = [](auto h) constexpr -> decltype(auto) {
  return h.template property<S>();
};

template <string_literal S>
static auto attr = []<typename H>(H h, typename H::indice step =
                                           0) constexpr -> decltype(auto) {
  using attr_n = attr_t<typename H::type, S>;
  return memory(step, mp::get<attr_n>(h.data().store()))[h.index().value()];
};

template <match::attribute T>
static constexpr decltype(auto) attr_memory(auto& data)
{
  return mp::get<T>(data.store());
};

template <match::item I, string_literal S>
static constexpr decltype(auto) attr_memory(auto& data)
{
  return mp::get<checked_attr_t<I, S>>(data.store());
};

template <string_literal S>
static constexpr auto is_identified_by =
    mp::is_a_model<[]<typename T>() constexpr {
      return match::tag<T, symbol<S>>;
    }>;

template <match::item I, string_literal S>
static auto prop_memory = [](auto& data) constexpr -> decltype(auto) {
  using info_t = get_info_t<decltype(data)>;
  constexpr auto tpl = mp::filter(
      mp::filter(typename info_t::all_properties_t{}, is_attached_storage<I>),
      is_identified_by<S>);

  if constexpr (mp::size(tpl) >= mp::size_c<1_c>) {
    using attached_storage_t = std::decay_t<decltype(tpl[0_c])>;
    return mp::get<attached_storage_t>(data.store());
  } else {
    []<bool flag = false>() {
      static_assert(flag, "attached storage not found");
    }();
  }
};

template <match::attribute T>
static constexpr decltype(auto) attr_values(auto& data, auto step)
{
  return memory(step, (attr_memory<T>(data)));
};

template <match::item I, string_literal S, typename D>
requires match::with_attribute<I, symbol<S>>
static constexpr decltype(auto) attr_values(D&& data, auto step)
{
  return memory(step, (attr_memory<I, S>(data)));
};

// Compile-time check: is there a property for Item with tag S?
template <typename Item, typename Data, string_literal S>
static constexpr bool has_property_v =
    mp::any_of(typename get_info_t<Data>::all_properties_t{},
               []<typename X>(X) -> bool {
                 return match::property<X> && match::tag<X, symbol<S>> &&
                        match::attached_storage<X, Item>;
               });

// Compile-time check: is there a dynamic_storage property for Item with tag S?
template <typename Item, typename Data, string_literal S>
static constexpr bool is_dynamic_storage_v =
    mp::any_of(typename get_info_t<Data>::all_properties_t{},
               []<typename X>(X) -> bool {
                 return match::property<X> &&
                        match::tag<X, symbol<S>> &&
                        match::attached_storage<X, Item> &&
                        requires { typename X::dynamic_storage_t; };
               });

template <typename Item, string_literal S>
static constexpr auto is_dynamic_storage_tagged =
    mp::is_a_model<[]<typename T>() constexpr {
      return match::property<T> && match::tag<T, symbol<S>> &&
             match::attached_storage<T, Item> &&
             requires { typename T::dynamic_storage_t; };
    }>;

// Runtime lookup: get or create the sparse_set for a dynamic property.
// The sparse_set is stored in the database-level _dynamic_properties map.
template <match::item I, string_literal S>
static auto prop_dynamic_memory = [](auto& data) -> decltype(auto) {
  using info_t = get_info_t<decltype(data)>;
  using env_t = typename info_t::template env<I>;
  using indice = typename env_t::indice;

  // Find the property to get its attribute type (Value)
  constexpr auto tpl = mp::filter(
      typename info_t::all_properties_t{}, is_dynamic_storage_tagged<I, S>);
  static_assert(mp::size(tpl) >= mp::size_c<1_c>,
                "dynamic property not found in topology");
  using prop_t = std::decay_t<decltype(tpl[0_c])>;
  using attr_t = typename prop_t::type;
  using value_t = typename traits::config<env_t>::template convert<attr_t>::type;

  // Build a unique key for this (Item, S) pair
  std::string key = std::string(typeid(I).name());
  key += ":";
  key += S.value;

  auto& map = data.store()._dynamic_properties;
  auto it = map.find(key);
  if (it == map.end()) {
    storage::sparse_set<indice, value_t> set;
    auto [new_it, _] = map.emplace(key, std::move(set));
    return std::any_cast<storage::sparse_set<indice, value_t>&>(new_it->second);
  }
  return std::any_cast<storage::sparse_set<indice, value_t>&>(it->second);
};

 template <match::item I, string_literal S>
 static auto prop_values = [](auto& data, auto step) -> decltype(auto) {
   if constexpr (is_dynamic_storage_v<I, std::decay_t<decltype(data)>, S>) {
     return prop_dynamic_memory<I, S>(data);
   } else {
     return memory(step, (prop_memory<I, S>(data)));
   }
 };

// fix H is a handle defined in storage
template <typename Hc>
static auto constexpr methods(Hc hc)
{
  using handle_t = typename decltype(+hc)::type;
  using data_t = typename handle_t::data_t;

  if constexpr (match::methods<handle_t>) {
    return handle_t{data_t{}, 0UL}.methods();
  }
  else {
    return gather<>{};
  }
}

template <typename Astor>
static auto constexpr attached_storage_name(Astor astor)
{
  using tag = typename Astor::tag;
  if constexpr (match::symbol<tag>) {
    return tag{};
  }
  else {
    []<bool flag = false>() {
      mp::type_trace<tag>();
      static_assert(flag, "tag is not a symbol!");
    }();
  }
}

template <typename Data, typename Attribute>
static auto constexpr get_storage_type(Data&& data, Attribute)
{
  using val_store_t =
      std::decay_t<decltype(mp::get<Attribute>(data.store()))>;
  using val_t = typename val_store_t::value_type::value_type;
  return val_t{};
}

template <typename System, typename Data, typename Attribute>
static auto constexpr convert_storage_type(System, Data&& data, Attribute)
{
  using info_t = get_info_t<Data>;

  return typename traits::config<typename info_t::template env<System>>::
      template convert<Attribute>::type{};
}

}  // namespace siconos::storage
