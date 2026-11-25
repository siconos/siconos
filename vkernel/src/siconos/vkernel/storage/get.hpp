// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/storage/info.hpp"
#include "siconos/vkernel/storage/memory.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/pattern/base.hpp"
#include "siconos/vkernel/storage/pattern/pattern.hpp"
#include "siconos/vkernel/storage/properties.hpp"
#include "siconos/vkernel/storage/sparse_set.hpp"
#include "siconos/vkernel/storage/traits/traits.hpp"

namespace siconos::vkernel::storage {

using namespace pattern;

template <typename A>
static auto get = mp::overload(
    // get<Attr>(data, step, handle)
    []<match::handle_attribute<A> Handle, match::store Data>(
        Data&& data, auto step, Handle&& handle) constexpr -> decltype(auto) {
      return memory(step,
                    mp::get<A>(static_cast<Data&&>(data).store()))[handle.index().value()];
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

template <string_literal S>
static auto param =
    [](auto h) constexpr -> decltype(auto) { return h.template env_param<S>(); };

/**
 * Where an attached storage is declared: the optional last template argument
 * of prop_values, prop and handle::property.
 *
 * - world_declared (the default): search the properties the world's handles
 *   know, i.e. its configuration (with_properties<...>);
 * - config_declared: the same;
 * - an item type: search only that item's own properties list.  This is how
 *   item code looks up the attached storages it declares (topology,
 *   system_kind, io, ...).
 */
struct world_declared {};
struct config_declared {};

template <typename Decl, typename Data>
struct declared_properties {
  using type = typename get_info_t<Data>::all_properties_t;
};

template <typename Decl, typename Data>
  requires match::item<Decl>
struct declared_properties<Decl, Data> {
  using type = std::decay_t<decltype(pattern::properties(Decl{}))>;
};

template <typename Decl, typename Data>
using declared_properties_t = typename declared_properties<Decl, Data>::type;

template <string_literal S, typename Decl = world_declared>
static auto prop =
    [](auto h) constexpr -> decltype(auto) { return h.template property<S, Decl>(); };

// the read access: no entry is created for an item without a dynamic property
template <string_literal S, typename Decl = world_declared>
static auto read_prop =
    [](auto h) constexpr -> decltype(auto) { return h.template read_property<S, Decl>(); };

template <string_literal S>
static auto attr =
    []<typename H>(H h, typename H::indice step = 0) constexpr -> decltype(auto) {
  using attr_n = attr_t<typename H::type, S>;
  return memory(step, mp::get<attr_n>(h.data().store()))[h.index().value()];
};

template <match::attribute T>
static constexpr decltype(auto) attr_memory(auto& data) {
  return mp::get<T>(data.store());
};

template <match::item I, string_literal S>
static constexpr decltype(auto) attr_memory(auto& data) {
  return mp::get<checked_attr_t<I, S>>(data.store());
};

template <string_literal S>
static constexpr auto is_identified_by =
    mp::is_a_model<[]<typename T>() constexpr { return match::tag<T, symbol<S>>; }>;

/// The (static) attached storage of item I with tag S, searched in the
/// properties declared by Decl (see world_declared).
template <match::item I, string_literal S, typename Decl, typename Data>
struct declared_attached_storage {
  static constexpr auto tpl =
      mp::filter(mp::filter(declared_properties_t<Decl, Data>{}, is_attached_storage<I>),
                 is_identified_by<S>);
  static constexpr bool found = mp::size(tpl) >= mp::size_c<1_c>;

  static constexpr auto first() {
    if constexpr (found) {
      return std::decay_t<decltype(tpl[0_c])>{};
    } else {
      []<bool flag = false>() { static_assert(flag, "attached storage not found"); }();
    }
  }

  using type = decltype(first());
};

template <match::item I, string_literal S, typename Decl = world_declared>
static auto prop_memory = [](auto& data) constexpr -> decltype(auto) {
  using data_t = std::decay_t<decltype(data)>;
  using attached_storage_t = typename declared_attached_storage<I, S, Decl, data_t>::type;
  return mp::get<attached_storage_t>(data.store());
};

template <match::attribute T>
static constexpr decltype(auto) attr_values(auto& data, auto step) {
  return memory(step, (attr_memory<T>(data)));
};

template <match::item I, string_literal S, typename D>
  requires match::with_attribute<I, symbol<S>>
static constexpr decltype(auto) attr_values(D&& data, auto step) {
  return memory(step, (attr_memory<I, S>(data)));
};

// Compile-time check: is there a property for Item with tag S?
template <typename Item, typename Data, string_literal S>
static constexpr bool has_property_v =
    mp::any_of(typename get_info_t<Data>::all_properties_t{}, []<typename X>(X) -> bool {
      return match::property<X> && match::tag<X, symbol<S>> &&
             match::attached_storage<X, Item>;
    });

// Compile-time check: is there a dynamic_storage property for Item with tag S
// (among the properties declared by Decl)?
template <typename Item, typename Data, string_literal S, typename Decl = world_declared>
static constexpr bool is_dynamic_storage_v =
    mp::any_of(declared_properties_t<Decl, Data>{}, []<typename X>(X) -> bool {
      return match::property<X> && match::tag<X, symbol<S>> &&
             match::attached_storage<X, Item> && requires { typename X::dynamic_storage_t; };
    });

template <typename Item, string_literal S>
static constexpr auto is_dynamic_storage_tagged = mp::is_a_model<[]<typename T>() constexpr {
  return match::property<T> && match::tag<T, symbol<S>> && match::attached_storage<T, Item> &&
         requires { typename T::dynamic_storage_t; };
}>;

template <typename T, typename Tag>
concept tagged_exactly = std::same_as<typename T::tag, Tag>;

// a dynamic property of Item whose tag is Tag, which may be any type (a
// symbol<S> is matched exactly, not as a base)
template <typename Item, typename Tag>
static constexpr auto is_dynamic_storage_with_tag = mp::is_a_model<[]<typename T>() constexpr {
  return match::property<T> && match::attached_storage<T, Item> &&
         requires { typename T::dynamic_storage_t; } && tagged_exactly<T, Tag>;
}>;

template <match::item I, string_literal S, typename Decl = world_declared>
static auto prop_dynamic_memory = [](auto& data) -> decltype(auto) {
  using data_t = std::decay_t<decltype(data)>;

  // Find the property to get its attribute type (Value)
  constexpr auto tpl =
      mp::filter(declared_properties_t<Decl, data_t>{}, is_dynamic_storage_tagged<I, S>);
  static_assert(mp::size(tpl) >= mp::size_c<1_c>, "dynamic property not found in topology");
  using prop_t = std::decay_t<decltype(tpl[0_c])>;
  // the sparse column of the property, created at registration and kept
  // consistent by storage::remove
  return mp::get<prop_t>(data.store());
};

/// prop_dynamic_memory for the dynamic property of I tagged by a type Tag
/// rather than by a name
template <match::item I, typename Tag, typename Decl = world_declared>
static auto prop_dynamic_memory_tagged = [](auto& data) -> decltype(auto) {
  using data_t = std::decay_t<decltype(data)>;

  constexpr auto tpl =
      mp::filter(declared_properties_t<Decl, data_t>{}, is_dynamic_storage_with_tag<I, Tag>);
  static_assert(mp::size(tpl) >= mp::size_c<1_c>, "dynamic property not found in declarer");
  using prop_t = std::decay_t<decltype(tpl[0_c])>;
  return mp::get<prop_t>(data.store());
};

template <match::item I, string_literal S, typename Decl = world_declared>
static auto prop_values = [](auto& data, auto step) -> decltype(auto) {
  if constexpr (is_dynamic_storage_v<I, std::decay_t<decltype(data)>, S, Decl>) {
    return prop_dynamic_memory<I, S, Decl>(data);
  } else {
    return memory(step, (prop_memory<I, S, Decl>(data)));
  }
};

// fix H is a handle defined in storage
template <typename Hc>
static auto constexpr methods(Hc hc) {
  using handle_t = typename decltype(+hc)::type;

  if constexpr (match::methods<handle_t>) {
    // the (name, &interface::member) table needs no handle state: call it on a
    // bare interface object.  Calling it on a handle_t{data_t{}, 0} would
    // build a whole temporary world (registry with all its columns) for every
    // bound item, at compile time in every TU and at import time.  methods()
    // bodies may use self() in unevaluated operands only.
    return typename handle_t::type::template interface<handle_t>{}.methods();
  } else {
    return gather<>{};
  }
}

template <typename Astor>
static auto constexpr attached_storage_name(Astor astor) {
  using tag = typename Astor::tag;
  if constexpr (match::symbol<tag>) {
    return tag{};
  } else {
    []<bool flag = false>() {
      mp::type_trace<tag>();
      static_assert(flag, "tag is not a symbol!");
    }();
  }
}

template <typename Data, typename Attribute>
static auto constexpr get_storage_type(Data&& data, Attribute) {
  using val_store_t = std::decay_t<decltype(mp::get<Attribute>(data.store()))>;
  using val_t = typename val_store_t::value_type::value_type;
  return val_t{};
}

template <typename System, typename Data, typename Attribute>
static auto constexpr convert_storage_type(System, Data&& data, Attribute) {
  using info_t = get_info_t<Data>;

  return typename traits::config<typename info_t::template env<System>>::template convert<
      Attribute>::type{};
}

}  // namespace siconos::vkernel::storage
