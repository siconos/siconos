#pragma once

#include <any>

#include "siconos/storage/info.hpp"
#include "siconos/storage/memory.hpp"
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

template <typename T>
struct access {
  static constexpr auto at = mp::overload(
      []<typename Data, typename U = T,
         match::handle<decltype(item_attribute<U>(typename get_info_t<Data>::all_items_t{}))>
             Handle>(Handle h, Data& data) -> decltype(auto) {
        return siconos::storage::get<U>(h, data);
      },
      []<typename U = T, typename FullHandle>(FullHandle h) -> decltype(auto) {
        return siconos::storage::get<U>(h.data(), h);
      },
      []<typename U = T, typename FullHandle>(FullHandle h, typename FullHandle::indice step)
          -> decltype(auto) { return siconos::storage::get<U>(h.data(), step, h); });
};

template <string_literal S>
static auto param =
    [](auto h) constexpr -> decltype(auto) { return h.template env_param<S>(); };

/**
 * Where an attached storage is declared: the optional last template argument
 * of prop_values, prop and handle::property.
 *
 * - world_declared (the default): search every property of the world, as
 *   before;
 * - config_declared: declared in the configuration (with_properties<...>);
 *   searches the world for now, and the configuration once a world no longer
 *   carries its item closure;
 * - an item type: search only that item's own properties list, without the
 *   world.  This is what lets item code look up the attached storages it
 *   declares (topology, moreau_jean_element, io, ...) without the world.
 *
 * Compile with VKERNEL_CHECK_DECLARER to check, at every call site, that the
 * declarer's list gives the same storage as the world search.
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

#ifdef VKERNEL_CHECK_DECLARER
  static_assert(
      std::is_same_v<Decl, world_declared> ||
          std::is_same_v<type,
                         typename declared_attached_storage<I, S, world_declared, Data>::type>,
      "the declarer's attached storage differs from the world's");
#endif
};

template <match::item I, string_literal S, typename Decl = world_declared>
static auto prop_memory = [](auto& data) constexpr -> decltype(auto) {
  using attached_storage_t =
      typename declared_attached_storage<I, S, Decl, std::decay_t<decltype(data)>>::type;
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

// Runtime lookup: get or create the sparse_set for a dynamic property.
// The sparse_set is stored in the database-level _dynamic_properties map.
// storage::make database: the sparse set lives in a map keyed by name
template <match::item I, string_literal S, typename prop_t, typename env_t, typename indice>
static decltype(auto) prop_dynamic_map(auto& data) {
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

template <match::item I, string_literal S, typename Decl = world_declared>
static auto prop_dynamic_memory = [](auto& data) -> decltype(auto) {
  using data_t = std::decay_t<decltype(data)>;
  using info_t = get_info_t<data_t>;
  using env_t = typename info_t::template env<I>;
  using indice = typename env_t::indice;

  // Find the property to get its attribute type (Value)
  constexpr auto tpl =
      mp::filter(declared_properties_t<Decl, data_t>{}, is_dynamic_storage_tagged<I, S>);
  static_assert(mp::size(tpl) >= mp::size_c<1_c>, "dynamic property not found in topology");
  using prop_t = std::decay_t<decltype(tpl[0_c])>;
#ifdef VKERNEL_CHECK_DECLARER
  constexpr auto world_tpl = mp::filter(declared_properties_t<world_declared, data_t>{},
                                        is_dynamic_storage_tagged<I, S>);
  static_assert(std::is_same_v<prop_t, std::decay_t<decltype(world_tpl[0_c])>>,
                "the declarer's dynamic property differs from the world's");
#endif
  if constexpr (requires { data.store()._registry; }) {
    // registry world: the sparse column of the property, created at
    // registration and kept consistent by storage::remove
    return mp::get<prop_t>(data.store());
  } else {
    return prop_dynamic_map<I, S, prop_t, env_t, indice>(data);
  }
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
    // bare interface object.  It used to be called on handle_t{data_t{}, 0},
    // which built a whole temporary world (database, or registry with all its
    // columns) for every bound item, at compile time in every TU and at
    // import time.  methods() bodies may use self() in unevaluated operands
    // only.
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

}  // namespace siconos::storage
