#pragma once

#include <functional>
#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <utility>

#include "boost/pfr.hpp"
#include "boost/pfr/core_name.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/some/some.hpp"

namespace siconos::storage::pattern {

using mp::append;
using mp::concat;
using mp::flatten;

// let rec
// https://stackoverflow.com/questions/2067988/recursive-lambda-functions-in-c11
template <typename F>
struct recursive {
  F f;
  template <typename... Ts>
  constexpr decltype(auto) operator()(Ts&&... ts) const {
    return f(std::ref(*this), std::forward<Ts>(ts)...);
  }

  template <typename... Ts>
  constexpr decltype(auto) operator()(Ts&&... ts) {
    return f(std::ref(*this), std::forward<Ts>(ts)...);
  }
};

template <typename F>
recursive(F) -> recursive<F>;

static auto const rec = [](auto f) constexpr { return recursive{std::move(f)}; };

static auto cons = []<typename A, typename... As>(A a, mp::tuple<As...> tpl) constexpr {
  return mp::insert(tpl, 0_c, a);
};

template <typename T, typename Tpl>
using cons_x = decltype(cons(T{}, Tpl{}));

using mp::transform;

template <typename T>
static constexpr auto instance = T{};

template <typename T>
struct hold {
  static constexpr auto value = T{};
};

template <typename U>
static auto contains_type = []<typename... Attrs>(gather<Attrs...>) constexpr -> bool {
  return (std::is_same_v<U, Attrs> || ...);
};

namespace must {
template <typename T, typename Tpl>
concept contains = contains_type<T>(instance<Tpl>);
}

namespace match {

template <typename T>
concept type_t = requires { typename T::type_t; };

template <typename T>
concept attribute = requires { typename T::attribute_t; };

template <typename T>
concept attribute_with_internal_type = attribute<T> && requires { typename T::type; };

template <typename T, typename I>
concept attached_storage = requires { typename T::tag; } && attribute<T> && item<I> &&
                           (std::derived_from<I, typename T::item> ||
                            std::derived_from<typename T::item, I>
                            // wrap case
                            || std::derived_from<typename I::type, typename T::item>);

template <typename T>
concept item_property = requires { typename T::item; } && std::derived_from<T, some::property>;

template <typename T, typename Tag>
concept tag = std::derived_from<typename T::tag, Tag>;

template <typename T>
concept item_ref = attribute<T> && item<typename T::type>;

template <typename T>
concept unbounded_storage = std::derived_from<T, some::unbounded_storage>;

template <typename T>
concept bounded_storage = std::derived_from<T, some::bounded_storage>;

template <typename T>
concept abstract_vector = std::derived_from<T, some::undefined_vector> ||
                          std::derived_from<T, some::undefined_unbounded_vector>;

template <typename T>
concept abstract_matrix = std::derived_from<T, some::undefined_matrix> ||
                          std::derived_from<T, some::undefined_unbounded_matrix>;

template <typename T>
concept property = std::derived_from<T, some::property>;

}  // namespace match

struct item {
  using item_t = void;

  using attributes = gather<>;

  template <typename H>
  struct interface {
    decltype(auto) self() {
      return static_cast<H*>(this);  // handle inherits from default_interface
    }
  };

  template <typename T>
  using methods = gather<>;
};

struct any_wrapper {};

template <template <typename... Ts> typename Wrapper, match::item Item, typename... Args>
struct wrap : Wrapper<Item, Args...>, Item, any_wrapper {
  using wrap_t = void;
  template <typename T>
  using wrapper = Wrapper<T, Args...>;
  //  using attributes = typename Item::attributes;
  using type = Item;
};

// Helper to build compile-time string from string_view
template <std::size_t N>
consteval auto string_view_to_fixed_string(std::string_view sv) {
  return [&sv]<std::size_t... Is>(std::index_sequence<Is...>) {
    const char arr[N] = {(Is < sv.size() ? sv[Is] : '\0')...};
    return string_literal<N>{arr};
  }(std::make_index_sequence<N>());
}

// Symbol type that stores PFR field name
template <typename T, std::size_t I>
struct pfr_symbol {
  static constexpr std::string_view name = boost::pfr::get_name<I, T>();
  // Convert to string_literal for compatibility
  static constexpr auto literal = [] {
    constexpr auto sv = name;
    constexpr std::size_t len = sv.size();
    // Build string_literal from constituent chars
    char arr[len + 1];
    for (std::size_t idx = 0; idx < len; ++idx) {
      arr[idx] = sv[idx];
    }
    arr[len] = '\0';
    return string_literal<len + 1>(arr);
  }();
};

// rename to attr_name ? (cf with_name)
template <string_literal Name, match::attribute A>
struct attribute : A, symbol<Name> {};

// Generate attribute from PFR field
namespace detail {
template <typename AttrStruct, std::size_t I>
struct pfr_field_attr {
  using field_type = boost::pfr::tuple_element_t<I, AttrStruct>;
  using symbol_info = pfr_symbol<AttrStruct, I>;

  // Create attribute with compile-time name from PFR
  using attribute_type = attribute<symbol_info::literal, field_type>;
};
}  // namespace detail

// Convert entire struct to gather<>
template <typename AttrStruct>
using struct_to_gather = decltype([]<std::size_t... Is>(std::index_sequence<Is...>) {
  return gather<typename detail::pfr_field_attr<AttrStruct, Is>::attribute_type...>{};
}(std::make_index_sequence<boost::pfr::tuple_size_v<AttrStruct>>{}));

// association for non nested type (should be the default now)
template <match::item Item, match::attribute A>
struct paired : A {
  using paired_t = void;
  using item = Item;
};

namespace match {
template <typename T>
concept paired = requires { typename T::paired_t; };
}  // namespace match

namespace match {
template <typename T>
concept named_item = item<T> && std::derived_from<T, any_symbol>;
}

template <match::item Item, match::attribute A>
static constexpr decltype(auto) named_attribute_maybe(Item, A) {
  // Always attach Item (or its underlying type) to make attributes unique per
  // item
  if constexpr (match::wrap<Item>) {
    return paired<typename Item::type, A>{};
  } else if constexpr (match::named_item<Item>) {
    return paired<typename Item::item, A>{};
  } else {
    return paired<Item, A>{};  // Always pair, never return A{}
  }
}

static auto attributes = []<match::item Item>(Item) constexpr -> decltype(auto) {
  if constexpr (std::is_aggregate_v<typename Item::attributes>) {
    auto unpaired_attrs = struct_to_gather<typename Item::attributes>{};
    return mp::transform(unpaired_attrs, [&]<typename A>(A) { return paired<Item, A>{}; });
  } else if constexpr (match::attributes<Item>) {
    return mp::transform(typename Item::attributes{}, [&]<match::attribute A>(A) {
      return named_attribute_maybe(Item{}, A{});
    });
  } else {
    return gather<>{};
  }
};

static auto properties = []<match::item Item>(Item) constexpr -> decltype(auto) {
  if constexpr (match::properties<Item>) {
    return typename Item::properties{};
  } else {
    return gather<>{};
  }
};

static auto is_a_ref =
    mp::is_a_model<[]<typename T>() consteval { return match::item_ref<T>; }>;

static auto is_a_poly_ref =
    mp::is_a_model<[]<typename T>() consteval { return match::polymorphic_type<T>; }>;

static auto all_items = rec([](auto&& all_items, match::item auto root_item) {
  using type_t = std::decay_t<decltype(root_item)>;

  auto items_ref = [&root_item]() {
    if constexpr (match::attributes<type_t>) {
      return transform(mp::filter(attributes(root_item), is_a_ref),
                       []<typename T>(T) { return typename T::type{}; });
    } else {
      return gather<>{};
    }
  };

  auto poly_ref = [&root_item]() {
    if constexpr (match::attributes<type_t>) {
      return transform(flatten(transform(mp::filter(attributes(root_item), is_a_poly_ref),
                                         []<typename T>(T) { return typename T::type{}; })),
                       []<typename T>(T) { return typename T::type{}; });
    } else {
      return gather<>{};
    }
  };

  if constexpr (match::items<type_t>) {
    return cons(root_item,
                flatten(transform(concat(items_ref(), typename type_t::items{}), all_items)));
  } else {
    return cons(root_item, flatten(transform(concat(items_ref(), poly_ref()), all_items)));
  }
});

static auto all_attributes = []<match::item Item>(Item) constexpr {
  return mp::tuple_unique(mp::concat_all(transform(all_items(Item{}), attributes)));
};

static auto all_properties = []<match::item Item>(Item) constexpr {
  return mp::tuple_unique(mp::concat_all(transform(all_items(Item{}), properties)));
};

namespace match {
template <typename I>
concept index = requires { typename I::index_t; };

template <typename H, typename A>
concept handle_attribute =
    attribute<A> && item<typename H::type> &&
    must::contains<A, decltype(siconos::storage::pattern::attributes(typename H::type{}))>;

template <typename T, typename I>
concept attribute_of = attribute<T> && item<I> &&
                       must::contains<T, decltype(siconos::storage::pattern::attributes(I{}))>;
}  // namespace match

template <string_literal S>
struct param : symbol<S> {};

template <auto V>
struct param_val {
  static constexpr auto value = V;
};

template <typename T>
struct param_type {
  using type = T;
};

namespace detail {
/**
 * @brief Index of the first item that owns Attr, as a plain constant-folded
 * scan over the pack.
 *
 * This is called once per storage column by storage::make, so it is on the
 * hot path of compilation. A recursive car/cdr walk costs O(N^2) compiler
 * memory per column (one fresh tuple type per step, plus the
 * reference_wrapper/__invoke_result trait train of the rec() combinator);
 * the bool array below instantiates nothing at all.
 */
template <typename Attr, typename... Items>
struct item_of {
  static constexpr std::size_t index = []() constexpr {
    constexpr bool found[] = {
        (match::attribute_of<Attr, Items> || match::attached_storage<Attr, Items>)...};
    for (std::size_t i = 0; i < sizeof...(Items); ++i) {
      if (found[i]) return i;
    }
    return sizeof...(Items);
  }();

  static_assert(index < sizeof...(Items), "item not found");

#if __has_builtin(__type_pack_element)
  using type = __type_pack_element<index, Items...>;
#else
  using type = std::tuple_element_t<index, std::tuple<Items...>>;
#endif
};
}  // namespace detail

template <match::attribute Attr>
static auto item_attribute = []<typename... Items>(mp::tuple<Items...>) constexpr {
  return typename detail::item_of<Attr, Items...>::type{};
};

static auto constexpr attribute_name(match::attribute auto a) { return a.str.value; };

template <typename... Ts>
std::tuple<Ts...> collect(Ts... ts) {
  return std::make_tuple(ts...);
};

template <typename F>
auto method(auto s, F f, auto doc) {
  return std::make_tuple(s, f, doc);
};

template <typename F>
auto method(auto s, F f) {
  return std::make_tuple(s, f);
};

namespace match {
template <typename T>
concept methods = requires(T h) { h.methods(); };
}  // namespace match

static auto constexpr method_name(auto m) { return std::get<0>(m); }

static auto constexpr method_def(auto m) { return std::get<1>(m); }

namespace match {
template <typename D>
concept store = requires(D d) { d.store(); };

template <typename Item>
concept without_attributes_bindings = requires { typename Item::without_attributes_bindings; };

template <typename Item>
concept without_attached_storages_bindings =
    requires { typename Item::without_attached_storages_bindings; };

template <typename Item>
concept batch_capable = requires { typename Item::batch_capable; };

template <typename Item>
concept without_batch_bindings = requires { typename Item::without_batch_bindings; };

template <typename Prop>
concept without_binding = requires { typename Prop::without_binding_t; };

template <typename Item>
concept empty_item = requires { typename Item::empty_item_t; };

template <typename Item, typename Symb>
concept with_attribute = requires {
  // Check that filter returns at least one element
  {
    mp::size(mp::filter(storage::pattern::attributes(Item{}), mp::derive_from<Symb>))
  } -> std::convertible_to<std::size_t>;
  requires(mp::size(mp::filter(storage::pattern::attributes(Item{}),
                               mp::derive_from<Symb>)) > mp::size_c<0_c>)
              .value;  // hana specific => hide under mp.hpp
};

}  // namespace match

struct empty_item : item {
  using empty_item_t = void;
};

template <match::item Item, string_literal S>
using attr_t =
    std::decay_t<decltype(mp::filter(attributes(Item{}), mp::derive_from<symbol<S>>)[0_c])>;

template <match::item Item, string_literal S>
  requires match::with_attribute<Item, symbol<S>>
using checked_attr_t = attr_t<Item, S>;

}  // namespace siconos::storage::pattern
