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
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/some/some.hpp"

namespace siconos::storage::pattern {

using mp::append;
using mp::concat;
using mp::flatten;
// namespace some = storage::some;
// namespace ground = storage::mp;

// let rec
// https://stackoverflow.com/questions/2067988/recursive-lambda-functions-in-c11
template <typename F>
struct recursive {
  F f;
  template <typename... Ts>
  constexpr decltype(auto) operator()(Ts&&... ts) const
  {
    return f(std::ref(*this), std::forward<Ts>(ts)...);
  }

  template <typename... Ts>
  constexpr decltype(auto) operator()(Ts&&... ts)
  {
    return f(std::ref(*this), std::forward<Ts>(ts)...);
  }
};

template <typename F>
recursive(F) -> recursive<F>;

static auto const rec = [](auto f) constexpr {
  return recursive{std::move(f)};
};

static auto proj = [](auto& data) constexpr -> decltype(auto) {
  return ([&data](auto&& fun) constexpr -> decltype(auto) {
    return ([&data,
             &fun]<typename... As>(As&&... args) constexpr -> decltype(auto) {
      return fun(std::forward<As>(args)..., data);
    });
  });
};

static auto compose = [](auto&& f, auto&& g) constexpr -> decltype(auto) {
  return [&f, &g]<typename Arg>(Arg arg) constexpr -> decltype(auto) {
    return f(g(std::forward<Arg>(arg)));
  };
};

static auto car = []<typename Tpl>(Tpl tpl) constexpr {
  static_assert(mp::size(Tpl{}) > mp::size_c<0_c>);
  return tpl[0_c];
};

// static_assert(std::is_same_v<decltype(car(gather<int, float, char>{})),
// int>);

static auto cdr =
    []<typename A0, typename... As>(mp::tuple<A0, As...> tpl) constexpr {
      return mp::drop_front(tpl, mp::size_c<1_c>);
    };

// static_assert(std::is_same_v<decltype(cdr(gather<int, float, char>{})),
//                              gather<float, char>>);

static auto cons = []<typename A, typename... As>(
                       A a, mp::tuple<As...> tpl) constexpr {
  return mp::insert(tpl, 0_c, a);
};

template <typename T, typename Tpl>
using cons_x = decltype(cons(T{}, Tpl{}));

// static_assert(std::is_same_v<decltype(cons(int{}, gather<float, char>{})),
//                              gather<int, float, char>>);

// static auto append = []<concepts::tuple_like... Tpls>(Tpls... tpls)
// constexpr
// {
//   return std::tuple_cat(tpls...);
// };

// static auto flatten = []<concepts::tuple_like Tpl>(Tpl tpl) constexpr {
//   if constexpr (std::tuple_size_v < Tpl >> 0) {
//     if constexpr (concepts::tuple_like<decltype(car(tpl))>) {
//       return std::apply(append, tpl);
//     }
//     else {
//       return tpl;
//     }
//   }
//   else {
//     return tpl;
//   }
// };

using mp::transform;

// template <typename F>
// static auto filter = []<typename Tpl>(const Tpl tpl) constexpr {
//   auto loop = rec([]<typename Itpl>(auto&& loop, Itpl) constexpr {
//     auto v = car(Itpl{});
//     using v_t = std::decay_t<decltype(v)>;
//     if constexpr (F{}.value(v_t{})) {
//       if constexpr (std::tuple_size_v < Itpl >> 1) {
//         return cons(v, loop(cdr(Itpl{})));
//       }
//       else {
//         return std::tuple<v_t>(v);
//       }
//     }
//     else if constexpr (std::tuple_size_v<Itpl> > 1) {
//       return loop(cdr(Itpl{}));
//     }
//     else {
//       return std::tuple<>{};
//     }
//   });

//   if constexpr (std::tuple_size_v < Tpl >> 0) {
//     return loop(Tpl{});
//   }
//   else {
//     return std::tuple<>{};
//   }
// };

template <std::size_t N>
struct size {
  using size_t = void;
  static constexpr std::size_t value = N;
};

template <std::size_t N>
struct degrees_of_freedom : size<N> {
  using degrees_of_freedom_t = void;
};

template <typename T>
static constexpr auto instance = T{};

template <typename T>
struct hold {
  static constexpr auto value = T{};
};
template <typename T>
static constexpr auto make_instance(T&&)
{
  return instance<T>;
};

template <typename U>
static auto contains_type =
    []<typename... Attrs>(gather<Attrs...>) constexpr -> bool {
  return (std::is_same_v<U, Attrs> || ...);
};

namespace must {
template <typename T, typename Tpl>
concept contains = contains_type<T>(instance<Tpl>);
}

// static_assert(match::size<std::vector<double>>);
// static_assert(match::size<std::array<double, 1>>);
// static_assert(match::push_back<std::vector<double>>);
// static_assert(!match::push_back<std::array<double, 3>>);

namespace match {

template <typename T>
concept type_t = requires { typename T::type_t; };

template <typename T>
concept attribute = requires { typename T::attribute_t; };

template <typename T>
concept attribute_with_internal_type =
    attribute<T> && requires { typename T::type; };

template <typename T>
concept attribute_or_item = attribute<T> || item<T>;

template <typename T, typename I>
concept attached_storage =
    requires { typename T::tag; } && attribute<T> && item<I> &&
    (std::derived_from<I, typename T::item> ||
     std::derived_from<typename T::item, I>
     // wrap case
     || std::derived_from<typename I::type, typename T::item>);

template <typename T>
concept item_property =
    requires { typename T::item; } && std::derived_from<T, some::property>;

template <typename T, typename Tag>
concept tag = std::derived_from<typename T::tag, Tag>;

template <typename T>
concept item_ref = attribute<T> && item<typename T::type>;

template <typename T>
concept unbounded_storage = std::derived_from<T, some::unbounded_storage>;

template <typename T>
concept bounded_storage = std::derived_from<T, some::bounded_storage>;

template <typename T>
concept abstract_vector =
    std::derived_from<T, some::undefined_vector> ||
    std::derived_from<T, some::undefined_unbounded_vector>;

template <typename T>
concept abstract_matrix =
    std::derived_from<T, some::undefined_matrix> ||
    std::derived_from<T, some::undefined_unbounded_matrix>;

template <typename T, std::size_t N>
concept cvector = requires(T a) { a[N - 1]; };

template <typename T>
concept property = std::derived_from<T, some::property>;

template <typename T, typename K>
concept property_of = property<T> && property<K> && std::derived_from<T, K>;

template <typename T, typename Ks>
concept any_of_property = mp::any_of(
    Ks{}, []<match::property K>(K) { return std::derived_from<T, K>; });

}  // namespace match

template <typename T>
struct tag {
  using type = T;
};

// concept item
template <match::item T>
struct use {
  using use_t = void;
  using type = T;
};

template <typename T>
struct attribute_p {
  static constexpr auto value = match::attribute<T>;
};

template <typename T>
struct degrees_of_freedom_p {
  static constexpr auto value = match::degrees_of_freedom<T>;
};

template <typename... Args>
struct frame {
  using args = gather<Args...>;

  static constexpr std::size_t dof =
      mp::find_if(args{},
                  []<typename T>(T) { return degrees_of_freedom_p<T>{}; })
          .value_or([]<bool flag = false>() {
            static_assert(flag, "need some dof");
          })
          .value;
};

struct item {
  using item_t = void;

  using attributes = gather<>;

  template <typename H>
  struct interface {
    decltype(auto) self()
    {
      return static_cast<H*>(this);  // handle inherits from default_interface
    }
  };

  template <typename T>
  using methods = gather<>;
};

struct any_wrapper {};

struct any_bounded_wrapper : any_wrapper {};
struct any_unbounded_wrapper : any_wrapper {};

template <template <typename... Ts> typename Wrapper, match::item Item,
          typename... Args>
struct wrap : Wrapper<Item, Args...>, Item, any_wrapper {
  using wrap_t = void;
  template <typename T>
  using wrapper = Wrapper<T, Args...>;
  //  using attributes = typename Item::attributes;
  using type = Item;
};

// namespace match {
//   template<typename T>
//   concept wrap = requires { typename T::wrap_t; };
// }

template <typename T>
struct place_holder {
  using type = std::array<T, 1>;
};

namespace concepts {
// T is a tag
template <typename T, typename Data>
concept vertex_item_t = requires(T t) {
  { static_cast<typename Data::vertex_items>(t) };
};
}  // namespace concepts

// Helper to build compile-time string from string_view
template <std::size_t N>
consteval auto string_view_to_fixed_string(std::string_view sv)
{
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
  using symbol_type = text<literal>;
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
using struct_to_gather = decltype([]<std::size_t... Is>(
                                      std::index_sequence<Is...>) {
  return gather<
      typename detail::pfr_field_attr<AttrStruct, Is>::attribute_type...>{};
}(std::make_index_sequence<boost::pfr::tuple_size_v<AttrStruct>>{}));

// association for non nested type (should be the default now)
template <match::item Item, match::attribute A>
struct paired : A {
  using paired_t = void;
  using item = Item;
};

namespace match {
template <typename PairedA, typename PairedB>
concept paired_similar =
    std::derived_from<typename PairedA::type, typename PairedB::type> ||
    std::derived_from<typename PairedB::type, typename PairedA::type>;

template <typename T>
concept paired = requires { typename T::paired_t; };
}  // namespace match

namespace must {

template <typename T, typename Tpl>
concept contains_similar_attribute =
    mp::any_of(Tpl{}, []<match::attribute A>(A) -> bool {
      if constexpr (match::paired<A> && match::paired<T>) {
        return std::derived_from<typename T::item, typename A::item> ||
               std::derived_from<typename A::item, typename T::item>;
      }
      else {
        return std::is_same_v<T, A>;
      }
    })();

}  // namespace must

namespace match {
template <typename T>
concept named_item = item<T> && std::derived_from<T, any_symbol>;
}

template <match::item Item, match::attribute A>
static constexpr decltype(auto) named_attribute_maybe(Item, A)
{
  // Always attach Item (or its underlying type) to make attributes unique per
  // item
  if constexpr (match::wrap<Item>) {
    return paired<typename Item::type, A>{};
  }
  else if constexpr (match::named_item<Item>) {
    return paired<typename Item::item, A>{};
  }
  else {
    return paired<Item, A>{};  // Always pair, never return A{}
  }
}

static auto attributes =
    []<match::item Item>(Item) constexpr -> decltype(auto) {
  if constexpr (std::is_aggregate_v<typename Item::attributes>) {
    auto unpaired_attrs = struct_to_gather<typename Item::attributes>{};
    return mp::transform(unpaired_attrs,
                         [&]<typename A>(A) { return paired<Item, A>{}; });
  }
  else if constexpr (match::attributes<Item>) {
    return mp::transform(typename Item::attributes{},
                         [&]<match::attribute A>(A) {
                           return named_attribute_maybe(Item{}, A{});
                         });
  }
  else {
    return gather<>{};
  }
};

template <typename Attrs, string_literal S>
using get_attr_t = std::decay_t<decltype(mp::filter(
    Attrs{}, mp::derive_from<symbol<S>>)[0_c])>;

template <match::item Item, string_literal S>
using attr_t = std::decay_t<decltype(mp::filter(
    attributes(Item{}), mp::derive_from<symbol<S>>)[0_c])>;

static auto properties =
    []<match::item Item>(Item) constexpr -> decltype(auto) {
  if constexpr (match::properties<Item>) {
    return typename Item::properties{};
  }
  else {
    return gather<>{};
  }
};

static auto is_a_ref =
    mp::is_a_model<[]<typename T>() consteval { return match::item_ref<T>; }>;

static auto is_a_poly_ref = mp::is_a_model<[]<typename T>() consteval {
  return match::polymorphic_type<T>;
}>;

static auto all_items = rec([](auto&& all_items, match::item auto root_item) {
  using type_t = std::decay_t<decltype(root_item)>;

  auto items_ref = [&root_item]() {
    if constexpr (match::attributes<type_t>) {
      return transform(mp::filter(attributes(root_item), is_a_ref),
                       []<typename T>(T) { return typename T::type{}; });
    }
    else {
      return gather<>{};
    }
  };

  auto poly_ref = [&root_item]() {
    if constexpr (match::attributes<type_t>) {
      return transform(flatten(transform(
                           mp::filter(attributes(root_item), is_a_poly_ref),
                           []<typename T>(T) { return typename T::type{}; })),
                       []<typename T>(T) { return typename T::type{}; });
    }
    else {
      return gather<>{};
    }
  };

  if constexpr (match::items<type_t>) {
    return cons(root_item, flatten(transform(
                               concat(items_ref(), typename type_t::items{}),
                               all_items)));
  }
  else {
    return cons(root_item, flatten(transform(concat(items_ref(), poly_ref()),
                                             all_items)));
    ;
  }
});

static auto all_attributes = []<match::item Item>(Item) constexpr {
  return mp::tuple_unique(
      mp::concat_all(transform(all_items(Item{}), attributes)));
};

static auto all_properties = []<match::item Item>(Item) constexpr {
  return mp::tuple_unique(
      mp::concat_all(transform(all_items(Item{}), properties)));
};

//  template<typename K>
//  static auto all_items_of_property = [](match::item auto&& t)
//    constexpr -> decltype(auto)
//  {
//    return filter<hold<decltype([]<typename T>(T) { return
//    match::property<T,K>; })>>
//      (all_items(t));
//  };

namespace match {
template <typename I>
concept index = requires { typename I::index_t; };

template <typename H, typename A>
concept handle_attribute =
    attribute<A> && item<typename H::type> &&
    must::contains<A, decltype(siconos::storage::pattern::attributes(
                          typename H::type{}))>;

template <typename H, typename A>
concept handle_attached_storage =
    item<typename H::type> &&
    std::tuple_size_v<
        std::decay_t<decltype(filter<hold<decltype([]<typename T>(T) {
          return (match::attached_storage<T, typename H::type> &&
                  match::tag<T, A>);
        })>>(typename H::info_t::all_properties_t{}))>> >=
        1;  // not an attached storage

template <typename T, typename I>
concept attribute_of =
    attribute<T> && item<I> &&
    must::contains<T, decltype(siconos::storage::pattern::attributes(I{}))>;
}  // namespace match
namespace types {
template <template <typename T> typename Transform, typename... Args>
using transform = decltype(transform(
    []<typename A>(A) { return Transform<A>{}; }, gather<Args...>{}));
template <match::attribute... Attrs>
using attributes = gather<Attrs...>;

template <match::item... Items>
using properties_of_items =
    decltype(mp::concat_all(typename Items::properties{}...));

}  // namespace types

template <match::item... Items>
using attributes_of_items = decltype(mp::concat_all(
    siconos::storage::pattern::attributes(Items{})...));

template <string_literal S>
struct indice_value : symbol<S> {
  std::size_t value;
};

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

template <match::attribute Attr>
static auto item_attribute = [](auto items) constexpr {
  using items_t = std::decay_t<decltype(items)>;

  auto loop = rec([]<typename Tpl>(auto&& loop, Tpl tpl) {
    using tpl_t = std::decay_t<Tpl>;
    using item_t = std::decay_t<decltype(car(tpl))>;

    if constexpr (match::attribute_of<Attr, item_t>) {
      return item_t{};
    }
    else if constexpr (match::attached_storage<Attr, item_t>) {
      return item_t{};
    }
    else if constexpr (mp::size(tpl_t{}) > mp::size_c<1>) {
      return loop(cdr(tpl));
    }
    else {
      []<typename Attribute = Attr, typename LastItem = item_t,
         typename Items = items_t, bool flag = false>() {
        static_assert(flag, "item not found");
      }();
    }
  });

  return loop(items);
};

static auto constexpr attribute_name(match::attribute auto a)
{
  return a.str.value;
};

template <typename... Ts>
std::tuple<Ts...> collect(Ts... ts)
{
  return std::make_tuple(ts...);
};

template <typename F>
auto method(auto s, F f, auto doc)
{
  return std::make_tuple(s, f, doc);
};

template <typename F>
auto method(auto s, F f)
{
  return std::make_tuple(s, f);
};

template <typename Tpl1, typename Tpl2>
auto concat_methods(Tpl1 m1, Tpl2 m2)
{
  return std::tuple_cat(m1, m2);
}

namespace match {
template <typename T>
concept def_method = requires { typename T::def_method_t; };

template <typename T>
concept methods = requires(T h) { h.methods(); };
}  // namespace match

static auto constexpr method_name(auto m) { return std::get<0>(m); }

static auto constexpr method_def(auto m) { return std::get<1>(m); }

namespace match {
template <typename T>
concept npy_format = (requires { typename T::value_type; } &&
                      std::is_scalar_v<typename T::value_type>) ||
                     requires { typename T; };

template <typename D>
concept store = requires(D d) { d.store(); };

template <typename Item>
concept without_attributes_bindings =
    requires { typename Item::without_attributes_bindings; };

template <typename Item>
concept without_attached_storages_bindings =
    requires { typename Item::without_attached_storages_bindings; };

template <typename Item>
concept batch_capable = requires { typename Item::batch_capable; };

template <typename Prop>
concept without_binding = requires { typename Prop::without_binding_t; };

template <typename Item>
concept empty_item = requires { typename Item::empty_item_t; };
}  // namespace match

struct empty_item : item {
  using empty_item_t = void;
};

}  // namespace siconos::storage::pattern
