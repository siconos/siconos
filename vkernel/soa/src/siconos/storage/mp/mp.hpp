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

#include <boost/hana/fwd/type.hpp>
#define BOOST_HANA_CONFIG_ENABLE_STRING_UDL 1
#include <algorithm>
#include <array>
#include <boost/hana.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <boost/hana/functional/overload_linearly.hpp>
#include <boost/hana/fwd/find_if.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/fwd/for_each.hpp>
#include <boost/hana/fwd/map.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/lazy.hpp>
#include <boost/hana/not_equal.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana/string.hpp>
#include <boost/type_index.hpp>
#include <cstdint>
#include <numeric>
#include <source_location>
#include <string>
#include <type_traits>

// #include "boost/hana/ext/boost/ublas.hpp"
//  cf
//  https://www.boost.org/doc/libs/1_80_0/libs/hana/doc/html/structboost_1_1hana_1_1string.html#ad77f7afff008c2ce15739ad16a8bf0a8

#include <string_view>

using namespace boost::hana::literals;

namespace siconos::storage::mp {

// debug (see_below for clang, gcc above...)
#if defined(__clang__)
template <typename... Ts>
constexpr bool type_trace() {
  std::tuple<Ts...> see_messages_below;
  return false;
};
#else
template <typename... Ts>
constexpr bool type_trace() {
  std::tuple<Ts...> see_messages_above;
  return false;
};
#endif

namespace hana = boost::hana;

namespace detail {

constexpr bool is_identifier_char(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') ||
         c == '_';
}

// the compiler's spelling of the function, with T in it:
// "... [with T = name; ...]" (gcc) or "... [T = name]" (clang)
template <typename T>
constexpr std::string_view pretty_function() {
  return std::source_location::current().function_name();
}

template <typename T>
constexpr std::string make_type_name() {
  const std::string_view f = pretty_function<T>();
  const auto b = f.find("T = ") + 4;
  const std::string_view full = f.substr(b, f.find_first_of(";]", b) - b);

  // without the namespaces: an identifier followed by :: is dropped
  std::string name;
  std::string word;
  for (std::size_t k = 0; k < full.size(); ++k) {
    const char c = full[k];
    if (is_identifier_char(c)) {
      word += c;
    } else if (c == ':' && k + 1 < full.size() && full[k + 1] == ':') {
      word.clear();
      ++k;
    } else {
      name += word;
      word.clear();
      name += c;
    }
  }
  name += word;

  // string_literal<N>{"x"} -> "x"
  for (auto p = name.find("string_literal<"); p != std::string::npos;
       p = name.find("string_literal<", p)) {
    const auto open = name.find('{', p);
    const auto close = name.find('}', open);
    name =
        name.substr(0, p) + name.substr(open + 1, close - open - 1) + name.substr(close + 1);
  }

  // gcc spells nested closing brackets "> >", clang ">>"
  for (auto p = name.find("> >"); p != std::string::npos; p = name.find("> >", p)) {
    name.erase(p + 1, 1);
  }
  return name;
}

// the characters of the name, computed once at compile time: a constexpr
// std::string cannot outlive the constant evaluation
template <typename T>
inline constexpr auto type_name_chars = [] {
  std::array<char, make_type_name<T>().size()> a{};
  std::ranges::copy(make_type_name<T>(), a.begin());
  return a;
}();

}  // namespace detail

// the name of a type, without the namespaces, the string literals of its
// template arguments as strings
template <typename T>
constexpr std::string_view type_name() {
  return {detail::type_name_chars<T>.data(), detail::type_name_chars<T>.size()};
}

using hana::append;
using hana::at;
using hana::bool_c;
using hana::equal;
using hana::flatten;
using hana::index_if;
using hana::is_valid;
using hana::pair;
using hana::prepend;
using hana::size;
using hana::size_c;
using hana::to;
using hana::to_set;
using hana::tuple;
using hana::unique;
template <typename... Ts, typename F>
constexpr auto unpack(hana::tuple<Ts...>, F&& f) {
  return f.template operator()<Ts...>(Ts{}...);
}

template <typename T, typename F>
constexpr auto unpack(T&& t, F&& f) {
  return hana::unpack(std::forward<T>(t), std::forward<F>(f));
}
static constexpr auto typeid_ = hana::typeid_;

template <template <typename... Ts> typename F>
static constexpr auto trait = hana::trait<F>;

using hana::compose;

static constexpr auto any_of = hana::any_of;

static constexpr auto make_tuple = hana::make_tuple;

using hana::range_c;

template <std::size_t N>
static constexpr auto range = hana::range_c<std::size_t, 0, N>;

static auto fold_left = []<typename Array, typename State, typename Fun>(
                            Array&& array, State&& initial_state,
                            Fun&& fun) constexpr -> decltype(auto) {
  using array_type = std::decay_t<Array>;

  /* ~ static */
  if constexpr (hana::Foldable<array_type>::value) {
    return hana::fold_left(array, initial_state, fun);
  }
  /* ~ dynamic (aka std::vector) */
  else if constexpr (std::copy_constructible<array_type> &&
                     std::equality_comparable<decltype(array.begin())> &&
                     std::input_iterator<decltype(array.begin())>) {
    return std::accumulate(array.begin(), array.end(), initial_state, fun);
  } else {
    // cf
    // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
    []<bool flag = false>() {
      static_assert(flag, "cannot fold_left with these parameters");
    }();
  }
};

static auto overload = hana::overload_linearly;

using hana::find_if;

using hana::for_each;

using hana::insert;

static auto first = hana::first;

static auto second = hana::second;

decltype(auto) tuple_first(auto&& tpl) { return tpl[0_c]; }
decltype(auto) tuple_second(auto&& tpl) { return tpl[1_c]; }
decltype(auto) tuple_third(auto&& tpl) { return tpl[2_c]; }

static auto transform = hana::transform;

static auto concat = hana::concat;

namespace detail {
// Concatenation of tuple types in one pass.  The fold of hana::concat below
// rebuilds a tuple at every step, which is quadratic in compile time over the
// item and property lists of a whole world.
template <typename... Tuples>
struct concat_types;
template <>
struct concat_types<> {
  using type = hana::tuple<>;
};
template <typename... Ts>
struct concat_types<hana::tuple<Ts...>> {
  using type = hana::tuple<Ts...>;
};
template <typename... As, typename... Bs, typename... Rest>
struct concat_types<hana::tuple<As...>, hana::tuple<Bs...>, Rest...>
    : concat_types<hana::tuple<As..., Bs...>, Rest...> {};

// tuples of empty types only (type lists): no value is lost by rebuilding them
template <typename T>
constexpr bool is_type_list = false;
template <typename... Ts>
constexpr bool is_type_list<hana::tuple<Ts...>> =
    ((std::is_empty_v<Ts> && std::is_default_constructible_v<Ts>) && ...);
}  // namespace detail

template <typename... Args>
auto constexpr concat_all(Args... args) {
  if constexpr ((detail::is_type_list<Args> && ...)) {
    return typename detail::concat_types<Args...>::type{};
  } else {
    return hana::fold_left(hana::make_tuple(args...), hana::make_tuple(), hana::concat);
  }
}

// static_assert(concat(make_tuple(1, 2, 3), make_tuple(4, 5, 6)) ==
//               make_tuple(1, 2, 3, 4, 5, 6));
// static_assert(concat_all(make_tuple(1), make_tuple(2), make_tuple(3, 4, 5))
// ==
//               make_tuple(1, 2, 3, 4, 5));

// using hana::pair;

using hana::type_c;

template <typename T>
static constexpr auto key = hana::type_c<T>;

template <typename First, typename Second>
using key_value = hana::pair<std::decay_t<decltype(key<First>)>, Second>;

using hana::make_pair;

template <typename... Pairs>
using map = hana::map<Pairs...>;

template <typename T, typename D>
static constexpr decltype(auto) get_m(D&& data) {
  return static_cast<D&&>(data)[key<T>];
};

// A store that resolves its own columns (storage::world): no key pack to scan,
// the store hands back the column for T directly.
template <typename T, typename S>
  requires requires(S&& s) { static_cast<S&&>(s).template column<T>(); }
decltype(auto) get(S&& s) {
  return static_cast<S&&>(s).template column<T>();
};

static auto make_type_c = []<typename T>(T) constexpr { return type_c<T>; };

static constexpr auto all_type_c(auto tpl) { return transform(tpl, make_type_c); };

namespace detail {
// First occurrence of each type, in order: what hana::to_set then to_tuple
// gives, without inserting into a set one element at a time (quadratic in
// compile time).  The comparisons are __is_same builtins inside constant
// evaluation, so no trait is instantiated per pair.
template <typename T, typename... Us>
constexpr std::array<bool, sizeof...(Us)> same_as_each = {__is_same(T, Us)...};

template <typename... Ts>
struct unique_types {
  static constexpr std::size_t n = sizeof...(Ts);
  static constexpr std::array<std::array<bool, n>, n> same = {same_as_each<Ts, Ts...>...};

  static constexpr std::array<bool, n> first = [] {
    std::array<bool, n> f{};
    for (std::size_t i = 0; i < n; ++i) {
      f[i] = true;
      for (std::size_t j = 0; j < i; ++j) {
        if (same[i][j]) {
          f[i] = false;
          break;
        }
      }
    }
    return f;
  }();

  static constexpr std::size_t count = [] {
    std::size_t c = 0;
    for (bool b : first) c += b;
    return c;
  }();

  static constexpr std::array<std::size_t, count> indices = [] {
    std::array<std::size_t, count> a{};
    for (std::size_t i = 0, c = 0; i < n; ++i) {
      if (first[i]) a[c++] = i;
    }
    return a;
  }();

  template <std::size_t... I>
  static auto build(std::index_sequence<I...>)
      -> hana::tuple<__type_pack_element<indices[I], Ts...>...>;

  using type = decltype(build(std::make_index_sequence<count>{}));
};

template <>
struct unique_types<> {
  using type = hana::tuple<>;
};
}  // namespace detail

template <typename... Ts>
static constexpr auto tuple_unique(hana::tuple<Ts...>) {
  return typename detail::unique_types<std::decay_t<Ts>...>::type{};
};

template <typename Xs, typename X>
static constexpr bool contains(Xs xs, X) {
  return hana::unpack(
      xs, [](auto... xs) { return (std::is_same_v<X, std::decay_t<decltype(xs)>> || ...); });
}

using hana::filter;

template <auto F, typename... Ts>
using check =
    std::conditional_t<F.template operator()<Ts...>(), std::true_type, std::false_type>;

template <auto F, typename... Ts>
struct on_concept {
  template <typename T2>
  struct is_a_model {
    using type = is_a_model<T2>;
    static constexpr bool value = check<F, T2, Ts...>::value;
  };
};

template <auto F, typename... Ts>
static constexpr auto is_a_model =
    compose(trait<on_concept<F, Ts...>::template is_a_model>, typeid_);

template <typename B>

static constexpr auto derive_from =
    is_a_model<[]<typename T>() consteval { return std::derived_from<T, B>; }>;

template <typename B>
static constexpr auto is_inside_type_parent =
    is_a_model<[]<typename T>() consteval { return std::derived_from<B, typename T::type>; }>;

template <typename... Ts>
constexpr auto std_tuple(const hana::tuple<Ts...>& htpl) {
  return hana::unpack(htpl,
                      []<typename... Elems>(Elems...) { return std::make_tuple(Elems{}...); });
}

// https://stackoverflow.com/questions/18063451/get-index-of-a-tuple-elements-type
template <class T, class... Ts>
constexpr std::size_t index_of(const std::tuple<Ts...>&) {
  int found{}, count{};
  ((!found ? (++count, found = std::is_same_v<T, Ts>) : 0), ...);
  return found ? count - 1 : count;
}

}  // namespace siconos::storage::mp
