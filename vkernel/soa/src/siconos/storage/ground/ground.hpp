#pragma once

#include <boost/hana/fwd/type.hpp>
#define BOOST_HANA_CONFIG_ENABLE_STRING_UDL 1
#include <algorithm>
#include <boost/hana.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/equal.hpp>
#ifdef __clang__
#include <boost/hana/experimental/type_name.hpp>
#endif
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
#include <numeric>
#include <type_traits>

// #include "boost/hana/ext/boost/ublas.hpp"
//  cf
//  https://www.boost.org/doc/libs/1_80_0/libs/hana/doc/html/structboost_1_1hana_1_1string.html#ad77f7afff008c2ce15739ad16a8bf0a8

#include <string_view>

using namespace boost::hana::literals;

namespace siconos::storage::ground {

template <typename Y>
constexpr decltype(auto) debug_type(Y)
{
  return boost::typeindex::type_id_with_cvr<Y>();
}

struct no_trace {};

template <typename T>
concept xdebug_type = std::derived_from<T, no_trace>;

// debug (see_below for clang, gcc above...)
#if defined(__clang__)
template <typename... Ts>
constexpr bool type_trace()
{
  std::tuple<Ts...> see_messages_below;
  return false;
};
#else
template <typename... Ts>
constexpr bool type_trace()
{
  std::tuple<Ts...> see_messages_above;
  return false;
};
#endif

namespace hana = boost::hana;

#ifdef __clang__
using hana::experimental::type_name;
#endif

using hana::take_while;

template <typename T>
using inner_type = typename std::remove_reference<T>::type;

using hana::append;
using hana::at;
using hana::equal;
using hana::eval;
using hana::flatten;
using hana::index_if;
using hana::integral_constant;
using hana::is_valid;
using hana::make_lazy;
using hana::pair;
using hana::prepend;
using hana::set;
using hana::size;
using hana::size_c;
using hana::to_map;
using hana::to_set;
using hana::tuple;
using hana::unique;
using hana::unpack;

static constexpr auto typeid_ = hana::typeid_;

template <template <typename... Ts> typename F>
static constexpr auto trait = hana::trait<F>;

using hana::compose;

static constexpr auto lockstep = hana::lockstep;

static constexpr auto any_of = hana::any_of;

static constexpr auto make_tuple = hana::make_tuple;

using hana::range_c;

template <std::size_t N>
static constexpr auto range = hana::range_c<std::size_t, 0, N>;

template <std::size_t N>
static constexpr auto iterate = hana::iterate<N>;

static auto fold_left = []<typename Array, typename State, typename Fun>(
                            Array &&array, State &&initial_state,
                            Fun &&fun) constexpr -> decltype(auto) {
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
  }
  else {
    // cf
    // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
    []<bool flag = false>() {
      static_assert(flag, "cannot fold_left with these parameters");
    }();
  }
};

static auto overload = hana::overload_linearly;

using hana::apply;

using hana::find_if;

using hana::for_each;

static auto insert = hana::insert;

static auto partial = hana::partial;

static auto first = hana::first;

static auto second = hana::second;

static auto front = hana::front;

static auto reverse = hana::reverse;

static auto scan_left = hana::scan_left;

static auto transform = hana::transform;

using hana::drop_front;
static auto concat = hana::concat;

template <typename... Args>
auto constexpr concat_all(Args... args)
{
  return hana::fold_left(hana::make_tuple(args...), hana::make_tuple(),
                         hana::concat);
}

// static_assert(concat(make_tuple(1, 2, 3), make_tuple(4, 5, 6)) ==
//               make_tuple(1, 2, 3, 4, 5, 6));
// static_assert(concat_all(make_tuple(1), make_tuple(2), make_tuple(3, 4, 5)) ==
//               make_tuple(1, 2, 3, 4, 5));

using hana::zip;

// f(T{}, ...) -> f<T>(...)
template <typename T>
static auto t_arg = []<typename F>(F &&f) { return ground::partial(f, T{}); };

// using hana::pair;

using hana::type_c;

template <typename T>
struct type_hash {
  using type = T;
  static constexpr int i{};
  static constexpr int const *value{&i};
};

template <typename T>
static constexpr auto type_hash_v = type_hash<T>::value;

template <auto hash>
struct hashed {};

template <typename T>
static constexpr auto key = hana::type_c<T>;

template <typename T>
static constexpr auto hashed_key = hana::type_c<hashed<type_hash_v<T>>>;

template <typename First, typename Second>
using key_value = hana::pair<std::decay_t<decltype(key<First>)>, Second>;

static auto make_key_value =
    []<typename First, typename Second>(
        First, Second &&second) constexpr -> decltype(auto) {
  return hana::make_pair(key<First>, static_cast<Second &&>(second));
};

static auto make_hashed_key_value =
    []<typename First, typename Second>(
        First, Second &&second) constexpr -> decltype(auto) {
  return hana::make_pair(hashed_key<typename First::type>,
                         static_cast<Second &&>(second));
};

template <typename Pair>
using hashed_pair_t =
    decltype(make_hashed_key_value(first(Pair{}), second(Pair{})));

using hana::make_pair;

template <typename... Pairs>
using pre_map = hana::tuple<Pairs...>;

template <typename... Pairs>
using map = hana::map<Pairs...>;

// The Hana map database takes a longer time to compile!
template <typename... Pairs>
struct database {
  using database_t = void;
  database() : store{} {};
  database(tuple<Pairs...> &&m)
      : store(to_map(static_cast<tuple<Pairs...> &&>(m))){};
  decltype(to_map(tuple<Pairs...>{})) store;
};

// Faster compilation but this error occurs with python bindings:
//    from ._nonos import __doc__, __version__, disks
// ImportError: vector::_M_realloc_insert
//
// template <typename... Pairs>
// struct database {
//   //  database() : store{} {};
//   //  database(tuple<Pairs...> &&m) : store(m){};
//   tuple<Pairs...> store;
// };

template <typename... Pairs>
auto to_database(tuple<Pairs...> &&data)
{
  return database<Pairs...>{static_cast<tuple<Pairs...> &&>(data)};
};

using hana::make_map;

template <typename Data, typename Key>
concept has_key = requires(Data m) { m[key<Key>]; };

template <typename T, typename D>
static constexpr decltype(auto) get_m(D &&data)
{
  return static_cast<D &&>(data)[key<T>];
};

template <typename T, typename D>
static constexpr decltype(auto) get_internal(D &&data)
{
  return static_cast<D &&>(data)[key<T>];
};

template <typename T, typename... Pairs>
static constexpr auto &get_internal(tuple<Pairs...> &&data)
{
  auto &&result =
      hana::find_if(static_cast<tuple<Pairs...> &&>(data),
                    []<typename P>(P) { return hana::first(P{}) == key<T>; });
  return hana::second(result.value());
};

template <typename T, typename... Pairs>
static constexpr auto &get_internal(tuple<Pairs...> &data)
{
  auto &&result = hana::find_if(
      data, []<typename P>(P) { return hana::first(P{}) == key<T>; });
  return hana::second(result.value());
};

template <typename T, typename... HPairs>
decltype(auto) get(ground::database<HPairs...> &&data)
{
  return get_internal<T>(
      static_cast<ground::database<HPairs...> &&>(data).store);
};

template <typename T, typename... HPairs>
decltype(auto) get(ground::database<HPairs...> &data)
{
  return get_internal<T>(data.store);
};

static auto make_type_c = []<typename T>(T) constexpr { return type_c<T>; };

static constexpr auto all_type_c(auto tpl)
{
  return transform(tpl, make_type_c);
};

static constexpr auto all_inside_types(auto tpl)
{
  return transform(tpl, []<typename T>(T) { return typename T::type{}; });
}

static constexpr auto tuple_unique(auto xs)
{
  return all_inside_types(hana::to_tuple(hana::to_set(all_type_c(xs))));
};

template <typename Xs, typename X>
static constexpr bool contains(Xs xs, X)
{
  return hana::contains(
      transform(xs, []<typename IX>(IX) { return type_c<IX>; }), type_c<X>);
}

using hana::filter;

static auto filter_t = []<typename Xs, typename Pred>(Xs &&xs, Pred &&pred) {
  return transform(hana::filter(all_type_c(static_cast<Xs &&>(xs)),
                                static_cast<Pred &&>(pred)),
                   []<typename Tc>(Tc) { return typename Tc::type{}; });
};

// map -> tuple -> tranform -> map
static auto map_transform = hana::demux(hana::to<hana::map_tag>)(
    compose(transform, hana::to<hana::tuple_tag>));

// dup(f)(x) = f(x, x)
static auto dup = []<typename F>(F &&f) constexpr -> decltype(auto) {
  return [&f]<typename X>(X &&x) {
    auto &&px = std::forward<X>(x);  // x must be forwarded once!!
    return std::forward<F>(f)(px, px);
  };
};

// static_assert(dup(hana::plus)(1) == 2);
// static_assert(dup(hana::mult)(2) == 4);

// map_transform pair(first, f(first, second)),
static auto map_value_transform =
    []<typename M, typename F>(M &&m, F &&f) constexpr -> decltype(auto) {
  return map_transform(
      std::forward<M>(m),
      dup(hana::lockstep(hana::make_pair)(
          hana::first, dup(hana::lockstep(std::forward<F>(f))(
                           hana::first, hana::second)))));
};

static auto pre_map_value_transform =
    []<typename M, typename F>(M &&m, F &&f) constexpr -> decltype(auto) {
  return transform(std::forward<M>(m),
                   dup(hana::lockstep(hana::make_pair)(
                       hana::first, dup(hana::lockstep(std::forward<F>(f))(
                                        hana::first, hana::second)))));
};

// compile-time itransform
static constexpr const auto itransform_ct(const auto &a, auto &&f)
{
  using array_type = std::decay_t<decltype(a)>;
  using size_type = typename array_type::size_type;
  array_type ta;  // 'a' passed as const ref is not a constant expression
  return [&f, &a]<size_type... I>(std::index_sequence<I...>) {
    return (array_type{f(I, a[I])...});
  }(std::make_integer_sequence<size_type, std::size(ta)>{});
}

static constexpr const auto itransform(const auto &array, auto &&func)
{
  using array_type = std::decay_t<decltype(array)>;
  if constexpr (hana::Foldable<array_type>::value) {
    return itransform_ct(array, func);
  }
  else if constexpr (std::equality_comparable<decltype(array.begin())> &&
                     std::input_iterator<decltype(array.begin())>) {
    array_type res;
    using size_type = typename array_type::size_type;
    std::transform(array.begin(), array.end(), std::back_inserter(res),
                   [&func, &res](const auto &x) {
                     size_type i = std::size(res);
                     return func(i, x);
                   });
    return res;
  }
  else {
    // cf
    // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
    []<bool flag = false>() {
      static_assert(flag, "cannot transform with these parameters");
    }();
  }
}
template <typename Base>
struct from {
  template <typename T>
  struct is_a_derivation {
    using type = is_a_derivation<T>;
    static constexpr bool value = std::derived_from<T, Base>;
  };
};

template <auto F, typename... Ts>
using check = std::conditional_t<F.template operator()<Ts...>(),
                                 std::true_type, std::false_type>;

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

static constexpr auto is_integral = is_a_model<[]<typename T>() consteval {
  return std::is_integral<T>::value;
}>;

template <typename B>
static constexpr auto derive_from = is_a_model<[]<typename T>() consteval {
  return std::derived_from<T, B>;
}>;

template <typename B>
static constexpr auto is_parent = is_a_model<[]<typename T>() consteval {
  return std::derived_from<B, T>;
}>;

template <typename B>
static constexpr auto is_inside_type_parent =
    is_a_model<[]<typename T>() consteval {
      return std::derived_from<B, typename T::type>;
    }>;

template <typename D>
static constexpr auto dump_keys(D, auto &&fun)
{
  for_each(D{}, [&fun]<typename KeyValue>(KeyValue kv) {
    fun(debug_type(inner_type<inner_type<decltype(first(kv))>>{})
            .pretty_name());
  });
}

template <auto I, typename R, typename F>
static constexpr R call_with_integral_constant_if_valid(R &&def_val, F &&fun)
{
  constexpr auto N = std::integral_constant<decltype(I), I>{};

  if constexpr (is_valid([](auto &&K) -> decltype(std::declval<F &&>()(K)) {
                })(N)) {
    return [&]() { return static_cast<F &&>(fun)(N); }();
  }
  else {
    return def_val;
  }
}

template <auto NumOfCases, typename ReturnType, typename F>
inline constexpr ReturnType call_with_index(auto index, ReturnType &&def_val,
                                            F &&f)
{
  constexpr auto fun_tab = []<std::size_t... I>(std::index_sequence<I...>) {
    return std::array{
        siconos::storage::ground::call_with_integral_constant_if_valid<
            I, ReturnType, F>...};
  }(std::make_index_sequence<NumOfCases>{});

  return fun_tab[index](static_cast<ReturnType &&>(def_val),
                        static_cast<F &&>(f));
}

template <typename... Ts>
auto std_tuple(const hana::tuple<Ts...> &htpl)
{
  return hana::unpack(htpl, []<typename... Elems>(Elems...) {
    return std::make_tuple(Elems{}...);
  });
}

// https://stackoverflow.com/questions/18063451/get-index-of-a-tuple-elements-type
template <class T, class... Ts>
constexpr std::size_t index_of(const std::tuple<Ts...> &)
{
  int found{}, count{};
  ((!found ? (++count, found = std::is_same_v<T, Ts>) : 0), ...);
  return found ? count - 1 : count;
}

}  // namespace siconos::storage::ground
