#pragma once

#include <cstddef>
#include <memory>
#include <tuple>
#include <type_traits>

#include "siconos/storage/mp/mp.hpp"

namespace siconos::storage::pattern {
using size_t = std::size_t;

template <typename... Ts>
using pointer = std::shared_ptr<Ts...>;

template <size_t N, typename tpl>
using nth_t = std::decay_t<decltype(tpl{}[mp::size_c<N>])>;

template <typename... Args>
using gather = mp::tuple<Args...>;

// template <typename T, typename Tpl>
// using cons_t = std::decay_t<decltype(std::tuple_cat<std::tuple<T>, Tpl>)>;

// https://ctrpeach.io/posts/cpp20-string-literal-template-parameters/
// https://stackoverflow.com/questions/62266052/c20-string-literal-template-argument-working-example
// https://www.cppstories.com/2021/constexpr-vecstr-cpp20/

template <std::size_t N>
struct string_literal {
  constexpr string_literal(const char (&str)[N])
  {
    std::copy_n(str, N, value);
  }
  char value[N];
};

template <string_literal a>
decltype(auto) operator""_x()
{
  return a;
};

template <typename T>
constexpr auto make_string_literal(T&& t)
{
  return string_literal(std::forward<T>(t));
};

template <string_literal Text>
struct text {
  static constexpr auto str = Text;
};

struct any_symbol {};
template <string_literal Symbol>
struct symbol : text<Symbol>, any_symbol {};

template <string_literal Symbol>
constexpr auto symb = symbol<Symbol>{};

template <string_literal Symbol>
struct name : text<Symbol> {};

template <size_t N>
constexpr auto make_string_literal(const string_literal<N>& a)
{
  return a;
}

template <size_t N>
consteval auto make_symbol(const string_literal<N>& a)
{
  symbol s = a;
  return s;
};

template <string_literal Descr>
struct description : text<Descr> {};

namespace match {
template <typename T>
concept scalar = std::is_scalar_v<T>;

template <typename T>
concept indice =
    std::is_scalar_v<T> && requires(T i) { std::array<double, 1>{}[i]; };

template <typename T>
concept symbol = std::derived_from<T, any_symbol>;

template <typename T>
concept ublas_matrix = requires(T m) { m(0, 0); };

template <typename T>
concept ublas_sparse_matrix = requires(T m) {
  m(0, 0);
  m.index1_data();
  m.index2_data();
};

template <typename T>
concept degrees_of_freedom = requires { typename T::degrees_of_freedom_t; };

template <typename T>
concept attributes = requires { typename T::attributes; };

template <typename T>
concept vdescriptor = requires { typename T::vdescriptor_t; };

template <typename T>
concept item = requires { typename T::item_t; };

template <typename T>
concept items = item<T> && requires { typename T::items; };

template <typename T>
concept properties = item<T> && requires { typename T::properties; };

template <typename T>
concept polymorphic_type = requires { typename T::polymorphic; };

template <typename T>
concept size = requires(T a) {
  { std::size(a) };
};

template <typename T>
concept push_back = requires(T a) {
  { a.push_back(typename T::value_type{}) };
};

template <typename T, typename I>
concept handle = std::derived_from<typename T::type, I> &&
                 requires { typename T::handle_t; };

template <typename T>
concept wrap = item<T> && requires { typename T::wrap_t; };

template <typename T>
concept database = requires { typename T::database_t; };

}  // namespace match
}  // namespace siconos::storage::pattern
