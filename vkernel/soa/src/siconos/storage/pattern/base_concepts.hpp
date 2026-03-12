#pragma once

#include <functional>
#include <iostream>
#include <tuple>

namespace siconos::storage::pattern {

template <typename Attrs>
static auto tags = []() constexpr {
  return transform(
      []<typename Attr>(Attr) { return std::tuple<typename Attr::tag>{}; },
      Attrs{});
};

namespace concepts {
// https://stackoverflow.com/questions/68443804/c20-concept-to-check-tuple-like-types
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2165r2.pdf
template <class T, std::size_t N>
concept has_tuple_element = requires(T t) {
  typename std::tuple_element_t<N, std::remove_const_t<T>>;
  {
    std::get<N>(t)
  } -> std::convertible_to<const std::tuple_element_t<N, T>&>;
};

template <class T>
concept tuple_like_ = !std::is_reference_v<T> && requires(T t) {
  typename std::tuple_size<T>::type;
  requires std::derived_from<
      std::tuple_size<T>,
      std::integral_constant<std::size_t, std::tuple_size_v<T>>>;
} && []<std::size_t... N>(std::index_sequence<N...>) {
  return (has_tuple_element<T, N> && ...);
}(std::make_index_sequence<std::tuple_size_v<T>>());

template <class T>
concept tuple_like = tuple_like_<std::decay_t<T>>;

template <typename T>
concept array_like = requires(T a) {
  typename T::size_type;
  typename T::value_type;
  {
    std::size(a)
  } -> std::convertible_to<std::size_t>;
} && []<std::size_t... N>(std::index_sequence<N...>) {
  return (has_tuple_element<T, N> && ...);
}(std::make_index_sequence<std::tuple_size_v<T>>());

template <typename T>
concept keep = requires(T t) {
  typename T::tag;
  {
    t.size
  };
};

template <typename T>
concept use = requires { typename T::use_t; };

template <typename T>
concept use_items = requires { typename T::use_items; };

template <typename T>
concept tag = requires {
  typename T::tag;
  typename T::tag_t;
};

template <typename T>
concept structure = requires { typename T::structure_t; };

template <typename T>
concept attribute = requires {
  typename T::tag;
  typename T::symbol;
};

template <typename T>
concept has_attributes = requires { typename T::attributes; };

template <typename T>
concept terminal_attribute =
    attribute<T> && requires { typename T::structure::type; };

template <typename... Ts>
concept at_least_one_terminal_attribute = (terminal_attribute<Ts> || ...);

template <typename T>
concept has_at_least_one_terminal_attribute =
    has_attributes<T> &&
    at_least_one_terminal_attribute<typename T::attributes>;

template <typename T>
concept item = has_attributes<typename T::definition>;

template <typename T>
concept vertex_item = requires {
  {
    T::definition::vertex_item == true
  };
};

}  // namespace concepts

}  // namespace siconos::pattern
