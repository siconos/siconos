#pragma once

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "siconos/storage/get.hpp"
#include "siconos/storage/sparse_set.hpp"

namespace siconos::storage::kernels {

/* Kernels: the functions that run over the items of a storage (the systems
   of an entity component system, the items being the entities and their
   attributes and properties the components).

   A kernel is a struct that declares
     - its iteration space: all the items of a type, the dense side of a
       sparse set, the items with an entry in a property...
     - the fields it reads (in), writes (out) or both (inout), at the current
       step or at the next one (the memory of keep<attr, 2>),
   and a static apply(p, fields...) that works on the fields of one item, p
   being the parameters of the run. apply knows nothing about the storage:
   run<K> finds the columns and loops.

     struct position : kernel<all<body>, in<"q">, in<"velocity">,
                              in<"velocity", next>, out<"q", next>> {
       static void apply(const auto& p, const auto& q, const auto& v,
                         const auto& v_next, auto& q_next) {
         q_next = q + p.h * (p.theta * v_next + (1.0 - p.theta) * v);
       }
     };

   The arguments of apply follow the declared fields, in order.

   fused<A, B> runs A then B on each item, in one loop. It is legal when A
   and B have the same iteration space: B then reads what A wrote for the
   same item only.

   The types in a kernel may be placeholders, resolved by the context of the
   run: a placeholder P defines P::resolved_in<Ctx>. An application defines
   its own placeholders and its own iteration spaces (see simul/kernels.hpp). */

// a field at the next step (the step + 1 slot of a keep<attr, 2> memory)
struct next {};

template <typename Ctx, typename T>
struct resolve {
  using type = T;
};
template <typename Ctx, typename T>
  requires requires { typename T::template resolved_in<Ctx>; }
struct resolve<Ctx, T> {
  using type = typename T::template resolved_in<Ctx>;
};
template <typename Ctx, typename T>
using resolve_t = typename resolve<Ctx, T>::type;

// the declarer among the options of a field: the first one that is not next,
// void for an attribute of the item
template <typename... Options>
struct declarer_of {
  using type = void;
};
template <typename Option, typename... Options>
struct declarer_of<Option, Options...> {
  using type = std::conditional_t<std::is_same_v<Option, next>,
                                  typename declarer_of<Options...>::type, Option>;
};

// a field of the iterated item: an attribute, or a property when a declarer
// is given
template <pattern::string_literal Name, typename... Options>
struct field {
  static constexpr std::size_t offset = (std::is_same_v<Options, next> || ...) ? 1 : 0;
  using declarer = typename declarer_of<Options...>::type;

  template <typename Ctx, typename Item>
  static decltype(auto) column(auto& data, auto step) {
    if constexpr (std::is_void_v<declarer>) {
      return storage::attr_values<Item, Name>(data, step + offset);
    } else {
      return storage::prop_values<Item, Name, resolve_t<Ctx, declarer>>(data, step + offset);
    }
  }
};

template <pattern::string_literal Name, typename... Options>
struct in : field<Name, Options...> {};

template <pattern::string_literal Name, typename... Options>
struct out : field<Name, Options...> {};

template <pattern::string_literal Name, typename... Options>
struct inout : field<Name, Options...> {};

// the value an iteration space gives with each item: the value of the item
// in a sparse set (set), its entry in a property (entries)
struct space_value {
  template <typename Ctx, typename Item>
  static space_value column(auto&, auto) {
    return {};
  }
};
struct entry : space_value {};

// the field of item i, v being the value given by the iteration space
template <typename Column>
decltype(auto) at(Column& column, auto i, auto&& v) {
  if constexpr (std::is_same_v<std::decay_t<Column>, space_value>) {
    return v;
  } else {
    return (column[i]);
  }
}

// iteration spaces: Space::for_each<Ctx>(data, step, n, f) calls f(i, v) for
// every item i, v being the value of the space (or nothing useful), n being
// the size of the columns of the kernel

// every item
template <typename Item>
struct all {
  using item = Item;
  template <typename Ctx>
  static void for_each(auto&, auto, std::size_t n, auto&& f) {
    for (std::size_t i = 0; i < n; ++i) {
      f(i, space_value{});
    }
  }
};

// the items of the sparse set tagged by Tag, declared by Declarer, with
// their value in the set
template <typename Item, typename Tag, typename Declarer>
struct set {
  using item = Item;
  template <typename Ctx>
  static void for_each(auto& data, auto, std::size_t, auto&& f) {
    auto& values = storage::prop_dynamic_memory_tagged<resolve_t<Ctx, Item>, Tag,
                                                       resolve_t<Ctx, Declarer>>(data);
    storage::for_each_entry(values, [&f](auto i, auto& value) { f(i, value); });
  }
};

// the items with an entry in the property Name declared by Declarer, with
// this entry: only these items when the property is a sparse set
template <typename Item, pattern::string_literal Name, typename Declarer>
struct entries {
  using item = Item;
  template <typename Ctx>
  static void for_each(auto& data, auto step, std::size_t, auto&& f) {
    auto& values =
        storage::prop_values<resolve_t<Ctx, Item>, Name, resolve_t<Ctx, Declarer>>(data, step);
    storage::for_each_entry(values, [&f](auto i, auto& value) { f(i, value); });
  }
};

template <typename Space, typename... Fields>
struct kernel {
  using space = Space;
  using fields = std::tuple<Fields...>;
};

// kernels run in one loop
template <typename... Kernels>
struct fused {
  using space = std::tuple_element_t<0, std::tuple<typename Kernels::space...>>;
  static_assert((std::is_same_v<typename Kernels::space, space> && ...),
                "fused kernels must have the same iteration space");
};

namespace detail {

// references on the columns of the fields (space values held by value)
template <typename Ctx, typename Item, typename... Fields>
auto columns(auto& data, auto step, std::tuple<Fields...>*) {
  return std::tuple<decltype(Fields::template column<Ctx, Item>(data, step))...>(
      Fields::template column<Ctx, Item>(data, step)...);
}

// the number of items: the size of the first column that is not a space value
template <typename Columns>
std::size_t count(Columns& columns) {
  std::size_t n = 0;
  bool found = false;
  std::apply(
      [&](auto&... column) {
        (
            [&](auto& c) {
              if constexpr (!std::is_same_v<std::decay_t<decltype(c)>, space_value>) {
                if (!found) {
                  n = std::size(c);
                  found = true;
                }
              }
            }(column),
            ...);
      },
      columns);
  return n;
}

}  // namespace detail

/// run the kernel K (or fused<...>) on the data of step, with the
/// parameters p. Ctx resolves the placeholders.
template <typename K, typename Ctx>
void run(auto& data, auto step, const auto& p) {
  if constexpr (requires { typename K::fields; }) {
    run<fused<K>, Ctx>(data, step, p);
  } else {
    [&]<typename... Kernels>(fused<Kernels...>*) {
      using item = resolve_t<Ctx, typename K::space::item>;
      auto cols = std::make_tuple(detail::columns<Ctx, item>(
          data, step, static_cast<typename Kernels::fields*>(nullptr))...);
      const std::size_t n = detail::count(std::get<0>(cols));

      K::space::template for_each<Ctx>(data, step, n, [&](auto i, auto&& v) {
        [&]<std::size_t... k>(std::index_sequence<k...>) {
          (std::apply([&](auto&... column) { Kernels::apply(p, at(column, i, v)...); },
                      std::get<k>(cols)),
           ...);
        }(std::index_sequence_for<Kernels...>{});
      });
    }(static_cast<K*>(nullptr));
  }
}

}  // namespace siconos::storage::kernels
