#pragma once

#include <cctype>
#include <cstddef>
#include <numeric>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
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

namespace detail {

// the name of a type, without the namespaces, the string literals of its
// template arguments as strings
template <typename T>
std::string type_name() {
  std::string_view f = __PRETTY_FUNCTION__;
  const auto b = f.find("T = ") + 4;
  const std::string full(f.substr(b, f.find_first_of(";]", b) - b));

  // without the namespaces: an identifier followed by :: is dropped
  std::string name;
  std::string word;
  for (std::size_t k = 0; k < full.size(); ++k) {
    const char c = full[k];
    if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') {
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
  return name;
}

// a placeholder or an item: its label if it has one
template <typename T>
std::string label() {
  if constexpr (requires { T::label; }) {
    return std::string(T::label);
  } else {
    return type_name<T>();
  }
}

}  // namespace detail

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

// the field Name of the item referenced by Ref, a reference of the iterated
// item: a gather. The reference is an attribute, or a property when a
// declarer is given in the options; Name is an attribute of the referenced
// item, or else its property with the same declarer.
template <pattern::string_literal Ref, pattern::string_literal Name, typename... Options>
struct through {
  template <typename Refs, typename Values>
  struct gathered {
    Refs& refs;
    Values& values;
    decltype(auto) operator[](std::size_t i) { return (values[refs[i].value()]); }
    std::size_t size() const { return std::size(refs); }
  };

  template <typename Ctx, typename Target>
  static decltype(auto) target_column(auto& data, auto step) {
    if constexpr (requires { storage::attr_values<Target, Name>(data, step); }) {
      return storage::attr_values<Target, Name>(data, step);
    } else {
      using declarer = typename declarer_of<Options...>::type;
      return storage::prop_values<Target, Name, resolve_t<Ctx, declarer>>(data, step);
    }
  }

  template <typename Ctx, typename Item>
  static auto column(auto& data, auto step) {
    using ref = field<Ref, Options...>;
    auto& refs = ref::template column<Ctx, Item>(data, step);
    using target = typename std::decay_t<decltype(refs[0])>::type;
    auto& values = target_column<Ctx, target>(data, step + ref::offset);
    return gathered<std::remove_reference_t<decltype(refs)>,
                    std::remove_reference_t<decltype(values)>>{refs, values};
  }
};

// the same, written: a scatter (several items may write the same target, the
// kernel must write the same value or be run in order)
template <pattern::string_literal Ref, pattern::string_literal Name, typename... Options>
struct out_through : through<Ref, Name, Options...> {};

// a handle on the item referenced by Ref, for the code that works on handles
// (the relations): host only
template <pattern::string_literal Ref, typename... Options>
struct referenced {
  template <typename Data, typename Refs>
  struct handles {
    Data& data;
    Refs& refs;
    auto operator[](std::size_t i) { return storage::make_handle(data, refs[i]); }
    std::size_t size() const { return std::size(refs); }
  };

  template <typename Ctx, typename Item>
  static auto column(auto& data, auto step) {
    auto& refs = field<Ref, Options...>::template column<Ctx, Item>(data, step);
    return handles<std::remove_reference_t<decltype(data)>,
                   std::remove_reference_t<decltype(refs)>>{data, refs};
  }
};

// a kernel that reads more than its fields (a cache shared by all the items,
// a whole column...): it is not fused with other kernels
struct reads_world {};

// a kernel that needs host objects (factorizations, handles, variants...):
// it never runs on a device
struct host_only {};

// a kernel whose items must be visited in order (side effects through its
// parameters, counters...): it never runs in parallel
struct ordered {};

// a kernel that writes into an assembled matrix or vector given in its
// parameters, each item in its own rows
struct assembles {};

// the value an iteration space gives with each item: the value of the item
// in a sparse set (set), its entry in a property (entries)
struct space_value {
  template <typename Ctx, typename Item>
  static space_value column(auto&, auto) {
    return {};
  }
};
struct entry : space_value {};

// the index of the item
struct item_index {
  template <typename Ctx, typename Item>
  static item_index column(auto&, auto) {
    return {};
  }
};

// the field of item i, v being the value given by the iteration space
template <typename Column>
decltype(auto) at(Column& column, auto i, auto&& v) {
  if constexpr (std::is_same_v<std::decay_t<Column>, space_value>) {
    return v;
  } else if constexpr (std::is_same_v<std::decay_t<Column>, item_index>) {
    return i;
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
  static std::string text() { return "all " + detail::label<Item>(); }
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
  static std::string text() {
    return detail::label<Item>() + " in " + detail::type_name<Tag>();
  }
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
  static std::string text() {
    return detail::label<Item>() + " with " + std::string(Name.value);
  }
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
  static_assert(sizeof...(Kernels) == 1 || (!std::is_base_of_v<reads_world, Kernels> && ...),
                "a kernel that reads the world cannot be fused");
};

/// the exclusive prefix sums of a column of counts or flags: offsets[k] is
/// the sum of the k first values, offsets[n] the total. For flags, the rank
/// of each flagged item among the flagged ones (a scan, the compaction of a
/// mark).
template <typename Index = std::size_t>
std::vector<Index> offsets(const auto& column) {
  std::vector<Index> o(std::size(column) + 1, 0);
  for (std::size_t k = 0; k < std::size(column); ++k) {
    o[k + 1] = o[k] + static_cast<Index>(column[k]);
  }
  return o;
}

/// v[k] = v[0] + ... + v[k], in place: with v[0] = 0 and v[k + 1] the size
/// of the output of item k, the offsets of the outputs of variable size
void inclusive_scan(auto& v) {
  std::inclusive_scan(std::begin(v), std::end(v), std::begin(v));
}

/// the number of items of an iteration space with a value (set, entries...)
template <typename Space, typename Ctx>
std::size_t count(auto& data, auto step) {
  std::size_t n = 0;
  Space::template for_each<Ctx>(data, step, 0, [&n](auto, auto&&) { ++n; });
  return n;
}

// what a kernel does, as text: its iteration space, the fields it reads and
// writes, its markers (see info<K>() and the dataflow observer of
// storage/schedule.hpp)
struct kernel_info {
  std::string name;
  std::string space;
  std::string item;  // the type of the items it ran on
  std::vector<std::string> reads;
  std::vector<std::string> writes;
  bool host_only = false;
  bool ordered = false;
  bool reads_world = false;
  bool assembles = false;
};

namespace detail {

template <typename S>
std::string space_text() {
  if constexpr (requires { S::text(); }) {
    return S::text();
  } else {
    return type_name<S>();
  }
}

template <typename Field>
struct access {
  static void add(kernel_info&) {}  // space values, item index
};

template <pattern::string_literal Name, typename... Options>
std::string field_label() {
  return std::string(Name.value) + (field<Name, Options...>::offset ? "_next" : "");
}

template <pattern::string_literal Name, typename... Options>
struct access<in<Name, Options...>> {
  static void add(kernel_info& i) { i.reads.push_back(field_label<Name, Options...>()); }
};
template <pattern::string_literal Name, typename... Options>
struct access<out<Name, Options...>> {
  static void add(kernel_info& i) { i.writes.push_back(field_label<Name, Options...>()); }
};
template <pattern::string_literal Name, typename... Options>
struct access<inout<Name, Options...>> {
  static void add(kernel_info& i) {
    i.reads.push_back(field_label<Name, Options...>());
    i.writes.push_back(field_label<Name, Options...>());
  }
};
template <pattern::string_literal Ref, pattern::string_literal Name, typename... Options>
struct access<through<Ref, Name, Options...>> {
  static void add(kernel_info& i) {
    i.reads.push_back(std::string(Ref.value) + "." + field_label<Name, Options...>());
  }
};
template <pattern::string_literal Ref, pattern::string_literal Name, typename... Options>
struct access<out_through<Ref, Name, Options...>> {
  static void add(kernel_info& i) {
    i.writes.push_back(std::string(Ref.value) + "." + field_label<Name, Options...>());
  }
};
template <pattern::string_literal Ref, typename... Options>
struct access<referenced<Ref, Options...>> {
  static void add(kernel_info& i) { i.reads.push_back(std::string(Ref.value) + " (handle)"); }
};

template <typename K, typename... Fields>
void add_accesses(kernel_info& i, std::tuple<Fields...>*) {
  (access<Fields>::add(i), ...);
}

}  // namespace detail

/// what the kernel K (or fused<...>) does
template <typename K>
kernel_info info() {
  kernel_info i;
  if constexpr (requires { typename K::fields; }) {
    i.name = detail::type_name<K>();
    i.space = detail::space_text<typename K::space>();
    detail::add_accesses<K>(i, static_cast<typename K::fields*>(nullptr));
    i.host_only = std::is_base_of_v<host_only, K>;
    i.ordered = std::is_base_of_v<ordered, K>;
    i.reads_world = std::is_base_of_v<reads_world, K>;
    i.assembles = std::is_base_of_v<assembles, K>;
  } else {
    [&]<typename... Kernels>(fused<Kernels...>*) {
      (
          [&](kernel_info k) {
            i.name += (i.name.empty() ? "" : " + ") + k.name;
            i.space = k.space;
            i.reads.insert(i.reads.end(), k.reads.begin(), k.reads.end());
            i.writes.insert(i.writes.end(), k.writes.begin(), k.writes.end());
            i.host_only = i.host_only || k.host_only;
            i.ordered = i.ordered || k.ordered;
            i.reads_world = i.reads_world || k.reads_world;
            i.assembles = i.assembles || k.assembles;
          }(info<Kernels>()),
          ...);
    }(static_cast<K*>(nullptr));
  }
  return i;
}

/// while not null, run appends the kernels it runs (see the dataflow
/// observer of storage/schedule.hpp)
inline thread_local std::vector<kernel_info>* recording = nullptr;

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
              using c_t = std::decay_t<decltype(c)>;
              if constexpr (!std::is_same_v<c_t, space_value> &&
                            !std::is_same_v<c_t, item_index>) {
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
void run(auto& data, auto step, auto&& p) {
  if constexpr (requires { typename K::fields; }) {
    run<fused<K>, Ctx>(data, step, p);
  } else {
    if (recording) {
      auto i = info<K>();
      i.item = detail::type_name<resolve_t<Ctx, typename K::space::item>>();
      recording->push_back(std::move(i));
    }
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
